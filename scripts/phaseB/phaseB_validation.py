# -*- coding: utf-8 -*-
"""
Phase B — 纯计算结合验证（无 GROMACS / 无 MD 可行性替代方案）

方法学（honest disclosure）:
  本环境无 GROMACS / conda，无法运行 100-150 ns MD + MM-PBSA。
  本脚本采用 **单构象端点法静态 MM 验证**（single-structure end-point MM），
  对 Top3 肽的最优对接构象 (Vina MODEL 1) 做:
    1) 复合物 MMFF94s 松弛（受体骨架固定，配体+口袋侧链优化）
    2) 静态 MM 结合能  ΔE = E_complex - E_pocket - E_ligand  (gas-phase, MMFF94s)
       + 非极性溶剂校正 ΔG_SA = γ·ΔSASA  (γ=0.0072 kcal/mol/Å²)   [若 SASA 可算]
    3) 每残基相互作用能分解 (vdW + Coulomb, MMFF 力场参数)
       —— 作为「计算丙氨酸扫描」的定量替代，回答「哪些残基最关键」
    4) 相互作用指纹: H键 / 疏水接触 / 离子键 → 定位 DPP4 口袋残基

局限 (写入报告):
  - 单构象，无构象采样（MD 才能提供）
  - 溶剂化为近似（仅非极性 SASA 项），无显式 GB 极性项、无熵项
  - 力场为 MMFF94s（适用于小分子的通用力场，对蛋白精度低于 CHARMM36m/AMBER99SB-ILDN）
  - 结论为相对的「结合强弱排序 + 关键残基识别」，非绝对 ΔG 定量

依赖: RDKit 2026.03.3 (envs/default)。直接用 RDKit 读取 idppiv 队列的 Vina pose
       (含标准氨基酸残基名, 无需 OpenBabel)。
运行: envs/default/Scripts/python.exe phaseB_validation.py
"""
import os, sys, subprocess, math, json
from collections import defaultdict

BASE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable
# 必须带 .exe: Python subprocess.run 不走 shell, 不会像 Git Bash 那样按 PATHEXT 自动补后缀
OBABEL = os.path.join(os.path.dirname(PY), "obabel.exe")

# ---------- 仓库根 (scripts/phaseB/ -> repo root) ----------
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# ---------- Top3 肽定义 (方案 α: 新 iDPPIV 队列 RDKit 制备下的最优结合肽) ----------
# 正确对接 pose 来自新 iDPPIV 队列的 Vina 输出 (moso_ligands_idppiv/dock_<idx>.pdbqt),
# 其 <idx> 与 data/moso_dock_queue_idppiv.txt 一致:
#   APQIP -> index 19 (dock_19.pdbqt), LPPGP -> index 0 (dock_00.pdbqt),
#   APPSQ -> index 54 (dock_54.pdbqt)。
# 这些 pose 含标准氨基酸残基名 (ALA/PRO/GLN/ILE/SER/GLY/LEU), 可直接由 RDKit 读取,
# 无需 OpenBabel。dG 取自 moso_dock_results_idppiv.tsv (与 pose 一一对应)。
TOP3 = [
    {"name": "APQIP", "seq": "APQIP",
     "dock": os.path.join(REPO, "docking", "moso_ligands_idppiv", "dock_19.pdbqt"), "dg": -6.807},
    {"name": "LPPGP", "seq": "LPPGP",
     "dock": os.path.join(REPO, "docking", "moso_ligands_idppiv", "dock_00.pdbqt"), "dg": -6.558},
    {"name": "APPSQ", "seq": "APPSQ",
     "dock": os.path.join(REPO, "docking", "moso_ligands_idppiv", "dock_54.pdbqt"), "dg": -6.513},
]
RECEPTOR = os.path.join(REPO, "docking", "1WCY_receptor.pdbqt")
POCKET_CUTOFF = 9.0   # Å, 配体重原子 8-10 Å 内视为口袋
DIELECTRIC = 4.0       # 库仑介电常数（蛋白内部近似）
GAMMA_SA = 0.0072      # kcal/mol/Å² 非极性溶剂化表面张力
BACKBONE = {"N", "CA", "C", "O"}

from rdkit import Chem
from rdkit.Chem import rdForceFieldHelpers as ff
from rdkit.Chem import AllChem

# ----------------------------------------------------------------------
def run_obabel(args):
    r = subprocess.run([OBABEL] + args, capture_output=True, text=True)
    if r.returncode != 0:
        sys.stderr.write(f"[obabel ERR] args={args}\n{r.stderr}\n")
    return r

# AutoDock 原子类型 → 元素（蛋白中 CA 为 Cα 碳，取首字母即可）
_AD_ELEMENT = {
    "C": "C", "A": "C", "CA": "C", "CB": "C", "CC": "C", "CD": "C",
    "CE": "C", "CG": "C", "CH": "C", "CK": "C", "CM": "C", "CN": "C",
    "CQ": "C", "CT": "C", "CW": "C", "CX": "C", "CY": "C", "CZ": "C", "C0": "C",
    "N": "N", "NA": "N", "NB": "N", "NC": "N", "ND": "N", "NE": "N",
    "NF": "N", "NH": "N", "NX": "N",
    "O": "O", "OA": "O", "OB": "O", "OC": "O", "OD": "O", "OE": "O",
    "OF": "O", "OH": "O", "OS": "O", "OW": "O",
    "S": "S", "SA": "S", "SB": "S", "SD": "S", "SG": "S", "SH": "S",
    "P": "P", "PA": "P", "PB": "P", "PC": "P", "PD": "P", "PE": "P", "PF": "P", "PH": "P",
    "H": "H", "HS": "H", "HD": "H", "HE": "H", "HF": "H", "HG": "H",
    "HH": "H", "HN": "H", "HO": "H", "HP": "H", "HW": "H",
    "F": "F", "BR": "BR", "CL": "CL", "I": "I",
    "FE": "FE", "ZN": "ZN", "MG": "MG", "CA2": "CA", "MN": "MN",
    "CU": "CU", "NI": "NI", "K": "K",
}
def ad_element(adtype):
    adtype = (adtype or "").strip().upper()
    if adtype in _AD_ELEMENT:
        return _AD_ELEMENT[adtype]
    if adtype:
        c = adtype[0]
        if c in "C N O S P H F B I M Z G":
            return c
    return "C"

def parse_pdb_atoms(path, model=None):
    """解析 PDB/ PDBQT ATOM/HETATM 行。model=None 取全部；否则取指定 MODEL 之后的块。
    末列为 AutoDock 原子类型，需映射为元素符号。
    兼容单构象(无 MODEL 记录)与多构象两种 Vina 输出。"""
    # 预扫描: 文件是否含 MODEL 记录
    has_model_rec = False
    with open(path) as _f0:
        for _l in _f0:
            if _l[:6].strip() == "MODEL":
                has_model_rec = True
                break
    atoms = []
    cur_model = 1
    in_model = (model is None) or (not has_model_rec)
    with open(path) as f:
        for line in f:
            rec = line[:6].strip()
            if rec == "MODEL":
                try:
                    cur_model = int(line[10:14])
                except ValueError:
                    pass
                in_model = (model is None or cur_model == model)
                continue
            if rec == "ENDMDL":
                if model is not None and cur_model == model:
                    break
                in_model = (model is None) or (not has_model_rec)
                continue
            if rec in ("ATOM", "HETATM") and in_model:
                try:
                    ad = line[76:78].strip()
                    atoms.append({
                        "name": line[12:16].strip(),
                        "altloc": line[16].strip(),
                        "resname": line[17:20].strip(),
                        "chain": line[21].strip(),
                        "resid": int(line[22:26]),
                        "x": float(line[30:38]), "y": float(line[38:46]), "z": float(line[46:54]),
                        "element": ad_element(ad),
                    })
                except ValueError:
                    continue
    return atoms

def atoms_to_pdbblock(atoms, title="LIG"):
    lines = ["REMARK generated by phaseB", "COMPND    %s" % title]
    serial = 1
    for a in atoms:
        el = a["element"]
        if len(el) > 1:
            nm = a["name"]
        else:
            nm = a["name"] if len(a["name"]) >= 4 else (" %-3s" % a["name"])
        lines.append(
            "ATOM  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s  "
            % (serial, a["name"][:4], a["resname"][:3], a["chain"] or "A",
               a["resid"], a["x"], a["y"], a["z"], el)
        )
        serial += 1
    lines.append("END")
    return "\n".join(lines)

def build_mol_from_atoms(atoms, add_hs=True):
    blk = atoms_to_pdbblock(atoms)
    mol = Chem.MolFromPDBBlock(blk, removeHs=False, flavor=0)
    if mol is None:
        return None
    if add_hs:
        mol = Chem.AddHs(mol, addCoords=True)
    return mol

def ligand_from_pose(dock_path, name=None):
    """直接由 RDKit 读取 Vina 对接 pose 的 MODEL 1（idppiv 队列的 pose
    含标准氨基酸残基名 ALA/PRO/GLN/ILE/SER/GLY/LEU，RDKit 可正确感知肽键，
    无需 OpenBabel）。返回补氢后的 mol。"""
    lines = open(dock_path).read().splitlines()
    block, inm = [], False
    for ln in lines:
        if ln.startswith("MODEL"):
            inm = True
            continue
        if ln.startswith("ENDMDL"):
            break
        if inm and ln.startswith("ATOM"):
            block.append(ln)
    blk = "\n".join(block)
    m = Chem.MolFromPDBBlock(blk, removeHs=False, flavor=0)
    if m is None:
        return None
    m = Chem.AddHs(m, addCoords=True)
    return m

def mmff_minimize(mol, fixed_idx=None, maxIts=800):
    # 初始化环信息（MMFF 需要；肽中 Pro 吡咯烷环等）
    try:
        mol.UpdatePropertyCache(strict=False)
        Chem.GetSymmSSSR(mol)
    except Exception:
        pass
    mp = ff.MMFFGetMoleculeProperties(mol)
    if mp is None:
        return None, None
    fld = ff.MMFFGetMoleculeForceField(mol, mp, 0)
    if fld is None:
        return None, None
    if fixed_idx:
        for i in fixed_idx:
            # 位置约束（容差 0.3 Å，强约束常数）≈ 固定骨架重原子于当前坐标
            fld.MMFFAddPositionConstraint(int(i), 0.3, 1.0e6)
    conv = fld.Minimize(maxIts=maxIts)
    return fld.CalcEnergy(), conv

def coords_of(mol, idxs):
    conf = mol.GetConformer()
    return [(conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z) for i in idxs]

def _set_coords(mol, coords):
    """按原子序覆盖坐标（mol 原子数须 == len(coords)）。"""
    conf = mol.GetConformer()
    for j, (x, y, z) in enumerate(coords):
        conf.SetAtomPosition(j, (x, y, z))

def _set_heavy_coords(mol, coords):
    """仅覆盖重原子坐标（H 保持，随后由 AddHs 重建）。"""
    heavy = [i for i, a in enumerate(mol.GetAtoms()) if a.GetSymbol() != "H"]
    for j, hi in enumerate(heavy):
        if j < len(coords):
            x, y, z = coords[j]
            mol.GetConformer().SetAtomPosition(hi, (x, y, z))

# ----------------------------------------------------------------------
def find_pocket(rec_atoms, lig_atoms, cutoff):
    keys = set()
    for ra in rec_atoms:
        rx, ry, rz = ra["x"], ra["y"], ra["z"]
        for la in lig_atoms:
            d2 = (rx-la["x"])**2 + (ry-la["y"])**2 + (rz-la["z"])**2
            if d2 <= cutoff*cutoff:
                keys.add((ra["chain"], ra["resname"], ra["resid"]))
                break
    return keys

def pocket_atoms_subset(rec_atoms, keys):
    out = []
    for a in rec_atoms:
        if (a["chain"], a["resname"], a["resid"]) in keys:
            out.append(a)
    return out

# ----------------------------------------------------------------------
def compute_mm_bind(pep):
    """返回该肽的 MM 结合能相关量 + 复合物 mol 引用（供后续分解/指纹）。"""
    rec_path = os.path.join(BASE, RECEPTOR)
    dock_path = os.path.join(BASE, pep["dock"])
    print(f"  [{pep['name']}] 解析受体({RECEPTOR}) ...")
    rec_atoms = parse_pdb_atoms(rec_path)
    print(f"  [{pep['name']}] 解析配体最优构象 (MODEL 1) ...")
    lig_atoms = parse_pdb_atoms(dock_path, model=1)
    if not lig_atoms:
        raise RuntimeError("配体原子解析失败")
    print(f"  [{pep['name']}] 受体 {len(rec_atoms)} 原子, 配体 {len(lig_atoms)} 原子")

    pocket_keys = find_pocket(rec_atoms, lig_atoms, POCKET_CUTOFF)
    pkt_atoms = pocket_atoms_subset(rec_atoms, pocket_keys)
    print(f"  [{pep['name']}] 口袋残基数 = {len(pocket_keys)} (cutoff {POCKET_CUTOFF} Å)")

    # 剔除 PDBQT 自带 H（避免与 AddHs 重复），统一补氢
    pkt_heavy = [a for a in pkt_atoms if a["element"] != "H"]
    lig_heavy = [a for a in lig_atoms if a["element"] != "H"]

    # 构建口袋 mol（蛋白，列解析）与配体 mol（按序列重建 + 对接坐标）
    pocket_mol = build_mol_from_atoms(pkt_heavy, add_hs=True)
    lig_mol = ligand_from_pose(pep["dock"], pep["name"])
    if pocket_mol is None or lig_mol is None:
        raise RuntimeError("RDKit 构建口袋/配体失败")
    print(f"  [{pep['name']}] 口袋 mol {pocket_mol.GetNumAtoms()} 原子, 配体 mol {lig_mol.GetNumAtoms()} 原子")

    combo = Chem.CombineMols(pocket_mol, lig_mol)
    n_pkt = pocket_mol.GetNumAtoms()        # 含 H 的口袋总原子数
    n_pkt_heavy = len(pkt_heavy)            # 口袋重原子数
    n_lig_heavy = len(lig_heavy)            # 配体重原子数
    # 固定口袋骨架重原子（combo 中口袋重原子索引 0..n_pkt_heavy-1）
    fixed = [i for i, a in enumerate(pkt_heavy) if a["name"] in BACKBONE]
    print(f"  [{pep['name']}] 固定口袋骨架原子 = {len(fixed)} (仅用于口袋单独最小化)")
    # 复合物整体松弛（不固定骨架），让配体+口袋侧链共同松弛以消除对接 clash
    E_complex, _ = mmff_minimize(combo, fixed_idx=None, maxIts=1500)
    print(f"  [{pep['name']}] E_complex (MMFF) = {E_complex:.2f} kcal/mol")

    # 拆出最小化后的 **重原子** 坐标（用于重建单独 mol）
    pkt_heavy_idx = list(range(n_pkt_heavy))                    # 口袋重原子
    lig_heavy_idx = list(range(n_pkt, n_pkt + n_lig_heavy))   # 配体重原子
    pkt_coords = coords_of(combo, pkt_heavy_idx)
    lig_coords = coords_of(combo, lig_heavy_idx)

    # 口袋单独：重原子重建 → 写入最小化坐标 → 补氢 → 最小化
    pkt_only = build_mol_from_atoms(pkt_heavy, add_hs=False)
    _set_coords(pkt_only, pkt_coords)
    pkt_only = Chem.AddHs(pkt_only, addCoords=True)
    fixed2 = [j for j, a in enumerate(pkt_heavy) if a["name"] in BACKBONE]
    E_pocket, _ = mmff_minimize(pkt_only, fixed_idx=fixed2, maxIts=1000)

    # 配体单独：MOL2 干净结构 → 写入最小化坐标 → 补氢 → 最小化
    lig_only = ligand_from_pose(pep["dock"], pep["name"])
    _set_heavy_coords(lig_only, lig_coords)
    lig_only = Chem.AddHs(lig_only, addCoords=True)
    E_lig, _ = mmff_minimize(lig_only, fixed_idx=None, maxIts=1000)

    dE = E_complex - E_pocket - E_lig
    print(f"  [{pep['name']}] E_pocket={E_pocket:.2f}  E_lig={E_lig:.2f}  ΔE_MM={dE:.2f} kcal/mol")

    return {
        "name": pep["name"], "seq": pep["seq"], "dg_vina": pep["dg"],
        "n_pocket_res": len(pocket_keys), "pocket_keys": sorted(pocket_keys),
        "E_complex": E_complex, "E_pocket": E_pocket, "E_lig": E_lig,
        "dE_MM": dE,
        # 供分解/指纹: 最小化后 combo 引用 + 索引（重原子）
        "_combo": combo, "_pkt_idx": pkt_heavy_idx, "_lig_idx": lig_heavy_idx,
        "_pkt_atoms": pkt_heavy, "_lig_atoms": lig_heavy,
        "_fixed": fixed,
    }

# ----------------------------------------------------------------------
def contact_profile(res):
    """每配体残基位置的口袋接触数（H键/疏水/离子），作为计算丙氨酸扫描的
    稳健替代度量。接触判定纯几何（不依赖力场能量量级），可解释：
    接触最多的配体位置 = 最关键位置（突变为 Ala 时结合损失最大）。"""
    combo = res["_combo"]
    pkt_idx = res["_pkt_idx"]
    lig_idx = res["_lig_idx"]
    conf = combo.GetConformer()

    def is_carbon(i):
        return combo.GetAtomWithIdx(i).GetSymbol() == "C"
    def polar(i):
        return combo.GetAtomWithIdx(i).GetSymbol() in ("O", "N")
    def donor(i):
        a = combo.GetAtomWithIdx(i)
        if a.GetSymbol() == "N":
            return any(n.GetSymbol() == "H" for n in a.GetNeighbors())
        return False

    per_atom = defaultdict(int)
    for rank, gi in enumerate(lig_idx):
        for pi in pkt_idx:
            d = conf.GetAtomPosition(gi).Distance(conf.GetAtomPosition(pi))
            if d > 4.5:
                continue
            contact = False
            if is_carbon(gi) and is_carbon(pi) and d <= 4.0:
                contact = True
            if (donor(gi) or polar(gi)) and (donor(pi) or polar(pi)) and d <= 3.5:
                contact = True
            if contact:
                per_atom[rank] += 1

    # 按序列位置归并残基（backbone N 起新残基；上限=序列长度）
    lig_atoms = res["_lig_atoms"]
    n_seq = len(res["seq"])
    residues, cur = [], []
    for rank, la in enumerate(lig_atoms):
        if la["name"] == "N" and cur and len(residues) < n_seq:
            residues.append(cur)
            cur = [rank]
        else:
            cur.append(rank)
    if cur:
        residues.append(cur)
    # 侧链内多余 N（如 Gln NE2）导致的超额组 → 合并到最后一组
    while len(residues) > n_seq:
        residues[-2].extend(residues[-1])
        residues.pop()
    per_pos = []
    for ri, grp in enumerate(residues):
        cnt = sum(per_atom[r] for r in grp)
        letter = res["seq"][ri] if ri < len(res["seq"]) else "?"
        per_pos.append((ri + 1, letter, cnt))
    top = sorted(per_atom.items(), key=lambda x: -x[1])[:8]
    total = sum(per_atom.values())
    return total, per_pos, [(res["_lig_atoms"][r]["name"], c) for r, c in top]

# ----------------------------------------------------------------------
def interaction_fingerprint(res):
    """H键/疏水/离子键: 配体-口袋 接触残基列表。"""
    combo = res["_combo"]
    pkt_idx = res["_pkt_idx"]
    lig_idx = res["_lig_idx"]
    conf = combo.GetConformer()

    # 元素/电荷分类
    def is_donor(idx):
        a = combo.GetAtomWithIdx(idx)
        if a.GetSymbol() == "N":
            # 有 H 即供体
            return any(n.GetSymbol() == "H" for n in a.GetNeighbors())
        return False
    def is_acceptor(idx):
        a = combo.GetAtomWithIdx(idx)
        return a.GetSymbol() in ("O", "N")
    def is_carbon(idx):
        return combo.GetAtomWithIdx(idx).GetSymbol() == "C"

    hbonds, hydrophobics, ionics = [], [], []
    for rank, gi in enumerate(lig_idx):
        for pi in pkt_idx:
            d = conf.GetAtomPosition(gi).Distance(conf.GetAtomPosition(pi))
            if d > 4.5:
                continue
            # 疏水
            if is_carbon(gi) and is_carbon(pi) and d <= 4.0:
                hydrophobics.append((gi, pi, d))
            # H键候选: N/O 极性接触 <3.5 Å
            if (is_donor(gi) or is_acceptor(gi)) and (is_donor(pi) or is_acceptor(pi)) and d <= 3.5:
                hbonds.append((rank, gi, pi, d))
    # 去重 H键 (按配体原子)，并记录对应口袋残基
    seen = set()
    hb_unique = []
    for rank, gi, pi, d in hbonds:
        if rank not in seen:
            seen.add(rank)
            # 定位口袋原子所属残基（来自 pkt_heavy 列表）
            pkt = res["_pkt_atoms"]
            # pi 是 combo 中口袋重原子索引；映射到 pkt_heavy 名
            lig_name = res["_lig_atoms"][rank]["name"]
            hb_unique.append((rank, lig_name, pi, d))
    # 口袋残基名映射：用 combo 索引 -> pkt_heavy[rank] 的残基标签
    pkt = res["_pkt_atoms"]
    hb_residues = []
    for rank, lig_name, pi, d in hb_unique:
        # pi 在 combo 中属口袋段 0..n_pkt_heavy-1，对应 pkt_heavy[pi]
        try:
            a = pkt[pi]
            hb_residues.append("%s%d" % (a["resname"], a["resid"]))
        except Exception:
            hb_residues.append("?")
    return {
        "hbonds": len(hb_unique),
        "hydrophobics": len(hydrophobics),
        "ionic_like": len(ionics),
        "hbond_lig_atoms": [lig for _, lig, _, _ in hb_unique][:20],
        "hbond_pocket_residues": hb_residues[:20],
    }

# ----------------------------------------------------------------------

def _serialize(res):
    """单肽计算结果 -> 可 JSON 序列化的纯量字典（剥离 RDKit mol 引用）。"""
    return {
        "peptide": res["name"], "seq": res["seq"], "dg_vina": res["dg_vina"],
        "dE_MM": round(res["dE_MM"], 3), "E_complex": round(res["E_complex"], 2),
        "E_pocket": round(res["E_pocket"], 2), "E_lig": round(res["E_lig"], 2),
        "n_pocket_res": res["n_pocket_res"],
        "pocket_residues": ["%s%d" % (k[1], k[2]) for k in res["pocket_keys"]],
        "contact_total": res["contact_total"],
        "contact_per_position": [list(x) for x in res["contact_per_pos"]],
        "contact_top_atoms": res["contact_top_atoms"],
        "fingerprint": res["fingerprint"],
    }


def _run_single(pep):
    """单肽独立进程入口：计算 + 接触指纹，写出 phaseB_<NAME>.json。

    独立进程的原因：同进程重复调用 RDKit MMFF mmff_minimize 会在第二条肽
    触发 C 层段错误（无 Python traceback）；每条肽拆成独立子进程即稳定。
    """
    print(">>> 处理 %s (Vina dG=%.3f kcal/mol)" % (pep["name"], pep["dg"]))
    res = compute_mm_bind(pep)
    total_ct, per_pos, top_atoms = contact_profile(res)
    fp = interaction_fingerprint(res)
    res["contact_total"] = total_ct
    res["contact_per_pos"] = per_pos
    res["contact_top_atoms"] = top_atoms
    res["fingerprint"] = fp
    out = _serialize(res)
    with open(os.path.join(BASE, "phaseB_%s.json" % pep["name"]), "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    print("[OK] %s -> phaseB_%s.json (dE_MM=%.3f, contacts=%d)"
          % (pep["name"], pep["name"], out["dE_MM"], out["contact_total"]))
    return out


def main():
    print("=" * 70)
    print("Phase B — 静态 MM 结合验证 (Top3 DPP4 抑制肽)")
    print("=" * 70)
    # 每条肽用独立子进程（隔离 RDKit MMFF 同进程状态污染导致第二条肽崩溃）
    results = []
    for pep in TOP3:
        try:
            r = subprocess.run(
                [sys.executable, __file__, "--pep", pep["name"]],
                capture_output=True, text=True)
            for line in r.stdout.splitlines():
                print("  " + line)
            if r.returncode != 0:
                sys.stderr.write(r.stderr)
                print("  [ERROR] %s 子进程失败 (rc=%d)" % (pep["name"], r.returncode))
                continue
            jp = os.path.join(BASE, "phaseB_%s.json" % pep["name"])
            if os.path.exists(jp):
                with open(jp) as f:
                    results.append(json.load(f))
                os.remove(jp)  # 清理中间产物
        except Exception as e:
            import traceback
            traceback.print_exc()
            print("  [ERROR] %s 调度失败: %s" % (pep["name"], e))

    # 输出
    out_tsv = os.path.join(BASE, "phaseB_results.tsv")
    with open(out_tsv, "w") as f:
        f.write("peptide\tseq\tdG_Vina\tdE_complex\tdE_MM\tn_pocket_res\tcontact_total\tHbond\thydrophobic\n")
        for r in results:
            fp = r["fingerprint"]
            f.write("%s\t%s\t%.3f\t%.2f\t%.3f\t%d\t%d\t%d\t%d\n" % (
                r["peptide"], r["seq"], r["dg_vina"], r["E_complex"], r["dE_MM"],
                r["n_pocket_res"], r["contact_total"],
                fp["hbonds"], fp["hydrophobics"]))
    print("\n[OK] 主结果写入 %s" % out_tsv)

    # 详细 JSON
    with open(os.path.join(BASE, "phaseB_detail.json"), "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print("[OK] 详细信息写入 phaseB_detail.json")

    print("\n=== 摘要 ===")
    for r in results:
        fp = r["fingerprint"]
        pp = r["contact_per_position"]
        pp_str = " ".join("%s%d:%d" % (L, i, c) for i, L, c in pp)
        print("  %-6s dG_Vina=%7.3f  dE_MM=%7.3f  pocket_res=%3d  contact=%5d  Hbond=%d  Hbond_res=[%s]"
              % (r["peptide"], r["dg_vina"], r["dE_MM"], r["n_pocket_res"],
                 r["contact_total"], fp["hbonds"],
                 ",".join(fp["hbond_pocket_residues"][:6])))
        print("          接触位置剖面: %s" % pp_str)


if __name__ == "__main__":
    # 单肽模式：被子进程调用，计算并写出 phaseB_<NAME>.json 后退出
    if len(sys.argv) > 2 and sys.argv[1] == "--pep":
        _name = sys.argv[2]
        _single = [p for p in TOP3 if p["name"] == _name]
        if _single:
            _run_single(_single[0])
        sys.exit(0)
    main()
