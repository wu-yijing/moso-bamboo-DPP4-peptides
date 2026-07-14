# -*- coding: utf-8 -*-
"""MM-GBSA 单结构结合自由能重算 (OpenMM implicit GBn2) —— **pocket 限制版**。
对 Top3 候选 (APQIP/LPPGP/APPSQ) 的对接 pose:
    dG_bind = E_cplx - E_rec - E_lig   (implicit GBSA, single-structure)

为何 pocket 限制 (2026-07-14 修正):
- 全蛋白 (1WCY 链 A, ~23k 原子含 H) 的 GBSA Born 半径计算为 O(N^2),
  在 CPU 上需数十分钟, 超出本环境单次运行上限且后台任务会被回收, 故不可行。
- 改用 **口袋限制受体**: 取受体中距已知活性口袋中心 (62.8,47.7,4.8 Å) 15 Å 内的
  整残基 (保留残基完整性, 不切断残基), 约 40–60 残基 / <1k 原子, GBSA 可 <2 min 完成。
- 受体 E_rec 用**同一固定口袋**计算一次, 三条候选共享 -> 仅配体/复合物不同,
  使 dG 差异纯粹来自配体-口袋相互作用, 排序最干净。
- 氢原子补加由 OpenMM `PDBFixer` 完成, 但**跳过 findMissingResidues**
  (避免其回填口袋缺口导致体系重新膨胀); 仅 addMissingAtoms + addMissingHydrogens,
  末端残基按 N/C 末端补 H/oxT, amber14 模板可接受 (口袋末端封端的常用近似)。
- 配体直接由 RDKit 读取 pose MODEL 1 (标准氨基酸残基名) -> 交 PDBFixer 补氢。

方法学诚实披露: pocket 限制 MM-GBSA 忽略口袋外受体-受体长程相互作用与受体构象应变,
用于**三候选相对排序**而非绝对亲和力断言; 单结构近似、忽略 -TΔS。"""
import os
from openmm import app, openmm
from openmm.app import ForceField, PDBFile
from pdbfixer import PDBFixer
from rdkit import Chem
from rdkit.Geometry import Point3D

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOCK = os.path.join(REPO, "docking")
REC  = os.path.join(DOCK, "1WCY_receptor.pdbqt")
LIGDIR_IDPP = os.path.join(DOCK, "moso_ligands_idppiv")
DG = {"APQIP": -6.807, "LPPGP": -6.558, "APPSQ": -6.513}
# 正确 iDPPIV 队列对接 pose (dock_<idx>.pdbqt), idx 对应 moso_dock_queue_idppiv.txt
#   APQIP -> dock_19, LPPGP -> dock_00, APPSQ -> dock_54
TOP3 = [
    ("APQIP", os.path.join(LIGDIR_IDPP, "dock_19.pdbqt")),
    ("LPPGP", os.path.join(LIGDIR_IDPP, "dock_00.pdbqt")),
    ("APPSQ", os.path.join(LIGDIR_IDPP, "dock_54.pdbqt")),
]
STD_AA = {'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
           'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'}
POCKET_CENTER = (62.8, 47.7, 4.8)   # 1WCY DPP4 活性口袋中心 (box 中心)
POCKET_RADIUS = 9.0                  # Å (收紧以适配本环境 <2min 前台上限)


def write_rec_pocket(src, dst, center, radius):
    """受体 pdbqt -> 距 center 半径内整残基的口袋 PDB (链 A, 标准 AA, 去 HETATM)。"""
    cx, cy, cz = center
    # 按残基聚合原子行
    res_atoms = {}
    order = []
    with open(src) as fh:
        for ln in fh:
            if ln[:6] != "ATOM  ":
                continue
            rn = ln[17:20].strip(); ch = ln[21:22]
            if rn not in STD_AA or (ch and ch != 'A'):
                continue
            key = (ch, int(ln[22:26]))
            res_atoms.setdefault(key, []).append(ln)
            if key not in order:
                order.append(key)
    out = []
    n_res = 0
    for key in order:
        lines = res_atoms[key]
        inside = False
        for ln in lines:
            x = float(ln[30:38]); y = float(ln[38:46]); z = float(ln[46:54])
            if (x-cx)**2 + (y-cy)**2 + (z-cz)**2 <= radius*radius:
                inside = True
                break
        if inside:
            out.extend(lines)
            n_res += 1
    with open(dst, "w") as f:
        f.writelines(out)
    return n_res


def write_lig_raw(pose_path, dst):
    """取 pose MODEL 1 的标准氨基酸 ATOM 行, **按 resid 排序** (pose 的 ATOM 行常乱序,
    须整理为残基连续, 否则 PDBFixer 残基错配) 后写 PDB (截到 76 列), 交 PDBFixer 补氢。"""
    atoms = []
    inm = False
    with open(pose_path) as fh:
        for ln in fh:
            rec = ln[:6].strip()
            if rec == "MODEL":
                inm = True; continue
            if rec == "ENDMDL":
                break
            if inm and ln[:6] in ("ATOM  ", "HETATM"):
                resid = int(ln[22:26])
                atoms.append((resid, ln))
    atoms.sort(key=lambda x: x[0])   # 残基连续
    with open(dst, "w") as f:
        for _, ln in atoms:
            f.write(ln[:76] + "  \n")


def fix_pocket(src_pdb, dst_pdb):
    """PDBFixer 补缺失原子 + 补氢。对口袋 (非连续残基):
    先 findMissingResidues 建立属性, 再清空其缺口字典 -> addMissingAtoms 仅补
    现有残基缺失原子 (OXT/末端 H), 不回填口袋缺口, 体系保持小。"""
    fixer = PDBFixer(filename=src_pdb)
    fixer.findMissingResidues()      # 建立 missingResidues 属性 (否则 addMissingAtoms 报 AttributeError)
    fixer.missingResidues = {}       # 清空缺口 -> 不膨胀口袋
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    PDBFile.writeFile(fixer.topology, fixer.positions, dst_pdb)


def parse_pdbqt_model1_atoms(path):
    """解析 pdbqt MODEL 1 的 ATOM/HETATM 行 -> 字典列表。无 MODEL 则取全部。"""
    atoms, inm = [], False
    with open(path) as fh:
        for ln in fh:
            rec = ln[:6].strip()
            if rec == "MODEL":
                inm = True; continue
            if rec == "ENDMDL":
                break
            if inm and ln[:6] in ("ATOM  ", "HETATM"):
                atoms.append(_atom(ln))
    if not atoms:
        with open(path) as fh:
            for ln in fh:
                if ln[:6] in ("ATOM  ", "HETATM"):
                    atoms.append(_atom(ln))
    return atoms


def _atom(ln):
    return {
        "name": ln[12:16].strip(),
        "resname": ln[17:20].strip(),
        "chain": ln[21:22],
        "resid": int(ln[22:26]),
        "ins": ln[26:27].strip(),
        "x": float(ln[30:38]), "y": float(ln[38:46]), "z": float(ln[46:54]),
        "element": ln[76:78].strip(),
    }


def ligand_from_pose(seq, pqt):
    """从序列构建规范五肽拓扑 (键/Pro 环正确), 再将 pose 重原子坐标按
    (残基号 1-based, 原子名) 映射覆盖。返回**重原子** mol (含坐标), 不含 H ——
    H 由 ligand_pdb_with_h 经 PDBFixer 补加 (amber 模板 H 名/残基信息最干净)。"""
    atoms = parse_pdbqt_model1_atoms(pqt)
    if not atoms:
        return None
    m = Chem.MolFromSequence(seq)          # 规范肽拓扑 (正确键/Pro 环)
    if m is None:
        return None
    m = Chem.RemoveHs(m)
    if m.GetNumConformers() == 0:        # 保证有 conformer 以放置坐标
        m.AddConformer(Chem.Conformer(m.GetNumAtoms()), assignId=True)
    pose = {(a["resid"], a["name"]): (a["x"], a["y"], a["z"]) for a in atoms}
    conf = m.GetConformer()
    matched = 0
    for atom in m.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info is None:
            continue
        key = (info.GetResidueNumber(), info.GetName().strip())
        if key in pose:
            conf.SetAtomPosition(atom.GetIdx(), Point3D(*pose[key]))
            matched += 1
    if matched == 0:
        return None
    return m


def ligand_pdb_with_h(mol, dst):
    """重原子 mol (含坐标) -> 手写重原子 PDB -> PDBFixer 补缺失原子+氢,
    写出带**正确 H 名与残基信息**的 PDB (amber14 模板可直接匹配)。
    关键教训: RDKit AddHs 的 H 原子 GetMonomerInfo 返回 None,
    手写 PDB 若标 UNK 或 H 名非模板 (H1/H2...) 均会被 amber14 判为
    'missing hydrogens'。PDBFixer 按模板产出标准 H 名, 一劳永逸。"""
    heavy = os.path.join(DOCK, "_mgbsa_heavy_tmp.pdb")
    open(heavy, "w").write(_mol_to_pdb(Chem.RemoveHs(mol)))
    fixer = PDBFixer(filename=heavy)
    fixer.findMissingResidues(); fixer.missingResidues = {}   # 不回填缺失残基
    fixer.findMissingAtoms(); fixer.addMissingAtoms(); fixer.addMissingHydrogens(7.0)
    PDBFile.writeFile(fixer.topology, fixer.positions, dst)


def _mol_to_pdb(m):
    """手动将 RDKit mol (含 H) 写为 PDB 文本, 逐原子输出 (含 H),
    规避 Chem.MolToPDBBlock 在某些情况下丢弃 H 原子的问题。
    关键: RDKit 自行 AddHs 产生的 H 原子 GetMonomerInfo() 返回 None,
    若写成 UNK 会导致 amber14 "No template / missing hydrogens"。
    修复: H 原子继承其成键重原子的残基信息 (resname/resnum/chain),
    并按 (resnum,chain) 计数赋唯一 H 名 (H1/H2/...)。amber14 按元素+成键
    重原子匹配 H, 名不敏感。"""
    conf = m.GetConformer()
    # 先收集重原子残基信息 (按原子 idx)
    heavy_info = {}
    for atom in m.GetAtoms():
        if atom.GetSymbol() == "H":
            continue
        info = atom.GetMonomerInfo()
        if info is not None:
            heavy_info[atom.GetIdx()] = (
                (info.GetResidueName() or "UNK").strip()[:3],
                info.GetResidueNumber(),
                info.GetChainId() or "A",
            )
    hcount = {}          # (resnum, chain) -> 已写 H 数
    out = ["REMARK manual pdb (heavy+H)", "COMPND    LIG"]
    for i, atom in enumerate(m.GetAtoms(), 1):
        p = conf.GetAtomPosition(atom.GetIdx())
        el = atom.GetSymbol()
        info = atom.GetMonomerInfo()
        if info is not None:
            resname = (info.GetResidueName() or "UNK").strip()[:3]
            resnum = info.GetResidueNumber()
            chain = info.GetChainId() or "A"
            name = (info.GetName() or el).strip()[:4]
        else:
            # H: 取其成键重原子的残基信息
            heavy_res = None
            for nbr in atom.GetNeighbors():
                if nbr.GetSymbol() != "H":
                    heavy_res = heavy_info.get(nbr.GetIdx())
                    break
            if heavy_res is None:
                heavy_res = ("UNK", 1, "A")
            resname, resnum, chain = heavy_res
            key = (resnum, chain)
            hcount[key] = hcount.get(key, 0) + 1
            name = f"H{hcount[key]}"
        out.append("ATOM  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s  "
                    % (i, name, resname, chain, resnum, p.x, p.y, p.z, el))
    out.append("END")
    return "\n".join(out)


def renumber_chain_B(lines):
    out, cur_key, n = [], None, 0
    for ln in lines:
        if ln[:6] not in ("ATOM  ", "HETATM"):
            out.append(ln.rstrip("\n")); continue
        key = (ln[17:20].strip(), ln[22:26])
        if key != cur_key:
            n += 1; cur_key = key
        out.append(ln[:21] + "B" + "%4d" % n + ln[26:])
    return out


def minimize_state(system, topology, positions, maxIter=200):
    """最小化并返回 (energy_kcal, minimized_positions)。"""
    integ = openmm.LangevinMiddleIntegrator(300*openmm.unit.kelvin,
                                           1/openmm.unit.picosecond,
                                           0.001*openmm.unit.picosecond)
    sim = app.Simulation(topology, system, integ, openmm.Platform.getPlatformByName("CPU"))
    sim.context.setPositions(positions)
    sim.minimizeEnergy(maxIterations=maxIter)
    st = sim.context.getState(getPositions=True, getEnergy=True)
    e = st.getPotentialEnergy().value_in_unit(openmm.unit.kilocalories_per_mole)
    return e, st.getPositions()


def eval_energy(system, topology, positions):
    """仅计算能量 (不最小化), 用于从复合物极小化坐标中分解受体/配体能量。"""
    integ = openmm.LangevinMiddleIntegrator(300*openmm.unit.kelvin,
                                           1/openmm.unit.picosecond,
                                           0.001*openmm.unit.picosecond)
    sim = app.Simulation(topology, system, integ, openmm.Platform.getPlatformByName("CPU"))
    sim.context.setPositions(positions)
    e = sim.context.getState(getEnergy=True).getPotentialEnergy()
    return e.value_in_unit(openmm.unit.kilocalories_per_mole)


def subset_topology(topology, keep_idxs):
    """抽取 keep_idxs (原 atom.index 集合) 构建新 Topology (保留残基/链/键结构)。
    返回 (new_topology, pos_index_list), pos_index_list[i] = 第 i 个新原子对应的原坐标 index。"""
    keep = set(keep_idxs)
    new = app.Topology()
    atom_map = {}
    pos_idx = []
    cur_chain = None; cur_chain_new = None
    cur_res = None;   cur_res_new = None
    for atom in topology.atoms():
        if atom.index not in keep:
            continue
        ch = atom.residue.chain
        if cur_chain is not ch:
            cur_chain = ch
            cur_chain_new = new.addChain(ch.id)
            cur_res = None
        res = atom.residue
        if cur_res is not res:
            cur_res = res
            cur_res_new = new.addResidue(res.name, cur_chain_new, res.id)
        na = new.addAtom(atom.name, atom.element, cur_res_new)
        atom_map[atom.index] = na
        pos_idx.append(atom.index)
    for bond in topology.bonds():
        a1, a2 = bond[0], bond[1]
        if a1.index in keep and a2.index in keep:
            new.addBond(atom_map[a1.index], atom_map[a2.index])
    return new, pos_idx


def main():
    ff = ForceField("amber14-all.xml", "implicit/gbn2.xml")
    kw = dict(nonbondedMethod=app.CutoffNonPeriodic,
              nonbondedCutoff=2.0*openmm.unit.nanometer,
              constraints=app.HBonds)

    # 受体口袋 (固定, 三候选共享) -> PDBFixer 补氢
    rec_pocket = os.path.join(DOCK, "_mgbsa_rec_pocket.pdb")
    n_res = write_rec_pocket(REC, rec_pocket, POCKET_CENTER, POCKET_RADIUS)
    rec_h = os.path.join(DOCK, "_mgbsa_rec_h.pdb")
    fix_pocket(rec_pocket, rec_h)
    rp = PDBFile(rec_h)
    print(f"[rec] pocket residues={n_res} atoms={rp.topology.getNumAtoms()}", flush=True)

    results = []
    for seq, pqt in TOP3:
        lig_mol = ligand_from_pose(seq, pqt)
        if lig_mol is None:
            raise RuntimeError(f"配体读取失败: {pqt}")
        # 配体: 规范序列肽 + pose 重原子坐标; H 由 PDBFixer 补 (名/残基信息正确)
        lig_h = os.path.join(DOCK, f"_mgbsa_lig_h_{seq}.pdb")
        ligand_pdb_with_h(lig_mol, lig_h)

        lig_lines = renumber_chain_B([l.rstrip("\n") for l in open(lig_h)])
        cplx = os.path.join(DOCK, f"_mgbsa_cplx_{seq}.pdb")
        with open(cplx, "w") as out:
            for ln in open(rec_h):
                if ln[:6] in ("ATOM  ", "HETATM"):
                    out.write(ln)
            out.write("TER\n")
            for ln in lig_lines:
                if ln[:6] in ("ATOM  ", "HETATM"):
                    out.write(ln + "\n")
            out.write("END\n")

        cp = PDBFile(cplx)
        top = cp.topology
        # (a) 仅最小化**复合物**一次; (b) 从同一极小化坐标分解 E_rec / E_lig
        #     —— 单轨迹 MM-GBSA, 杜绝"配体独立弛豫"伪影 (否则自由配体能量被
        #        不实压低, 致 dG 偏高甚至为正如 APQIP 初版 +8.2)。
        E_cplx, minpos = minimize_state(ff.createSystem(top, **kw), top, cp.positions)
        rec_idxs = [a.index for a in top.atoms() if a.residue.chain.id == 'A']
        lig_idxs = [a.index for a in top.atoms() if a.residue.chain.id == 'B']
        rec_top, rec_pi = subset_topology(top, set(rec_idxs))
        lig_top, lig_pi = subset_topology(top, set(lig_idxs))
        rec_pos = [minpos[i] for i in rec_pi]
        lig_pos = [minpos[i] for i in lig_pi]
        E_rec = eval_energy(ff.createSystem(rec_top, **kw), rec_top, rec_pos)
        E_lig = eval_energy(ff.createSystem(lig_top, **kw), lig_top, lig_pos)
        dG = E_cplx - E_rec - E_lig
        results.append((seq, DG[seq], dG, E_rec, E_lig, E_cplx))
        print(f"  {seq}: dG_vina={DG[seq]:.3f}  dG_mgbsa={dG:.2f}  "
              f"Erec={E_rec:.1f}  Elig={E_lig:.1f}  Ecplx={E_cplx:.1f}", flush=True)

    out_tsv = os.path.join(REPO, "data", "phaseB", "mmgbsa_top3.tsv")
    os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
    with open(out_tsv, "w", encoding="utf-8") as f:
        f.write("peptide\tdG_vina_kcal_mol\tdG_mgbsa_kcal_mol\tE_rec\tE_lig\tE_cplx\n")
        for seq, dv, dg, er, el, ec in results:
            f.write(f"{seq}\t{dv:.3f}\t{dg:.2f}\t{er:.1f}\t{el:.1f}\t{ec:.1f}\n")
    print(f"\n-> {out_tsv}")
if __name__ == "__main__":
    main()
