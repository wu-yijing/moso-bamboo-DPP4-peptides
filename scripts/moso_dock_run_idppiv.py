# -*- coding: utf-8 -*-
"""
[配体准备 + 批量对接] 毛竹 iDPPIV-SCM 优先化 Top-60 队列 -> Vina 对接 DPP4 (1WCY)
--------------------------------------------------------------------------
制备方法(关键修复):
  旧 openbabel make3D 在本环境因 MMFF94/ring-fragments 数据文件缺失(且 wheel 内
  文件名拼写错误 rigid-*) 而损坏, 生成的 PDBQT 扭转树畸形导致 Vina 静默崩溃。
  现改用: RDKit MolFromSequence(自带 MMFF94, 无需外部数据) 生成合规 3D
        -> 写 SDF -> openbabel 仅做 SDF->PDBQT 格式转换(不调 make3D)。
  产出 PDBQT 电荷为 0.00、原子类型为 C/OA/N/NA/HD 等退化型,
  与既有 moso_dock_ranking.txt(旧代理队列) 的制备惯例一致, 可比。

鲁棒性改进(本次重跑):
  - 断点续跑: 已完成肽写入结果与时跳过, 中途被杀可重跑补齐。
  - 即时落盘: 每肽对接完立即 append 到结果 TSV, 不依赖最终一次性写盘。
  - vina 超时保护: 单肽子进程超时 300s, 超时记 NA 不阻塞整批。
"""
import os, subprocess, sys
from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import pybel

HERE = "E:/workbuddy/Claw"
QUEUE   = os.path.join(HERE, "moso_dock_queue_idppiv.txt")
REC     = os.path.join(HERE, "1WCY_receptor.pdbqt")
BOX     = os.path.join(HERE, "moso_box.txt")
VINA    = os.path.join(HERE, "vina.exe")
LIGDIR  = os.path.join(HERE, "moso_ligands_idppiv")
RESULTS = os.path.join(HERE, "moso_dock_results_idppiv.tsv")
os.makedirs(LIGDIR, exist_ok=True)
# 抑制 openbabel 缺失 MMFF94 数据时的报错噪声(转换步骤已不依赖之)
os.environ["BABEL_DATADIR"] = r"C:/Users/Administrator/.workbuddy/binaries/python/envs/default/Lib/site-packages/openbabel/share/openbabel/3.2.1"

def prep_pdbqt(seq, idx):
    """RDKit 生成 3D(SDF) -> openbabel 转 PDBQT(不调 make3D)"""
    m = Chem.MolFromSequence(seq)
    if m is None:
        return None
    m = Chem.AddHs(m)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    if not AllChem.EmbedMultipleConfs(m, numConfs=1, params=params):
        return None
    AllChem.MMFFOptimizeMoleculeConfs(m)   # RDKit 自带 MMFF94
    sdf = os.path.join(LIGDIR, f"{idx:02d}_{seq}.sdf")
    w = Chem.SDWriter(sdf); w.write(m); w.close()
    mb = pybel.readfile("sdf", sdf).__next__()
    pqt = os.path.join(LIGDIR, f"{idx:02d}_{seq}.pdbqt")
    mb.write("pdbqt", pqt, overwrite=True)
    return pqt if os.path.getsize(pqt) > 50 else None

def dock(pqt, idx):
    out = os.path.join(LIGDIR, f"dock_{idx:02d}.pdbqt")
    try:
        r = subprocess.run([VINA, "--receptor", REC, "--ligand", pqt,
                            "--config", BOX, "--out", out, "--cpu", "4"],
                           capture_output=True, text=True, timeout=300)
    except subprocess.TimeoutExpired:
        return None
    for line in r.stdout.splitlines():
        if line.strip().startswith("1 "):
            try:
                return float(line.split()[1])
            except (ValueError, IndexError):
                pass
    return None

# ---- 断点续跑: 加载已完成 ----
done = set()
if os.path.exists(RESULTS):
    with open(RESULTS, encoding="utf-8") as f:
        for ln in f:
            ln = ln.strip()
            if ln and not ln.startswith("peptide"):
                done.add(ln.split("\t")[0])

rows = [l.rstrip("\n").split("\t") for l in open(QUEUE, encoding="utf-8") if l.strip()]
print(f"iDPPIV 队列: {len(rows)} 条, 已完成 {len(done)} 条", flush=True)

# 若是首次, 写表头
if not os.path.exists(RESULTS):
    with open(RESULTS, "w", encoding="utf-8") as f:
        f.write("peptide\tiDPPIV_score\tdG_kcal_mol\n")

count = 0
for i, (seq, sc, _tier) in enumerate(rows):
    if seq in done:
        continue
    try:
        pqt = prep_pdbqt(seq, i)
        aff = None
        if pqt is not None:
            aff = dock(pqt, i)
        # 即时落盘
        with open(RESULTS, "a", encoding="utf-8") as f:
            f.write(f"{seq}\t{float(sc):.3f}\t{aff if aff is not None else 'NA'}\n")
        done.add(seq)
        count += 1
        print(f"  {seq:7s} iDPPIV={float(sc):>6.2f}  dG={aff}", flush=True)
    except Exception as e:
        with open(RESULTS, "a", encoding="utf-8") as f:
            f.write(f"{seq}\t{float(sc):.3f}\tERR\n")
        done.add(seq)
        print(f"  {seq:7s} 异常: {e}", flush=True)

print(f"\n本轮新增 {count} 条; 累计 {len(done)}/{len(rows)} 条 -> {RESULTS}", flush=True)
