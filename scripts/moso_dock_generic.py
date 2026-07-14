# -*- coding: utf-8 -*-
"""
通用 Vina 对接(同一 RDKit 制备管线) —— 可用于任意队列, 断点续跑 + 即时落盘.
用法: python moso_dock_generic.py <QUEUE> <RESULTS_TSV> <LIGDIR>
"""
import os, sys, subprocess
from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import pybel

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # scripts/ -> 仓库根
DATA = os.path.join(REPO, "data")
DOCK = os.path.join(REPO, "docking")

def _resolve_in(name, *dirs):
    """按给定目录顺序解析已存在的输入文件; 绝对路径/已存在则原样返回;
    都找不到时默认落在第一个目录(用于输出)。"""
    if os.path.isabs(name) or os.path.exists(name):
        return name
    for d in dirs:
        cand = os.path.join(d, name)
        if os.path.exists(cand):
            return cand
    return os.path.join(dirs[0], name)

QUEUE  = _resolve_in(sys.argv[1], DATA, DOCK)                         # 队列多在 data/
RESULTS= sys.argv[2] if os.path.isabs(sys.argv[2]) else os.path.join(DOCK, sys.argv[2])
LIGDIR = sys.argv[3] if os.path.isabs(sys.argv[3]) else os.path.join(DOCK, sys.argv[3])
REC    = os.path.join(DOCK, "1WCY_receptor.pdbqt")
BOX    = os.path.join(DOCK, "moso_box.txt")
# AutoDock Vina 第三方二进制(不纳入版本库): 默认走 PATH; 可用 VINA_EXE 覆盖。
VINA   = os.environ.get("VINA_EXE", "vina")
os.makedirs(LIGDIR, exist_ok=True)

def _autolocate_babel_datadir():
    if os.environ.get("BABEL_DATADIR"):
        return
    try:
        import openbabel as _ob
        share = os.path.join(os.path.dirname(_ob.__file__), "share", "openbabel")
        if os.path.isdir(share):
            vers = [d for d in os.listdir(share) if os.path.isdir(os.path.join(share, d))]
            if vers:
                os.environ["BABEL_DATADIR"] = os.path.join(share, sorted(vers)[-1])
    except Exception:
        pass

_autolocate_babel_datadir()

def prep_pdbqt(seq, idx):
    m = Chem.MolFromSequence(seq)
    if m is None: return None
    m = Chem.AddHs(m)
    p = AllChem.ETKDGv3(); p.randomSeed = 42
    if not AllChem.EmbedMultipleConfs(m, numConfs=1, params=p): return None
    AllChem.MMFFOptimizeMoleculeConfs(m)
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
            try: return float(line.split()[1])
            except (ValueError, IndexError): pass
    return None

done = set()
if os.path.exists(RESULTS):
    for ln in open(RESULTS, encoding="utf-8"):
        ln = ln.strip()
        if ln and not ln.startswith("peptide"):
            done.add(ln.split("\t")[0])

rows = [l.rstrip("\n").split("\t") for l in open(QUEUE, encoding="utf-8") if l.strip()]
print(f"队列 {sys.argv[1]}: {len(rows)} 条, 已完成 {len(done)} 条", flush=True)
if not os.path.exists(RESULTS):
    with open(RESULTS, "w", encoding="utf-8") as f:
        f.write("peptide\tiDPPIV_score\tdG_kcal_mol\n")

count = 0
for i, parts in enumerate(rows):
    seq = parts[0]; sc = parts[1] if len(parts) > 1 else "0"
    if seq in done: continue
    try:
        pqt = prep_pdbqt(seq, i); aff = None
        if pqt is not None: aff = dock(pqt, i)
        with open(RESULTS, "a", encoding="utf-8") as f:
            f.write(f"{seq}\t{float(sc):.3f}\t{aff if aff is not None else 'NA'}\n")
        done.add(seq); count += 1
        print(f"  {seq:7s} dG={aff}", flush=True)
    except Exception as e:
        with open(RESULTS, "a", encoding="utf-8") as f:
            f.write(f"{seq}\t{float(sc):.3f}\tERR\n")
        done.add(seq)
        print(f"  {seq:7s} 异常: {e}", flush=True)

print(f"\n本轮新增 {count} 条; 累计 {len(done)}/{len(rows)} -> {RESULTS}", flush=True)
