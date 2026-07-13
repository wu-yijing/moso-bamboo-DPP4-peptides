# -*- coding: utf-8 -*-
"""
[配体准备 + 批量对接] 毛竹候选肽 -> Vina 对接 DPP4 (1WCY)
依赖(用户本机): rdkit, openbabel, autodock-vina
读取 moso_dock_queue.txt -> 为每个肽:
  1) RDKit 生成 3D 构象 -> SDF
  2) obabel SDF -> PDBQT
  3) vina --receptor 1WCY_receptor.pdbqt --ligand X.pdbqt --config box.txt
  4) 解析最佳 binding affinity (kcal/mol)
输出: moso_dock_results.tsv (肽, PR评分, 亲和力, 等级)
"""
import os, subprocess, sys

QUEUE = "E:/workbuddy/Claw/moso_dock_queue.txt"
REC   = "E:/workbuddy/Claw/1WCY_receptor.pdbqt"
BOX   = "E:/workbuddy/Claw/moso_box.txt"
OUTDIR= "E:/workbuddy/Claw/dock_out"
N_CONF = 20          # 每个肽生成构象数
TOP_N  = 50          # 实际跑前 N 条(默认60队列全部)

# ---- 1) 口袋盒子 ----
open(BOX,"w").write(
    "center_x = 62.8\ncenter_y = 47.7\ncenter_z = 4.8\n"
    "size_x = 30\nsize_y = 30\nsize_z = 30\n"
)
os.makedirs(OUTDIR, exist_ok=True)

from rdkit import Chem
from rdkit.Chem import AllChem

def pep_to_pdbqt(seq, idx):
    m = Chem.MolFromSequence(seq)
    if m is None: return None
    m = Chem.AddHs(m)
    AllChem.EmbedMultipleConfs(m, numConfs=N_CONF, randomSeed=42)
    sdf = f"{OUTDIR}/lig_{idx}.sdf"
    w = Chem.SDWriter(sdf); w.write(m); w.close()
    pqt = f"{OUTDIR}/lig_{idx}.pdbqt"
    subprocess.run(["obabel", sdf, "-O", pqt, "-xr", "-p", "7.4"], check=True)
    return pqt

def dock(pqt, idx):
    out = f"{OUTDIR}/dock_{idx}.pdbqt"
    r = subprocess.run(["vina","--receptor",REC,"--ligand",pqt,
                        "--config",BOX,"--out",out,"--cpu","4"],
                       capture_output=True, text=True)
    best=None
    for line in r.stdout.splitlines():
        if line.strip().startswith("1 "):
            best=float(line.split()[1]); break
    return best

rows=[l.rstrip("\n").split("\t") for l in open(QUEUE) if l.strip()]
results=[]
for i,(seq,score,_tier) in enumerate(rows[:TOP_N]):
    try:
        pqt = pep_to_pdbqt(seq, i)
        aff = dock(pqt, i) if pqt else None
        results.append((seq, float(score), aff))
        print(f"  {seq:6s} PR={score}  dG={aff}")
    except Exception as e:
        print(f"  {seq:6s} 失败: {e}")

with open("E:/workbuddy/Claw/moso_dock_results.tsv","w") as f:
    f.write("peptide\tPR_score\tdG_kcal_mol\n")
    for s,sc,a in sorted(results, key=lambda x:(x[2] if x[2] is not None else 9e9)):
        f.write(f"{s}\t{sc:.3f}\t{a if a is not None else 'NA'}\n")
print(f"\n对接完成 -> E:/workbuddy/Claw/moso_dock_results.tsv ({len(results)} 条)")
