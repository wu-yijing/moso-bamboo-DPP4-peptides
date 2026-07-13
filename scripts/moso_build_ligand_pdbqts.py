# -*- coding: utf-8 -*-
"""
毛竹候选肽 → 批量生成 3D PDBQT 文件（用 openbabel）
输出: moso_ligands/ 目录下所有 <序号>_<肽名>.pdbqt
同时生成 batch_dock.cmd (Windows 一键 Vina 对接)
"""
import os, subprocess, glob
from openbabel import pybel

# 氨基酸三字母 SMILES 片段（简化版：N→C 端直接拼接）
AA = {'A':'C[C@@H](N)C(=O)O','C':'SC[C@@H](N)C(=O)O','D':'OC(=O)C[C@@H](N)C(=O)O',
      'E':'OC(=O)CC[C@@H](N)C(=O)O','F':'C1=CC=CC=C1C[C@@H](N)C(=O)O',
      'G':'NCC(=O)O','H':'C1=C[N+]=CN1C[C@@H](N)C(=O)O',
      'I':'CC[C@H](C)[C@@H](N)C(=O)O','K':'NCCCC[C@@H](N)C(=O)O',
      'L':'CC(C)C[C@@H](N)C(=O)O','M':'CSCC[C@@H](N)C(=O)O',
      'N':'NC(=O)C[C@@H](N)C(=O)O','P':'C1C[C@H](N)C1C(=O)O',
      'Q':'NC(=O)CC[C@@H](N)C(=O)O','R':'N=C(N)NCCC[C@@H](N)C(=O)O',
      'S':'OC[C@@H](N)C(=O)O','T':'C[C@H](O)[C@@H](N)C(=O)O',
      'V':'CC(C)[C@@H](N)C(=O)O','W':'C1=CC=C2C(=C1)C(=CN2)C[C@@H](N)C(=O)O',
      'Y':'C1=CC(=CC=C1C[C@@H](N)C(=O)O)O'}

def seq_to_smiles(seq):
    """Build proper peptide SMILES with peptide bonds"""
    if len(seq)==0 or any(a not in AA for a in seq): return None
    parts = [AA[a] for a in seq]
    # N-term: keep NH2; C-term: keep COOH
    # Between: replace COOH of i with CONH of i+1
    # S M I L E S: concatenate with proper NC(=O) linkage
    result = parts[0]  # N-term
    for p in parts[1:]:
        # Remove OH from previous C-term, add NC(=O) then attach next amino (remove NH2)
        result = result[:-1] + 'NC(=O)' + p[1:]
    return result

# 读队列
queue = [l.strip().split('\t') for l in open("moso_dock_queue.txt") if l.strip()]
print(f"对接队列: {len(queue[:60])} 条")

outdir = "moso_ligands"
os.makedirs(outdir, exist_ok=True)

batch = ["@echo off", "REM 批量 Vina 对接 - 用 Docker 或本地 Vina"]
batch.append(f"set REC=1WCY_receptor.pdbqt")
batch.append(f"set BOX=moso_box.txt")

for i, row in enumerate(queue[:60]):
    seq = row[0]
    pr_score = row[1]
    smi = seq_to_smiles(seq)
    if smi is None:
        print(f"  {i:02d} {seq:6s} SKIP (unknown AA)")
        continue
    try:
        mol = pybel.readstring('smi', smi)
        mol.make3D()
        fname = f"{i:02d}_{seq}.pdbqt"
        fpath = os.path.join(outdir, fname)
        mol.write('pdbqt', fpath, overwrite=True)
        size = os.path.getsize(fpath)
        batch.append(
            f"echo --- Docking {seq} ---"
        )
        batch.append(
            f"vina --receptor %REC% --ligand {outdir}\\{fname} "
            f"--config %BOX% --out {outdir}\\dock_{i:02d}_{seq}.pdbqt --cpu 4"
        )
        print(f"  ✅ {i:02d} {seq:6s} -> {fname} ({size} B)")
    except Exception as e:
        print(f"  ❌ {i:02d} {seq:6s} error: {e}")

# Write batch script
with open("batch_dock.cmd", "w") as f:
    f.write("\n".join(batch))
    f.write("\necho ALL DOCKING COMPLETE\npause\n")

# Write Linux bash version too
batch_sh = ["#!/bin/bash",
    "REC=1WCY_receptor.pdbqt",
    "BOX=moso_box.txt"]
for i, row in enumerate(queue[:60]):
    seq = row[0]
    batch_sh.append(
        f"vina --receptor $REC --ligand {outdir}/{i:02d}_{seq}.pdbqt "
        f"--config $BOX --out {outdir}/dock_{i:02d}_{seq}.pdbqt --cpu 4"
    )
with open("batch_dock.sh", "w") as f:
    f.write("\n".join(batch_sh) + "\n")

n_ok = len([f for f in os.listdir(outdir) if f.endswith('.pdbqt')])
print(f"\n成功生成 {n_ok}/{len(queue[:60])} 个配体 PDBQT")
print(f"Windows:   batch_dock.cmd")
print(f"Linux/WSL: bash batch_dock.sh")
print(f"出文件目录: {outdir}/")
