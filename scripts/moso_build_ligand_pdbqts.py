# -*- coding: utf-8 -*-
"""
Moso-bamboo candidate peptides -> batch-generate 3D PDBQT files (using openbabel)
Output: all <index>_<peptide>.pdbqt files under moso_ligands/
Also generate batch_dock.cmd (one-click Windows Vina docking)
"""
import os, subprocess, glob
from openbabel import pybel

# amino-acid three-letter SMILES fragments (simplified: direct N->C concatenation)
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
    if len(seq) == 0 or any(a not in AA for a in seq): return None
    parts = [AA[a] for a in seq]
    # N-term: keep NH2; C-term: keep COOH
    # Between: replace COOH of i with CONH of i+1
    # SMILES: concatenate with proper NC(=O) linkage
    result = parts[0]  # N-term
    for p in parts[1:]:
        # Remove OH from previous C-term, add NC(=O) then attach next amino (remove NH2)
        result = result[:-1] + 'NC(=O)' + p[1:]
    return result

# read queue
queue = [l.strip().split('\t') for l in open("moso_dock_queue.txt") if l.strip()]
print(f"Docking queue: {len(queue[:60])} entries")

outdir = "moso_ligands"
os.makedirs(outdir, exist_ok=True)

batch = ["@echo off", "REM batch Vina docking - using Docker or local Vina"]
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
        print(f"  [OK] {i:02d} {seq:6s} -> {fname} ({size} B)")
    except Exception as e:
        print(f"  [ERR] {i:02d} {seq:6s} error: {e}")

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
print(f"\nSuccessfully generated {n_ok}/{len(queue[:60])} ligand PDBQT files")
print(f"Windows:   batch_dock.cmd")
print(f"Linux/WSL: bash batch_dock.sh")
print(f"Output directory: {outdir}/")
