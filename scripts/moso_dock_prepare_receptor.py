# -*- coding: utf-8 -*-
"""
[Receptor prep] 1WCY.pdb -> 1WCY_receptor.pdbqt
Steps: strip water / native ligand (sitagliptin A1201) -> ADT prepare_receptor4.py to pdbqt
Deps: AutoDockTools (pythonsh + MGLTools) or openbabel (obabel -xr)
This machine has no ADT; here we emit a clean receptor pdb for the user to convert locally.
"""
import re, os
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # scripts/ -> repo root
DOCK = os.path.join(REPO, "docking")
# 1WCY.pdb must be downloaded separately from RCSB PDB and placed in docking/ (not under version control)
SRC=os.path.join(DOCK, "1WCY.pdb")
OUT=os.path.join(DOCK, "1WCY_clean.pdb")
keep=[]
for line in open(SRC):
    rec=line[:6]
    if rec in ("TER   ","END   ","ENDMDL"):
        keep.append(line); continue
    if rec=="HETATM":
        # remove ligand A1201 (sitagliptin) and crystal water (HOH); ions/glycerol etc. may be kept at discretion
        name=line[17:20].strip()
        if name=="A1201" or name=="HOH":
            continue
    if rec in ("ATOM  ","HETATM"):
        keep.append(line); continue
    keep.append(line)
open(OUT,"w").writelines(keep)
print(f"Wrote clean receptor (ligand/water stripped): {OUT}  ({len(keep)} lines)")

print("\n--- run on user's machine (after installing ADT/MGLTools) ---")
print("  pythonsh $MGLTOOLS/prepare_receptor4.py -r 1WCY_clean.pdb -o 1WCY_receptor.pdbqt -A checkhydrogens")
print("  or: obabel 1WCY_clean.pdb -O 1WCY_receptor.pdbqt -xr")
print(f"\nPocket (grid) center (measured from sitagliptin A1201): center_x=62.8 center_y=47.7 center_z=4.8  size=30 30 30")
print("The template paper gives (54,62,37) as an approximate same-pocket value; the measured values are more reliable.")
