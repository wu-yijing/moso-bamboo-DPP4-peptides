# -*- coding: utf-8 -*-
"""
Benchmark validation: known DPP4 inhibitors (peptides) -> reuse this pipeline to dock 1WCY.
Validates two points:
  (a) whether known inhibitors dock into the DPP4 active pocket (overlapping the literature binding site);
  (b) whether Vina dG ranking agrees with literature relative inhibitory activity (stronger inhibitors -> more negative dG).
Fully reuses moso_dock_run_idppiv.py's prep_pdbqt + dock logic to guarantee method consistency.
Vina binary: specified via env var VINA_EXE (consistent with current pipeline convention); if unset, uses vina on PATH.
"""
import os, subprocess, math
from collections import defaultdict
from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import pybel

REPO   = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOCK   = os.path.join(REPO, "docking")
REC    = os.path.join(DOCK, "1WCY_receptor.pdbqt")
BOX    = os.path.join(DOCK, "moso_box.txt")
VINA   = os.environ.get("VINA_EXE", "vina")
LIGDIR = os.path.join(DOCK, "benchmark_ligands")
os.makedirs(LIGDIR, exist_ok=True)

# (seq, name, ic50_note, role)
# role: positive=classic human DPP4 inhibitor; food=food-source control; negative=non-inhibitory negative control
KNOWN = [
    ("IPI", "Diprotin A",          "classic competitive human DPP4 inhibitor (Streptomyces), IC50 ~ 20-50 uM", "positive"),
    ("VPL", "Diprotin B",          "classic human DPP4 inhibitor, IC50 ~ 30-60 uM",                       "positive"),
    ("IPP", "Ile-Pro-Pro (milk-derived)", "milk-protein-source DPP4 inhibitory peptide, IC50 ~ tens-to-hundreds uM", "food"),
    ("VPP", "Val-Pro-Pro (milk-derived)", "milk-protein-source DPP4 inhibitory peptide, IC50 ~ tens-to-hundreds uM", "food"),
    ("AA",  "Ala-Ala (negative control)", "non DPP4 inhibitor, expected weak/no activity",                   "negative"),
]

def prep_pdbqt(seq, idx):
    m = Chem.MolFromSequence(seq)
    if m is None:
        return None
    m = Chem.AddHs(m)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    if not AllChem.EmbedMultipleConfs(m, numConfs=1, params=params):
        return None
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
        return None, None
    best = None
    for line in r.stdout.splitlines():
        if line.strip().startswith("1 "):
            try:
                best = float(line.split()[1])
            except (ValueError, IndexError):
                pass
            break
    return best, (out if os.path.exists(out) else None)

def parse_atoms(path, model=0):
    """Parse pdbqt: model=0 collects all (receptor has no MODEL header); model=N collects only the N-th MODEL (ligand pose)."""
    atoms = []; cur = 0
    with open(path) as fh:
        for ln in fh:
            s = ln[:6]
            if s == "MODEL ":
                try:
                    cur = int(ln.split()[1])
                except Exception:
                    cur = 0
                continue
            if s == "ENDMDL":
                continue
            if s in ("ATOM  ", "HETATM"):
                if model == 0 or cur == model:
                    atoms.append(dict(
                        x=float(ln[30:38]), y=float(ln[38:46]), z=float(ln[46:54]),
                        name=ln[12:16].strip(), resName=ln[17:20].strip(),
                        chain=ln[21:22], resSeq=int(ln[22:26])))
    return atoms

def analyze(out_pqt, receptor_pdbqt):
    lig = parse_atoms(out_pqt, model=1)
    rec = parse_atoms(receptor_pdbqt, model=0)
    if not lig:
        return None, []
    cx = sum(a["x"] for a in lig)/len(lig)
    cy = sum(a["y"] for a in lig)/len(lig)
    cz = sum(a["z"] for a in lig)/len(lig)
    res_atoms = defaultdict(list)
    for a in rec:
        rn = a["resName"]
        if not rn or len(rn) < 2:
            continue
        res_atoms[(a["chain"], rn, a["resSeq"])].append(a)
    def dist(a, b):
        return math.sqrt((a["x"]-b["x"])**2 + (a["y"]-b["y"])**2 + (a["z"]-b["z"])**2)
    res_min = {}
    for res, ats in res_atoms.items():
        res_min[res] = min(dist(la, ra) for la in lig for ra in ats)
    near = sorted(res_min.items(), key=lambda kv: kv[1])[:12]
    # safe label: avoid [1] indexing on short/empty resName
    near_disp = []
    for res, d in near:
        chain, rn, resSeq = (list(res) + [None, None, None])[:3]
        rn = rn or "?"
        near_disp.append((f"{rn[:3]:>3}{resSeq}", d))
    return (cx, cy, cz), near_disp

def main():
    rows = []
    for i, (seq, name, note, role) in enumerate(KNOWN):
        pqt = prep_pdbqt(seq, i)
        dg, out = (None, None)
        if pqt:
            dg, out = dock(pqt, i)
        center, near = (None, None)
        if out:
            center, near = analyze(out, REC)
        rows.append((seq, name, role, note, dg, center, near))
        near_str = ", ".join(f"{lab}:{d:.1f}" for lab, d in (near or [])[:6])
        print(f"  {seq:5s} {name:22s} dG={dg} center={tuple(round(v,1) for v in center) if center else None} near=[{near_str}]", flush=True)
    out_tsv = os.path.join(DOCK, "benchmark_dpp4_inhibitors.tsv")
    with open(out_tsv, "w", encoding="utf-8") as f:
        f.write("seq\tname\trole\tic50_note\tdG_kcal_mol\tpose_center\tnearest_residues\n")
        for seq, name, role, note, dg, center, near in rows:
            cstr = ",".join(f"{v:.1f}" for v in center) if center else ""
            nstr = "; ".join(f"{lab}({d:.1f}A)" for lab, d in (near or [])[:10])
            f.write(f"{seq}\t{name}\t{role}\t{note}\t{dg if dg is not None else 'NA'}\t{cstr}\t{nstr}\n")
    print(f"\n-> benchmark result: {out_tsv}")

if __name__ == "__main__":
    main()
