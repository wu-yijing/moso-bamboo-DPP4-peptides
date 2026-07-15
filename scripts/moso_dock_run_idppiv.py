# -*- coding: utf-8 -*-
"""
[Ligand prep + batch docking] moso-bamboo iDPPIV-SCM-prioritized Top-60 queue -> Vina docking vs DPP4 (1WCY)
--------------------------------------------------------------------------
Preparation fix (key):
  the old openbabel make3D was broken in this environment (missing MMFF94/ring-fragments
  data files, plus a filename typo rigid-* inside the wheel), producing malformed PDBQTs
  whose torsion trees crashed Vina silently.
  Now we use: RDKit MolFromSequence (ships its own MMFF94, no external data) to make a
        valid 3D -> write SDF -> openbabel only does SDF->PDBQT format conversion (no make3D).
  Produced PDBQTs have charge 0.00 and degenerate atom types C/OA/N/NA/HD,
  consistent with the existing moso_dock_ranking.txt (old proxy queue) preparation
  convention, hence comparable.

Robustness improvements (this re-run):
  - resume from checkpoint: finished peptides are written to results and skipped, so an interrupted run can resume.
  - immediate flush: each peptide is appended to the result TSV right after docking, not relying on a final one-shot write.
  - Vina timeout guard: per-peptide subprocess timeout 300s; on timeout record NA without blocking the batch.
"""
import os, subprocess, sys
from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import pybel

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # scripts/ -> repo root
DATA = os.path.join(REPO, "data")
DOCK = os.path.join(REPO, "docking")
QUEUE   = os.path.join(DATA, "moso_dock_queue_idppiv.txt")
REC     = os.path.join(DOCK, "1WCY_receptor.pdbqt")
BOX     = os.path.join(DOCK, "moso_box.txt")
# AutoDock Vina is a third-party binary (not under version control): defaults to `vina` on PATH,
# or set the absolute path via env var VINA_EXE (e.g. vina.exe under Windows).
VINA    = os.environ.get("VINA_EXE", "vina")
LIGDIR  = os.path.join(DOCK, "moso_ligands_idppiv")
RESULTS = os.path.join(DOCK, "moso_dock_results_idppiv.tsv")
os.makedirs(LIGDIR, exist_ok=True)

def _autolocate_babel_datadir():
    """Auto-derive BABEL_DATADIR from the installed openbabel package location
    (portable, avoids machine-specific absolute paths).
    This pipeline only uses openbabel for SDF->PDBQT format conversion, not MMFF94;
    it is set only to silence noise."""
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
    """RDKit makes 3D (SDF) -> openbabel converts to PDBQT (no make3D)"""
    m = Chem.MolFromSequence(seq)
    if m is None:
        return None
    m = Chem.AddHs(m)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    if not AllChem.EmbedMultipleConfs(m, numConfs=1, params=params):
        return None
    AllChem.MMFFOptimizeMoleculeConfs(m)   # RDKit ships MMFF94
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

# ---- resume from checkpoint: load finished ----
done = set()
if os.path.exists(RESULTS):
    with open(RESULTS, encoding="utf-8") as f:
        for ln in f:
            ln = ln.strip()
            if ln and not ln.startswith("peptide"):
                done.add(ln.split("\t")[0])

rows = [l.rstrip("\n").split("\t") for l in open(QUEUE, encoding="utf-8") if l.strip()]
print(f"iDPPIV queue: {len(rows)} entries, {len(done)} already done", flush=True)

# if first run, write header
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
        # immediate flush
        with open(RESULTS, "a", encoding="utf-8") as f:
            f.write(f"{seq}\t{float(sc):.3f}\t{aff if aff is not None else 'NA'}\n")
        done.add(seq)
        count += 1
        print(f"  {seq:7s} iDPPIV={float(sc):>6.2f}  dG={aff}", flush=True)
    except Exception as e:
        with open(RESULTS, "a", encoding="utf-8") as f:
            f.write(f"{seq}\t{float(sc):.3f}\tERR\n")
        done.add(seq)
        print(f"  {seq:7s} exception: {e}", flush=True)

print(f"\nthis run added {count} entries; cumulative {len(done)}/{len(rows)} -> {RESULTS}", flush=True)
