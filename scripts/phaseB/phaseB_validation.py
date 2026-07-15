# -*- coding: utf-8 -*-
"""
Phase B — purely computational binding validation (no GROMACS / no MD feasibility alternative)

Methodology (honest disclosure):
  This environment has no GROMACS / conda, so 100-150 ns MD + MM-PBSA cannot run.
  This script uses a **single-structure end-point MM validation**
  (single-structure end-point MM) on the Top3 peptides' best docking conformation (Vina MODEL 1):
    1) complex MMFF94s relaxation (receptor backbone fixed, ligand + pocket side-chains optimized)
    2) static MM binding energy ΔE = E_complex - E_pocket - E_ligand (gas-phase, MMFF94s)
       + nonpolar solvation correction ΔG_SA = γ·ΔSASA  (γ=0.0072 kcal/mol/Å²)  [if SASA computable]
    3) per-residue interaction-energy decomposition (vdW + Coulomb, MMFF force field)
       -- as a quantitative substitute for "computational alanine scanning", answering "which residues are most critical"
    4) interaction fingerprint: H-bonds / hydrophobic contacts / ionic bonds -> locate DPP4 pocket residues

Limitations (written into the report):
  - single conformation, no conformational sampling (only MD provides that)
  - solvation is approximate (nonpolar SASA term only), no explicit GB polar term, no entropy term
  - force field is MMFF94s (general-purpose small-molecule FF, lower protein accuracy than CHARMM36m/AMBER99SB-ILDN)
  - conclusion is a relative "binding-strength ranking + key-residue identification", not absolute ΔG quantification

Dependencies: RDKit 2026.03.3 (envs/default). Reads the idppiv-queue Vina poses directly via RDKit
       (standard AA residue names, no OpenBabel needed).
Run: envs/default/Scripts/python.exe phaseB_validation.py
"""
import os, sys, subprocess, math, json
from collections import defaultdict

BASE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable
# must include .exe: Python subprocess.run does not go through a shell, so it will not
# auto-append a suffix by PATHEXT the way Git Bash does
OBABEL = os.path.join(os.path.dirname(PY), "obabel.exe")

# ---------- repo root (scripts/phaseB/ -> repo root) ----------
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# ---------- Top3 peptide definition (Plan α: best-binding peptides under the new iDPPIV-queue RDKit prep) ----------
# Correct docking poses come from the new iDPPIV-queue Vina output (moso_ligands_idppiv/dock_<idx>.pdbqt),
# where <idx> matches data/moso_dock_queue_idppiv.txt:
#   APQIP -> index 19 (dock_19.pdbqt), LPPGP -> index 0 (dock_00.pdbqt),
#   APPSQ -> index 54 (dock_54.pdbqt).
# These poses carry standard AA residue names (ALA/PRO/GLN/ILE/SER/GLY/LEU) and are read directly by RDKit,
# no OpenBabel needed. dG taken from moso_dock_results_idppiv.tsv (one-to-one with the pose).
TOP3 = [
    {"name": "APQIP", "seq": "APQIP",
     "dock": os.path.join(REPO, "docking", "moso_ligands_idppiv", "dock_19.pdbqt"), "dg": -6.807},
    {"name": "LPPGP", "seq": "LPPGP",
     "dock": os.path.join(REPO, "docking", "moso_ligands_idppiv", "dock_00.pdbqt"), "dg": -6.558},
    {"name": "APPSQ", "seq": "APPSQ",
     "dock": os.path.join(REPO, "docking", "moso_ligands_idppiv", "dock_54.pdbqt"), "dg": -6.513},
]
RECEPTOR = os.path.join(REPO, "docking", "1WCY_receptor.pdbqt")
POCKET_CUTOFF = 9.0   # Å, ligand heavy atoms within 8-10 Å count as pocket
DIELECTRIC = 4.0       # Coulomb dielectric constant (protein-interior approximation)
GAMMA_SA = 0.0072      # kcal/mol/Å² nonpolar solvation surface tension
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

# AutoDock atom type -> element (in proteins CA is the Cα carbon; first letter suffices)
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
    """Parse PDB / PDBQT ATOM/HETATM lines. model=None -> all; else the block after the specified MODEL.
    The last column is the AutoDock atom type, mapped to an element symbol.
    Compatible with both single-conformation (no MODEL record) and multi-conformation Vina output."""
    # pre-scan: does the file contain a MODEL record?
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
    """Read the Vina docking pose's MODEL 1 directly via RDKit (the idppiv-queue pose
    carries standard AA residue names ALA/PRO/GLN/ILE/SER/GLY/LEU, which RDKit
    perceives correctly as peptide bonds, no OpenBabel needed). Returns the H-added mol."""
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
    # initialize ring info (MMFF needs it; e.g. Pro pyrrolidine ring in peptides)
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
            # position constraint (tolerance 0.3 Å, strong constant) ~= fix backbone heavy atoms at current coords
            fld.MMFFAddPositionConstraint(int(i), 0.3, 1.0e6)
    conv = fld.Minimize(maxIts=maxIts)
    return fld.CalcEnergy(), conv

def coords_of(mol, idxs):
    conf = mol.GetConformer()
    return [(conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z) for i in idxs]

def _set_coords(mol, coords):
    """Overwrite coordinates by atom order (mol atom count must == len(coords))."""
    conf = mol.GetConformer()
    for j, (x, y, z) in enumerate(coords):
        conf.SetAtomPosition(j, (x, y, z))

def _set_heavy_coords(mol, coords):
    """Overwrite heavy-atom coordinates only (H preserved, rebuilt later by AddHs)."""
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
    """Return the peptide's MM-binding-energy-related quantities + complex mol reference
    (for subsequent decomposition / fingerprint)."""
    rec_path = os.path.join(BASE, RECEPTOR)
    dock_path = os.path.join(BASE, pep["dock"])
    print(f"  [{pep['name']}] parsing receptor ({RECEPTOR}) ...")
    rec_atoms = parse_pdb_atoms(rec_path)
    print(f"  [{pep['name']}] parsing ligand best conformation (MODEL 1) ...")
    lig_atoms = parse_pdb_atoms(dock_path, model=1)
    if not lig_atoms:
        raise RuntimeError("ligand atom parsing failed")
    print(f"  [{pep['name']}] receptor {len(rec_atoms)} atoms, ligand {len(lig_atoms)} atoms")

    pocket_keys = find_pocket(rec_atoms, lig_atoms, POCKET_CUTOFF)
    pkt_atoms = pocket_atoms_subset(rec_atoms, pocket_keys)
    print(f"  [{pep['name']}] pocket residue count = {len(pocket_keys)} (cutoff {POCKET_CUTOFF} Å)")

    # strip PDBQT's own H (to avoid duplication with AddHs), then add H uniformly
    pkt_heavy = [a for a in pkt_atoms if a["element"] != "H"]
    lig_heavy = [a for a in lig_atoms if a["element"] != "H"]

    # build pocket mol (protein, parsed by column) and ligand mol (rebuilt from sequence + dock coords)
    pocket_mol = build_mol_from_atoms(pkt_heavy, add_hs=True)
    lig_mol = ligand_from_pose(pep["dock"], pep["name"])
    if pocket_mol is None or lig_mol is None:
        raise RuntimeError("RDKit pocket/ligand construction failed")
    print(f"  [{pep['name']}] pocket mol {pocket_mol.GetNumAtoms()} atoms, ligand mol {lig_mol.GetNumAtoms()} atoms")

    combo = Chem.CombineMols(pocket_mol, lig_mol)
    n_pkt = pocket_mol.GetNumAtoms()        # pocket total atoms incl. H
    n_pkt_heavy = len(pkt_heavy)            # pocket heavy atoms
    n_lig_heavy = len(lig_heavy)            # ligand heavy atoms
    # fix pocket backbone heavy atoms (pocket heavy-atom indices 0..n_pkt_heavy-1 in combo)
    fixed = [i for i, a in enumerate(pkt_heavy) if a["name"] in BACKBONE]
    print(f"  [{pep['name']}] fixed pocket backbone atoms = {len(fixed)} (used for pocket-only minimization)")
    # relax the whole complex (no backbone fix) so ligand + pocket side-chains relax together, removing dock clashes
    E_complex, _ = mmff_minimize(combo, fixed_idx=None, maxIts=1500)
    print(f"  [{pep['name']}] E_complex (MMFF) = {E_complex:.2f} kcal/mol")

    # extract minimized **heavy-atom** coords (to rebuild separate mols)
    pkt_heavy_idx = list(range(n_pkt_heavy))                    # pocket heavy atoms
    lig_heavy_idx = list(range(n_pkt, n_pkt + n_lig_heavy))   # ligand heavy atoms
    pkt_coords = coords_of(combo, pkt_heavy_idx)
    lig_coords = coords_of(combo, lig_heavy_idx)

    # pocket alone: rebuild heavy atoms -> write minimized coords -> add H -> minimize
    pkt_only = build_mol_from_atoms(pkt_heavy, add_hs=False)
    _set_coords(pkt_only, pkt_coords)
    pkt_only = Chem.AddHs(pkt_only, addCoords=True)
    fixed2 = [j for j, a in enumerate(pkt_heavy) if a["name"] in BACKBONE]
    E_pocket, _ = mmff_minimize(pkt_only, fixed_idx=fixed2, maxIts=1000)

    # ligand alone: clean MOL2 structure -> write minimized coords -> add H -> minimize
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
        # for decomposition / fingerprint: minimized combo reference + indices (heavy atoms)
        "_combo": combo, "_pkt_idx": pkt_heavy_idx, "_lig_idx": lig_heavy_idx,
        "_pkt_atoms": pkt_heavy, "_lig_atoms": lig_heavy,
        "_fixed": fixed,
    }

# ----------------------------------------------------------------------
def contact_profile(res):
    """Per-ligand-residue pocket contact count (H-bond / hydrophobic / ionic), as a
    robust substitute metric for computational alanine scanning. Contact judgment is
    purely geometric (independent of force-field energy magnitude), interpretable:
    the most-contacted ligand position = the most critical position (largest binding
    loss when mutated to Ala)."""
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

    # merge atoms into residues by sequence position (new residue at backbone N; cap = sequence length)
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
    # extra N groups within a side chain (e.g. Gln NE2) cause surplus groups -> merge into last group
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
    """H-bond / hydrophobic / ionic: the list of ligand-pocket contacting residues."""
    combo = res["_combo"]
    pkt_idx = res["_pkt_idx"]
    lig_idx = res["_lig_idx"]
    conf = combo.GetConformer()

    # element / charge classification
    def is_donor(idx):
        a = combo.GetAtomWithIdx(idx)
        if a.GetSymbol() == "N":
            # has H -> donor
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
            # hydrophobic
            if is_carbon(gi) and is_carbon(pi) and d <= 4.0:
                hydrophobics.append((gi, pi, d))
            # H-bond candidate: N/O polar contact <3.5 Å
            if (is_donor(gi) or is_acceptor(gi)) and (is_donor(pi) or is_acceptor(pi)) and d <= 3.5:
                hbonds.append((rank, gi, pi, d))
    # de-duplicate H-bonds (by ligand atom), and record the corresponding pocket residue
    seen = set()
    hb_unique = []
    for rank, gi, pi, d in hbonds:
        if rank not in seen:
            seen.add(rank)
            # locate the residue the pocket atom belongs to (from pkt_heavy list)
            pkt = res["_pkt_atoms"]
            # pi is the pocket heavy-atom index in combo; map to pkt_heavy name
            lig_name = res["_lig_atoms"][rank]["name"]
            hb_unique.append((rank, lig_name, pi, d))
    # pocket residue-name map: combo index -> pkt_heavy[rank] residue label
    pkt = res["_pkt_atoms"]
    hb_residues = []
    for rank, lig_name, pi, d in hb_unique:
        # pi is in the pocket segment 0..n_pkt_heavy-1 of combo, corresponding to pkt_heavy[pi]
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
    """Single-peptide result -> JSON-serializable pure-quantity dict (strip RDKit mol refs)."""
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
    """Single-peptide independent-process entry: compute + contact fingerprint,
    write phaseB_<NAME>.json.

    Independent-process reason: repeated RDKit MMFF mmff_minimize calls in the same
    process trigger a C-layer segfault (no Python traceback) on the second peptide;
    splitting each peptide into its own subprocess is stable.
    """
    print(">>> processing %s (Vina dG=%.3f kcal/mol)" % (pep["name"], pep["dg"]))
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
    print("Phase B — static MM binding validation (Top3 DPP4 inhibitory peptides)")
    print("=" * 70)
    # each peptide in its own subprocess (isolate RDKit MMFF same-process state corruption that crashed peptide #2)
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
                print("  [ERROR] %s subprocess failed (rc=%d)" % (pep["name"], r.returncode))
                continue
            jp = os.path.join(BASE, "phaseB_%s.json" % pep["name"])
            if os.path.exists(jp):
                with open(jp) as f:
                    results.append(json.load(f))
                os.remove(jp)  # clean up intermediate
        except Exception as e:
            import traceback
            traceback.print_exc()
            print("  [ERROR] %s scheduling failed: %s" % (pep["name"], e))
    # output
    out_tsv = os.path.join(BASE, "phaseB_results.tsv")
    with open(out_tsv, "w") as f:
        f.write("peptide\tseq\tdG_Vina\tdE_complex\tdE_MM\tn_pocket_res\tcontact_total\tHbond\thydrophobic\n")
        for r in results:
            fp = r["fingerprint"]
            f.write("%s\t%s\t%.3f\t%.2f\t%.3f\t%d\t%d\t%d\t%d\n" % (
                r["peptide"], r["seq"], r["dg_vina"], r["E_complex"], r["dE_MM"],
                r["n_pocket_res"], r["contact_total"],
                fp["hbonds"], fp["hydrophobics"]))
    print("\n[OK] main result written to %s" % out_tsv)
    # detailed JSON
    with open(os.path.join(BASE, "phaseB_detail.json"), "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print("[OK] detailed info written to phaseB_detail.json")

    print("\n=== summary ===")
    for r in results:
        fp = r["fingerprint"]
        pp = r["contact_per_position"]
        pp_str = " ".join("%s%d:%d" % (L, i, c) for i, L, c in pp)
        print("  %-6s dG_Vina=%7.3f  dE_MM=%7.3f  pocket_res=%3d  contact=%5d  Hbond=%d  Hbond_res=[%s]"
              % (r["peptide"], r["dg_vina"], r["dE_MM"], r["n_pocket_res"],
                 r["contact_total"], fp["hbonds"],
                 ",".join(fp["hbond_pocket_residues"][:6])))
        print("          contact-position profile: %s" % pp_str)


if __name__ == "__main__":
    # single-peptide mode: called by subprocess, compute and write phaseB_<NAME>.json then exit
    if len(sys.argv) > 2 and sys.argv[1] == "--pep":
        _name = sys.argv[2]
        _single = [p for p in TOP3 if p["name"] == _name]
        if _single:
            _run_single(_single[0])
        sys.exit(0)
    main()
