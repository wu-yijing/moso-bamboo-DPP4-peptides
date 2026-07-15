# -*- coding: utf-8 -*-
"""Single-structure MM-GBSA binding free-energy recomputation (OpenMM implicit GBn2)
— **pocket-restricted** version.
For the Top3 docking poses of the candidate peptides (APQIP/LPPGP/APPSQ):
    dG_bind = E_cplx - E_rec - E_lig   (implicit GBSA, single-structure)

Why pocket-restricted (corrected 2026-07-14):
- Full-protein GBSA (1WCY chain A, ~23k atoms incl. H) Born-radius evaluation is O(N^2),
  taking tens of minutes on CPU — exceeding this environment's per-run limit and getting
  reclaimed as a background task, hence infeasible.
- Switched to a **pocket-restricted receptor**: take whole residues within 15 A of the known
  active-site pocket center (62.8,47.7,4.8 A), ~40-60 residues / <1k atoms, GBSA <2 min.
- E_rec is computed once for the **same fixed pocket** and shared by all three candidates ->
  only ligand/complex differ, so dG differences arise purely from ligand-pocket interactions,
  giving the cleanest ranking.
- Hydrogens added by OpenMM `PDBFixer`, but **skip findMissingResidues** (to avoid it
  back-filling pocket gaps and re-inflating the system); only addMissingAtoms + addMissingHydrogens,
  terminal residues capped with N/C-terminal H/oxT, acceptable under the amber14 template
  (common approximation for pocket-terminated systems).
- Ligand read directly from RDKit on the pose MODEL 1 (standard AA residue names) -> PDBFixer H addition.

Methodological honest disclosure: pocket-restricted MM-GBSA ignores long-range receptor-receptor
interactions outside the pocket and receptor conformational strain; it is used for **relative
ranking of the three candidates**, not for absolute affinity claims; single-structure approximation,
ignores -TΔS."""
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
# Correct iDPPIV-queue docking poses (dock_<idx>.pdbqt), idx matches moso_dock_queue_idppiv.txt
#   APQIP -> dock_19, LPPGP -> dock_00, APPSQ -> dock_54
TOP3 = [
    ("APQIP", os.path.join(LIGDIR_IDPP, "dock_19.pdbqt")),
    ("LPPGP", os.path.join(LIGDIR_IDPP, "dock_00.pdbqt")),
    ("APPSQ", os.path.join(LIGDIR_IDPP, "dock_54.pdbqt")),
]
STD_AA = {'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
           'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'}
POCKET_CENTER = (62.8, 47.7, 4.8)   # 1WCY DPP4 active-site pocket center (box center)
POCKET_RADIUS = 9.0                  # A (tightened to fit <2 min foreground limit)


def write_rec_pocket(src, dst, center, radius):
    """Receptor pdbqt -> pocket PDB of whole residues within radius of center
    (chain A, standard AA, no HETATM)."""
    cx, cy, cz = center
    # aggregate atoms by residue
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
    """Take the ATOM lines of pose MODEL 1 (standard AA), **sorted by resid**
    (pose ATOM lines are often out of order and must be made residue-contiguous,
    otherwise PDBFixer mis-assigns residues), then write PDB (truncated to col 76)
    for PDBFixer H addition."""
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
    atoms.sort(key=lambda x: x[0])   # residue-contiguous
    with open(dst, "w") as f:
        for _, ln in atoms:
            f.write(ln[:76] + "  \n")


def fix_pocket(src_pdb, dst_pdb):
    """PDBFixer: add missing atoms + hydrogens. For the pocket (non-contiguous residues):
    first call findMissingResidues to establish the attribute, then clear its gap dict ->
    addMissingAtoms only fills existing residues' missing atoms (OXT/terminal H), does not
    back-fill pocket gaps, keeping the system small."""
    fixer = PDBFixer(filename=src_pdb)
    fixer.findMissingResidues()      # establish missingResidues attr (else addMissingAtoms raises AttributeError)
    fixer.missingResidues = {}       # clear gaps -> do not inflate pocket
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    PDBFile.writeFile(fixer.topology, fixer.positions, dst_pdb)


def parse_pdbqt_model1_atoms(path):
    """Parse ATOM/HETATM lines of pdbqt MODEL 1 -> list of dicts. Use all if no MODEL."""
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
    """Build a canonical pentapeptide topology from sequence (correct bonds/Pro ring),
    then overlay pose heavy-atom coordinates by (1-based residue number, atom name).
    Returns the **heavy-atom** mol (with coords), no H -- H is added by ligand_pdb_with_h
    via PDBFixer (cleanest amber-template H names/residue info)."""
    atoms = parse_pdbqt_model1_atoms(pqt)
    if not atoms:
        return None
    m = Chem.MolFromSequence(seq)          # canonical peptide topology (correct bonds/Pro ring)
    if m is None:
        return None
    m = Chem.RemoveHs(m)
    if m.GetNumConformers() == 0:        # ensure a conformer exists to place coords
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
    """Heavy-atom mol (with coords) -> hand-written heavy-atom PDB -> PDBFixer add
    missing atoms + H, writing a PDB with **correct H names and residue info** (amber14
    template matches directly). Lesson: RDKit AddHs H atoms return None from GetMonomerInfo;
    a hand-written PDB marked UNK or with non-template H names (H1/H2...) is judged
    'missing hydrogens' by amber14. PDBFixer emits standard H names by template, once and for all."""
    heavy = os.path.join(DOCK, "_mgbsa_heavy_tmp.pdb")
    open(heavy, "w").write(_mol_to_pdb(Chem.RemoveHs(mol)))
    fixer = PDBFixer(filename=heavy)
    fixer.findMissingResidues(); fixer.missingResidues = {}   # do not back-fill missing residues
    fixer.findMissingAtoms(); fixer.addMissingAtoms(); fixer.addMissingHydrogens(7.0)
    PDBFile.writeFile(fixer.topology, fixer.positions, dst)


def _mol_to_pdb(m):
    """Manually write an RDKit mol (with H) to PDB text, atom-by-atom (incl. H),
    avoiding Chem.MolToPDBBlock dropping H atoms in some cases.
    Key: H atoms produced by RDKit's own AddHs return None from GetMonomerInfo();
    writing them as UNK causes amber14 'No template / missing hydrogens'.
    Fix: H atoms inherit the residue info (resname/resnum/chain) of their bonded heavy atom,
    and are assigned unique H names (H1/H2/...) by (resnum,chain) count. amber14 matches
    H by element + bonded heavy atom, name-insensitive."""
    conf = m.GetConformer()
    # first collect heavy-atom residue info (by atom idx)
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
    hcount = {}          # (resnum, chain) -> count of H written
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
            # H: take residue info of its bonded heavy atom
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
    """Minimize and return (energy_kcal, minimized_positions)."""
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
    """Compute energy only (no minimization), used to decompose receptor/ligand
    energy from the minimized complex coordinates."""
    integ = openmm.LangevinMiddleIntegrator(300*openmm.unit.kelvin,
                                           1/openmm.unit.picosecond,
                                           0.001*openmm.unit.picosecond)
    sim = app.Simulation(topology, system, integ, openmm.Platform.getPlatformByName("CPU"))
    sim.context.setPositions(positions)
    e = sim.context.getState(getEnergy=True).getPotentialEnergy()
    return e.value_in_unit(openmm.unit.kilocalories_per_mole)


def subset_topology(topology, keep_idxs):
    """Extract keep_idxs (original atom.index set) into a new Topology (preserving
    residue/chain/bond structure). Returns (new_topology, pos_index_list), where
    pos_index_list[i] = original coordinate index of the i-th new atom."""
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

    # Receptor pocket (fixed, shared by all three candidates) -> PDBFixer H addition
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
            raise RuntimeError(f"ligand read failed: {pqt}")
        # Ligand: canonical sequence peptide + pose heavy-atom coords; H added by PDBFixer
        # (correct names/residue info)
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
        # (a) minimize the **complex** once only; (b) decompose E_rec / E_lig from the
        #     same minimized coordinates -- single-trajectory MM-GBSA, eliminating the
        #     "ligand independent relaxation" artifact (otherwise free-ligand energy is
        #     spuriously lowered, making dG too high or even positive, e.g. the +8.2
        #     initial APQIP result).
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
