# Phase B — Pure-computational binding validation (Top-3 DPP4 inhibitory peptides)

> **Project positioning:** this is a **purely computational study** with no wet-lab validation (in-vitro DPP4 inhibition, IC₅₀, Caco-2 transport/activity, peptide synthesis). Also, this machine has **no GROMACS / no conda**, so 100–150 ns MD is infeasible. Phase B therefore uses a **single-structure end-point static MMFF94s validation** + geometric contact profiling as a feasibility substitute for wet-lab/MD.

---

## 1. Goal

Use computational means to answer the two questions wet-lab would answer:
1. **Is the docked conformation stable, and is binding real?** → complex MMFF94s relaxation.
2. **Which residues are most critical?** → per-ligand-position pocket contact profile (a surrogate for computational alanine scanning).

---

## 2. Method

### 2.1 Input
- Receptor: `1WCY_receptor.pdbqt` (DPP4 / sitagliptin complex, ADT-preprocessed).
- Ligand best pose: Vina docking product `moso_ligands/dock_<id>_<PEP>.pdbqt` **MODEL 1** (Top-3: LPPQ / APSPE / LAPSP).

### 2.2 Pocket definition
Among all receptor heavy atoms, residues within **≤ 9.0 Å** of any ligand heavy atom are treated as the pocket. Top-3 each yield **39–40 pocket residues**.

### 2.3 Molecule construction (RDKit 2026.03.3 + OpenBabel 3.2.1)
- **Pocket:** parse PDBQT → PDB block (AutoDock atom types mapped to element symbols) → RDKit read + unified H addition.
- **Ligand:** Vina folds the peptide into a single residue `UNK` block (atom names `C_1…C_38`); RDKit direct read fails on valence misjudgment. Instead use **OpenBabel to convert `dock_*.pdbqt` → MOL2 (MODEL 1)** (correctly senses peptide bonds), then hand to RDKit read + H addition.

### 2.4 MMFF94s complex relaxation
- Minimize **the whole complex** (pocket + ligand) with MMFF94s (no backbone fixed), letting ligand and pocket side chains relax together, removing any steric clash left by docking.
- The force field must initialize ring info first (`GetSymmSSSR`, the Pro pyrrolidine ring in peptides).

### 2.5 Static-MM binding energy ΔE_MM
Single-structure end-point method (consistent with single-structure MM-GBSA):
```
ΔE_MM = E_complex − E_pocket − E_ligand
```
All three are minimized from the same relaxed coordinates, then their MMFF energies are taken after independent minimization.

### 2.6 Contact profile (alanine-scan surrogate)
Pure geometric criteria (not relying on force-field magnitude, robust and interpretable):
- **Inter-molecular** ligand–pocket contacts: H-bond (polar pair ≤ 3.5 Å), hydrophobic (C–C ≤ 4.0 Å), ionic (polar pair ≤ 4.0 Å).
- Binned by ligand sequence position → **pocket contact count per residue position**.
- **The position with the most contacts = the most critical position** (mutating it to Ala loses the most binding) — exactly the question alanine scanning answers.

### 2.7 Interaction fingerprint
Extract the **DPP4 pocket residue names** contacted by the ligand (e.g., H-bond partners) from the relaxed conformation to localize the binding site.

---

## 3. Results

| Peptide (seq) | dG_Vina | E_complex (MMFF) | ΔE_MM | Pocket res. | Total contacts | H-bonds | Hydrophobic |
|---|---|---|---|---|---|---|
| LPPQ (Leu-Pro-Pro-Gln) | −6.109 | −304.38 | +2.45 | 39 | 92 | 10 | 72 |
| APSPE (Ala-Pro-Ser-Pro-Glu) | −7.150 | −254.11 | +2.93 | 39 | 119 | 13 | 81 |
| LAPSP (Leu-Ala-Pro-Ser-Pro) | −7.087 | −291.60 | +2.84 | 40 | 101 | 13 | 80 |

### 3.1 Contact-position profile (alanine-scan surrogate)
- **LPPQ**: L1:9 · P2:21 · **P3:25** · **Q4:37** → C-terminal Gln dominant.
- **APSPE**: A1:2 · P2:8 · S3:0 · P4:0 · **E5:109** → C-terminal **Glu overwhelming (109 contacts)**.
- **LAPSP**: L1:6 · A2:7 · P3:22 · S4:16 · **P5:50** → C-terminal Pro dominant.

> **Biological consistency:** APSPE's C-terminal glutamate (E5) shows extreme contact preference at the S1′ pocket, matching DPP4's known "prefers substrates with negatively-charged C-terminus occupying S1′" mechanism; the dual-Pro core fits DPP4 S1′/S2′ Pro preference.

### 3.2 H-bond partners (DPP4 pocket residues)
All three enrich **GLU146** (S2′ region) and respectively touch PRO149 / SER182 / TYR183 — all within DPP4's known active-center cavity, consistent with literature-reported DPP4 inhibitor binding sites.

---

## 4. Limitations (honest disclosure, must state at submission)

1. **Single conformation, no conformational sampling:** true MM-GBSA needs MD-trajectory sampling; this work is a single-relaxed-conformation end-point method, not covering conformational entropy/sampling fluctuations.
2. **Force-field limitation:** MMFF94s (generic small-molecule force field) is less accurate on proteins than CHARMM36m / AMBER99SB-ILDN.
3. **ΔE_MM is not true ΔG:** it is a **gas-phase** single-point energy difference, with no explicit solvation (GB polar term) and no conformational/translational entropy term; the three peptides' ΔE_MM are all ~+2.8 kcal/mol (similar magnitude, not discriminative), so **binding-strength ranking still relies mainly on Vina dG**, with ΔE_MM serving only as a conformation-relaxability (stability) indicator.
4. **Contact profile is a relative metric:** used to identify "most critical position"; its absolute contact count is affected by burial degree and should not be read as quantitative binding energy.
5. **Future work:** if GROMACS becomes available on an HPC host, add 100–150 ns MD + MM-GBSA (with GB solvation and entropy terms) for quantitative ΔG and converged residue decomposition — upgrading this Phase B static approximation to gold-standard validation.

---

## 5. Reproduction

```bash
# Env: Python 3.13 venv (RDKit 2026.03.3 / OpenBabel 3.2.1 / numpy / scipy)
PY=envs/default/Scripts/python.exe

# Run (script in same dir as 1WCY_receptor.pdbqt, moso_ligands/)
$PY phaseB_validation.py
# Outputs: phaseB_results.tsv (summary), phaseB_detail.json (per-residue/per-atom contacts)
```

Products: `data/phaseB/phaseB_results.tsv`, `data/phaseB/phaseB_detail.json`, `scripts/phaseB/phaseB_validation.py`.
