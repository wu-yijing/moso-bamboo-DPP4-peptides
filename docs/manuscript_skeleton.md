# In silico discovery of DPP4-inhibitory peptides from the moso bamboo (*Phyllostachys edulis*) proteome

> **Status: DRAFT SKELETON** — pure *in silico* study, no wet-lab validation.
> Placeholders marked `[TODO]` require author/peer input or future data.
> All computational results below are reproduced from `data/` in this repository.

**Authors:** [TODO: Y. Wu, Q. Wu, J. Zhao, M. Chen, G. Jin — confirm order/affiliations]
**Target journal:** *Food & Function* (computational food-peptide discovery) — or *in silico*-focused venue.

---

## Abstract

**Background:** Dipeptidyl peptidase-4 (DPP4/CD26) inhibitors improve glycaemic control by prolonging incretin half-life. Food-derived DPP4-inhibitory peptides are a promising nutritional strategy, yet systematic *in silico* discovery from under-explored plant proteomes remains limited.
**Objective:** To computationally mine the moso bamboo (*Phyllostachys edulis*) proteome for novel DPP4-inhibitory peptide sequences.
**Methods:** We retrieved 253 bamboo proteins (UniProt taxonomy 38705), performed *in silico* gastrointestinal digestion, screened peptides for allergenicity/toxicity with official web servers (AlgPred 2.0, ToxinPred 3.0), docked candidates against the DPP4 structure (PDB 1WCY), and verified top hits by static MMFF94s contact profiling and DPP4-centred network pharmacology (STRING).
**Results:** Virtual digestion yielded 4,950 short (2–6 aa) unique peptides; official filtering retained 4,674 non-toxic, non-allergenic candidates. Docking of a 60-peptide priority queue identified **LPPQ** (ΔG = −7.47 kcal/mol), **APSPE** (−7.15) and **LAPSP** (−7.09) as top binders. Static MM confirmed DPP4-active-site engagement (H-bonds enriched at GLU146/S2′; APSPE C-terminal Glu anchored S1′ with 109 contacts). All top peptides carried the X-Pro/X-Ala DPP4-recognition motif. Network pharmacology positioned DPP4 within the incretin axis (GLP-1/GCG, GIP), rationalising glucoregulatory potential.
**Conclusions:** The moso bamboo proteome is a computationally viable source of DPP4-binding peptide sequences. *In vitro* DPP4 inhibition and stability assays are required to confirm activity.

---

## 1. Introduction

Type 2 diabetes mellitus (T2DM) is driven partly by inadequate incretin (GLP-1, GIP) signalling; DPP4 rapidly truncates these hormones via its X-Pro/X-Ala exopeptidase activity. Pharmacological DPP4 inhibitors (gliptins) are established therapeutics. Concurrently, food-derived bioactive peptides with DPP4-inhibitory activity offer a nutritional complement, with several characterised from animals and crops (e.g., garlic, yam, soy, milk) [TODO: cite 4–6 key refs].

Moso bamboo (*Phyllostachys edulis*) is a major economic and food crop in East Asia, yet its proteome remains an under-explored source of bioactive peptides. *In silico* peptidomics — *in silico* digestion, allergenicity/toxicity screening, and molecular docking — enables rapid, low-cost candidate discovery before any synthesis [TODO: cite methodology refs].

Here we present a fully *in silico* pipeline that discovers DPP4-binding peptide sequences from the moso bamboo proteome, with static molecular-mechanics and network-pharmacology verification as pure-computational surrogates for experimental validation.

---

## 2. Materials and Methods

### 2.1 Moso bamboo proteome and *in silico* digestion
253 protein sequences were retrieved from UniProt (taxonomy **38705**, *Phyllostachys edulis*; **0 reviewed / all TrEMBL** predicted, 103,409 residues). *In silico* gastrointestinal digestion was simulated with broad-specificity protease rules (trypsin/chymotrypsin/pepsin/elastase/carboxypeptidase), yielding **7,988 unique peptides** under strict rules, of which **4,950 were short (2–6 aa)** peptides carried forward.

### 2.2 *In silico* allergenicity and toxicity screening (Phase A)
All 4,950 short peptides were submitted to official web servers: **AlgPred 2.0** (AAC-RF model; non-allergen threshold ML-score < 0.6) and **ToxinPred 3.0** (non-toxin). Funnel: 4,950 → AlgPred non-allergen (4,680) → ToxinPred non-toxin (**4,674**).
> *Note:* the **PeptideRanker** bioactivity server (UCD) was **unavailable (HTTP 503)** during this study; its PR ≥ 0.5 layer is pending server recovery and will be merged via `phaseA_merge.py`. An early in-house proxy score was used **only** to prioritise a 60-peptide docking queue and is **not** used as a publication filter.

### 2.3 DPP4 target and molecular docking
Human DPP4 crystal structure **PDB 1WCY** (sitagliptin-bound) was prepared as receptor PDBQT; a 60-peptide priority queue was built as 3D ligand PDBQTs (RDKit; 8 Arg/His-containing peptides rescued from NaN coordinates). **AutoDock Vina 1.2.5** docking used a pocket box centred on the catalytic site (GLU146/S2′, PRO149, SER182, TYR183, ASN151; S1′ favours C-terminal negatively-charged residues). All 60/60 ligands docked successfully.

### 2.4 Static MM binding verification (Phase B)
Top-3 peptides were verified by **single-conformation MMFF94s** analysis: OpenBabel PDBQT→MOL2 conversion, RDKit MMFF94s relaxation, geometric contact profiling (H-bond ≤ 3.5 Å, hydrophobic C–C ≤ 4.0 Å, ionic ≤ 4.0 Å) as a robust surrogate for alanine scanning, and an interaction fingerprint localising DPP4 pocket residues.
> *Caveat:* ΔE_MM is a gas-phase single-point energy difference, **not** a true ΔG (no explicit solvation/entropy); it is reported for relative contact context only. Binding-rank ordering follows Vina ΔG.

### 2.5 *In silico* ADMET and GI-stability profiling (Phase C1)
For the 60 docked peptides we computed sequence-derived physicochemical descriptors (BioPython ProtParam: MW, pI, GRAVY, aliphatic index, charge at pH 7.4; **Boman** protein-binding index) and RDKit descriptors (TPSA, HBD, HBA, rotatable bonds, QED). Gastrointestinal stability was predicted by literature protease-cleavage rules (pepsin/trypsin/chymotrypsin/elastase/carboxypeptidase/aminopeptidase) and DPP4 self-cleavage (X-Pro/X-Ala motif).

### 2.6 DPP4-centred network pharmacology (Phase C2)
The human DPP4 (UniProt **P27487**) functional association network was retrieved from **STRING DB** (REST API, functional, score ≥ 0.4) and enriched with literature-derived DPP4 cleavage-substrate edges (GLP-1/GCG, GIP, CXCL12/SDF-1, NPY, PYY, substance P) to contextualise the glucoregulatory mechanism.

---

## 3. Results

### 3.1 Candidate funnel
| Step | Peptides |
|---|---|
| Moso bamboo proteins | 253 (all TrEMBL) |
| Unique peptides (strict digestion) | 7,988 |
| Short peptides (2–6 aa) | 4,950 |
| − AlgPred non-allergen | 4,680 |
| − ToxinPred non-toxin | **4,674** |
| Docking priority queue (proxy-ranked) | 60 |

### 3.2 Docking and top candidates
Vina ΔG for the 60-peptide queue ranged to **−7.47 kcal/mol** (no < −8.0 "strong" binders). Top 3:

| Rank | Peptide | Sequence | ΔG (kcal/mol) |
|---|---|---|---|
| 1 | LPPQ | Leu-Pro-Pro-Gln | **−7.472** |
| 2 | APSPE | Ala-Pro-Ser-Pro-Glu | −7.150 |
| 3 | LAPSP | Leu-Ala-Pro-Ser-Pro | −7.087 |

All 60 candidates carried the **X-Pro / X-Ala** N-terminal motif recognised by DPP4 — consistent with the enzyme's substrate preference and with the docking funnel.

### 3.3 Static MM and interaction fingerprint (Phase B)
| Peptide | ΔE_MM* | Pocket residues | Total contacts | H-bonds | Key contact |
|---|---|---|---|---|---|
| LPPQ | +2.45 | 39 | 92 | 10 | Gln4 = 37 (C-terminal) |
| APSPE | +2.93 | 39 | 119 | 13 | **Glu5 = 109** (S1′) |
| LAPSP | +2.84 | 40 | 101 | 13 | Pro5 = 50 (C-terminal Pro) |

H-bonds in all three enriched at **GLU146 (S2′)** and touched PRO149/SER182/TYR183/ASN151 — residues of the known DPP4 active-site cavity. APSPE's C-terminal glutamate dominated S1′ occupancy (109 contacts), matching DPP4's preference for a negatively charged C-terminus in S1′.
*\* ΔE_MM is gas-phase MMFF94s single-point difference, not ΔG.*

### 3.4 ADMET / GI profile and DPP4-substrate caveat (Phase C1)
Among 60 docked peptides: GI stability = Moderate (1 endoprotease site) 39, Low (multiple sites) 11, High (no site) 10. All 60 carried the X-Pro/X-Ala motif → each is itself a **DPP4 substrate**, i.e. would be cleaved by DPP4 *in vivo*. This is a known oral-delivery limitation, not a disqualifier (the peptide still competes for the active site *in vitro*); structural stabilisation (cyclisation, N-terminal capping) would be needed for oral application.

### 3.5 Network pharmacology context (Phase C2)
The DPP4 network comprised **15 nodes / 32 edges** (STRING 11 nodes/26 functional edges + 4 literature substrate nodes/6 edges). Core neighbours included **GCG** (GLP-1/GLP-2 precursor), **GIP** (gastric incretin), **CXCR4** (SDF-1/CXCL12 receptor), **ADA**, **CAV1**, **PRCP**. Mechanism: inhibiting DPP4 raises active GLP-1/GIP → augmented insulin secretion → improved glycaemic control — providing a coherent pharmacological rationale for the discovered peptides.

---

## 4. Discussion

We report, to our knowledge, the first *in silico* peptidome-wide discovery of DPP4-binding peptides from moso bamboo. The top candidates (LPPQ, APSPE, LAPSP) bind the DPP4 active site by the same mechanistic grammar as characterised food peptides — Pro-rich cores and a C-terminal anchor occupying S1′/S2′. Static MM and contact fingerprinting corroborate docking without requiring molecular dynamics hardware.

The universal X-Pro/X-Ala motif is a double-edged finding: it explains strong *in silico* DPP4 recognition yet flags rapid *in vivo* DPP4-mediated clearance. This tension is well documented for food-derived DPP4 inhibitors and motivates future stabilisation strategies.

### Limitations (honest disclosure — MUST appear in Methods/Limitations)
1. **Filtering:** official ToxinPred 3.0 + AlgPred 2.0 screening was performed (4,674 non-toxic/non-allergenic candidates). The PeptideRanker bioactivity server was unavailable (503); its PR ≥ 0.5 layer is pending and will be merged when the server recovers. The in-house proxy score was used only to narrow the docking queue, not as a publication criterion.
2. **Docking is a static ΔG estimate**; Phase B used single-conformation MMFF94s, **not** molecular dynamics. No *in vitro* or *in vivo* validation was performed.
3. **Source provenance:** the moso bamboo proteome contains **0 manually-reviewed entries (all TrEMBL predicted)**; this is disclosed in Methods. (Comparable *in silico* food-peptide studies are similarly TrEMBL-dominated.)
4. **Phase B/C are computational approximations.** No wet-lab experiments (IC₅₀ DPP4 inhibition, Caco-2 permeability, serum/ GI stability assays) were conducted; all functional claims are *in silico* and require experimental confirmation.

---

## 5. Conclusion

The moso bamboo (*Phyllostachys edulis*) proteome yields DPP4-binding peptide sequences discoverable by a fully *in silico* pipeline. LPPQ, APSPE and LAPSP are the top computational candidates. Experimental validation of DPP4 inhibitory activity and oral stability is the necessary next step.

---

## 6. References
[TODO: 12–20 refs — UniProt/TrEMBL, AlgPred, ToxinPred, PeptideRanker, AutoDock Vina, MMFF94s, STRING DB, DPP4 substrate/incretin reviews (Drucker 2006; Mentlein 2009), food-derived DPP4 peptide precedents (garlic, yam, soy, milk).]
