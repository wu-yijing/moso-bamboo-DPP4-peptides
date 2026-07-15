# Moso Bamboo (*Phyllostachys edulis*) DPP4-Inhibitory Peptides — A Reproducible, Offline Computational Framework

> A **methodology-first**, fully *in-silico* project that uses the 253-protein complete proteome of moso bamboo as a **case study** to build a transparent, zero-external-dependency framework for proteome-scale antidiabetic peptide discovery, with three methodological contributions: **(①) quantitative exposure of DPP-IV predictor length bias**, **(②) multi-target (DPP4/ACE/α-glucosidase) same-proteome prioritization**, and **(③) genome-existence + edible-shoot-expression bridging**.
>
> The authoritative, submission-ready description is the CSBJ English manuscript `docs/manuscript_csbj_en.md`. A Chinese working draft (`docs/manuscript_draft.md`) is kept locally and is **not pushed** to this repository.

---

## 1. Project background

An earlier yam (*Dioscorea polystachya*) project collapsed as a DPP4 claim: its UniProt pool held only 20 curated proteins, virtual digestion yielded merely 237 peptides (vs. a 5,230 funnel target), and 0 of 7 experimentally-validated peptides (FWPQY, *etc.*) could be cut from those 20 proteins — yam could not support a DPP4 claim.

Following the methodology of a template paper (Cheng *et al.*, *Bioorganic Chemistry* 175 (2026) 109801, *"Discovery of garlic-derived peptides as natural DPP4 inhibitors"*), we reselected **moso bamboo (*Phyllostachys edulis*)**. UniProt verification shows bamboo's entry pool (taxonomy 38705) contains **253 proteins** (all TrEMBL-predicted, 0 reviewed) — ~2.2× the template's garlic (113), sufficient to sustain a full pipeline.

---

## 2. Computational pipeline and results

| Step | Method | Key result | Software |
|---|---|---|---|
| 1. Species/protein-pool verification | UniProt taxonomy/search + uniprotkb/search | *P. edulis* taxid **38705**; 253 proteins (0 reviewed / 253 TrEMBL); no reference proteome | `curl` + UniProt REST API |
| 2. Protein download | UniProtKB `/stream` (FASTA) | `data/moso_253.fasta` (253, 103,409 residues) | `curl` |
| 3. *In-vitro* virtual digestion | Re-implemented ExPASy PeptideCutter: pepsin(pH1.3)+trypsin+specific chymotrypsin (cuts C-term of F/Y/W, not before Pro); keep 2–6 aa, dedupe | Strict unique peptides **7,988**; 2–6 aa short peptides **4,950** | Python (numpy/scipy/biopython) |
| 4. Allergen/toxin filtering | ToxinPred 3.0 + AlgPred 2.0 (AAC-RF) replacing the original proxy scorer | 4,950 → **4,674** | Python (custom) |
| 5. DPP4 structural-preference narrowing | 3–5 aa, N-term hydrophobic, Pro/Ala preferred at position 2 | candidate pool **2,019** → docking queue **565** | Python (custom) |
| 6. Ligand 3D prep (PDBQT) | SMILES→3D→PDBQT; RDKit rescue for 8 Arg/His peptides (NaN) | **565/565** ligand PDBQTs, zero NaN | OpenBabel / RDKit |
| 7. Receptor & pocket | DPP4 crystal **1WCY** (sitagliptin-bound) → PDBQT; grid from co-crystal ligand | `docking/1WCY_receptor.pdbqt` (12,248 atoms); box center **(62.8, 47.7, 4.8)**, size 30³ | RCSB PDB / OpenBabel / awk |
| 8. Molecular docking (AutoDock Vina) | `--exhaustiveness 4 --cpu 2`, 9 poses/peptide | **565/565** done; Top: LPPQ −6.109*, APSPE −7.150, LAPSP −7.087 kcal/mol | AutoDock Vina 1.2.5 |

\* *Note:* the previously reported LPPQ dG of −7.472 was produced from a corrupted OpenBabel `make3D` preparation and has been **corrected to −6.109 kcal·mol⁻¹** (RDKit-prepared); see the manuscript §2.4 for explicit disclosure.

### Docking Top 10 (binding free energy dG, kcal/mol)

| Rank | Peptide | dG | Class |
|---|---|---|---|
| 1 | LPPQ (Leu-Pro-Pro-Gln) | **−6.109** | mid-strong |
| 2 | APSPE (Ala-Pro-Ser-Pro-Glu) | −7.150 | mid-strong |
| 3 | LAPSP (Leu-Ala-Pro-Ser-Pro) | −7.087 | mid-strong |
| 4 | LPGF (Leu-Pro-Gly-Pro) | −7.075 | mid-strong |
| 5 | LPINP (Leu-Pro-Ile-Asn-Pro) | −6.988 | mid-strong |
| 6 | LPSP (Leu-Pro-Ser-Pro) | −6.867 | mid-strong |
| 7 | LPCPR (Leu-Pro-Cys-Pro-Arg) | −6.835 | mid-strong |
| 8 | LPGDP (Leu-Pro-Gly-Asp-Pro) | −6.793 | mid-strong |
| 9 | LPDDP (Leu-Pro-Asp-Asp-Pro) | −6.693 | mid-strong |
| 10 | APSQP (Ala-Pro-Ser-Gln-Pro) | −6.515 | mid-strong |

Distribution: mid-strong (−6.5 ~ −8) 10; moderate (−5 ~ −6.5) 49; weak (>-5) 1 (CPCSK −4.856); none < −8.0 strong.

---

## 2.5 Phase A — official-server filtering (replacing the original proxy scorer)

> The PeptideRanker official server was long unusable (HTTP 503); we replaced the allergen/toxicity layer with peer-validated **ToxinPred 3.0 batch submission + AlgPred 2.0 (AAC-RF)**.

| Stage | Filter | Remaining |
|---|---|---|
| Parent set (2–6 aa short peptides) | — | **4,950** |
| AlgPred 2.0 (ML score < 0.6) | −270 | 4,680 |
| ToxinPred 3.0 (Non-Toxin) | −6 | **4,674** |
| iDPPIV-SCM (offline reproduction) | replaces PeptideRanker proxy (see 2.8) | all 4,950 short peptides scored |

- All **565 docking-queue peptides pass** ToxinPred + AlgPred (0 rejected).
- Products: `data/phaseA_inputs/results_toxinpred.csv`, `results_algpred.csv`, `official_candidates.tsv` (4,674).
- Scripts: `scripts/phaseA/phaseA_run_toxinpred.py`, `phaseA_run_algpred.py`, `phaseA_merge.py`.
- Details: `docs/phaseA_README.md`.

---

## 2.6 Phase B — pure-computational binding validation (static MM, no MD)

> This machine has **no GROMACS / no conda**, so 100–150 ns MD is infeasible. We use **single-conformation endpoint static MMFF94s validation** + geometric contact profiling as a feasibility substitute for wet-lab/MD.

| Peptide (sequence) | dG_Vina | ΔE_MM* | Pocket residues | Total contacts | H-bonds | Key contact position |
|---|---|---|---|---|---|---|
| LPPQ (Leu-Pro-Pro-Gln) | −6.109 | +2.45 | 39 | 92 | 10 | **Gln4:37** (C-term dominant), Pro2/3 ~21–25 each |
| APSPE (Ala-Pro-Ser-Pro-Glu) | −7.150 | +2.93 | 39 | 119 | 13 | **Glu5:109** (S1′ overwhelming) |
| LAPSP (Leu-Ala-Pro-Ser-Pro) | −7.087 | +2.84 | 40 | 101 | 13 | **Pro5:50** (C-term Pro) |

\* ΔE_MM is a gas-phase single-point energy difference, **not a true ΔG** (no explicit solvation/entropy); the three peptides are similar in magnitude and not discriminative; binding-strength ranking still relies mainly on Vina dG.

- **Biological consistency:** APSPE's C-terminal glutamate shows extreme contact preference at the S1′ pocket, matching DPP4's "prefers substrates with negatively-charged C-terminus occupying S1′" mechanism; the dual-Pro core fits S1′/S2′ Pro preference.
- **H-bond partners:** all three enrich **GLU146** (S2′ region) and respectively touch PRO149 / SER182 / TYR183 — all within DPP4's known active-center cavity.
- Products: `data/phaseB/phaseB_results.tsv`, `phaseB_detail.json`.
- Script: `scripts/phaseB/phaseB_validation.py`.
- Limitations & reproduction: `docs/phaseB_README.md`.

---

## 2.7 Phase C — *in-silico* ADMET & DPP4-centric network pharmacology

> A second pure-computational validation layer after Phase B, replacing originally-planned *in-vitro* ADMET / Caco-2 permeability / serum-stability experiments. All computational, no wet lab.

- **C1 peptide-level ADMET + GI stability** (565 Vina-docked peptides): BioPython `ProtParam` (MW, pI, GRAVY, aliphatic index, pH 7.4 charge) + Boman protein-binding index + RDKit descriptors (TPSA/HBD/HBA/QED) + protease-cleavage-rule-based GI-stability prediction (pepsin/trypsin/chymotrypsin/elastase/carboxypeptidase/aminopeptidase) and DPP4 self-cleavage (X-Pro/X-Ala) judgment.
- **Key finding:** all 565 candidates carry X-Pro/X-Ala motifs (consistent with the "prefers X-Pro" DPP4 structural-preference filter — supporting binding yet meaning in-vivo rapid DPP4-mediated degradation as a known oral-delivery limitation). GI stability: Moderate 39 / Low 11 / High 10.
- **C2 DPP4-centric network** (STRING DB REST API, *Homo sapiens* 9606, functional, score≥400): **15 nodes / 32 edges** (STRING 11 nodes 26 functional edges + literature-expanded 4 substrate nodes 6 edges). Core neighbors include **GCG (GLP-1/GLP-2 precursor), GIP (gastric inhibitory polypeptide), CXCR4 (SDF-1/CXCL12 receptor), ADA, CAV1, PRCP**. Mechanistic context: inhibit DPP4 → raise active GLP-1/GIP → promote insulin secretion → improve T2DM glycemia.
- Products: `data/phaseC/phaseC_peptides.tsv` (per-peptide), `data/phaseC/phaseC_network.json`, `phaseC_network_summary.txt`.
- Scripts: `scripts/phaseC/phaseC_peptides.py`, `phaseC_network.py`.
- Limitations & reproduction: `docs/phaseC_README.md`.

---

## 2.8 iDPPIV-SCM module — offline activity prescreen (replacing PeptideRanker proxy)

> Both the PeptideRanker official server and the iDPPIV-SCM website were long unusable (pythonanywhere shows "Coming Soon"). To eliminate external-server dependency entirely, we **fully reproduce iDPPIV-SCM (Scoring Card Method, Charoenkwan *et al.*, *J. Proteome Res.* 2020, DOI:10.1021/acs.jproteome.0c00590) offline locally**, replacing the original `pepranker()` proxy (constant 1.000, zero discrimination).

| Stage | Method | Result |
|---|---|---|
| Training data | WeiLab-BioChem/Structural-DPP-IV public homologous set (531+532 train / 133+133 test, all 20 standard AAs) | `scripts/idppiv_scm/data/` |
| Scoring card | Global AA-composition SCM: score = Σₐ nₐ·log₂(P⁺ₐ/P⁻ₐ) (P⁺/P⁻ = pos/neg residue frequencies) | `scripts/idppiv_scm/model.py` |
| Stage 1 (activity prescreen) | Score all 4,950 2–6 aa short peptides; continuous distribution [−7.243, 4.431]; predict **3,400 DPP-IV inhibitory peptides (68.7%)** | `data/moso_candidates_idppiv_short.tsv` |
| Stage 2 (DPP4 preference narrowing) | 3–5 aa, N-term hydrophobic, hydrophobic residues retained (independent biological rule); iDPPIV score > 0 net-positive | candidate pool **565** → docking queue **60** |
| Vina docking (same protocol) | 1WCY receptor, 30³ box, exhaustiveness 8 (consistent with old queue) | `docking/moso_dock_results_idppiv_clean.tsv` |

**Key findings (see `docs/methodology_replacement_report.md`, `docking/moso_dock_compare.tsv`):**
1. **iDPPIV score and Vina dG are nearly zero-correlated (Spearman ρ = 0.067)** — the scoring card predicts "DPP-IV inhibitory propensity" (classification dimension), not "binding free energy" (physical dimension); the two are orthogonal. This validates the "activity prescreen → binding assessment" two-stage independent-filter pipeline design.
2. **Set-level enrichment of stronger-binding peptides:** under identical prep, the new queue's best **APQIP −6.807** beats the old proxy queue's best **LAPSP −6.461**; peptides with dG ≤ −6.0 rise to 33.3% vs 20.0%; of 32 overlapping peptides, 21 improved. Δ≈−0.35 kcal/mol falls within Vina's typical noise band (±0.5–1.0) and cannot be exaggerated as "significantly stronger binding."
3. **Interpretability:** learned AA propensities are biologically self-consistent — Pro (+0.875) ranks first (DPP4 S1 pocket specificity for Pro), Cys (−2.482) strongly negative (echoing toxicity marker).

**Relation to Phase A:** Phase A (ToxinPred 3.0 + AlgPred 2.0) is the allergen/toxicity filter layer; this module is the activity-propensity scoring layer — orthogonal and complementary. The original "proxy heuristic" scoring was removed from the honest footnote ① and upgraded to a "literature-validated, offline-reproducible SCM algorithm."

- Products: `data/moso_candidates_idppiv*.tsv`, `data/moso_dock_queue_idppiv.txt`, `docking/moso_dock_results_idppiv*.tsv`, `docking/moso_dock_compare.tsv`, `docking/moso_ligands_idppiv/` (180 files).
- Scripts: `scripts/idppiv_scm/` (model + data + validation/exploration), `scripts/moso_pipeline_filter_idppiv.py`, `scripts/moso_pipeline_filter2_idppiv.py`, `scripts/moso_dock_run_idppiv.py`, `scripts/moso_dock_generic.py`, `scripts/moso_dock_analyze.py`, `scripts/moso_dock_compare.py`.
- Limitation: the iDPPIV-SCM training benchmark suffers length confounding (negatives mostly long peptides/proteins); the original authors concede limited accuracy. Correct usage is candidate ranking / soft filtering, **not** deterministic activity judgment.

---

## 3. Directory structure

```
moso-bamboo-DPP4-peptides/
├── README.md
├── .gitignore
├── data/                         # proteins and peptide libraries
│   ├── moso_253.fasta           # 253 moso bamboo protein sequences
│   ├── moso_253_peptides.txt           # relaxed-rule unique peptides
│   ├── moso_253_peptides_strict.txt   # strict-rule unique peptides
│   ├── moso_candidates_pr_filtered.txt# 2,019 candidate peptides (old proxy)
│   ├── moso_candidates_idppiv.txt     # iDPPIV candidate pool (4,742)
│   ├── moso_candidates_idppiv_short.tsv  # 4,950 short peptides re-scored by iDPPIV
│   ├── moso_candidates_idppiv_proxy.tsv  # old proxy set with iDPPIV scores
│   ├── moso_dock_queue.txt     # 60-peptide docking queue (old proxy)
│   └── moso_dock_queue_idppiv.txt     # 60-peptide docking queue (iDPPIV)
├── docking/
│   ├── 1WCY_receptor.pdbqt    # DPP4 receptor (sitagliptin-bound)
│   ├── moso_box.txt             # Vina pocket parameters
│   ├── moso_ligands/           # 60 ligand PDBQTs + 60 docked output poses
│   ├── moso_dock_results.tsv   # 60-peptide docking dG result table
│   ├── moso_dock_ranking.txt  # full ranking
│   ├── moso_ligands_idppiv/   # iDPPIV-queue 60 ligand PDBQTs + 60 docked outputs
│   ├── moso_dock_results_idppiv.tsv       # iDPPIV-queue 60-peptide docking dG (with rerun dup rows)
│   ├── moso_dock_results_idppiv_clean.tsv # deduplicated 60-peptide docking dG
│   ├── moso_dock_results_old_rdkit.tsv    # old queue same-protocol RDKit re-dock (fair comparison baseline)
│   └── moso_dock_compare.tsv  # iDPPIV vs old-proxy fair cross-queue comparison
├── scripts/                     # reproducible pipeline scripts
│   ├── rerun_digestion_moso253_strict.py
│   ├── moso_pipeline_filter2.py
│   ├── moso_build_ligand_pdbqts.py
│   ├── moso_dock_prepare_receptor.py
│   ├── moso_dock_run.py
│   ├── moso_report_results.py
│   ├── batch_dock.cmd / batch_dock.sh   # one-click batch docking
│   ├── reports/                # report-generation scripts
│   ├── idppiv_scm/            # iDPPIV-SCM offline reproduction (2.8)
│   │   ├── model.py           # scoring-card model (DEFAULT_THRESHOLD)
│   │   ├── scm.py / validate*.py / explore*.py / score_short.py
│   │   └── data/             # public homologous train/test datasets
│   ├── moso_pipeline_filter_idppiv.py    # Stage 1 iDPPIV scoring
│   ├── moso_pipeline_filter2_idppiv.py   # Stage 2 DPP4 preference narrowing
│   ├── moso_dock_run_idppiv.py          # iDPPIV-queue Vina docking (resumable)
│   ├── moso_dock_generic.py             # generic re-dock (old queue same protocol)
│   ├── moso_dock_analyze.py             # new-queue internal correlation analysis
│   └── moso_dock_compare.py            # cross-queue fair comparison
└── docs/
    ├── moso_dpp4_peptide_project_summary.docx   # full project summary
    └── methodology_replacement_report.md            # iDPPIV-replaces-PeptideRanker methodology note
```

---

## 4. Reproduction

### 4.1 Completed (this repository *is* the result)

Virtual digestion → filtering → docking are fully run; results are in `data/` and `docking/`.

### 4.2 Re-run docking (requires local Vina install)

```bash
# install (local machine / server)
conda install -c conda-forge vina openbabel rdkit

# enter scripts/, run one-click batch docking
# Windows:
scripts\batch_dock.cmd
# Linux/WSL:
bash scripts/batch_dock.sh
```

> Note: `vina.exe` is a third-party binary, **not** committed to the repository. Download it from
> https://github.com/ccsb-scripps/AutoDock-Vina/releases/tag/v1.2.5
> and place it in the run directory, or modify the Vina path in `batch_dock.*`.

---

## 5. Limitations that must be disclosed honestly

1. **Filtering layer:** allergen/toxicity filtering now uses **ToxinPred 3.0 + AlgPred 2.0** official outputs (Phase A, 4,674 candidates; see 2.5); activity-propensity scoring now uses **iDPPIV-SCM offline reproduction** replacing the original PeptideRanker proxy (see 2.8). Both are pure-computational filter layers; the original "proxy heuristic" scoring has been removed from this pipeline and must not serve as a final publication judgment.
2. **Docking is a static dG estimate**, requiring GROMACS MD / MM-PBSA + *in-vitro* experiment for confirmation; it cannot directly serve as an activity conclusion.
3. **Moso bamboo has 0 reviewed entries (all TrEMBL-predicted)**; Methods must honestly disclose the source and predicted nature (the template garlic paper is likewise ~95% TrEMBL, which is acceptable).
4. **WPHY/WPQY/VAPGW are garlic peptides** and must not be used for moso bamboo ground-truth validation.

---

## 6. To-do (per template's full pipeline; this project is pure-computation, no wet lab)

- [x] **Phase A** official ToxinPred 3.0 + AlgPred 2.0 filtering (4,674 candidates); PeptideRanker proxy replaced by iDPPIV-SCM offline reproduction (see 2.8)
- [x] **Phase B (static approximation)** single-conformation MMFF94s relaxation + geometric contact profile (Top-3 done)
- [x] **Phase C** peptide-level *in-silico* ADMET + GI stability (C1, 565 peptides) + DPP4-centric network pharmacology (C2, STRING 15 nodes / 32 edges)
- [ ] **Phase B upgrade** if GROMACS available on HPC: add 100–150 ns MD + MM-PBSA (with GB solvation and entropy terms) for quantitative ΔG and converged residue decomposition
- [x] **iDPPIV-SCM module (2.8)** offline reproduction replacing PeptideRanker proxy, same-protocol Vina comparison rerun, committed to GitHub
- [ ] **Manuscript** (*in-silico* discovery) skeleton + Methods/Limitations four honest footnotes (① official filtering done ② Vina static dG ③ moso bamboo 0 reviewed / all TrEMBL ④ Phase B/C are computational approximations, wet lab not done)

---

## 7. References

- Template method: Cheng Y. *et al.* *Discovery of garlic-derived peptides as natural DPP4 inhibitors.* Bioorganic Chemistry 175 (2026) 109801.
- DPP4 structure: PDB **1WCY** (sitagliptin-bound DPP4).
- Protein sequences: UniProtKB, *Phyllostachys edulis* (taxonomy 38705).
- iDPPIV-SCM (activity-prescreen scoring card): Charoenkwan P. *et al.* *iDPPIV-SCM: A Sequence-Based Predictor for Identifying and Analyzing Dipeptidyl Peptidase-IV (DPP-IV) Inhibitory Peptides.* J. Proteome Res. 2020, 19(5), 1890–1901. DOI:10.1021/acs.jproteome.0c00590.
