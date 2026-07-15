# Moso-Bamboo Protein-Derived DPP4 Inhibitory Peptides (Pure-Computation Study): Literature Support and Improvement Plan

> Search date: 2026-07-15 | Sources: PubMed / PMC / Nature (*Scientific Reports*) / Web of Science-indexed journals (confirmed via journal official sites and indexes) | Language: mostly Chinese; proper nouns / numbers retained as in source.

## 0. Search notes and honesty statement

- This search covers PubMed (incl. PMC open full text), Nature's *Scientific Reports*, and journals in the Web of Science core collection (several cited in WoS, e.g. iDPPIV-SCM cited 106 times). **This machine has no direct access to the Web of Science subscription database**; all listed WoS-journal literature was confirmed for bibliographic entry and DOI via the publisher's official site / index page, with citation status subject to the journal's WoS coverage.
- Literature is tagged as "support" or "improve" for direct mapping onto manuscript sections.
- All values (IC₅₀, ACC, ΔG, etc.) are taken from retrieved abstracts / text; before formal submission, the original PDFs must be double-checked for page numbers and exact values.

---

## 1. Core literature list (grouped by topic)

### Group A · Computational methodology (our Phase-1 screener and its upgraded version)

| ID | Citation | Key conclusion | Tag |
|---|---|---|---|
| **A1** | Charoenkwan P, Kanthawong S, Nantasenamat C, Hasan MM, Shoombuatong W. **iDPPIV-SCM: A Sequence-Based Predictor for Identifying and Analyzing Dipeptidyl Peptidase IV (DPP-IV) Inhibitory Peptides Using a Scoring Card Method.** *J Proteome Res.* 2020, 19(10): 4125–4136. DOI: 10.1021/acs.jproteome.0c00590 | PMID: 32897718 | First sequence-type DPP4 inhibitory-peptide predictor (Scoring Card Method); independent test ACC≈0.797, better than kNN/LR/DT, comparable to SVM; **authors themselves state the benchmark has length confounding and limited model accuracy**; public web server (camt.pythonanywhere.com, now "Coming Soon" offline). | **Support** (direct methodological source) |
| **A2** | Zou H. **iDPPIV-SI: identifying dipeptidyl peptidase IV inhibitory peptides by using multiple sequence information.** *J Biomol Struct Dyn.* 2024, 42(4): 2144–2152. DOI: 10.1080/07391102.2023.2203257 | PMID: 37125813 | Next-gen predictor: 50 physicochemical properties + first/second-order correlations + discrete wavelet transform + LASSO feature selection + SVM; **train ACC 91.26%, independent ACC 98.12%**, significantly better than prior predictors; dataset and MATLAB code open-sourced (figshare). | **Improve** (Phase-1 upgrade / cross-validation candidate) |

### Group B · DPP4 substrate specificity and structural mechanism (support contact-fingerprint / MM-GBSA interpretation)

| ID | Citation | Key conclusion | Tag |
|---|---|---|---|
| **B1** | Hirabayashi K, et al. **Structural Basis of Proline-Specific Exopeptidase Activity as Observed in Human Dipeptidyl Peptidase-IV.** *Structure.* 2003, 11(8): 947–956. DOI: 10.1016/S0969-2126(03)00160-6 | Human DPP-IV co-crystal with Diprotin A (Ile-Pro-Ile): reveals **S1 pocket** (Tyr662/Tyr666/Val656/Val711/Trp659/Tyr631), **dual-Glu motif Glu205/Glu206** (anchors substrate N-terminus), catalytic triad Ser630/Asp708/His740; establishes "Pro at position 2 required to enter S1 pocket." | **Support** (catalytic-core anchor + Pro@P1 mechanism) |
| **B2** | Engel M, et al. **Contribution of amino acids in the active site of dipeptidyl peptidase 4 to the catalytic action of the enzyme.** *PLoS One.* 2023. DOI: 10.1371/journal.pone.0289239 | Systematically validates active-center residues (Glu205/Glu206/Asn710/Arg125) contribution to hydrolysis vs inhibition binding; reaffirms S1-specific pocket, P1 Pro/Ala preference, and **length dependence**. | **Support** (short-peptide fit, P1=Pro) |
| **B3** | The Medical Biochemistry Page. **Dipeptidyl Peptidase 4 (DPP4).** (review entry) | Textbook confirmation: DPP4 only cleaves N-terminal "X-Pro / X-Ala" dipeptides, S1/S2/S3 pocket division, monomer inactive requiring homodimer. | **Support** (background and boundary) |

### Group C · Food-derived (esp. plant-source) DPP4 inhibitory peptides (support topic choice + provide benchmark active peptides)

| ID | Citation | Key conclusion | Tag |
|---|---|---|---|
| **C1** | Xie Peng, Yuan Shaofei, Zhang Jian, Huang Qi, Wang Hongyan, He Liang. **Virtual screening and activity of DPP-IV inhibitory peptides from moso bamboo shoot proteins.** *World Bamboo and Rattan.* 2026, 24(1): 25–32. DOI: 10.12168/sjzttx.2025.11.25.001 | **Most directly comparable study**: 8 moso-bamboo protein sequences → PeptideCutter 4-enzyme virtual hydrolysis → 2597 peptides → PeptideRanker score>0.5 yields 270 2–5-mers → ADMET (AdmetSAR/Innovagen/ToxinPred/AllerTOP) screen → 6 → CDOCKER docking selects EGF/RY/LR → **in-vitro EGF tripeptide IC₅₀ = 366.44 ± 2.15 μg/mL**. | **Support + comparison anchor** |
| **C2** | **Identification and Characterization of Dipeptidyl Peptidase-IV Inhibitory Peptides from Oat Proteins.** (PMC9141920, 2022) | Oat-source 10 new peptides; **IPQHY / VPQHY / VAVVPF / VPLGGF IC₅₀ < 50 μM** (mixed-inhibition type), docked to active-site residues. | **Support** (plant short peptides truly potent) |
| **C3** | **Identification of Novel Dipeptidyl Peptidase-IV Inhibitory Peptides in Chickpea Protein Hydrolysates.** *J Agric Food Chem.* 2023, 71(21): 8211–8219. DOI: 10.1021/acs.jafc.3c00603 | Chickpea-source; molecular docking validates AAWPGHPEF/LAFP/IAIPPGIPYW/PPGIPYW binding active center; **IAIPPGIPYW IC₅₀ = 12.43 μM** (contains repeated Pro motif), effective in Caco-2 cells. | **Support** (Pro-enriched peptide potent + cellular validation paradigm) |
| **C4** | **Deciphering the molecular mechanism of DPP-IV inhibition by quinoa-derived peptides…** *Food Chemistry.* 2026. DOI: 10.1016/j.foodchem.2026.013270 (in-text S0308814626013270) | Quinoa pentapeptides **WLAFR / LLPFR IC₅₀ = 78 / 130 μM**; kinetics + spectroscopy + in silico, hydrophobicity-dominated; MD validates conformational stability. | **Support** (pentapeptide activity + docking/MD paradigm) |
| **C5** | **Screening and evaluation of novel DPP-IV inhibitory peptides in goat milk…** *ScienceDirect* 2025 (S259015752500063X) | Explicitly states: **"Vina score correlates only partially with DPP4 inhibition, molecular docking cannot distinguish substrates from inhibitors"**; lists food-derived peptide docking-energy range (−5.5~−10.8 kcal/mol). | **Support** (directly supports our "Vina/contact-count ≠ activity" argument) |
| **C6** | **Virtual Hydrolysis-Based Screening of Wheat-Derived DPP-IV Inhibitory Peptides…** *J Agric Food Chem.* 2025. DOI: 10.1021/acs.jafc.5c03006 | PMID: 40623964 | ConPLex deep learning + virtual hydrolysis + cellular experiment + MD: 4 peptides IC₅₀ 2.89–4.96 mM (competitive), tau-RaMD computes binding residence time. | **Support + wet-experiment route reference** |
| **C7** | **Bee pollen-derived peptide with dual DPP-IV inhibition and glucose transport modulation.** *Nat Sci Rep.* 2026. DOI: 10.1038/s41598-026-39009-1 | Bee-pollen peptide ATHALLA (IC₅₀ 52.63 μM) + molecular docking + Caco-2 + in-silico ADMET; natural-protein-source functional-peptide mining paradigm. | **Support** (recent Nature direct analogue) |
| **C8** | Nongonierma AB, FitzGerald RJ series reviews (*Food Protein Hydrolysates as Source of DPP-IV Inhibitory Peptides*, Proc Nutr Soc 2013, etc.) | Summarizes food-derived DPP4 inhibitory peptides (milk IPI/Ile-Pro-Ile positive control, soy IAVPTGVA 223 μM, lupin LTFPGSAED 228 μM, etc.) IC₅₀ table. | **Support** (expand literature-comparison table) |

### Group D · MM-GBSA reliability and protocols (support §2.6/§3.5 methodology and honest disclosure)

| ID | Citation | Key conclusion | Tag |
|---|---|---|---|
| **D1** | **Comparison of affinity ranking using AutoDock-GPU and MM-GBSA scores for BACE-1 inhibitors in the D3R Grand Challenge 4.** *J Comput Aided Mol Des.* 2020 (PMC7027993) | **MM-GBSA rescoring did NOT improve correlation with experiment**; results depend on initial conformation, protonation, ligand charge. | **Support** (honest disclosure + defense of "why still use MM-GBSA") |
| **D2** | Karaman B, Sippl W. **Docking and binding free energy calculations of sirtuin inhibitors.** *Eur J Med Chem.* 2015. DOI: 10.1016/j.ejmech.2015.01.014 | MM-GBSA significantly correlates with activity data, better than docking score; used as "post-docking filter." igb=1 performs best, igb=8 also good. | **Support** (rationale for MM-GBSA as ranking/filter) |
| **D3** | **Improving predictive performance of MM/PBSA(GBSA) for protein–cyclic peptide complexes.** *Brief Bioinform.* 2025 | Systematically optimizes force field / dielectric constant, correlation Rp improved to −0.732 (vs traditional +131.6%). | **Improve** (MM-GBSA protocol tuning direction) |
| **D4** | **GPU-Accelerated Virtual Screening and MD Simulations for Identification of Novel DPP-4 Inhibitors.** *ACS Omega.* 2025. DOI: 10.1021/acsomega.5c08231 | On DPP-4, MM-GBSA ΔG range −58 ~ −5.8 kcal/mol; **includes negative control acetaminophen (−8.50) validating discrimination**. | **Support + Improve** (DPP-4 system usability + negative-control practice) |

---

## 2. Section-by-section support mapping

| Manuscript section | Citable literature | Use |
|---|---|---|
| Introduction (topic rationale) | C1, C2, C7, C8 | Plant/natural protein sources truly release DPP4 inhibitory peptides; moso bamboo already has experimental validation (C1) → well-founded topic |
| Phase 1 (iDPPIV-SCM) | A1 (method source), A2 (upgrade), D1 (limitation) | Method origin + honest disclosure "length confounding / limited accuracy" directly quoted from A1 authors |
| Docking protocol benchmark (①) | B1, D4, C5 | Known inhibitor / active-pocket residue controls; negative-control idea from D4 |
| §3.5 MM-GBSA | D1, D2, D3, D4 | Single-structure approximation, relative ranking, protocol-tuning basis |
| §3.6/§3.7 Contact profile and joint interpretation | B1, B2, C5 | "Glu146 anchors catalytic core" "P1=Pro fits S1" "docking cannot distinguish substrate/inhibitor" |
| §3.8 Contact fingerprint | B1, C3 | Structural corroboration of Pro-enriched peptide strong anchoring |
| §4 Discussion / limitations | D1, A1, C5 | Forward citations of MM-GBSA limitation, SCM limitation, docking limitation |
| Conclusion / outlook (wet-experiment route) | C1, C6, C3, C4 | "computational screening → synthesis → IC₅₀ → Caco-2 → MD mechanism" full paradigm |

---

## 3. Key arguments directly writable into the paper (with literature anchors)

1. **"Moso bamboo (bamboo shoot) endogenous proteins truly release DPP4 inhibitory peptides, and this is already in-vitro validated"** — cite C1 (EGF IC₅₀ 366 μg/mL) and patent CN118702773A (bamboo peptide IC₅₀ 305–541 μM). This provides **experimental-level circumstantial evidence** for our "moso-bamboo full-proteome as source" pure-computation topic (note: C1 uses PeptideRanker+CDOCKER, complementary not duplicative with our iDPPIV-SCM+Vina+MM-GBSA).
2. **"DPP4 only recognizes substrates with Pro/Ala at N-terminal position 2; the S1 hydrophobic pocket is the primary specificity determinant"** — cite B1/B2. Our three candidates (LPPGP / APPSQ / APQIP) **all have Pro at position 2 (P1=Pro)**, which should be explicitly stated in Results/Discussion as "structural fit," corroborating the Pro propensity +0.875 in iDPPIV-SCM, forming strong activity-supporting evidence.
3. **"More contacts ≠ higher-quality binding; docking cannot distinguish substrate from inhibitor"** — cite C5 (original statement) + B1 (Diprotin A is actually a slow-hydrolysis substrate, not a strong inhibitor). This is the theoretical basis for our using **MM-GBSA + contact profile (anchoring catalytic core Glu146 rather than glycosylation bias)** to distinguish "true inhibitor vs false positive," and should be front-cited in §4.2 to strengthen methodology.
4. **"MM-GBSA is best as relative ranking / post-docking filter, not an absolute-affinity claim"** — cite D2 (supportive usage) + D1 (BACE-1 rescoring no better than docking, hinting limitation). Our "Vina + MM-GBSA + iDPPIV-SCM three-signal joint, divergence = information" strategy itself aligns with literature's robust understanding.
5. **"Short peptides (2–6 aa) are the mainstream length of food-derived DPP4 inhibitory peptides"** — cite C2/C3/C4/C8 (oat, chickpea, quinoa, milk peptides are all 2–6 mers). Consistent with our candidate length range, supporting the "short-peptide-priority" screening preference.

---

## 4. Improvement plan (by priority)

### Plan 1 (high priority · easy) · Introduce iDPPIV-SI cross-validation in Phase 1
- **Basis**: A2 independent ACC 98.12%, dataset/code open-source (figshare).
- **Action**: Use A2's SVM+LASSO model to rescore our 7,988 candidates, intersect/union with iDPPIV-SCM; report two-model agreement (Kappa/overlap rate) in the manuscript to compensate SCM's length-confounding weakness.
- **Output**: new `data/phaseA/idppiv_si_crosscheck.tsv` + a methodology section.

### Plan 2 (high priority · easy) · Expand literature-comparison table (incl. IC₅₀ benchmarks)
- **Basis**: C2/C3/C4/C8 provide experimentally validated food-derived peptide IC₅₀.
- **Action**: Append ≥10 experimentally-IC₅₀ peptides to `data/literature_dpp4_peptides.tsv` (IPI positive control, oat IPQHY, chickpea IAIPPGIPYW, quinoa WLAFR/LLPFR, soy IAVPTGVA, lupin LTFPGSAED, wheat 4 peptides, bee ATHALLA, bamboo EGF), letting reviewers compare "our candidates vs validated active peptides" horizontally.
- **Note**: Pure-computation candidates have no IC₅₀; clearly distinguish in-table as "in silico only / to-be-validated."

### Plan 3 (high priority · argument strengthening) · Explicitly discuss P1=Pro mechanism fit
- **Basis**: B1/B2 + our candidates all have Pro at position 2.
- **Action**: Add 1 paragraph in §3 Results or §4 Discussion linking "sequence Pro preference (SCM propensity) + structural S1-pocket Pro specificity + candidate P1=Pro" into the **core structural argument supporting activity**.

### Plan 4 (medium priority · methodology uplift) · MM-GBSA protocol tuning + negative control
- **Basis**: D3 (force-field/dielectric tuning significantly improves), D4 (negative control acetaminophen validates discrimination).
- **Action**: For pocket-constrained MM-GBSA, try GB model `igb=1` or `igb=8` (per D2), scan internal dielectric εin ∈ {1.0, 2.0, 4.0}; add 1–2 **negative-control pentapeptides** (random inactive sequence / known weak-binding peptide) to validate ranking discrimination; multiple sampling (e.g. 3× independent single-structure minimization averaged) to reduce single-structure bias.
- **Output**: methodology subsection + sensitivity appendix table.

### Plan 5 (medium priority · forward-looking) · Introduce MD sampling / binding residence time
- **Basis**: C6 (wheat, tau-RaMD residence time), C4 (quinoa, 100 ns MD validation).
- **Action**: Run 50–100 ns MD on Top-2 (LPPGP/APPSQ) (local compute needs assessment), replace single-structure approximation with trajectory-averaged MM-GBSA, report binding free energy RMSD/RMSF stability over time.
- **Risk**: This environment kills foreground processes at 2 min and reclaims background tasks; must split into short segments or offline batch; **feasibility to be assessed**, listed as future work.

### Plan 6 (low priority · boundary defense) · Forward-cite MM-GBSA limitation to defend against reviewer doubt
- **Basis**: D1 (BACE-1 rescoring no better than docking).
- **Action**: In §4.2 state clearly: "We are aware MM-GBSA rescoring may not beat docking (D1), hence we position it as a relative-ranking metric within the three-signal joint rather than an independent claim, and interpret APQIP's divergence accordingly."

### Plan 7 (low priority · benchmarking) · Method-comparison table with moso-bamboo homologue study (C1)
- **Basis**: C1 is the most direct prior (same organism, same target) but different method (PeptideRanker+CDOCKER vs iDPPIV-SCM+Vina+MM-GBSA+contact fingerprint).
- **Action**: Add a small section/table comparing the two workflows' inputs (8 hand-picked proteins vs 253 UniProt full proteome), screeners, docking tools, validation modes (C1 has wet experiment vs our pure computation), highlighting our incremental contribution of "omics-scale + multi-signal joint + contact fingerprint."

### Plan 8 (forward-looking · route) · Provide an actionable wet-experiment validation route
- **Basis**: C1 (EGF IC₅₀), C6 (wheat full paradigm), C3 (Caco-2 cells).
- **Action**: In Conclusion/Outlook state: "If wet-experiment conditions are available, we recommend per C6 paradigm: chemically synthesize LPPGP/APPSQ → DPP4 enzyme-inhibition IC₅₀ assay → Caco-2 uptake/stability → MD mechanism resolution," forming a benchmarking anchor with C1's moso-bamboo EGF result.

---

## 5. Reference summary (for pasting into manuscript References)

1. Charoenkwan P, Kanthawong S, Nantasenamat C, Hasan MM, Shoombuatong W. iDPPIV-SCM: A Sequence-Based Predictor for Identifying and Analyzing Dipeptidyl Peptidase IV (DPP-IV) Inhibitory Peptides Using a Scoring Card Method. *J Proteome Res.* 2020;19(10):4125–4136. doi:10.1021/acs.jproteome.0c00590. PMID:32897718.
2. Zou H. iDPPIV-SI: identifying dipeptidyl peptidase IV inhibitory peptides by using multiple sequence information. *J Biomol Struct Dyn.* 2024;42(4):2144–2152. doi:10.1080/07391102.2023.2203257. PMID:37125813.
3. Hirabayashi K, et al. Structural Basis of Proline-Specific Exopeptidase Activity as Observed in Human Dipeptidyl Peptidase-IV. *Structure.* 2003;11(8):947–956. doi:10.1016/S0969-2126(03)00160-6.
4. Engel M, et al. Contribution of amino acids in the active site of dipeptidyl peptidase 4 to the catalytic action of the enzyme. *PLoS One.* 2023. doi:10.1371/journal.pone.0289239.
5. Xie Peng, Yuan Shaofei, Zhang Jian, Huang Qi, Wang Hongyan, He Liang. Virtual screening and activity of DPP-IV inhibitory peptides from moso bamboo shoot proteins. *World Bamboo and Rattan.* 2026;24(1):25–32. doi:10.12168/sjzttx.2025.11.25.001.
6. Identification and Characterization of Dipeptidyl Peptidase-IV Inhibitory Peptides from Oat Proteins. (PMC9141920, 2022).
7. Identification of Novel Dipeptidyl Peptidase-IV Inhibitory Peptides in Chickpea Protein Hydrolysates. *J Agric Food Chem.* 2023;71(21):8211–8219. doi:10.1021/acs.jafc.3c00603.
8. Deciphering the molecular mechanism of DPP-IV inhibition by quinoa-derived peptides. *Food Chemistry.* 2026. doi:10.1016/j.foodchem.2026.013270.
9. Screening and evaluation of novel DPP-IV inhibitory peptides in goat milk based on molecular docking and molecular dynamics simulation. *ScienceDirect* 2025 (S259015752500063X).
10. Virtual Hydrolysis-Based Screening of Wheat-Derived DPP-IV Inhibitory Peptides. *J Agric Food Chem.* 2025. doi:10.1021/acs.jafc.5c03006. PMID:40623964.
11. Bee pollen-derived peptide with dual DPP-IV inhibition and glucose transport modulation. *Nat Sci Rep.* 2026. doi:10.1038/s41598-026-39009-1.
12. Comparison of affinity ranking using AutoDock-GPU and MM-GBSA scores for BACE-1 inhibitors in the D3R Grand Challenge 4. *J Comput Aided Mol Des.* 2020 (PMC7027993).
13. Karaman B, Sippl W. Docking and binding free energy calculations of sirtuin inhibitors. *Eur J Med Chem.* 2015. doi:10.1016/j.ejmech.2015.01.014.
14. Improving predictive performance of MM/PBSA(GBSA) for protein–cyclic peptide complexes. *Brief Bioinform.* 2025.
15. GPU-Accelerated Virtual Screening and MD Simulations for Identification of Novel DPP-4 Inhibitors. *ACS Omega.* 2025. doi:10.1021/acsomega.5c08231.
16. Nongonierma AB, FitzGerald RJ. Food protein hydrolysates as a source of dipeptidyl peptidase IV inhibitory peptides. *Proc Nutr Soc.* 2013;72(1):1–12.

---

## 6. One-sentence overview

> Literature consistently shows: **plant/natural protein sources truly release experimentally-active DPP4 inhibitory short peptides** (especially characterized by P1=Pro), and **docking score ≠ activity, MM-GBSA best as relative ranking** — this both strongly supports our pure-computation topic and the "three-signal joint + contact fingerprint" methodology, and reveals 8 immediately-actionable improvements (highest priority: iDPPIV-SI cross-validation, expanded IC₅₀ literature table, explicit P1=Pro mechanism-fit discussion). The most direct prior is Xie Peng et al. (2026) moso-bamboo homologue study (already in-vitro validated EGF IC₅₀≈366 μg/mL), usable as a method-comparison and benchmarking anchor.

---

## 7. Implementation status (2026-07-15)

- **Plan ⑦ (method-comparison table with moso-bamboo homologue study) — IMPLEMENTED.** Written into manuscript `docs/manuscript_draft.md` **§3.9** (Table 4: eight-dimension methodological comparison, highlighting this study's omics-scale + dual-independent-predictor + contact-fingerprint increment).
- **Plan ① (iDPPIV-SI cross-validation) — IMPLEMENTED (methodologically-consistent offline reproduction).**
  - Script `scripts/phaseA/idppiv_si_crosscheck.py`: reproduces iDPPIV-SI's LASSO-SVM paradigm on the homologous 665+665 benchmark (composition features 467-dim → LASSO selects 114 → RBF-SVM).
  - Results `data/phaseA/idppiv_si_crosscheck.tsv` / `_summary.txt`, written into manuscript **§3.10**.
  - Key values: benchmark test ACC=0.669 (lower than SCM ~0.797, reflecting simplified composition features); the three finalists independently locate **LPPGP (si_prob 0.892, rank 34/4742) consistent top, APQIP (si_prob 0.249, rank 2738/4742) independently downgraded (corroborating MM-GBSA/contact-fingerprint downgrade conclusion), APPSQ neutral (0.563)**; two-model continuous-score Spearman ρ=0.290.
  - **Honest disclosure**: Zou H (2024) original-model MATLAB weights are on figshare; this environment has no MATLAB, so it is a methodologically-consistent Python reproduction rather than byte-level reimplementation; A2 citation (Zou H, *J Biomol Struct Dyn* 42(4):2144–2152, DOI 10.1080/07391102.2023.2203257, PMID 37125813) was verified via search to **truly exist**.
- **Plan ② (expand IC₅₀ literature table) — IMPLEMENTED** (previous round): `data/literature_dpp4_peptides.tsv` appended with 14 experimentally-IC₅₀ food-derived peptides.
- **Plan ③ (P1=Pro mechanism fit) — PARTIALLY implemented**: §3.9/§4.1 already state candidates all have Pro at position 2, fitting S1-pocket specificity; can be made an explicit paragraph in §3.2 or §4 (to-be-added).
- **Plan ④⑤⑧ (MM-GBSA tuning / MD sampling / wet-experiment route) — NOT implemented**: Plans ④⑤ need extra compute (this environment kills foreground at 2 min, reclaims background, listed as future work); Plan ⑧ route already given in §4.3.
