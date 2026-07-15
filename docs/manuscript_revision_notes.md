# Manuscript Revision Notes (review of the reconstructed `manuscript_draft.md`)

> Review date: 2026-07-15. Overall assessment: the three-pillar structure (① length bias ② multi-target ③ genome/transcriptome) is complete, the honesty boundaries are in place, and the differentiation from Xie Peng 2026 is clear. The items below are graded by **severity**: Class A are hard errors that must be fixed before submission; Class B are structure and numbering; Class C are typos/factual issues; Class D are submission-readiness items.

---

## A. Critical consistency errors (must fix before submission)

### A1. [Most severe] The DPP4-axis definition in the multi-target framework is self-contradictory
- **Symptom**: In §3.11 Table 8, the DPP4 row shows `corr(score, length) = +0.155` (positive correlation: longer → higher score); but §3.9 explicitly states that the SCM score ~ length ρ = **−0.188** (negative correlation: shorter → higher score), and §3.11 line 254 itself says "opposite in direction to the §3.9 DPP4 benchmark … SCM score vs length is **negatively correlated**".
- **Root cause**: The "DPP4 normalised score" in the multi-target matrix (Table 9: LPPGP 0.438 / APQIP 0.330 / APPSQ 0.242) and the "iDPPIV-SCM probability 0.490 / rank 77/4742" in §3.13 **cannot be the same score** — if LPPGP ranked 77/4742 in SCM (top 1.6%, normalised score ≈0.98), its "DPP4 normalised score" of only 0.438 (43.8th percentile) is impossible. This shows the multi-target DPP4 axis is actually the **"reference-set signature method" homologous to ACE/αG (cosine + 3-mer Jaccard)**, not the "benchmark-validated predictor iDPPIV-SCM" claimed in §2.9's table.
- **Fix (recommended option A, minimal change and more honest)**:
  1. In §2.9's three-axis table, change the DPP4-axis source to "**reference-set signature method (homologous to ACE/αG, not the SCM predictor)**", and delete the "benchmark-validated predictor (reliable)" label;
  2. In §3.11 state explicitly: "the three axes of the multi-target framework **uniformly adopt the reference-set signature method** for comparability; the independently validated iDPPIV-SCM predictor score is reported separately in §3.9/§3.13";
  3. Keep the Table 8 DPP4 `corr=+0.155` (self-consistent with the signature method's positive correlation), but change §3.11 line 254 from "opposite in direction to the §3.9 negative correlation" to "all three axes are **positively correlated** under the signature method, whereas the §3.9 SCM predictor itself is **negatively correlated** — the two DPP4 scores exhibit different length behaviours, which exactly shows the predictor and the signature method capture different signals".
  - This change actually **strengthens the honesty narrative**: ACE/αG have no validated predictor, so a single method across all three axes is the most transparent; SCM is listed separately as the validated DPP4-specific score.

### A2. §3.7 "iDPPIV double-first" conflicts with the absolute ranking in §3.13
- **Symptom**: §3.7 table/line 156 calls LPPGP "(MM-GBSA, iDPPIV) **double-first**"; but §3.13 shows LPPGP's iDPPIV-SCM **rank 77 / 4,742** (top 1.6%, not "first").
- **Fix**: change "double-first" to "**highest iDPPIV probability among the three finalist peptides** (0.490)", or "**top overall priority**", to avoid colliding with the absolute rank 77.

### A3. Candidate-pool count 4,742 vs 4,717 unexplained
- **Symptom**: §3.1 / §3.13 state the SCM-positive pool is 4,742 peptides; §3.11's multi-target matrix analyses **4,717** (a gap of 25).
- **Fix**: add a sentence explaining (expected: 25 peptides containing non-standard characters or with all three-axis 3-mer overlaps equal to 0 were removed), or confirm the removal rule in the script and state it.

### A4. APPSQ "hit-target count 0" yet called "DPP4-focused"
- **Symptom**: In Table 9, APPSQ has `hit-target count = 0` (MULTI tag none), but §3.11 line 271 calls "APPSQ / APQIP **DPP4-focused** candidates". Under the framework's own rule "top 10% = a hit", APPSQ did not enter the top 10% on the DPP4 signature axis (0.242, the lowest of the three finalists), so 0 hits is self-consistent, but the "DPP4-focused" wording contradicts its own 0 hits.
- **Fix**: restrict "DPP4-focused" to "**by the SCM and docking priority in §3.7/§3.9, APPSQ ranks second among DPP4 single-target candidates**", and clearly separate this from the multi-target signature's 0 hits in two distinct sentences.

---

## B. Structure and numbering

### B1. Table numbering jumps abruptly from "Table 5"
- The inline tables in §3.2–§3.7 are **unnumbered**; §3.9 suddenly refers to "Table 5", with Tables 1–4 missing in between.
- **Fix**: number all tables sequentially (suggest: §3.2 candidate table = Table 1, §3.4 benchmark = Table 2, §3.5 MM-GBSA = Table 3, §3.7 priority master table = Table 4, then continue), and update the in-text references accordingly.

### B2. Figure numbering has Fig 3/Fig 4 with no Fig 1/Fig 2 before them
- §3.8 references Fig 3 and Fig 4 (contact fingerprints), but §3.1–§3.7 contain no figures.
- **Fix**: either turn the screening funnel (§3.1), the priority master table (§3.7), or the docking pose into Fig 1/Fig 2 to make the numbering continuous; or relabel the contact-fingerprint figures as Fig 1/Fig 2 and state in the preamble "this paper contains only two figures".

### B3. Inconsistent figure paths
- §3.8 body text writes `../figures/contact_fingerprint_LPPGP_APPSQ.svg`, while the data index (line 358) writes `figures/...` (missing `../`). Unify to one relative path.

---

## C. Typos and facts

### C1. Journal name spelling "Gigascience" → "GigaScience"
- Appears in §3.10 line 210 and reference line 345. Correct is **GigaScience** (Peng et al. 2018, GigaScience).

### C2. Table 9 tag "MULTI" typo
- Line 267 `MULTI（DPP4+ACE）` is a typo for **MULTI** (multi-target).

### C3. §2.8 length-baseline thresholds T=6/8/10/12 but only T=10 reported
- §2.8 line 62 lists 4 thresholds, but §3.9 and the abstract report only len≤10. State "primary threshold T=10", with the others as sensitivity scans (add a table if present; otherwise keep only T=10).

### C4. Benchmark-set size off-by-one: 1,063 vs 1,062
- §3.9: train 1,063 / test 266; §3.13: train 531/531, test 133/133 (total 1,328). The homologous set used by SCM and SI reproduction differs by 1 entry. Verify the actual row counts in `scripts/idppiv_scm/data/` and the SI script, and unify to a single number.

### C5. Methods describe in present tense docking already done in another environment
- §2.4–§2.6 describe DPP4 docking and MM-GBSA in present tense ("adopt / recompute with …"), but this environment no longer has rdkit/openbabel/vina (see the idea② exploration log).
- **Fix**: add a sentence "the above DPP4 docking and MM-GBSA calculations were completed in a separate computing environment; all input structures, receptor PDBQT, scripts, and outputs are archived with the manuscript to ensure reproducibility", to avoid reviewers questioning "why are the tools described in Methods currently unavailable".

### C6. §1.1 example peptides uncited
- "Soy-derived DYPAY/IAVPTGVA, camel-milk-derived LPVP/MPVQA" etc. should be cited (can be merged into the `reference` column of `data/literature_dpp4_peptides.tsv`).

---

## D. Submission readiness

### D1. References "to be added"
- Line 342 "References (to be added)". **This is a hard blocker**: before submission you must complete the formal entries with DOIs (Charoenkwan 2020, Zou 2024, Peng 2018/2013, Xie Peng 2026, 1WCY PDB, original food-derived peptide literature). At minimum, complete the 6 entries already named in the text.

### D2. Language version vs target journal
- Currently a Chinese manuscript (mixed Chinese-English). If submitting to a bioinformatics venue per the established positioning (*BMC Bioinformatics / PLOS Computational Biology / Scientific Reports* etc.), **a complete English version is required**; if submitting to a Chinese journal (e.g. the parent system of *World Bamboo and Rattan*), it may be retained. It is recommended to state the target journal and language version in the status line at the top of the manuscript.

### D3. Missing Data Availability statement
- The manuscript repeatedly says "scripts released with the manuscript", but there is no formal Data Availability / Code Availability section. It is recommended to add a sentence pointing to this GitHub repository (`git@github.com:wu-yijing/moso-bamboo-DPP4-peptides.git`) and the source of the benchmark sets used.

---

## Priority suggestions
1. Fix **A1 (multi-target DPP4-axis definition)** first — this is the core of methodological self-consistency and is what reviewers will catch;
2. Fix **A2/A3/A4** in parallel (numeric and wording conflicts);
3. Then fix **B1–B3** (numbering makes the manuscript look "finished" rather than a draft);
4. Finally complete **C/D** (typos + references + data statement).
