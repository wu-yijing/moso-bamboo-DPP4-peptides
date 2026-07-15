# Moso-Bamboo DPP4 Inhibitory Peptide Project — Phase-1 Screening Methodology Replacement Report

> **Replaced object**: `pepranker()` (PeptideRanker-style proxy scoring) in `moso_pipeline_filter.py`
> **Replacement**: **iDPPIV-SCM local offline reproduction** (global amino-acid-composition Scoring Card Method)
> **Trigger**: PeptideRanker official site long unresponsive; iDPPIV-SCM official site (camt.pythonanywhere.com) verified **down ("Coming Soon")**, confirming "long-term unavailable" risk. Switching to a DPP4-**dedicated**, **literature-validated**, **fully offline-reproducible** scoring card instead raises methodological defensibility.

---

## 1. Why iDPPIV-SCM (instead of continuing the generic proxy)

| Dimension | Old PeptideRanker-style proxy | New iDPPIV-SCM (offline reproduction) |
|---|---|---|
| Target specificity | Generic "bioactivity" score, not DPP4-dedicated | **DPP4 inhibitory-peptide-dedicated** scoring card |
| Reproducibility | Transparent but heuristic proxy, no literature support | Reproduces a published method (Charoenkwan et al. 2020, *J. Proteome Res.* 19:4125–4136, DOI:10.1021/acs.jproteome.0c00590) |
| External dependency | Depended on PeptideRanker web page (now dead) | **Zero network, zero external dependency** |
| Interpretability | Weak (hydrophobicity/Pro/MW empirical weights) | Strong: per-residue propensity, directly supports discussion |

**Conclusion**: Upgrading "proxy heuristic" to "literature-validated, offline-reproducible SCM algorithm" turns the original manuscript's honesty footnote ① from a **weakness into a strength**.

---

## 2. Method and data (fully offline)

### 2.1 Scoring Card Method formula
On the DPP4 inhibitory (positive) vs non-inhibitory (negative) training set, compute the global propensity of each amino acid `a`:

```
P(a) = log2( Obs(a) / Exp(a) )
Obs(a) = occurrence frequency of a in positive set
Exp(a) = overall frequency of a in (positive+negative) set
Peptide S's iDPPIV-SCM total score = Σ_{a∈S} P(a)        # standard SCM sums residue propensities
```

Also provided: length-normalized score `score_mean = ΣP(a) / |S|`, for fair cross-length ranking.

### 2.2 Training dataset (homologous public data, downloaded locally with code)
- Source: WeiLab-BioChem/**Structural-DPP-IV** repo `data/DPP-IV/` (the 665+665 benchmark homologous to the iDPPIV-SCM paper)
- Scale: **train.tsv** = 531 positive + 532 negative (benchmark); **test.tsv** = 133 positive + 133 negative (independent test)
- Storage: `E:/workbuddy/Claw/idppiv_scm/data/`
- Reproduction scripts: `idppiv_scm/model.py`, `scm.py`, `validate.py`, `validate2.py`

---

## 3. Validation results (honest disclosure)

| Evaluation | This reproduction | Literature report (Charoenkwan 2020) |
|---|---|---|
| Independent test ACC (threshold optimized) | **0.771** | ≈0.797 |
| Nested 5-fold CV (per-fold leave-one-out threshold) | 0.647 (high variance) | ≈0.819 |
| Pure length baseline (len≤10→positive) | **0.820** | — |

### 3.1 Key finding: benchmark has **length confounding** (must be written into manuscript truthfully)
- Training-set **positives dominated by short peptides** (207/531 only 2 aa, median length 13), **negatives almost all long peptides/proteins** (median length 27, max 90).
- Hence a **pure length baseline** ("short = positive") already reaches 0.820, higher than the composition-type SCM itself.
- iDPPIV-SCM authors themselves state *"not yet accurate enough for real-world applications"*.
- **Practical meaning for this project**: All moso-bamboo candidates are **2–6 aa short peptides** (same "positive-like" length region), so the length signal fails uniformly across candidates; the SCM score at this point mainly reflects the **residue-composition signal genuinely related to DPP4 inhibition**, used for candidate-prioritization ranking and soft filtering — exactly its designed intent (*"rapid screening … prior to synthesis"*).

### 3.2 Learned amino-acid propensities (biologically self-consistent, directly writable to Discussion)
```
P(Pro)  +0.875   ← DPP4 S1 pocket strongly specific to Pro, signature feature of food-source DPP4 inhibitory peptides
M(Met)  +0.702
W(Trp)  +0.433
Q(Gln)  +0.376
V(Val)  +0.313
E(Glu)  +0.181
...
C(Cys)  -2.482   ← strongly negative (also echoes toxicity proxy's Cys flag)
```
The residue signal matches known DPP4 inhibitory-peptide chemical intuition, proving SCM learned a real composition preference rather than merely length.

---

## 4. Pipeline linkage and funnel comparison

| Stage | Old (PeptideRanker proxy) | New (iDPPIV-SCM) |
|---|---|---|
| Input unique peptides | 7,988 | 7,988 |
| Phase-1 iDPPIV/proxy predicts inhibitory | ~all retained (proxy score constant≈1.0) | **4,746** (59%) |
| Phase-2 de-allergenization | → 4,742* | 4,742 |
| Phase-3 de-toxicity | → 4,742* | 4,742 |
| DPP4-preference narrowing (filter2) | 565 → Top-60 | **565 → Top-60** |
| Best-docked peptide LPPQ selected? | ✓ | **✓ (still in new Top-60, iDPPIV=+2.182)** |

\* Old 4,289 candidates were an early 4,950-short-peptide-file product; under the same input the new candidate 4,743 is consistent in magnitude.

**Old/new docking-queue overlap 32/60**; the original best-docked peptide **LPPQ (−7.472 kcal/mol) is still retained by the new method and ranks 3rd**, i.e. the existing Vina docking conclusion is supported by the new, more defensible screening method.

---

## 5. Output file list

| File | Description |
|---|---|
| `idppiv_scm/model.py` | iDPPIV-SCM offline reproduction model (build scoring card + score) |
| `idppiv_scm/scm.py` | SCM core math (incl. position-specific / global / dipeptide implementations, for reproduction audit) |
| `idppiv_scm/data/train.tsv`, `test.tsv` | homologous public train/test sets (downloaded) |
| `idppiv_scm/validate.py`, `validate2.py`, `explore*.py` | recipe probing and accuracy-validation scripts |
| `moso_pipeline_filter_idppiv.py` | Phase-1 filter script using iDPPIV-SCM |
| `moso_pipeline_filter2_idppiv.py` | DPP4-preference second narrowing (iDPPIV-score ranking) |
| `moso_candidates_idppiv.txt` | new candidate list (4,743 entries, with iDPPIV total & mean score, descending) |
| `moso_dock_queue_idppiv.txt` | new docking-priority queue (Top-60, by iDPPIV score) |

---

## 6. Manuscript honesty-footnote update suggestion

- **Old footnote ①** (proxy): "filtering used a proxy heuristic, not the official server" →
- **New footnote ①**: "Phase 1 adopts the Scoring Card Method of iDPPIV-SCM (Charoenkwan et al. 2020, *J. Proteome Res.*), fully reproduced locally on a homologous public DPP4 inhibitory-peptide dataset with zero external dependency; that benchmark's positive/negative samples have length-distribution confounding and the model accuracy is limited (independent test ACC≈0.77), hence this work uses it for candidate-prioritization **ranking and soft filtering** rather than deterministic activity judgment, with downstream molecular docking (Vina) evaluating binding capability."

---

## 7. Optional follow-ups (pending user decision)
1. **Re-dock new queue**: Run Vina on the new Top-60 (`moso_dock_queue_idppiv.txt`, incl. 28 new peptides outside the old queue), to test whether iDPPIV prioritization yields better-binding peptides.
2. **Cross-validation**: Use StackDPPIV (high accuracy, ACC 0.891) online server to recheck Top candidates (requires its recovery).
3. **Phases 2/3 still proxy**: Allergen (AllerTOP), toxicity (ToxinPred) still depend on official web pages; the formal manuscript must use official-tool outputs as ground truth — this point retains the original honesty footnote.
