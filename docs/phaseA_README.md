# Phase A — Official-server filtering replacing the proxy scorer

## Purpose

Previously, the PeptideRanker-style scores in `moso_candidates_pr_filtered.txt` were a **proxy heuristic** (all values 1.000, zero discrimination) and could not serve as a screening basis for a computational paper. This phase uses **official online servers** to run real predictions on the **4,950 2–6 aa short peptides** obtained from moso bamboo virtual digestion, producing a publishable "official-filter" candidate set that replaces the proxy results.

## Inputs (in `data/phaseA_inputs/`)

| File | Description |
|---|---|
| `moso_short_2to6.fasta` | Main input, all 4,950, with `>pep_00001`-style IDs |
| `moso_short_2to6.txt` | Plain list, one per line |
| `peptideranker/batch_XX.fasta` | 500 per batch, 10 batches (manual submission) |

## Three-stage filter funnel

```
parent (4,950)
  └─ [optional] PeptideRanker ≥ 0.5 — server temporarily unusable
  └─ AlgPred 2.0 ML Score < 0.6 — ✅ auto-submit
  └─ ToxinPred3 Non-Toxin         — ✅ auto-submit
  └─ official_candidates.tsv
```

## How to run

### 1️⃣ ToxinPred (auto — ~5–10 min)
```bash
python scripts/phaseA/phaseA_run_toxinpred.py
```
- Auto-submits 4,950 peptides to the ToxinPred batch server in batches of **50**.
- Parses the HTML result table, extracting SVM Score + Prediction.
- Output → `data/phaseA_inputs/results_toxinpred.csv`

### 2️⃣ AlgPred 2.0 (auto — ~5–10 min)
```bash
python scripts/phaseA/phaseA_run_algpred.py
```
- Auto-submits to the AlgPred 2.0 Batch server (AAC-RF mode).
- Parses the result table, extracting ML Score + Prediction.
- Output → `data/phaseA_inputs/results_algpred.csv`

### 3️⃣ PeptideRanker (manual — server temporarily unusable)
- URL: `https://peptide.ucd.ie/peptideranker/` 👈 currently returns **503 Service Unavailable**
- After recovery: open the URL, paste `data/phaseA_inputs/moso_short_2to6.fasta`
  (or each `data/phaseA_inputs/peptideranker/batch_XX.fasta`)
- Save results separately as `data/phaseA_inputs/results_peptideranker.tsv`

### 4️⃣ Merge all results
```bash
python scripts/phaseA/phaseA_merge.py
```
The script prints the **funnel counts** and `official_candidates.tsv`. It auto-skips any missing result file.

## Server details

| Server | URL | Returns | Pass condition |
|---|---|---|---|
| **PeptideRanker** (UCD) | https://peptide.ucd.ie/peptideranker/ ⚠ 503 | 0–1 bioactivity probability | ≥ 0.5 |
| **AlgPred 2.0** | https://webs.iiitd.edu.in/raghava/algpred2/batch.html | ML Score + Allergen/Non-allergen | ML Score < 0.6 |
| **ToxinPred3 Batch** | http://www.raghavagps.net/raghava/toxinpred/multi_submit.php | SVM Score + Toxin/Non-Toxin | Non-Toxin |

## Config parameters (top CONFIG block of `scripts/phaseA/phaseA_merge.py`)
- `THR_ALGPRED = 0.6` — AlgPred ML Score threshold (below = non-allergen)
- `TOXIN_NEG_LABELS` — ToxinPred negative-label set
- `THR_PR = 0.5` — PeptideRanker threshold (enable after server recovery)

## Honesty statement (to write into Methods)
- ToxinPred and AlgPred predictions were completed via official-server HTTP batch submission;
- PeptideRanker requires manual submission (server currently unusable);
- The moso bamboo proteome consists of UniProt TrEMBL predicted entries (0 Swiss-Prot reviewed), disclosed honestly in Methods.
