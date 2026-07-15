# -*- coding: utf-8 -*-
"""
Re-score ALL 4,950 moso-bamboo DPP4 short peptides (2-6 aa) with the locally
reproduced iDPPIV-SCM, replacing the PeptideRanker-style proxy scores
(original moso_candidates_pr_filtered.txt was all 1.000).

Output:
  moso_candidates_idppiv_short.tsv   -- iDPPIV-SCM scores for all 4,950 short peptides
  moso_candidates_idppiv_proxy.tsv   -- the old proxy candidate set (4,289) with iDPPIV scores
"""
import os, sys
_HERE = os.path.dirname(os.path.abspath(__file__))          # scripts/idppiv_scm
_SCRIPTS = os.path.dirname(_HERE)                            # scripts/
REPO = os.path.dirname(_SCRIPTS)                             # repo root
sys.path.insert(0, _SCRIPTS)                                 # make `import idppiv_scm.model` work
from idppiv_scm.model import build_propensity, score, score_mean, predict, DEFAULT_THRESHOLD

AA = set("ACDEFGHIKLMNPQRSTVWY")
HERE = os.path.join(REPO, "data")                           # input/output data dir

# 1) train the scoring card
P, stats = build_propensity()
Tp, Tn, npos, nneg = stats
print(f"[scoring card] train pos={npos} neg={nneg}; threshold(predict inhibitory)= {DEFAULT_THRESHOLD:+.3f}")

def load_short(path):
    out = []
    for l in open(path, encoding="utf-8"):
        s = l.strip()
        if 2 <= len(s) <= 6 and all(c in AA for c in s):
            out.append(s)
    return out

def score_file(in_path, out_path, label):
    seqs = load_short(in_path) if "strict" in in_path else [
        l.strip() for l in open(in_path, encoding="utf-8") if l.strip()
    ]
    rows = []
    for s in seqs:
        sc = score(s, P)
        sm = score_mean(s, P)
        pr = predict(s, P)
        rows.append((s, len(s), sc, sm, pr))
    rows.sort(key=lambda r: -r[2])
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("peptide\tlength\tiDPPIV_score\tiDPPIV_mean\tpredicted_DPP4_inhibitory\n")
        for s, ln, sc, sm, pr in rows:
            f.write(f"{s}\t{ln}\t{sc:.4f}\t{sm:.4f}\t{pr}\n")
    npos = sum(1 for r in rows if r[4] == 1)
    print(f"[{label}] total {len(rows)} entries -> predicted DPP-IV inhibitory peptide: {npos} ({100*npos/len(rows):.1f}%)")
    print(f"        iDPPIV_score range: [{min(r[2] for r in rows):.3f}, {max(r[2] for r in rows):.3f}]")
    print(f"        file: {out_path}")
    return rows

# 2a) all 4,950 short peptides
all_rows = score_file(
    os.path.join(HERE, "moso_253_peptides_strict.txt"),
    os.path.join(HERE, "moso_candidates_idppiv_short.tsv"),
    "all short peptides (2-6 aa)",
)

# 2b) old proxy candidate set (4,289) -> iDPPIV scores
proxy_path = os.path.join(HERE, "moso_candidates_pr_filtered.txt")
if os.path.exists(proxy_path):
    score_file(
        proxy_path,
        os.path.join(HERE, "moso_candidates_idppiv_proxy.tsv"),
        "old proxy candidate set (orig. PeptideRanker)",
    )

# 3) show Top-20
print("\n=== Top-20 (by iDPPIV_score descending) ===")
print(f"{'peptide':<10}{'len':>4}{'iDPPIV_score':>15}{'iDPPIV_mean':>14}{'pred':>6}")
for s, ln, sc, sm, pr in all_rows[:20]:
    print(f"{s:<10}{ln:>4}{sc:>15.4f}{sm:>14.4f}{pr:>6}")

# 4) compare the old proxy scores ("all 1.000") with the true iDPPIV distribution
import statistics
vals = [r[2] for r in all_rows]
print(f"\n=== distribution comparison ===")
print(f"old PeptideRanker proxy score: constant = 1.000 (no discriminative power)")
print(f"new iDPPIV-SCM score:      mean={statistics.mean(vals):.3f}  median={statistics.median(vals):.3f}  stdev={statistics.pstdev(vals):.3f}")
