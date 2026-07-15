# -*- coding: utf-8 -*-
"""Internal analysis of the new iDPPIV queue: de-dup + iDPPIV score vs measured dG correlation (free of preparation-method confound)"""
import collections, statistics, json, os

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # scripts/ -> repo root
DOCK = os.path.join(REPO, "docking")

RAW = os.path.join(DOCK, "moso_dock_results_idppiv.tsv")
# aggregate all docking values per peptide (multiple re-runs)
agg = collections.defaultdict(list)
for l in open(RAW, encoding="utf-8"):
    l = l.rstrip("\n")
    if not l or l.startswith("peptide"):
        continue
    seq, sc, dg = l.split("\t")
    if dg in ("NA", "ERR"):
        continue
    agg[seq].append((float(sc), float(dg)))

# clean 60-peptide table: take each peptide's best (most negative) dG
clean = []
for seq, vals in agg.items():
    sc = vals[0][0]
    dgs = [v[1] for v in vals]
    best = min(dgs)
    spread = max(dgs) - min(dgs)
    clean.append((seq, sc, best, spread, len(vals)))

clean.sort(key=lambda x: x[2])   # ascending by dG (most negative first)
print(f"unique peptides: {len(clean)}  peptides with >1 re-run: {sum(1 for c in clean if c[4]>1)}")
print("\n=== new iDPPIV queue Top-15 (best dG) ===")
print(f"{'rank':<5}{'peptide':<8}{'iDPPIV':>8}{'dG_best':>9}{'spread':>8}{'n_dock':>7}")
for i, (seq, sc, best, spread, n) in enumerate(clean[:15], 1):
    print(f"{i:<5}{seq:<8}{sc:>8.3f}{best:>9.3f}{spread:>8.3f}{n:>7}")

# Spearman correlation between iDPPIV ranking score and dG
def spearman(xs, ys):
    n = len(xs)
    rx = sorted(range(n), key=lambda i: xs[i])
    ry = sorted(range(n), key=lambda i: ys[i])
    rankx = [0]*n; ranky = [0]*n
    for r, i in enumerate(rx): rankx[i] = r
    for r, i in enumerate(ry): ranky[i] = r
    d2 = sum((rankx[i]-ranky[i])**2 for i in range(n))
    return 1 - 6*d2/(n*(n**2-1))

scs = [c[1] for c in clean]
dgs = [c[2] for c in clean]
rho = spearman(scs, dgs)
print(f"\n=== correlation (n={len(clean)}) ===")
print(f"iDPPIV_score vs dG  Spearman rho = {rho:.3f}")
# bin: dG median of iDPPIV high/low groups
hi = [d for s, d in zip(scs, dgs) if s >= statistics.median(scs)]
lo = [d for s, d in zip(scs, dgs) if s < statistics.median(scs)]
print(f"iDPPIV high-group dG median: {statistics.median(hi):.3f}  (n={len(hi)})")
print(f"iDPPIV low-group  dG median: {statistics.median(lo):.3f}  (n={len(lo)})")

# save clean table
_clean_out = os.path.join(DOCK, "moso_dock_results_idppiv_clean.tsv")
with open(_clean_out, "w", encoding="utf-8") as f:
    f.write("peptide\tiDPPIV_score\tdG_best\tdG_spread\tn_dock\n")
    for seq, sc, best, spread, n in clean:
        f.write(f"{seq}\t{sc:.3f}\t{best:.3f}\t{spread:.3f}\t{n}\n")
print(f"\nclean table written -> {_clean_out}")
