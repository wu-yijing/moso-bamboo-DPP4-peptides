# -*- coding: utf-8 -*-
"""
Fair cross-queue comparison: new iDPPIV queue vs old proxy queue
(both docked against 1WCY with the same RDKit pipeline)
------------------------------------------------------------------------
Read:
  moso_dock_results_idppiv_clean.tsv  (new queue 60, de-dup best kept)
  moso_dock_results_old_rdkit.tsv    (old queue 60, re-docked same way, de-dup best kept)
Output: console comparison + moso_dock_compare.tsv
"""
import collections, statistics, os

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # scripts/ -> repo root
DOCK = os.path.join(REPO, "docking")

def load_best(path, has_score=True):
    agg = collections.defaultdict(list)
    for l in open(path, encoding="utf-8"):
        l = l.rstrip("\n")
        if not l or l.startswith("peptide"):
            continue
        parts = l.split("\t")
        seq = parts[0]
        sc = float(parts[1]) if len(parts) > 1 else 0.0
        if len(parts) < 3 or parts[2] in ("NA", "ERR"):
            continue
        agg[seq].append((sc, float(parts[2])))
    out = {}
    for seq, vals in agg.items():
        out[seq] = (vals[0][0], min(v[1] for v in vals))  # (score, best_dG)
    return out

new = load_best(os.path.join(DOCK, "moso_dock_results_idppiv_clean.tsv"))
old = load_best(os.path.join(DOCK, "moso_dock_results_old_rdkit.tsv"))

print(f"new queue unique peptides: {len(new)}   old queue unique peptides: {len(old)}")

new_dg = [v[1] for v in new.values()]
old_dg = [v[1] for v in old.values()]

def stats(name, xs):
    xs_sorted = sorted(xs)
    print(f"\n[{name}] n={len(xs)}")
    print(f"  best dG (most negative): {min(xs):.3f}")
    print(f"  median dG:        {statistics.median(xs):.3f}")
    print(f"  mean dG:        {statistics.mean(xs):.3f}")
    print(f"  frac dG <= -6.5: {sum(1 for x in xs if x<=-6.5)/len(xs)*100:.1f}%")
    print(f"  frac dG <= -6.0: {sum(1 for x in xs if x<=-6.0)/len(xs)*100:.1f}%")

stats("new iDPPIV queue (same RDKit method)", new_dg)
stats("old proxy queue (same RDKit method)", old_dg)

# directly pair-compare overlapping peptides (shared by both queues)
overlap = set(new) & set(old)
print(f"\n=== paired comparison of overlapping peptides (n={len(overlap)}) ===")
better_new = sum(1 for p in overlap if new[p][1] < old[p][1])   # more negative = better
better_old = sum(1 for p in overlap if old[p][1] < new[p][1])
print(f"  peptides with better (more negative) dG in new queue: {better_new}")
print(f"  peptides with better (more negative) dG in old queue: {better_old}")

new_best = min(new_dg); old_best = min(old_dg)
print(f"\n=== conclusive comparison ===")
print(f"  new-queue best-binding peptide dG = {new_best:.3f} ({min(new, key=lambda p: new[p][1])})")
print(f"  old-queue best-binding peptide dG = {old_best:.3f} ({min(old, key=lambda p: old[p][1])})")
if new_best < old_best:
    print(f"  >>> under the same preparation, the new iDPPIV prioritization yields a better-binding peptide (Δ={new_best-old_best:+.3f} kcal/mol)")
else:
    print(f"  >>> under the same preparation, the new queue's best does not surpass the old queue (Δ={new_best-old_best:+.3f} kcal/mol)")
    print(f"      i.e. iDPPIV prioritization brought no stronger-binding peptide on this receptor/box; its value lies in the activity-screening dimension (see analysis).")

# save
with open(os.path.join(DOCK, "moso_dock_compare.tsv"), "w", encoding="utf-8") as f:
    f.write("metric\tnew_idppiv\told_proxy\n")
    f.write(f"n\t{len(new)}\t{len(old)}\n")
    f.write(f"best_dG\t{new_best:.3f}\t{old_best:.3f}\n")
    f.write(f"median_dG\t{statistics.median(new_dg):.3f}\t{statistics.median(old_dg):.3f}\n")
    f.write(f"mean_dG\t{statistics.mean(new_dg):.3f}\t{statistics.mean(old_dg):.3f}\n")
    f.write(f"frac<=-6.5\t{sum(1 for x in new_dg if x<=-6.5)/len(new_dg):.3f}\t{sum(1 for x in old_dg if x<=-6.5)/len(old_dg):.3f}\n")
    f.write(f"overlap_n\t{len(overlap)}\t\n")
    f.write(f"overlap_better_new\t{better_new}\t\n")
    f.write(f"overlap_better_old\t{better_old}\t\n")
print(f"\ncomparison written -> {os.path.join(DOCK, 'moso_dock_compare.tsv')}")
