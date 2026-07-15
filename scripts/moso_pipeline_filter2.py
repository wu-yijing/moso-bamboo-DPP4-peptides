# -*- coding: utf-8 -*-
"""Moso-bamboo DPP4 candidate peptides -- second narrowing (known structural
preferences of DPP4 inhibitory peptides)
Based on literature: N-terminal hydrophobic (I/L/V/F/W/M); length 2-5 aa optimal;
position-2 Pro/Ala preferred; hydrophobic-rich.
On top of the 4289 PeptideRanker>0.5 and non-allergenic / non-toxic candidates,
layer the DPP4 preference rules -> docking-scale set.
"""
import os
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # scripts/ -> repo root
DATA = os.path.join(REPO, "data")

HYDRO = set("ILVFMWACY")          # hydrophobic
def dpp4_rules(seq, score):
    if len(seq) < 3 or len(seq) > 5: return False   # 3-5aa is typical DPP4 inhibitory-peptide length
    if seq[0] not in HYDRO: return False          # N-term hydrophobic
    if score < 0.55: return False                # slightly raised threshold
    if not any(a in HYDRO for a in seq): return False
    return True

rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(DATA, "moso_candidates_pr_filtered.txt")) if l.strip()]
cands = [(p, float(s)) for p, s in rows]
narrow = [(p, s) for p, s in cands if dpp4_rules(p, s)]
print(f"4289 candidates -> DPP4 preference narrowed -> {len(narrow)} entries")

# secondary priority: position-2 Pro/Ala marker
def p2(p): return len(p) >= 2 and p[1] in "PA"
tier1 = [(p, s) for p, s in narrow if p2(p)]      # position-2 P/A (optimal)
tier2 = [(p, s) for p, s in narrow if not p2(p)]
print(f"  of which position-2 P/A (optimal): {len(tier1)} | rest: {len(tier2)}")

# output: docking-priority set (tier1 first)
n_top = min(60, len(narrow))
top = (sorted(tier1, key=lambda x:-x[1]) + sorted(tier2, key=lambda x:-x[1]))[:n_top]
_queue = os.path.join(DATA, "moso_dock_queue.txt")
with open(_queue, "w") as f:
    for p, s in top:
        f.write(f"{p}\t{s:.3f}\t{'P2' if p2(p) else ''}\n")
print(f"\ndocking queue (top {n_top}) -> {_queue}")
print("samples (first 15):")
for p, s in top[:15]:
    print(f"  {p:6s} score={s:.3f} {'[P2]' if p2(p) else ''}")

print("\n=== full funnel (template garlic 1442->249->34 synthesized) ===")
print(f"moso bamboo: 7988 -> 4333 -> 4333 -> 4289 -> {len(narrow)} (DPP4 preference) -> docking queue {len(top)}")
