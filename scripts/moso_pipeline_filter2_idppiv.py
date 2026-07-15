# -*- coding: utf-8 -*-
"""
Moso-bamboo DPP4 candidate peptides -- second narrowing (known structural
preferences of DPP4 inhibitory peptides) -- iDPPIV-SCM version
==========================================================================
Based on literature: N-terminal hydrophobic (I/L/V/F/W/M); length 2-5 aa
(this rule uses 3-5) is optimal; hydrophobic-rich;
        DPP4 S1 pocket specificity for Pro (position-2 Pro/Ala is optimal).
On top of the iDPPIV-SCM candidates (which already carry DPP4-inhibition propensity),
layer the DPP4 preference rules -> docking-scale set.
"""
import os
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # scripts/ -> repo root
DATA = os.path.join(REPO, "data")

HYDRO = set("ILVFMWACY")          # hydrophobic

def dpp4_rules(seq, score):
    if len(seq) < 3 or len(seq) > 5: return False   # 3-5aa is typical DPP4 inhibitory-peptide length
    if seq[0] not in HYDRO: return False              # N-term hydrophobic
    if score <= 0.0: return False                   # iDPPIV net-positive (more inhibitory than random)
    if not any(a in HYDRO for a in seq): return False
    return True

# read iDPPIV candidate file (first row is header, skip)
rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(DATA, "moso_candidates_idppiv.txt")) if l.strip()]
if rows and rows[0][0].lower().startswith("peptide"):
    rows = rows[1:]
cands = [(p, float(s)) for p, s, *_ in rows]
narrow = [(p, s) for p, s in cands if dpp4_rules(p, s)]
print(f"iDPPIV candidates {len(cands)} -> DPP4 preference narrowed -> {len(narrow)} entries")

# secondary priority: position-2 Pro/Ala marker (optimal for DPP4 S1 pocket)
def p2(p): return len(p) >= 2 and p[1] in "PA"
tier1 = [(p, s) for p, s in narrow if p2(p)]
tier2 = [(p, s) for p, s in narrow if not p2(p)]
print(f"  of which position-2 P/A (optimal): {len(tier1)} | rest: {len(tier2)}")

# output: docking-priority set (tier1 first, both sorted by iDPPIV score descending)
n_top = min(60, len(narrow))
top = (sorted(tier1, key=lambda x: -x[1]) + sorted(tier2, key=lambda x: -x[1]))[:n_top]
_queue = os.path.join(DATA, "moso_dock_queue_idppiv.txt")
with open(_queue, "w") as f:
    for p, s in top:
        f.write(f"{p}\t{s:.3f}\t{'P2' if p2(p) else ''}\n")
print(f"\ndocking queue (top {n_top}) -> {_queue}")
print("samples (first 15):")
for p, s in top[:15]:
    print(f"  {p:6s} iDPPIV={s:+.3f} {'[P2]' if p2(p) else ''}")

# ---- cross-check with old queue: is the best docked peptide LPPQ still selected, overlap how much ----
old = set(l.split("\t")[0].strip() for l in open(os.path.join(DATA, "moso_dock_queue.txt")) if l.strip())
new = set(p for p, _ in top)
print(f"\nold docking queue ({len(old)}) vs new docking queue ({len(new)}): overlap {len(old & new)} entries")
for key in ["LPPQ", "APSPE", "LAPSP", "LPGP"]:
    print(f"  old best {key}: {'present in new queue' if key in new else 'absent from new queue'}")
