# -*- coding: utf-8 -*-
"""
Moso-bamboo DPP4 inhibitory-peptide pipeline -- filtering stage
(iDPPIV-SCM offline-reproduction version)
====================================================================
Replaces the original PeptideRanker-style proxy with the literature-validated,
fully offline-reproducible iDPPIV-SCM (global amino-acid composition
Scoring Card Method) score.

Stage 1: iDPPIV-SCM DPP-IV inhibition-propensity scoring (continuous;
         higher = more DPP-IV-inhibitory-peptide-like), and a soft filter
         that predicts "inhibitory peptide" at the training-optimized threshold.
Stage 2: AllerTOP-style allergen prediction (drop high risk)  [still proxy]
Stage 3: ToxinPred-style toxicity prediction (drop high risk)   [still proxy]

Output: funnel counts per stage + final candidate list (with iDPPIV total
        score and length-normalized score, sorted by score descending).
Note: Stage 1 is DPP-IV-specific, literature-validated, offline-reproducible;
      Stages 2/3 still depend on official web tools, and the formal
      manuscript must report the official tool/server outputs.
"""
import sys, os
_HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(_HERE)                      # scripts/ -> repo root
DATA = os.path.join(REPO, "data")
sys.path.insert(0, os.path.join(_HERE, "idppiv_scm"))
from model import build_propensity, score, score_mean, predict, DEFAULT_THRESHOLD

# ---- physicochemical params (stage 2/3 proxy, same as original) ----
KD = {'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,
       'H':-3.2,'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,
       'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
MW = {'A':89,'R':174,'N':132,'D':133,'C':121,'Q':146,'E':147,'G':75,'H':155,
       'I':131,'L':131,'K':146,'M':149,'F':165,'P':115,'S':105,'T':119,'W':204,'Y':181,'V':117}

def aller_top(seq):
    rk = (seq.count('R')+seq.count('K'))/len(seq)
    n  = seq.count('N')/len(seq)
    return (rk*0.6 + n*0.4) > 0.45   # True=high risk (drop)

def toxin_pred(seq):
    if seq.count('C') >= 2: return True
    if (seq.count('R')+seq.count('K')) >= 3: return True
    return False

# ---- load ----
peps = [l.strip() for l in open(os.path.join(DATA, "moso_253_peptides_strict.txt")) if l.strip()]
print(f"input unique peptides = {len(peps)}")

# ---- build iDPPIV-SCM scoring card (from public dataset, offline) ----
P, stats = build_propensity()
print(f"[iDPPIV-SCM] scoring card built (train pos={stats[2]}, neg={stats[3]}; threshold={DEFAULT_THRESHOLD})")

# ---- stage 1 iDPPIV-SCM ----
scored = [(p, score(p, P), score_mean(p, P)) for p in peps]
s1 = [(p, sc, mn) for p, sc, mn in scored if predict(p, P) == 1]
print(f"stage1 iDPPIV-SCM predicted inhibitory peptide : {len(s1)}")

# ---- stage 2 AllerTOP (proxy) ----
s2 = [(p, sc, mn) for p, sc, mn in s1 if not aller_top(p)]
print(f"stage2 drop allergen (AllerTOP): {len(s2)}")

# ---- stage 3 ToxinPred (proxy) ----
s3 = [(p, sc, mn) for p, sc, mn in s2 if not toxin_pred(p)]
print(f"stage3 drop toxicity (ToxinPred) : {len(s3)}")

# ---- sanity check: template known-active peptides should survive ----
refs = ["WPHY","WPQY","VAPGW","WPH","WPQ","VAP"]
print("\nsanity check (template active-peptide fragments, if present in pool):")
for r in refs:
    hit = [p for p,_,_ in s3 if r in p]
    print(f"  {r}: hit {len(hit)} entries" + (f"  e.g. {hit[0]}" if hit else " (not in pool)"))

# ---- output (sorted by iDPPIV total score descending) ----
s3_sorted = sorted(s3, key=lambda x: -x[1])
_cand = os.path.join(DATA, "moso_candidates_idppiv.txt")
with open(_cand,"w") as f:
    f.write("peptide\tiDPPIV_score\tiDPPIV_mean\tpredicted_DPP4_inhibitory\n")
    for p, sc, mn in s3_sorted:
        f.write(f"{p}\t{sc:.3f}\t{mn:.3f}\t1\n")
print(f"\nfinal candidates -> {_cand} ({len(s3)} entries, sorted by iDPPIV score descending)")

# ---- docking queue Top-60 (by iDPPIV total score) ----
top60 = s3_sorted[:60]
_queue = os.path.join(DATA, "moso_dock_queue_idppiv.txt")
with open(_queue,"w") as f:
    for p, sc, mn in top60:
        f.write(f"{p}\t{sc:.3f}\n")
print(f"docking queue Top-60 -> {_queue}")

# ---- funnel cross-check ----
print("\n=== funnel cross-check ===")
print(f"old (PeptideRanker proxy): 4950 -> 2019 -> ... -> 60 (docking)")
print(f"new (iDPPIV-SCM)     : {len(peps)} -> {len(s1)} -> {len(s2)} -> {len(s3)}  (Top-60 into docking)")
