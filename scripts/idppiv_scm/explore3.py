# -*- coding: utf-8 -*-
"""Recipe probe v3: length-confounding analysis + length-normalized SCM + odds-ratio variant"""
import sys, os, math, random
from collections import defaultdict
sys.path.insert(0, os.path.dirname(__file__))
from scm import load_tsv
AA = "ACDEFGHIKLMNPQRSTVWY"; AAS = set(AA)

trS, trL = load_tsv("idppiv_scm/data/train.tsv")
teS, teL = load_tsv("idppiv_scm/data/test.tsv")
pos = [s for s, l in zip(trS, trL) if l == 1]; neg = [s for s, l in zip(trS, trL) if l == 0]

def acc_at(seqs, labs, preds):
    return sum(1 for y, p in zip(labs, preds) if p == y) / len(labs)

# ---- (1) pure-length baseline ----
print("=== (1) pure-length baseline (predict short = positive) ===")
for T in [5, 6, 8, 10, 12, 15]:
    p_tr = [1 if len(s) <= T else 0 for s in trS]; p_te = [1 if len(s) <= T else 0 for s in teS]
    print(f"  len<={T}: trainACC={acc_at(trS, trL, p_tr):.3f}  testACC={acc_at(teS, teL, p_te):.3f}")

# ---- build two global cards: Obs/Exp and Odds(pos/neg) ----
def build_global(pos, neg, mode):
    Np = defaultdict(int); Nn = defaultdict(int)
    Tp = sum(len(s) for s in pos); Tn = sum(len(s) for s in neg)
    for s in pos:
        for a in s: Np[a] += 1
    for s in neg:
        for a in s: Nn[a] += 1
    P = {}
    for a in AA:
        fp = Np[a] / Tp if Tp else 0; fn = Nn[a] / Tn if Tn else 0
        if mode == "obs/exp":
            exp = ((Np[a] + Nn[a]) / (Tp + Tn)) if (Tp + Tn) else 0
            P[a] = math.log2(fp / exp) if fp > 0 and exp > 0 else 0.0
        else:  # odds
            P[a] = math.log2(fp / fn) if fp > 0 and fn > 0 else 0.0
    return P

def sc_sum(seq, P):
    return sum(P.get(a, 0.0) for a in seq)

def sc_mean(seq, P):
    return sc_sum(seq, P) / len(seq) if seq else 0.0

def tune(seqs, labs, scores):
    cand = sorted(set(scores)); best = (0, 0)
    ts = [min(scores) - 1e-9] + [(cand[i] + cand[i + 1]) / 2 for i in range(len(cand) - 1)] + [max(scores) + 1e-9]
    for t in ts:
        tp = fp = tn = fn = 0
        for sc, y in zip(scores, labs):
            p = 1 if sc > t else 0
            tp += (y == 1 and p == 1); fp += (y == 0 and p == 1); tn += (y == 0 and p == 0); fn += (y == 1 and p == 0)
        den = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
        mcc = ((tp * tn) - (fp * fn)) / den if den else 0
        if mcc > best[0]:
            best = (mcc, t)
    return best[1]

for mode in ["obs/exp", "odds"]:
    for agg in ["sum", "mean"]:
        P = build_global(pos, neg, mode)
        sc = lambda s: (sc_sum if agg == "sum" else sc_mean)(s, P)
        s_tr = [sc(s) for s in trS]; s_te = [sc(s) for s in teS]
        t = tune(trS, trL, s_tr)
        print(f"[global {mode} {agg}] thr={t:.3f} trainACC={acc_at(trS, trL, [1 if x > t else 0 for x in s_tr]):.3f} testACC={acc_at(teS, teL, [1 if x > t else 0 for x in s_te]):.3f}")
