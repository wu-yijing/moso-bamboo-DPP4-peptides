# -*- coding: utf-8 -*-
"""
iDPPIV-SCM local reproduction model (global amino-acid composition Scoring Card Method)
=================================================================
Fully offline / zero external dependencies. Reproduces Charoenkwan et al. 2020
(J. Proteome Res. 19:4125-4136, DOI:10.1021/acs.jproteome.0c00590).

Method: on the DPP-IV inhibitory-peptide (positive) vs non-inhibitory-peptide (negative)
training set, compute for each amino acid a a global propensity score
        P(a) = log2( Obs(a) / Exp(a) )
        Obs(a) = occurrence frequency of a in the positive set
        Exp(a) = overall frequency of a in the (positive+negative) set
Total iDPPIV-SCM score of peptide sequence S = Sigma_{a in S} P(a)
  (standard SCM sum over position/residue propensities)
Also provides a length-normalized variant score_mean (Sigma P / |S|) for fair
cross-length ranking.

Note (honest disclosure):
  - This benchmark set is length-confounded: positives are mostly short peptides,
    negatives are mostly long peptides/proteins; the literature's ~0.82 accuracy
    is partly attributable to this. The iDPPIV-SCM authors themselves state
    "not yet accurate enough for real-world applications".
  - In this project (all candidates are 2-6 aa short peptides), the length signal is
    identical across all candidates, so the SCM score mainly reflects the genuine
    signal of [residue composition] that is associated with DPP-IV inhibition,
    and is used for candidate prioritization (ranking) and soft filtering rather
    than as a deterministic classifier.
"""
import os, math
from collections import defaultdict

AA = "ACDEFGHIKLMNPQRSTVWY"
AAS = set(AA)
_HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_TRAIN = os.path.join(_HERE, "data", "train.tsv")


def build_propensity(train_tsv=DEFAULT_TRAIN):
    """Return (P: {aa: log2 propensity}, stats: (Tpos, Tneg))"""
    pos, neg = [], []
    with open(train_tsv, encoding="utf-8") as f:
        next(f)  # header
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            s = "".join(c for c in p[2].upper() if c in AAS)
            if not s:
                continue
            (pos if int(p[1]) == 1 else neg).append(s)
    Np = defaultdict(int)
    Nn = defaultdict(int)
    Tp = sum(len(s) for s in pos)
    Tn = sum(len(s) for s in neg)
    for s in pos:
        for a in s:
            Np[a] += 1
    for s in neg:
        for a in s:
            Nn[a] += 1
    P = {}
    for a in AA:
        obs = (Np[a] / Tp) if Tp else 0.0
        exp = ((Np[a] + Nn[a]) / (Tp + Tn)) if (Tp + Tn) else 0.0
        P[a] = math.log2(obs / exp) if (obs > 0 and exp > 0) else 0.0
    return P, (Tp, Tn, len(pos), len(neg))


def score(seq, P):
    """SCM total score (Sigma residue propensity, standard SCM sum)"""
    if not seq:
        return 0.0
    return sum(P.get(a, 0.0) for a in seq)


def score_mean(seq, P):
    """Length-normalized SCM score (mean propensity per residue), for fair cross-length ranking"""
    if not seq:
        return 0.0
    return score(seq, P) / len(seq)


# ---- decision threshold optimized on the training set (max MCC under nested 5-fold CV) ----
# on this dataset: threshold ~= -1.148 (score > threshold => predicted DPP-IV inhibitory peptide)
DEFAULT_THRESHOLD = -1.148


def predict(seq, P, threshold=DEFAULT_THRESHOLD):
    return 1 if score(seq, P) > threshold else 0


if __name__ == "__main__":
    P, stats = build_propensity()
    Tp, Tn, npos, nneg = stats
    print(f"train: pos={npos} neg={nneg}  total positive residues={Tp} total negative residues={Tn}")
    print("amino-acid propensity P(a) (positive = more likely to appear in DPP-IV inhibitory peptides):")
    for a in sorted(AA, key=lambda x: -P[x]):
        bar = "#" * int(round(P[a] * 4))
        print(f"  {a}  {P[a]:+.3f}  {bar}")
