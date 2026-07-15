# -*- coding: utf-8 -*-
"""
iDPPIV-SCM local reproduction (Scoring Card Method, SCM)
=====================================================
Fully offline, zero external dependencies. Reproduces Charoenkwan et al. 2020
(J. Proteome Res. 19:4125-4136, DOI:10.1021/acs.jproteome.0c00590)
proposed iDPPIV-SCM scoring-card method.

SCM core:
  for position i (0-based) and amino acid a, propensity score
      P(i,a) = log2( Obs(i,a) / Exp(i,a) )
  where (counted on positive / negative training sets):
      Obs(i,a) = Npos(i,a) / Npos(i)       # frequency of a at position i in positive set
      Exp(i,a) = (Npos(i,a)+Nneg(i,a)) / (Npos(i)+Nneg(i))  # expected frequency
  peptide sequence S=(a1..aL) total SCM score:
      Score(S) = Sigma_{i=0}^{L-1} P(i, a_{i+1})
  decision: Score > 0 predicts DPP-IV inhibitory peptide (positive).
  (threshold 0 comes from log2(1)=0; a single-site contribution is positive only
   when Obs>Exp, consistent with the original SCM paper definition.)

Variable-length handling: only peptides with length >= i+1 contribute the i-th term;
short peptides naturally accumulate only their valid positions.
"""
import math
from collections import defaultdict

AA = "ACDEFGHIKLMNPQRSTVWY"
AA_SET = set(AA)


def load_tsv(path):
    """Read index\tlabel\tsequence format, return (seqs, labels)"""
    seqs, labels = [], []
    with open(path, encoding="utf-8") as f:
        header = f.readline()  # skip header
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                # degenerate: sequence only
                parts = (None, None, line)
            _, label, seq = parts[0], parts[1], parts[2].strip()
            seq = "".join(c for c in seq.upper() if c in AA_SET)
            if not seq:
                continue
            seqs.append(seq)
            labels.append(int(label))
    return seqs, labels


def build_scm(pos_seqs, neg_seqs, L_max=None):
    """Compute SCM scoring card P[i][a] from positive/negative sequence sets"""
    if L_max is None:
        L_max = max([len(s) for s in pos_seqs + neg_seqs] + [1])

    Npos = defaultdict(lambda: defaultdict(int))   # Npos[i][a]
    Nneg = defaultdict(lambda: defaultdict(int))
    Npos_i = defaultdict(int)                       # num positive-set peptides at position i (length>=i+1)
    Nneg_i = defaultdict(int)

    for s in pos_seqs:
        for i, a in enumerate(s):
            Npos[i][a] += 1
            Npos_i[i] += 1
    for s in neg_seqs:
        for i, a in enumerate(s):
            Nneg[i][a] += 1
            Nneg_i[i] += 1

    P = {}
    for i in range(L_max):
        P[i] = {}
        denom = Npos_i[i] + Nneg_i[i]
        for a in AA:
            obs = (Npos[i][a] / Npos_i[i]) if Npos_i[i] > 0 else 0.0
            exp = ((Npos[i][a] + Nneg[i][a]) / denom) if denom > 0 else 0.0
            if obs > 0 and exp > 0:
                P[i][a] = math.log2(obs / exp)
            else:
                P[i][a] = 0.0  # unobserved -> neutral
    return P, L_max


def score(seq, P, L_max):
    """SCM total score (sum of log2 propensities)"""
    s = 0.0
    for i, a in enumerate(seq):
        if i >= L_max:
            break
        s += P[i].get(a, 0.0)
    return s


def predict(seq, P, L_max, threshold=0.0):
    return 1 if score(seq, P, L_max) > threshold else 0


def evaluate(seqs, labels, P, L_max, threshold=0.0):
    """Return (ACC, MCC, prediction_list)"""
    tp = fp = tn = fn = 0
    preds = []
    for s, y in zip(seqs, labels):
        p = predict(s, P, L_max, threshold)
        preds.append(p)
        if y == 1 and p == 1: tp += 1
        elif y == 0 and p == 1: fp += 1
        elif y == 0 and p == 0: tn += 1
        else: fn += 1
    acc = (tp + tn) / (tp + fp + tn + fn) if (tp+fp+tn+fn) else 0.0
    denom = math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
    mcc = ((tp*tn) - (fp*fn)) / denom if denom > 0 else 0.0
    return acc, mcc, preds
