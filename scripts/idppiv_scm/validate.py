# -*- coding: utf-8 -*-
"""Validate whether the local iDPPIV-SCM reproduces the literature accuracy
   Literature: CV~0.819, independent test~0.797 (Charoenkwan 2020)
"""
import sys, os, random
sys.path.insert(0, os.path.dirname(__file__))
from scm import load_tsv, build_scm, score, evaluate

random.seed(42)

tr_seqs, tr_lab = load_tsv("idppiv_scm/data/train.tsv")
te_seqs, te_lab = load_tsv("idppiv_scm/data/test.tsv")
print(f"Training set: {len(tr_seqs)} entries (pos={sum(tr_lab)}, neg={len(tr_lab)-sum(tr_lab)})")
print(f"Test set:     {len(te_seqs)} entries (pos={sum(te_lab)}, neg={len(te_lab)-sum(te_lab)})")

pos = [s for s, y in zip(tr_seqs, tr_lab) if y == 1]
neg = [s for s, y in zip(tr_seqs, tr_lab) if y == 0]
P, L = build_scm(pos, neg)
print(f"SCM scoring-card dims: positions L_max={L}, {len(AA:=set('ACDEFGHIKLMNPQRSTVWY'))} amino acids per position")

# ---- evaluate on the independent test set ----
acc_te, mcc_te, _ = evaluate(te_seqs, te_lab, P, L)
print(f"\n[independent test] ACC={acc_te:.3f}  MCC={mcc_te:.3f}  (literature~0.797)")

# ---- 5-fold cross-validation (within training set) ----
idx = list(range(len(tr_seqs)))
random.shuffle(idx)
folds = [idx[i::5] for i in range(5)]
accs = []
for k in range(5):
    test_idx = set(folds[k])
    tr_idx = [i for i in idx if i not in test_idx]
    p = [tr_seqs[i] for i in tr_idx if tr_lab[i] == 1]
    n = [tr_seqs[i] for i in tr_idx if tr_lab[i] == 0]
    Pk, Lk = build_scm(p, n)
    fold_seqs = [tr_seqs[i] for i in folds[k]]
    fold_lab = [tr_lab[i] for i in folds[k]]
    a, _, _ = evaluate(fold_seqs, fold_lab, Pk, Lk)
    accs.append(a)
cv = sum(accs) / len(accs)
print(f"[5-fold CV] mean ACC={cv:.3f}  (folds: {[f'{a:.3f}' for a in accs]})  (literature~0.819)")

# ---- threshold sensitivity (full model swept on training set, reference only) ----
acc_tr, mcc_tr, _ = evaluate(tr_seqs, tr_lab, P, L)
print(f"[full training set] ACC={acc_tr:.3f}  MCC={mcc_tr:.3f}")

print("\n=> If independent-test ACC falls in 0.75-0.82 and CV in 0.78-0.85, treat it as a faithful reproduction.")
