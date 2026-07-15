#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plan (1): iDPPIV-SI-style second independent predictor cross-validation
================================================
Purpose: use a sequence classifier trained independently and differently from iDPPIV-SCM
      to re-score our candidate pool, and report the two-model agreement
      (Spearman correlation of continuous scores + Top-K overlap rate) to compensate
      for SCM's benchmark length-confound / limited-accuracy weakness.

Methodology (faithfully reproduce iDPPIV-SI's "feature selection + SVM" paradigm,
offline zero-external-dependency):
  - features: amino-acid composition (20) + dipeptide composition (400) + 1st/2nd-order autocorrelation (40)
          = 460 dims (iDPPIV-SI uses 50 physicochemical properties + 1/2-order correlation + DWT;
          here composition-type features serve as the equivalent reproducible physicochemical-property carrier, see manuscript honest disclosure)
  - feature selection: LASSO (LassoCV, 5-fold)
  - classifier: RBF-SVM (class_weight=balanced, small grid for C)
  - train/eval: homologous public 665+665 benchmark (train/test.tsv), report test-set ACC/Sens/Spec/MCC

Honest statement: Zou H (2024)'s original model weights are stored as figshare MATLAB .mat;
this environment has no MATLAB, so this is a methodology-consistent Python offline
reproduction (iDPPIV-SI-style), not a byte-exact re-implementation.
"""
import csv, sys, itertools
import numpy as np
from scipy.stats import spearmanr
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LassoCV
from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.metrics import accuracy_score, recall_score, precision_score, matthews_corrcoef

AA = "ACDEFGHIKLMNPQRSTVWY"
AA_IDX = {a: i for i, a in enumerate(AA)}
PAIRS = ["".join(p) for p in itertools.product(AA, repeat=2)]

def aa_comp(seq):
    v = np.zeros(20)
    L = len(seq)
    if L == 0:
        return v
    for c in seq:
        if c in AA_IDX:
            v[AA_IDX[c]] += 1
    return v / L

def dipep_comp(seq):
    v = np.zeros(400)
    if len(seq) < 2:
        return v
    for i in range(len(seq) - 1):
        d = seq[i:i+2]
        if d in PAIRS_MAP:
            v[PAIRS_MAP[d]] += 1
    return v / (len(seq) - 1)

def autocorr(seq, lag):
    # for each amino-acid type a, compute the indicator-property Moreau-Broto normalized autocorrelation:
    #   AC_a(lag) = sum_i [res_i==a and res_{i+lag}==a] / (L - lag)
    # returns 20 dims (equivalent carrier of iDPPIV-SI's '1-order / 2-order correlation')
    L = len(seq)
    v = np.zeros(20)
    if L <= lag:
        return v
    denom = L - lag
    for i in range(L - lag):
        a = AA_IDX.get(seq[i]); b = AA_IDX.get(seq[i + lag])
        if a is not None and b is not None and a == b:
            v[a] += 1.0
    return v / denom

def aa_comp_single(c):
    v = np.zeros(20)
    if c in AA_IDX:
        v[AA_IDX[c]] = 1
    return v

def pc_props(seq):
    # precisely computable physicochemical properties (equivalent reproducible subset of iDPPIV-SI's '50 PC properties')
    L = len(seq) or 1
    net = sum(1 for c in seq if c in "KRH") - sum(1 for c in seq if c in "DE")
    hyd = sum(1 for c in seq if c in "AVLIMFWY") / L
    pro = seq.count("P") / L
    chg = sum(1 for c in seq if c in "KRHDE") / L
    aro = sum(1 for c in seq if c in "FWY") / L
    pol = sum(1 for c in seq if c in "STNQCY") / L
    ali = sum(1 for c in seq if c in "AVLIM") / L
    return np.array([net, hyd, pro, chg, aro, pol, ali], dtype=float)

def featurize(seq):
    return np.concatenate([aa_comp(seq), dipep_comp(seq),
                           autocorr(seq, 1), autocorr(seq, 2), pc_props(seq)])

PAIRS_MAP = {p: i for i, p in enumerate(PAIRS)}

def load_bench(path):
    seqs, labels = [], []
    with open(path) as f:
        r = csv.reader(f, delimiter="\t")
        next(r)
        for row in r:
            if len(row) < 3:
                continue
            seqs.append(row[2].strip())
            labels.append(int(row[1]))
    return seqs, np.array(labels)

def load_candidates(path):
    rows = []
    with open(path) as f:
        r = csv.reader(f, delimiter="\t")
        header = next(r)
        for row in r:
            if len(row) < 4:
                continue
            rows.append(row)
    return header, rows

def main():
    base = "scripts/idppiv_scm/data"
    tr_seq, tr_lab = load_bench(f"{base}/train.tsv")
    te_seq, te_lab = load_bench(f"{base}/test.tsv")
    print(f"[bench] train={len(tr_seq)} ({tr_lab.sum()} pos)  test={len(te_seq)} ({te_lab.sum()} pos)", flush=True)

    Xtr = np.array([featurize(s) for s in tr_seq])
    Xte = np.array([featurize(s) for s in te_seq])
    print(f"[feat] dim={Xtr.shape[1]}", flush=True)

    scaler = StandardScaler().fit(Xtr)
    Xtr_s = scaler.transform(Xtr)
    Xte_s = scaler.transform(Xte)

    # LASSO feature selection
    lasso = LassoCV(cv=5, max_iter=5000).fit(Xtr_s, tr_lab)
    mask = lasso.coef_ != 0
    if mask.sum() == 0:
        mask = np.ones(Xtr_s.shape[1], dtype=bool)
    print(f"[lasso] selected {mask.sum()}/{mask.shape[0]} features", flush=True)
    Xtr_l = Xtr_s[:, mask]
    Xte_l = Xte_s[:, mask]

    # RBF-SVM small grid for C
    svc = GridSearchCV(SVC(kernel="rbf", class_weight="balanced", probability=True),
                       {"C": [0.1, 1.0, 10.0]}, cv=StratifiedKFold(3, shuffle=True, random_state=0))
    svc.fit(Xtr_l, tr_lab)
    clf = svc.best_estimator_
    print(f"[svm] best C={svc.best_params_['C']}", flush=True)

    te_pred = clf.predict(Xte_l)
    te_prob = clf.predict_proba(Xte_l)[:, 1]
    acc = accuracy_score(te_lab, te_pred)
    sens = recall_score(te_lab, te_pred, pos_label=1)
    spec = recall_score(te_lab, te_pred, pos_label=0)
    mcc = matthews_corrcoef(te_lab, te_pred)
    print(f"[eval] test ACC={acc:.3f} Sens={sens:.3f} Spec={spec:.3f} MCC={mcc:.3f}", flush=True)

    # re-score the candidate pool
    header, cands = load_candidates("data/moso_candidates_idppiv.txt")
    c_seq = [r[0] for r in cands]
    scm_score = np.array([float(r[1]) for r in cands])
    scm_mean = np.array([float(r[2]) for r in cands])
    scm_pred = np.array([int(r[3]) for r in cands])
    Xc = scaler.transform(np.array([featurize(s) for s in c_seq]))[:, mask]
    si_score = clf.decision_function(Xc)
    si_prob = clf.predict_proba(Xc)[:, 1]
    si_pred = (si_prob >= 0.5).astype(int)
    print(f"[score] candidates={len(c_seq)}  si_prob>0.5: {si_pred.sum()}", flush=True)

    # agreement: continuous-score Spearman
    rho, prho = spearmanr(scm_score, si_score)
    print(f"[agree] Spearman(scm_score, si_score) rho={rho:.3f} (p={prho:.1e})", flush=True)

    # Top-K overlap
    def topk(idx, k):
        return set(np.argsort(idx)[::-1][:k])
    for k in (50, 100, 200):
        a = topk(scm_score, k); b = topk(si_score, k)
        ov = len(a & b)
        print(f"[agree] Top-{k} overlap: {ov}/{k} ({ov/k:.0%})", flush=True)

    # locate the three finalist peptides in both models
    si_rank = {c_seq[idx]: r + 1 for r, idx in enumerate(np.argsort(si_score)[::-1])}
    scm_rank = {c_seq[idx]: r + 1 for r, idx in enumerate(np.argsort(scm_score)[::-1])}
    finals = {"LPPGP": None, "APPSQ": None, "APQIP": None}
    for i, s in enumerate(c_seq):
        if s in finals:
            finals[s] = (scm_score[i], si_score[i], si_prob[i],
                         si_rank[s], scm_rank[s])
    print("[finals]", flush=True)
    for name in finals:
        if finals[name]:
            sc, ss, sp, rk_si, rk_scm = finals[name]
            print(f"  {name}: scm={sc:.3f}  si_score={ss:+.3f}  si_prob={sp:.3f}  rankSI={rk_si}/{len(c_seq)}  rankSCM={rk_scm}/{len(c_seq)}", flush=True)

    # write TSV
    out = "data/phaseA/idppiv_si_crosscheck.tsv"
    with open(out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["peptide", "iDPPIV_score", "iDPPIV_mean", "scm_pred",
                    "si_score", "si_prob", "si_pred"])
        order = np.argsort(si_score)[::-1]
        for i in order:
            w.writerow([c_seq[i], f"{scm_score[i]:.4f}", f"{scm_mean[i]:.4f}",
                        scm_pred[i], f"{si_score[i]:.4f}", f"{si_prob[i]:.4f}", si_pred[i]])
    # summary line (for manuscript citation)
    with open("data/phaseA/idppiv_si_crosscheck_summary.txt", "w") as f:
        f.write(f"test_ACC\t{acc:.3f}\n")
        f.write(f"test_Sens\t{sens:.3f}\n")
        f.write(f"test_Spec\t{spec:.3f}\n")
        f.write(f"test_MCC\t{mcc:.3f}\n")
        f.write(f"spearman_rho\t{rho:.3f}\n")
        f.write(f"spearman_p\t{prho:.2e}\n")
        f.write(f"candidates_n\t{len(c_seq)}\n")
        f.write(f"si_pos_n\t{int(si_pred.sum())}\n")
        for k in (50, 100, 200):
            a = topk(scm_score, k); b = topk(si_score, k)
            f.write(f"top{k}_overlap\t{len(a & b)}\n")
        for name in finals:
            if finals[name]:
                sc, ss, sp, rk_si, rk_scm = finals[name]
                f.write(f"final_{name}\tscm={sc:.3f}\tsi_score={ss:+.3f}\tsi_prob={sp:.3f}\trankSI={rk_si}/{len(c_seq)}\trankSCM={rk_scm}/{len(c_seq)}\n")
    print(f"[done] wrote {out} + summary", flush=True)

if __name__ == "__main__":
    main()
