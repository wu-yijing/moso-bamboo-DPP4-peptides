# -*- coding: utf-8 -*-
"""
Tier0 length-debiased benchmark evaluation of iDPPIV-SCM
========================================================
Goal: separate the *length* signal from the *residue-composition* (true DPP-IV)
signal in the homologous iDPPIV-SCM benchmark, to test whether SCM retains
genuine discriminative power once the length confound is removed.

Design
------
1. Build a *length-matched balanced pool*: for every peptide length L present in
   BOTH the positive and negative sets, keep min(n_pos_L, n_neg_L) peptides from
   each class (random down-sampling). The resulting POS and NEG subsets have an
   *identical* length distribution, so peptide length carries zero class
   information (length baseline -> ~0.5 ACC / ~0.5 AUC).
2. Repeat the random matched sampling R times (different seeds) to average out
   the down-sampling randomness.
3. For each repeat, run stratified k-fold cross-validation: fit the SCM scoring
   card on the training folds, score the held-out fold, and pool out-of-fold
   predictions. Report ACC (threshold 0), MCC and ROC-AUC (threshold-free).
4. Compare four settings:
     (a) confounded benchmark   -- SCM         (CV on the full pooled set)
     (b) confounded benchmark   -- length baseline
     (c) length-matched pool    -- SCM         (CV, averaged over R repeats)
     (d) length-matched pool    -- length baseline
   Expectation: (b) high, (d) ~0.5; (c) > (d) proves SCM captures real
   composition signal beyond length.

Fully offline, standard library only. Deterministic (fixed master seed).
"""
import os
import random
from collections import Counter, defaultdict
from statistics import mean, pstdev

import scm  # local module (build_scm, score, load_tsv, evaluate)

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
OUT = os.path.join(HERE, "tier0_results.tsv")

MASTER_SEED = 20260717
R_REPEATS = 50
K_FOLDS = 5


# ----------------------------- data loading ----------------------------------
def load_pool():
    seqs, labels = [], []
    for fn in ("train.tsv", "test.tsv"):
        s, y = scm.load_tsv(os.path.join(DATA, fn))
        seqs += s
        labels += y
    return seqs, labels


# --------------------------- metric utilities --------------------------------
def confusion(preds, labels):
    tp = fp = tn = fn = 0
    for p, y in zip(preds, labels):
        if y == 1 and p == 1: tp += 1
        elif y == 0 and p == 1: fp += 1
        elif y == 0 and p == 0: tn += 1
        else: fn += 1
    return tp, fp, tn, fn


def acc_mcc(preds, labels):
    import math
    tp, fp, tn, fn = confusion(preds, labels)
    n = tp + fp + tn + fn
    acc = (tp + tn) / n if n else 0.0
    denom = math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
    mcc = ((tp*tn) - (fp*fn)) / denom if denom > 0 else 0.0
    return acc, mcc


def roc_auc(scores, labels):
    """Rank-based AUC (Mann-Whitney U). Ties get average rank."""
    paired = sorted(zip(scores, labels), key=lambda t: t[0])
    # assign average ranks
    ranks = [0.0] * len(paired)
    i = 0
    while i < len(paired):
        j = i
        while j + 1 < len(paired) and paired[j+1][0] == paired[i][0]:
            j += 1
        avg = (i + j) / 2.0 + 1.0  # 1-based average rank
        for k in range(i, j + 1):
            ranks[k] = avg
        i = j + 1
    n_pos = sum(1 for _, y in paired if y == 1)
    n_neg = len(paired) - n_pos
    if n_pos == 0 or n_neg == 0:
        return 0.5
    sum_ranks_pos = sum(r for r, (_, y) in zip(ranks, paired) if y == 1)
    u = sum_ranks_pos - n_pos * (n_pos + 1) / 2.0
    return u / (n_pos * n_neg)


# ------------------------- length-matched sampling ---------------------------
def length_matched_pool(pos, neg, rng):
    """Down-sample each length bin so POS and NEG share an identical length
    distribution. Returns (mseqs, mlabels)."""
    pos_by_len = defaultdict(list)
    neg_by_len = defaultdict(list)
    for s in pos:
        pos_by_len[len(s)].append(s)
    for s in neg:
        neg_by_len[len(s)].append(s)
    mseqs, mlabels = [], []
    for L in sorted(set(pos_by_len) & set(neg_by_len)):
        m = min(len(pos_by_len[L]), len(neg_by_len[L]))
        ps = rng.sample(pos_by_len[L], m)
        ns = rng.sample(neg_by_len[L], m)
        mseqs += ps + ns
        mlabels += [1]*m + [0]*m
    return mseqs, mlabels


# ------------------------------ CV harness -----------------------------------
def stratified_folds(labels, k, rng):
    idx_pos = [i for i, y in enumerate(labels) if y == 1]
    idx_neg = [i for i, y in enumerate(labels) if y == 0]
    rng.shuffle(idx_pos)
    rng.shuffle(idx_neg)
    folds = [[] for _ in range(k)]
    for f, i in enumerate(idx_pos):
        folds[f % k].append(i)
    for f, i in enumerate(idx_neg):
        folds[f % k].append(i)
    return folds


def cv_scm(seqs, labels, k, rng):
    """Stratified k-fold CV of SCM. Returns pooled out-of-fold (acc, mcc, auc)."""
    folds = stratified_folds(labels, k, rng)
    oof_scores = [0.0] * len(seqs)
    oof_preds = [0] * len(seqs)
    for f in range(k):
        test_idx = set(folds[f])
        tr_pos = [seqs[i] for i in range(len(seqs)) if i not in test_idx and labels[i] == 1]
        tr_neg = [seqs[i] for i in range(len(seqs)) if i not in test_idx and labels[i] == 0]
        if not tr_pos or not tr_neg:
            continue
        P, L_max = scm.build_scm(tr_pos, tr_neg)
        for i in test_idx:
            sc = scm.score(seqs[i], P, L_max)
            oof_scores[i] = sc
            oof_preds[i] = 1 if sc > 0.0 else 0
    acc, mcc = acc_mcc(oof_preds, labels)
    auc = roc_auc(oof_scores, labels)
    return acc, mcc, auc


def cv_length_baseline(seqs, labels, k, rng):
    """Length baseline under CV: on each training fold pick the length threshold
    k* (len<=k* -> positive) maximizing training ACC, apply to test fold.
    AUC uses -length as the score (shorter = more positive)."""
    folds = stratified_folds(labels, k, rng)
    oof_preds = [0] * len(seqs)
    oof_scores = [0.0] * len(seqs)
    max_len = max(len(s) for s in seqs)
    for f in range(k):
        test_idx = set(folds[f])
        tr = [(seqs[i], labels[i]) for i in range(len(seqs)) if i not in test_idx]
        best_k, best_acc = 1, -1.0
        for kk in range(1, max_len + 1):
            tp = sum(1 for s, y in tr if y == 1 and len(s) <= kk)
            fp = sum(1 for s, y in tr if y == 0 and len(s) <= kk)
            tn = sum(1 for s, y in tr if y == 0) - fp
            fn = sum(1 for s, y in tr if y == 1) - tp
            acc = (tp + tn) / len(tr)
            if acc > best_acc:
                best_acc, best_k = acc, kk
        for i in test_idx:
            oof_preds[i] = 1 if len(seqs[i]) <= best_k else 0
            oof_scores[i] = -len(seqs[i])
    acc, mcc = acc_mcc(oof_preds, labels)
    auc = roc_auc(oof_scores, labels)
    return acc, mcc, auc


# --------------------------------- main --------------------------------------
def main():
    seqs, labels = load_pool()
    pos = [s for s, y in zip(seqs, labels) if y == 1]
    neg = [s for s, y in zip(seqs, labels) if y == 0]
    print(f"Pooled benchmark: n={len(seqs)} pos={len(pos)} neg={len(neg)}")

    master = random.Random(MASTER_SEED)

    rows = []  # (setting, metric, mean, sd, n_pool)

    # (a)/(b) confounded benchmark, single CV (repeat R times with reshuffle for CI)
    a_acc, a_mcc, a_auc = [], [], []
    b_acc, b_mcc, b_auc = [], [], []
    for r in range(R_REPEATS):
        rng = random.Random(master.randint(0, 10**9))
        acc, mcc, auc = cv_scm(seqs, labels, K_FOLDS, rng)
        a_acc.append(acc); a_mcc.append(mcc); a_auc.append(auc)
        rng2 = random.Random(master.randint(0, 10**9))
        acc, mcc, auc = cv_length_baseline(seqs, labels, K_FOLDS, rng2)
        b_acc.append(acc); b_mcc.append(mcc); b_auc.append(auc)

    # (c)/(d) length-matched pool, R repeats of matched sampling + CV
    c_acc, c_mcc, c_auc = [], [], []
    d_acc, d_mcc, d_auc = [], [], []
    pool_sizes = []
    for r in range(R_REPEATS):
        rng = random.Random(master.randint(0, 10**9))
        mseqs, mlabels = length_matched_pool(pos, neg, rng)
        pool_sizes.append(len(mseqs))
        acc, mcc, auc = cv_scm(mseqs, mlabels, K_FOLDS, rng)
        c_acc.append(acc); c_mcc.append(mcc); c_auc.append(auc)
        acc, mcc, auc = cv_length_baseline(mseqs, mlabels, K_FOLDS, rng)
        d_acc.append(acc); d_mcc.append(mcc); d_auc.append(auc)

    def summ(name, pool, accs, mccs, aucs):
        print(f"\n[{name}]  (pool n={pool})")
        print(f"  ACC = {mean(accs):.4f} +/- {pstdev(accs):.4f}")
        print(f"  MCC = {mean(mccs):.4f} +/- {pstdev(mccs):.4f}")
        print(f"  AUC = {mean(aucs):.4f} +/- {pstdev(aucs):.4f}")
        rows.append((name, pool, mean(accs), pstdev(accs),
                     mean(mccs), pstdev(mccs), mean(aucs), pstdev(aucs)))

    summ("Confounded / SCM",            len(seqs), a_acc, a_mcc, a_auc)
    summ("Confounded / length-baseline", len(seqs), b_acc, b_mcc, b_auc)
    summ("Length-matched / SCM",         int(mean(pool_sizes)), c_acc, c_mcc, c_auc)
    summ("Length-matched / length-baseline", int(mean(pool_sizes)), d_acc, d_mcc, d_auc)

    # (e) label-permutation null on the length-matched pool: shuffle labels and
    #     re-run SCM CV. A valid model's real AUC should exceed this null.
    perm_auc = []
    for r in range(R_REPEATS):
        rng = random.Random(master.randint(0, 10**9))
        mseqs, mlabels = length_matched_pool(pos, neg, rng)
        shuffled = mlabels[:]
        rng.shuffle(shuffled)
        _, _, auc = cv_scm(mseqs, shuffled, K_FOLDS, rng)
        perm_auc.append(auc)
    print("\n[Length-matched / SCM label-permutation null]")
    print(f"  AUC = {mean(perm_auc):.4f} +/- {pstdev(perm_auc):.4f}  (target ~0.5)")
    obs = mean(c_auc)
    n_ge = sum(1 for a in perm_auc if a >= obs)
    p_emp = (n_ge + 1) / (len(perm_auc) + 1)
    print(f"  observed matched-SCM AUC = {obs:.4f};  empirical p (null>=obs) = {p_emp:.4f}")
    rows.append(("Length-matched / SCM permutation-null", int(mean(pool_sizes)),
                 0.0, 0.0, 0.0, 0.0, mean(perm_auc), pstdev(perm_auc)))

    # write TSV
    with open(OUT, "w", encoding="utf-8") as f:
        f.write("setting\tpool_n\tACC_mean\tACC_sd\tMCC_mean\tMCC_sd\tAUC_mean\tAUC_sd\n")
        for row in rows:
            f.write("\t".join(
                (row[0], str(row[1]),
                 f"{row[2]:.4f}", f"{row[3]:.4f}",
                 f"{row[4]:.4f}", f"{row[5]:.4f}",
                 f"{row[6]:.4f}", f"{row[7]:.4f}")) + "\n")
    print(f"\nWritten: {OUT}")

    # headline interpretation
    print("\n[Interpretation]")
    print(f"  Length signal removed?  matched length-baseline ACC = "
          f"{mean(d_acc):.3f} (target ~0.5), AUC = {mean(d_auc):.3f} (target ~0.5)")
    print(f"  SCM residual (true) signal?  matched SCM AUC = {mean(c_auc):.3f} "
          f"vs matched baseline AUC = {mean(d_auc):.3f}")
    lift = mean(c_auc) - mean(d_auc)
    print(f"  De-confounded AUC lift over length baseline = {lift:+.3f}")


if __name__ == "__main__":
    main()
