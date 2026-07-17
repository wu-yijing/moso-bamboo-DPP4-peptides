# -*- coding: utf-8 -*-
"""
Tier0 debias exploration
=========================
Quantify the length confound in the iDPPIV-SCM homologous benchmark and
determine a feasible length-matched debiasing strategy.

Steps:
  1. Merge train.tsv + test.tsv.
  2. Report POS/NEG length distributions (median, mean, min, max, histogram).
  3. Report per-length-bin POS vs NEG counts to find the overlap region where
     length-matched subsampling is possible.
  4. Report the naive length baseline (len<=k -> positive) accuracy on the full
     pooled set to reconfirm the confound.
"""
import os
from collections import Counter, defaultdict
from statistics import median, mean

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")

AA_SET = set("ACDEFGHIKLMNPQRSTVWY")


def load_tsv(path):
    seqs, labels = [], []
    with open(path, encoding="utf-8") as f:
        f.readline()
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            label, seq = parts[1], parts[2].strip()
            seq = "".join(c for c in seq.upper() if c in AA_SET)
            if not seq:
                continue
            seqs.append(seq)
            labels.append(int(label))
    return seqs, labels


def describe(name, lengths):
    lengths = sorted(lengths)
    print(f"  {name}: n={len(lengths)}  min={lengths[0]}  max={lengths[-1]}  "
          f"median={median(lengths)}  mean={mean(lengths):.2f}")


def main():
    seqs, labels = [], []
    for fn in ("train.tsv", "test.tsv"):
        s, y = load_tsv(os.path.join(DATA, fn))
        print(f"{fn}: n={len(s)}  pos={sum(y)}  neg={len(y)-sum(y)}")
        seqs += s
        labels += y

    pos = [s for s, y in zip(seqs, labels) if y == 1]
    neg = [s for s, y in zip(seqs, labels) if y == 0]
    print(f"\nPOOLED: n={len(seqs)}  pos={len(pos)}  neg={len(neg)}")

    print("\n[Length distribution]")
    describe("POS", [len(s) for s in pos])
    describe("NEG", [len(s) for s in neg])

    pos_len = Counter(len(s) for s in pos)
    neg_len = Counter(len(s) for s in neg)
    all_lens = sorted(set(pos_len) | set(neg_len))

    print("\n[Per-length-bin counts]  len : POS / NEG / matched(min)")
    total_matched = 0
    overlap_lens = []
    for L in all_lens:
        p = pos_len.get(L, 0)
        n = neg_len.get(L, 0)
        m = min(p, n)
        total_matched += m
        flag = ""
        if p > 0 and n > 0:
            overlap_lens.append(L)
            flag = "  <-- overlap"
        print(f"  {L:>3} : {p:>4} / {n:>4} / {m:>4}{flag}")

    print(f"\nOverlap lengths (POS&NEG both present): {overlap_lens}")
    print(f"Max length-matched balanced pool: 2 x {total_matched} = {2*total_matched} peptides")

    # Naive length baseline on pooled set: choose best threshold k (len<=k -> pos)
    print("\n[Naive length baseline: len<=k -> positive]")
    best = None
    for k in range(1, max(all_lens) + 1):
        tp = sum(1 for s in pos if len(s) <= k)
        fp = sum(1 for s in neg if len(s) <= k)
        tn = len(neg) - fp
        fn = len(pos) - tp
        acc = (tp + tn) / len(seqs)
        if best is None or acc > best[1]:
            best = (k, acc, tp, fp, tn, fn)
    k, acc, tp, fp, tn, fn = best
    print(f"  best k={k}  ACC={acc:.4f}  tp={tp} fp={fp} tn={tn} fn={fn}")
    print("  => confound reconfirmed: length alone separates the classes.")


if __name__ == "__main__":
    main()
