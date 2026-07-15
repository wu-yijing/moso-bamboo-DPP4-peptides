#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Phase A - Step 1: extract the 2-6 aa short peptides produced by the moso-bamboo
in silico digestion, to serve as input for the official web servers
(PeptideRanker / AllerTOP v2 / ToxinPred3).

Input : data/moso_253_peptides_strict.txt  (one peptide per line)
Filter: length 2-6, standard 20 amino acids only, deduplicated
Output: data/phaseA_inputs/
  - moso_short_2to6.fasta      (main input, all 2-6 aa, with IDs)
  - moso_short_2to6.txt        (plain list, one per line, for some servers)
  - <server>/batch_XX.fasta     (500 per batch, one copy per server)
"""
import os
import math

# --- repo-relative paths (scripts/phaseA/ -> repo root) ---
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SRC = os.path.join(REPO, "data", "moso_253_peptides_strict.txt")
OUT_DIR = os.path.join(REPO, "data", "phaseA_inputs")
BATCH = 500
STANDARD = set("ACDEFGHIKLMNPQRSTVWY")
SERVERS = ["peptideranker", "allertop", "toxinpred"]


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    peps = []
    with open(SRC, encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if len(s) < 2 or len(s) > 6:
                continue
            if any(c not in STANDARD for c in s):
                continue
            peps.append(s)

    seen = set()
    uniq = []
    for p in peps:
        if p in seen:
            continue
        seen.add(p)
        uniq.append(p)

    print(f"[parsed] total 2-6aa std peptides : {len(peps)}")
    print(f"[unique] deduplicated            : {len(uniq)}")

    # master FASTA
    master_fa = os.path.join(OUT_DIR, "moso_short_2to6.fasta")
    with open(master_fa, "w") as f:
        for i, p in enumerate(uniq, 1):
            f.write(f">pep_{i:05d}\n{p}\n")

    # simple list
    simple = os.path.join(OUT_DIR, "moso_short_2to6.txt")
    with open(simple, "w") as f:
        for p in uniq:
            f.write(p + "\n")

    # batches per server
    n = len(uniq)
    nb = math.ceil(n / BATCH)
    for srv in SERVERS:
        d = os.path.join(OUT_DIR, srv)
        os.makedirs(d, exist_ok=True)
        for b in range(nb):
            chunk = uniq[b * BATCH:(b + 1) * BATCH]
            bf = os.path.join(d, f"batch_{b + 1:02d}.fasta")
            with open(bf, "w") as f:
                for i, p in enumerate(chunk, b * BATCH + 1):
                    f.write(f">pep_{i:05d}\n{p}\n")

    print(f"[done] master + simple + {nb} batches x {len(SERVERS)} servers written to:")
    print(f"       {OUT_DIR}")


if __name__ == "__main__":
    main()
