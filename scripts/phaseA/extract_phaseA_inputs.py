#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Phase A - Step 1: 提取毛竹虚拟消化产生的 2-6aa 短肽，作为官方服务器
(PeptideRanker / AllerTOP v2 / ToxinPred3) 的输入。

输入: E:/workbuddy/Claw/moso_253_peptides_strict.txt  (每行一条肽)
过滤: 长度 2-6, 仅标准 20 氨基酸, 去重
输出: data/phaseA_inputs/
  - moso_short_2to6.fasta      (主输入, 全部 2-6aa, 带 ID)
  - moso_short_2to6.txt        (简单列表, 每行一条, 供部分服务器)
  - <server>/batch_XX.fasta     (每批 500, 三服务器各一份)
"""
import os
import math

SRC = r"E:/workbuddy/Claw/moso_253_peptides_strict.txt"
OUT_DIR = r"E:/workbuddy/moso-bamboo-DPP4-peptides/data/phaseA_inputs"
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
