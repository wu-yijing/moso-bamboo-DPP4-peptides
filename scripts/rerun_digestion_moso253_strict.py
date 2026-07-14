#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
毛竹 模拟消化 — 严格按模板论文(Chenget al. Bioorg Chem 175,2026,大蒜DPP4)的 PeptideCutter 规则:
  pepsin(pH1.3):     切 F/Y/W/L 的 N 端 (这些残基成为片段起点)
  trypsin:            切 K/R 的 C 端 (K,R 后断开)
  chymotrypsin(spec): 切 F/Y/W 的 C 端 (L 不切；Pro 前不停)
合并三酶切点 -> 片段；保留 2..20 aa 且不含 X 的唯一肽。
"""
import os
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # scripts/ -> 仓库根
DATA = os.path.join(REPO, "data")

def parse_fasta(path):
    seqs={}; name=None; buf=[]
    for line in open(path):
        line=line.rstrip("\n")
        if line.startswith(">"):
            if name: seqs[name]="".join(buf)
            name=line[1:].split()[0]; buf=[]
        else: buf.append(line.strip())
    if name: seqs[name]="".join(buf)
    return seqs

seqs=parse_fasta(os.path.join(DATA, "moso_253.fasta"))
n=len(seqs); tot=sum(len(s) for s in seqs.values())
print(f"蛋白={n}  总残基={tot:,}  均={tot/n:.0f}aa")

pep=set("FYWL")   # N端前切
tryp=set("KR")    # C端后切
chym=set("FYW")   # C端后切(不含L)

def cut_points(seq):
    pts={0,len(seq)}
    for i,a in enumerate(seq):
        if a in pep:  pts.add(i)      # a前切 -> a为起点
        if a in tryp: pts.add(i+1)    # a后切 -> i+1为起点
        if a in chym: pts.add(i+1)
    return sorted(pts)

ALL=set()
for s in seqs.values():
    if not s: continue
    pts=cut_points(s)
    for a,b in zip(pts[:-1],pts[1:]):
        fr=s[a:b]
        if fr and "X" not in fr: ALL.add(fr)

from collections import Counter
lens=Counter(len(p) for p in ALL)
p26=sum(1 for p in ALL if 2<=len(p)<=6)
p220=sum(1 for p in ALL if 2<=len(p)<=20)
print(f"\n唯一肽(全部,无X) = {len(ALL)}")
print(f"  2..6 aa  = {p26}")
print(f"  2..20 aa = {p220}")
print("长度分布(2..12):",{L:lens.get(L,0) for L in range(2,13)})

print("\n=== 对比 ===")
print(f"{'指标':<22}{'毛竹(严格)':>14}{'大蒜模板':>14}")
print(f"{'蛋白':<22}{n:>14}{113:>14}")
print(f"{'唯一肽全部':<22}{len(ALL):>14}{5672:>14}")
print(f"{'唯一 2..6aa':<22}{p26:>14}{1442:>14}")
print(f"{'PR>0.5候选(估40%)':<22}{int(p26*0.4):>14}{249:>14}")

_out = os.path.join(DATA, "moso_253_peptides_strict.txt")
with open(_out,"w") as f:
    for p in sorted(ALL,key=lambda x:(len(x),x)): f.write(p+"\n")
print(f"\n-> {_out} ({len(ALL)} 条)")
