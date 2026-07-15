#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Moso-bamboo simulated digestion -- strictly following the template paper
(Chenget al. Bioorg Chem 175, 2026, garlic DPP4) PeptideCutter rules:
  pepsin(pH1.3):     cut at N-terminus of F/Y/W/L (these residues start fragments)
  trypsin:            cut after C-terminus of K/R (break after K,R)
  chymotrypsin(spec): cut after C-terminus of F/Y/W (L not cut; pause before Pro)
Merge the three enzyme cut points -> fragments; keep unique peptides of 2..20 aa
that contain no X.
"""
import os
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # scripts/ -> repo root
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
print(f"proteins={n}  total residues={tot:,}  mean={tot/n:.0f}aa")

pep=set("FYWL")   # N-term pre-cut
tryp=set("KR")    # C-term post-cut
chym=set("FYW")   # C-term post-cut (excl. L)

def cut_points(seq):
    pts={0,len(seq)}
    for i,a in enumerate(seq):
        if a in pep:  pts.add(i)      # pre-cut before a -> a is start
        if a in tryp: pts.add(i+1)    # post-cut after a -> i+1 is start
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
print(f"\nunique peptides (all, no X) = {len(ALL)}")
print(f"  2..6 aa  = {p26}")
print(f"  2..20 aa = {p220}")
print("length distribution (2..12):",{L:lens.get(L,0) for L in range(2,13)})

print("\n=== comparison ===")
print(f"{'metric':<22}{'moso (strict)':>14}{'garlic template':>14}")
print(f"{'proteins':<22}{n:>14}{113:>14}")
print(f"{'unique peptides all':<22}{len(ALL):>14}{5672:>14}")
print(f"{'unique 2..6aa':<22}{p26:>14}{1442:>14}")
print(f"{'PR>0.5 candidates (est 40%)':<22}{int(p26*0.4):>14}{249:>14}")

_out = os.path.join(DATA, "moso_253_peptides_strict.txt")
with open(_out,"w") as f:
    for p in sorted(ALL,key=lambda x:(len(x),x)): f.write(p+"\n")
print(f"\n-> {_out} ({len(ALL)} entries)")
