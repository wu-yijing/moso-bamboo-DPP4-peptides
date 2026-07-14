# -*- coding: utf-8 -*-
"""iDPPIV-SCM 配方探测: 测试多种 SCM 变体, 找能复现~0.80的配方"""
import sys, os, math
from collections import defaultdict
sys.path.insert(0, os.path.dirname(__file__))
from scm import load_tsv

AA = "ACDEFGHIKLMNPQRSTVWY"
AAS = set(AA)

def load(path):
    seqs, labs = [], []
    first = True
    for line in open(path, encoding="utf-8"):
        if first: first = False; continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 3: continue
        s = "".join(c for c in p[2].upper() if c in AAS)
        if s: seqs.append(s); labs.append(int(p[1]))
    return seqs, labs

trS, trL = load("idppiv_scm/data/train.tsv")
teS, teL = load("idppiv_scm/data/test.tsv")
pos = [s for s,l in zip(trS,trL) if l==1]
neg = [s for s,l in zip(trS,trL) if l==0]

def build_mono(pos, neg, Lmax):
    Np=defaultdict(lambda:defaultdict(int)); Nn=defaultdict(lambda:defaultdict(int))
    Npi=defaultdict(int); Nni=defaultdict(int)
    for s in pos:
        for i,a in enumerate(s): Np[i][a]+=1; Npi[i]+=1
    for s in neg:
        for i,a in enumerate(s): Nn[i][a]+=1; Nni[i]+=1
    P={}
    for i in range(Lmax):
        P[i]={}
        denom=Npi[i]+Nni[i]
        for a in AA:
            obs=Np[i][a]/Npi[i] if Npi[i] else 0
            exp=(Np[i][a]+Nn[i][a])/denom if denom else 0
            P[i][a]=math.log2(obs/exp) if obs>0 and exp>0 else 0.0
    return P, Lmax

def build_mono_global(pos, neg):
    """全局组成型: P(a) 不依赖位置"""
    Np=defaultdict(int); Nn=defaultdict(int)
    Tp=sum(len(s) for s in pos); Tn=sum(len(s) for s in neg)
    for s in pos:
        for a in s: Np[a]+=1
    for s in neg:
        for a in s: Nn[a]+=1
    P={}
    for a in AA:
        obs=Np[a]/Tp if Tp else 0
        exp=((Np[a]+Nn[a])/(Tp+Tn)) if (Tp+Tn) else 0
        P[a]=math.log2(obs/exp) if obs>0 and exp>0 else 0.0
    return P

def build_di(pos, neg, Lmax):
    """二肽位置特异性"""
    Np=defaultdict(lambda:defaultdict(int)); Nn=defaultdict(lambda:defaultdict(int))
    Npi=defaultdict(int); Nni=defaultdict(int)
    for s in pos:
        for i in range(len(s)-1):
            d=s[i:i+2]; Np[i][d]+=1; Npi[i]+=1
    for s in neg:
        for i in range(len(s)-1):
            d=s[i:i+2]; Nn[i][d]+=1; Nni[i]+=1
    P={}
    for i in range(Lmax):
        P[i]={}
        denom=Npi[i]+Nni[i]
        # 遍历出现过的二肽
        keys=set(Np[i]) | set(Nn[i])
        for d in keys:
            obs=Np[i][d]/Npi[i] if Npi[i] else 0
            exp=(Np[i][d]+Nn[i][d])/denom if denom else 0
            P[i][d]=math.log2(obs/exp) if obs>0 and exp>0 else 0.0
    return P, Lmax

Lmax = max(len(s) for s in pos+neg)

def score_mono(seq, P, Lmax, globalP=None):
    if globalP:
        return sum(globalP.get(a,0.0) for a in seq)
    return sum(P[i].get(a,0.0) for i,a in enumerate(seq) if i<Lmax)

def score_di(seq, P, Lmax):
    return sum(P[i].get(seq[i:i+2],0.0) for i in range(len(seq)-1) if i<Lmax)

def best_thresh(seqs, labs, scores):
    """在训练集上扫描阈值, 取 MCC 最大"""
    cand=sorted(set(scores))
    best=(-1,0,None)
    for t in [min(scores)-1e-9]+[(cand[i]+cand[i+1])/2 for i in range(len(cand)-1)]+[max(scores)+1e-9]:
        tp=fp=tn=fn=0
        for sc,y in zip(scores,labs):
            p=1 if sc>t else 0
            if y==1 and p==1: tp+=1
            elif y==0 and p==1: fp+=1
            elif y==0 and p==0: tn+=1
            else: fn+=1
        den=math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
        mcc=((tp*tn)-(fp*fn))/den if den else 0
        if mcc>best[0]: best=(mcc,t,sc)
    return best[1]

def acc_at(seqs,labs,scores,th):
    ok=sum(1 for sc,y in zip(scores,labs) if (1 if sc>th else 0)==y)
    return ok/len(labs)

# ---- 变体1: 位置特异性单残基 ----
Pm,Lm=build_mono(pos,neg,Lmax)
sm=[score_mono(s,Pm,Lm) for s in trS]; se=[score_mono(s,Pm,Lm) for s in teS]
tm=best_thresh(trS,trL,sm); print(f"[单残基-位置] 最优阈值={tm:.3f}  trainACC={acc_at(trS,trL,sm,tm):.3f}  testACC={acc_at(teS,teL,se,tm):.3f}")

# ---- 变体2: 全局组成单残基 ----
Pmg=build_mono_global(pos,neg)
smg=[score_mono(s,None,Lmax,Pmg) for s in trS]; seg=[score_mono(s,None,Lmax,Pmg) for s in teS]
tmg=best_thresh(trS,trL,smg); print(f"[单残基-全局] 最优阈值={tmg:.3f}  trainACC={acc_at(trS,trL,smg,tmg):.3f}  testACC={acc_at(teS,teL,seg,tmg):.3f}")

# ---- 变体3: 二肽位置特异性 ----
Pd,Ld=build_di(pos,neg,Lmax)
sd=[score_di(s,Pd,Ld) for s in trS]; de=[score_di(s,Pd,Ld) for s in teS]
td=best_thresh(trS,trL,sd); print(f"[二肽-位置]   最优阈值={td:.3f}  trainACC={acc_at(trS,trL,sd,td):.3f}  testACC={acc_at(teS,teL,de,td):.3f}")

# ---- 变体4: 混合 (单残基位置 + 二肽位置) ----
sh=[sm[i]+sd[i] for i in range(len(sm))]
dh=[se[i]+de[i] for i in range(len(se))]
th=best_thresh(trS,trL,sh); print(f"[混合 单+二]  最优阈值={th:.3f}  trainACC={acc_at(trS,trL,sh,th):.3f}  testACC={acc_at(teS,teL,dh,th):.3f}")

# ---- 变体5: 混合 (全局单残基 + 二肽位置) ----
sh2=[smg[i]+sd[i] for i in range(len(smg))]
dh2=[seg[i]+de[i] for i in range(len(seg))]
th2=best_thresh(trS,trL,sh2); print(f"[混合 全局单+二] 最优阈值={th2:.3f}  trainACC={acc_at(trS,trL,sh2,th2):.3f}  testACC={acc_at(teS,teL,dh2,th2):.3f}")
