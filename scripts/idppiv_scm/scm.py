# -*- coding: utf-8 -*-
"""
iDPPIV-SCM 本地复现版 (Scoring Card Method, SCM)
=====================================================
完全离线、零外部依赖。复现 Charoenkwan et al. 2020
(J. Proteome Res. 19:4125-4136, DOI:10.1021/acs.jproteome.0c00590)
提出的 iDPPIV-SCM 评分卡方法。

SCM 核心:
  对位置 i (0-based) 与氨基酸 a, 倾向性得分
      P(i,a) = log2( Obs(i,a) / Exp(i,a) )
  其中 (在正/负训练集上统计):
      Obs(i,a) = Npos(i,a) / Npos(i)       # 正集中 a 出现在位置 i 的频率
      Exp(i,a) = (Npos(i,a)+Nneg(i,a)) / (Npos(i)+Nneg(i))  # 期望频率
  肽序列 S=(a1..aL) 的 SCM 总分:
      Score(S) = Σ_{i=0}^{L-1} P(i, a_{i+1})
  判定: Score > 0 预测为 DPP-IV 抑制肽(阳性)。
  (阈值 0 来自 log2(1)=0; 仅当 Obs>Exp 时单点贡献为正,
   与 SCM 原始论文定义一致。)

变长处理: 长度 >= i+1 的肽才贡献第 i 项; 短肽自然只累加其有效位置。
"""
import math
from collections import defaultdict

AA = "ACDEFGHIKLMNPQRSTVWY"
AA_SET = set(AA)


def load_tsv(path):
    """读取 index\tlabel\tsequence 格式, 返回 (seqs, labels)"""
    seqs, labels = [], []
    with open(path, encoding="utf-8") as f:
        header = f.readline()  # 跳过表头
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                # 退化: 仅 sequence
                parts = (None, None, line)
            _, label, seq = parts[0], parts[1], parts[2].strip()
            seq = "".join(c for c in seq.upper() if c in AA_SET)
            if not seq:
                continue
            seqs.append(seq)
            labels.append(int(label))
    return seqs, labels


def build_scm(pos_seqs, neg_seqs, L_max=None):
    """从正/负序列集合计算 SCM 评分卡 P[i][a]"""
    if L_max is None:
        L_max = max([len(s) for s in pos_seqs + neg_seqs] + [1])

    Npos = defaultdict(lambda: defaultdict(int))   # Npos[i][a]
    Nneg = defaultdict(lambda: defaultdict(int))
    Npos_i = defaultdict(int)                       # 位置 i 上正集肽数(长度>=i+1)
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
                P[i][a] = 0.0  # 未观测到 → 中性
    return P, L_max


def score(seq, P, L_max):
    """SCM 总分 (log2 倾向性之和)"""
    s = 0.0
    for i, a in enumerate(seq):
        if i >= L_max:
            break
        s += P[i].get(a, 0.0)
    return s


def predict(seq, P, L_max, threshold=0.0):
    return 1 if score(seq, P, L_max) > threshold else 0


def evaluate(seqs, labels, P, L_max, threshold=0.0):
    """返回 (ACC, MCC, 预测列表)"""
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
