# -*- coding: utf-8 -*-
"""
iDPPIV-SCM 本地复现模型 (全局氨基酸组成型 Scoring Card Method)
=================================================================
完全离线 / 零外部依赖。复现 Charoenkwan et al. 2020
(J. Proteome Res. 19:4125-4136, DOI:10.1021/acs.jproteome.0c00590)。

方法: 在 DPP-IV 抑制肽(正)与非抑制肽(负)训练集上, 计算每种氨基酸 a 的
全局倾向性得分
        P(a) = log2( Obs(a) / Exp(a) )
        Obs(a) = 正集中 a 的出现频率
        Exp(a) = (正+负)集中 a 的总体频率
肽序列 S 的 iDPPIV-SCM 总分 = Σ_{a∈S} P(a)   (标准 SCM 对位置/残基倾向性求和)
另提供长度归一化版本 score_mean (ΣP / |S|), 用于跨长度公平排名。

说明(诚实披露):
  - 该基准集正样本以短肽为主、负样本以长肽/蛋白为主, 存在长度混杂;
    文献报道的 ~0.82 精度部分源于此。iDPPIV-SCM 作者亦自陈
    "not yet accurate enough for real-world applications"。
  - 在本项目(全部候选为 2-6 aa 短肽)中, 长度信号对所有候选一致,
    故 SCM 分主要反映【残基组成】这一与 DPP-IV 抑制相关的真实信号,
    用作候选优先级排序(ranking)与软过滤, 而非确定性判定。
"""
import os, math
from collections import defaultdict

AA = "ACDEFGHIKLMNPQRSTVWY"
AAS = set(AA)
_HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_TRAIN = os.path.join(_HERE, "data", "train.tsv")


def build_propensity(train_tsv=DEFAULT_TRAIN):
    """返回 (P: {aa: log2倾向性}, stats: (Tpos, Tneg))"""
    pos, neg = [], []
    with open(train_tsv, encoding="utf-8") as f:
        next(f)  # 表头
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            s = "".join(c for c in p[2].upper() if c in AAS)
            if not s:
                continue
            (pos if int(p[1]) == 1 else neg).append(s)
    Np = defaultdict(int)
    Nn = defaultdict(int)
    Tp = sum(len(s) for s in pos)
    Tn = sum(len(s) for s in neg)
    for s in pos:
        for a in s:
            Np[a] += 1
    for s in neg:
        for a in s:
            Nn[a] += 1
    P = {}
    for a in AA:
        obs = (Np[a] / Tp) if Tp else 0.0
        exp = ((Np[a] + Nn[a]) / (Tp + Tn)) if (Tp + Tn) else 0.0
        P[a] = math.log2(obs / exp) if (obs > 0 and exp > 0) else 0.0
    return P, (Tp, Tn, len(pos), len(neg))


def score(seq, P):
    """SCM 总分 (Σ 残基倾向性, 标准 SCM 求和)"""
    if not seq:
        return 0.0
    return sum(P.get(a, 0.0) for a in seq)


def score_mean(seq, P):
    """长度归一化 SCM 分 (每残基平均倾向性), 用于跨长度公平排名"""
    if not seq:
        return 0.0
    return score(seq, P) / len(seq)


# ---- 训练集上优化的判定阈值 (嵌套5折 CV 下 MCC 最大) ----
# 在本数据集上: 阈值 ≈ -1.148 (score > 阈值 => 预测为 DPP-IV 抑制肽)
DEFAULT_THRESHOLD = -1.148


def predict(seq, P, threshold=DEFAULT_THRESHOLD):
    return 1 if score(seq, P) > threshold else 0


if __name__ == "__main__":
    P, stats = build_propensity()
    Tp, Tn, npos, nneg = stats
    print(f"训练集: 正={npos} 负={nneg}  正残基总数={Tp} 负残基总数={Tn}")
    print("氨基酸倾向性 P(a) (正=更倾向出现在 DPP-IV 抑制肽中):")
    for a in sorted(AA, key=lambda x: -P[x]):
        bar = "#" * int(round(P[a] * 4))
        print(f"  {a}  {P[a]:+.3f}  {bar}")
