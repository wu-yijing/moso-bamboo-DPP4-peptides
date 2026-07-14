# -*- coding: utf-8 -*-
"""
毛竹 DPP4 抑制肽 全流程之 过滤阶段（iDPPIV-SCM 离线复现版）
====================================================================
替换原 PeptideRanker-style 代理为 文献验证、完全离线可复现的
iDPPIV-SCM (全局氨基酸组成型 Scoring Card Method) 评分。

阶段1 iDPPIV-SCM DPP-IV 抑制倾向性评分 (连续分, 越高越像 DPP-IV 抑制肽;
        并以训练集优化阈值做"预测为抑制肽"的软过滤)
阶段2 AllerTOP-style 过敏原预测 (剔除高风险)  [仍代理]
阶段3 ToxinPred-style 毒性预测 (剔除高风险)    [仍代理]

输出: 各阶段漏斗 + 终选候选清单(含 iDPPIV 总分与长度归一化分, 按分降序)
说明: 阶段1 已为 DPP-IV 专用、文献验证、离线可复现; 阶段2/3 仍依赖
      官方网页工具, 正式稿件须以官方工具/服务器输出为准。
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "idppiv_scm"))
from model import build_propensity, score, score_mean, predict, DEFAULT_THRESHOLD

# ---- 理化参数(阶段2/3 代理, 沿用原版) ----
KD = {'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,
       'H':-3.2,'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,
       'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
MW = {'A':89,'R':174,'N':132,'D':133,'C':121,'Q':146,'E':147,'G':75,'H':155,
       'I':131,'L':131,'K':146,'M':149,'F':165,'P':115,'S':105,'T':119,'W':204,'Y':181,'V':117}

def aller_top(seq):
    rk = (seq.count('R')+seq.count('K'))/len(seq)
    n  = seq.count('N')/len(seq)
    return (rk*0.6 + n*0.4) > 0.45   # True=高风险(剔除)

def toxin_pred(seq):
    if seq.count('C') >= 2: return True
    if (seq.count('R')+seq.count('K')) >= 3: return True
    return False

# ---- 读入 ----
peps = [l.strip() for l in open("E:/workbuddy/Claw/moso_253_peptides_strict.txt") if l.strip()]
print(f"输入唯一肽 = {len(peps)}")

# ---- 构建 iDPPIV-SCM 评分卡 (基于公开数据集, 离线) ----
P, stats = build_propensity()
print(f"[iDPPIV-SCM] 评分卡已构建 (训练正={stats[2]}, 负={stats[3]}; 阈值={DEFAULT_THRESHOLD})")

# ---- 阶段1 iDPPIV-SCM ----
scored = [(p, score(p, P), score_mean(p, P)) for p in peps]
s1 = [(p, sc, mn) for p, sc, mn in scored if predict(p, P) == 1]
print(f"阶段1 iDPPIV-SCM 预测为抑制肽 : {len(s1)}")

# ---- 阶段2 AllerTOP (代理) ----
s2 = [(p, sc, mn) for p, sc, mn in s1 if not aller_top(p)]
print(f"阶段2 去过敏原(AllerTOP): {len(s2)}")

# ---- 阶段3 ToxinPred (代理) ----
s3 = [(p, sc, mn) for p, sc, mn in s2 if not toxin_pred(p)]
print(f"阶段3 去毒性(ToxinPred) : {len(s3)}")

# ---- 校验: 模板已知活性肽应被保留 ----
refs = ["WPHY","WPQY","VAPGW","WPH","WPQ","VAP"]
print("\n校验(模板活性肽片段, 若存在于池中):")
for r in refs:
    hit = [p for p,_,_ in s3 if r in p]
    print(f"  {r}: 命中 {len(hit)} 条" + (f" 例 {hit[0]}" if hit else " (未在池中)"))

# ---- 输出 (按 iDPPIV 总分降序) ----
s3_sorted = sorted(s3, key=lambda x: -x[1])
with open("E:/workbuddy/Claw/moso_candidates_idppiv.txt","w") as f:
    f.write("peptide\tiDPPIV_score\tiDPPIV_mean\tpredicted_DPP4_inhibitory\n")
    for p, sc, mn in s3_sorted:
        f.write(f"{p}\t{sc:.3f}\t{mn:.3f}\t1\n")
print(f"\n终选候选 -> E:/workbuddy/Claw/moso_candidates_idppiv.txt ({len(s3)} 条, 按 iDPPIV 分降序)")

# ---- 对接队列 Top-60 (按 iDPPIV 总分) ----
top60 = s3_sorted[:60]
with open("E:/workbuddy/Claw/moso_dock_queue_idppiv.txt","w") as f:
    for p, sc, mn in top60:
        f.write(f"{p}\t{sc:.3f}\n")
print(f"对接队列 Top-60 -> E:/workbuddy/Claw/moso_dock_queue_idppiv.txt")

# ---- 漏斗对照 ----
print("\n=== 漏斗对照 ===")
print(f"旧(PeptideRanker代理): 4950 -> 2019 -> ... -> 60(对接)")
print(f"新(iDPPIV-SCM)     : {len(peps)} -> {len(s1)} -> {len(s2)} -> {len(s3)}  (Top-60 进对接)")
