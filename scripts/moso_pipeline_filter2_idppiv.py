# -*- coding: utf-8 -*-
"""
毛竹 DPP4 候选肽 二次收窄 (DPP4 抑制肽已知结构偏好) —— iDPPIV-SCM 版
==========================================================================
基于文献: N端疏水(I/L/V/F/W/M); 长度2-5aa(本规则取3-5)最佳; 富疏水;
        DPP4 S1 口袋对 Pro 特异性(第2位 Pro/Ala 最优)。
在 iDPPIV-SCM 候选(已含 DPP4 抑制倾向性)基础上, 叠加 DPP4 偏好规则 -> 可对接量级。
"""
import os
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # scripts/ -> 仓库根
DATA = os.path.join(REPO, "data")

HYDRO = set("ILVFMWACY")          # 疏水

def dpp4_rules(seq, score):
    if len(seq) < 3 or len(seq) > 5: return False   # 3-5aa 为 DPP4 抑制肽典型长度
    if seq[0] not in HYDRO: return False              # N端疏水
    if score <= 0.0: return False                   # iDPPIV 净正向(比随机更抑制样)
    if not any(a in HYDRO for a in seq): return False
    return True

# 读 iDPPIV 候选文件 (首行为表头, 跳过)
rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(DATA, "moso_candidates_idppiv.txt")) if l.strip()]
if rows and rows[0][0].lower().startswith("peptide"):
    rows = rows[1:]
cands = [(p, float(s)) for p, s, *_ in rows]
narrow = [(p, s) for p, s in cands if dpp4_rules(p, s)]
print(f"iDPPIV 候选 {len(cands)} -> DPP4偏好收窄 -> {len(narrow)} 条")

# 二级优先: 第2位 Pro/Ala 标记 (DPP4 S1 口袋最优)
def p2(p): return len(p) >= 2 and p[1] in "PA"
tier1 = [(p, s) for p, s in narrow if p2(p)]
tier2 = [(p, s) for p, s in narrow if not p2(p)]
print(f"  其中 第2位P/A(最优): {len(tier1)} | 其余: {len(tier2)}")

# 输出: 对接优先集(tier1 在前, 均按 iDPPIV 分降序)
n_top = min(60, len(narrow))
top = (sorted(tier1, key=lambda x: -x[1]) + sorted(tier2, key=lambda x: -x[1]))[:n_top]
_queue = os.path.join(DATA, "moso_dock_queue_idppiv.txt")
with open(_queue, "w") as f:
    for p, s in top:
        f.write(f"{p}\t{s:.3f}\t{'P2' if p2(p) else ''}\n")
print(f"\n对接队列(优先前{n_top}) -> {_queue}")
print("样例(前15):")
for p, s in top[:15]:
    print(f"  {p:6s} iDPPIV={s:+.3f} {'[P2]' if p2(p) else ''}")

# ---- 与旧队列对照: 最佳对接肽 LPPQ 是否仍入选, 重叠多少 ----
old = set(l.split("\t")[0].strip() for l in open(os.path.join(DATA, "moso_dock_queue.txt")) if l.strip())
new = set(p for p, _ in top)
print(f"\n旧对接队列({len(old)}) vs 新对接队列({len(new)}): 重叠 {len(old & new)} 条")
for key in ["LPPQ", "APSPE", "LAPSP", "LPGP"]:
    print(f"  旧最佳 {key}: {'在新队列中 ✓' if key in new else '不在新队列 ✗'}")
