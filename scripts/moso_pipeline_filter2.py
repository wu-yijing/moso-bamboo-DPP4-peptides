# -*- coding: utf-8 -*-
"""毛竹 DPP4 候选肽 二次收窄 (DPP4 抑制肽已知结构偏好)
基于文献: N端疏水(I/L/V/F/W/M); 长度2-5aa最佳; 第2位Pro/Ala优先; 富疏水。
在 4289 条 PeptideRanker>0.5 且非过敏/非毒 的基础上, 叠加 DPP4 偏好规则 -> 可对接量级。
"""
HYDRO = set("ILVFMWACY")          # 疏水
def dpp4_rules(seq, score):
    if len(seq) < 3 or len(seq) > 5: return False   # 3-5aa 为 DPP4 抑制肽典型长度
    if seq[0] not in HYDRO: return False          # N端疏水
    if score < 0.55: return False                # 略提阈值
    if not any(a in HYDRO for a in seq): return False
    return True

rows = [l.rstrip("\n").split("\t") for l in open("E:/workbuddy/Claw/moso_candidates_pr_filtered.txt") if l.strip()]
cands = [(p, float(s)) for p, s in rows]
narrow = [(p, s) for p, s in cands if dpp4_rules(p, s)]
print(f"4289 候选 -> DPP4偏好收窄 -> {len(narrow)} 条")

# 二级优先: 第2位 Pro/Ala 标记
def p2(p): return len(p) >= 2 and p[1] in "PA"
tier1 = [(p, s) for p, s in narrow if p2(p)]      # 第2位P/A (最优)
tier2 = [(p, s) for p, s in narrow if not p2(p)]
print(f"  其中 第2位P/A(最优): {len(tier1)} | 其余: {len(tier2)}")

# 输出: 对接优先集(tier1 在前)
n_top = min(60, len(narrow))
top = (sorted(tier1, key=lambda x:-x[1]) + sorted(tier2, key=lambda x:-x[1]))[:n_top]
with open("E:/workbuddy/Claw/moso_dock_queue.txt", "w") as f:
    for p, s in top:
        f.write(f"{p}\t{s:.3f}\t{'P2' if p2(p) else ''}\n")
print(f"\n对接队列(优先前{n_top}) -> E:/workbuddy/Claw/moso_dock_queue.txt")
print("样例(前15):")
for p, s in top[:15]:
    print(f"  {p:6s} score={s:.3f} {'[P2]' if p2(p) else ''}")

print("\n=== 完整漏斗 (模板大蒜 1442->249->34合成) ===")
print(f"毛竹: 7988 -> 4333 -> 4333 -> 4289 -> {len(narrow)} (DPP4偏好) -> 对接队列{len(top)}")
