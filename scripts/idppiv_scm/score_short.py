# -*- coding: utf-8 -*-
"""
对毛竹 DPP4 项目【全部 4,950 条 2-6 aa 短肽】用本地复现的 iDPPIV-SCM 重新打分,
替换 PeptideRanker-style 代理分 (原 moso_candidates_pr_filtered.txt 全为 1.000)。

输出:
  moso_candidates_idppiv_short.tsv   —— 全部 4,950 条短肽的 iDPPIV-SCM 评分
  moso_candidates_idppiv_proxy.tsv   —— 旧代理候选集(4,289 条)换成 iDPPIV 分
"""
import os, sys
sys.path.insert(0, "E:/workbuddy/Claw")
from idppiv_scm.model import build_propensity, score, score_mean, predict, DEFAULT_THRESHOLD

AA = set("ACDEFGHIKLMNPQRSTVWY")
HERE = "E:/workbuddy/Claw"

# 1) 训练评分卡
P, stats = build_propensity()
Tp, Tn, npos, nneg = stats
print(f"[评分卡] 训练集 正={npos} 负={nneg}; 阈值(预测抑制)= {DEFAULT_THRESHOLD:+.3f}")

def load_short(path):
    out = []
    for l in open(path, encoding="utf-8"):
        s = l.strip()
        if 2 <= len(s) <= 6 and all(c in AA for c in s):
            out.append(s)
    return out

def score_file(in_path, out_path, label):
    seqs = load_short(in_path) if "strict" in in_path else [
        l.strip() for l in open(in_path, encoding="utf-8") if l.strip()
    ]
    rows = []
    for s in seqs:
        sc = score(s, P)
        sm = score_mean(s, P)
        pr = predict(s, P)
        rows.append((s, len(s), sc, sm, pr))
    rows.sort(key=lambda r: -r[2])
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("peptide\tlength\tiDPPIV_score\tiDPPIV_mean\tpredicted_DPP4_inhibitory\n")
        for s, ln, sc, sm, pr in rows:
            f.write(f"{s}\t{ln}\t{sc:.4f}\t{sm:.4f}\t{pr}\n")
    npos = sum(1 for r in rows if r[4] == 1)
    print(f"[{label}] 共 {len(rows)} 条 -> 预测为 DPP-IV 抑制肽: {npos} 条 ({100*npos/len(rows):.1f}%)")
    print(f"        iDPPIV_score 范围: [{min(r[2] for r in rows):.3f}, {max(r[2] for r in rows):.3f}]")
    print(f"        文件: {out_path}")
    return rows

# 2a) 全部 4,950 条短肽
all_rows = score_file(
    os.path.join(HERE, "moso_253_peptides_strict.txt"),
    os.path.join(HERE, "moso_candidates_idppiv_short.tsv"),
    "全部短肽(2-6aa)",
)

# 2b) 旧代理候选集 (4,289) 换成 iDPPIV 分
proxy_path = os.path.join(HERE, "moso_candidates_pr_filtered.txt")
if os.path.exists(proxy_path):
    score_file(
        proxy_path,
        os.path.join(HERE, "moso_candidates_idppiv_proxy.tsv"),
        "旧代理候选集(原PeptideRanker)",
    )

# 3) 展示 Top-20
print("\n=== Top-20 (按 iDPPIV_score 降序) ===")
print(f"{'肽':<10}{'len':>4}{'iDPPIV_score':>15}{'iDPPIV_mean':>14}{'pred':>6}")
for s, ln, sc, sm, pr in all_rows[:20]:
    print(f"{s:<10}{ln:>4}{sc:>15.4f}{sm:>14.4f}{pr:>6}")

# 4) 对比旧代理分的"全 1.000" 与 iDPPIV 真实分布
import statistics
vals = [r[2] for r in all_rows]
print(f"\n=== 分布对比 ===")
print(f"旧 PeptideRanker 代理分: 恒定 = 1.000 (无区分力)")
print(f"新 iDPPIV-SCM 分:     mean={statistics.mean(vals):.3f}  median={statistics.median(vals):.3f}  stdev={statistics.pstdev(vals):.3f}")
