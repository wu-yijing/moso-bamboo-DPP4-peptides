# -*- coding: utf-8 -*-
"""
公平跨队列对比: 新 iDPPIV 队列 vs 旧代理队列 (两者均用同一 RDKit 管线对接 1WCY)
------------------------------------------------------------------------
读取:
  moso_dock_results_idppiv_clean.tsv  (新队列 60, 已去重取最优)
  moso_dock_results_old_rdkit.tsv    (旧队列 60, 同法重对接, 去重取最优)
输出: 控制台对比 + moso_dock_compare.tsv
"""
import collections, statistics, os

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # scripts/ -> 仓库根
DOCK = os.path.join(REPO, "docking")

def load_best(path, has_score=True):
    agg = collections.defaultdict(list)
    for l in open(path, encoding="utf-8"):
        l = l.rstrip("\n")
        if not l or l.startswith("peptide"):
            continue
        parts = l.split("\t")
        seq = parts[0]
        sc = float(parts[1]) if len(parts) > 1 else 0.0
        if len(parts) < 3 or parts[2] in ("NA", "ERR"):
            continue
        agg[seq].append((sc, float(parts[2])))
    out = {}
    for seq, vals in agg.items():
        out[seq] = (vals[0][0], min(v[1] for v in vals))  # (score, best_dG)
    return out

new = load_best(os.path.join(DOCK, "moso_dock_results_idppiv_clean.tsv"))
old = load_best(os.path.join(DOCK, "moso_dock_results_old_rdkit.tsv"))

print(f"新队列唯一肽: {len(new)}   旧队列唯一肽: {len(old)}")

new_dg = [v[1] for v in new.values()]
old_dg = [v[1] for v in old.values()]

def stats(name, xs):
    xs_sorted = sorted(xs)
    print(f"\n[{name}] n={len(xs)}")
    print(f"  最佳 dG (最负): {min(xs):.3f}")
    print(f"  中位 dG:        {statistics.median(xs):.3f}")
    print(f"  均值 dG:        {statistics.mean(xs):.3f}")
    print(f"  dG <= -6.5 占比: {sum(1 for x in xs if x<=-6.5)/len(xs)*100:.1f}%")
    print(f"  dG <= -6.0 占比: {sum(1 for x in xs if x<=-6.0)/len(xs)*100:.1f}%")

stats("新 iDPPIV 队列 (RDKit 同法)", new_dg)
stats("旧 代理队列 (RDKit 同法)", old_dg)

# 重叠肽(两队列共有)直接配对比较
overlap = set(new) & set(old)
print(f"\n=== 重叠肽配对比较 (n={len(overlap)}) ===")
better_new = sum(1 for p in overlap if new[p][1] < old[p][1])   # 更负=更好
better_old = sum(1 for p in overlap if old[p][1] < new[p][1])
print(f"  新队列 dG 更优(更负)的肽: {better_new}")
print(f"  旧队列 dG 更优(更负)的肽: {better_old}")

new_best = min(new_dg); old_best = min(old_dg)
print(f"\n=== 结论性对比 ===")
print(f"  新队列最佳结合肽 dG = {new_best:.3f} ({min(new, key=lambda p: new[p][1])})")
print(f"  旧队列最佳结合肽 dG = {old_best:.3f} ({min(old, key=lambda p: old[p][1])})")
if new_best < old_best:
    print(f"  >>> 同法制备下, 新 iDPPIV 优先化给出更优结合肽 (Δ={new_best-old_best:+.3f} kcal/mol)")
else:
    print(f"  >>> 同法制备下, 新队列最佳未超越旧队列 (Δ={new_best-old_best:+.3f} kcal/mol)")
    print(f"      即 iDPPIV 优先化在本受体/盒子上未带来结合更强肽; 其价值在活性筛选维度(见分析)。")

# 保存
with open(os.path.join(DOCK, "moso_dock_compare.tsv"), "w", encoding="utf-8") as f:
    f.write("metric\tnew_idppiv\told_proxy\n")
    f.write(f"n\t{len(new)}\t{len(old)}\n")
    f.write(f"best_dG\t{new_best:.3f}\t{old_best:.3f}\n")
    f.write(f"median_dG\t{statistics.median(new_dg):.3f}\t{statistics.median(old_dg):.3f}\n")
    f.write(f"mean_dG\t{statistics.mean(new_dg):.3f}\t{statistics.mean(old_dg):.3f}\n")
    f.write(f"frac<=-6.5\t{sum(1 for x in new_dg if x<=-6.5)/len(new_dg):.3f}\t{sum(1 for x in old_dg if x<=-6.5)/len(old_dg):.3f}\n")
    f.write(f"overlap_n\t{len(overlap)}\t\n")
    f.write(f"overlap_better_new\t{better_new}\t\n")
    f.write(f"overlap_better_old\t{better_old}\t\n")
print(f"\n对比已写 -> {os.path.join(DOCK, 'moso_dock_compare.tsv')}")
