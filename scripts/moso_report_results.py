# -*- coding: utf-8 -*-
"""解析 Vina 对接结果，输出排名报告"""
import os, re

# 读结果
if not os.path.exists("moso_dock_results.tsv"):
    print("结果文件不存在，检查后台任务是否完成。")
    print("查看进度: wc -l moso_dock_results.tsv && ls moso_ligands/dock_*.pdbqt | wc -l")
    exit()

rows = []
skipped = 0
with open("moso_dock_results.tsv") as f:
    header = f.readline()  # skip header
    for line in f:
        line = line.strip()
        if not line or line.startswith("peptide"): continue
        parts = line.split()
        if len(parts) < 2:
            skipped += 1
            continue
        try:
            pep = parts[0]
            dg = float(parts[1])
            rows.append((pep, dg))
        except:
            skipped += 1

print(f"\n====== 毛竹 DPP4 对接排名 ======")
print(f"成功对接: {len(rows)}/{len(rows)+skipped}" + (f" (跳过的: {skipped})" if skipped else ""))
print(f"{'排名':>4} {'肽':>8} {'结合能(dG)':>12} {'分类':>10}")
print("-"*40)

# 按 dG 排序（最负=最好）
ranked = sorted(rows, key=lambda x: x[1])

# 结合标签
def label(dg):
    if dg <= -8.0: return "强结合"
    elif dg <= -6.5: return "中-强"
    elif dg <= -5.5: return "中等"
    elif dg <= -4.5: return "弱-中"
    else: return "弱结合"

for i, (pep, dg) in enumerate(ranked):
    tag = label(dg)
    print(f"{i+1:>4} {pep:>8} {dg:>10.3f}  {tag}")

# Top 候选
print(f"\n===== Top 10 候选肽 (推荐优先验证) =====")
top10 = ranked[:10]
print(f"{'序':>3} {'肽':>8} {'结合能':>10} {'特性':>15}")
for i, (pep, dg) in enumerate(top10):
    print(f"{i+1:>3} {pep:>8} {dg:>8.3f}  {label(dg)}")

# 聚类分析
cats = {"强结合":0,"中-强":0,"中等":0,"弱-中":0,"弱结合":0}
for _, dg in rows:
    cats[label(dg)] += 1
print(f"\n===== 结合能力分布 =====")
for k,v in sorted(cats.items()):
    bar = "#" * (v//2)
    print(f"  {k}: {v:>3}  {bar}")

# 保存报告
with open("moso_dock_ranking.txt", "w") as f:
    f.write("peptide\tdG_kcal_mol\tcategory\trmsd_best\n")
    for i, (pep, dg) in enumerate(ranked):
        # Try to get RMSD from output file
        rmsd = ""
        fpath = f"moso_ligands/dock_{i:02d}_{pep}.pdbqt"
        # Find correct file
        import glob
        matches = glob.glob(f"moso_ligands/dock_*_{pep}.pdbqt")
        if matches:
            with open(matches[0]) as df:
                for dl in df:
                    if dl.startswith("   1 "):
                        rmsd = dl.split()[2] if len(dl.split())>2 else ""
                        break
        f.write(f"{pep}\t{dg:.3f}\t{label(dg)}\t{rmsd}\n")

print(f"\n详细排名 -> moso_dock_ranking.txt")
print(f"对接输出 -> moso_ligands/dock_*.pdbqt")
print(f"\n下一步建议:")
print(f"  • dG < -7.0 kcal/mol 的候选 → Gly-Pro-pNA 体外 DPP4 抑制验证")
print(f"  • 最佳结合肽 → MD 模拟 + MM/PBSA (GROMACS, 50-150ns)")
print(f"  • 网络药理学 SwissTargetPrediction + STRING + DAVID")
print(f"  • 细胞活性 Caco-2 原位 DPP4 抑制 (复刻模板论文)")
