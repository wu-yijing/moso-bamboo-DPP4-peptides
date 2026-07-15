# -*- coding: utf-8 -*-
"""Parse Vina docking results and emit a ranking report"""
import os, re

# read results
if not os.path.exists("moso_dock_results.tsv"):
    print("Result file not found; check whether the background task finished.")
    print("Check progress: wc -l moso_dock_results.tsv && ls moso_ligands/dock_*.pdbqt | wc -l")
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

print(f"\n====== Moso-bamboo DPP4 docking ranking ======")
print(f"successful docking: {len(rows)}/{len(rows)+skipped}" + (f" (skipped: {skipped})" if skipped else ""))
print(f"{'rank':>4} {'peptide':>8} {'bind energy(dG)':>12} {'class':>10}")
print("-"*40)

# sort by dG (most negative = best)
ranked = sorted(rows, key=lambda x: x[1])

# binding label
def label(dg):
    if dg <= -8.0: return "strong"
    elif dg <= -6.5: return "med-strong"
    elif dg <= -5.5: return "medium"
    elif dg <= -4.5: return "weak-med"
    else: return "weak"

for i, (pep, dg) in enumerate(ranked):
    tag = label(dg)
    print(f"{i+1:>4} {pep:>8} {dg:>10.3f}  {tag}")

# Top candidates
print(f"\n===== Top 10 candidate peptides (recommended for priority validation) =====")
top10 = ranked[:10]
print(f"{'no.':>3} {'peptide':>8} {'bind energy':>10} {'property':>15}")
for i, (pep, dg) in enumerate(top10):
    print(f"{i+1:>3} {pep:>8} {dg:>8.3f}  {label(dg)}")

# clustering
cats = {"strong":0,"med-strong":0,"medium":0,"weak-med":0,"weak":0}
for _, dg in rows:
    cats[label(dg)] += 1
print(f"\n===== binding-ability distribution =====")
for k,v in sorted(cats.items()):
    bar = "#" * (v//2)
    print(f"  {k}: {v:>3}  {bar}")

# save report
with open("moso_dock_ranking.txt", "w") as f:
    f.write("peptide\tdG_kcal_mol\tcategory\trmsd_best\n")
    for i, (pep, dg) in enumerate(ranked):
        # try to get RMSD from output file
        rmsd = ""
        fpath = f"moso_ligands/dock_{i:02d}_{pep}.pdbqt"
        # find correct file
        import glob
        matches = glob.glob(f"moso_ligands/dock_*_{pep}.pdbqt")
        if matches:
            with open(matches[0]) as df:
                for dl in df:
                    if dl.startswith("   1 "):
                        rmsd = dl.split()[2] if len(dl.split())>2 else ""
                        break
        f.write(f"{pep}\t{dg:.3f}\t{label(dg)}\t{rmsd}\n")

print(f"\ndetailed ranking -> moso_dock_ranking.txt")
print(f"docking output -> moso_ligands/dock_*.pdbqt")
print(f"\nNext-step suggestions:")
print(f"  - candidates with dG < -7.0 kcal/mol -> in-vitro DPP4 inhibition assay with Gly-Pro-pNA")
print(f"  - best-binding peptide -> MD simulation + MM/PBSA (GROMACS, 50-150ns)")
print(f"  - network pharmacology: SwissTargetPrediction + STRING + DAVID")
print(f"  - cellular activity: Caco-2 in-situ DPP4 inhibition (reproduce template paper)")
