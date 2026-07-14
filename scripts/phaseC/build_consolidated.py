# -*- coding: utf-8 -*-
"""
生成 Top 候选 consolidated 表 (data/phaseC/top_candidates_consolidated.tsv)

合并两个来源，按肽序列 (peptide) 连接：
  1) docking/moso_dock_results_idppiv_clean.tsv
       -> iDPPIV_score (Phase A 离线活性预筛) + dG_best/dG_spread/n_dock (Vina 对接)
  2) data/phaseC/phaseC_peptides.tsv
       -> ADMET (MW/pI/charge/GRAVY/Boman/TPSA/HBD/HBA/QED) + GI 稳定性 + DPP4 底物基序

排序：dG_best 升序（结合最强优先，与稿件对接小节一致）；取 Top 20。
iDPPIV 已作为构建对接队列的预筛门槛，故此处以对接 dG 为主序。
纯计算研究，无任何湿实验验证。
"""
import os, csv

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DOCK = os.path.join(REPO, "docking", "moso_dock_results_idppiv_clean.tsv")
PHASEC = os.path.join(REPO, "data", "phaseC", "phaseC_peptides.tsv")
OUT = os.path.join(REPO, "data", "phaseC", "top_candidates_consolidated.tsv")
TOPN = 20


def read_tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def main():
    dock = {r["peptide"]: r for r in read_tsv(DOCK)}
    phasec = {r["peptide"]: r for r in read_tsv(PHASEC)}

    merged = []
    for pep, d in dock.items():
        if pep not in phasec:
            continue
        p = phasec[pep]
        merged.append({
            "peptide": pep,
            "length": p["length"],
            "iDPPIV_score": d["iDPPIV_score"],
            "dG_Vina": d["dG_best"],
            "dG_spread": d["dG_spread"],
            "n_dock": d["n_dock"],
            "MW": p["rdkit_mw"],
            "pI": p["pI"],
            "net_charge_7p4": p["net_charge_7p4"],
            "GRAVY": p["GRAVY"],
            "Boman": p["boman_index"],
            "TPSA": p["TPSA"],
            "HBD": p["HBD"],
            "HBA": p["HBA"],
            "QED": p["QED"],
            "GI_stability": p["GI_stability_class"],
            "DPP4_substrate_motif": p["DPP4_substrate_motif"],
            "is_Top3": p["is_Top3"],
            # 排序键
            "_dg": float(d["dG_best"]),
        })

    merged.sort(key=lambda x: x["_dg"])  # 升序 = 最负 = 结合最强
    top = merged[:TOPN]

    cols = ["peptide", "length", "iDPPIV_score", "dG_Vina", "dG_spread", "n_dock",
            "MW", "pI", "net_charge_7p4", "GRAVY", "Boman", "TPSA", "HBD", "HBA",
            "QED", "GI_stability", "DPP4_substrate_motif", "is_Top3"]
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in top:
            w.writerow(r)
    print(f"[OK] 写出 {len(top)} 条 -> {OUT}")
    # 控制台摘要
    print(f"{'pep':7s} {'iDPPIV':>7s} {'dG':>8s} {'MW':>7s} {'pI':>5s} {'Boman':>6s} {'QED':>6s}  GI")
    for r in top:
        print(f"{r['peptide']:7s} {float(r['iDPPIV_score']):7.3f} {float(r['dG_Vina']):8.3f} "
              f"{float(r['MW']):7.1f} {float(r['pI']):5.2f} {float(r['Boman']):6.2f} "
              f"{float(r['QED']):6.3f}  {r['GI_stability']}")


if __name__ == "__main__":
    main()
