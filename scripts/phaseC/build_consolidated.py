# -*- coding: utf-8 -*-
"""
Generate a Top-candidate consolidated table (data/phaseC/top_candidates_consolidated.tsv)

Merge two sources and join by peptide sequence (peptide):
  1) docking/moso_dock_results_idppiv_clean.tsv
       -> iDPPIV_score (Phase A offline activity pre-screen) + dG_best/dG_spread/n_dock (Vina docking)
  2) data/phaseC/phaseC_peptides.tsv
       -> ADMET (MW/pI/charge/GRAVY/Boman/TPSA/HBD/HBA/QED) + GI stability + DPP4 substrate motif
Order by dG_best ascending (strongest binding first, consistent with the docking
subsection of the manuscript); take Top 20.
iDPPIV already serves as the pre-screen threshold when building the docking queue,
so here docking dG is the primary sort key.
Pure computational study; no wet-lab validation of any kind.
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
            # sort key
            "_dg": float(d["dG_best"]),
        })

    merged.sort(key=lambda x: x["_dg"])  # ascending = most negative = strongest binding
    top = merged[:TOPN]

    cols = ["peptide", "length", "iDPPIV_score", "dG_Vina", "dG_spread", "n_dock",
            "MW", "pI", "net_charge_7p4", "GRAVY", "Boman", "TPSA", "HBD", "HBA",
            "QED", "GI_stability", "DPP4_substrate_motif", "is_Top3"]
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in top:
            w.writerow(r)
    print(f"[OK] wrote {len(top)} entries -> {OUT}")
    # console summary
    print(f"{'pep':7s} {'iDPPIV':>7s} {'dG':>8s} {'MW':>7s} {'pI':>5s} {'Boman':>6s} {'QED':>6s}  GI")
    for r in top:
        print(f"{r['peptide']:7s} {float(r['iDPPIV_score']):7.3f} {float(r['dG_Vina']):8.3f} "
              f"{float(r['MW']):7.1f} {float(r['pI']):5.2f} {float(r['Boman']):6.2f} "
              f"{float(r['QED']):6.3f}  {r['GI_stability']}")


if __name__ == "__main__":
    main()
