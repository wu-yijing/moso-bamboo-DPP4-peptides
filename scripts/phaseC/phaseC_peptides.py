#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Phase C1 — Peptide-level in silico ADMET & GI-stability profiling
===============================================================
Pure-computational surrogate for wet-lab ADMET / Caco-2 / serum-stability
validation. NO external dependencies, NO network.

For each of the 60 Vina-docked moso-bamboo peptides we compute:
  1. Sequence-derived physicochemical descriptors (BioPython ProtParam +
     Boman protein-binding index).
  2. RDKit molecular descriptors (MW, TPSA, HBD, HBA, rotatable bonds, QED)
     where the peptide can be built cleanly.
  3. GI-tract stability via literature protease-cleavage rules
     (pepsin / trypsin / chymotrypsin / elastase / carboxypeptidase /
     aminopeptidase) and DPP4 self-cleavage (X-Pro / X-Ala motif).

Honest caveat (see docs/phaseC_README.md):
  * These are SEQUENCE-DERIVED DESCRIPTORS, not measured ADMET values.
  * GI-stability rules are simplified specificity models, not ex-vivo assays.
  * A positive DPP4 self-cleavage flag means the peptide is itself a DPP4
    substrate motif (rapid in-vivo degradation) — a stability consideration,
    not a disqualifier (it still competes for the active site).
"""

import os
import csv
import sys

import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    HAVE_RDKIT = True
except Exception as e:
    HAVE_RDKIT = False
    sys.stderr.write(f"[WARN] RDKit unavailable: {e}\n")

# ---- Boman (1995) free-energy-of-transfer table (water -> cyclohexane) ----
BOMAN_DF = {
    'A': 0.5, 'R': 1.8, 'N': -1.5, 'D': -2.5, 'C': -1.0,
    'Q': -0.7, 'E': -3.0, 'G': 0.0, 'H': -1.0, 'I': 3.0,
    'L': 2.8, 'K': 1.8, 'M': 2.0, 'F': 3.7, 'P': -0.3,
    'S': -1.0, 'T': -0.3, 'W': 3.2, 'Y': 2.5, 'V': 2.4,
}
HYDROPHOBIC = set('AVLIFMWCY')
CHARGED = set('DEKRH')
PROLINE = set('P')

# --- repo-relative paths (scripts/phaseC/ -> repo root) ---
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DOCK = os.path.join(REPO, "docking")
DATA = os.path.join(REPO, "data")

# Docking result file (peptide column = single-letter sequence, plus dG)
# 方案 α: 用 RDKit 制备的 iDPPIV 队列干净结果（含真实 iDPPIV 分 + 校正后 dG）
DOCK_TSV = os.path.join(DOCK, "moso_dock_results_idppiv_clean.tsv")
OUT_TSV = os.path.join(DATA, "phaseC", "phaseC_peptides.tsv")

# Top3 (方案 α): 新 iDPPIV 队列 RDKit 制备下的最优结合肽
TOP3 = {"APQIP": -6.807, "LPPGP": -6.558, "APPSQ": -6.513}


def load_dock():
    """Parse the docking TSV: col0 = peptide sequence, dG col by header
    keyword (dg/dock/affinity/best) with last-numeric fallback."""
    peptides = []
    with open(DOCK_TSV, "r", encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        hl = [h.lower() for h in header]
        dG_idx = None
        for i, h in enumerate(hl[1:], start=1):
            if any(k in h for k in ("dg", "dock", "affinity", "best")):
                dG_idx = i
                break
        if dG_idx is None:  # fallback: last column
            dG_idx = len(hl) - 1
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            f = line.split("\t")
            seq = f[0].strip()
            try:
                dG = float(f[dG_idx])
            except (ValueError, IndexError):
                dG = None
            peptides.append((seq, dG))
    return peptides


def boman_index(seq):
    """Boman protein-binding potential = (1/n) * sum(Δf).
    Negative -> higher propensity to bind proteins."""
    total = 0.0
    n = 0
    for aa in seq:
        if aa in BOMAN_DF:
            total += BOMAN_DF[aa]
            n += 1
    return (total / n) if n else float("nan")


def rdkit_descriptors(seq):
    """Build peptide with RDKit and compute drug-likeness descriptors.
    Returns dict; empty on failure."""
    out = {"rdkit_mw": None, "tpsa": None, "hbd": None,
            "hba": None, "rot_bonds": None, "qed": None}
    if not HAVE_RDKIT:
        return out
    try:
        mol = Chem.MolFromSequence(seq)
        if mol is None:
            return out
        out["rdkit_mw"] = Descriptors.ExactMolWt(mol)
        out["tpsa"] = Descriptors.TPSA(mol)
        out["hbd"] = Descriptors.NumHDonors(mol)
        out["hba"] = Descriptors.NumHAcceptors(mol)
        out["rot_bonds"] = Descriptors.NumRotatableBonds(mol)
        try:
            out["qed"] = Descriptors.qed(mol)
        except Exception:
            out["qed"] = None
    except Exception as e:
        sys.stderr.write(f"[WARN] RDKit fail {seq}: {e}\n")
    return out


def gi_cleavage_sites(seq):
    """Count internal protease cleavage sites + terminal vulnerabilities
    using simplified specificity rules."""
    n = len(seq)
    sites = 0
    detail = {}

    # Trypsin: after K/R unless next is P
    tr = [i for i in range(n - 1) if seq[i] in "KR" and seq[i + 1] != "P"]
    # Chymotrypsin: after F/Y/W (aromatic)
    ch = [i for i in range(n - 1) if seq[i] in "FYW"]
    # Pepsin: after hydrophobic/aromatic (broad) F/Y/W/L/M/A
    pe = [i for i in range(n - 1) if seq[i] in "FYWLM"]
    # Elastase: after small residues A/G/S/V
    el = [i for i in range(n - 1) if seq[i] in "AGSV"]

    detail = {"trypsin": len(tr), "chymotrypsin": len(ch),
              "pepsin": len(pe), "elastase": len(el)}
    internal = len(tr) + len(ch) + len(pe) + len(el)

    # Carboxypeptidase: C-terminal vulnerable unless last is Pro
    cpa_blocked = (n >= 1 and seq[-1] == "P")
    cpa_vuln = not cpa_blocked
    # Aminopeptidase: free N-terminus vulnerable unless N-terminal Pro (blocks it)
    nterm_blocked = (n >= 1 and seq[0] == "P")
    nterm_vuln = not nterm_blocked

    # DPP4 self-cleavage: X-Pro or X-Ala motif (position 2 is P or A)
    dpp4_substrate = (n >= 2 and seq[1] in "PA")

    # Endoprotease internal cuts are the primary GI-stability determinant;
    # exopeptidase vulnerabilities are near-universal for free linear peptides.
    total = internal + (1 if cpa_vuln else 0) + (1 if nterm_vuln else 0)
    exo_vuln = (cpa_vuln or nterm_vuln)
    return {
        "internal_sites": internal,
        "total_sites": total,
        "cpa_vuln": int(cpa_vuln),
        "nterm_vuln": int(nterm_vuln),
        "cpa_blocked": int(cpa_blocked),
        "nterm_blocked": int(nterm_blocked),
        "exo_vuln": int(exo_vuln),
        "dpp4_substrate": dpp4_substrate,
        "detail": detail,
    }


def gi_class(seq, g):
    """Assign a GI-stability class based on endoprotease (internal) sites."""
    ins = g["internal_sites"]
    if ins == 0 and g["nterm_blocked"]:
        return "Very High (no endoprotease site; N-term Pro-blocked)"
    if ins == 0:
        return "High (no endoprotease site)"
    if ins == 1:
        return "Moderate (1 endoprotease site)"
    return "Low (multiple endoprotease sites)"


def main():
    peptides = load_dock()
    sys.stderr.write(f"[INFO] Loaded {len(peptides)} docked peptides\n")

    rows = []
    for seq, dG in peptides:
        try:
            pa = ProteinAnalysis(seq)
            mw = pa.molecular_weight()
            pI = pa.isoelectric_point()
            aa_pct = pa.amino_acids_percent
            gravy = pa.gravy()
            instab = pa.instability_index()
            # aliphatic index (BioPython >=1.80 removed the built-in);
            # AI = A + 2.9*V + 3.9*(I+L)  [X = mole percent]
            aliphatic = (aa_pct.get("A", 0.0)
                         + 2.9 * aa_pct.get("V", 0.0)
                         + 3.9 * (aa_pct.get("I", 0.0) + aa_pct.get("L", 0.0)))
            charge = pa.charge_at_pH(7.4)
        except Exception as e:
            sys.stderr.write(f"[WARN] ProtParam fail {seq}: {e}\n")
            continue

        boman = boman_index(seq)
        frac_hyd = sum(aa_pct.get(a, 0.0) for a in HYDROPHOBIC) / 100.0
        frac_chg = sum(aa_pct.get(a, 0.0) for a in CHARGED) / 100.0
        frac_pro = aa_pct.get("P", 0.0) / 100.0

        rd = rdkit_descriptors(seq)
        g = gi_cleavage_sites(seq)
        gclass = gi_class(seq, g)

        rows.append({
            "peptide": seq,
            "length": len(seq),
            "biopy_mw": round(mw, 2),
            "rdkit_mw": (round(rd["rdkit_mw"], 2) if rd["rdkit_mw"] else ""),
            "pI": round(pI, 2),
            "net_charge_7p4": round(charge, 2),
            "GRAVY": round(gravy, 3),
            "instability_index": round(instab, 1),
            "aliphatic_index": round(aliphatic, 1),
            "boman_index": round(boman, 3),
            "frac_hydrophobic": round(frac_hyd, 3),
            "frac_charged": round(frac_chg, 3),
            "frac_proline": round(frac_pro, 3),
            "TPSA": (round(rd["tpsa"], 1) if rd["tpsa"] else ""),
            "HBD": (rd["hbd"] if rd["hbd"] is not None else ""),
            "HBA": (rd["hba"] if rd["hba"] is not None else ""),
            "rot_bonds": (rd["rot_bonds"] if rd["rot_bonds"] is not None else ""),
            "QED": (round(rd["qed"], 3) if rd["qed"] is not None else ""),
            "GI_internal_sites": g["internal_sites"],
            "GI_total_sites": g["total_sites"],
            "GI_cpa_vuln": g["cpa_vuln"],
            "GI_nterm_vuln": g["nterm_vuln"],
            "GI_cpa_blocked": g["cpa_blocked"],
            "GI_nterm_blocked": g["nterm_blocked"],
            "GI_exo_vuln": g["exo_vuln"],
            "DPP4_substrate_motif": int(g["dpp4_substrate"]),
            "GI_stability_class": gclass,
            "dG_Vina": (round(dG, 3) if dG is not None else ""),
            "is_Top3": int(seq in TOP3),
        })

    # sort by dG (most negative first)
    rows.sort(key=lambda r: r["dG_Vina"] if isinstance(r["dG_Vina"], float) else 0.0)

    fields = list(rows[0].keys())
    os.makedirs(os.path.dirname(OUT_TSV), exist_ok=True)
    with open(OUT_TSV, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    sys.stderr.write(f"[INFO] Wrote {len(rows)} rows -> {OUT_TSV}\n")

    # ---- console summary ----
    print(f"Peptides profiled : {len(rows)}")
    n_top = sum(r["is_Top3"] for r in rows)
    print(f"Top3 present     : {n_top}/3")
    # GI class distribution
    from collections import Counter
    gc = Counter(r["GI_stability_class"] for r in rows)
    print("GI-stability distribution:")
    for k, v in gc.most_common():
        print(f"   {k:35s} {v}")
    n_sub = sum(r["DPP4_substrate_motif"] for r in rows)
    print(f"DPP4-substrate-motif peptides : {n_sub}/{len(rows)}")
    print("\nTop3 detail:")
    for r in rows:
        if r["is_Top3"]:
            print(f"   {r['peptide']:6s} dG={r['dG_Vina']:>7} pI={r['pI']:>5} "
                  f"boman={r['boman_index']:>7} GRAVY={r['GRAVY']:>6} "
                  f"GI={r['GI_stability_class']:30s} DPP4sub={r['DPP4_substrate_motif']}")


if __name__ == "__main__":
    main()
