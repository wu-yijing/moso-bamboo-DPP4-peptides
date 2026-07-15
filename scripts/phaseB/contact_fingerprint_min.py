#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LPPGP / APPSQ contact fingerprint (receptor-residue interface) -- geometry of the
*minimized complex* fully consistent with Section 3.6 / 3.7
==================================================================================
Reuses scripts/phaseB/phaseB_validation.py's compute_mm_bind (MMFF94s relaxation:
receptor backbone fixed, ligand + pocket side chains optimized; pocket = 9A cutoff);
from the minimized complex, compute per-receptor-residue contact counts
(total / hydrophobic / polar-HBond) -> contact fingerprint.
Output:
  data/phaseB/contact_fingerprint_LPPGP_APPSQ.tsv
  data/phaseB/contact_fingerprint_LPPGP_APPSQ.json
"""
import os, json, sys
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import phaseB_validation as pb

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT  = os.path.join(ROOT, "data", "phaseB")
os.makedirs(OUT, exist_ok=True)

PEPS = [p for p in pb.TOP3 if p["name"] in ("LPPGP", "APPSQ")]

CATALYTIC = {"GLU146", "ARG147"}
S1S2 = {"TYR166", "ASN150", "ASN151", "TYR132", "TYR173", "TYR183", "TRP168", "SER209"}
INTERFACE = {"ARG184", "GLU145", "GLU177", "LYS175", "ILE143", "ILE148", "ILE172",
             "ILE176", "ILE185", "LEU164", "LEU180", "VAL167", "VAL174", "THR144",
             "THR152", "THR186", "PRO149", "PRO178", "PRO181", "ALA130", "SER131",
             "ASP133", "ASP171", "ASN119", "ASN169", "ASN179", "TRP187", "TYR118",
             "PHE", "HIS"}
GLYCAN = {"NAG1201", "NAG1", "FUC"}
def role_of(rn, ri):
    k = "%s%d" % (rn, ri)
    if k in CATALYTIC: return "catalytic"
    if k in GLYCAN:    return "glycan"
    if k in S1S2:      return "S1S2"
    if k in INTERFACE: return "interface"
    return "pocket"

def receptor_contact_map(res):
    """Compute per-receptor-residue contacts from the minimized complex (res['_combo'])."""
    combo = res["_combo"]
    conf  = combo.GetConformer()
    pkt_idx  = res["_pkt_idx"]     # pocket heavy-atom indices within the combo
    lig_idx  = res["_lig_idx"]     # ligand heavy-atom indices within the combo
    pkt_atoms = res["_pkt_atoms"]  # metadata [(resname,resid,name,chain,...)]

    def sym(i): return combo.GetAtomWithIdx(i).GetSymbol()
    def is_donor(i):
        a = combo.GetAtomWithIdx(i)
        return a.GetSymbol() == "N" and any(n.GetSymbol() == "H" for n in a.GetNeighbors())
    def is_acc(i): return combo.GetAtomWithIdx(i).GetSymbol() in ("O", "N")

    rec_cnt = defaultdict(lambda: {"total":0,"hb":0,"hp":0})
    rec_hb  = set()
    for gi in lig_idx:
        gi_sym = sym(gi)
        for pi in pkt_idx:
            d = conf.GetAtomPosition(gi).Distance(conf.GetAtomPosition(pi))
            if d > 4.5: continue
            hydro = (gi_sym == "C" and sym(pi) == "C" and d <= 4.0)
            polar = ((is_donor(gi) or is_acc(gi)) and (is_donor(pi) or is_acc(pi)) and d <= 3.5)
            if not (hydro or polar): continue
            a = pkt_atoms[pi]
            key = (a["resname"], a["resid"])
            rec_cnt[key]["total"] += 1
            if hydro: rec_cnt[key]["hp"] += 1
            if polar:
                rec_cnt[key]["hb"] += 1
                rec_hb.add("%s%d" % (a["resname"], a["resid"]))
    return rec_cnt, rec_hb

def main():
    per_pep = {}
    for pep in PEPS:
        print(f"=== {pep['name']} : MMFF minimization + fingerprint ===", flush=True)
        res = pb.compute_mm_bind(pep)
        rec_cnt, rec_hb = receptor_contact_map(res)
        per_pep[pep["name"]] = {
            "rec": rec_cnt, "hb": rec_hb,
            "n_total": sum(c["total"] for c in rec_cnt.values()),
        }
        print(f"  [{pep['name']}] total_contacts={per_pep[pep['name']]['n_total']} "
              f"H-bond_residues={sorted(rec_hb)}", flush=True)

    all_keys = set()
    for d in per_pep.values():
        all_keys.update(d["rec"].keys())
    ordered = sorted(all_keys, key=lambda k: (k[1], k[0]))

    rows = []
    for (rn, ri) in ordered:
        role = role_of(rn, ri)
        lp = per_pep["LPPGP"]["rec"].get((rn, ri), {"total":0,"hb":0,"hp":0})
        sq = per_pep["APPSQ"]["rec"].get((rn, ri), {"total":0,"hb":0,"hp":0})
        rows.append({
            "res": "%s%d" % (rn, ri), "role": role,
            "LPPGP_contacts": lp["total"], "LPPGP_hbond": "Y" if lp["hb"]>0 else "",
            "LPPGP_hydrophobic": lp["hp"],
            "APPSQ_contacts": sq["total"], "APPSQ_hbond": "Y" if sq["hb"]>0 else "",
            "APPSQ_hydrophobic": sq["hp"],
        })

    tsv = os.path.join(OUT, "contact_fingerprint_LPPGP_APPSQ.tsv")
    with open(tsv, "w") as f:
        f.write("res\trole\tLPPGP_contacts\tLPPGP_hbond\tLPPGP_hydrophobic\t"
                "APPSQ_contacts\tAPPSQ_hbond\tAPPSQ_hydrophobic\n")
        for r in rows:
            f.write(f"{r['res']}\t{r['role']}\t{r['LPPGP_contacts']}\t{r['LPPGP_hbond']}\t"
                    f"{r['LPPGP_hydrophobic']}\t{r['APPSQ_contacts']}\t{r['APPSQ_hbond']}\t"
                    f"{r['APPSQ_hydrophobic']}\n")
    jsn = os.path.join(OUT, "contact_fingerprint_LPPGP_APPSQ.json")
    with open(jsn, "w") as f:
        json.dump({"rows": rows,
                   "summary": {n: {"n_total": per_pep[n]["n_total"],
                                   "hbonds": sorted(per_pep[n]["hb"])}
                               for n in per_pep}}, f, indent=2)
    print(f"[OK] wrote {tsv}\n[OK] wrote {jsn}", flush=True)

if __name__ == "__main__":
    main()
