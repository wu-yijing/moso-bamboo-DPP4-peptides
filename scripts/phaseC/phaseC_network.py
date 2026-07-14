#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Phase C2 — DPP4-centred network-pharmacology graph
====================================================
Builds a DPP4-centred pharmacology knowledge graph as the pure-computational
surrogate for wet-lab target/mechanism validation.

Sources (clearly labelled in the output):
  * STRING DB functional association network for DPP4 (Homo sapiens, 9606),
    fetched via REST API (network_type=functional, required_score=400).
    -> edges typed "functional_string" with STRING combined score.
  * Curated DPP4 cleavage-SUBSTRATE edges from the primary literature
    (DPP4 inactivates incretin & peptide hormones by X-Pro/X-Ala trimming).
    -> edges typed "cleavage_substrate_literature" with reference tag.

The graph contextualises WHY a DPP4 inhibitor is glucoregulatory:
inhibiting DPP4 raises active GLP-1/GIP -> insulin secretion -> T2DM benefit.
"""

import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "..", "data", "phaseC")
RAW = os.path.join(DATA, "dpp4_network.json")
OUT_JSON = os.path.join(DATA, "phaseC_network.json")
OUT_TXT = os.path.join(DATA, "phaseC_network_summary.txt")

# ---------------------------------------------------------------- annotations
# role/pathway knowledge for nodes appearing in the STRING network + a few
# well-established DPP4 peptide substrates added as literature nodes.
NODE_ANNOT = {
    "DPP4": {
        "gene": "DPP4", "uniprot": "P27487",
        "role": "Dipeptidyl peptidase 4 / CD26 — target protease. Trims "
                "N-terminal X-Pro/X-Ala dipeptides, inactivating incretin and "
                "peptide hormones.",
        "pathway": "Incretin / glucose homeostasis / immune co-stimulation",
    },
    "GCG": {
        "gene": "GCG", "uniprot": "P01275",
        "role": "Preproglucagon — source of GLP-1 and GLP-2, the principal "
                "DPP4 substrates secreted from intestinal L-cells.",
        "pathway": "Incretin (GLP-1) signalling",
    },
    "GIP": {
        "gene": "GIP", "uniprot": "P09681",
        "role": "Glucose-dependent insulinotropic polypeptide (gastric incretin); "
                "a major DPP4 substrate.",
        "pathway": "Incretin (GIP) signalling",
    },
    "ADA": {
        "gene": "ADA", "uniprot": "P00813",
        "role": "Adenosine deaminase — forms a covalent signalling complex with "
                "DPP4/CD26; immune modulation.",
        "pathway": "T-cell co-stimulation",
    },
    "CAV1": {
        "gene": "CAV1", "uniprot": "Q03135",
        "role": "Caveolin-1 — scaffolds DPP4 at the plasma membrane; "
                "co-stimulation of T-cells.",
        "pathway": "Membrane scaffolding / immune",
    },
    "CXCR4": {
        "gene": "CXCR4", "uniprot": "P61073",
        "role": "Receptor for CXCL12/SDF-1; DPP4 cleaves CXCL12, modulating "
                "haemopoietic stem-cell homing and chemotaxis.",
        "pathway": "Chemokine / stem-cell homing",
    },
    "PTPRC": {
        "gene": "PTPRC", "uniprot": "P08575",
        "role": "CD45 tyrosine phosphatase — DPP4-CD45 axis in T-cell signalling.",
        "pathway": "T-cell signalling",
    },
    "FN1": {
        "gene": "FN1", "uniprot": "P02751",
        "role": "Fibronectin — DPP4 adhesion/interaction partner.",
        "pathway": "Cell adhesion",
    },
    "ITGB1": {
        "gene": "ITGB1", "uniprot": "P05556",
        "role": "Integrin beta-1 — adhesion/migration interaction partner.",
        "pathway": "Cell adhesion / migration",
    },
    "PRCP": {
        "gene": "PRCP", "uniprot": "Q9H3Z3",
        "role": "Prolyl carboxypeptidase — DPP4-family serine protease homolog.",
        "pathway": "Protease family",
    },
    "ACE2": {
        "gene": "ACE2", "uniprot": "Q9BYF1",
        "role": "Angiotensin-converting enzyme 2 — renin-angiotensin system; "
                "relevant to T2DM cardiovascular comorbidity.",
        "pathway": "Renin-angiotensin / CVD",
    },
    # ---- literature-only DPP4 peptide substrates (added as enriched nodes) ----
    "CXCL12": {
        "gene": "CXCL12", "uniprot": "P48061",
        "role": "SDF-1 / CXCL12 — DPP4 substrate; chemotaxis & stem-cell homing.",
        "pathway": "Chemokine (literature substrate)",
        "literature_only": True,
    },
    "NPY": {
        "gene": "NPY", "uniprot": "P01303",
        "role": "Neuropeptide Y — DPP4 substrate (appetite / autonomic control).",
        "pathway": "Peptide hormone (literature substrate)",
        "literature_only": True,
    },
    "PYY": {
        "gene": "PYY", "uniprot": "P10082",
        "role": "Peptide YY — DPP4 substrate; satiety signalling.",
        "pathway": "Peptide hormone (literature substrate)",
        "literature_only": True,
    },
    "TAC1": {
        "gene": "TAC1", "uniprot": "P20366",
        "role": "Substance P — DPP4 substrate; sensory/inflammatory peptide.",
        "pathway": "Peptide hormone (literature substrate)",
        "literature_only": True,
    },
}

# Curated DPP4 cleavage-substrate edges (from primary literature).
# ref tags map to a short citation list printed in the summary.
LIT_SUBSTRATE_EDGES = [
    ("DPP4", "GCG",   "Mentlein 2009; Drucker 2006"),
    ("DPP4", "GIP",   "Mentlein 2009; Drucker 2006"),
    ("DPP4", "CXCL12","Christopherson 2002"),
    ("DPP4", "NPY",   "Mentlein 2009"),
    ("DPP4", "PYY",   "Mentlein 2009"),
    ("DPP4", "TAC1",  "Mentlein 2009"),
]


def main():
    raw = json.load(open(RAW, encoding="utf-8"))

    nodes = {}
    edges = []

    # STRING functional edges
    for e in raw:
        a, b = e["preferredName_A"], e["preferredName_B"]
        for n in (a, b):
            if n not in nodes:
                ann = NODE_ANNOT.get(n, {"gene": n, "role": "STRING neighbour",
                                          "pathway": "—"})
                nodes[n] = {"id": n, "name": n,
                            "gene": ann.get("gene", n),
                            "uniprot": ann.get("uniprot", ""),
                            "role": ann.get("role", ""),
                            "pathway": ann.get("pathway", ""),
                            "source": "string"}
        edges.append({
            "source": a, "target": b,
            "type": "functional_string",
            "score": round(float(e["score"]), 3),
            "nscore": e.get("nscore", 0), "tscore": e.get("tscore", 0),
            "ascore": e.get("ascore", 0), "escore": e.get("escore", 0),
        })

    # Enrich with literature substrate nodes + edges
    for (src, tgt, ref) in LIT_SUBSTRATE_EDGES:
        if tgt not in nodes:
            ann = NODE_ANNOT.get(tgt, {"gene": tgt, "role": "DPP4 substrate",
                                        "pathway": "—", "literature_only": True})
            nodes[tgt] = {"id": tgt, "name": tgt,
                          "gene": ann.get("gene", tgt),
                          "uniprot": ann.get("uniprot", ""),
                          "role": ann.get("role", ""),
                          "pathway": ann.get("pathway", ""),
                          "source": "literature_substrate"}
        edges.append({
            "source": src, "target": tgt,
            "type": "cleavage_substrate_literature",
            "score": None, "ref": ref,
        })

    graph = {
        "central_target": "DPP4",
        "central_target_uniprot": "P27487",
        "string_params": {
            "species": 9606, "network_type": "functional",
            "required_score": 400,
            "note": "fetched via STRING REST API (network endpoint)",
        },
        "nodes": list(nodes.values()),
        "edges": edges,
        "legend": {
            "functional_string": "STRING DB functional association (combined score 0-1)",
            "cleavage_substrate_literature": "DPP4 cleaves/inactivates this peptide hormone (primary literature)",
        },
    }
    json.dump(graph, open(OUT_JSON, "w", encoding="utf-8"),
              indent=2, ensure_ascii=False)

    # ---- readable summary ----
    lines = []
    lines.append("DPP4-centred network-pharmacology summary")
    lines.append("=" * 50)
    lines.append(f"Central target : DPP4 (UniProt P27487)")
    lines.append(f"Nodes          : {len(graph['nodes'])}")
    lines.append(f"  STRING      : {sum(1 for n in graph['nodes'] if n['source']=='string')}")
    lines.append(f"  Literature  : {sum(1 for n in graph['nodes'] if n['source']!='string')}")
    lines.append(f"Edges          : {len(graph['edges'])}")
    lines.append(f"  functional  : {sum(1 for e in edges if e['type']=='functional_string')}")
    lines.append(f"  substrate   : {sum(1 for e in edges if e['type']=='cleavage_substrate_literature')}")
    lines.append("")
    lines.append("Key pharmacology context (why DPP4 inhibition is glucoregulatory):")
    lines.append("  DPP4 trims N-terminal X-Pro/X-Ala from peptide hormones.")

    def role(name):
        return NODE_ANNOT.get(name, {}).get("role", "")
    for sub in ["GCG", "GIP", "CXCL12", "NPY", "PYY", "TAC1"]:
        if sub in nodes:
            lines.append(f"   - {sub:8s}: {role(sub)}")
    lines.append("")
    lines.append("Representative references:")
    lines.append("  Drucker DJ. Nat Rev Drug Discov 2006 (incretin/DPP4).")
    lines.append("  Mentlein R. Regul Pept 2009 (DPP4 substrates review).")
    lines.append("  Christopherson K et al. PNAS 2002 (CXCL12/SDF-1 cleavage).")
    lines.append("  Szabo et al. / Kameoka et al. (ADA, CAV1 complexes).")
    lines.append("")
    lines.append("HONEST LIMITATION: STRING functional edges are human PPI/")
    lines.append("association background, NOT proof that the bamboo peptides act")
    lines.append("on these nodes. Substrate edges are literature-derived context.")
    lines.append("No wet-lab target-validation was performed (pure in silico).")

    open(OUT_TXT, "w", encoding="utf-8").write("\n".join(lines))

    print("\n".join(lines))
    print(f"\n[INFO] Wrote {OUT_JSON}")
    print(f"[INFO] Wrote {OUT_TXT}")


if __name__ == "__main__":
    main()
