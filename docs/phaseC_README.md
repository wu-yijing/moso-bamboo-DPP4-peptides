# Phase C — In silico ADMET & DPP4-centred Network Pharmacology

**Project positioning**: This is a purely computational study (no wet-lab). This phase is the **second purely computational validation layer** after the static MM validation of Phase B, replacing the functional validation that would otherwise depend on in vitro ADMET / Caco-2 permeability / serum-stability experiments. All outputs are computationally derived; no wet experiments were performed.

---

## 1. Objectives

1. **Peptide-level in silico ADMET profile** (C1): compute physicochemical and drug-likeness descriptors for the 60 Vina-docked peptides, and predict gastrointestinal (GI) stability based on protease cleavage-site rules.
2. **DPP4-centred network pharmacology** (C2): retrieve the human DPP4 protein–protein interaction / functional-association network from STRING DB, overlay a literature-curated DPP4 substrate layer, and build a pharmacological knowledge graph that explains the mechanistic context of "inhibit DPP4 → raise GLP-1/GIP → improve glycaemia".

---

## 2. Methods

### 2.1 C1 — Peptide descriptors (purely computational)
- **Sequence-derived descriptors** (BioPython `ProtParam`, v1.87): molecular weight (MW), isoelectric point (pI), net charge at pH 7.4, GRAVY (grand average of hydropathy), instability index, aliphatic index (`AI = A + 2.9·V + 3.9·(I+L)`, X in mol%), charge at pH 7.4.
- **Boman protein-binding index** (Boman 1995): `Boman = (1/n)·ΣΔf`, where Δf is the amino-acid water→cyclohexane transfer free energy. **A negative value → higher protein-binding propensity**.
- **RDKit molecular descriptors**: ExactMolWt, TPSA, HBD (NumHDonors), HBA (NumHAcceptors), rotatable bonds, QED drug-likeness. Peptides are built via `Chem.MolFromSequence` (2D descriptors, no conformational sampling needed).

### 2.2 C1 — GI-tract stability (protease cleavage rules)
For each peptide, count endopeptidase cleavage sites (simplified specificity model):
- **Trypsin**: after K/R (unless followed by P)
- **Chymotrypsin**: after F/Y/W
- **Pepsin**: after F/Y/W/L/M (broad hydrophobic)
- **Elastase**: after A/G/S/V
- **Carboxypeptidase / aminopeptidase**: C/N-terminal exopeptidases (an N-terminal Pro cap resists aminopeptidases; a C-terminal Pro resists carboxypeptidases)
- **DPP4 self-cleavage**: Pro/Ala at N-terminal position 2 (X-Pro / X-Ala motif) → the peptide itself is a DPP4 substrate

GI-stability grading centres on **endopeptidase sites**:
`0 sites AND N-terminal Pro cap → Very High`; `0 sites → High`; `1 site → Moderate`; `≥2 sites → Low`.

### 2.3 C2 — DPP4 network (STRING + literature)
- **STRING DB REST API** (`network`, species=9606, functional, required_score=400): retrieve the DPP4 (UniProt P27487) functional-association network → edge type `functional_string` (combined score 0–1).
- **Literature-curated substrate layer**: overlay peptide hormones cleaved/inactivated by DPP4 (GCG/GLP-1, GIP, CXCL12/SDF-1, NPY, PYY, TAC1/substance P) → edge type `cleavage_substrate_literature` (with citation tags).
- Script: `scripts/phaseC/phaseC_network.py`

---

## 3. Results

### 3.1 C1 — Profile of the 60 docked peptides
- All 60 carry an **X-Pro / X-Ala motif** (Pro or Ala at positions 1–2), consistent with the "X-Pro preference" DPP4 structural-preference filter applied during candidate selection — this supports binding but also means they will be rapidly degraded by DPP4 in vivo as substrates (see Limitations).
- **GI-stability distribution**: Moderate (1 endopeptidase site) 39 / Low (multiple sites) 11 / High (no endopeptidase site) 10. Free linear short peptides are generally sensitive to exopeptidases.
- **Top-3 descriptors**:
  | Peptide | dG (kcal/mol) | pI | Boman | GRAVY | GI grade |
  |---|---|---|---|---|---|
  | LPPQ | −6.109 | 5.53 | +0.375 | −0.725 | Moderate |
  | APSPE | −7.150 | 4.60 | −0.820 | −1.140 | Low |
  | LAPSP | −7.087 | 5.53 | +0.340 | +0.320 | Low |
  - APSPE has the most negative Boman index (−0.82), indicating a relatively higher protein-binding propensity.
- Full per-peptide data: `data/phaseC/phaseC_peptides.tsv`

### 3.2 C2 — DPP4-centred network
- **15 nodes / 32 edges** (STRING 11 nodes + 26 functional edges; literature 4 nodes + 6 substrate edges).
- Core STRING neighbours: **GCG** (GLP-1/GLP-2 precursor), **GIP** (gastric inhibitory polypeptide), **CXCR4** (CXCL12/SDF-1 receptor), **ADA** (adenosine deaminase, covalent complex with DPP4/CD26), **CAV1** (caveolin-1), **PRCP** (prolyl carboxypeptidase, DPP4-family homologue enzyme), **ACE2**, **FN1**, **ITGB1**, **PTPRC/CD45**.
- Mechanistic context: DPP4 trims and inactivates incretins such as GLP-1 and GIP via X-Pro/X-Ala → inhibiting DPP4 raises active GLP-1/GIP → promotes insulin secretion → improves glycaemia in type 2 diabetes.
- Outputs: `data/phaseC/phaseC_network.json` (machine-readable), `data/phaseC/phaseC_network_summary.txt` (human-readable summary).

---

## 4. Honest limitations (must be disclosed in the submitted Methods / Limitations)

1. **Descriptors are sequence-derived, not measured ADMET**: MW/pI/GRAVY/Boman/TPSA etc. are computed from sequence/structure and were not validated by any in vitro permeability, metabolism, or serum-stability experiment.
2. **GI stability is a simplified rule**: based on published protease specificities, not ex vivo digestion or in vivo animal stability assays; short peptides are generally sensitive to exopeptidases, so the conclusion is only a rough screen.
3. **DPP4 self-cleavage tension**: all 60 candidates bear the X-Pro/X-Ala motif (both the DPP4 preferential-binding motif and its substrate motif). This indicates that oral delivery requires structural stabilisation (e.g. cyclisation, N-terminal capping), otherwise they will be rapidly degraded by DPP4 in vivo — a known limitation that does not negate their in vitro / computational binding activity.
4. **The STRING network is a human interaction background**: it demonstrates that DPP4 sits within the incretin / glucose-metabolism context, but **does not prove** that moso-bamboo peptides act on these nodes; substrate edges are literature-curated context.
5. **The instability index has limited reference value for short peptides** (it was designed for full-length proteins); recorded only for completeness.
6. This Phase C is a **computational approximation** with no wet experiments (Caco-2, serum stability, and in vitro DPP4 inhibition IC₅₀ were all not performed).

---

## 5. Reproduction

```bash
# C1 peptide-level ADMET + GI stability
python scripts/phaseC/phaseC_peptides.py
#   -> data/phaseC/phaseC_peptides.tsv

# C2 DPP4 network (requires network access to fetch STRING; or reuse the stored raw_dpp4_network.json)
curl "https://string-db.org/api/json/network?identifiers=9606.ENSP00000353731&species=9606&required_score=400&network_type=functional" -o data/phaseC/dpp4_network.json
python scripts/phaseC/phaseC_network.py
#   -> data/phaseC/phaseC_network.json / phaseC_network_summary.txt
```

Dependencies: Python ≥3.10, `rdkit`, `biopython`, `numpy`, `requests` (for STRING fetch).
