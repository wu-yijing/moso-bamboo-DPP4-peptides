# Idea ③ Transcriptome Expression Layer (edible-shoot expression evidence) — Section Draft

> Positioning: Second pillar of the bioinformatics methods paper *Reproducible Offline Framework* (genome presence → edible-part expression).
> This layer upgrades "candidates are genome-encoded, naturally existing" (local part of Idea ③) further to
> "candidate source proteins belong to functional families **highly expressed in edible shoots**," as qualitative evidence,
> to defuse the soft spot of "253 entries include non-food proteins, weak food-transformability."

---

## 3.x Edible-shoot expression evidence layer (transcriptome bridging)

### 3.x.1 Motivation and boundary
A purely computational study cannot provide wet-experiment "organ-level quantitative expression (TPM/FPKM)." However, published moso-bamboo transcriptome resources offer a **literature-anchored qualitative criterion** for "whether a candidate source protein belongs to a highly expressed functional family in edible shoots."
The goal of this layer is: without relying on wet experiments and without relying on BLAST / large-file downloads (unrunnable in this environment),
to deliver a reproducible, auditable "shoot relevance" grading, and to explicitly declare the remaining steps of the quantitative linkage.

### 3.x.2 Method
**Step 1 — Resource audit.** Systematically review public moso-bamboo transcriptomes (Table 1), distinguishing
(i) shoot-specific resources (Peng 2013 de novo shoot transcriptome; MDPI 2020 'Pachyloen' shoot developmental-stage transcriptome),
(ii) non-shoot controls (GA-treated culm/seedling, rhizome, culm circRNA).
Key fact: neither shoot-specific resource uses **UniProt accession as key**
(Peng 2013 = de novo unigene; MDPI 2020 = reference-genome gene model `PH0100…G…`).

**Step 2 — Gene-model linkage break verification.** Extract a real UniProt TrEMBL entry (`A0AA96NL62`) from our 253 source proteins and verify via JSON:
its cross-references contain only EMBL / AlphaFoldDB / InterPro / Pfam etc., **not** the `PH0100…` genome gene-model ID.
I.e., the direct mapping "UniProt accession → genome gene model" is **broken** in current public annotations.

**Step 3 — Qualitative "shoot relevance" grading of source proteins.** For the 253 source proteins, overlap their functional description with
published shoot-highly-expressed gene families (Peng 2013 Fig.3: CYCA, EXP, FTK, BGL, ARF, MYB, MYC,
DOF, SAUR, AUX1, GID/GID1; MDPI 2020: transcription factors, hormone signaling, cell-wall metabolism enrichment)
via keyword overlap, assigning a five-level grade high / medium-high / medium / low-medium / low
(grading logic in `scripts/phaseA/transcriptome_expression_layer.py`'s `shoot_tier`).
Candidate peptides are set-attributed by their source-protein set (consistent with presence mapping, avoiding first-hit bias).

### 3.x.3 Results

**Table 1 Public moso-bamboo transcriptome resource audit**

| Dataset | Tissue | Stage-resolved | ID scheme | Shoot-specific | UniProt key | Currently linkable |
|---|---|---|---|---|---|---|
| Peng 2013 (PMC3820679) | shoot (6 heights)+culm | height gradient | de novo unigene | yes | no | no (needs BLAST) |
| MDPI 2020 Forests 11(8):861 | shoot (H1–12 young/M1–12 mature) | developmental stage | genome model PH0100… | yes | no | indirect (needs acc↔PH0100 mapping) |
| GSE104596 (Zhang 2018) | seedling/culm (Mock/GA) | GA treatment | platform/transcript | no (culm) | no | partial |
| GSE90517 (Wang 2017) | rhizome | — | transcript/AS | no | no | no |
| GSE104951 | culm | growth | circRNA | no | no | no |

**Table 2 Source-protein "shoot relevance" grade distribution (n=253)**

| Grade | # proteins | share |
|---|---|---|
| high | 124 | 49.0% |
| medium-high | 31 | 12.3% |
| medium | 76 | 30.0% |
| low-medium | 21 | 8.3% |
| low | 1 | 0.4% |

Candidate pool (4,950 2–6 aa short peptides) set-attributed by source-protein set:
**high touches 65.1%**, medium-high 38.8%, medium 53.4%, low-medium 19.4%, low 3.6%.
By "most favorable source protein," **3,837 / 4,950 candidates (77.5%) have their best source protein in high/medium-high shoot relevance**.

**Shoot relevance of finalist candidates' source proteins:**
- **LPPGP** → B3VN36 cytochrome P450 73A33 (metabolic enzyme) → medium-high; P450 family broadly expressed in young-tissue secondary metabolism.
- **APPSQ** → A0A3Q8AYS5 Squamosa-promoter binding protein-like (SBP transcription factor) → high; SBP/SPL regulate shoot apical meristem and vegetative growth, highly expressed in young shoots.
- **APQIP** → X2F5C1 MADS-box protein 4 (transcription factor) → high; MADS-box core-regulates shoot apical meristem and shoot development.

### 3.x.4 Honesty boundary (methodological statement)
This layer **does not provide per-candidate TPM quantification**; the root causes are two, both explicitly declared within the method:
1. UniProt TrEMBL (tax 38705) entries lack `PH0100…` gene-model cross-references, so direct mapping is broken;
2. No public "UniProt-keyed" shoot TPM matrix exists; and this execution environment has no BLAST / large-file download capability.

The remaining quantitative linkage (per-candidate shoot expression magnitude) is encapsulated as a **one-click reproducible protocol script**
`scripts/phaseA/transcriptome_join_protocol.py`, to be run on a host with BLAST/MMseqs2:
- Path A: use our 253 proteins to align against the de novo shoot transcriptome, taking best-hit identity/coverage as evidence the source protein is detectable in shoots;
- Path B: join the MDPI 2020 genome gene-model expression matrix, requiring a first derivation of accession↔PH0100 mapping from the genome GFF (the same alignment can also yield it in reverse).

That is: **this layer's current deliverable = resource audit + source-protein qualitative shoot grading (literature-anchored) + quantitative-linkage protocol (fillable)**,
the three together forming a transparent, reproducible response to the "253 include non-food proteins" soft spot.

### 3.x.5 Summary
Candidates not only back-trace 100% to moso-bamboo chromosome-level genome annotated proteome (presence), but **their source proteins' main body
(77.5% of candidates) falls within published shoot high/medium-high expression functional families** (transcription factors, cell-wall synthesis, cell cycle, hormone signaling).
This advances "naturally existing" to "**naturally existing within the edible-shoot expression lineage**," providing computational-level edible-part evidence for candidates' food-transformability without wet experiments — forming a methodological contrast with Xie Peng 2026, who only reported "hand-picked edible proteins."
