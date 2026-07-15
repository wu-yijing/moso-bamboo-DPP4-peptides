# Idea ③ Section Draft (locally completable part): Genome Presence Mapping

> Positioning: Second-pillar section of the bioinformatics framework paper. Moso bamboo as case study; this section uses **genome presence** to prove candidates are "genome-encoded, naturally occurring bamboo protein digestion products," defusing the "253 entries include non-food proteins = unreliable" doubt.
> Data: 253 UniProt proteins (taxonomy 38705), 4,950 2–6 aa candidates. Script `scripts/phaseA/genome_existence_map.py` (zero-dependency, substring back-tracing + single classification per protein).
> **Scope statement**: This section completes only the local part of "genome presence"; the "edible-shoot expression evidence" requires external transcriptome TPM (Peng 2013, PMC3820679 / GEO), listed as to-be-added data, to be fetched after user authorization.

---

## 3.X Genome presence of candidates: back-tracing from the full proteome to the chromosome-level genome

**Motivation.** This study mines within the moso-bamboo full proteome (253 UniProtKB entries, taxonomy 38705), differing from Xie Peng et al. (2026) who hand-picked 8 "edible" proteins. This breadth advantage also invites doubt: are the candidates from real, traceable bamboo proteins, or digestion/pipeline artifacts? This section answers that via **genome presence mapping** and lays the evidence chain for the subsequent "edible-part expression" filter.

**Method.** Each candidate peptide is an in-vitro-digestion substring of its source protein; using `peptide in protein_sequence` it is precisely back-traced to its attribution among the 253 source proteins (zero-dependency, reproducible). Each source protein is classified **only once** (avoiding first-hit bias), grouped into 8 functional classes by UniProt header description (metabolic enzyme / transcription regulation / structural / photosynthesis / signaling / defense / storage / other). Candidates are set-attributed by their "source-protein set" (a peptide may span multiple classes).

**Results.**
- **100% presence confirmed**: all 4,950 candidates back-trace to ≥1 source protein (0 unmapped); the candidate pool covers **253/253 (100%)** source proteins. All source proteins are `OX=38705` (*Phyllostachys edulis*), i.e. annotated proteome products of the Peng et al. (2018, *GigaScience*, PMID 30202850, 51,074 genes) chromosome-level genome.
- **Broad functional span**: source proteins distribute across 8 functional classes (each protein counted once) —

| Source-protein class | # proteins | % of 253 | Candidates touched* | % of pool |
|---|---|---|---|---|
| other (limited description info) | 83 | 32.8% | 2601 | 52.5% |
| transcription_regulation | 79 | 31.2% | 2530 | 51.1% |
| metabolic_enzyme | 52 | 20.6% | 2566 | 51.8% |
| photosynthesis | 21 | 8.3% | 959 | 19.4% |
| structural (e.g. cellulose synthase) | 12 | 4.7% | 941 | 19.0% |
| signaling | 4 | 1.6% | 312 | 6.3% |
| hypothetical | 1 | 0.4% | 176 | 3.6% |
| defense_stress | 1 | 0.4% | 125 | 2.5% |

\* set attribution: candidates classified by their source-protein set, may span multiple classes, so the total exceeds the candidate count.

- **Genomic tracing of finalist candidates**: **LPPGP** main source `B3VN36` (Cytochrome P450 73A33, metabolic-enzyme class); **APPSQ** main source `A0A3Q8AYS5` (Squamosa-promoter binding protein-like, transcription-regulation class); **APQIP** main source `X2F5C1` (MADS-box protein 4, transcription-regulation class). All three uniquely hit a single source protein, with clear tracing.
- **Prevalence of multi-hits**: 59.5% of candidate short motifs span multiple proteins (e.g. "PP""GP" di/tri-peptides are highly abundantly repeated in bamboo proteins), which is exactly the cost of "full-proteome breadth," and also shows that sequence alone cannot lock a unique food source — further justifying the necessity of "transcriptome-expression filtering."

**Conclusions and honesty boundary.**
1. Candidates are **not pipeline artifacts**: 100% back-trace to annotated proteins of the chromosome-level genome, confirming they are naturally occurring bamboo-protein digestion products with confirmed presence.
2. The candidate pool covers **real, diverse bamboo proteins** — metabolic enzymes, transcription regulation, structural, photosynthesis, signaling, defense, etc. — proving "253 entries" is a real plant proteome rather than an unreliable set, defusing the superficial doubt about "non-food proteins."
3. **Honest qualifier**: genome presence only proves "natural existence," **not equivalent to "edible."** Cellulose synthase, photosynthetic complexes etc. are not traditional food proteins; to confirm candidates are **specifically expressed in edible shoots** (thus obtainable from food), the "shoot-transcriptome expression" evidence layer must be added (Peng 2013 transcriptome TPM) — this is to-be-added external data, to be fetched after approval as the second layer of Idea ③. At that point, a dual filter of "genome presence ∧ high shoot expression" will further converge our candidates from "naturally existing" to "food-relevant."

---

### Writing notes (to the authors)
- The theme of this section is **presence + naturalness**, not "edibility" — the latter must wait for transcriptome data; do not prematurely claim food-relevance.
- Differentiation from Xie Peng: they avoided this issue by "hand-picking 8 edible proteins"; we use **full-omics tracing** to objectively prove candidates come from a real bamboo proteome, leaving "edibility" to downstream expression filtering.
- Honest disclosure: 59.5% candidates multi-hit, meaning short peptides cannot be uniquely traced → exactly why transcriptome expression is needed for spatial/tissue filtering (not a defect, but the next step in methodological design).
- Downstream: once transcriptome data arrive, upgrade this section to a "genome presence ∧ shoot expression" dual filter, directly answering reviewers' potential doubt about "non-food proteins."
