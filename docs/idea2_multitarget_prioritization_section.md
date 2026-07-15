# Idea ② Section Draft: Multi-target Antidiabetic Peptide Prioritization within the Same Proteome (second pillar of the bioinformatics methods paper)

> Positioning: This paper is a **reproducible offline framework + predictor-bias quantification** methods paper; moso bamboo is only a case study.
> Idea ② sits between Idea ① (DPP4 predictor length bias) and Idea ③ (genome/transcriptome presence), upgrading "single-target screening" to "**multi-target (DPP4 + ACE + α-glucosidase) prioritization within the same proteome**," and for the first time upgrading Idea ①'s bias argument into a **cross-target general problem**.

---

## 2.5 Multi-target prioritization framework

### 2.5.1 Motivation and gap
Postprandial hyperglycemia in type 2 diabetes is driven cooperatively by multiple enzymes: DPP4 degrades incretins, ACE raises vasoconstriction, and α-glucosidase accelerates carbohydrate absorption. A single DPP4-inhibitory peptide (such as those screened by Xie Peng et al. 2026) can hardly cover the whole process; peptides **derived from the same proteome that can act on multiple targets simultaneously** have greater functional ingredient value. To our knowledge, however, **no moso bamboo (*Phyllostachys edulis*) proteome-scale multi-target antidiabetic peptide prioritization study exists**. Under zero-experiment, fully-offline constraints, this section uniformly maps 4,717 DPP4 pre-screen positive candidates onto three target axes and builds a reproducible prioritization matrix.

### 2.5.2 Definition of the three target-scoring axes (method and honesty boundary)
| Axis | Source | Method | Nature |
|---|---|---|---|
| **DPP4** | Project Phase A `iDPPIV-SCM` (665+665 benchmark, offline reproduction, independent ACC≈0.77) | Composition-type log-odds score | **Benchmark-validated predictor** (reliable) |
| **ACE** | Literature-anchored known-active-peptide reference set (N=55; classic cheese/fermented-food ACE peptides such as IPP/VPP/LHLPLP/RYLGY, multi-source citations) | (a) 20-dim amino-acid-composition cosine of candidate vs reference set; (b) max per-peptide 3-mer Jaccard of candidate vs reference set; equal-weight combination | **knowledge-type composition + sequence-similarity signature** (not a benchmark predictor) |
| **α-glucosidase** | Literature-anchored known-active-peptide reference set (N=16, compiled from 4 papers: black soybean FLKEAFGV/RADLPGVK/NNNPFKF; large black-pig ham FELLKHQK etc. 5; longan seed SSYYPFKGFA/VKGPGLYSDI; yellow mealworm 6) | same as ACE | same as above (smaller reference set, lower confidence) |

**Honesty boundary (must write):** The ACE / α-glucosidase axes' scores are "similarity signatures to known active peptides," derived from *cited literature*, **not benchmark-validated machine-learning predictors**; they only indicate that the candidate is "similar to confirmed active peptides" at the composition and motif level, **not equivalent to predicted activity**. The reference sets are small and short-peptide biased; the length-bias caveat revealed in Idea ① applies equally — see §2.5.4.

### 2.5.3 Prioritization matrix and multi-target identification
For each candidate, compute the three-axis normalized score (min-max to 0–1); rule: "each axis's top 10% = a hit on that target":
- **Multi-target candidates (≥2 axes simultaneously in top 10%): 218** (4.7% of pool).
- Single-target: DPP4 310 / ACE 256 / α-glucosidase 412.
- **Multi-target top candidates are highly enriched in Pro** (e.g. `LPPQGHIPEK`, `VVAPPER`, `EPPVK`, `VPPNPTPPPS`, `LPPMPAPAPVH`), consistent with the known structure–activity relationship that DPP4/ACE both prefer N/C-terminal Pro.

### 2.5.4 Cross-target length bias (core methodological contribution, echoing Idea ①)
Length statistics on the three targets' literature reference sets, plus "score–length" correlation on candidates:

| Reference set | N | Median length | ≤6 aa share | Length range | Mean | corr(score, length) |
|---|---|---|---|---|---|---|
| DPP4 (known active) | 37 | 6.0 | 70.3% | 3–12 | 5.97 | **+0.155** |
| ACE (known active) | 55 | 6.0 | 65.5% | 3–15 | 6.53 | **+0.414** |
| α-glucosidase (known active) | 16 | 7.5 | 37.5% | 2–10 | 5.88 | **+0.258** |

**Conclusion:** The known-active-peptide reference sets for all three targets are **dominated by short peptides** (median 6–7.5 aa, ≤6 aa share 37–70%), and the candidates' ACE/αG similarity scores are **positively correlated with length** (longer → more "like" known active peptides — because longer peptides contain more matchable composition mass and 3-mer motifs). This is **opposite in direction but homologous in origin** to the DPP4 benchmark in Idea ① ("short peptides dominate, SCM score negatively correlated with length"): regardless of the scoring scheme, **the length confounding of benchmark/reference sets seeps into candidate ranking**. → This upgrades Idea ①'s argument into a general claim: "**length bias is a cross-target universal problem in antidiabetic-peptide (multi-target) computational discovery; any training/signature based on short-peptide-enriched benchmarks requires explicit quantification and length-robust correction.**"

### 2.5.5 Positioning of the finalist candidates in the three-target matrix (linking Ideas ①/③)
| Finalist | Length | DPP4 norm | ACE norm | αG norm | # targets hit | Tag |
|---|---|---|---|---|---|---|
| **LPPGP** | 5 | 0.438 | 0.553 | 0.193 | **2** | MULTI (DPP4+ACE) |
| APPSQ | 5 | 0.242 | 0.394 | 0.236 | 0 | none |
| APQIP | 5 | 0.330 | 0.406 | 0.151 | 1 | SINGLE-DPP4 |

**Interpretation:**
- **LPPGP** is further highlighted as the **lead candidate** under this multi-target framework — beyond DPP4 (Idea ① priority, best MM-GBSA), it also enters the top 10% on the ACE axis, showing the broadest antidiabetic-peptide profile; consistent with its origin from a cytochrome P450 metabolic-enzyme family that is highly expressed in edible shoots (Idea ③).
- APPSQ / APQIP are **DPP4-focused** candidates (consistent with Idea ① SCM/docking conclusions), not reaching top 10% on ACE/αG axes — presented truthfully, without overclaiming multi-target activity.

---

## Differentiation from published work (usable in Discussion)
1. **First moso-bamboo multi-target antidiabetic peptide prioritization**: Xie Peng et al. 2026 only did single-target DPP4; this work is the first to incorporate ACE and α-glucosidase into a same-proteome framework.
2. **Methodological contribution, not another screening paper**: upgrade "length bias" from single-point DPP4 evidence to a cross-target general methodological caveat (§2.5.4), and provide a reproducible, zero-external-dependency scoring scheme of "composition + sequence similarity."
3. **Fully offline and reproducible**: the DPP4 axis is the benchmark-validated iDPPIV-SCM; the ACE/αG axes are literature-anchored signatures (scripts and reference sets released with the manuscript), recomputable by anyone.
4. **Front-loaded honesty boundary**: explicitly state ACE/αG are "similarity signatures" not predicted activity, and that reference-set short-peptide bias has been quantified — reviewer-friendly, avoiding the risk of "overclaiming multi-target activity."

## Outputs
- Script: `scripts/phaseB/multitarget_prioritization.py` (zero-dependency)
- Data: `data/phaseC_multitarget/multitarget_priority_matrix.tsv` (4,717×12 full matrix), `multitarget_summary.txt`
- Reference sets (literature-anchored): DPP4=37 / ACE=55 / αG=16, all provided with the manuscript and cited.

## To-add / future work
- If conditions allow, replace the ACE/αG axes with **benchmark-validated ML predictors** (e.g. pLM4ACE, PepBench `ace_inhibitory`) to upgrade this section from "signature" to "prediction"; currently the sandbox cannot download large benchmarks, so literature-anchored signatures are used, declared in the boundary.
- Wet-experiment cross-validation of multi-target candidates (triple IC₅₀ assay) is listed as future work.
