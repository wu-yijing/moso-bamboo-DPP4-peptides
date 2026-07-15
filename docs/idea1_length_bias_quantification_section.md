# Idea ① Section Draft: Quantifying the Length Bias of DPP-IV Predictors and Its Impact on Proteome-Scale Mining

> Positioning: The **methodological core contribution** of this study (primary bioinformatics angle). Moso bamboo serves only as a case study; the theme is "expose and quantify the length bias of black-box predictors and how it distorts candidate ranking."
> Data source: iDPPIV-SCM homologous public benchmark (train 1063 / test 266, balanced positive/negative); moso-bamboo candidate short-peptide pool (n=4950, 2–6 aa). All analyses are zero-external-dependency and reproducible (`scripts/phaseA/length_bias_analysis.py`).

---

## 3.X Quantifying predictor length bias and its impact on proteome-scale mining

**Background and motivation.** The reliability of computational DPP-IV inhibitory peptide predictors has long been questioned; their accuracy is believed to suffer from "length confounding" — positive samples are mostly short peptides while negative samples are mostly long peptides or full-length proteins (flagged repeatedly as an open problem by PepBenchmark, StackDPPIV, etc.). However, quantitative evidence is lacking on the extent to which length bias is misread as an "activity signal," and how it systematically distorts the candidate prioritization of "proteome-scale virtual mining." Using iDPPIV-SCM (Charoenkwan et al. 2020, locally reproduced offline) as a case object, this study provides the first quantitative results on the above questions and proposes a length-robust candidate-ranking protocol.

**Methods.**
1. Reproduce the **global-composition SCM** on the iDPPIV-SCM homologous public benchmark (score = Σ_a∈S P(a); decision threshold −1.148 optimized via nested 5-fold cross-validation); report independent-test ACC and MCC.
2. Define a **trivial length baseline**: peptide length ≤ T (T = 6/8/10/12) → classified as DPP-IV inhibitory positive; compare its ACC/MCC with SCM.
3. Tabulate the length distributions of benchmark positive/negative samples and compare SCM vs length-baseline ACC stratified by length (≤6 aa / >6 aa) to localize the bias source.
4. Within the moso-bamboo candidate short-peptide pool (n=4950, 2–6 aa), compute the SCM global score, length-normalized score (per-residue mean), independent SI-style probability, and their Spearman correlation with peptide length; propose a **length-robust residual reranking** (regress SCM score on length, take residuals, rerank) and compare the overlap between the original Top20 and the robust Top20 to assess the practical impact of bias on ranking.

**Results.**
- Global-composition SCM independent test **ACC = 0.771 (MCC = 0.555)**, close to the literature report (≈0.797); the position-specific variant shows lower ACC (0.620).
- However, a **trivial baseline relying only on peptide length** (len ≤ 10 → positive) already reaches **ACC = 0.820 (MCC = 0.648), higher than SCM itself**. In the long-peptide stratum (>6 aa) SCM is 0.692 vs length baseline 0.833; in the short-peptide stratum (≤6 aa) SCM is only 0.518 vs length baseline 0.800 — i.e., SCM is not superior to the pure-length heuristic in either length stratum.
- **Root cause**: benchmark positive samples have median length **4** (mean 5.31, range 2–18), negative samples median length **15** (mean 15.53, range 5–75); the length distributions are severely confounded (Table X1).
- Within the candidate pool, the length bias persists as "**shorter → higher score**": SCM global score ~ length ρ = **−0.188**, independent SI probability ~ length ρ = **−0.212**; whereas the length-normalized score (per-residue mean) reduces the correlation to ρ = **−0.056**, proving length normalization partially mitigates the bias.
- **Length-robust residual reranking** replaces 3 of the original Top20 (overlap 17/20); our three finalist candidates **LPPGP / APPSQ / APQIP** keep their ranks after residual reranking (16→16, 420→420, 107→107), indicating their prioritization is **robust to length bias**.

*Table X1. Benchmark length-distribution confounding and length-baseline performance (core evidence for Idea ①)*

| Metric | Value |
|---|---|
| Benchmark train / test | 1063 / 266 (half positive, half negative) |
| Global SCM independent ACC / MCC | 0.771 / 0.555 |
| Pure length baseline (len≤10→pos) ACC / MCC | 0.820 / 0.648 |
| Positive-sample length median / mean | 4 / 5.31 |
| Negative-sample length median / mean | 15 / 15.53 |
| Candidate-pool SCM score ~ length ρ | −0.188 |
| Candidate-pool SI probability ~ length ρ | −0.212 |
| Candidate-pool length-normalized score ~ length ρ | −0.056 |
| Top20 overlap after residual reranking | 17/20 |

**Conclusions and implications for mining.** A considerable part of the accuracy reported by iDPPIV-SCM originates from length confounding rather than a genuine residue-composition signal; in proteome-scale mining, ranking directly by SCM score systematically biases the candidate pool toward extremely short peptides. We recommend: (a) use **length-normalized score or residual reranking** for length-robust ranking; (b) when reporting any "high-activity candidate," **explicitly disclose the length distribution and perform stratified validation**. This quantitative result explains at the methodological level why "simply re-screening once more" struggles to yield reliable candidates, and also provides a bias-correction basis for the candidate-prioritization conclusion in §3.7 of this study — our finalist candidates are proven robust to length bias, so their prioritization does not depend on this confounded signal.

---

### Writing notes (to the authors)
- This is a "methodological contribution," not a "screening result": moso bamboo is a case study; do not write it as "we screened moso bamboo again."
- The differentiation point from Xie Peng 2026 is: they treated PeptideRanker as a black box; we **exposed and quantified the black box's bias**.
- Honesty boundary: **the predictor itself has limited accuracy (ACC≈0.77) and suffers length confounding**, so the SCM score in this study is used only for candidate-prioritization ranking and soft filtering, not for deterministic activity judgment.
- Downstream §3.9 (multi-target) and §3.x (genome/transcriptome bridging) will continue to follow this "length-robust" principle.
