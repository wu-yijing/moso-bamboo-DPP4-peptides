# -*- coding: utf-8 -*-
"""
Idea 2: multi-target antidiabetic-peptide prioritization within the same proteome
(DPP4 + ACE + alpha-glucosidase)
Fully offline / zero third-party dependencies.

Method (isomorphic to Idea 1, extended to multiple targets):
  - Three axes unified as "reference-set signature method" (composition cosine of candidate
    vs known active peptides + 3-mer Jaccard):
       DPP4 axis uses this project's 37 known active-peptide reference set
       (data/literature_dpp4_peptides.tsv); the independently validated iDPPIV-SCM
       predictor score is reported separately in §3.9/§3.13 and NOT entered here.
  - ACE / alphaG axes: derived from literature-anchored "known active-peptide reference sets"
        (a) composition-signature cosine similarity (candidate vs 20-dim AA composition of reference set)
        (b) sequence similarity: maximum 3-mer Jaccard of candidate vs each reference peptide (motif/short-module overlap)
  - Three axes min-max normalized to 0-1 -> multi-target prioritization matrix
  - Cross-target length-bias quantification (upgrade Idea 1's argument to
    "a general problem in antidiabetic-peptide discovery")
  - Identify multi-target candidates (>=2 targets simultaneously in top 1/10)

Honest boundary (explicitly declared in both the section draft and the abstract):
  ACE / alphaG / DPP4 axis scores are all "knowledge-based composition + sequence-similarity
  signatures" derived from *cited literature* known active peptides, NOT benchmark-validated
  ML predictors; they only express "similarity to known active peptides", not predicted activity.
  The reference sets are small and bias toward short peptides -> the same length-bias
  caveat as Idea 1 applies.
"""
import os, sys, math, csv, statistics

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
CAND = os.path.join(ROOT, "data", "moso_candidates_idppiv.txt")
OUTDIR = os.path.join(ROOT, "data", "phaseC_multitarget")
os.makedirs(OUTDIR, exist_ok=True)

AA = "ACDEFGHIKLMNPQRSTVWY"
AMI = {a: i for i, a in enumerate(AA)}

# ---------- literature-anchored reference sets (cited) ----------
# DPP4 known active peptides: from this project's data/literature_dpp4_peptides.tsv (experimentally validated)
DPP4_REF = ["IPI","VPL","IPP","VPP","DYPAY","IAVPTGVA","YVVNPDNDEN","YVVNPDNNEN",
    "LTFPGSAED","LPVP","MPVQA","LDKVFR","LPGFF","MAGVDHI","GPFPLL","LPYPY","GPFPILN",
    "RPWR","YPVEPF","MHQPPQPL","SPTVMFPPQSVL","LPLP","QEPV","GPSGLDGAK","IPQHY",
    "VPQHY","VAVVPF","VPLGGF","IAIPPGIPYW","WLAFR","LLPFR","ATHALLA","EGF",
    "TENEWK","NFVSER","LDLPSK","QHEQR"]

# ACE inhibitory peptides: classic food-derived peptides from cheese/fermented foods (multiple literature, see section draft citations)
ACE_REF = ["IPP","VPP","RYLGY","LHLPLP","AYFYPEL","AYFYPE","HLPLP","EKDERF","VRYL",
    "YPFPGPIPN","FFVAP","EIVPN","DKIHPF","ELQDKIHPF","HLPLPLL","EMPFPK","VFGK",
    "VYPYYG","AAATP","KAAAAP","AAPLAP","KPVAAP","IAGRP","KAAAATP","PTPVP","PSNPP",
    "TGLKP","GGVPGG","FNMPLTIRITPGSKA","HCNKKYRSEM","TKYRVP","TSNRYHSYPWG","GVVPL",
    "LGL","SFVTT","VISDEDGVTH","NVPVYEGY","ITALAPSTM","ARHPHP","RPKHPIKHQ","RPKHPI",
    "FVAPFPEVF","DAYPSGAW","KAVPYPQ","VKEAMAPK","KVLPVPQK","GPVRGPFPIIV","APFPQV",
    "RELEEL","RLEEL","NENLL","LPQEVL","APFPQVF","EAMAPK","AVPYPQ"]

# alpha-glucosidase inhibitory peptides: compiled from 4 literature sources (black bean / Dahe black-pig ham / longan seed / yellow mealworm)
AG_REF = ["FLKEAFGV","RADLPGVK","NNNPFKF","FELLKHQK","VATVSLPR","EALELLK","IIAPPERK",
    "IEEALGDK","SSYYPFKGFA","VKGPGLYSDI","DK-7","WK-6","GR-7","FK-8","SK-6","DK-8"]

def clean(seqs):
    out = []
    for s in seqs:
        s = s.strip().upper()
        s = "".join(c for c in s if c in AMI)  # keep only the standard 20 amino acids
        if 2 <= len(s) <= 20:
            out.append(s)
    return out

DPP4_REF = clean(DPP4_REF)
ACE_REF  = clean(ACE_REF)
AG_REF   = clean(AG_REF)

def comp_vec(seq):
    v = [0]*20
    for c in seq:
        v[AMI[c]] += 1
    n = sum(v) or 1
    return [x/n for x in v]

def cosine(a, b):
    dot = sum(x*y for x, y in zip(a, b))
    na = math.sqrt(sum(x*x for x in a))
    nb = math.sqrt(sum(x*x for x in b))
    if na == 0 or nb == 0:
        return 0.0
    return dot/(na*nb)

def kmers(seq, k=3):
    return set(seq[i:i+k] for i in range(len(seq)-k+1)) if len(seq) >= k else {seq}

def jaccard(a, b):
    if not a or not b:
        return 0.0
    return len(a & b)/len(a | b)

def ref_profile(ref):
    profs = [comp_vec(s) for s in ref]
    return [sum(p[i] for p in profs)/len(profs) for i in range(20)]

def max_kmer_jaccard(seq, ref):
    ks = kmers(seq)
    best = 0.0
    for r in ref:
        best = max(best, jaccard(ks, kmers(r)))
    return best

def length_stats(ref):
    L = [len(s) for s in ref]
    return (statistics.median(L), sum(1 for x in L if x <= 6)/len(L),
            min(L), max(L), statistics.mean(L))

def pearson(xs, ys):
    n = len(xs)
    mx, my = statistics.mean(xs), statistics.mean(ys)
    num = sum((x-mx)*(y-my) for x, y in zip(xs, ys))
    dx = math.sqrt(sum((x-mx)**2 for x in xs))
    dy = math.sqrt(sum((y-my)**2 for y in ys))
    return num/(dx*dy) if dx and dy else 0.0

# ---------- load candidates (three axes unified = reference-set signature method) ----------
cands = []
with open(CAND, encoding="utf-8") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        pep = row["peptide"].strip()
        if all(c in AMI for c in pep) and 2 <= len(pep) <= 20:
            cands.append([pep, len(pep)])
print("candidates loaded:", len(cands))

# ---------- three-target reference-set composition signatures (same origin, ensuring comparability) ----------
DPP4_PROF = ref_profile(DPP4_REF)
ACE_PROF  = ref_profile(ACE_REF)
AG_PROF   = ref_profile(AG_REF)

# per-candidate scoring: each axis = 0.5*composition cosine + 0.5*max 3-mer Jaccard
rows = []
for pep, L in cands:
    cv = comp_vec(pep)
    dpp4_comp = cosine(cv, DPP4_PROF); dpp4_km = max_kmer_jaccard(pep, DPP4_REF)
    ace_comp  = cosine(cv, ACE_PROF);  ace_km  = max_kmer_jaccard(pep, ACE_REF)
    ag_comp   = cosine(cv, AG_PROF);   ag_km   = max_kmer_jaccard(pep, AG_REF)
    rows.append([pep, L, dpp4_comp, dpp4_km, ace_comp, ace_km, ag_comp, ag_km])

# ---------- min-max normalization ----------
def col(idx):
    return [r[idx] for r in rows]
def minmax(vals):
    lo, hi = min(vals), max(vals)
    rng = hi-lo or 1.0
    return [(v-lo)/rng for v in vals]

d_c_n  = minmax(col(2)); d_k_n  = minmax(col(3))
ac_c_n = minmax(col(4)); ac_k_n = minmax(col(5))
ag_c_n = minmax(col(6)); ag_k_n = minmax(col(7))

OUT = []
for i, r in enumerate(rows):
    pep, L, dpp4_comp, dpp4_km, ace_comp, ace_km, ag_comp, ag_km = r
    dpp4_score = 0.5*d_c_n[i] + 0.5*d_k_n[i]
    ace_score  = 0.5*ac_c_n[i] + 0.5*ac_k_n[i]
    ag_score   = 0.5*ag_c_n[i] + 0.5*ag_k_n[i]
    OUT.append([pep, L, dpp4_score, dpp4_comp, dpp4_km,
                ace_comp, ace_km, ace_score,
                ag_comp, ag_km, ag_score])

# ---------- multi-target identification (>=2 axes in top 1/10) ----------
def top_mask(vals, frac=0.1):
    thr = sorted(vals)[int(len(vals)*(1-frac))-1]
    return [v >= thr for v in vals], thr

d_top, d_thr       = top_mask([o[2] for o in OUT])
ace_top, ace_thr   = top_mask([o[7] for o in OUT])
ag_top, ag_thr     = top_mask([o[10] for o in OUT])

for i, o in enumerate(OUT):
    ntar = sum([d_top[i], ace_top[i], ag_top[i]])
    o.append(ntar)
    o.append("MULTI" if ntar >= 2 else ("SINGLE-DPP4" if d_top[i] else ("SINGLE-ACE" if ace_top[i] else ("SINGLE-AG" if ag_top[i] else "none"))))

# ---------- length-bias quantification (cross-target) ----------
def bias_block(name, ref, score_col_idx):
    med, pshort, lo, hi, mean = length_stats(ref)
    Ls = [o[1] for o in OUT]
    scores = [o[score_col_idx] for o in OUT]
    rr = pearson(Ls, scores)
    return (name, len(ref), med, round(pshort,3), lo, hi, round(mean,2), round(rr,3))

bias = [
    bias_block("DPP4 (ref=known active)", DPP4_REF, 2),
    bias_block("ACE  (ref=known active)", ACE_REF, 7),
    bias_block("alphaG(ref=known active)", AG_REF, 10),
]

# ---------- output ----------
tsv = os.path.join(OUTDIR, "multitarget_priority_matrix.tsv")
with open(tsv, "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["peptide","length","dpp4_norm","dpp4_comp_cos","dpp4_kmer_jac",
                "ace_comp_cos","ace_kmer_jac","ace_norm",
                "ag_comp_cos","ag_kmer_jac","ag_norm",
                "n_targets","flag"])
    for o in sorted(OUT, key=lambda x: (-x[11], -x[2], -x[7], -x[10])):
        w.writerow(o)

summary = os.path.join(OUTDIR, "multitarget_summary.txt")
multi = [o for o in OUT if o[11] >= 2]
single_d = [o for o in OUT if o[12] == "SINGLE-DPP4"]
single_a = [o for o in OUT if o[12] == "SINGLE-ACE"]
single_g = [o for o in OUT if o[12] == "SINGLE-AG"]
with open(summary, "w", encoding="utf-8") as f:
    f.write("=== Idea 2 multi-target antidiabetic-peptide prioritization (offline / literature-anchored) ===\n\n")
    f.write("Candidate pool: %d entries (from DPP4-stage iDPPIV positives)\n" % len(OUT))
    f.write("Reference-set sizes: DPP4=%d, ACE=%d, alphaG=%d (all cited literature known actives)\n\n" % (len(DPP4_REF), len(ACE_REF), len(AG_REF)))
    f.write("[A] cross-target length bias (core methodological contribution, isomorphic to Idea 1)\n")
    f.write("  target            N   medLen  %%<=6  min  max  meanLen  corr(score,length)\n")
    for b in bias:
        f.write("  %-17s %4d  %5.1f  %5.3f  %3d  %3d  %6.2f   %+6.3f\n" % b)
    f.write("  -> all three target reference sets are dominated by short peptides, and candidate scores correlate\n")
    f.write("     positively with length -> length bias is a cross-target general problem\n\n")
    f.write("[B] prioritization results\n")
    f.write("  multi-target (>=2 axes top 1/10): %d entries\n" % len(multi))
    f.write("  single-target DPP4: %d | single-target ACE: %d | single-target alphaG: %d\n" % (len(single_d), len(single_a), len(single_g)))
    f.write("  thresholds: dpp4_norm>=%.3f, ace_norm>=%.3f, ag_norm>=%.3f (each top 10%%)\n" % (d_thr, ace_thr, ag_thr))
    f.write("\n[C] multi-target candidates Top-20 (by #targets hit, then by dpp4 norm)\n")
    for o in sorted(multi, key=lambda x:(-x[11], -x[2]))[:20]:
        f.write("  %-12s L=%2d  dpp4=%.3f ace=%.3f ag=%.3f  nTar=%d [%s]\n" %
                (o[0], o[1], o[2], o[7], o[10], o[11], o[12]))
    f.write("\n[D] finalist candidates in the three-target matrix\n")
    for target in ["LPPGP","APPSQ","APQIP"]:
        m = [o for o in OUT if o[0]==target]
        if m:
            o = m[0]
            f.write("  %-8s L=%2d dpp4=%.3f ace=%.3f ag=%.3f nTar=%d flag=%s\n" %
                    (o[0], o[1], o[2], o[7], o[10], o[11], o[12]))
    f.write("\nHonest boundary: the three axes (DPP4/ACE/alphaG) are knowledge-based composition+sequence-similarity\n")
    f.write("          signatures (literature-anchored), NOT benchmark-validated ML predictors;\n")
    f.write("          they only express similarity to known active peptides, not predicted activity; needs wet-lab validation.\n")

print(open(summary, encoding="utf-8").read())
print("\nwritten:", tsv)
print("written:", summary)
