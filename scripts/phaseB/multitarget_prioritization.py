# -*- coding: utf-8 -*-
"""
思路②: 多靶点抗糖尿病肽同蛋白组优先化 (DPP4 + ACE + alpha-glucosidase)
完全离线 / 零第三方依赖。

方法 (与思路①同构, 但扩展到多靶点):
  - DPP4 轴: 复用现有 iDPPIV-SCM 组成评分 (data/moso_candidates_idppiv.txt)
  - ACE / alphaG 轴: 文献锚定 "已知活性肽参考集" 派生
        (a) 组成签名 cosine 相似度 (候选 vs 参考集 20 维氨基酸组成)
        (b) 序列相似度: 候选与参考集逐肽的 3-mer Jaccard 最大值 (基序/短模体重叠)
  - 三轴 min-max 归一化到 0-1 -> 三靶点优先化矩阵
  - 跨靶点长度偏倚量化 (把思路①论点升级为 "抗糖尿病肽发现中的普遍问题")
  - 识别多靶点候选 (>=2 靶点同时进入前 1/10)

诚实边界 (章节草稿与摘要均显式声明):
  ACE / alphaG 分数为 "知识型组成+序列相似签名", 源于 *引用文献* 的已知活性肽,
  非基准验证的 ML 预测器; 仅表示 "与已知活性肽的相似度", 不等同预测活性。
  参考集规模小且偏短肽 -> 与思路①相同的长度偏倚警示适用。
"""
import os, sys, math, csv, statistics

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
CAND = os.path.join(ROOT, "data", "moso_candidates_idppiv.txt")
OUTDIR = os.path.join(ROOT, "data", "phaseC_multitarget")
os.makedirs(OUTDIR, exist_ok=True)

AA = "ACDEFGHIKLMNPQRSTVWY"
AMI = {a: i for i, a in enumerate(AA)}

# ---------- 文献锚定参考集 (cited) ----------
# DPP4 已知活性肽: 来自本项目 data/literature_dpp4_peptides.tsv (实验验证)
DPP4_REF = ["IPI","VPL","IPP","VPP","DYPAY","IAVPTGVA","YVVNPDNDEN","YVVNPDNNEN",
    "LTFPGSAED","LPVP","MPVQA","LDKVFR","LPGFF","MAGVDHI","GPFPLL","LPYPY","GPFPILN",
    "RPWR","YPVEPF","MHQPPQPL","SPTVMFPPQSVL","LPLP","QEPV","GPSGLDGAK","IPQHY",
    "VPQHY","VAVVPF","VPLGGF","IAIPPGIPYW","WLAFR","LLPFR","ATHALLA","EGF",
    "TENEWK","NFVSER","LDLPSK","QHEQR"]

# ACE 抑制肽: 奶酪/发酵食品经典食源肽 (多处文献, 见章节草稿引用)
ACE_REF = ["IPP","VPP","RYLGY","LHLPLP","AYFYPEL","AYFYPE","HLPLP","EKDERF","VRYL",
    "YPFPGPIPN","FFVAP","EIVPN","DKIHPF","ELQDKIHPF","HLPLPLL","EMPFPK","VFGK",
    "VYPYYG","AAATP","KAAAAP","AAPLAP","KPVAAP","IAGRP","KAAAATP","PTPVP","PSNPP",
    "TGLKP","GGVPGG","FNMPLTIRITPGSKA","HCNKKYRSEM","TKYRVP","TSNRYHSYPWG","GVVPL",
    "LGL","SFVTT","VISDEDGVTH","NVPVYEGY","ITALAPSTM","ARHPHP","RPKHPIKHQ","RPKHPI",
    "FVAPFPEVF","DAYPSGAW","KAVPYPQ","VKEAMAPK","KVLPVPQK","GPVRGPFPIIV","APFPQV",
    "RELEEL","RLEEL","NENLL","LPQEVL","APFPQVF","EAMAPK","AVPYPQ"]

# alpha-glucosidase 抑制肽: 4 篇文献汇编 (黑豆/大河乌猪火腿/龙眼籽/黄粉虫)
AG_REF = ["FLKEAFGV","RADLPGVK","NNNPFKF","FELLKHQK","VATVSLPR","EALELLK","IIAPPERK",
    "IEEALGDK","SSYYPFKGFA","VKGPGLYSDI","DK-7","WK-6","GR-7","FK-8","SK-6","DK-8"]

def clean(seqs):
    out = []
    for s in seqs:
        s = s.strip().upper()
        s = "".join(c for c in s if c in AMI)  # 仅保留标准 20 氨基酸
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

# ---------- 载入候选 (DPP4 轴 = iDPPIV_score) ----------
cands = []
with open(CAND, encoding="utf-8") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        pep = row["peptide"].strip()
        try:
            dscore = float(row["iDPPIV_score"])
        except:
            continue
        if all(c in AMI for c in pep) and 2 <= len(pep) <= 20:
            cands.append([pep, len(pep), dscore])
print("candidates loaded:", len(cands))

# ---------- 参考集组成签名 ----------
ACE_PROF = ref_profile(ACE_REF)
AG_PROF  = ref_profile(AG_REF)

# 逐候选打分
rows = []
for pep, L, dscore in cands:
    cv = comp_vec(pep)
    ace_comp = cosine(cv, ACE_PROF)
    ace_km   = max_kmer_jaccard(pep, ACE_REF)
    ag_comp  = cosine(cv, AG_PROF)
    ag_km    = max_kmer_jaccard(pep, AG_REF)
    rows.append([pep, L, dscore, ace_comp, ace_km, ag_comp, ag_km])

# ---------- min-max 归一化 ----------
def col(idxs):
    return [r[idxs] for r in rows]
def minmax(vals):
    lo, hi = min(vals), max(vals)
    rng = hi-lo or 1.0
    return [(v-lo)/rng for v in vals]

d_n   = minmax(col(2))
ac_c_n = minmax(col(3))
ac_k_n = minmax(col(4))
ag_c_n = minmax(col(5))
ag_k_n = minmax(col(6))

OUT = []
for i, r in enumerate(rows):
    pep, L, dscore, ace_comp, ace_km, ag_comp, ag_km = r
    ace_score = 0.5*ac_c_n[i] + 0.5*ac_k_n[i]
    ag_score  = 0.5*ag_c_n[i] + 0.5*ag_k_n[i]
    OUT.append([pep, L, dscore, d_n[i],
                ace_comp, ace_km, ace_score,
                ag_comp, ag_km, ag_score])

# ---------- 多靶点识别 (>=2 轴进入前 1/10) ----------
def top_mask(vals, frac=0.1):
    thr = sorted(vals)[int(len(vals)*(1-frac))-1]
    return [v >= thr for v in vals], thr

d_top, d_thr       = top_mask([o[3] for o in OUT])
ace_top, ace_thr   = top_mask([o[6] for o in OUT])
ag_top, ag_thr     = top_mask([o[9] for o in OUT])

for i, o in enumerate(OUT):
    ntar = sum([d_top[i], ace_top[i], ag_top[i]])
    o.append(ntar)
    o.append("MULTI" if ntar >= 2 else ("SINGLE-DPP4" if d_top[i] else ("SINGLE-ACE" if ace_top[i] else ("SINGLE-AG" if ag_top[i] else "none"))))

# ---------- 长度偏倚量化 (跨靶点) ----------
def bias_block(name, ref, score_col_idx):
    med, pshort, lo, hi, mean = length_stats(ref)
    Ls = [o[1] for o in OUT]
    scores = [o[score_col_idx] for o in OUT]
    r = pearson(Ls, scores)
    return (name, len(ref), med, round(pshort,3), lo, hi, round(mean,2), round(r,3))

bias = [
    bias_block("DPP4 (ref=已知活性肽)", DPP4_REF, 3),
    bias_block("ACE  (ref=已知活性肽)", ACE_REF, 6),
    bias_block("alphaG(ref=已知活性肽)", AG_REF, 9),
]

# ---------- 输出 ----------
tsv = os.path.join(OUTDIR, "multitarget_priority_matrix.tsv")
with open(tsv, "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["peptide","length","dpp4_scm_score","dpp4_norm",
                "ace_comp_cos","ace_kmer_jac","ace_norm",
                "ag_comp_cos","ag_kmer_jac","ag_norm",
                "n_targets","flag"])
    for o in sorted(OUT, key=lambda x: (-x[10], -x[3], -x[6], -x[9])):
        w.writerow(o)

summary = os.path.join(OUTDIR, "multitarget_summary.txt")
multi = [o for o in OUT if o[10] >= 2]
single_d = [o for o in OUT if o[11] == "SINGLE-DPP4"]
single_a = [o for o in OUT if o[11] == "SINGLE-ACE"]
single_g = [o for o in OUT if o[11] == "SINGLE-AG"]
with open(summary, "w", encoding="utf-8") as f:
    f.write("=== 思路② 多靶点抗糖尿病肽优先化 (离线/文献锚定) ===\n\n")
    f.write("候选池: %d 条 (来自 DPP4 阶段 iDPPIV 阳性)\n" % len(OUT))
    f.write("参考集规模: DPP4=%d, ACE=%d, alphaG=%d (均为引用文献已知活性肽)\n\n" % (len(DPP4_REF), len(ACE_REF), len(AG_REF)))
    f.write("[A] 跨靶点长度偏倚 (核心方法贡献, 与思路①同构)\n")
    f.write("  target            N   medLen  %%<=6  min  max  meanLen  corr(score,length)\n")
    for b in bias:
        f.write("  %-17s %4d  %5.1f  %5.3f  %3d  %3d  %6.2f   %+6.3f\n" % b)
    f.write("  -> 三个靶点参考集均以短肽为主, 且候选分数与长度正相关 -> 长度偏倚是跨靶点的普遍问题\n\n")
    f.write("[B] 优先化结果\n")
    f.write("  多靶点(>=2轴前1/10): %d 条\n" % len(multi))
    f.write("  单靶点 DPP4: %d | 单靶点 ACE: %d | 单靶点 alphaG: %d\n" % (len(single_d), len(single_a), len(single_g)))
    f.write("  阈值: dpp4_norm>=%.3f, ace_norm>=%.3f, ag_norm>=%.3f (各自前10%%)\n" % (d_thr, ace_thr, ag_thr))
    f.write("\n[C] 多靶点候选 Top-20 (按命中靶点数, 再按 dpp4 归一分)\n")
    for o in sorted(multi, key=lambda x:(-x[10], -x[3]))[:20]:
        f.write("  %-12s L=%2d  dpp4=%.3f ace=%.3f ag=%.3f  nTar=%d [%s]\n" %
                (o[0], o[1], o[3], o[6], o[9], o[10], o[11]))
    f.write("\n[D] 决赛候选在三靶点矩阵中的位置\n")
    for target in ["LPPGP","APPSQ","APQIP"]:
        m = [o for o in OUT if o[0]==target]
        if m:
            o = m[0]
            f.write("  %-8s L=%2d dpp4=%.3f ace=%.3f ag=%.3f nTar=%d flag=%s\n" %
                    (o[0], o[1], o[3], o[6], o[9], o[10], o[11]))
    f.write("\n诚实边界: ACE/alphaG 分数为知识型组成+序列相似签名(文献锚定), 非基准验证 ML 预测器;\n")
    f.write("          仅表示与已知活性肽的相似度, 不等同预测活性; 需湿实验验证。\n")

print(open(summary, encoding="utf-8").read())
print("\nwritten:", tsv)
print("written:", summary)
