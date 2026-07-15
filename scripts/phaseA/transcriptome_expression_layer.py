# -*- coding: utf-8 -*-
"""
思路③（需外部数据部分）：转录组表达层（定性 + 可复现连接协议）
======================================================
把候选的"基因组存在性"升级为"可食嫩笋表达"证据。

现实约束（已用真实数据核验）:
  * UniProt TrEMBL 条目 (taxonomy 38705) 不含 PH0100... 基因模型交叉引用
    -> "UniProt accession -> 基因组基因模型" 直接映射在本环境断链。
  * 公开毛竹笋转录组(Peng 2013 de novo; MDPI 2020 'Pachyloen' 发育阶段) 均
    不以 UniProt accession 为键; 且无本地 BLAST / 大文件下载能力
    -> 逐候选 TPM 定量连接在沙箱内不可行。
  * 因此本层交付: (a) 资源审计; (b) 源蛋白"笋相关性分级"(基于与已发表
    笋高表达基因家族的重叠, 文献锚定, 非 TPM); (c) 定量连接协议脚本(留给有
    BLAST 主机的环境一键补齐)。

输入:
  data/moso_253.fasta                    # 253 条毛竹 UniProt 蛋白
  data/moso_candidates_idppiv_short.tsv  # 4950 条 2-6 aa 候选(含 SCM 分)
输出:
  data/phaseA/transcriptome_expression_audit.tsv
  data/phaseA/source_protein_shoot_relevance.tsv
  data/phaseA/transcriptome_expression_summary.txt
  docs/思路③_转录组表达层_章节草稿.md  (由单独脚本/人工整合)

零外部依赖（仅 stdlib）。
"""
import os
from collections import defaultdict, Counter

HERE = os.path.dirname(os.path.abspath(__file__))
FASTA = os.path.join(HERE, "..", "..", "data", "moso_253.fasta")
SHORT = os.path.join(HERE, "..", "..", "data", "moso_candidates_idppiv_short.tsv")
OUT_AUDIT = os.path.join(HERE, "..", "..", "data", "phaseA", "transcriptome_expression_audit.tsv")
OUT_REL = os.path.join(HERE, "..", "..", "data", "phaseA", "source_protein_shoot_relevance.tsv")
OUT_SUM = os.path.join(HERE, "..", "..", "data", "phaseA", "transcriptome_expression_summary.txt")

FINALS = ("LPPGP", "APPSQ", "APQIP")

# ---- 已发表毛竹笋/幼组织高表达基因家族 (Peng 2013 PLoS One e78944 Fig.3;
#      MDPI 2020 Forests 11(8):861 发育阶段笋转录组) ----
SHOOT_DOC_FAMILIES = (
    "cyclin(CYCA)", "expansin(EXP)", "fructokinase(FTK)", "beta-glucanase(BGL)",
    "auxin response factor(ARF)", "MYB", "MYC", "DOF", "SAUR",
    "auxin influx(AUX1)", "gibberellin receptor(GID/GID1)",
    "MADS-box", "Squamosa/SBP", "bZIP", "WRKY",
    "cellulose synthase(CesA)", "xylan biosynthesis", "cell-wall metabolism",
    "CDK/cyclin-dependent kinase", "hormone signaling",
)


def parse_fasta(path):
    recs = {}
    acc = None; seq = []; hdr = ""
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if acc:
                    recs[acc] = ("".join(seq), hdr)
                hdr = line[1:].strip()
                acc = hdr.split("|")[1] if "|" in hdr else hdr.split()[0]
                seq = []
            else:
                seq.append(line.strip())
        if acc:
            recs[acc] = ("".join(seq), hdr)
    return recs


def gene_name(hdr):
    for tok in hdr.split():
        if tok.startswith("GN="):
            return tok[3:]
    return ""


def categorize(hdr):
    d = hdr.lower()
    if any(k in d for k in ["storage protein", "globulin", "albumin", "11s albumin",
                             "seed storage", "ferritin", "vicilin", "legumin"]):
        return "storage_reserve"
    if any(k in d for k in ["photosystem", "chlorophyll", "rubisco", "psii", "psi ",
                             "light-harvesting", "thylakoid", "photosynth"]):
        return "photosynthesis"
    if any(k in d for k in ["cellulose", "xylan", "lignin", "cell wall", "tubulin",
                             "actin", "structural", "cytoskeleton", "pectin", "hemicell"]):
        return "structural"
    if any(k in d for k in ["defensin", "pathogen", "chitinase", "stress", "heat shock",
                             "disease", "resistance", "pr protein", "antimicrobial"]):
        return "defense_stress"
    if any(k in d for k in ["transcription", "ribosomal", "histone", "translation",
                             "wrky", "myb", "mads", "squamosa", "binding protein",
                             "transcription factor", "tf ", "nucle", "chromatin", "factor"]):
        return "transcription_regulation"
    if any(k in d for k in ["receptor", "auxin", "hormone", "signal", "responsive"]):
        return "signaling"
    if any(k in d for k in ["synthase", "synthase", "oxidase", "dehydrogenase", "transferase",
                             "hydrolase", "peptidase", "reductase", "lyase", "mutase",
                             "isomerase", "kinase", "phosphatase", "enzyme", "synth",
                             "cytochrome", "p450"]):
        return "metabolic_enzyme"
    if any(k in d for k in ["hypothetical", "uncharacterized", "predicted protein", "unknown protein"]):
        return "hypothetical"
    return "other"


def shoot_tier(cat, desc, gn):
    """基于与已发表笋高表达家族的重叠给定性分级 (high/medium/low)。"""
    d = (desc + " " + gn).lower()
    # 显式高相关家族关键词 (直接命中 Peng2013/MDPI2020 报道的笋高表达家族)
    high_kw = ["cyclin", "cdk", "expansin", "mads", "squamosa", "sbpl", "sbpa",
                "myb", "myc", "dof", "saur", "arf", "aux1", "gid", "gibberellin",
                "cellulose synthase", "cesa", "bzip", "wrky", "transcription factor",
                "xylan", "cell wall", "histone", "chromatin"]
    medhigh_kw = ["cytochrome", "p450", "glycosyltransferase", "transferase",
                   "fructokinase", "kinase", "signal", "hormone", "receptor"]
    if any(k in d for k in high_kw):
        return "high", "命中笋高表达家族关键词(转录因子/细胞周期/细胞壁/激素信号)"
    if cat == "transcription_regulation":
        return "high", "转录调控类(TF)为笋转录组最富集的功能类之一"
    if cat == "structural":
        return "high", "结构/细胞壁合成类(CesA等), 笋快速生长与可食质地核心"
    if cat == "metabolic_enzyme":
        if any(k in d for k in medhigh_kw):
            return "medium-high", "代谢酶(含P450/激酶/糖基转移酶), 笋幼嫩组织次生代谢活跃"
        return "medium", "一般代谢酶, 部分属组成型表达"
    if cat == "signaling":
        return "medium", "信号/激素类, 笋发育激素信号活跃"
    if cat == "photosynthesis":
        return "low-medium", "光合类, 嫩笋光合弱于成熟叶, 但幼组织仍有一定表达"
    if cat == "other":
        return "medium", "未归类(other), 需逐条功能注释判定, 保守计为中等"
    if cat in ("hypothetical", "defense_stress", "storage_reserve"):
        return "low", f"{cat} 类, 与笋特异性高表达关联弱"
    return "medium", "默认中等"


def load_cands(path):
    out = []
    with open(path, encoding="utf-8") as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 5:
                continue
            out.append((p[0], int(p[1]), float(p[2]), float(p[3])))
    return out


def main():
    recs = parse_fasta(FASTA)
    accs = list(recs.keys())
    seqs = [recs[a][0] for a in accs]
    n_prot = len(accs)

    prot_cat = {a: categorize(recs[a][1]) for a in accs}
    prot_tier = {}
    prot_reason = {}
    for a in accs:
        t, r = shoot_tier(prot_cat[a], recs[a][1], gene_name(recs[a][1]))
        prot_tier[a] = t
        prot_reason[a] = r
    tier_count = Counter(prot_tier.values())

    # 候选 -> 源蛋白集合 -> 笋分级(set 归属)
    cands = load_cands(SHORT)
    cand_tier_set = Counter()
    tier_rank = {"high": 3, "medium-high": 2, "medium": 2, "low-medium": 1, "low": 0}
    cand_best_tier = {}
    finals_info = {}
    for pep, L, sc, mn in cands:
        hits = [a for a, s in zip(accs, seqs) if pep in s]
        if not hits:
            continue
        ts = set(prot_tier[a] for a in hits)
        for t in ts:
            cand_tier_set[t] += 1
        best = max(ts, key=lambda x: tier_rank.get(x, 1))
        cand_best_tier[pep] = best
        if pep in FINALS:
            a0 = hits[0]
            finals_info[pep] = (len(hits), a0, recs[a0][1].split(" OS=")[0].split("|")[-1],
                                 gene_name(recs[a0][1]), prot_cat[a0], prot_tier[a0])

    n_high_cand = sum(1 for p in cand_best_tier if cand_best_tier[p] in ("high", "medium-high"))

    # ---------- 资源审计 ----------
    audit = [
        ["dataset", "tissue", "stage_resolution", "id_scheme", "shoot_specific",
         "uniprot_keyed", "joinable_now", "notes"],
        ["Peng 2013 (PMC3820679 / PLoS One e78944)", "shoots (6 heights) + culms",
         "height-graded", "de novo unigenes", "YES", "NO",
         "NO (needs BLAST vs our 253 proteins)",
         "6076 up / 4613 down DEGs; Suppl Table S4 putative genes (XLS); 81% reads to genome"],
        ["MDPI 2020 Forests 11(8):861 ('Pachyloen' shoots)", "shoots (H1-H12 young, M1-M12 mature)",
         "developmental stage", "reference-genome gene models (PH0100...)", "YES", "NO",
         "INDIRECT (needs acc<->PH0100 map, currently blocked)",
         "mapping rates Table 2; GO of shoot-specific genes Table 3; Suppl likely has expr matrix"],
        ["GSE104596 (Zhang 2018 BMC Plant Biol 18:125)", "seedlings / culms (Mock vs GA)",
         "GA treatment", "platform/transcript IDs", "NO (culm/seedling)", "NO",
         "PARTIAL", "GA-response; 6 samples; aligned to P.heterocycla-v1.0; not shoot"],
        ["GSE90517 / PRJNA354950 (Wang 2017 Plant J 91:684)", "rhizome",
         "rhizome-associated", "transcript/AS", "NO (rhizome)", "NO", "NO",
         "alternative splicing / polyadenylation; 16 samples"],
        ["GSE104951 (circRNA in culms)", "culms", "culm growth", "circRNA",
         "NO (culm)", "NO", "NO", "circular RNAs in rapidly growing culms"],
    ]

    # ---------- 输出 ----------
    with open(OUT_AUDIT, "w", encoding="utf-8") as f:
        for row in audit:
            f.write("\t".join(row) + "\n")

    with open(OUT_REL, "w", encoding="utf-8") as f:
        f.write("accession\tgene_name\tcategory\tshoot_tier\trationale\tdescription\n")
        for a in accs:
            f.write("\t".join([a, gene_name(recs[a][1]), prot_cat[a], prot_tier[a],
                                prot_reason[a], recs[a][1].split(" OS=")[0].split("|")[-1]]) + "\n")

    lines = []
    lines.append("=== 思路③ 转录组表达层（定性 + 可复现连接协议）— 摘要 ===")
    lines.append(f"源蛋白总数: {n_prot}")
    lines.append("")
    lines.append("源蛋白【笋相关性分级】分布 (基于与已发表笋高表达家族的重叠, 文献锚定):")
    for t, c in tier_count.most_common():
        lines.append(f"  {t:14s} {c:4d} 条源蛋白 ({100*c/n_prot:.1f}%)")
    lines.append("")
    lines.append(f"候选池 (4950) 按源蛋白集合做【笋分级 set 归属】(可跨多级, 合计>候选数):")
    for t, c in sorted(cand_tier_set.items(), key=lambda x: -x[1]):
        lines.append(f"  {t:14s} 触及 {c:5d} 条候选 ({100*c/len(cands):.1f}%)")
    lines.append("")
    n_best = len(cand_best_tier)
    lines.append(f"以'最有利源蛋白'判定: {n_high_cand}/{n_best} 条候选 ({100*n_high_cand/n_best:.1f}%) "
                 "其最优源蛋白属 high/medium-high 笋相关性 -> 候选主要来自笋高/中高表达的功能家族")
    lines.append("")
    lines.append("决赛候选的源蛋白 (笋相关性):")
    for pep in FINALS:
        if pep in finals_info:
            nh, a, nm, gn, cat, tier = finals_info[pep]
            lines.append(f"  {pep}: 主源={a} ({nm}, GN={gn}) 类别={cat} 笋分级={tier} 命中源蛋白数={nh}")
    lines.append("")
    lines.append("资源审计 (见 transcriptome_expression_audit.tsv): 5 个公开毛竹转录组资源; "
                 "仅 Peng2013 与 MDPI2020 为【笋特异性】, 但二者均不以 UniProt accession 为键。")
    lines.append("")
    lines.append("量化连接缺口 (诚实边界):")
    lines.append("  (1) UniProt TrEMBL (tax 38705) 条目不含 PH0100... 基因模型交叉引用 -> 直接映射断链;")
    lines.append("  (2) 无公开【以 UniProt 为键】的笋 TPM 矩阵; 沙箱无 BLAST / 大文件下载能力;")
    lines.append("  -> 逐候选 TPM 定量连接需: 在有 BLAST 主机上运行 scripts/phaseA/transcriptome_join_protocol.py")
    lines.append("     (对 de novo 笋转录组 vs 本 253 蛋白做序列比对, 或对接 MDPI2020 基因组基因模型表达矩阵)。")
    lines.append("  本层当前交付 = 资源审计 + 源蛋白定性笋分级(文献锚定) + 定量连接协议(可一键补齐)。")

    with open(OUT_SUM, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print("\nAUDIT  ->", OUT_AUDIT)
    print("REL     ->", OUT_REL)
    print("SUMMARY ->", OUT_SUM)


if __name__ == "__main__":
    main()
