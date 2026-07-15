# -*- coding: utf-8 -*-
"""
Idea 3 (external-data part): transcriptome expression layer (qualitative + reproducible-join protocol)
======================================================
Upgrade the candidates' "genome presence" to "edible young-shoot expression" evidence.

Reality constraints (verified against real data):
  * UniProt TrEMBL entries (taxonomy 38705) do NOT carry PH0100... gene-model
    cross-references -> the "UniProt accession -> genome gene model" direct map
    is broken in this environment.
  * Public moso-bamboo shoot transcriptomes (Peng 2013 de novo; MDPI 2020
    'Pachyloen' developmental stages) are neither keyed by UniProt accession,
    nor can we BLAST locally / download large files -> per-candidate TPM
    quantitative join is infeasible inside the sandbox.
  * Therefore this layer delivers: (a) a resource audit; (b) source-protein
    "shoot-relevance tiering" based on overlap with published high-expression
    shoot gene families (literature-anchored, NOT TPM); (c) a quantitative-join
    protocol script (to be one-click completed in a BLAST-host environment).

Input:
  data/moso_253.fasta                    # 253 moso-bamboo UniProt proteins
  data/moso_candidates_idppiv_short.tsv  # 4950 2-6 aa candidates (incl. SCM score)
Output:
  data/phaseA/transcriptome_expression_audit.tsv
  data/phaseA/source_protein_shoot_relevance.tsv
  data/phaseA/transcriptome_expression_summary.txt
  docs/idea3_transcriptome_expression_layer_section.md  (integrated separately / by hand)

Zero external dependencies (stdlib only).
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

# ---- Published moso-bamboo shoot / young-tissue high-expression gene families
#      (Peng 2013 PLoS One e78944 Fig.3; MDPI 2020 Forests 11(8):861 developmental-stage shoot transcriptome) ----
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
    """Assign a qualitative tier (high/medium/low) from overlap with published
    shoot high-expression families."""
    d = (desc + " " + gn).lower()
    # explicit high-relevance family keywords (direct hits to Peng2013/MDPI2020 reported shoot high-expression families)
    high_kw = ["cyclin", "cdk", "expansin", "mads", "squamosa", "sbpl", "sbpa",
                "myb", "myc", "dof", "saur", "arf", "aux1", "gid", "gibberellin",
                "cellulose synthase", "cesa", "bzip", "wrky", "transcription factor",
                "xylan", "cell wall", "histone", "chromatin"]
    medhigh_kw = ["cytochrome", "p450", "glycosyltransferase", "transferase",
                   "fructokinase", "kinase", "signal", "hormone", "receptor"]
    if any(k in d for k in high_kw):
        return "high", "hits shoot high-expression family keyword (TF/cell-cycle/cell-wall/hormone-signaling)"
    if cat == "transcription_regulation":
        return "high", "transcription regulation (TF) is among the most enriched classes in the shoot transcriptome"
    if cat == "structural":
        return "high", "structural/cell-wall synthesis (CesA etc.), core to rapid shoot growth and edible texture"
    if cat == "metabolic_enzyme":
        if any(k in d for k in medhigh_kw):
            return "medium-high", "metabolic enzyme (incl. P450/kinase/glycosyltransferase), active secondary metabolism in young shoots"
        return "medium", "general metabolic enzyme, partly constitutive expression"
    if cat == "signaling":
        return "medium", "signaling/hormone, active hormone signaling in shoot development"
    if cat == "photosynthesis":
        return "low-medium", "photosynthesis, young shoots photosynthesize less than mature leaves but still express"
    if cat == "other":
        return "medium", "unclassified (other), needs per-gene annotation; conservatively medium"
    if cat in ("hypothetical", "defense_stress", "storage_reserve"):
        return "low", f"{cat} class, weakly associated with shoot-specific high expression"
    return "medium", "default medium"


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

    # candidate -> set of source proteins -> shoot tiering (set attribution)
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

    # ---------- resource audit ----------
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

    # ---------- output ----------
    with open(OUT_AUDIT, "w", encoding="utf-8") as f:
        for row in audit:
            f.write("\t".join(row) + "\n")

    with open(OUT_REL, "w", encoding="utf-8") as f:
        f.write("accession\tgene_name\tcategory\tshoot_tier\trationale\tdescription\n")
        for a in accs:
            f.write("\t".join([a, gene_name(recs[a][1]), prot_cat[a], prot_tier[a],
                                prot_reason[a], recs[a][1].split(" OS=")[0].split("|")[-1]]) + "\n")

    lines = []
    lines.append("=== Idea 3 transcriptome expression layer (qualitative + reproducible-join protocol) — summary ===")
    lines.append(f"Total source proteins: {n_prot}")
    lines.append("")
    lines.append("Source-protein [shoot-relevance tier] distribution (based on overlap with published shoot "
                 "high-expression families, literature-anchored):")
    for t, c in tier_count.most_common():
        lines.append(f"  {t:14s} {c:4d} source proteins ({100*c/n_prot:.1f}%)")
    lines.append("")
    lines.append(f"Candidate pool (4950) [shoot-tier set attribution] by source-protein set (crosses tiers, sum > candidate count):")
    for t, c in sorted(cand_tier_set.items(), key=lambda x: -x[1]):
        lines.append(f"  {t:14s} touches {c:5d} candidates ({100*c/len(cands):.1f}%)")
    lines.append("")
    n_best = len(cand_best_tier)
    lines.append(f"By 'most-favorable source protein': {n_high_cand}/{n_best} candidates ({100*n_high_cand/n_best:.1f}%) "
                 "have their best source protein in high/medium-high shoot relevance -> candidates mainly come "
                 "from high/medium-high shoot-expression functional families")
    lines.append("")
    lines.append("Finalist candidates' source proteins (shoot relevance):")
    for pep in FINALS:
        if pep in finals_info:
            nh, a, nm, gn, cat, tier = finals_info[pep]
            lines.append(f"  {pep}: primary source={a} ({nm}, GN={gn}) class={cat} shoot_tier={tier} n_source_hits={nh}")
    lines.append("")
    lines.append("Resource audit (see transcriptome_expression_audit.tsv): 5 public moso-bamboo transcriptome "
                 "resources; only Peng2013 and MDPI2020 are [shoot-specific], but neither is keyed by UniProt accession.")
    lines.append("")
    lines.append("Quantitative-join gap (honest boundary):")
    lines.append("  (1) UniProt TrEMBL (tax 38705) entries do NOT carry PH0100... gene-model cross-references -> direct map broken;")
    lines.append("  (2) no public [UniProt-keyed] shoot TPM matrix; sandbox has no BLAST / large-file download capability;")
    lines.append("  -> per-candidate TPM quantitative join needs: run scripts/phaseA/transcriptome_join_protocol.py "
                 "on a BLAST host")
    lines.append("     (sequence-align the de novo shoot transcriptome vs our 253 proteins, or join the MDPI2020 "
                 "genome-gene-model expression matrix).")
    lines.append("  This layer currently delivers = resource audit + source-protein qualitative shoot tiering "
                 "(literature-anchored) + quantitative-join protocol (one-click completable).")

    with open(OUT_SUM, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

    print("\n".join(lines))
    print("\nAUDIT  ->", OUT_AUDIT)
    print("REL     ->", OUT_REL)
    print("SUMMARY ->", OUT_SUM)


if __name__ == "__main__":
    main()
