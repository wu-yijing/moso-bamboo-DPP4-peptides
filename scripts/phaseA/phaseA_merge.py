#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Phase A - Step 3: 合并三个官方服务器的返回结果，得到"官方过滤"候选集。

流程:
  母集 (4,950 条 2-6aa) ── PeptideRanker (PR >= THR_PR)
                         ── AllerTOP v2  (非过敏原)
                         ── ToxinPred3   (Non-Toxin)
        └─> 三者交集 = 官方过滤候选集  ->  data/phaseA_inputs/official_candidates.tsv

用法:
  1. 把三个服务器的结果另存为 CSV/TSV（列名见 CONFIG 区注释）。
  2. 在 CONFIG 区填好文件路径与列名（或保持自动探测默认值）。
  3. 运行: python phaseA_merge.py
  4. 查看输出的漏斗计数与 official_candidates.tsv。

说明:
  - 服务器返回格式可能略有差异，本脚本对每类结果用"序列列 + 分值/标签列"
    做通用加载；列名可用 CONFIG 显式指定，未指定时按关键词自动探测。
  - 任一结果文件缺失时，该阶段会被跳过并在漏斗中标注 (skipped)，便于分阶段回填。
  - 仅用标准库 (csv)，无需额外依赖。
"""
import os
import csv

BASE = r"E:/workbuddy/moso-bamboo-DPP4-peptides/data/phaseA_inputs"
MASTER = os.path.join(BASE, "moso_short_2to6.txt")
OUT = os.path.join(BASE, "official_candidates.tsv")

# ============================== CONFIG ==============================
# 三个服务器结果文件路径（把网页/邮件结果另存为 csv 或 tsv 后填这里）
PR_FILE = os.path.join(BASE, "results_peptideranker.tsv")        # PeptideRanker 结果
ALLER_FILE = os.path.join(BASE, "results_allertop.tsv")         # AllerTOP v2 结果
TOX_FILE = os.path.join(BASE, "results_toxinpred.csv")           # ToxinPred3 结果

# 阈值 / 标签
THR_PR = 0.5                       # PeptideRanker: 生物活性概率 >= 0.5 视为有活性
NON_ALLERGEN_LABELS = ("non-allergen", "non_allergen", "nonallergen",
                         "non allergen", "0", "no", "negative")   # AllerTOP 阴性标签（小写匹配）
TOXIN_NEG_LABELS = ("non-toxin", "nontoxin", "non_toxin",
                      "0", "no", "negative")                      # ToxinPred 阴性标签（小写匹配）

# 列名自动探测关键词
SEQ_HINTS = ("sequence", "peptide", "seq", "pep")
PR_SCORE_HINTS = ("score", "prob", "probability", "rank", "pr")
ALLER_HINTS = ("allerg", "class", "prediction", "label", "result")
TOX_HINTS = ("prediction", "label", "class", "result", "toxin")

# 若自动探测不准，可在下方显式指定（留空 "" 则用自动探测）
PR_SEQ_COL = ""
PR_SCORE_COL = ""
ALLER_SEQ_COL = ""
ALLER_LABEL_COL = ""
TOX_SEQ_COL = ""
TOX_LABEL_COL = ""
# ====================================================================


def _find_col(header, hints, explicit):
    if explicit:
        return explicit
    hl = [h.lower() for h in header]
    for h in hints:
        for i, name in enumerate(hl):
            if h in name:
                return header[i]
    return None


def load_table(path):
    if not os.path.exists(path):
        return None, None
    with open(path, encoding="utf-8-sig", newline="") as f:
        sample = f.read(4096)
        f.seek(0)
        delim = "\t" if sample.count("\t") >= sample.count(",") else ","
        reader = csv.DictReader(f, delimiter=delim)
        rows = list(reader)
        return reader.fieldnames, rows


def load_pr(path, universe):
    fn, rows = load_table(path)
    if fn is None:
        return None
    seqc = _find_col(fn, SEQ_HINTS, PR_SEQ_COL)
    scrc = _find_col(fn, PR_SCORE_HINTS, PR_SCORE_COL)
    if seqc is None or scrc is None:
        print(f"  [PR] 列探测失败 header={fn}")
        return None
    out = {}
    for r in rows:
        seq = (r.get(seqc) or "").strip().upper()
        try:
            val = float(r.get(scrc))
        except (TypeError, ValueError):
            continue
        if seq:
            out[seq] = val
    print(f"  [PR] 读入 {len(out)} 条; 阈值 >= {THR_PR}")
    return out


def load_label(path, hints, explicit, neg_labels, name):
    fn, rows = load_table(path)
    if fn is None:
        return None
    seqc = _find_col(fn, SEQ_HINTS, explicit[0])
    labc = _find_col(fn, hints, explicit[1])
    if seqc is None or labc is None:
        print(f"  [{name}] 列探测失败 header={fn}")
        return None
    out = {}
    for r in rows:
        seq = (r.get(seqc) or "").strip().upper()
        lab = (r.get(labc) or "").strip().lower()
        if seq:
            out[seq] = lab
    print(f"  [{name}] 读入 {len(out)} 条")
    return out


def main():
    with open(MASTER, encoding="utf-8") as f:
        universe = [l.strip() for l in f if l.strip()]
    print(f"[母集] 2-6aa 短肽共 {len(universe)} 条\n")

    print("[加载结果]")
    pr = load_pr(PR_FILE, universe)
    aller = load_label(ALLER_FILE, ALLER_HINTS,
                       (ALLER_SEQ_COL, ALLER_LABEL_COL), NON_ALLERGEN_LABELS, "AllerTOP")
    tox = load_label(TOX_FILE, TOX_HINTS,
                     (TOX_SEQ_COL, TOX_LABEL_COL), TOXIN_NEG_LABELS, "ToxinPred")

    survivors = set(universe)
    funnel = [("母集 (2-6aa)", len(universe))]

    # PeptideRanker
    if pr is None:
        print("  [PR] 文件缺失 -> 跳过该阶段 (skipped)\n")
    else:
        passed = {s for s in survivors if pr.get(s, -1) >= THR_PR}
        survivors = passed
        funnel.append((f"PeptideRanker >= {THR_PR}", len(survivors)))

    # AllerTOP
    if aller is None:
        print("  [AllerTOP] 文件缺失 -> 跳过该阶段 (skipped)\n")
    else:
        passed = {s for s in survivors if aller.get(s, "?").lower() in NON_ALLERGEN_LABELS}
        removed = len(survivors) - len(passed)
        survivors = passed
        funnel.append((f"AllerTOP 非过敏原 (-{removed})", len(survivors)))

    # ToxinPred
    if tox is None:
        print("  [ToxinPred] 文件缺失 -> 跳过该阶段 (skipped)\n")
    else:
        passed = {s for s in survivors if tox.get(s, "?").lower() in TOXIN_NEG_LABELS}
        removed = len(survivors) - len(passed)
        survivors = passed
        funnel.append((f"ToxinPred Non-Toxin (-{removed})", len(survivors)))

    print("\n[漏斗]")
    for stage, n in funnel:
        print(f"  {stage:32s} : {n}")

    # 输出
    survivors = sorted(survivors)
    with open(OUT, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["peptide"])
        for s in survivors:
            w.writerow([s])
    print(f"\n[完成] 官方过滤候选集 {len(survivors)} 条 -> {OUT}")


if __name__ == "__main__":
    main()
