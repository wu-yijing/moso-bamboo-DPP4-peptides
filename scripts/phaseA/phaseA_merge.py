#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Phase A — Step 3: 合并三个官方服务器的返回结果，得到"官方过滤"候选集。

流程:
  母集 (4,950 条 2-6aa)
    ── [可选] PeptideRanker (PR >= THR_PR)  — 服务器暂不可用，手动提交后回填
    ── AlgPred 2.0  (ML Score < THR_ALGPRED) — 自动提交
    ── ToxinPred3   (Non-Toxin)               — 自动提交
    └─> 交集 = 官方过滤候选集

用法:
  1. 确保以下文件存在：
     data/phaseA_inputs/results_toxinpred.csv   (由 phaseA_run_toxinpred.py 生成)
     data/phaseA_inputs/results_algpred.csv     (由 phaseA_run_algpred.py 生成)
  2. (可选) 放 results_peptideranker.tsv 到同一目录
  3. 运行: python phaseA_merge.py
  4. 查看漏斗计数与 official_candidates.tsv

注意:
  - PeptideRanker 服务器当前 503（peptide.ucd.ie），恢复后手动提交并在下方
    取消掉 PR_FILE 路径前的注释即可启用。
  - AlgPred 对短肽（2-6aa）的过敏原预测阈值需结合 ML Score 定量使用，
    默认 THR_ALGPRED=0.6（ML Score < 0.6 视为非过敏原）。
    若所有肽均被判为 Allergen，可调高阈值；若希望收紧，调低阈值。
  - ToxinPred 返回的预测标签是 "Non-Toxin" / "Toxin"。
"""
import os
import csv

BASE = os.path.normpath(r"E:/workbuddy/moso-bamboo-DPP4-peptides/data/phaseA_inputs")
MASTER = os.path.join(BASE, "moso_short_2to6.txt")
OUT = os.path.join(BASE, "official_candidates.tsv")

# ============================== CONFIG ==============================
# --- 结果文件路径 ---
# PeptideRanker（手动提交后取消注释）
# PR_FILE = os.path.join(BASE, "results_peptideranker.tsv")

# AlgPred 2.0（自动生成）
ALG_FILE = os.path.join(BASE, "results_algpred.csv")

# ToxinPred3（自动生成）
TOX_FILE = os.path.join(BASE, "results_toxinpred.csv")

# --- 阈值 ---
THR_PR = 0.5                     # PeptideRanker bioactivity threshold
THR_ALGPRED = 0.6               # AlgPred ML Score cutoff (< 0.6 = non-allergen)
TOXIN_NEG_LABELS = ("non-toxin", "nontoxin", "non_toxin",
                    "0", "no", "negative")
# ====================================================================


def load_csv(path):
    """Load CSV file, return (header, rows) or (None, None) if missing."""
    if not os.path.exists(path):
        return None, None
    with open(path, encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        return reader.fieldnames, rows


def build_seq_map(header, rows, seq_col_hints, value_col, convert=None):
    """Build {sequence_upper: value} map from CSV rows."""
    if header is None:
        return None

    hl = [h.lower() for h in header]
    seq_col = None
    for hint in seq_col_hints:
        for i, name in enumerate(hl):
            if hint in name:
                seq_col = header[i]
                break
        if seq_col:
            break

    if seq_col is None:
        print(f"  [WARN] 找不到序列列 (hints={seq_col_hints}, header={header})")
        return None

    out = {}
    for r in rows:
        seq = (r.get(seq_col) or "").strip().upper()
        if not seq:
            continue
        # Also try peptide_id column if seq is short (e.g., from numbering)
        if len(seq) < 2 and r.get("sequence", "").strip():
            seq = r["sequence"].strip().upper()
        val_str = r.get(value_col, "").strip()
        if convert:
            try:
                val = convert(val_str)
            except (ValueError, TypeError):
                continue
        else:
            val = val_str.lower()
        out[seq] = val

    return out


def build_toxinpred_map():
    """Load ToxinPred results. Returns {sequence: prediction_label}."""
    fn, rows = load_csv(TOX_FILE)
    if fn is None:
        return None

    # Columns: peptide_id, sequence, svm_score, prediction
    out = {}
    for r in rows:
        seq = (r.get("sequence") or "").strip().upper()
        pred = (r.get("prediction") or "").strip()
        if seq and pred:
            out[seq] = pred.lower()
    print(f"  [ToxinPred] 读入 {len(out)} 条")
    return out


def build_algpred_map():
    """
    Load AlgPred results.
    AlgPred doesn't return sequence in results, only peptide_id and ML score.
    We need to map pep_XXXXX IDs back to the sequences from the master FASTA.
    Returns {sequence: ml_score_float}
    """
    fn, rows = load_csv(ALG_FILE)
    if fn is None:
        return None

    # Build ID→seq map from master FASTA
    id_to_seq = {}
    with open(os.path.join(BASE, "moso_short_2to6.fasta")) as f:
        lines = f.read().strip().split("\n")
    cur_id = None
    for line in lines:
        if line.startswith(">"):
            cur_id = line[1:].strip()
        elif cur_id:
            id_to_seq[cur_id] = line.strip().upper()
            cur_id = None

    out = {}
    for r in rows:
        pid = (r.get("peptide_id") or "").strip()
        score_str = (r.get("ml_score") or "").strip()
        seq = id_to_seq.get(pid, "")
        if seq and score_str:
            try:
                out[seq] = float(score_str)
            except ValueError:
                pass
    print(f"  [AlgPred] 读入 {len(out)} 条, 阈値 ML Score < {THR_ALGPRED}")
    return out


def main():
    with open(MASTER) as f:
        universe = [l.strip() for l in f if l.strip()]
    print(f"[母集] 2-6aa 短肽共 {len(universe)} 条\n")

    print("[加载结果]")

    # ToxinPred (automated)
    tox = build_toxinpred_map()
    if tox is None:
        print("  [ToxinPred] 文件缺失 -> 跳过")

    # AlgPred (automated)
    alg = build_algpred_map()
    if alg is None:
        print("  [AlgPred] 文件缺失 -> 跳过")

    # PeptideRanker (manual, optional)
    pr = None
    pr_file_path = None
    for candidate in ["results_peptideranker.tsv", "results_peptideranker.csv"]:
        p = os.path.join(BASE, candidate)
        if os.path.exists(p):
            pr_file_path = p
            break

    if pr_file_path:
        print(f"  [PR] 找到文件: {pr_file_path}")
        fn, rows = load_csv(pr_file_path)
        if fn:
            seq_col = None
            score_col = None
            hl = [h.lower() for h in fn]
            for h, name in zip(hl, fn):
                if any(x in h for x in ("sequence", "peptide", "seq", "pep")):
                    seq_col = name
                if any(x in h for x in ("score", "prob", "probability", "rank")):
                    score_col = name
            if seq_col and score_col:
                pr = {}
                for r in rows:
                    seq = (r.get(seq_col) or "").strip().upper()
                    try:
                        val = float(r.get(score_col))
                    except (TypeError, ValueError):
                        continue
                    if seq:
                        pr[seq] = val
                print(f"  [PR] 读入 {len(pr)} 条, 阈值 >= {THR_PR}")
    else:
        print("  [PR] 文件缺失 -> 跳过 (PeptideRanker 服务器暂不可用)")

    # ── Funnel ──
    survivors = set(universe)
    funnel = [("母集 (2-6aa)", len(universe))]

    if pr is not None:
        passed = {s for s in survivors if pr.get(s, -1) >= THR_PR}
        survivors = passed
        funnel.append((f"PeptideRanker >= {THR_PR}", len(survivors)))

    if alg is not None:
        passed = {s for s in survivors if alg.get(s, 999) < THR_ALGPRED}
        removed = len(survivors) - len(passed)
        survivors = passed
        funnel.append((f"AlgPred ML < {THR_ALGPRED} (-{removed})", len(survivors)))

    if tox is not None:
        passed = {s for s in survivors if tox.get(s, "?") in TOXIN_NEG_LABELS}
        removed = len(survivors) - len(passed)
        survivors = passed
        funnel.append((f"ToxinPred Non-Toxin (-{removed})", len(survivors)))

    print("\n[漏斗]")
    for stage, n in funnel:
        print(f"  {stage:40s} : {n}")

    survivors = sorted(survivors)
    with open(OUT, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["peptide"])
        for s in survivors:
            w.writerow([s])
    print(f"\n[完成] 官方过滤候选集 {len(survivors)} 条 -> {OUT}")


if __name__ == "__main__":
    main()
