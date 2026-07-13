# Phase A — 用官方服务器替代代理过滤

## 目的
此前 `moso_candidates_pr_filtered.txt` 中的 PeptideRanker 类打分是**代理启发式**（分值全为 1.000，无区分度），
不能作为计算论文的筛选依据。本阶段用官方在线服务器对毛竹虚拟消化得到的 **4,950 条 2–6aa 短肽** 做真实预测，
得到可发表的"官方过滤"候选集，替换代理结果。

## 输入（位于 `data/phaseA_inputs/`）
| 文件 | 说明 |
|---|---|
| `moso_short_2to6.fasta` | 主输入，全部 4,950 条，带 `>pep_00001` 式 ID |
| `moso_short_2to6.txt`   | 简单列表，每行一条 |
| `peptideranker/batch_XX.fasta` | 每批 500，共 10 批（手动提交用） |

## 三阶段过滤漏斗

```
母集 (4,950)
  └─ [可选] PeptideRanker ≥ 0.5  — 服务器暂不可用
  └─ AlgPred 2.0 ML Score < 0.6  — ✅ 自动提交
  └─ ToxinPred3 Non-Toxin         — ✅ 自动提交
  └─ official_candidates.tsv
```

## 运行方式

### 1️⃣ ToxinPred（自动 — 约需 5–10 分钟）
```bash
python scripts/phaseA/phaseA_run_toxinpred.py
```
- 自动将 4,950 条按每批 **50 条**提交到 ToxinPred 批量服务器
- 解析 HTML 结果表，提取 SVM Score + Prediction
- 输出 → `data/phaseA_inputs/results_toxinpred.csv`

### 2️⃣ AlgPred 2.0（自动 — 约需 5–10 分钟）
```bash
python scripts/phaseA/phaseA_run_algpred.py
```
- 自动提交到 AlgPred 2.0 Batch 服务器（AAC-RF 模式）
- 解析结果表，提取 ML Score + Prediction
- 输出 → `data/phaseA_inputs/results_algpred.csv`

### 3️⃣ PeptideRanker（手动 — 服务器暂不可用）
- URL: `https://peptide.ucd.ie/peptideranker/` 👈 当前返回 **503 Service Unavailable**
- 恢复后：打开 URL，粘贴 `data/phaseA_inputs/moso_short_2to6.fasta`
  （或逐个 `data/phaseA_inputs/peptideranker/batch_XX.fasta`）
- 将结果另存为 `data/phaseA_inputs/results_peptideranker.tsv`

### 4️⃣ 合并所有结果
```bash
python scripts/phaseA/phaseA_merge.py
```
脚本输出**漏斗计数**与 `official_candidates.tsv`。任一结果文件暂缺时自动跳过。

## 服务器详情
| 服务器 | URL | 返回 | 通过条件 |
|---|---|---|---|
| **PeptideRanker** (UCD) | https://peptide.ucd.ie/peptideranker/ ⚠ 503 | 0–1 生物活性概率 | ≥ 0.5 |
| **AlgPred 2.0** | https://webs.iiitd.edu.in/raghava/algpred2/batch.html | ML Score + Allergen/Non-allergen | ML Score < 0.6 |
| **ToxinPred3 Batch** | http://www.raghavagps.net/raghava/toxinpred/multi_submit.php | SVM Score + Toxin/Non-Toxin | Non-Toxin |

## 配置参数（`scripts/phaseA/phaseA_merge.py` 顶部 CONFIG 区）
- `THR_ALGPRED = 0.6` — AlgPred ML Score 阈值（低于此值视为非过敏原）
- `TOXIN_NEG_LABELS` — ToxinPred 阴性标签集合
- `THR_PR = 0.5` — PeptideRanker 阈值（待服务器恢复后启用）

## 诚信声明（写入稿件 Methods）
- ToxinPred 与 AlgPred 的预测经由官方服务器 HTTP 批量提交完成；
- PeptideRanker 需手动提交（服务器当前不可用）；
- 毛竹蛋白组为 UniProt TrEMBL 预测条目（0 条 Swiss-Prot 人工审阅），Methods 如实披露。
