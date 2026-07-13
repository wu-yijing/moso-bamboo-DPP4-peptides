# Phase A — 用官方服务器替代代理过滤（PeptideRanker / AllerTOP v2 / ToxinPred3）

## 目的
此前 `moso_candidates_pr_filtered.txt` 中的 PeptideRanker 类打分是**代理启发式**（分值全为 1.000，无区分度），
不能作为计算论文的筛选依据。本阶段用三个官方在线服务器对毛竹虚拟消化得到的 **4,950 条 2–6aa 短肽** 做真实预测，
得到可发表的"官方过滤"候选集，替换代理结果。

## 输入（已生成，位于 `data/phaseA_inputs/`）
| 文件 | 说明 |
|---|---|
| `moso_short_2to6.fasta` | 主输入，全部 4,950 条，带 `>pep_00001` 式 ID |
| `moso_short_2to6.txt`   | 简单列表，每行一条（部分服务器用） |
| `peptideranker/batch_XX.fasta` | 每批 500，共 10 批 |
| `allertop/batch_XX.fasta`     | 同上 |
| `toxinpred/batch_XX.fasta`    | 同上 |

> 想一次提交全部：直接用 `moso_short_2to6.fasta`；若服务器报错/超时，改传 `batch_XX.fasta`。

## 三个服务器
| 服务器 | URL | 返回 | 通过阈值 |
|---|---|---|---|
| **PeptideRanker** (UCD) | http://bioware.ucd.ie/~compass/biowareweb/Server_pages/ （或 http://distilldeep.ucd.ie/PeptideRanker/） | 每条肽的生物活性概率 0–1，按概率排序 | **≥ 0.5** 为有活性 |
| **AllerTOP v2** | https://www.ddg-pharmfac.net/AllerTOP/ | 过敏原预测（Allergen / Non-allergen）+ 评分 | 保留 **Non-allergen** |
| **ToxinPred3** | https://webs.iiitd.edu.in/raghava/toxinpred3/ （Batch submission） | 毒性预测（Toxin / Non-Toxin）+ 理化性质 | 保留 **Non-Toxin** |

> 诚实注记：ToxinPred3 有命令行版（`pip install toxinpred3`），但**强制依赖 scikit-learn==1.0.2**，
> 而本机仅 Python 3.13，该旧版 sklearn 无法在 3.13 上安装/反序列化模型。故 ToxinPred 同样走网页批量提交。
> 三个服务器均为网页交互，无法由本机脚本自动调用——以下为人工提交步骤。

## 逐步提交
### 1. PeptideRanker
1. 打开上面的 URL，粘贴 `moso_short_2to6.fasta`（或逐个 `batch_XX.fasta`）。
2. 选择 short-peptide 预测器（肽长 <20aa 会自动用短肽模型）。
3. 提交，等待返回**按概率排序**的列表。
4. 将结果另存为 `data/phaseA_inputs/results_peptideranker.tsv`（列含：序列 + 概率分）。

### 2. AllerTOP v2
1. 打开 URL，提交同样的 FASTA（或批次）。
2. 取回过敏原判定表，另存为 `data/phaseA_inputs/results_allertop.tsv`
   （列含：序列 + Allergen/Non-allergen 标签）。

### 3. ToxinPred3（Batch）
1. 打开 URL 的 **Batch submission** 页，提交同样的 FASTA。
2. 取回结果 CSV，另存为 `data/phaseA_inputs/results_toxinpred.csv`
   （默认列：Sequence ID, Sequence, ML Score, Hybrid Score, Prediction, PPV）。

## 合并（回填后自动过滤）
结果就位后，编辑 `scripts/phaseA/phaseA_merge.py` 顶部的 `CONFIG` 区：
- 确认三个结果文件路径正确；
- 若服务器返回列名与自动探测关键词不符，显式填写 `*_COL` 列名。

然后运行：
```bash
python scripts/phaseA/phaseA_merge.py
```
脚本输出**漏斗计数**（母集 → PR 通过 → 非过敏原 → Non-Toxin）与
`data/phaseA_inputs/official_candidates.tsv`（官方过滤候选集，供后续对接队列重排）。

> 任一结果文件暂缺时，该阶段会被跳过（标注 skipped），可分批回填后重跑。

## 与下游的衔接
官方候选集将**替换** `moso_candidates_pr_filtered.txt` 与 `moso_dock_queue.txt` 中的 60 条对接队列来源：
重新按官方分数排序，取 Top 60（或更多）进入 AutoDock Vina 对接（Phase B 之前的数据准备）。
当前 60 条对接结果（LPPQ −7.472 等）暂不废止，待官方候选集生成后做一致性比对。

## 诚信声明（写入稿件 Methods）
- 短肽生物活性/过敏原/毒性预测分别由 PeptideRanker、AllerTOP v2、ToxinPred3 官方服务器完成；
- 原始代理打分（PeptideRanker-style、AllerTOP/ToxinPred 代理）已弃用；
- 毛竹蛋白组为 UniProt TrEMBL 预测条目（0 条 Swiss-Prot 人工审阅），Methods 如实披露。
