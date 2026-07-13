# Moso Bamboo (Phyllostachys edulis) DPP4-Inhibitory Peptides

> 基于 *in silico* 虚拟酶解 + 分子对接发现毛竹源 DPP4 抑制肽的计算项目。
> 方法学参考模板论文：Cheng et al., *Bioorganic Chemistry* 175 (2026) 109801 —
> 《Discovery of garlic-derived peptides as natural DPP4 inhibitors》.

---

## 1. 项目背景

原山药（*Dioscorea polystachya*）项目因 UniProt 仅 20 条 curated 蛋白、虚拟消化仅产 237 条肽（崩溃于原 5,230 漏斗），且 7 条实验验证肽（FWPQY 等）0/7 可由该 20 条蛋白切出，**山药无法支撑 DPP4 声明**。

遂参照上述模板论文方法，重选 **毛竹（Moso bamboo, *Phyllostachys edulis*）** 为研究对象。经 UniProt 核实，毛竹条目池（taxonomy 38705）共 **253 条**蛋白（全 TrEMBL 预测，0 条 Swiss-Prot 人工审阅），约为模板大蒜（113 条）的 2.2 倍，足以撑起完整管线。

---

## 2. 计算管线与结果

| 步骤 | 方法 | 主要结果 | 软件 |
|---|---|---|---|
| 1. 物种/蛋白池核实 | UniProt taxonomy/search + uniprotkb/search | *P. edulis* taxid **38705**；253 条（0 reviewed / 253 all）；无参考蛋白质组 | `curl` + UniProt REST API |
| 2. 蛋白下载 | UniProtKB `/stream` (FASTA) | `data/moso_253.fasta`（253 条，103,409 残基，均 409 aa） | `curl` |
| 3. 虚拟胃肠消化 | 复刻 ExPASy PeptideCutter：pepsin(pH1.3)+trypsin+特异 chymotrypsin（切 F/Y/W 的 C 端，Pro 前不停），去重保留 2–6 aa | 严格规则唯一肽 **7,988**；2–6 aa 短肽 **4,950**（大蒜模板 1,442 的 3.4 倍） | Python（numpy/scipy/biopython） |
| 4. PeptideRanker 风格打分 + 去过敏/毒 | 代理启发式（N 端疏水、Pro 含量、分子量）→ AllerTOP/ToxinPred 风格过滤 | >0.5：**4,333** → 过滤后 **4,289** | Python（自写） |
| 5. DPP4 结构偏好收窄 | 3–5 aa、N 端疏水、第 2 位 Pro/Ala 优先 | 候选池 **2,019** → 对接队列 **60** | Python（自写） |
| 6. 配体 3D 制备 (PDBQT) | SMILES→3D→PDBQT；8 条含 Arg/His 肽 NaN 用 RDKit 补救 | **60/60** 配体 PDBQT，零 NaN | OpenBabel / RDKit |
| 7. 受体与口袋 | DPP4 晶体 **1WCY**（西格列汀结合态）→ PDBQT；实测配体坐标定网格 | `docking/1WCY_receptor.pdbqt`（12,248 原子）；盒中心 **(62.8, 47.7, 4.8)** size 30³ | RCSB PDB / OpenBabel / awk |
| 8. 分子对接 (AutoDock Vina) | `--exhaustiveness 4 --cpu 2`，输出 9 构象/肽 | **60/60** 完成；Top：**LPPQ -7.472**、APSPE -7.150、LAPSP -7.087 kcal/mol | AutoDock Vina 1.2.5 |

### 对接 Top 10（结合自由能 dG, kcal/mol）

| 排名 | 肽 | dG | 类别 |
|---|---|---|---|
| 1 | LPPQ (Leu-Pro-Pro-Gln) | **-7.472** | 中-强 |
| 2 | APSPE (Ala-Pro-Ser-Pro-Glu) | -7.150 | 中-强 |
| 3 | LAPSP (Leu-Ala-Pro-Ser-Pro) | -7.087 | 中-强 |
| 4 | LPGF (Leu-Pro-Gly-Pro) | -7.075 | 中-强 |
| 5 | LPINP (Leu-Pro-Ile-Asn-Pro) | -6.988 | 中-强 |
| 6 | LPSP (Leu-Pro-Ser-Pro) | -6.867 | 中-强 |
| 7 | LPCPR (Leu-Pro-Cys-Pro-Arg) | -6.835 | 中-强 |
| 8 | LPGDP (Leu-Pro-Gly-Asp-Pro) | -6.793 | 中-强 |
| 9 | LPDDP (Leu-Pro-Asp-Asp-Pro) | -6.693 | 中-强 |
| 10 | APSQP (Ala-Pro-Ser-Gln-Pro) | -6.515 | 中-强 |

分布：中-强（-6.5 ~ -8）10 条 / 中等（-5 ~ -6.5）49 条 / 弱（>-5）1 条（CPPSK -4.856）；无 < -8.0 强结合。

---

## 3. 目录结构

```
moso-bamboo-DPP4-peptides/
├── README.md
├── .gitignore
├── data/                         # 蛋白与肽库
│   ├── moso_253.fasta           # 253 条毛竹蛋白序列
│   ├── moso_253_peptides.txt           # 宽松规则唯一肽
│   ├── moso_253_peptides_strict.txt   # 严格规则唯一肽
│   ├── moso_candidates_pr_filtered.txt# 2,019 候选肽
│   └── moso_dock_queue.txt     # 60 条对接队列
├── docking/
│   ├── 1WCY_receptor.pdbqt    # DPP4 受体（西格列汀结合态）
│   ├── moso_box.txt             # Vina 口袋参数
│   ├── moso_ligands/           # 60 配体 PDBQT + 60 对接输出构象
│   ├── moso_dock_results.tsv   # 60 肽对接 dG 结果表
│   └── moso_dock_ranking.txt  # 完整排名
├── scripts/                     # 可复现管线脚本
│   ├── rerun_digestion_moso253_strict.py
│   ├── moso_pipeline_filter2.py
│   ├── moso_build_ligand_pdbqts.py
│   ├── moso_dock_prepare_receptor.py
│   ├── moso_dock_run.py
│   ├── moso_report_results.py
│   ├── batch_dock.cmd / batch_dock.sh   # 一键批量对接
│   └── aux/                   # 报告生成脚本
└── docs/
    └── 毛竹DPP4抑制肽_项目总结.docx   # 项目完整总结
```

---

## 4. 复现方法

### 4.1 已完成（本仓库即结果）
虚拟消化 → 过滤 → 对接 全流程已跑通，结果见 `data/` 与 `docking/`。

### 4.2 重跑对接（需本地安装 Vina）
```bash
# 安装（本机/服务器）
conda install -c conda-forge vina openbabel rdkit

# 进入 scripts/，运行一键批量对接
# Windows:
scripts\batch_dock.cmd
# Linux/WSL:
bash scripts/batch_dock.sh
```
> 注：`vina.exe` 为第三方二进制，未纳入版本库。请从
> https://github.com/ccsb-scripps/AutoDock-Vina/releases/tag/v1.2.5
> 下载并置于运行目录，或修改 `batch_dock.*` 中的 vina 路径。

---

## 5. 必须如实披露的局限

1. **过滤打分是代理启发式**（基于文献理化特征），非官方 PeptideRanker / AllerTOP / ToxinPred 服务器输出，正式稿件须替换并注明。
2. **对接为静态 dG 估计**，需经 GROMACS MD / MM-PBSA + 体外实验确认，不能直接作为活性结论。
3. **毛竹 0 条人工审阅条目（全 TrEMBL 预测）**，Methods 须如实披露来源与预测性质（模板大蒜论文同样约 95% 为 TrEMBL，可接受）。
4. **WPHY/WPQY/VAPGW 是大蒜肽**，不应用于毛竹真值校验。

---

## 6. 后续待办（按模板全程）

- [ ] 官方 PeptideRanker / AllerTOP / ToxinPred 服务器验证 Top 候选
- [ ] GROMACS MD（50→150 ns）+ MM/PBSA 精算 Top 肽结合自由能
- [ ] 网络药理学：SwissTargetPrediction → STRING → DAVID
- [ ] 体外活性验证：Gly-Pro-pNA DPP4 抑制（IC₅₀）+ Caco-2 原位抑制

---

## 7. 引用

- 模板方法：Cheng Y. et al. *Discovery of garlic-derived peptides as natural DPP4 inhibitors.* Bioorganic Chemistry 175 (2026) 109801.
- DPP4 结构：PDB **1WCY**（sitagliptin-bound DPP4）.
- 蛋白序列：UniProtKB，*Phyllostachys edulis* (taxonomy 38705).
