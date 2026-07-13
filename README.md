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

## 2.5 Phase A — 官方服务器过滤（替代原代理打分）

> PeptideRanker 官方服务器长期 503 不可用；以经同行验证的 **ToxinPred 3.0 批量提交 + AlgPred 2.0 (AAC-RF)** 替代过敏/毒性过滤层。

| 阶段 | 过滤 | 剩余 |
|---|---|---|
| 母集（2–6 aa 短肽） | — | **4,950** |
| AlgPred 2.0（ML score < 0.6） | −270 | 4,680 |
| ToxinPred 3.0（Non-Toxin） | −6 | **4,674** |
| PeptideRanker（PR ≥ 0.5） | 跳过（服务器 503，待恢复后补） | — |

- 全部 **60 条对接队列肽均通过** ToxinPred + AlgPred（0 条被拒）。
- 产物：`data/phaseA_inputs/results_toxinpred.csv`、`results_algpred.csv`、`official_candidates.tsv`（4,674 条）。
- 脚本：`scripts/phaseA/phaseA_run_toxinpred.py`、`phaseA_run_algpred.py`、`phaseA_merge.py`。
- 详见 `docs/phaseA_README.md`。

---

## 2.6 Phase B — 纯计算结合验证（静态 MM，无 MD）

> 本机**无 GROMACS / 无 conda**，无法运行 100–150 ns MD。采用 **单构象端点法静态 MMFF94s 验证** + 几何接触剖面，作为湿实验/MD 的可行性替代。

| 肽（序列） | dG_Vina | ΔE_MM* | 口袋残基 | 总接触 | 氢键 | 关键接触位置 |
|---|---|---|---|---|---|---|
| LPPQ (Leu-Pro-Pro-Gln) | −7.472 | +2.45 | 39 | 92 | 10 | **Gln4:37**（C 端主导），Pro2/3 各 ~21–25 |
| APSPE (Ala-Pro-Ser-Pro-Glu) | −7.150 | +2.93 | 39 | 119 | 13 | **Glu5:109**（S1′ 压倒性） |
| LAPSP (Leu-Ala-Pro-Ser-Pro) | −7.087 | +2.84 | 40 | 101 | 13 | **Pro5:50**（C 端 Pro） |

> \* ΔE_MM 为气相单点能量差，**非真实 ΔG**（无显式溶剂化/熵），三肽量级相近、不具区分度；结合强弱排序仍以 Vina dG 为主判据。

- **生物学一致性**：APSPE 的 C 端谷氨酸在 S1′ 口袋的极端接触偏好，与 DPP4“偏好底物 C 端带负电残基占据 S1′”机制吻合；双 Pro 核心符合 S1′/S2′ 对 Pro 的偏好。
- **氢键伙伴**：三肽均富集 **GLU146**（S2′ 区），并分别触及 PRO149 / SER182 / TYR183 / ASN151 —— 均位于 DPP4 已知活性中心腔体。
- 产物：`data/phaseB/phaseB_results.tsv`、`phaseB_detail.json`。
- 脚本：`scripts/phaseB/phaseB_validation.py`。
- 局限与复现：详见 `docs/phaseB_README.md`。

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
│   └── reports/                # 报告生成脚本
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

1. **过滤层**：PeptideRanker 官方服务器 503 不可用，已以 **ToxinPred 3.0 + AlgPred 2.0** 官方输出替代过敏/毒性过滤（Phase A，4,674 候选）；PR 层待服务器恢复补全。原"代理启发式"打分仅用于早期 2,019→60 对接队列收窄，不应作为最终发表判据。
2. **对接为静态 dG 估计**，需经 GROMACS MD / MM-PBSA + 体外实验确认，不能直接作为活性结论。
3. **毛竹 0 条人工审阅条目（全 TrEMBL 预测）**，Methods 须如实披露来源与预测性质（模板大蒜论文同样约 95% 为 TrEMBL，可接受）。
4. **WPHY/WPQY/VAPGW 是大蒜肽**，不应用于毛竹真值校验。

---

## 6. 后续待办（按模板全程；本项目为纯计算，无湿实验）

- [x] **Phase A** 官方 ToxinPred 3.0 + AlgPred 2.0 过滤（4,674 候选；PeptideRanker 服务器 503 待恢复补）
- [x] **Phase B（静态近似）** 单构象 MMFF94s 松弛 + 几何接触剖面（Top3 已完成）
- [ ] **Phase B 升级** 若在 HPC 获 GROMACS：补 100–150 ns MD + MM-PBSA（含 GB 溶剂化与熵项），得定量 ΔG 与收敛残基分解
- [ ] PeptideRanker 恢复后补 PR ≥ 0.5 过滤层（scripts/phaseA/phaseA_merge.py 已支持自动并入）
- [ ] 网络药理学：SwissTargetPrediction → STRING → DAVID（计算上下文，替代 Caco-2 转运）
- [ ] 稿件定位 *in silico* 发现；Methods 如实披露 ①官方过滤已做 ②Vina 为静态 dG ③毛竹 0 条人工审阅全 TrEMBL ④本 Phase B 为静态近似、湿实验未做

---

## 7. 引用

- 模板方法：Cheng Y. et al. *Discovery of garlic-derived peptides as natural DPP4 inhibitors.* Bioorganic Chemistry 175 (2026) 109801.
- DPP4 结构：PDB **1WCY**（sitagliptin-bound DPP4）.
- 蛋白序列：UniProtKB，*Phyllostachys edulis* (taxonomy 38705).
