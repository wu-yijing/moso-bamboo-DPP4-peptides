# Phase C — In silico ADMET & DPP4-centred Network Pharmacology

**项目定位**：纯计算研究（no wet-lab）。本阶段是 Phase B 静态 MM 验证之后的
**第二轮纯计算验证层**，替代原本依赖体外 ADMET / Caco-2 渗透 / 血清稳定性
实验的功能性验证。全部产物均为计算推导，无任何湿实验。

---

## 1. 目标

1. **肽级 in silico ADMET 概貌**（C1）：对 60 条 Vina 对接肽计算
   物化与类药描述符，并基于蛋白酶裂解位点规则预测胃肠道（GI）稳定性。
2. **DPP4 中心网络药理学**（C2）：通过 STRING DB 获取 DPP4 的人类
   蛋白互作/功能关联网络，叠加文献策展的 DPP4 底物层，构建药理学知识图谱，
   解释"抑制 DPP4 → 升高 GLP-1/GIP → 改善血糖"的机制语境。

---

## 2. 方法

### 2.1 C1 — 肽描述符（纯计算）
- **序列衍生描述符**（BioPython `ProtParam`，v1.87）：分子量（MW）、
  等电点（pI）、pH 7.4 净电荷、GRAVY（平均疏水性）、不稳定指数、
  脂肪族指数（按 `AI = A + 2.9·V + 3.9·(I+L)`，X 为摩百分比）、
  pH 7.4 电荷。
- **Boman 蛋白结合指数**（Boman 1995）：`Boman = (1/n)·ΣΔf`，
  Δf 为氨基酸水→环己烷转移自由能。**负值 → 更高的蛋白结合倾向**。
- **RDKit 分子描述符**：ExactMolWt、TPSA、HBD（NumHDonors）、
  HBA（NumHAcceptors）、可旋转键、QED 类药性。肽由 `Chem.MolFromSequence`
  构建后计算（2D 描述符，无需构象采样）。

### 2.2 C1 — GI 道稳定性（蛋白酶裂解规则）
对每条肽统计内切酶裂解位点（简化特异性模型）：
- **胰蛋白酶 trypsin**：K/R 后（除非后是 P）
- **胰凝乳蛋白酶 chymotrypsin**：F/Y/W 后
- **胃蛋白酶 pepsin**：F/Y/W/L/M 后（广谱疏水）
- **弹性蛋白酶 elastase**：A/G/S/V 后
- **羧肽酶 / 氨肽酶**：C/N 端外切（N 端 Pro 封端可抗氨肽酶，C 端 Pro 抗羧肽酶）
- **DPP4 自身裂解**：N 端第 2 位为 Pro/Ala（X-Pro / X-Ala 基序）→ 该肽本身为 DPP4 底物

GI 稳定性分级以**内切酶位点**为核心：
`0 位点且 N 端 Pro 封端 → Very High`；`0 → High`；`1 → Moderate`；`≥2 → Low`。

### 2.3 C2 — DPP4 网络（STRING + 文献）
- **STRING DB REST API**（`network`，species=9606，functional，
  required_score=400）：获取 DPP4（UniProt P27487）功能关联网络
  → 边类型 `functional_string`（combined score 0–1）。
- **文献策展底物层**：叠加 DPP4 裂解灭活的肽类激素
  （GCG/GLP-1、GIP、CXCL12/SDF-1、NPY、PYY、TAC1/substance P）
  → 边类型 `cleavage_substrate_literature`（带引用标签）。
- 脚本：`scripts/phaseC/phaseC_network.py`

---

## 3. 结果

### 3.1 C1 — 60 条对接肽概貌
- 全部 60 条均携带 **X-Pro / X-Ala 基序**（第 1–2 位为 Pro 或 Ala），
  与候选筛选阶段"偏好 X-Pro"的 DPP4 结构偏好过滤一致——既支撑结合，
  也意味着它们在体内同样会被 DPP4 当作底物快速降解（见局限性）。
- **GI 稳定性分布**：Moderate（1 内切位点）39 条 / Low（多内切位点）11 条 /
  High（无内切位点）10 条。游离线性短肽普遍对外切酶敏感。
- **Top3 描述符**：
  | 肽 | dG (kcal/mol) | pI | Boman | GRAVY | GI 级 |
  |---|---|---|---|---|---|
  | LPPQ | −7.472 | 5.53 | +0.375 | −0.725 | Moderate |
  | APSPE | −7.150 | 4.60 | −0.820 | −1.140 | Low |
  | LAPSP | −7.087 | 5.53 | +0.340 | +0.320 | Low |
  - APSPE 的 Boman 指数最负（−0.82），提示相对更高的蛋白结合倾向。
- 完整逐肽数据：`data/phaseC/phaseC_peptides.tsv`

### 3.2 C2 — DPP4 中心网络
- **15 节点 / 32 边**（STRING 11 节点 26 功能边 + 文献 4 节点 6 底物边）。
- STRING 网络核心邻居：**GCG**（GLP-1/GLP-2 前体）、**GIP**（胃增泌素）、
  **CXCR4**（CXCL12/SDF-1 受体）、**ADA**（腺苷脱氨酶，DPP4/CD26 共价复合物）、
  **CAV1**（窖蛋白-1）、**PRCP**（脯氨酰羧肽酶，DPP4 家族同源酶）、
  **ACE2**、**FN1**、**ITGB1**、**PTPRC/CD45**。
- 机制语境：DPP4 通过 X-Pro/X-Ala 修剪灭活 GLP-1、GIP 等增泌素 →
  抑制 DPP4 可升高活性 GLP-1/GIP → 促进胰岛素分泌 → 改善 2 型糖尿病血糖。
- 产物：`data/phaseC/phaseC_network.json`（机读）、
  `data/phaseC/phaseC_network_summary.txt`（可读摘要）。

---

## 4. 诚实局限（投稿 Methods / Limitations 必披露）

1. **描述符为序列衍生，非实测 ADMET**：MW/pI/GRAVY/Boman/TPSA 等由
   序列/结构推算，未经任何体外渗透、代谢、血清稳定性实验验证。
2. **GI 稳定性为简化规则**：基于 published 蛋白酶特异性，非 ex-vivo 消化或
   动物体内稳定性测定；短肽对外切酶普遍敏感，结论仅为粗筛。
3. **DPP4 自身裂解张力**：全部 60 候选均为 X-Pro/X-Ala 基序（既是 DPP4
   偏好结合基序，也是其底物基序）。提示口服递送需结构稳定化
   （如环化、N 端封端），否则体内会被 DPP4 快速降解——属已知局限，
   不否定其体外/计算结合活性。
4. **STRING 网络为人类互作背景**：证明 DPP4 处于增泌素/糖代谢语境，
   **不等于**证明毛竹肽作用于这些节点；底物边为文献策展语境。
5. **不稳定指数对短肽参考价值有限**（该指标为全长蛋白设计），仅作记录。
6. 本 Phase C 为**计算近似**，无任何湿实验（Caco-2、血清稳定性、
   体外 DPP4 抑制 IC₅₀ 均未做）。

---

## 5. 复现

```bash
# C1 肽级 ADMET + GI 稳定性
python scripts/phaseC/phaseC_peptides.py
#   -> data/phaseC/phaseC_peptides.tsv

# C2 DPP4 网络（需联网拉取 STRING；亦可复用已存 raw_dpp4_network.json）
curl "https://string-db.org/api/json/network?identifiers=9606.ENSP00000353731&species=9606&required_score=400&network_type=functional" -o data/phaseC/dpp4_network.json
python scripts/phaseC/phaseC_network.py
#   -> data/phaseC/phaseC_network.json / phaseC_network_summary.txt
```

依赖：Python ≥3.10、`rdkit`、`biopython`、`numpy`、`requests`（STRING 拉取）。
