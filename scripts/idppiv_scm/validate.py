# -*- coding: utf-8 -*-
"""验证本地 iDPPIV-SCM 是否复现文献精度
   文献: CV≈0.819, 独立测试≈0.797 (Charoenkwan 2020)
"""
import sys, os, random
sys.path.insert(0, os.path.dirname(__file__))
from scm import load_tsv, build_scm, score, evaluate

random.seed(42)

tr_seqs, tr_lab = load_tsv("idppiv_scm/data/train.tsv")
te_seqs, te_lab = load_tsv("idppiv_scm/data/test.tsv")
print(f"训练集: {len(tr_seqs)} 条 (正={sum(tr_lab)}, 负={len(tr_lab)-sum(tr_lab)})")
print(f"测试集: {len(te_seqs)} 条 (正={sum(te_lab)}, 负={len(te_lab)-sum(te_lab)})")

pos = [s for s, y in zip(tr_seqs, tr_lab) if y == 1]
neg = [s for s, y in zip(tr_seqs, tr_lab) if y == 0]
P, L = build_scm(pos, neg)
print(f"SCM 评分卡维度: 位置数 L_max={L}, 每位置 {len(AA:=set('ACDEFGHIKLMNPQRSTVWY'))} 个氨基酸")

# ---- 在独立测试集上评估 ----
acc_te, mcc_te, _ = evaluate(te_seqs, te_lab, P, L)
print(f"\n[独立测试集] ACC={acc_te:.3f}  MCC={mcc_te:.3f}  (文献≈0.797)")

# ---- 5 折交叉验证 (训练集内) ----
idx = list(range(len(tr_seqs)))
random.shuffle(idx)
folds = [idx[i::5] for i in range(5)]
accs = []
for k in range(5):
    test_idx = set(folds[k])
    tr_idx = [i for i in idx if i not in test_idx]
    p = [tr_seqs[i] for i in tr_idx if tr_lab[i] == 1]
    n = [tr_seqs[i] for i in tr_idx if tr_lab[i] == 0]
    Pk, Lk = build_scm(p, n)
    fold_seqs = [tr_seqs[i] for i in folds[k]]
    fold_lab = [tr_lab[i] for i in folds[k]]
    a, _, _ = evaluate(fold_seqs, fold_lab, Pk, Lk)
    accs.append(a)
cv = sum(accs) / len(accs)
print(f"[5折CV] 平均 ACC={cv:.3f}  (各折: {[f'{a:.3f}' for a in accs]})  (文献≈0.819)")

# ---- 阈值敏感性 (在训练集上用全量模型走一遍, 仅供参考) ----
acc_tr, mcc_tr, _ = evaluate(tr_seqs, tr_lab, P, L)
print(f"[训练集全量] ACC={acc_tr:.3f}  MCC={mcc_tr:.3f}")

print("\n=> 若 独立测试 ACC 落在 0.75-0.82、CV 落在 0.78-0.85, 即视为忠实复现。")
