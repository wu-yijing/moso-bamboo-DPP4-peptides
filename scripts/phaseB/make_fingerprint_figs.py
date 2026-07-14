#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
由 contact_fingerprint_LPPGP_APPSQ.json 与 phaseB_detail.json 生成接触指纹可视化:
  figures/contact_fingerprint_LPPGP_APPSQ.svg  (残基层面 barcode 指纹, 核心)
  figures/contact_per_position.svg             (逐配体位置接触剖面)
  figures/contact_fingerprint.html             (合并浏览器视图)
"""
import os, json

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FIG  = os.path.join(ROOT, "figures")
os.makedirs(FIG, exist_ok=True)
DATA = os.path.join(ROOT, "data", "phaseB")

with open(os.path.join(DATA, "contact_fingerprint_LPPGP_APPSQ.json")) as f:
    FP = json.load(f)
with open(os.path.join(DATA, "phaseB_detail.json")) as f:
    DET = json.load(f)

ROLE_COLOR = {
    "catalytic": "#c0392b",   # 催化核心 (红)
    "S1S2":      "#e67e22",   # S1/S2 口袋 (橙)
    "glycan":    "#7f8c8d",   # 糖基化偏倚 (灰)
    "interface": "#2980b9",   # 界面/口袋壁 (蓝)
    "pocket":    "#16a085",   # 其他口袋 (青)
}
ROLE_LABEL = {
    "catalytic": "催化核心 (Glu146/Arg147)",
    "S1S2":      "S1/S2 口袋",
    "glycan":    "糖基化偏倚 (NAG/FUC)",
    "interface": "界面/口袋壁",
    "pocket":    "口袋残基",
}

def lerp(a, b, t):
    return tuple(int(round(a[i] + (b[i]-a[i])*t)) for i in range(3))
# 接触强度配色: 浅黄 -> 深红
C0 = (255, 247, 236); C1 = (179, 0, 0)
def intensity_color(n, nmax):
    if nmax <= 0: return "#ffffff"
    t = min(1.0, n / nmax)
    r, g, b = lerp(C0, C1, t)
    return f"#{r:02x}{g:02x}{b:02x}"

# ---------------- 残基层指纹 barcode ----------------
rows = FP["rows"]
nmax = max(max(r["LPPGP_contacts"], r["APPSQ_contacts"]) for r in rows)
colw = 46
left = 70; top = 60; lane_h = 56; lane_gap = 26
W = max(560, left + len(rows)*colw + 30)
H = top + 2*lane_h + lane_gap + 190

svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
       f'font-family="DejaVu Sans, Arial, sans-serif" font-size="11">']
svg.append(f'<rect width="{W}" height="{H}" fill="#ffffff"/>')
svg.append(f'<text x="{left}" y="26" font-size="15" font-weight="bold">LPPGP / APPSQ 残基层接触指纹 (最小化复合体, ≤4.5 Å)</text>')
svg.append(f'<text x="{left}" y="44" font-size="11" fill="#555">细胞填充=接触数 (浅→深); ★=极性/H-bond; 边框色=残基角色</text>')

def lane_y(pep_i):
    return top + pep_i*(lane_h + lane_gap)
labels = ["LPPGP", "APPSQ"]
for pi in (0, 1):
    y0 = lane_y(pi)
    svg.append(f'<rect x="{left-6}" y="{y0}" width="{len(rows)*colw+4}" height="{lane_h}" fill="#fafafa" stroke="#ddd"/>')
    svg.append(f'<text x="{left-12}" y="{y0+lane_h/2+4}" text-anchor="end" font-weight="bold">{labels[pi]}</text>')

for ci, r in enumerate(rows):
    x = left + ci*colw
    # 列分隔
    svg.append(f'<line x1="{x}" y1="{top}" x2="{x}" y2="{top+2*lane_h+lane_gap}" stroke="#eee"/>')
    vals = [r["LPPGP_contacts"], r["APPSQ_contacts"]]
    hbs  = [r["LPPGP_hbond"]=="Y", r["APPSQ_hbond"]=="Y"]
    for pi in (0,1):
        y0 = lane_y(pi)
        v = vals[pi]
        cx = x + colw/2; cy = y0 + lane_h/2
        if v > 0:
            rad = 6 + 13*(v/nmax)
            col = intensity_color(v, nmax)
            rc = ROLE_COLOR.get(r["role"], "#16a085")
            svg.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{rad:.1f}" fill="{col}" '
                       f'stroke="{rc}" stroke-width="2"/>')
            svg.append(f'<text x="{cx:.1f}" y="{cy+4:.1f}" text-anchor="middle" font-size="10" '
                       f'fill="#222">{v}</text>')
            if hbs[pi]:
                svg.append(f'<text x="{cx:.1f}" y="{y0+10:.1f}" text-anchor="middle" font-size="12" '
                           f'fill="#c0392b" font-weight="bold">★</text>')
    # x 标签 (残基短名)
    short = r["res"]
    svg.append(f'<text x="{x+colw/2:.1f}" y="{top+2*lane_h+lane_gap+14}" text-anchor="middle" '
               f'font-size="9" fill="{ROLE_COLOR.get(r["role"],"#16a085")}">{short}</text>')
# 斜线引导
svg.append(f'<line x1="{left}" y1="{top+2*lane_h+lane_gap}" x2="{left}" y2="{top+2*lane_h+lane_gap+4}" stroke="#999"/>')

# 图例 (竖排, 避免超出画布)
ly = top + 2*lane_h + lane_gap + 40
svg.append(f'<text x="{left}" y="{ly}" font-weight="bold">残基角色:</text>')
# 两列布局
items = list(ROLE_COLOR.items())
col_x = [left, left + 250]
for idx, (role, col) in enumerate(items):
    cx = col_x[idx % 2]
    cy = ly + 20 + (idx // 2) * 20
    svg.append(f'<rect x="{cx}" y="{cy-10}" width="12" height="12" fill="#fff" stroke="{col}" stroke-width="2"/>')
    svg.append(f'<text x="{cx+16}" y="{cy}" font-size="10">{ROLE_LABEL[role]}</text>')
foot_y = ly + 20 + ((len(items)+1)//2) * 20 + 8
svg.append(f'<text x="{left}" y="{foot_y}" font-size="10" fill="#777">★ 极性/H-bond 接触; 圆面积按接触原子对数着色 (浅→深)。</text>')
svg.append('</svg>')

svg_path = os.path.join(FIG, "contact_fingerprint_LPPGP_APPSQ.svg")
with open(svg_path, "w") as f:
    f.write("\n".join(svg))
print(f"[OK] {svg_path}")

# ---------------- 逐位置接触剖面 (bar) ----------------
seqs = {d["peptide"]: d["seq"] for d in DET}
ppos = {d["peptide"]: {p[0]: p[2] for p in d["contact_per_position"]} for d in DET}
targets = ["LPPGP", "APPSQ"]
barw = 34; group_gap = 26; left2 = 60; top2 = 60
W2 = left2 + len(targets)*5*barw + len(targets)*group_gap + 40
H2 = 320
svg2 = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W2}" height="{H2}" '
        f'font-family="DejaVu Sans, Arial, sans-serif" font-size="11">']
svg2.append(f'<rect width="{W2}" height="{H2}" fill="#fff"/>')
svg2.append(f'<text x="{left2}" y="26" font-size="15" font-weight="bold">逐配体位置接触数 (最小化复合体)</text>')
ymax = max(max(ppos[t].values()) for t in targets)
y_top = top2; y_bot = H2-50
def ymap(v): return y_bot - (v/ymax)*(y_bot-y_top)
svg2.append(f'<line x1="{left2}" y1="{y_bot}" x2="{W2-20}" y2="{y_bot}" stroke="#333"/>')
for v in range(0, ymax+1, max(1, ymax//4)):
    yy = ymap(v)
    svg2.append(f'<line x1="{left2}" y1="{yy:.1f}" x2="{W2-20}" y2="{yy:.1f}" stroke="#eee"/>')
    svg2.append(f'<text x="{left2-6}" y="{yy+4:.1f}" text-anchor="end" font-size="9" fill="#777">{v}</text>')
tcol = {"LPPGP": "#27ae60", "APPSQ": "#8e44ad"}
for ti, t in enumerate(targets):
    seq = seqs[t]; gx = left2 + ti*(5*barw + group_gap)
    svg2.append(f'<text x="{gx+5*barw/2:.1f}" y="{y_bot+18}" text-anchor="middle" font-weight="bold" fill="{tcol[t]}">{t}</text>')
    for pi0 in range(5):
        pos = pi0+1
        v = ppos[t].get(pos, 0)
        bx = gx + pi0*barw
        by = ymap(v)
        svg2.append(f'<rect x="{bx+2}" y="{by:.1f}" width="{barw-4}" height="{y_bot-by:.1f}" '
                    f'fill="{tcol[t]}" opacity="0.85"/>')
        svg2.append(f'<text x="{bx+barw/2:.1f}" y="{by-4:.1f}" text-anchor="middle" font-size="9">{v}</text>')
        # 位置字母 (seq)
        svg2.append(f'<text x="{bx+barw/2:.1f}" y="{y_bot+32}" text-anchor="middle" font-size="10" font-weight="bold">{seq[pi0]}{pos}</text>')
svg2.append('</svg>')
svg2_path = os.path.join(FIG, "contact_per_position.svg")
with open(svg2_path, "w") as f:
    f.write("\n".join(svg2))
print(f"[OK] {svg2_path}")

# ---------------- HTML 合并视图 ----------------
html = f"""<!DOCTYPE html><html lang="zh"><head><meta charset="utf-8">
<title>LPPGP / APPSQ 接触指纹</title>
<style>body{{font-family:system-ui,Arial,sans-serif;margin:24px;color:#222}}
h1{{font-size:20px}} .fig{{border:1px solid #ddd;border-radius:8px;padding:12px;margin:16px 0}}
.cap{{color:#555;font-size:13px;margin-top:8px}}</style></head>
<body><h1>LPPGP / APPSQ 接触指纹可视化</h1>
<div class="fig"><img src="contact_fingerprint_LPPGP_APPSQ.svg" width="100%">
<div class="cap">图1. 残基层接触指纹（最小化复合体几何）。每一列为一个与任一肽发生接触的受体残基，
按残基号排序；上/下两条分别为 LPPGP 与 APPSQ。圆面积∝接触原子对数，颜色深浅同；★标记极性/H-bond 接触；
外框颜色=残基在 DPP4 口袋中的角色（红=催化核心 Glu146/Arg147，橙=S1/S2 口袋，灰=糖基化偏倚，蓝=界面/口袋壁）。</div></div>
<div class="fig"><img src="contact_per_position.svg" width="70%">
<div class="cap">图2. 逐配体位置接触数。横轴为肽的 5 个残基位置（字母+序号），柱高=该位置与口袋的接触原子对数。
LPPGP 接触集中于 Gly4 + Pro2；APPSQ 集中于 Ser4 + Pro3，足迹更广。</div></div>
<p style="color:#777;font-size:12px">接触判据：总接触 ≤4.5 Å；疏水（双方碳）≤4.0 Å；极性/H-bond（双方 O/N）≤3.5 Å。
几何来自 MMFF94s 松弛后的复合体（与 §3.6/§3.7 一致）。</p>
</body></html>"""
html_path = os.path.join(FIG, "contact_fingerprint.html")
with open(html_path, "w") as f:
    f.write(html)
print(f"[OK] {html_path}")
