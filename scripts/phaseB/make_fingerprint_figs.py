#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate contact-fingerprint visualizations from contact_fingerprint_LPPGP_APPSQ.json
and phaseB_detail.json:
  figures/contact_fingerprint_LPPGP_APPSQ.svg  (residue-level barcode fingerprint, core)
  figures/contact_per_position.svg             (per-ligand-position contact profile)
  figures/contact_fingerprint.html             (combined browser view)
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
    "catalytic": "#c0392b",   # catalytic core (red)
    "S1S2":      "#e67e22",   # S1/S2 pocket (orange)
    "glycan":    "#7f8c8d",   # glycan bias (gray)
    "interface": "#2980b9",   # interface / pocket wall (blue)
    "pocket":    "#16a085",   # other pocket (teal)
}
ROLE_LABEL = {
    "catalytic": "catalytic core (Glu146/Arg147)",
    "S1S2":      "S1/S2 pocket",
    "glycan":    "glycan bias (NAG/FUC)",
    "interface": "interface / pocket wall",
    "pocket":    "pocket residue",
}

def lerp(a, b, t):
    return tuple(int(round(a[i] + (b[i]-a[i])*t)) for i in range(3))
# contact-intensity color ramp: pale yellow -> deep red
C0 = (255, 247, 236); C1 = (179, 0, 0)
def intensity_color(n, nmax):
    if nmax <= 0: return "#ffffff"
    t = min(1.0, n / nmax)
    r, g, b = lerp(C0, C1, t)
    return f"#{r:02x}{g:02x}{b:02x}"

# ---------------- residue-level fingerprint barcode ----------------
rows = FP["rows"]
nmax = max(max(r["LPPGP_contacts"], r["APPSQ_contacts"]) for r in rows)
colw = 46
left = 70; top = 60; lane_h = 56; lane_gap = 26
W = max(560, left + len(rows)*colw + 30)
H = top + 2*lane_h + lane_gap + 190

svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
       f'font-family="DejaVu Sans, Arial, sans-serif" font-size="11">']
svg.append(f'<rect width="{W}" height="{H}" fill="#ffffff"/>')
svg.append(f'<text x="{left}" y="26" font-size="15" font-weight="bold">LPPGP / APPSQ residue-level contact fingerprint (minimized complex, ≤4.5 Å)</text>')
svg.append(f'<text x="{left}" y="44" font-size="11" fill="#555">cell fill = contact count (pale→deep); ★ = polar/H-bond; border color = residue role</text>')

def lane_y(pep_i):
    return top + pep_i*(lane_h + lane_gap)
labels = ["LPPGP", "APPSQ"]
for pi in (0, 1):
    y0 = lane_y(pi)
    svg.append(f'<rect x="{left-6}" y="{y0}" width="{len(rows)*colw+4}" height="{lane_h}" fill="#fafafa" stroke="#ddd"/>')
    svg.append(f'<text x="{left-12}" y="{y0+lane_h/2+4}" text-anchor="end" font-weight="bold">{labels[pi]}</text>')

for ci, r in enumerate(rows):
    x = left + ci*colw
    # column divider
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
    # x-axis label (short residue name)
    short = r["res"]
    svg.append(f'<text x="{x+colw/2:.1f}" y="{top+2*lane_h+lane_gap+14}" text-anchor="middle" '
               f'font-size="9" fill="{ROLE_COLOR.get(r["role"],"#16a085")}">{short}</text>')
# slanted guide
svg.append(f'<line x1="{left}" y1="{top+2*lane_h+lane_gap}" x2="{left}" y2="{top+2*lane_h+lane_gap+4}" stroke="#999"/>')

# legend (vertical, to avoid exceeding canvas)
ly = top + 2*lane_h + lane_gap + 40
svg.append(f'<text x="{left}" y="{ly}" font-weight="bold">Residue role:</text>')
# two-column layout
items = list(ROLE_COLOR.items())
col_x = [left, left + 250]
for idx, (role, col) in enumerate(items):
    cx = col_x[idx % 2]
    cy = ly + 20 + (idx // 2) * 20
    svg.append(f'<rect x="{cx}" y="{cy-10}" width="12" height="12" fill="#fff" stroke="{col}" stroke-width="2"/>')
    svg.append(f'<text x="{cx+16}" y="{cy}" font-size="10">{ROLE_LABEL[role]}</text>')
foot_y = ly + 20 + ((len(items)+1)//2) * 20 + 8
svg.append(f'<text x="{left}" y="{foot_y}" font-size="10" fill="#777">★ polar/H-bond contact; circle area encodes number of contacting atom pairs (pale→deep).</text>')
svg.append('</svg>')

svg_path = os.path.join(FIG, "contact_fingerprint_LPPGP_APPSQ.svg")
with open(svg_path, "w") as f:
    f.write("\n".join(svg))
print(f"[OK] {svg_path}")

# ---------------- per-position contact profile (bar) ----------------
seqs = {d["peptide"]: d["seq"] for d in DET}
ppos = {d["peptide"]: {p[0]: p[2] for p in d["contact_per_position"]} for d in DET}
targets = ["LPPGP", "APPSQ"]
barw = 34; group_gap = 26; left2 = 60; top2 = 60
W2 = left2 + len(targets)*5*barw + len(targets)*group_gap + 40
H2 = 320
svg2 = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W2}" height="{H2}" '
        f'font-family="DejaVu Sans, Arial, sans-serif" font-size="11">']
svg2.append(f'<rect width="{W2}" height="{H2}" fill="#fff"/>')
svg2.append(f'<text x="{left2}" y="26" font-size="15" font-weight="bold">Per-ligand-position contact count (minimized complex)</text>')
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
        # position letter (from seq)
        svg2.append(f'<text x="{bx+barw/2:.1f}" y="{y_bot+32}" text-anchor="middle" font-size="10" font-weight="bold">{seq[pi0]}{pos}</text>')
svg2.append('</svg>')
svg2_path = os.path.join(FIG, "contact_per_position.svg")
with open(svg2_path, "w") as f:
    f.write("\n".join(svg2))
print(f"[OK] {svg2_path}")

# ---------------- HTML combined view ----------------
html = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<title>LPPGP / APPSQ contact fingerprint</title>
<style>body{{font-family:system-ui,Arial,sans-serif;margin:24px;color:#222}}
h1{{font-size:20px}} .fig{{border:1px solid #ddd;border-radius:8px;padding:12px;margin:16px 0}}
.cap{{color:#555;font-size:13px;margin-top:8px}}</style></head>
<body><h1>LPPGP / APPSQ contact fingerprint visualization</h1>
<div class="fig"><img src="contact_fingerprint_LPPGP_APPSQ.svg" width="100%">
<div class="cap">Fig.1. Residue-level contact fingerprint (minimized-complex geometry). Each column is a receptor
residue that contacts either peptide, sorted by residue number; the upper/lower lanes are LPPGP and APPSQ
respectively. Circle area ∝ number of contacting atom pairs (same color shading); ★ marks polar/H-bond contacts;
border color = the residue's role in the DPP4 pocket (red = catalytic core Glu146/Arg147, orange = S1/S2 pocket,
gray = glycan bias, blue = interface / pocket wall).</div></div>
<div class="fig"><img src="contact_per_position.svg" width="70%">
<div class="cap">Fig.2. Per-ligand-position contact count. Horizontal axis = the peptide's 5 residue positions
(letter + index); bar height = number of pocket-contacting atom pairs at that position.
LPPGP contacts concentrate at Gly4 + Pro2; APPSQ concentrates at Ser4 + Pro3, with a broader footprint.</div></div>
<p style="color:#777;font-size:12px">Contact criteria: total contact ≤4.5 Å; hydrophobic (both carbons) ≤4.0 Å;
polar/H-bond (both O/N) ≤3.5 Å. Geometry from the MMFF94s-relaxed complex (consistent with §3.6/§3.7).</p>
</body></html>"""
html_path = os.path.join(FIG, "contact_fingerprint.html")
with open(html_path, "w") as f:
    f.write(html)
print(f"[OK] {html_path}")
