# -*- coding: utf-8 -*-
"""
Render the Tier0 length-debiasing result as an SVG bar chart matching the
project figure style (white background, Arial). Zero external dependencies.

Reads tier0_results.tsv and emits figures/phaseA_tier0_debias.svg (AUC, with SD
error bars) comparing four settings + permutation null.
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "tier0_results.tsv")
FIG = os.path.abspath(os.path.join(HERE, "..", "..", "figures", "phaseA_tier0_debias.svg"))

# short labels + colors (red = positive/SCM signal per project convention)
BARS = [
    ("Confounded\nlength baseline", "0.8845", "0.0000", "#7f8c8d"),
    ("Confounded\nSCM",             "0.6906", "0.0065", "#e67e22"),
    ("Length-matched\nlength baseline", "0.5000", "0.0000", "#95a5a6"),
    ("Length-matched\nSCM",         "0.7106", "0.0174", "#c0392b"),
    ("Length-matched\nSCM (label-perm null)", "0.4901", "0.0336", "#bdc3c7"),
]


def load_from_tsv():
    """Prefer live values from the TSV when present."""
    if not os.path.exists(RES):
        return
    m = {}
    with open(RES, encoding="utf-8") as f:
        f.readline()
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 8:
                m[p[0]] = (p[6], p[7])  # AUC_mean, AUC_sd
    mapping = {
        "Confounded / length-baseline": 0,
        "Confounded / SCM": 1,
        "Length-matched / length-baseline": 2,
        "Length-matched / SCM": 3,
        "Length-matched / SCM permutation-null": 4,
    }
    for k, idx in mapping.items():
        if k in m:
            lab, _, _, col = BARS[idx]
            BARS[idx] = (lab, m[k][0], m[k][1], col)


def main():
    load_from_tsv()
    W, H = 760, 460
    x0, y0 = 90, 360         # axis origin
    plot_w, plot_h = 620, 280
    ymin, ymax = 0.4, 0.95   # AUC range for visual contrast
    n = len(BARS)
    slot = plot_w / n
    bw = slot * 0.5

    def yv(v):
        return y0 - (v - ymin) / (ymax - ymin) * plot_h

    s = []
    s.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
             f'viewBox="0 0 {W} {H}" font-family="Arial, Helvetica, sans-serif">')
    s.append(f'<rect width="{W}" height="{H}" fill="#ffffff"/>')
    s.append(f'<text x="{W/2}" y="30" text-anchor="middle" font-size="18" '
             f'font-weight="bold" fill="#1a1a1a">Tier0 length-debiased benchmark validation (ROC-AUC)</text>')
    s.append(f'<text x="{W/2}" y="52" text-anchor="middle" font-size="12" fill="#555">'
             f'50x repeated length-matched sampling + stratified 5-fold CV (mean +/- SD)</text>')

    # gridlines + y labels
    for gv in [0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
        yy = yv(gv)
        s.append(f'<line x1="{x0}" y1="{yy:.1f}" x2="{x0+plot_w}" y2="{yy:.1f}" '
                 f'stroke="#e0e0e0" stroke-width="1"/>')
        s.append(f'<text x="{x0-8}" y="{yy+4:.1f}" text-anchor="end" font-size="11" '
                 f'fill="#555">{gv:.2f}</text>')
    # 0.5 chance line emphasized
    yy5 = yv(0.5)
    s.append(f'<line x1="{x0}" y1="{yy5:.1f}" x2="{x0+plot_w}" y2="{yy5:.1f}" '
             f'stroke="#333" stroke-width="1.2" stroke-dasharray="5,4"/>')
    s.append(f'<text x="{x0+plot_w-4}" y="{yy5-6:.1f}" text-anchor="end" font-size="10" '
             f'fill="#333">AUC = 0.5 (chance)</text>')

    # axes
    s.append(f'<line x1="{x0}" y1="{y0}" x2="{x0+plot_w}" y2="{y0}" stroke="#333" stroke-width="1.5"/>')
    s.append(f'<line x1="{x0}" y1="{y0-plot_h}" x2="{x0}" y2="{y0}" stroke="#333" stroke-width="1.5"/>')
    s.append(f'<text x="34" y="{y0-plot_h/2}" text-anchor="middle" font-size="12" '
             f'fill="#333" transform="rotate(-90 34 {y0-plot_h/2})">ROC-AUC</text>')

    for i, (lab, mean_s, sd_s, col) in enumerate(BARS):
        v = float(mean_s); sd = float(sd_s)
        cx = x0 + slot * (i + 0.5)
        bx = cx - bw/2
        by = yv(v)
        s.append(f'<rect x="{bx:.1f}" y="{by:.1f}" width="{bw:.1f}" height="{y0-by:.1f}" '
                 f'fill="{col}"/>')
        # value label
        s.append(f'<text x="{cx:.1f}" y="{by-8:.1f}" text-anchor="middle" font-size="12" '
                 f'font-weight="bold" fill="#1a1a1a">{v:.3f}</text>')
        # error bar
        if sd > 0:
            e_hi = yv(v + sd); e_lo = yv(v - sd)
            s.append(f'<line x1="{cx:.1f}" y1="{e_hi:.1f}" x2="{cx:.1f}" y2="{e_lo:.1f}" '
                     f'stroke="#333" stroke-width="1.2"/>')
            s.append(f'<line x1="{cx-5:.1f}" y1="{e_hi:.1f}" x2="{cx+5:.1f}" y2="{e_hi:.1f}" stroke="#333" stroke-width="1.2"/>')
            s.append(f'<line x1="{cx-5:.1f}" y1="{e_lo:.1f}" x2="{cx+5:.1f}" y2="{e_lo:.1f}" stroke="#333" stroke-width="1.2"/>')
        # x label (two lines)
        for j, part in enumerate(lab.split("\n")):
            s.append(f'<text x="{cx:.1f}" y="{y0+18+j*13:.1f}" text-anchor="middle" '
                     f'font-size="10.5" fill="#333">{part}</text>')

    s.append('</svg>')
    with open(FIG, "w", encoding="utf-8") as f:
        f.write("\n".join(s))
    print(f"Written figure: {FIG}")


if __name__ == "__main__":
    main()
