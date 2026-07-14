#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Render phaseC_network.json -> a static SVG network diagram."""
import json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "..", "data", "phaseC")
NET = os.path.join(DATA, "phaseC_network.json")
OUT = os.path.join(DATA, "phaseC_network.svg")

net = json.load(open(NET, encoding="utf-8"))
nodes = {n["id"]: n for n in net["nodes"]}
edges = net["edges"]

CX, CY = 340, 340

# separate DPP4, STRING neighbours, literature substrates
dpp4 = net["central_target"]
string_n = [n["id"] for n in net["nodes"]
            if n["source"] == "string" and n["id"] != dpp4]
lit_n = [n["id"] for n in net["nodes"] if n["source"] != "string"]

# layout: STRING neighbours on inner ring r=190, literature on outer ring r=300
def ring(ids, r, start=0.0):
    out = {}
    n = len(ids)
    for i, k in enumerate(ids):
        a = start + 2 * math.pi * i / max(n, 1)
        out[k] = (CX + r * math.cos(a), CY + r * math.sin(a))
    return out

pos = {dpp4: (CX, CY)}
pos.update(ring(string_n, 190, -math.pi / 2))
pos.update(ring(lit_n, 300, 0.0))

# node colors by source
def fill(nid):
    n = nodes[nid]
    if nid == dpp4:
        return "#e6194b"          # central target red
    if n["source"] == "string":
        return "#4361ee"          # STRING neighbour blue
    return "#f4a261"              # literature substrate orange

svg = []
svg.append(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 680 680">')
svg.append('<rect width="680" height="680" fill="#0f1420"/>')
svg.append(f'<text x="340" y="28" fill="#e8eef7" font-size="17" '
           f'text-anchor="middle" font-family="Segoe UI,Arial">'
           f'DPP4-centred network pharmacology (15 nodes / {len(edges)} edges)</text>')

# edges
for e in edges:
    a, b = e["source"], e["target"]
    if a not in pos or b not in pos:
        continue
    x1, y1 = pos[a]; x2, y2 = pos[b]
    if e["type"] == "functional_string":
        col, dash = "#4361ee", "0"
        w = max(0.6, float(e.get("score") or 0) * 2.2)
    else:
        col, dash = "#f4a261", "5 4"
        w = 1.6
    svg.append(f'<line x1="{x1:.0f}" y1="{y1:.0f}" x2="{x2:.0f}" y2="{y2:.0f}" '
               f'stroke="{col}" stroke-width="{w:.2f}" stroke-dasharray="{dash}" opacity="0.75"/>')

# nodes
for nid, (x, y) in pos.items():
    n = nodes[nid]
    c = fill(nid)
    r = 16 if nid == dpp4 else (13 if n["source"] == "string" else 11)
    svg.append(f'<circle cx="{x:.0f}" cy="{y:.0f}" r="{r}" fill="{c}" '
               f'stroke="#0f1420" stroke-width="2"/>')
    # label
    svg.append(f'<text x="{x:.0f}" y="{y+r+13:.0f}" fill="#cdd6e4" font-size="11" '
               f'text-anchor="middle" font-family="Segoe UI,Arial">{nid}</text>')

# legend
lx, ly = 20, 620
svg.append(f'<circle cx="{lx}" cy="{ly}" r="7" fill="#e6194b"/>'
           f'<text x="{lx+14}" y="{ly+4}" fill="#cdd6e4" font-size="11" '
           f'font-family="Segoe UI,Arial">DPP4 target</text>')
svg.append(f'<circle cx="{lx}" cy="{ly+20}" r="7" fill="#4361ee"/>'
           f'<text x="{lx+14}" y="{ly+24}" fill="#cdd6e4" font-size="11" '
           f'font-family="Segoe UI,Arial">STRING functional neighbour</text>')
svg.append(f'<circle cx="{lx}" cy="{ly+40}" r="7" fill="#f4a261"/>'
           f'<text x="{lx+14}" y="{ly+44}" fill="#cdd6e4" font-size="11" '
           f'font-family="Segoe UI,Arial">Literature DPP4 substrate</text>')
svg.append(f'<line x1="{lx+150}" y1="{ly}" x2="{lx+185}" y2="{ly}" '
           f'stroke="#4361ee" stroke-width="2"/>'
           f'<text x="{lx+192}" y="{ly+4}" fill="#cdd6e4" font-size="11" '
           f'font-family="Segoe UI,Arial">functional (STRING)</text>')
svg.append(f'<line x1="{lx+150}" y1="{ly+20}" x2="{lx+185}" y2="{ly+20}" '
           f'stroke="#f4a261" stroke-width="2" stroke-dasharray="5 4"/>'
           f'<text x="{lx+192}" y="{ly+24}" fill="#cdd6e4" font-size="11" '
           f'font-family="Segoe UI,Arial">cleavage substrate (lit.)</text>')

svg.append('</svg>')
open(OUT, "w", encoding="utf-8").write("\n".join(svg))
print(f"[INFO] wrote {OUT}")
