# -*- coding: utf-8 -*-
"""
[受体准备] 1WCY.pdb -> 1WCY_receptor.pdbqt
步骤: 去水/去原配体(sitagliptin A1201) -> ADT prepare_receptor4.py 转 pdbqt
依赖: AutoDockTools (pythonsh + MGLTools) 或 openbabel(obabel -xr)
本机无 ADT, 此处生成干净受体 pdb 供用户本机转换。
"""
import re
SRC="E:/workbuddy/Claw/1WCY.pdb"
OUT="E:/workbuddy/Claw/1WCY_clean.pdb"
keep=[]
for line in open(SRC):
    rec=line[:6]
    if rec in ("TER   ","END   ","ENDMDL"):
        keep.append(line); continue
    if rec=="HETATM":
        # 去除配体 A1201(sitagliptin) 与结晶水(HOH); 保留离子/甘油等可酌情
        name=line[17:20].strip()
        if name=="A1201" or name=="HOH":
            continue
    if rec in ("ATOM  ","HETATM"):
        keep.append(line); continue
    keep.append(line)
open(OUT,"w").writelines(keep)
print(f"已写干净受体(去配体/水): {OUT}  ({len(keep)} 行)")

print("\n--- 用户本机执行(装好 ADT/MGLTools 后) ---")
print("  pythonsh $MGLTOOLS/prepare_receptor4.py -r 1WCY_clean.pdb -o 1WCY_receptor.pdbqt -A checkhydrogens")
print("  或: obabel 1WCY_clean.pdb -O 1WCY_receptor.pdbqt -xr")
print(f"\n口袋(grid)中心(实测 sitagliptin A1201): center_x=62.8 center_y=47.7 center_z=4.8  size=30 30 30")
print("模板论文给 (54,62,37) 为同口袋近似, 以实测为准更稳。")
