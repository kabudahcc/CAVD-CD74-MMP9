# -*- coding: utf-8 -*-
"""选题B 升级9b：分子对接 —— AutoDock Vina 1.2.7 + MMP9 催化域 (PDB 1GKC, chain A)
中文路径会导致 openbabel/vina 读取失败，全部在 C:\\dock_tmp 纯英文暂存目录中运行，结果拷回工作区。
"""
import os, sys, json, shutil, subprocess
import numpy as np
import pandas as pd
import requests

WS = r"C:\Users\98757\Desktop\生信分析"
DK = os.path.join(WS, "data", "docking")
OUT = os.path.join(WS, "论文_选题B")
TMP = r"C:\dock_tmp"
os.makedirs(TMP, exist_ok=True)
for f in ("vina.exe", "1GKC.pdb"):
    shutil.copy(os.path.join(DK, f), os.path.join(TMP, f))
VINA = os.path.join(TMP, "vina.exe")

# ---------- 1) 受体准备（chain A 蛋白 + ZN/CA） ----------
rec_lines, nf = [], []
for line in open(os.path.join(DK, "1GKC.pdb")):
    rec, chain, resn = line[:6].strip(), line[21], line[17:20].strip()
    if rec == "ATOM" and chain == "A":
        rec_lines.append(line)
    elif rec == "HETATM" and chain == "A" and resn in ("ZN", "CA"):
        rec_lines.append(line)
    elif rec == "HETATM" and resn == "NFH" and chain == "A":
        nf.append(line)
rec_lines.append("END\n")
open(os.path.join(TMP, "receptor.pdb"), "w").writelines(rec_lines)

xyz = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in nf])
center = xyz.mean(axis=0)
print("配体 NFH 原子数:", len(nf), "| 对接盒中心:", np.round(center, 2))

def obabel(infile, in_fmt, out_fmt, outfile, add_h=False):
    code = (f"from openbabel import openbabel as ob\n"
            f"c = ob.OBConversion(); c.SetInAndOutFormats('{in_fmt}','{out_fmt}')\n"
            f"m = ob.OBMol(); assert c.ReadFile(m, r'{infile}')\n"
            + ("m.AddHydrogens()\n" if add_h else "")
            + f"c.WriteFile(m, r'{outfile}'); print(m.NumAtoms())")
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, cwd=TMP)
    if r.returncode != 0:
        raise RuntimeError(f"obabel failed for {infile}: {r.stderr[-300:]}")
    return r.stdout.strip()

print("receptor pdbqt atoms:", obabel(os.path.join(TMP, "receptor.pdb"), "pdb", "pdbqt",
                                      os.path.join(TMP, "receptor.pdbqt")))
# obabel 会写入配体风格 ROOT/BRANCH 标签，刚性受体需剥除
bad = ("ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF")
with open(os.path.join(TMP, "receptor.pdbqt")) as f:
    lines = [l for l in f if not l.startswith(bad)]
open(os.path.join(TMP, "receptor.pdbqt"), "w").writelines(lines)
print("receptor pdbqt 清理后行数:", len(lines))

# ---------- 2) 配体：PubChem 3D SDF → PDBQT ----------
ligands = ["SB-3CT", "celecoxib", "toremifene", "dapagliflozin", "rolipram"]
for name in ligands:
    sdf = os.path.join(TMP, name + ".sdf")
    pdbqt = os.path.join(TMP, name + ".pdbqt")
    if not os.path.exists(pdbqt):
        r = requests.get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/SDF",
                         params={"record_type": "3d"}, timeout=60)
        r.raise_for_status()
        open(sdf, "wb").write(r.content)
        n = obabel(sdf, "sdf", "pdbqt", pdbqt, add_h=True)
        print(name, "pdbqt atoms:", n)

# ---------- 3) 逐个对接 ----------
results = {}
for name in ligands:
    outp = os.path.join(TMP, name + "_out.pdbqt")
    cmd = [VINA, "--receptor", os.path.join(TMP, "receptor.pdbqt"),
           "--ligand", os.path.join(TMP, name + ".pdbqt"),
           "--center_x", f"{center[0]:.3f}", "--center_y", f"{center[1]:.3f}", "--center_z", f"{center[2]:.3f}",
           "--size_x", "22", "--size_y", "22", "--size_z", "22",
           "--exhaustiveness", "16", "--num_modes", "9",
           "--out", outp]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=280, cwd=TMP)
    if r.returncode != 0:
        raise RuntimeError(f"vina failed for {name}: {r.stdout[-300:]} {r.stderr[-300:]}")
    open(os.path.join(TMP, name + "_log.txt"), "w").write(r.stdout)
    aff = None
    for line in r.stdout.splitlines():
        p = line.split()
        if len(p) >= 2 and p[0] == "1":
            try:
                aff = float(p[1]); break
            except ValueError:
                pass
    results[name] = aff
    print(f"{name:15s} best affinity = {aff} kcal/mol")

df = pd.DataFrame([{"compound": k,
                    "role": "positive control (known MMP9 inhibitor)" if k == "SB-3CT" else "LINCS predicted reverser",
                    "best_affinity_kcal_mol": v} for k, v in results.items()])
df = df.sort_values("best_affinity_kcal_mol")
df.to_csv(os.path.join(OUT, "TableS_docking.csv"), index=False, encoding="utf-8-sig")
json.dump({"pdb": "1GKC chain A (MMP9 catalytic domain, residues 107-215/391-443)",
           "box_center": [round(x, 3) for x in center.tolist()], "box_size": [22, 22, 22],
           "exhaustiveness": 16, "affinities": results},
          open(os.path.join(OUT, "results_upgrade9_docking.json"), "w", encoding="utf-8"), indent=1)
# 拷回全部中间文件供复现
for f in os.listdir(TMP):
    if f not in ("vina.exe", "1GKC.pdb"):
        shutil.copy(os.path.join(TMP, f), os.path.join(DK, f))
print("\nsaved TableS_docking.csv + results_upgrade9_docking.json; 中间文件已拷回 data/docking")
