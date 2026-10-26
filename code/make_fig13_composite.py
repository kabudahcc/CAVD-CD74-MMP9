# -*- coding: utf-8 -*-
"""Fig. 13：A LINCS reversers | B MMP9 对接亲和力 | C ceRNA 网络"""
import json, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx

OUT = r"C:\Users\98757\Desktop\生信分析\论文_选题B"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5,
                     "axes.linewidth": 0.8, "axes.spines.top": False, "axes.spines.right": False})

cmap = pd.read_csv(os.path.join(OUT, "TableS_cMAP_LINCS.csv"))
dock = json.load(open(os.path.join(OUT, "results_upgrade9_docking.json")))
cerna = json.load(open(os.path.join(OUT, "results_upgrade9_cerna.json"), encoding="utf-8"))

fig = plt.figure(figsize=(10.5, 8.2))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.35], hspace=0.32, wspace=0.30,
                      left=0.16, right=0.97, top=0.94, bottom=0.08)

# ---- A: top reversers ----
axA = fig.add_subplot(gs[0, 0])
rev = cmap[cmap["zsum"] < 0].nsmallest(10, "zsum").iloc[::-1]
named = [d for d in rev["drug"] if not d.startswith("BRD-")]
colors = ["#d6604d" if d in named else "#bababa" for d in rev["drug"]]
axA.barh(rev["drug"] + " (" + rev["cell"] + ")", rev["zsum"], color=colors, edgecolor="black", lw=0.4)
axA.set_xlabel("z-sum score (negative = signature reverser)")
axA.set_title("A  LINCS L1000 predicted signature reversers", loc="left", fontweight="bold")
axA.axvline(0, color="black", lw=0.6)

# ---- B: docking ----
axB = fig.add_subplot(gs[0, 1])
aff = pd.Series(dock["affinities"]).sort_values()
cols = ["#2166ac" if k == "SB-3CT" else "#4393c3" for k in aff.index]
axB.bar(aff.index, aff.values, color=cols, edgecolor="black", lw=0.4)
for i, v in enumerate(aff.values):
    axB.text(i, v - 0.25, f"{v:.2f}", ha="center", va="top", fontsize=7.5)
axB.axhline(aff["SB-3CT"], color="#2166ac", ls="--", lw=0.8)
axB.set_ylabel("Best affinity (kcal/mol)")
axB.set_ylim(0, min(aff.values) - 1.5)
axB.set_title("B  Docking to MMP9 catalytic domain (1GKC)", loc="left", fontweight="bold")
axB.tick_params(axis="x", rotation=25)
axB.text(0.02, 0.03, "dashed line: SB-3CT (known MMP9 inhibitor)", transform=axB.transAxes, fontsize=7)

# ---- C: ceRNA network ----
axC = fig.add_subplot(gs[1, :])
G = nx.Graph()
for e in cerna["lncrna_mirna"]:
    G.add_edge(e["lncRNA"], e["miRNA"])
for mi in {e["miRNA"] for e in cerna["lncrna_mirna"]}:
    G.add_edge(mi, "CD74"); G.add_edge(mi, "MMP9")
mrna = ["CD74", "MMP9"]
mir = sorted({e["miRNA"] for e in cerna["lncrna_mirna"]})
lnc = sorted({e["lncRNA"] for e in cerna["lncrna_mirna"]})

pos = {}
pos["CD74"] = (-0.6, 0); pos["MMP9"] = (0.6, 0)
# miRNA 中间椭圆环
for i, m in enumerate(mir):
    a = 2 * np.pi * i / len(mir)
    pos[m] = (2.6 * np.cos(a), 1.5 * np.sin(a))
# lncRNA 外环，靠近其 miRNA 的角度
lnc_ang = {}
for e in cerna["lncrna_mirna"]:
    lnc_ang.setdefault(e["miRNA"], []).append(e["lncRNA"])
for i, m in enumerate(mir):
    a0 = 2 * np.pi * i / len(mir)
    for j, l in enumerate(lnc_ang[m]):
        a = a0 + (j - (len(lnc_ang[m]) - 1) / 2) * 0.09
        pos[l] = (4.6 * np.cos(a), 2.8 * np.sin(a))

nx.draw_networkx_edges(G, pos, ax=axC, edge_color="#cccccc", width=0.6)
nx.draw_networkx_nodes(G, pos, nodelist=lnc, node_color="#92c5de", node_size=28, ax=axC, linewidths=0)
nx.draw_networkx_nodes(G, pos, nodelist=mir, node_color="#f4a582", node_size=260, ax=axC, linewidths=0.4, edgecolors="black")
nx.draw_networkx_nodes(G, pos, nodelist=mrna, node_color="#b2182b", node_size=900, ax=axC, linewidths=0.6, edgecolors="black")
for n in mrna:
    axC.annotate(n, pos[n], ha="center", va="center", fontsize=8, fontweight="bold", color="white")
for m in mir:
    x, y = pos[m]
    r = np.hypot(x / 2.6, y / 1.5)
    axC.annotate(m.replace("hsa-", ""), (x * 1.13, y * 1.16), ha="center", va="center", fontsize=5.6)
axC.set_title("C  ceRNA network: lncRNA (blue) - shared miRNA (orange) - CD74/MMP9 (red)",
              loc="left", fontweight="bold")
axC.axis("off")
axC.set_xlim(-6.4, 6.4); axC.set_ylim(-3.9, 3.9)

fig.savefig(os.path.join(OUT, "Fig13_drug_docking_ceRNA.png"), dpi=300)
print("saved Fig13_drug_docking_ceRNA.png")
