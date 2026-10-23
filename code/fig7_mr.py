# -*- coding: utf-8 -*-
"""Fig7: MR 森林图（CD74/MMP9 对各终点的 OR）"""
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Arial"
MR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "mr")
OUT = r"C:\Users\98757\Desktop\生信分析\论文_选题B"

df = pd.read_csv(os.path.join(MR, "Table_MR_results.csv"))
df["label"] = (df["gene"] + " | " + df["exposure"].str.replace("_", " ")
               + " -> " + df["outcome"].str.replace("_", " ")
               + "\n(IVs=" + df["n_instruments"].astype(str)
               + ", F=" + df["median_F"].round(1).astype(str) + ")")
df = df.iloc[::-1].reset_index(drop=True)

fig, ax = plt.subplots(figsize=(8.2, 4.6), dpi=300)
y = np.arange(len(df))
# 对数坐标，CI 截断显示
lo_clip, hi_clip = 0.05, 20
for i, r in df.iterrows():
    lo, hi = max(r["OR_lo"], lo_clip), min(r["OR_hi"], hi_clip)
    ax.plot([lo, hi], [i, i], color="#4472C4", lw=1.6, zorder=2)
    ax.scatter(r["OR"], i, s=90, color="#C00000", zorder=3, marker="s")
    if r["OR_hi"] > hi_clip:
        ax.plot(hi_clip, i, marker=">", color="#4472C4", markersize=6)
    ax.text(21.5, i, f"OR {r['OR']:.2f} ({r['OR_lo']:.2f}-{r['OR_hi']:.2f}), P={r['p']:.2g}",
            va="center", fontsize=8.5)
ax.axvline(1, color="grey", ls="--", lw=1)
ax.set_xscale("log")
ax.set_xlim(lo_clip, 60)
ax.set_xticks([0.1, 0.25, 0.5, 1, 2, 4, 10])
ax.set_xticklabels(["0.1", "0.25", "0.5", "1", "2", "4", "10"])
ax.set_yticks(y)
ax.set_yticklabels(df["label"], fontsize=8.5)
ax.set_xlabel("Odds ratio for calcific aortic valve stenosis per 1-SD higher genetically predicted expression (log scale)", fontsize=9)
ax.set_title("Exploratory Mendelian randomisation: cis-eQTL instruments for CD74 and MMP9", fontsize=10)
ax.spines[["top", "right"]].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "Fig7_MR_forest.png"), bbox_inches="tight")
print("saved Fig7_MR_forest.png")
