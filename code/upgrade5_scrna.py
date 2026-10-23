# -*- coding: utf-8 -*-
"""选题B 升级5a：GSE220774 单细胞 —— CD74/MMP9 细胞来源定位 + 疾病阶段趋势
Fig11: (A) dotplot 各细胞群表达 (B) 供体疾病阶段趋势
"""
import os, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Arial"
OUT = r"C:\Users\98757\Desktop\生信分析\论文_选题B"

d = json.load(open("data/gse220774_hub_counts.json"))
counts = {g: np.array(v, dtype=float) for g, v in d["genes"].items()}
cells = d["header"]
print({g: (c > 0).mean() for g, c in counts.items()})

meta = pd.read_csv("data/GSE220774_metadata.csv.gz", compression="gzip",
                   usecols=["Unnamed: 0", "sample", "donor", "layer", "nCount_RNA", "new_cluster"])
meta = meta.rename(columns={"Unnamed: 0": "cell"})
meta = meta.set_index("cell").loc[cells]  # 对齐顺序
assert (meta.index == pd.Index(cells)).all()

# 归一化：count / nCount_RNA * 1e4, log1p
norm = {}
for g, c in counts.items():
    norm[g] = np.log1p(c / meta["nCount_RNA"].values * 1e4)

meta["stage"] = meta["donor"].str.extract(r"#(\d)").astype(int)  # 1=非病变 .. 5=晚期钙化
meta["group"] = meta["new_cluster"].str.replace(r"\d+$", "", regex=True)
print(meta["group"].value_counts())

# ---- 每 cluster 统计 ----
order = ["VEC1","VEC2","VEC3","VEC4","VEC5","VEC6","VEC7",
         "VIC1","VIC2","VIC3","VIC4","VIC5","VIC6","VIC7","VIC8","VIC9",
         "Trans1","Trans2","Trans3","Mo1","Mo2","Mo3","Mo4","Mo5","Mo6",
         "T-Cell","Stromal"]
stats_rows = []
for cl in order:
    m = (meta["new_cluster"] == cl).values
    row = {"cluster": cl, "n": int(m.sum())}
    for g in ["CD74", "MMP9"]:
        row[f"{g}_mean"] = float(norm[g][m].mean())
        row[f"{g}_pct"] = float((counts[g][m] > 0).mean() * 100)
    stats_rows.append(row)
sdf = pd.DataFrame(stats_rows)
sdf.to_csv(os.path.join(OUT, "TableS_scRNA_cluster_expression.csv"), index=False)
print(sdf[["cluster","CD74_mean","CD74_pct","MMP9_mean","MMP9_pct"]].round(2).to_string(index=False))

# ---- 供体 × 阶段（每供体所有细胞 & 仅免疫细胞）----
donor_rows = []
for donor, sub in meta.groupby("donor"):
    stage = int(sub["stage"].iloc[0])
    row = {"donor": donor, "stage": stage, "n": len(sub)}
    for g in ["CD74", "MMP9"]:
        idx = sub.index
        pos = meta.index.get_indexer(idx)
        row[f"{g}_all"] = float(norm[g][pos].mean())
        imm = sub["new_cluster"].str.startswith(("Mo", "T-Cell")).values
        if imm.sum() > 10:
            row[f"{g}_immune"] = float(norm[g][pos[imm]].mean())
    donor_rows.append(row)
ddf = pd.DataFrame(donor_rows).sort_values("stage")
print(ddf.round(3).to_string(index=False))

# ---- Fig11 ----
fig = plt.figure(figsize=(9.2, 7.6), dpi=300)
gs = fig.add_gridspec(2, 2, height_ratios=[1.6, 1], width_ratios=[2.6, 1],
                      hspace=0.32, wspace=0.25)
# (A) dotplot
ax = fig.add_subplot(gs[0, :])
xs, ys, sizes, colors = [], [], [], []
genes = ["CD74", "MMP9"]
vmax = max(sdf[[f"{g}_mean" for g in genes]].values.max(), 0.01)
for yi, cl in enumerate(order[::-1]):
    r = sdf[sdf["cluster"] == cl].iloc[0]
    for xi, g in enumerate(genes):
        xs.append(xi); ys.append(yi)
        sizes.append(max(r[f"{g}_pct"], 0.5) * 6)
        colors.append(r[f"{g}_mean"] / vmax)
sc = ax.scatter(xs, ys, s=sizes, c=colors, cmap="Reds", vmin=0, vmax=1, edgecolors="grey", lw=0.4)
ax.set_xticks(range(len(genes))); ax.set_xticklabels(genes, fontsize=10)
ax.set_yticks(range(len(order))); ax.set_yticklabels(order[::-1], fontsize=8)
ax.set_title("A  CD74 and MMP9 expression across 27 aortic valve cell clusters (GSE220774)", fontsize=10, loc="left")
cb = fig.colorbar(sc, ax=ax, fraction=0.025, pad=0.01)
cb.set_label("Mean normalised expression (scaled)", fontsize=8)
cb.ax.tick_params(labelsize=7)
# 大小组例（放在点阵右侧空白区）
for pct, ypos in [(10, 20.5), (30, 16.5), (60, 12)]:
    ax.scatter([2.15], [ypos], s=pct * 6, c="white", edgecolors="grey")
    ax.text(2.45, ypos, f"{pct}%", va="center", fontsize=7.5)
ax.text(2.15, 23.2, "% expressing", fontsize=8, ha="left")
ax.set_xlim(-0.5, 3.1)
ax.set_ylim(-1, len(order))

# (B) 供体阶段趋势
for i, (g, c) in enumerate([("CD74", "#C0504D"), ("MMP9", "#4F81BD")]):
    axb = fig.add_subplot(gs[1, i])
    x = ddf["stage"].values
    y_all = ddf[f"{g}_all"].values
    y_imm = ddf[f"{g}_immune"].values
    axb.plot(x, y_all, "o-", c=c, lw=1.8, label="All cells")
    axb.plot(x, y_imm, "s--", c="#E8862E", lw=1.5, label="Immune (Mo+T)")
    rho_all = np.corrcoef(x, y_all)[0, 1]
    axb.set_xlabel("Disease stage (donor rank)", fontsize=9)
    axb.set_ylabel("Mean normalised expression", fontsize=9)
    axb.set_title(f"{'B' if i==0 else 'C'}  {g} vs fibrocalcification stage (r={rho_all:.2f})", fontsize=9.5, loc="left")
    axb.tick_params(labelsize=8)
    axb.legend(fontsize=7.5)
    axb.set_xticks([1, 2, 3, 4, 5])
    axb.set_xticklabels(["1\nnon-diseased", "2", "3", "4", "5\nadvanced"], fontsize=7.5)
fig.savefig(os.path.join(OUT, "Fig11_scRNA_GSE220774.png"), bbox_inches="tight")
plt.close(fig)

json.dump({"cluster_stats": sdf.round(4).to_dict("records"),
           "donor_stats": ddf.round(4).to_dict("records")},
          open(os.path.join(OUT, "results_upgrade5.json"), "w"), indent=1)
print("saved Fig11 + tables")
