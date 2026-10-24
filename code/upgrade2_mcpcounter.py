# -*- coding: utf-8 -*-
"""选题B 升级2：MCPcounter 标准反卷积（官方 markers，ebecht/MCPcounter）
评分 = 各细胞群标志基因表达均值（log 尺度），在 GSE83453 与 GSE153555 上计算
"""
import os, sys, json
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pipe_common import load_gse83453, load_gse153555, bh_fdr

plt.rcParams["font.family"] = "Arial"
OUT = r"C:\Users\98757\Desktop\生信分析\论文_选题B"

mk = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), "mcpcounter_genes.txt"), sep="\t")
markers = {c: g["HUGO symbols"].tolist() for c, g in mk.groupby("Cell population")}
print("cell types:", {k: len(v) for k, v in markers.items()})


def mcp_scores(G):
    """G: genes x samples (log 尺度). 返回 celltypes x samples"""
    avail = set(G.index)
    sc = {}
    for ct, genes in markers.items():
        gg = [g for g in genes if g in avail]
        if len(gg) >= 1:  # 与官方 MCPcounter.estimate 一致：单 marker 也纳入
            sc[ct] = G.loc[gg].mean(axis=0)
    return pd.DataFrame(sc).T, {ct: len([g for g in genes if g in avail]) for ct, genes in markers.items()}


def compare(scores, gr, case, ctrl):
    a = scores.loc[:, gr == case]; b = scores.loc[:, gr == ctrl]
    p = np.array([stats.mannwhitneyu(a.iloc[i].values, b.iloc[i].values,
                                     alternative="two-sided").pvalue for i in range(a.shape[0])])
    r = pd.DataFrame({"cell": scores.index, "mean_case": a.mean(axis=1),
                      "mean_ctrl": b.mean(axis=1), "P": p})
    r["FDR"] = bh_fdr(r["P"].values)
    return r.sort_values("P")


G1, gr1 = load_gse83453()
print("GSE83453 range:", float(G1.values.min()), float(G1.values.max()))
G2, gr2 = load_gse153555()

s1, nmark1 = mcp_scores(G1)
s2, nmark2 = mcp_scores(G2)
print("markers used:", nmark1)

r1 = compare(s1, gr1, "Calcified", "Normal")
mask2 = gr2.isin(["AS", "NC"])
r2 = compare(s2.loc[:, mask2], gr2[mask2], "AS", "NC")
r1.to_csv(os.path.join(OUT, "TableS_MCPcounter_GSE83453.csv"), index=False)
r2.to_csv(os.path.join(OUT, "TableS_MCPcounter_GSE153555.csv"), index=False)
print("\nGSE83453 MCPcounter:")
print(r1.to_string(index=False))
print("\nGSE153555 MCPcounter:")
print(r2.to_string(index=False))

# hub 基因与 MCPcounter 评分的 Spearman 相关（GSE83453）
corr_rows = []
for gene in ["CD74", "MMP9"]:
    x = G1.loc[gene]
    for ct in s1.index:
        rho, p = stats.spearmanr(x, s1.loc[ct])
        corr_rows.append(dict(gene=gene, cell=ct, rho=rho, P=p))
corr = pd.DataFrame(corr_rows)
corr["FDR"] = bh_fdr(corr["P"].values)
corr.to_csv(os.path.join(OUT, "TableS_MCPcounter_hubcorr.csv"), index=False)

# ---------- Fig5b：MCPcounter 分组比较（双队列两行） ----------
sig = r1[r1["FDR"] < 0.05]["cell"].tolist()
order = [c for c in r1["cell"].tolist() if c in set(r2["cell"])]  # 统一细胞群顺序便于跨队列比较

mask2 = gr2.isin(["AS", "NC"])


def draw_panel(ax, scores, gr, r, ct, case, ctrl):
    a = scores.loc[ct, gr == case]; b = scores.loc[ct, gr == ctrl]
    row = r[r["cell"] == ct].iloc[0]
    bp = ax.boxplot([b, a], tick_labels=[ctrl, case], widths=0.55,
                    patch_artist=True, medianprops=dict(color="black"))
    for patch, c in zip(bp["boxes"], ["#8EC3E6", "#E8A198"]):
        patch.set_facecolor(c)
    for i, v in enumerate([b, a]):
        ax.scatter(np.full(len(v), i + 1) + np.random.uniform(-0.08, 0.08, len(v)),
                   v, s=9, c="black", alpha=0.55, zorder=3)
    star = "***" if row["FDR"] < 0.001 else "**" if row["FDR"] < 0.01 else "*" if row["FDR"] < 0.05 else "ns"
    fdr_txt = "FDR<0.001" if row["FDR"] < 0.001 else f"FDR={row['FDR']:.3f}"
    ax.set_title(f"{ct}\n{fdr_txt} {star}", fontsize=8)
    ax.tick_params(labelsize=7.5)


ncols = 5
nrows = int(np.ceil(len(order) / ncols))
fig, axes = plt.subplots(nrows * 2, ncols, figsize=(12.5, 3.1 * nrows * 2), dpi=300)
for j, ct in enumerate(order):
    draw_panel(axes[j // ncols, j % ncols], s1, gr1, r1, ct, "Calcified", "Normal")
for j, ct in enumerate(order):
    draw_panel(axes[nrows + j // ncols, j % ncols], s2.loc[:, mask2], gr2[mask2], r2, ct, "AS", "NC")
for a in axes.ravel()[len(order):nrows * ncols]:
    a.axis("off")
for a in axes.ravel()[nrows * ncols + len(order):]:
    a.axis("off")
axes[0, 0].annotate("GSE83453 (microarray)", xy=(-0.32, 0.5), xycoords="axes fraction",
                    rotation=90, va="center", ha="center", fontsize=10, fontweight="bold")
axes[nrows, 0].annotate("GSE153555 (RNA-seq)", xy=(-0.32, 0.5), xycoords="axes fraction",
                        rotation=90, va="center", ha="center", fontsize=10, fontweight="bold")
fig.suptitle("MCPcounter deconvolution scores in two independent CAVD cohorts", fontsize=11)
fig.tight_layout(rect=[0.02, 0, 1, 0.96])
fig.savefig(os.path.join(OUT, "Fig5_MCPcounter.png"), bbox_inches="tight")
plt.close(fig)

json.dump({"markers_used_GSE83453": nmark1,
           "sig_cells_GSE83453": sig,
           "r1": r1.round(5).to_dict("records"),
           "r2": r2.round(5).to_dict("records")},
          open(os.path.join(OUT, "results_upgrade2.json"), "w"), indent=1)
print("\nsaved Fig5_MCPcounter.png + tables")
