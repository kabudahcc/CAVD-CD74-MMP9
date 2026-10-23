# -*- coding: utf-8 -*-
"""选题B 升级4：GSE12644 第三独立队列验证（GPL570，20 例人主动脉瓣，Case=钙化主动脉瓣）
输出：DE 结果、hub 基因统计、双基因模型在本队列的 AUC、Fig6b 表达图
"""
import os, sys, json
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pipe_common import read_series_matrix, collapse_probes, de, bh_fdr, DATA, load_gse83453

plt.rcParams["font.family"] = "Arial"
OUT = r"C:\Users\98757\Desktop\生信分析\论文_选题B"

expr, titles = read_series_matrix(os.path.join(DATA, "GSE12644_series_matrix.txt.gz"))
gr = pd.Series(["Case" if "Case" in t else "Ctrl" for t in titles], index=expr.columns)
print("samples:", gr.value_counts().to_dict())

p2s = json.load(open(os.path.join(DATA, "gpl570_probe2symbol.json")))
p2s = {k: v for k, v in p2s.items() if v}
G = collapse_probes(expr, p2s)
print("genes:", G.shape, "range:", float(G.values.min()), float(G.values.max()))
# 若非 log 尺度则 log2 转换
if G.values.max() > 100:
    G = np.log2(G + 1)
    print("log2 transformed")

res = de(G, gr, "Case", "Ctrl")
res.to_csv(os.path.join(OUT, "TableS_DE_GSE12644.csv"))

hub_stats = []
for g in ["CD74", "MMP9"]:
    if g in res.index:
        r = res.loc[g]
        hub_stats.append(dict(gene=g, log2FC=r["log2FC"], P=r["P.Value"], FDR=r["FDR"]))
        print(g, dict(hub_stats[-1]))
    else:
        print(g, "NOT in matrix")

# 双基因模型 AUC（队列内 z-score + 训练系数）
G1, gr1 = load_gse83453()
from sklearn.linear_model import LogisticRegression
trX = G1.T[["CD74", "MMP9"]].values
ytr = (gr1 == "Calcified").astype(int).values
sc = StandardScaler().fit(trX)
clf = LogisticRegression(penalty="l2", C=1.0, max_iter=5000).fit(sc.transform(trX), ytr)
teX = G.loc[["CD74", "MMP9"]].T.values
sc3 = StandardScaler().fit(teX)
score = (sc3.transform(teX) * clf.coef_[0]).sum(axis=1)
yte = (gr == "Case").astype(int).values
auc3 = roc_auc_score(yte, score)
print("GSE12644 2-gene AUC:", auc3)

# Fig6b：hub 基因在第三队列的表达
fig, axes = plt.subplots(1, 2, figsize=(6.4, 3.4), dpi=300)
for ax, g in zip(axes, ["CD74", "MMP9"]):
    a = G.loc[g, gr == "Case"]; b = G.loc[g, gr == "Ctrl"]
    bp = ax.boxplot([b, a], tick_labels=["Control", "CAVD"], widths=0.55, patch_artist=True,
                    medianprops=dict(color="black"))
    for patch, c in zip(bp["boxes"], ["#8EC3E6", "#E8A198"]):
        patch.set_facecolor(c)
    for i, v in enumerate([b, a]):
        ax.scatter(np.full(len(v), i + 1) + np.random.uniform(-0.09, 0.09, len(v)),
                   v, s=12, c="black", alpha=0.6, zorder=3)
    r = res.loc[g]
    ax.set_title(f"{g}\nlog2FC={r['log2FC']:.2f}, FDR={r['FDR']:.3g}", fontsize=9)
    ax.set_ylabel("Expression (log2)")
    ax.tick_params(labelsize=8)
fig.suptitle("Third-cohort validation: GSE12644 (human aortic valves, GPL570)", fontsize=10)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(os.path.join(OUT, "Fig6b_GSE12644.png"), bbox_inches="tight")
plt.close(fig)

json.dump(dict(hub_stats=hub_stats, auc_2gene=auc3, n_case=int((gr == "Case").sum()),
               n_ctrl=int((gr == "Ctrl").sum()), n_genes=int(G.shape[0])),
          open(os.path.join(OUT, "results_upgrade4.json"), "w"), indent=1)
print("saved Fig6b + TableS_DE_GSE12644")
