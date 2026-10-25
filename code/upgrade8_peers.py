# -*- coding: utf-8 -*-
"""选题B 升级8（学习同类文章）：
1) 单基因 ROC：CD74 / MMP9 各自在三个队列的 AUC + bootstrap 95% CI
2) STRING PPI 网络：28 个复制候选基因的互作网络（score>=0.4）
输出：TableS_pergene_ROC.csv、results_upgrade8.json、Fig4 面板 C 素材（ppi_network.png）
"""
import os, sys, json
import numpy as np
import pandas as pd
import requests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
from sklearn.metrics import roc_auc_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pipe_common import load_gse83453, load_gse153555, read_series_matrix, collapse_probes, DATA

OUT = r"C:\Users\98757\Desktop\生信分析\论文_选题B"
HUB = ["CD74", "MMP9"]
rng = np.random.RandomState(42)

G1, gr1 = load_gse83453()
G2, gr2 = load_gse153555()
expr, titles = read_series_matrix(os.path.join(DATA, "GSE12644_series_matrix.txt.gz"))
gr3 = pd.Series(["Case" if "Case" in t else "Ctrl" for t in titles], index=expr.columns)
p2s = {k: v for k, v in json.load(open(os.path.join(DATA, "gpl570_probe2symbol.json"))).items() if v}
G3 = collapse_probes(expr, p2s)
if G3.values.max() > 100:
    G3 = np.log2(G3 + 1)

cohorts = [("GSE83453", G1, (gr1 == "Calcified").astype(int).values),
           ("GSE153555", G2, (gr2 == "AS").astype(int).values),
           ("GSE12644", G3, (gr3 == "Case").astype(int).values)]

def boot_auc(y, s, B=2000):
    aucs = []
    n1, n0 = (y == 1).sum(), (y == 0).sum()
    i1, i0 = np.where(y == 1)[0], np.where(y == 0)[0]
    for _ in range(B):
        idx = np.concatenate([rng.choice(i1, n1, replace=True), rng.choice(i0, n0, replace=True)])
        aucs.append(roc_auc_score(y[idx], s[idx]))
    return np.percentile(aucs, [2.5, 97.5])

rows = []
for name, G, y in cohorts:
    for g in HUB:
        s = G.loc[g].values
        auc = roc_auc_score(y, s)
        lo, hi = boot_auc(y, s)
        rows.append(dict(cohort=name, gene=g, AUC=round(float(auc), 3),
                         CI_low=round(float(lo), 3), CI_high=round(float(hi), 3)))
        print(f"{name} {g}: AUC={auc:.3f} ({lo:.3f}-{hi:.3f})")
roc_df = pd.DataFrame(rows)
roc_df.to_csv(os.path.join(OUT, "TableS_pergene_ROC.csv"), index=False, encoding="utf-8-sig")

# ---------- STRING PPI（28 个候选） ----------
cand28 = json.load(open(os.path.join(OUT, "results_B.json"), encoding="utf-8"))["cand_genes"]
url = "https://string-db.org/api/tsv/network"
r = requests.post(url, data={"identifiers": "%0d".join(cand28), "species": 9606,
                             "required_score": 400}, timeout=60)
r.raise_for_status()
from io import StringIO
edges = pd.read_csv(StringIO(r.text), sep="\t")
edges = edges[["preferredName_A", "preferredName_B", "score"]]
edges.columns = ["gene_A", "gene_B", "combined_score"]
edges.to_csv(os.path.join(OUT, "TableS_PPI_STRING.csv"), index=False, encoding="utf-8-sig")

G_net = nx.Graph()
for _, e in edges.iterrows():
    G_net.add_edge(e["gene_A"], e["gene_B"], weight=e["combined_score"])
deg = dict(G_net.degree())
deg_s = pd.Series(deg).sort_values(ascending=False)
print("PPI 节点数:", G_net.number_of_nodes(), "边数:", G_net.number_of_edges())
print("degree 前五:", deg_s.head(5).to_dict())
print("CD74 degree:", deg.get("CD74"), "| MMP9 degree:", deg.get("MMP9"))

fig, ax = plt.subplots(figsize=(7.2, 6.4), dpi=300)
pos = nx.spring_layout(G_net, k=0.55, seed=42)
node_sizes = [40 + 60 * deg[n] for n in G_net.nodes]
colors = ["#C0504D" if n in HUB else "#78909C" for n in G_net.nodes]
nx.draw_networkx_edges(G_net, pos, ax=ax, alpha=0.35, edge_color="#90A4AE", width=1.0)
nx.draw_networkx_nodes(G_net, pos, ax=ax, node_size=node_sizes, node_color=colors,
                       edgecolors="white", linewidths=0.6)
for n in G_net.nodes:
    fs = 9 if n in HUB else 6.5
    fw = "bold" if n in HUB else "normal"
    ax.annotate(n, pos[n], textcoords="offset points", xytext=(4, 3), fontsize=fs,
                fontweight=fw, color="#1F3238")
ax.set_title("STRING PPI network of the 28 replicated candidates\n(node size ∝ degree; hub genes in red)", fontsize=10)
ax.axis("off")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "Fig4_PPI_network.png"), bbox_inches="tight")
plt.close(fig)

json.dump({"pergene_roc": rows,
           "ppi": {"n_nodes": G_net.number_of_nodes(), "n_edges": G_net.number_of_edges(),
                   "degree_top5": {k: int(v) for k, v in deg_s.head(5).items()},
                   "CD74_degree": int(deg.get("CD74", 0)), "MMP9_degree": int(deg.get("MMP9", 0)),
                   "string_min_score": 0.4}},
          open(os.path.join(OUT, "results_upgrade8.json"), "w", encoding="utf-8"), indent=1)
print("saved TableS_pergene_ROC.csv + TableS_PPI_STRING.csv + Fig4_PPI_network.png + results_upgrade8.json")
