# -*- coding: utf-8 -*-
"""选题B 投稿前图片升级：
1) Fig_workflow.png 研究流程图
2) Fig_heatmap.png 三队列 top DEG 热图（GSE83453）
3) 多面板合成：Fig2_expression = 火山图(A)+热图(B)；Fig5_model = ROC(A)+校准(B)+DCA(C)+列线图(D)
"""
import os, sys, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from PIL import Image, ImageDraw, ImageFont

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pipe_common import load_gse83453, de, bh_fdr

plt.rcParams["font.family"] = "Arial"
OUT = r"C:\Users\98757\Desktop\生信分析\论文_选题B"

# ============ 1) 流程图（自上而下流向） ============
fig, ax = plt.subplots(figsize=(9.6, 8.0), dpi=300)
ax.set_xlim(0, 100); ax.set_ylim(-14, 100); ax.axis("off")

def box(x, y, w, h, text, fc="#EAF1F3", ec="#2E6E6A", fs=8.2, bold=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.6",
                                fc=fc, ec=ec, lw=1.3))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
            fontweight="bold" if bold else "normal", color="#1F3238")

def arrow(x1, y1, x2, y2, rad=0.0):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                 mutation_scale=13, color="#78909C", lw=1.3,
                                 connectionstyle=f"arc3,rad={rad}"))

ax.text(50, 98.5, "Integrated analysis framework", ha="center", fontsize=11.5, fontweight="bold", color="#1F3238")

# Row A 数据层
box(3, 84, 36, 11, "GEO bulk transcriptomes\nGSE83453 (19 vs 8) · GSE153555 (5 AS vs 5 NC)\nGSE12644 (10 vs 10, untouched)", "#E3F0E8")
box(43, 84, 26, 11, "Single-cell RNA-seq\nGSE220774\n82,356 cells, 5 donors", "#E3F0E8")
box(73, 84, 25, 11, "Genetics\nFinnGen R12 (12,418 cases)\neQTL Catalogue / GTEx v8", "#E3F0E8")

# Row B/C 左列主流程
box(3, 63, 36, 13, "Discovery & replication\n561 DEGs ∩ 1,174 immune ∩ 2,132 secreted\n= 37 candidates → 28 replicated\n(GSE153555 + DESeq2)", "#FFFFFF")
box(3, 42, 36, 13, "Machine-learning consensus\nLASSO ∩ SVM-RFE ∩ Random forest\nHub genes: CD74, MMP9", "#FFFFFF")
arrow(21, 84, 21, 77.5)
arrow(21, 63, 21, 56.5)

# Row D/E 下部分支
box(3, 21, 36, 13, "Two-gene diagnostic model\nNested CV AUC 0.868 (0.691–1.000)\nReplication 0.920 · Untouched 0.700\nCalibration · DCA · Nomogram", "#FDF3E7", ec="#C08030")
box(43, 42, 26, 13, "Single-cell localization\nCD74: macrophages + diseased\nendothelium · MMP9: Mo5", "#EFE7F5", ec="#6A4A8A")
box(72, 42, 26, 13, "Exploratory MR\n3 cis-eQTL instruments\nNo evidence of causal effect", "#EFE7F5", ec="#6A4A8A")
box(43, 21, 55, 13, "Immune landscape\nMarker-gene z-scores + MCP-counter deconvolution\nHub gene–infiltration correlations", "#FDF3E7", ec="#C08030")
arrow(21, 42, 21, 35.5)            # ML → model
arrow(39.5, 48.5, 42.5, 48.5)      # ML → scRNA characterization
arrow(56, 84, 56, 56.5)            # scRNA data → scRNA analysis
arrow(85.5, 84, 85.5, 56.5)        # genetics → MR
arrow(56, 42, 60, 35.5)            # scRNA → immune 汇入（单细胞解析免疫景观）
arrow(39.5, 45, 42.5, 34.5)        # ML/hub genes → immune landscape
arrow(90, 42, 80, -1.2, rad=-0.38)  # MR → 结论（绕开 immune 框右侧）

# 探索性转化层
box(43, 3, 55, 12, "Exploratory translation (hypothesis-generating)\nLINCS L1000 signature reversers (dapagliflozin, celecoxib)\nMMP9 docking · ceRNA network (14 shared miRNAs)", "#EFE7F5", ec="#6A4A8A")
arrow(39.5, 45, 50, 16.5, rad=-0.25)  # hub genes → 探索性转化
arrow(70, 21, 70, 16.5)               # immune → 探索性转化
arrow(70, 3, 62, -1.0)                # 探索性转化 → 结论

# 结论
box(18, -12.5, 64, 10, "CD74 + MMP9: candidate tissue transcript-level immune-related\nfeatures of calcific aortic valve disease", "#D9E8E7", ec="#2E6E6A", fs=9.5, bold=True)
arrow(21, 21, 34, -1.0)            # model → conclusion
plt.tight_layout()
plt.savefig(os.path.join(OUT, "Fig_workflow.png"), bbox_inches="tight")
plt.close(fig)
print("saved Fig_workflow.png")

# ============ 2) 热图（GSE83453 top DEGs） ============
G1, gr1 = load_gse83453()
res = de(G1, gr1, "Calcified", "Normal")
sig = res[res["FDR"] < 0.05]
top = pd.concat([sig.nlargest(25, "log2FC"), sig.nsmallest(25, "log2FC")]).index.tolist()
top = [g for g in top if g in G1.index]
M = G1.loc[top]
Z = M.sub(M.mean(axis=1), axis=0).div(M.std(axis=1).replace(0, np.nan), axis=0).fillna(0)
order = np.argsort((gr1 == "Calcined").values.astype(int)) if False else np.argsort((gr1 == "Calcified").values.astype(int))
Z = Z.iloc[:, order]
