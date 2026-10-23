# -*- coding: utf-8 -*-
"""重画手稿B Fig6：GSE153555 面板标注改为 DESeq2 统计量"""
import sys, os
sys.path.insert(0, r'C:\Users\98757\Desktop\生信分析\data')
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pipe_common import load_gse83453, load_gse153555, DATA

OUT = r'C:\Users\98757\Desktop\生信分析\论文_选题B'
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.unicode_minus'] = False

G1, gr1 = load_gse83453()
G2, gr2 = load_gse153555()
d2 = pd.read_csv(os.path.join(DATA, 'r_crosscheck', 'R_DESeq2_GSE153555.csv')).set_index('gene')

hub = ['CD74', 'MMP9']
fig, axes = plt.subplots(1, 4, figsize=(12.8, 4.2), dpi=300)
for j, g in enumerate(hub):
    # GSE83453 面板保持原样（Mann-Whitney）
    ax = axes[j * 2]
    d = [G1.loc[g, gr1 == 'Normal'].values, G1.loc[g, gr1 == 'Calcified'].values]
    bp = ax.boxplot(d, tick_labels=['Normal', 'Calcified'], patch_artist=True, widths=0.55, medianprops=dict(color='black'))
    bp['boxes'][0].set_facecolor('#7B9E89'); bp['boxes'][1].set_facecolor('#C0504D')
    p = stats.mannwhitneyu(d[0], d[1]).pvalue
    ax.set_title(f'{g} (GSE83453)\nP={p:.1e}', fontsize=9)
    ax.set_ylabel('Expression')
    # GSE153555 面板：DESeq2 统计量
    ax = axes[j * 2 + 1]
    d = [G2.loc[g, gr2 == 'NC'].values, G2.loc[g, gr2 == 'AS'].values]
    bp = ax.boxplot(d, tick_labels=['NC', 'AS'], patch_artist=True, widths=0.55, medianprops=dict(color='black'))
    bp['boxes'][0].set_facecolor('#7B9E89'); bp['boxes'][1].set_facecolor('#C0504D')
    row = d2.loc[g]
    ax.set_title(f'{g} (GSE153555)\nDESeq2 FDR={row["FDR"]:.1e}', fontsize=9)
    ax.set_ylabel('Expression')
fig.suptitle('Hub gene expression in discovery and validation cohorts', fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(os.path.join(OUT, 'Fig6_hub_expression.png'))
plt.close(fig)
for g in hub:
    r = d2.loc[g]
    print(g, 'DESeq2 log2FC=%.3f P=%.2e FDR=%.2e' % (r['log2FC'], r['P.Value'], r['FDR']))
print('Fig6 已更新')
