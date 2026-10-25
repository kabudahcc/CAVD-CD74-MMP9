# -*- coding: utf-8 -*-
"""选题B完整流水线：六图 + results_B.json（供手稿引用）"""
import sys, os, json
sys.path.insert(0, r'C:\Users\98757\Desktop\生信分析\data')
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegressionCV, LogisticRegression
from sklearn.svm import SVC
from sklearn.feature_selection import RFE
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_predict, StratifiedKFold
from collections import Counter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pipe_common import bh_fdr, de, load_gse83453, load_gse153555, enrichr, immune_scores, IMMUNE_MARKERS, DATA

OUT = r'C:\Users\98757\Desktop\生信分析\论文_选题B'
os.makedirs(OUT, exist_ok=True)
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.unicode_minus'] = False
R = {}  # numbers for manuscript

lists = json.load(open(os.path.join(DATA, 'gene_lists.json'), encoding='utf-8'))
IMMUNE = set(lists['immune_kegg']); SECRETED = set(lists['secreted_uniprot'])
R['n_immune_list'] = len(IMMUNE); R['n_secreted_list'] = len(SECRETED)

G1, gr1 = load_gse83453()
G2, gr2 = load_gse153555()
R['gse83453_groups'] = gr1.value_counts().to_dict()
R['gse153555_groups'] = gr2.value_counts().to_dict()
R['gse83453_genes'] = G1.shape[0]; R['gse153555_genes'] = G2.shape[0]

r1 = de(G1, gr1, 'Calcified', 'Normal')
r2 = de(G2, gr2, 'AS', 'NC')
R['deg_83453_fdr'] = int((r1['FDR'] < 0.05).sum())
R['deg_83453_strong'] = int(((r1['FDR'] < 0.05) & (r1['log2FC'].abs() > 0.585)).sum())
R['deg_153555_fdr'] = int((r2['FDR'] < 0.05).sum())

# ---- Fig1: volcano ----
fig, ax = plt.subplots(figsize=(7, 5.5), dpi=300)
nl = -np.log10(r1['P.Value'].clip(lower=1e-300))
sg = (r1['P.Value'] < 0.05) & (r1['log2FC'].abs() > 0.585)
ax.scatter(r1.loc[~sg, 'log2FC'], nl[~sg], s=4, c='#B0BEC5', alpha=0.5, label='NS')
ax.scatter(r1.loc[sg & (r1['log2FC'] > 0), 'log2FC'], nl[sg & (r1['log2FC'] > 0)], s=5, c='#C0504D', alpha=0.6, label='Up in calcified')
ax.scatter(r1.loc[sg & (r1['log2FC'] < 0), 'log2FC'], nl[sg & (r1['log2FC'] < 0)], s=5, c='#4F81BD', alpha=0.6, label='Down in calcified')
for g in ['SPP1', 'MMP9', 'CD74', 'MMP12', 'CXCL12', 'CCL19']:
    if g in r1.index:
        ax.annotate(g, (r1.loc[g, 'log2FC'], nl[g]), textcoords='offset points', xytext=(6, 4), fontsize=8)
ax.axvline(0.585, ls='--', c='grey', lw=0.8); ax.axvline(-0.585, ls='--', c='grey', lw=0.8)
ax.axhline(-np.log10(0.05), ls='--', c='grey', lw=0.8)
ax.set_xlabel('log2 fold change (calcified / normal)'); ax.set_ylabel('-log10(P value)')
ax.set_title('GSE83453: calcified vs normal aortic valves')
ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'Fig1_volcano.png')); plt.close(fig)

# ---- 候选基因 ----
deg1 = set(r1[(r1['P.Value'] < 0.05) & (r1['log2FC'].abs() > 0.585)].index)
cand = deg1 & IMMUNE & SECRETED & set(r2.index)
both = r1.loc[list(cand), ['log2FC', 'P.Value']].join(r2.loc[list(cand), ['log2FC', 'P.Value']], lsuffix='_1', rsuffix='_2')
cand2 = both[(np.sign(both['log2FC_1']) == np.sign(both['log2FC_2'])) & (both['P.Value_2'] < 0.05)]
# 固定候选基因顺序：集合运算顺序随进程随机，排序保证 LASSO/SVM-RFE/RF 结果可重复
genes = sorted(cand2.index.tolist())
R['n_deg1'] = len(deg1); R['n_cand'] = len(cand); R['n_cand2'] = len(cand2)
R['cand_genes'] = sorted(genes)

# ---- Fig2: 富集分析 ----
enr = enrichr(genes)
go = enr['GO_Biological_Process_2023'].head(10).iloc[::-1]
fig, ax = plt.subplots(figsize=(7.5, 5), dpi=300)
ax.barh(range(len(go)), -np.log10(go['adjP']), color='#2E6E6A')
ax.set_yticks(range(len(go)))
ax.set_yticklabels([t.split(' (')[0][:55] for t in go['term']], fontsize=8)
ax.set_xlabel('-log10(adjusted P)')
ax.set_title('GO Biological Process enrichment\n(28 replicated candidate genes)', fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'Fig2_GO_enrichment.png')); plt.close(fig)
R['go_top'] = enr['GO_Biological_Process_2023'].head(5)[['term', 'adjP']].to_dict('records')
R['kegg_top'] = enr['KEGG_2021_Human'].head(5)[['term', 'adjP']].to_dict('records')
enr['GO_Biological_Process_2023'].head(30).to_csv(os.path.join(OUT, 'TableS_GO.csv'), index=False, encoding='utf-8-sig')
enr['KEGG_2021_Human'].head(30).to_csv(os.path.join(OUT, 'TableS_KEGG.csv'), index=False, encoding='utf-8-sig')

# ---- 机器学习 ----
X = G1.loc[genes].T.values
y = (gr1 == 'Calcified').astype(int).values
Xs = StandardScaler().fit_transform(X)
n_sel = min(5, len(genes))
lr = LogisticRegressionCV(Cs=20, penalty='l1', solver='liblinear', cv=StratifiedKFold(5), max_iter=5000, random_state=42).fit(Xs, y)
lasso_genes = [g for g, c in zip(genes, lr.coef_[0]) if abs(c) > 1e-6]
rfe = RFE(SVC(kernel='linear', random_state=42), n_features_to_select=n_sel).fit(Xs, y)
svm_genes = [g for g, s in zip(genes, rfe.support_) if s]
rf = RandomForestClassifier(n_estimators=500, random_state=42).fit(Xs, y)
rf_imp = pd.Series(rf.feature_importances_, index=genes).sort_values(ascending=False)
rf_genes = rf_imp.head(n_sel).index.tolist()
votes = Counter(lasso_genes + svm_genes + rf_genes)
hub = sorted(set(lasso_genes) & set(svm_genes) & set(rf_genes)) or sorted([g for g, v in votes.items() if v >= 2])
R['lasso_genes'] = lasso_genes; R['svm_genes'] = svm_genes; R['rf_genes'] = rf_genes; R['hub_genes'] = hub

# Fig3: 三算法选择概览
fig, axes = plt.subplots(1, 3, figsize=(13, 4.6), dpi=300)
ax = axes[0]
lc = pd.Series(lr.coef_[0], index=genes)
lc = lc[lc.abs() > 1e-6].sort_values()
ax.barh(lc.index, lc.values, color=['#C0504D' if v > 0 else '#4F81BD' for v in lc])
ax.set_title(f'LASSO ({len(lc)} genes)', fontsize=10); ax.set_xlabel('Coefficient'); ax.tick_params(labelsize=7)
ax.axvline(0, c='black', lw=0.8)
ax = axes[1]
rk = pd.Series({g: i + 1 for i, g in enumerate([g for g, _ in sorted(zip(genes, rfe.ranking_), key=lambda t: t[1])])})
sel = pd.Series(1 / rfe.ranking_.astype(float), index=genes).sort_values().tail(12)
ax.barh(sel.index, sel.values, color='#2E6E6A')
ax.set_title('SVM-RFE (1/rank, top 12)', fontsize=10); ax.set_xlabel('1 / rank'); ax.tick_params(labelsize=7)
ax = axes[2]
ri = rf_imp.head(12).sort_values()
ax.barh(ri.index, ri.values, color='#E8862E')
ax.set_title('Random Forest importance (top 12)', fontsize=10); ax.set_xlabel('Importance'); ax.tick_params(labelsize=7)
fig.suptitle('Three-algorithm feature selection', fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.93]); fig.savefig(os.path.join(OUT, 'Fig3_ML_selection.png')); plt.close(fig)

# ---- 模型 + ROC ----
trX = G1.T[hub].values
teX = G2.loc[hub].T.values
ytr = y; yte = (gr2 == 'AS').astype(int).values
sc = StandardScaler().fit(trX)
clf = LogisticRegression(penalty='l2', C=1.0, max_iter=5000).fit(sc.transform(trX), ytr)
cv_n = min(5, sum(ytr), sum(1 - ytr))
# 与 upgrade6_stats.py 口径一致：打乱的分层 5 折，固定种子，保证可重复
cv_split = StratifiedKFold(cv_n, shuffle=True, random_state=42)
loo = cross_val_predict(LogisticRegression(penalty='l2', C=1.0, max_iter=5000), sc.transform(trX), ytr, cv=cv_split, method='predict_proba')[:, 1]
R['auc_train'] = float(roc_auc_score(ytr, clf.predict_proba(sc.transform(trX))[:, 1]))
R['auc_cv'] = float(roc_auc_score(ytr, loo))
sc2 = StandardScaler().fit(teX)
score_te = (sc2.transform(teX) * clf.coef_[0]).sum(axis=1)
R['auc_ext'] = float(roc_auc_score(yte, score_te))
R['coef'] = dict(zip(hub, clf.coef_[0].round(3)))

fig, ax = plt.subplots(figsize=(5.6, 5.2), dpi=300)
fpr, tpr, _ = roc_curve(ytr, clf.predict_proba(sc.transform(trX))[:, 1])
ax.plot(fpr, tpr, c='#C0504D', lw=2, label=f"GSE83453 training AUC = {R['auc_train']:.3f}")
fpr, tpr, _ = roc_curve(ytr, loo)
ax.plot(fpr, tpr, c='#E8862E', lw=1.6, ls='-.', label=f"5-fold CV AUC = {R['auc_cv']:.3f}")
fpr, tpr, _ = roc_curve(yte, score_te)
ax.plot(fpr, tpr, c='#4F81BD', lw=2, label=f"GSE153555 replication AUC = {R['auc_ext']:.3f}")
ax.plot([0, 1], [0, 1], ls='--', c='grey')
ax.set_xlabel('1 - specificity'); ax.set_ylabel('Sensitivity')
ax.set_title(f"Diagnostic model: {' + '.join(hub)}")
ax.legend(fontsize=9, loc='lower right')
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'Fig4_ROC.png')); plt.close(fig)

# ---- 免疫浸润 ----
sc1 = immune_scores(G1)
imm_diff = []
for cell in sc1.columns:
    a = sc1.loc[gr1[gr1 == 'Calcified'].index, cell]
    b = sc1.loc[gr1[gr1 == 'Normal'].index, cell]
    p = stats.mannwhitneyu(a, b).pvalue
    imm_diff.append({'cell': cell, 'mean_calc': a.mean(), 'mean_norm': b.mean(), 'P': p})
imm_diff = pd.DataFrame(imm_diff)
imm_diff['FDR'] = bh_fdr(imm_diff['P'].values)
R['immune_infiltration'] = imm_diff.round(4).to_dict('records')
imm_diff.round(4).to_csv(os.path.join(OUT, 'TableS_immune.csv'), index=False, encoding='utf-8-sig')

fig, ax = plt.subplots(figsize=(10, 4.6), dpi=300)
cells = imm_diff.sort_values('P')['cell'].tolist()
data = [[sc1.loc[gr1[gr1 == 'Normal'].index, c], sc1.loc[gr1[gr1 == 'Calcified'].index, c]] for c in cells]
pos = np.arange(len(cells)) * 2
for i, (n, c_) in enumerate(data):
    bp1 = ax.boxplot(n, positions=[pos[i] - 0.4], widths=0.7, patch_artist=True, medianprops=dict(color='black'))
    bp2 = ax.boxplot(c_, positions=[pos[i] + 0.4], widths=0.7, patch_artist=True, medianprops=dict(color='black'))
    bp1['boxes'][0].set_facecolor('#7B9E89'); bp2['boxes'][0].set_facecolor('#C0504D')
    star = '*' if imm_diff.set_index('cell').loc[cells[i], 'FDR'] < 0.05 else ''
    ax.text(pos[i], max(np.max(n), np.max(c_)) + 0.15, star, ha='center', fontsize=12)
ax.set_xticks(pos); ax.set_xticklabels(cells, rotation=35, ha='right', fontsize=8)
ax.set_ylabel('Immune score (z)')
ax.set_title('Immune cell infiltration scores: calcified vs normal valves (GSE83453)')
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color='#7B9E89', label='Normal'), Patch(color='#C0504D', label='Calcified')], fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'Fig5_immune.png')); plt.close(fig)

# hub 与免疫细胞相关性
corr = {}
for g in hub:
    ge = G1.loc[g]
    corr[g] = {c: stats.spearmanr(ge.values, sc1[c].values).statistic for c in sc1.columns}
corr = pd.DataFrame(corr)
R['hub_immune_corr'] = corr.round(3).to_dict()

# ---- Fig6: hub 基因表达 ----
fig, axes = plt.subplots(1, len(hub) * 2, figsize=(3.2 * len(hub) * 2, 4.2), dpi=300)
for j, g in enumerate(hub):
    for k, (G_, gr_, case, ctrl, title) in enumerate([(G1, gr1, 'Calcified', 'Normal', 'GSE83453'), (G2, gr2, 'AS', 'NC', 'GSE153555')]):
        ax = axes[j * 2 + k]
        d = [G_.loc[g, gr_ == ctrl].values, G_.loc[g, gr_ == case].values]
        bp = ax.boxplot(d, tick_labels=[ctrl, case], patch_artist=True, widths=0.55, medianprops=dict(color='black'))
        bp['boxes'][0].set_facecolor('#7B9E89'); bp['boxes'][1].set_facecolor('#C0504D')
        p = stats.mannwhitneyu(d[0], d[1]).pvalue
        ax.set_title(f'{g} ({title})\nP={p:.1e}', fontsize=9)
        ax.set_ylabel('Expression')
fig.suptitle('Hub gene expression in discovery and validation cohorts', fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.93]); fig.savefig(os.path.join(OUT, 'Fig6_hub_expression.png')); plt.close(fig)

# hub 基因统计明细
hubstats = []
for g in hub:
    hubstats.append({'gene': g,
                     'log2FC_GSE83453': r1.loc[g, 'log2FC'], 'P_GSE83453': r1.loc[g, 'P.Value'], 'FDR_GSE83453': r1.loc[g, 'FDR'],
                     'log2FC_GSE153555': r2.loc[g, 'log2FC'], 'P_GSE153555': r2.loc[g, 'P.Value'], 'FDR_GSE153555': r2.loc[g, 'FDR'],
                     'coef': clf.coef_[0][hub.index(g)]})
hubstats = pd.DataFrame(hubstats).round(5)
hubstats.to_csv(os.path.join(OUT, 'Table_hub_genes.csv'), index=False, encoding='utf-8-sig')
R['hubstats'] = hubstats.to_dict('records')

json.dump(R, open(os.path.join(OUT, 'results_B.json'), 'w', encoding='utf-8'), indent=1, ensure_ascii=False, default=str)
print(json.dumps({k: v for k, v in R.items() if not isinstance(v, (list, dict))}, indent=1, ensure_ascii=False))
print('B 流水线完成:', OUT)
