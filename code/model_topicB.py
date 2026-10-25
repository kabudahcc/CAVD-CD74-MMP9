# -*- coding: utf-8 -*-
"""选题B正式建模预演：免疫相关分泌蛋白 ∩ DEG → LASSO+SVM-RFE+RF → 诊断模型 → GSE153555 外部验证"""
import gzip, io, os, json, re
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
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA = r'C:\Users\98757\Desktop\生信分析\data'
OUT = r'C:\Users\98757\Desktop\生信分析\预分析_选题B'
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

def bh_fdr(p):
    p = np.asarray(p, float); n = len(p)
    o = np.argsort(p); r = p[o] * n / (np.arange(n) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    f = np.empty(n); f[o] = np.clip(r, 0, 1); return f

def de(g, gr, case, ctrl):
    a = g.loc[:, gr == case]; b = g.loc[:, gr == ctrl]
    logfc = a.mean(axis=1) - b.mean(axis=1)
    p = stats.ttest_ind(a.values, b.values, axis=1, equal_var=False, nan_policy='omit').pvalue
    r = pd.DataFrame({'log2FC': logfc, 'P.Value': p}).dropna()
    r['FDR'] = bh_fdr(r['P.Value'].values)
    return r

lists = json.load(open(os.path.join(DATA, 'gene_lists.json'), encoding='utf-8'))
IMMUNE = set(lists['immune_kegg']); SECRETED = set(lists['secreted_uniprot'])

# ---------- 重建两个基因矩阵 ----------
# GSE83453
titles_line = None
with gzip.open(os.path.join(DATA, 'GSE83453_series_matrix.txt.gz'), 'rt', encoding='utf-8', errors='ignore') as f:
    buf = []
    for line in f:
        if line.startswith('!Sample_title'):
            titles_line = line.strip()
        if line.startswith('"ID_REF"') or line.startswith('ID_REF'):
            buf.append(line); break
    for line in f:
        if line.startswith('!series_matrix_table_end'):
            break
        buf.append(line)
expr = pd.read_csv(io.StringIO(''.join(buf)), sep='\t', index_col=0)
expr.index = expr.index.astype(str)
titles = [t.strip().strip('"') for t in titles_line.split('\t')[1:]]
g1_gr = pd.Series(['Calcified' if ('bicuspid' in t or 'stenotic' in t) else 'Normal' for t in titles], index=expr.columns)
tbl = []
hdr = False
with gzip.open(os.path.join(DATA, 'GPL10558.annot.gz'), 'rt', encoding='utf-8', errors='ignore') as f:
    for l in f:
        if not hdr:
            if l.count('\t') > 10 and l.split('\t')[0].strip() == 'ID':
                hdr = True; tbl.append(l)
        elif l.count('\t') > 5:
            tbl.append(l)
annot = pd.read_csv(io.StringIO(''.join(tbl)), sep='\t', low_memory=False)
p2s = dict(zip(annot['ID'].astype(str), annot['Gene symbol'].astype(str)))
p2s = {k: v for k, v in p2s.items() if v not in ('nan', '', '---')}
expr = expr.loc[expr.index.isin(p2s)]
expr['symbol'] = expr.index.map(p2s)
expr['rowmean'] = expr.drop(columns=['symbol']).mean(axis=1)
expr = expr.sort_values('rowmean', ascending=False).drop_duplicates('symbol').drop(columns='rowmean')
G1 = expr.set_index('symbol')

# GSE153555
rc = pd.read_csv(os.path.join(DATA, 'GSE153555_readcount_genename.txt.gz'), sep='\t', compression='gzip')
samp = [c for c in rc.columns if c.startswith('AVLF_')]
rc = rc[rc['Gene_name'].notna() & (rc['Gene_name'].astype(str).str.strip() != '')]
bio = {}
for c in samp:
    bio.setdefault(re.sub(r'[ab]$', '', c), []).append(c)
m2 = pd.DataFrame({b: rc[cols].sum(axis=1).values for b, cols in bio.items()}, index=rc['Gene_name'].astype(str))
g2_gr = pd.Series(['NC' if c.startswith('AVLF_D') else 'AS' if c.startswith('AVLF_S') else 'AI' for c in m2.columns], index=m2.columns)
cpm = m2.div(m2.sum(axis=0), axis=1) * 1e6
G2 = np.log2(cpm.loc[(cpm >= 1).sum(axis=1) >= 4] + 0.5).groupby(level=0).mean()

# ---------- 差异分析 ----------
r1 = de(G1, g1_gr, 'Calcified', 'Normal')
r2 = de(G2, g2_gr, 'AS', 'NC')

# ---------- 候选基因：DEG ∩ 免疫 ∩ 分泌 ∩ 双数据集重复 ----------
deg1 = set(r1[(r1['P.Value'] < 0.05) & (r1['log2FC'].abs() > 0.585)].index)
cand = deg1 & IMMUNE & SECRETED & set(r2.index)
both = r1.loc[list(cand), ['log2FC','P.Value']].join(r2.loc[list(cand), ['log2FC','P.Value']], lsuffix='_1', rsuffix='_2')
cand2 = both[(np.sign(both['log2FC_1']) == np.sign(both['log2FC_2'])) & (both['P.Value_2'] < 0.05)]
print(f'DEG {len(deg1)} ∩ 免疫 {len(IMMUNE)} ∩ 分泌 {len(SECRETED)} = {len(cand)}；双数据集重复 = {len(cand2)}')
print('候选基因:', ', '.join(cand2.index.tolist()))

# ---------- 机器学习三算法 ----------
# 固定候选基因顺序：集合运算顺序随进程随机，排序保证三算法结果可重复
genes = sorted(cand2.index.tolist())
X = G1.loc[genes].T.values
y = (g1_gr == 'Calcified').astype(int).values
scaler = StandardScaler().fit(X)
Xs = scaler.transform(X)
n_sel = min(5, len(genes))

lr = LogisticRegressionCV(Cs=20, penalty='l1', solver='liblinear', cv=StratifiedKFold(5), max_iter=5000, random_state=42).fit(Xs, y)
lasso_genes = [g for g, c in zip(genes, lr.coef_[0]) if abs(c) > 1e-6]
print('LASSO 选中:', lasso_genes)

rfe = RFE(SVC(kernel='linear', random_state=42), n_features_to_select=n_sel).fit(Xs, y)
svm_genes = [g for g, s in zip(genes, rfe.support_) if s]
print('SVM-RFE 选中:', svm_genes)

rf = RandomForestClassifier(n_estimators=500, random_state=42).fit(Xs, y)
rf_genes = [g for g, _ in sorted(zip(genes, rf.feature_importances_), key=lambda t: -t[1])[:n_sel]]
print('RF 重要性前5:', rf_genes)

hub = sorted(set(lasso_genes) & set(svm_genes) & set(rf_genes))
if not hub:
    # 票数>=2
    from collections import Counter
    votes = Counter(lasso_genes + svm_genes + rf_genes)
    hub = sorted([g for g, v in votes.items() if v >= 2])
print('Hub 基因（三算法交集/多数票）:', hub)

# ---------- 诊断模型与 AUC ----------
def eval_model(train_g, train_Xdf, train_y, test_g, test_Xdf, test_y, label):
    clf = LogisticRegression(penalty='l2', C=1.0, max_iter=5000)
    sc = StandardScaler().fit(train_Xdf[train_g].values)
    clf.fit(sc.transform(train_Xdf[train_g].values), train_y)
    # 内部 LOOCV
    loo = cross_val_predict(LogisticRegression(penalty='l2', C=1.0, max_iter=5000),
                            sc.transform(train_Xdf[train_g].values), train_y,
                            cv=min(len(train_y), sum(train_y), sum(1-train_y)), method='predict_proba')[:, 1]
    auc_cv = roc_auc_score(train_y, loo)
    # 外部验证（在验证集上独立标准化）
    common = [g for g in train_g if g in test_Xdf.index]
    sc2 = StandardScaler().fit(test_Xdf.loc[common].T.values)
    # 用训练集系数在验证集上打分（方向对齐：各自 z 标准化）
    coef = pd.Series(clf.coef_[0], index=train_g)
    score = (sc2.transform(test_Xdf.loc[common].T.values) * coef[common].values).sum(axis=1) + 0
    auc_ext = roc_auc_score(test_y, score)
    return auc_cv, auc_ext, coef

trX = G1.T; trY = (g1_gr == 'Calcified').astype(int).values
teX = G2; teY = (g2_gr == 'AS').astype(int).values
auc_cv, auc_ext, coef = eval_model(hub, trX, trY, hub, teX, teY, 'B')
print(f'模型 hub 基因 {len(hub)} 个: 内部LOOCV AUC={auc_cv:.3f}, 外部验证(GSE153555) AUC={auc_ext:.3f}')

# ---------- 图 ----------
fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), dpi=150)
# ROC
ax = axes[0]
sc = StandardScaler().fit(trX[hub].values)
clf = LogisticRegression(penalty='l2', C=1.0, max_iter=5000).fit(sc.transform(trX[hub].values), trY)
p_tr = clf.predict_proba(sc.transform(trX[hub].values))[:, 1]
fpr, tpr, _ = roc_curve(trY, p_tr)
ax.plot(fpr, tpr, c='#C0504D', lw=2, label=f'GSE83453 训练集 AUC={roc_auc_score(trY, p_tr):.3f}')
common = [g for g in hub if g in teX.index]
sc2 = StandardScaler().fit(teX.loc[common].T.values)
score = (sc2.transform(teX.loc[common].T.values) * coef[common].values).sum(axis=1)
fpr2, tpr2, _ = roc_curve(teY, score)
ax.plot(fpr2, tpr2, c='#4F81BD', lw=2, label=f'GSE153555 外部验证 AUC={auc_ext:.3f}')
ax.plot([0,1],[0,1], ls='--', c='grey')
ax.set_xlabel('1 - 特异度'); ax.set_ylabel('灵敏度'); ax.set_title('诊断模型 ROC 曲线')
ax.legend(fontsize=9)
# 系数
ax = axes[1]
cc = coef.sort_values()
colors = ['#C0504D' if v > 0 else '#4F81BD' for v in cc]
ax.barh(cc.index, cc.values, color=colors)
ax.axvline(0, c='black', lw=0.8)
ax.set_xlabel('Logistic 回归系数'); ax.set_title(f'Hub 基因系数（{len(hub)} 个）')
fig.tight_layout(); fig.savefig(os.path.join(OUT, '图B3_机器学习模型ROC与系数.png')); plt.close(fig)

with pd.ExcelWriter(os.path.join(OUT, '选题B建模结果.xlsx'), engine='openpyxl') as w:
    cand2.round(6).to_excel(w, sheet_name='候选基因_免疫分泌双重复')
    pd.DataFrame({'gene': hub, 'coef': coef[hub].values}).to_excel(w, sheet_name='Hub基因与系数', index=False)
    pd.DataFrame({'LASSO': pd.Series(lasso_genes), 'SVM_RFE': pd.Series(svm_genes), 'RF': pd.Series(rf_genes)}).to_excel(w, sheet_name='三算法选中基因', index=False)
print('输出目录:', OUT)
