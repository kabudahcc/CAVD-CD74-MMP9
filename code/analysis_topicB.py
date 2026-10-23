# -*- coding: utf-8 -*-
"""选题B预分析：CAVD（钙化性主动脉瓣疾病）差异分析 + 跨数据集验证
发现集 GSE83453（芯片：钙化瓣19 vs 正常瓣8）；验证集 GSE153555（RNA-seq：AS6 vs NC5）"""
import gzip, io, os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA = r'C:\Users\98757\Desktop\生信分析\data'
OUT = r'C:\Users\98757\Desktop\生信分析\预分析_选题B'
os.makedirs(OUT, exist_ok=True)
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
    return r.sort_values('FDR')

# ---------- 1. GSE83453 发现集 ----------
print('=== GSE83453 发现集 ===')
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
groups1 = pd.Series(['Calcified' if ('bicuspid' in t or 'stenotic' in t) else 'Normal' for t in titles], index=expr.columns)
print('分组:', groups1.value_counts().to_dict())
if np.nanmax(expr.values) > 50:
    expr = np.log2(expr + 1); print('已 log2 转换')
else:
    print('已为 log 尺度, max=%.1f' % np.nanmax(expr.values))

# GPL10558 注释
tbl_lines = []
header_found = False
with gzip.open(os.path.join(DATA, 'GPL10558.annot.gz'), 'rt', encoding='utf-8', errors='ignore') as f:
    for l in f:
        if not header_found:
            if l.count('\t') > 10 and l.split('\t')[0].strip() == 'ID':
                header_found = True
                tbl_lines.append(l)
        else:
            if l.count('\t') > 5:
                tbl_lines.append(l)
annot = pd.read_csv(io.StringIO(''.join(tbl_lines)), sep='\t', low_memory=False)
symcol = None
for c in annot.columns:
    if c.lower() in ('symbol', 'gene_symbol', 'ilmn_gene', 'gene symbol'):
        symcol = c if annot[c].astype(str).str.contains(r'^[A-Z0-9\-\.]+$', na=False).mean() > 0.5 else symcol
    if c.lower() == 'symbol':
        symcol = 'Symbol'
print('注释列:', annot.columns.tolist()[:12], '... 使用:', symcol)
idcol = annot.columns[0]
p2s = dict(zip(annot[idcol].astype(str), annot[symcol].astype(str)))
p2s = {k: v for k, v in p2s.items() if v not in ('nan', '', '---')}

expr = expr.loc[expr.index.isin(p2s)]
expr['symbol'] = expr.index.map(p2s)
expr['rowmean'] = expr.drop(columns=['symbol']).mean(axis=1)
expr = expr.sort_values('rowmean', ascending=False).drop_duplicates('symbol').drop(columns='rowmean')
g1 = expr.set_index('symbol')
print('基因水平矩阵:', g1.shape)

res1 = de(g1, groups1, 'Calcified', 'Normal')
print(f'GSE83453 钙化 vs 正常: FDR<0.05 {(res1["FDR"]<0.05).sum()} 个；P<0.05 {(res1["P.Value"]<0.05).sum()} 个；'
      f'强差异(FDR<0.05且|FC|>0.585) {((res1["FDR"]<0.05)&(res1["log2FC"].abs()>0.585)).sum()} 个')
res1.to_csv(os.path.join(OUT, 'DEG_钙化vs正常_GSE83453.csv'), encoding='utf-8-sig')

# ---------- 2. GSE153555 验证集 ----------
print('\n=== GSE153555 验证集 ===')
rc = pd.read_csv(os.path.join(DATA, 'GSE153555_readcount_genename.txt.gz'), sep='\t', compression='gzip')
samp_cols = [c for c in rc.columns if c.startswith('AVLF_')]
rc = rc[rc['Gene_name'].notna() & (rc['Gene_name'].astype(str).str.strip() != '')]
# 合并技术重复（D1a+D1b）
import re
def bio_id(c):
    return re.sub(r'[ab]$', '', c)
bio = {}
for c in samp_cols:
    bio.setdefault(bio_id(c), []).append(c)
mat2 = pd.DataFrame({b: rc[cols].sum(axis=1).values for b, cols in bio.items()}, index=rc['Gene_name'].astype(str))
groups2 = pd.Series(['NC' if c.startswith('AVLF_D') else 'AS' if c.startswith('AVLF_S') else 'AI' for c in mat2.columns], index=mat2.columns)
print('生物学样本分组:', groups2.value_counts().to_dict())
# CPM + log2
lib = mat2.sum(axis=0)
cpm = mat2.div(lib, axis=1) * 1e6
keep = (cpm >= 1).sum(axis=1) >= 4
g2 = np.log2(cpm.loc[keep] + 0.5)
g2 = g2.groupby(g2.index).mean()
print('过滤后基因数:', g2.shape[0])

res2 = de(g2, groups2, 'AS', 'NC')
print(f'GSE153555 AS vs NC: FDR<0.05 {(res2["FDR"]<0.05).sum()} 个；P<0.05 {(res2["P.Value"]<0.05).sum()} 个；'
      f'强差异 {((res2["FDR"]<0.05)&(res2["log2FC"].abs()>0.585)).sum()} 个')
res2.to_csv(os.path.join(OUT, 'DEG_ASvsNC_GSE153555.csv'), encoding='utf-8-sig')

# ---------- 3. 跨数据集一致性 ----------
both = res1[['log2FC','P.Value','FDR']].join(res2[['log2FC','P.Value','FDR']], lsuffix='_83453', rsuffix='_153555', how='inner')
print('\n共有基因:', len(both))
# 小样本用 P<0.05 作为名义显著
nom = both[(both['P.Value_83453'] < 0.05) & (both['log2FC_83453'].abs() > 0.585)]
conc = (np.sign(nom['log2FC_83453']) == np.sign(nom['log2FC_153555'])).mean()
rep = nom[(np.sign(nom['log2FC_83453'])==np.sign(nom['log2FC_153555'])) & (nom['P.Value_153555']<0.05)]
print(f'发现集名义显著且强差异基因: {len(nom)} 个；方向一致率 {conc*100:.1f}%；双数据集重复 {len(rep)} 个')

# ---------- 4. 免疫分泌蛋白富集度（可行性代理指标）----------
pat = r'^(IL\d|CCL\d|CXCL|CCL|CX3CL1|TNF|IFN|TGFB|C5|C3|C1Q|SFRP|POSTN|SPP1|MMP|ADAMTS|BMP|NOG|CHI3L1|VCAM1|ICAM1|SELE|SELL|SELP)'
immune = [g for g in both.index if re.match(pat, g)]
imm_nom = nom.loc[nom.index.isin(immune)]
imm_rep = rep.loc[rep.index.isin(immune)]
print(f'免疫/分泌相关基因(模式匹配): 全集 {len(immune)} 个；发现集显著强差异 {len(imm_nom)} 个；双数据集重复 {len(imm_rep)} 个')
print('双重复的免疫分泌基因:', ', '.join(imm_rep.index.tolist()[:40]))

# ---------- 5. 图 ----------
fig, ax = plt.subplots(figsize=(8.5, 6.5), dpi=150)
r = res1
nl = -np.log10(r['P.Value'].clip(lower=1e-300))
sg = (r['P.Value']<0.05)&(r['log2FC'].abs()>0.585)
ax.scatter(r.loc[~sg,'log2FC'], nl[~sg], s=4, c='#B0BEC5', alpha=0.5, label='不显著')
ax.scatter(r.loc[sg&(r['log2FC']>0),'log2FC'], nl[sg&(r['log2FC']>0)], s=5, c='#C0504D', alpha=0.6, label='上调')
ax.scatter(r.loc[sg&(r['log2FC']<0),'log2FC'], nl[sg&(r['log2FC']<0)], s=5, c='#4F81BD', alpha=0.6, label='下调')
for g in imm_rep.index[:12]:
    ax.annotate(g, (r.loc[g,'log2FC'], nl[g]), textcoords='offset points', xytext=(6,4), fontsize=8, color='#8a4b00')
ax.axvline(0.585, ls='--', c='grey', lw=0.8); ax.axvline(-0.585, ls='--', c='grey', lw=0.8)
ax.axhline(-np.log10(0.05), ls='--', c='grey', lw=0.8)
ax.set_xlabel('log2 Fold Change（钙化瓣 / 正常瓣）'); ax.set_ylabel('-log10(P)')
ax.set_title('GSE83453：钙化性主动脉瓣 vs 正常瓣膜 差异分析火山图\n（标注＝双数据集重复的免疫分泌相关基因）')
ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(OUT, '图B1_火山图_GSE83453.png')); plt.close(fig)

fig, ax = plt.subplots(figsize=(8, 7), dpi=150)
ax.scatter(nom['log2FC_83453'], nom['log2FC_153555'], s=8, c='#78909C', alpha=0.5, label=f'发现集显著基因 (n={len(nom)})')
ax.scatter(imm_nom['log2FC_83453'], imm_nom['log2FC_153555'], s=45, c='#E8862E', edgecolors='black', linewidths=0.5, zorder=5, label='免疫分泌相关基因')
for g in imm_rep.index[:12]:
    ax.annotate(g, (nom.loc[g,'log2FC_83453'], nom.loc[g,'log2FC_153555']), textcoords='offset points', xytext=(6,4), fontsize=8)
ax.axhline(0,c='grey',lw=0.8); ax.axvline(0,c='grey',lw=0.8)
ax.set_xlabel('log2FC — GSE83453 芯片（钙化19 vs 正常8）')
ax.set_ylabel('log2FC — GSE153555 RNA-seq（AS6 vs NC5）')
ax.set_title(f'选题B跨数据集一致性（钙化性主动脉瓣疾病）\n方向一致率 {conc*100:.1f}%，双重复基因 {len(rep)} 个')
ax.legend(fontsize=9, loc='upper left')
fig.tight_layout(); fig.savefig(os.path.join(OUT, '图B2_跨数据集一致性散点图.png')); plt.close(fig)

rep.sort_values('log2FC_83453').to_csv(os.path.join(OUT, '双数据集重复基因.csv'), encoding='utf-8-sig')
with pd.ExcelWriter(os.path.join(OUT, '选题B预分析结果.xlsx')) as w:
    rep.round(6).to_excel(w, sheet_name='双数据集重复基因')
    imm_rep.round(6).to_excel(w, sheet_name='免疫分泌相关重复基因')
    res1.head(500).round(6).to_excel(w, sheet_name='Top500_GSE83453')
print('输出目录:', OUT)
