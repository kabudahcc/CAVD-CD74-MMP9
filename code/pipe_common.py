# -*- coding: utf-8 -*-
"""公共分析模块：数据加载、差异分析、富集分析、免疫浸润评分"""
import gzip, io, os, json, re, urllib.request
import numpy as np
import pandas as pd
from scipy import stats

DATA = r'C:\Users\98757\Desktop\生信分析\data'

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

def read_series_matrix(path):
    titles_line = None
    with gzip.open(path, 'rt', encoding='utf-8', errors='ignore') as f:
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
    return expr, titles

def collapse_probes(expr, p2s):
    expr = expr.loc[expr.index.isin(p2s)].copy()
    expr['symbol'] = expr.index.map(p2s)
    expr['rowmean'] = expr.drop(columns=['symbol']).mean(axis=1)
    expr = expr.sort_values('rowmean', ascending=False).drop_duplicates('symbol').drop(columns='rowmean')
    return expr.set_index('symbol')

def load_gse83453():
    expr, titles = read_series_matrix(os.path.join(DATA, 'GSE83453_series_matrix.txt.gz'))
    gr = pd.Series(['Calcified' if ('bicuspid' in t or 'stenotic' in t) else 'Normal' for t in titles], index=expr.columns)
    tbl = []; hdr = False
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
    G = collapse_probes(expr, p2s)
    return G, gr  # Calcified / Normal

def load_gse153555():
    rc = pd.read_csv(os.path.join(DATA, 'GSE153555_readcount_genename.txt.gz'), sep='\t', compression='gzip')
    samp = [c for c in rc.columns if c.startswith('AVLF_')]
    rc = rc[rc['Gene_name'].notna() & (rc['Gene_name'].astype(str).str.strip() != '')]
    bio = {}
    for c in samp:
        bio.setdefault(re.sub(r'[ab]$', '', c), []).append(c)
    m2 = pd.DataFrame({b: rc[cols].sum(axis=1).values for b, cols in bio.items()}, index=rc['Gene_name'].astype(str))
    gr = pd.Series(['NC' if c.startswith('AVLF_D') else 'AS' if c.startswith('AVLF_S') else 'AI' for c in m2.columns], index=m2.columns)
    cpm = m2.div(m2.sum(axis=0), axis=1) * 1e6
    G = np.log2(cpm.loc[(cpm >= 1).sum(axis=1) >= 4] + 0.5).groupby(level=0).mean()
    return G, gr  # NC / AS / AI

def load_gse66360():
    expr, titles = read_series_matrix(os.path.join(DATA, 'GSE66360_series_matrix.txt.gz'))
    groups = pd.Series(['AMI' if 'Myocardial Infarction' in t else 'Control' for t in titles], index=expr.columns)
    cohorts = pd.Series(['DISCOVERY' if 'DISCOVERY' in t else 'VALIDATION' for t in titles], index=expr.columns)
    if np.nanmax(expr.values) > 50:
        expr = np.log2(expr + 1)
    p2s = json.load(open(os.path.join(DATA, 'gpl570_probe2symbol.json'), encoding='utf-8'))
    G = collapse_probes(expr, p2s)
    return G, groups, cohorts

def enrichr(genes, libs=('GO_Biological_Process_2023', 'KEGG_2021_Human')):
    boundary = '----kimiboundary'
    body = (f'--{boundary}\r\nContent-Disposition: form-data; name="list"\r\n\r\n' + '\n'.join(genes) +
            f'\r\n--{boundary}\r\nContent-Disposition: form-data; name="description"\r\n\r\nanalysis\r\n--{boundary}--\r\n').encode()
    req = urllib.request.Request('https://maayanlab.cloud/Enrichr/addList', data=body,
                                 headers={'Content-Type': f'multipart/form-data; boundary={boundary}', 'User-Agent': 'Mozilla/5.0'})
    uid = json.loads(urllib.request.urlopen(req, timeout=60).read())['userListId']
    out = {}
    for lib in libs:
        url = f'https://maayanlab.cloud/Enrichr/enrich?userListId={uid}&backgroundType={lib}'
        res = json.loads(urllib.request.urlopen(urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'}), timeout=60).read())
        rows = [{'term': r[1], 'p': r[2], 'adjP': r[6], 'genes': ','.join(r[5])} for r in res.get(lib, [])]
        out[lib] = pd.DataFrame(rows)
    return out

# 免疫细胞标志物（主流谱系，文献通用标志基因）
IMMUNE_MARKERS = {
    'CD8+ T cells': ['CD8A', 'CD8B', 'GZMB', 'PRF1'],
    'CD4+ T cells': ['CD4', 'IL7R', 'LTB'],
    'Tregs': ['FOXP3', 'IL2RA', 'CTLA4'],
    'B cells': ['CD19', 'MS4A1', 'CD79A'],
    'NK cells': ['NCAM1', 'NKG7', 'KLRD1', 'GNLY'],
    'Macrophages': ['CD68', 'CD163', 'MRC1'],
    'Monocytes': ['CD14', 'FCGR3A', 'VCAN'],
    'Neutrophils': ['S100A8', 'S100A9', 'FCGR3B', 'CSF3R'],
    'Dendritic cells': ['CD209', 'HLA-DRA', 'CLEC9A'],
    'Mast cells': ['TPSAB1', 'TPSB2', 'KIT'],
    'Th1 cells': ['IFNG', 'TBX21'],
    'Th2 cells': ['GATA3', 'IL4', 'IL5'],
    'Th17 cells': ['IL17A', 'RORC'],
}

def immune_scores(G):
    """每种免疫细胞：成员基因 z 分数均值（ssGSEA 简化替代）"""
    Z = (G.T - G.mean(axis=1)) / G.std(axis=1).replace(0, np.nan)
    Z = Z.T
    scores = {}
    for cell, genes in IMMUNE_MARKERS.items():
        present = [g for g in genes if g in Z.index]
        if present:
            scores[cell] = Z.loc[present].mean(axis=0)
    return pd.DataFrame(scores)
