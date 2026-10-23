# -*- coding: utf-8 -*-
"""cis-MR：GTEx eQTL (CD74/MMP9) → FinnGen CAVS(I9_CAVS_OPERATED)
数据源：EBI eQTL Catalogue GTEx v8（含 beta/se/rsid）+ FinnGen R12 tabix 远程区域抓取
"""
import sys, os, json, math
sys.path.insert(0, r'C:\Users\98757\Desktop\生信分析\data')
import numpy as np
import pandas as pd
from scipy import stats
from tabix_remote import tabix_query, fetch_header
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'Arial'
OUT = r'C:\Users\98757\Desktop\生信分析\论文_选题B'
os.makedirs(OUT, exist_ok=True)

GTEX = {
    'Whole blood': 'https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge/Whole_Blood.tsv.gz',
    'Artery aorta': 'https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge/Artery_Aorta.tsv.gz',
}
FINNGEN = 'https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_I9_CAVS_OPERATED.gz'
GENES = {'CD74': ('ENSG00000019582', '5', 149_400_000, 151_420_000),
         'MMP9': ('ENSG00000100985', '20', 43_600_000, 45_660_000)}

# ---------- 1. 取 GTEx 表头并查询 ----------
def gtex_rows(url, ens, chrom, beg, end):
    hdr = fetch_header(url)
    idx = {c: i for i, c in enumerate(hdr)}
    rows = []
    for ln in tabix_query(url, chrom, beg, end):
        f = ln.split('\t')
        if len(f) < len(hdr):
            continue
        if f[idx['molecular_trait_id']] == ens:
            rows.append({'rsid': f[idx['rsid']], 'chrom': f[idx['chromosome']], 'pos': int(f[idx['position']]),
                         'ref': f[idx['ref']].upper(), 'alt': f[idx['alt']].upper(),
                         'beta_e': float(f[idx['beta']]), 'se_e': float(f[idx['se']]),
                         'p_e': float(f[idx['pvalue']]), 'maf_e': float(f[idx['maf']])})
    return pd.DataFrame(rows).drop_duplicates(subset=['chrom', 'pos', 'ref', 'alt'])

exp_all = {}
for gene, (ens, chrom, beg, end) in GENES.items():
    for tis, url in GTEX.items():
        try:
            df = gtex_rows(url, ens, chrom, beg, end)
        except Exception as e:
            print('GTEx 查询失败', gene, tis, e); df = pd.DataFrame()
        exp_all[(gene, tis)] = df
        print(f'{gene} / {tis}: {len(df)} eQTL rows, min P = {df.p_e.min() if len(df) else float("nan"):.2e}')

# ---------- 2. 工具变量选择（P 阈值 + 距离剪枝）----------
def select_iv(df, pth=5e-8, dist=250_000, maxn=12):
    for p in [pth, 1e-6, 1e-5]:
        d = df[df['p_e'] < p].sort_values('p_e')
        if len(d) >= 2 or p == 1e-5:
            break
    kept = []
    for _, r in d.iterrows():
        if all(abs(r['pos'] - k['pos']) > dist for k in kept):
            kept.append(r)
        if len(kept) >= maxn:
            break
    return pd.DataFrame(kept)

ivs = {k: select_iv(v) for k, v in exp_all.items()}
for k, v in ivs.items():
    print(k, '->', len(v), 'instruments')

# ---------- 3. FinnGen 结局数据 ----------
fg_hdr = fetch_header(FINNGEN)
print('FinnGen cols:', fg_hdr)
fidx = {c: i for i, c in enumerate(fg_hdr)}

def fg_rows(chrom, beg, end):
    rows = {}
    for ln in tabix_query(FINNGEN, chrom, beg, end):
        f = ln.split('\t')
        if len(f) < 10:
            continue
        rsids = f[fidx['rsids']] if 'rsids' in fidx else ''
        key2 = (f[fidx['#chrom']].replace('chr', ''), int(f[fidx['pos']]))
        rows[key2] = {'rsids': rsids, 'ref': f[fidx['ref']].upper(), 'alt': f[fidx['alt']].upper(),
                      'beta_o': float(f[fidx['beta']]), 'se_o': float(f[fidx['sebeta']]),
                      'p_o': float(f[fidx['pval']])}
    return rows

outcome = {}
for gene, (ens, chrom, beg, end) in GENES.items():
    outcome[gene] = fg_rows(chrom, beg, end)
    print(gene, 'FinnGen rows:', len(outcome[gene]))

# ---------- 4. 谐和 + MR ----------
def comp(a, b):
    m = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(m[x] for x in a) == b

def harmonize(iv, outmap):
    recs = []
    for _, r in iv.iterrows():
        c = outmap.get((r['chrom'].replace('chr', ''), r['pos']))
        if c is None:
            continue
        ea, oa = r['alt'], None
        if r['alt'] == c['alt'] and r['ref'] == c['ref']:
            beta_o = c['beta_o']
        elif r['alt'] == c['ref'] and r['ref'] == c['alt']:
            beta_o = -c['beta_o']
        elif comp(r['alt'], c['alt']) and comp(r['ref'], c['ref']):
            beta_o = c['beta_o']
        elif comp(r['alt'], c['ref']) and comp(r['ref'], c['alt']):
            beta_o = -c['beta_o']
        else:
            continue
        # 回文 SNP 且频率接近 0.5 → 丢弃
        if {r['ref'], r['alt']} in ({'A', 'T'}, {'C', 'G'}) and 0.4 < r['maf_e'] < 0.6:
            continue
        recs.append({**r.to_dict(), 'beta_o': beta_o, 'se_o': c['se_o'], 'p_o': c['p_o']})
    return pd.DataFrame(recs)

mr_results = []
for (gene, tis), iv in ivs.items():
    if not len(iv):
        continue
    h = harmonize(iv, outcome[gene])
    if not len(h):
        print(gene, tis, 'harmonize 后无 SNP'); continue
    h['wald'] = h['beta_o'] / h['beta_e']
    h['wald_se'] = np.sqrt(h['se_o'] ** 2 / h['beta_e'] ** 2 + h['beta_o'] ** 2 * h['se_e'] ** 2 / h['beta_e'] ** 4)
    if len(h) == 1:
        b, se = h['wald'].iloc[0], h['wald_se'].iloc[0]
    else:  # IVW 固定效应
        w = 1 / h['wald_se'] ** 2
        b = (w * h['wald']).sum() / w.sum()
        se = math.sqrt(1 / w.sum())
    z = b / se
    p = 2 * stats.norm.sf(abs(z))
    mr_results.append({'gene': gene, 'tissue': tis, 'n_snp': len(h), 'beta': b, 'se': se,
                       'OR': math.exp(b), 'OR_lo': math.exp(b - 1.96 * se), 'OR_hi': math.exp(b + 1.96 * se), 'P': p})
    h.to_csv(os.path.join(OUT, f'TableS_MR_SNPs_{gene}_{tis.replace(" ", "_")}.csv'), index=False, encoding='utf-8-sig')
    print(gene, tis, f'n={len(h)} OR={math.exp(b):.3f} ({math.exp(b-1.96*se):.3f}-{math.exp(b+1.96*se):.3f}) P={p:.2e}')

mr = pd.DataFrame(mr_results)
mr.to_csv(os.path.join(OUT, 'Table_MR_results.csv'), index=False, encoding='utf-8-sig')
json.dump(mr_results, open(os.path.join(OUT, 'mr_results.json'), 'w'), indent=1)

# ---------- 5. 森林图 ----------
if len(mr):
    fig, ax = plt.subplots(figsize=(7.5, 3.6), dpi=300)
    mr = mr.iloc[::-1].reset_index(drop=True)
    for i, r in mr.iterrows():
        col = '#C0504D' if r['P'] < 0.05 else '#7F7F7F'
        ax.errorbar(r['OR'], i, xerr=[[r['OR'] - r['OR_lo']], [r['OR_hi'] - r['OR']]], fmt='o', color=col, capsize=3, ms=6)
    ax.axvline(1, ls='--', c='grey')
    ax.set_xscale('log')
    ax.set_yticks(range(len(mr)))
    ax.set_yticklabels([f"{r['gene']} ({r['tissue']}, n SNP={r['n_snp']})" for _, r in mr.iterrows()], fontsize=9)
    ax.set_xlabel('Odds ratio for CAVS per 1-SD higher gene expression (95% CI, log scale)')
    ax.set_title('cis-MR: genetically proxied expression and calcific aortic valve stenosis')
    for i, r in mr.iterrows():
        ax.text(ax.get_xlim()[1] * 1.05, i, f"OR={r['OR']:.2f}, P={r['P']:.1e}", va='center', fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'Fig7_MR_forest.png'), bbox_inches='tight')
    plt.close(fig)
    print('Fig7 已保存')
