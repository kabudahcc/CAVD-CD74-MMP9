# -*- coding: utf-8 -*-
"""cis-MR v2：GTEx eQTL (CD74/MMP9) → FinnGen CAVS"""
import sys, os, json, math
sys.path.insert(0, r'C:\Users\98757\Desktop\生信分析\data')
import numpy as np
import pandas as pd
from scipy import stats
from tabix_remote import tabix_query, fetch_header

GTEX = {
    'Whole blood': 'https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge/Whole_Blood.tsv.gz',
    'Artery aorta': 'https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge/Artery_Aorta.tsv.gz',
}
FINNGEN = 'https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_I9_CAVS_OPERATED.gz'
GENES = {'CD74': ('ENSG00000019582', '5', 150_150_000, 150_670_000),
         'MMP9': ('ENSG00000100985', '20', 44_350_000, 44_910_000)}
CACHE = r'C:\Users\98757\Desktop\生信分析\data\mr'
os.makedirs(CACHE, exist_ok=True)

def parse_rows(url, chrom, beg, end, colmap, filter_fn, cache_name):
    """colmap: {字段名: 语义} 依表头定位；filter_fn(f)->bool"""
    cf = os.path.join(CACHE, cache_name)
    if os.path.exists(cf):
        return open(cf, encoding='utf-8').read().splitlines()
    hdr = fetch_header(url)
    idx = {c: i for i, c in enumerate(hdr)}
    need = {k: idx[v] for k, v in colmap.items()}
    pos_i = need['pos']
    lines = tabix_query(url, chrom, beg, end, nfields=len(hdr))
    keep = []
    for ln in lines:
        f = ln.split('\t')
        try:
            p = int(f[pos_i])
        except Exception:
            continue
        if not (beg <= p <= end):
            continue
        if filter_fn(f, need):
            keep.append(ln)
    open(cf, 'w', encoding='utf-8').write('\n'.join(keep))
    return keep

GTEX_COLS = {'gene': 'molecular_trait_id', 'chrom': 'chromosome', 'pos': 'position', 'ref': 'ref',
             'alt': 'alt', 'beta': 'beta', 'se': 'se', 'p': 'pvalue', 'maf': 'maf', 'rsid': 'rsid'}

def gtex(tis, gene, ens, chrom, beg, end):
    url = GTEX[tis]
    hdr = fetch_header(url)
    idx = {c: i for i, c in enumerate(hdr)}
    need = {k: idx[v] for k, v in GTEX_COLS.items()}
    lines = parse_rows(url, chrom, beg, end, GTEX_COLS, lambda f, n: f[n['gene']] == ens,
                       f'gtex_{tis.replace(" ","_")}_{gene}.txt')
    rows = []
    for ln in lines:
        f = ln.split('\t')
        rows.append({'rsid': f[need['rsid']], 'chrom': f[need['chrom']].replace('chr', ''), 'pos': int(f[need['pos']]),
                     'ref': f[need['ref']].upper(), 'alt': f[need['alt']].upper(),
                     'beta_e': float(f[need['beta']]), 'se_e': float(f[need['se']]),
                     'p_e': float(f[need['p']]), 'maf_e': float(f[need['maf']])})
    return pd.DataFrame(rows).drop_duplicates(subset=['chrom', 'pos', 'ref', 'alt'])

FG_COLS = {'chrom': '#chrom', 'pos': 'pos', 'ref': 'ref', 'alt': 'alt', 'rsids': 'rsids',
           'beta': 'beta', 'se': 'sebeta', 'p': 'pval'}

def finngen(chrom, beg, end):
    hdr = fetch_header(FINNGEN)
    idx = {c: i for i, c in enumerate(hdr)}
    need = {k: idx[v] for k, v in FG_COLS.items()}
    lines = parse_rows(FINNGEN, chrom, beg, end, FG_COLS, lambda f, n: True,
                       f'finngen_{chrom}_{beg}.txt')
    out = {}
    for ln in lines:
        f = ln.split('\t')
        out[(f[need['chrom']].replace('chr', ''), int(f[need['pos']]))] = {
            'rsids': f[need['rsids']], 'ref': f[need['ref']].upper(), 'alt': f[need['alt']].upper(),
            'beta_o': float(f[need['beta']]), 'se_o': float(f[need['se']]), 'p_o': float(f[need['p']])}
    return out

if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if which == 'gtex':
        for gene, (ens, chrom, beg, end) in GENES.items():
            for tis in GTEX:
                df = gtex(tis, gene, ens, chrom, beg, end)
                print(gene, tis, len(df), 'rows; minP =', f"{df.p_e.min():.2e}" if len(df) else 'NA')
    elif which == 'finngen':
        for gene, (ens, chrom, beg, end) in GENES.items():
            o = finngen(chrom, beg, end)
            print(gene, len(o), 'outcome rows')
