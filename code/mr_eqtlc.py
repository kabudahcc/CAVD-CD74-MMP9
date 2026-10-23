# -*- coding: utf-8 -*-
"""从 eQTL Catalogue (r7) 与 GTEx v8 抓取 CD74/MMP9 cis-eQTL 窗口
控制请求频率避免 EBI 封 IP：单线程、每次请求间隔 sleep"""
import os, sys, time, gzip, io
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tabix_remote import tabix_query, fetch_header

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "mr")
os.makedirs(OUT, exist_ok=True)

GENES = {
    "CD74": dict(ens="ENSG00000019582", chrom="5", beg=150_150_000, end=150_660_000),
    "MMP9": dict(ens="ENSG00000100985", chrom="20", beg=45_760_000, end=46_270_000),
}

EQTLC = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/{qts}/{qtd}/{qtd}.all.tsv.gz"
GTEX = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge/{tis}.tsv.gz"

# (标签, url模板参数, 适用基因)
DATASETS = [
    # eQTL Catalogue r7
    ("BLUEPRINT_monocyte",  EQTLC.format(qts="QTS000002", qtd="QTD000021"), "eqtlc", ["MMP9", "CD74"]),
    ("BLUEPRINT_neutrophil",EQTLC.format(qts="QTS000002", qtd="QTD000026"), "eqtlc", ["MMP9"]),
    ("Quach_monocyte",      EQTLC.format(qts="QTS000024", qtd="QTD000409"), "eqtlc", ["MMP9", "CD74"]),
    ("Lepik_blood",         EQTLC.format(qts="QTS000019", qtd="QTD000373"), "eqtlc", ["CD74", "MMP9"]),
    ("OneK1K_CD14Mono",     EQTLC.format(qts="QTS000038", qtd="QTD000609"), "eqtlc", ["CD74", "MMP9"]),
    ("OneK1K_NK",           EQTLC.format(qts="QTS000038", qtd="QTD000620"), "eqtlc", ["CD74"]),
    # GTEx v8 (EBI imported)
    ("GTEx_Artery_Aorta",   GTEX.format(tis="Artery_Aorta"), "gtex", ["MMP9", "CD74"]),
]


def parse_eqtlc(lines, gene):
    g = GENES[gene]
    rows = []
    for ln in lines:
        f = ln.split("\t")
        if len(f) < 19:
            continue
        # molecular_trait_id,chromosome,position,ref,alt,variant,ma_samples,maf,pvalue,beta,se,type,ac,an,r2,molecular_trait_object_id,gene_id,median_tpm,rsid
        if f[0] != g["ens"]:
            continue
        try:
            pos = int(f[2]); p = float(f[8]); beta = float(f[9]); se = float(f[10])
        except ValueError:
            continue
        if not (g["beg"] <= pos <= g["end"]):
            continue
        rows.append(dict(rsid=f[18].strip(), pos=pos, ref=f[3], alt=f[4], beta=beta, se=se,
                         p=p, maf=f[7], an=f[13]))
    return rows


def parse_gtex(lines, gene):
    g = GENES[gene]
    rows = []
    for ln in lines:
        f = ln.split("\t")
        if len(f) != 18:
            continue
        # variant,r2,pvalue,molecular_trait_object_id,molecular_trait_id,maf,gene_id,median_tpm,beta,se,an,ac,chromosome,position,ref,alt,type,rsid
        mtid = f[4].split(".")[0]
        if mtid != g["ens"]:
            continue
        try:
            pos = int(f[13]); p = float(f[2]); beta = float(f[8]); se = float(f[9])
        except ValueError:
            continue
        if not (g["beg"] <= pos <= g["end"]):
            continue
        rows.append(dict(rsid=f[17], pos=pos, ref=f[14], alt=f[15], beta=beta, se=se,
                         p=p, maf=f[5], an=f[10]))
    return rows


def main():
    import csv
    summary = []
    for label, url, kind, genes in DATASETS:
        for gene in genes:
            out = os.path.join(OUT, f"eqtl_{label}_{gene}.csv")
            if os.path.exists(out):
                print(f"[skip] {out}"); continue
            g = GENES[gene]
            print(f"[fetch] {label} {gene} chr{g['chrom']}:{g['beg']}-{g['end']}", flush=True)
            try:
                nf = 19 if kind == "eqtlc" else 18
                lines = tabix_query(url, g["chrom"], g["beg"], g["end"], nfields=nf, threads=1)
                rows = parse_eqtlc(lines, gene) if kind == "eqtlc" else parse_gtex(lines, gene)
            except Exception as e:
                print(f"  ERROR {type(e).__name__}: {str(e)[:200]}", flush=True)
                summary.append((label, gene, -1, None))
                time.sleep(3)
                continue
            if rows:
                with open(out, "w", newline="") as fh:
                    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
                    w.writeheader(); w.writerows(rows)
                minp = min(r["p"] for r in rows)
                print(f"  rows={len(rows)} minP={minp:.2e}", flush=True)
                summary.append((label, gene, len(rows), minp))
            else:
                print("  rows=0", flush=True)
                summary.append((label, gene, 0, None))
            time.sleep(2)  # 防 DoS 限速
    print("\n=== summary ===")
    for s in summary:
        print(s)


if __name__ == "__main__":
    main()
