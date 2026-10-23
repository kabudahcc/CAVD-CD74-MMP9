# -*- coding: utf-8 -*-
"""MetaLCL 远程 parquet 行组裁剪：只取基因 cis 窗口内的变异"""
import sys
import fsspec
import pyarrow.parquet as pq
import pandas as pd

BASE = ("https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/MetaLCL/"
        "Liang_2013_ALSPAC_CAP_eur_GENCORD_GEUVADIS_GTEx_TwinsUK_CoLausR2_"
        "naive_LCL_metaanalysis_{ens}_id.parquet")

GENES = {
    "CD74": dict(ens="ENSG00000019582", chrom="5", beg=150_300_000, end=150_520_000),
    "MMP9": dict(ens="ENSG00000100985", chrom="20", beg=45_900_000, end=46_120_000),
}


def fetch(gene):
    g = GENES[gene]
    url = BASE.format(ens=g["ens"])
    print(f"[{gene}] open {url.split('/')[-1]}", flush=True)
    with fsspec.open(url, "rb") as f:
        pf = pq.ParquetFile(f)
        md = pf.metadata
        print(f"[{gene}] rows={md.num_rows} rowgroups={md.num_row_groups}", flush=True)
        print(f"[{gene}] columns={pf.schema_arrow.names}", flush=True)
        keep = []
        for i in range(md.num_row_groups):
            rg = md.row_group(i)
            stats = {}
            for j in range(rg.num_columns):
                col = rg.column(j)
                name = col.path_in_schema
                if col.is_stats_set:
                    s = col.statistics
                    if s.has_min_max:
                        stats[name] = (s.min, s.max)
            # 找染色体与位置列
            chrom_col = next((c for c in ("chromosome", "chrom", "chr") if c in stats), None)
            pos_col = next((c for c in ("position", "pos", "base_pair_location") if c in stats), None)
            if chrom_col is None or pos_col is None:
                if i == 0:
                    print(f"[{gene}] rg0 stats keys: {list(stats.keys())}", flush=True)
                continue
            cmin, cmax = stats[chrom_col]
            pmin, pmax = stats[pos_col]
            cmin, cmax = str(cmin), str(cmax)
            # 行组染色体为单值时才可裁剪
            if cmin == cmax and cmin.lstrip("chr") != g["chrom"]:
                continue
            if pmax < g["beg"] or pmin > g["end"]:
                continue
            keep.append(i)
        print(f"[{gene}] keep rowgroups: {keep}", flush=True)
        if not keep:
            return None
        tbl = pf.read_row_groups(keep)
        df = tbl.to_pandas()
        pos_col = next(c for c in ("position", "pos", "base_pair_location") if c in df.columns)
        chrom_col = next(c for c in ("chromosome", "chrom", "chr") if c in df.columns)
        df = df[(df[chrom_col].astype(str).str.lstrip("chr") == g["chrom"])
                & (df[pos_col] >= g["beg"]) & (df[pos_col] <= g["end"])]
        out = f"data/mr/metalcl_{gene}.csv"
        df.to_csv(out, index=False)
        print(f"[{gene}] saved {len(df)} rows -> {out}", flush=True)
        print(df.head(3).to_string(), flush=True)
        return df


if __name__ == "__main__":
    for gene in (sys.argv[1:] or ["CD74", "MMP9"]):
        try:
            fetch(gene)
        except Exception as e:
            print(f"[{gene}] ERROR {type(e).__name__}: {e}", flush=True)
