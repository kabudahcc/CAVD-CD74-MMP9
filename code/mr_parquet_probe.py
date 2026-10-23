# -*- coding: utf-8 -*-
"""远程探测 MetaLCL parquet schema 与行组统计量"""
import fsspec
import pyarrow.parquet as pq

URL_CD74 = ("https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/MetaLCL/"
            "Liang_2013_eQTLGen_Huang_2018_LCL_twinsuk_ensembl_gene_id/"
            "Liang_2013_WholeBlood_eQTLGen_Huang_2018_LCL_twinsuk."
            "qtl_nominal_ENSG00000019582_id.parquet")

# 先确认真实 URL（目录列表）
import requests
base = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/MetaLCL/"
r = requests.get(base, timeout=60)
print(r.status_code)
import re
links = re.findall(r'href="([^"]+)"', r.text)
for l in links[:20]:
    print(l)
print('total links:', len(links))
