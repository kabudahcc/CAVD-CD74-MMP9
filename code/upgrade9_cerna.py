# -*- coding: utf-8 -*-
"""选题B 升级9c-2：ceRNA 网络组装
miRNA-mRNA 层来自 multiMiR（mirna_shared.csv：miRTarBase 验证/双库预测的共享 miRNA）
lncRNA-miRNA 层来自 ENCORI/starBase CLIP 数据（clipExpNum>=1，每 miRNA 取 top 5）
输出：TableS_ceRNA.csv + results_upgrade9_cerna.json
"""
import os, io, json, time
import pandas as pd
import requests

WS = r"C:\Users\98757\Desktop\生信分析"
OUT = os.path.join(WS, "论文_选题B")
API = "https://rnasysu.com/encori/api/miRNATarget/"

shared = pd.read_csv(os.path.join(WS, "data", "mirna_shared.csv"))
mirnas = shared["miRNA"].tolist()
print("共享 miRNA:", len(mirnas))

def encori_lnc(mi, retries=3):
    for i in range(retries):
        try:
            r = requests.get(API, params={"assembly": "hg38", "geneType": "lncRNA", "miRNA": mi}, timeout=90)
            r.raise_for_status()
            lines = [l for l in r.text.splitlines() if not l.startswith("#")]
            if len(lines) < 2:
                return pd.DataFrame()
            return pd.read_csv(io.StringIO("\n".join(lines)), sep="\t")
        except Exception as e:
            print("  retry", i + 1, mi, e); time.sleep(5)
    return pd.DataFrame()

lnc_rows = []
for mi in mirnas:
    df = encori_lnc(mi)
    if df.empty:
        print(mi, "lncRNA 查询为空"); continue
    df = df[df["clipExpNum"] >= 1]
    top = df.groupby("geneName")["clipExpNum"].max().sort_values(ascending=False).head(5)
    for g, c in top.items():
        lnc_rows.append({"lncRNA": g, "miRNA": mi, "clipExpNum": int(c)})
    print(mi, "CLIP+ lncRNA:", df['geneName'].nunique(), "取", len(top))
    time.sleep(1)

lnc = pd.DataFrame(lnc_rows)
lnc = lnc.drop_duplicates(["lncRNA", "miRNA"])
print("\nlncRNA-miRNA 边:", len(lnc), "| 独立 lncRNA:", lnc["lncRNA"].nunique())

rows = []
for _, e in lnc.iterrows():
    for gene in ("CD74", "MMP9"):
        rows.append({"lncRNA": e["lncRNA"], "miRNA": e["miRNA"], "mRNA": gene,
                     "lnc_clipExpNum": e["clipExpNum"]})
ce = pd.DataFrame(rows)
ce.to_csv(os.path.join(OUT, "TableS_ceRNA.csv"), index=False, encoding="utf-8-sig")

json.dump({"shared_miRNAs": shared.to_dict("records"),
           "lncrna_mirna": lnc.to_dict("records"),
           "note": "miRNA-mRNA edges: miRTarBase validated or TargetScan+miRDB co-predicted (multiMiR); "
                   "lncRNA-miRNA edges: ENCORI CLIP-supported (clipExpNum>=1), top5 per miRNA"},
          open(os.path.join(OUT, "results_upgrade9_cerna.json"), "w", encoding="utf-8"),
          indent=1, ensure_ascii=False)
print("saved TableS_ceRNA.csv + results_upgrade9_cerna.json")
print("网络规模:", lnc["lncRNA"].nunique(), "lncRNA -", len(mirnas), "miRNA - 2 mRNA")
