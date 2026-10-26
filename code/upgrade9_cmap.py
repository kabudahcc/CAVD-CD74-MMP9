# -*- coding: utf-8 -*-
"""选题B 升级9a：类 cMAP 药物预测 —— SigCom LINCS（LINCS L1000 化学扰动签名）
流程（官方 issue #70 示例）：基因 symbol → 实体 UUID → enrich/ranktwosided → 签名元数据
z-sum > 0 = mimicker（与疾病同向）；z-sum < 0 = reverser（逆转签名 = 候选药）
输出：TableS_cMAP_LINCS.csv + results_upgrade9.json + 原始响应 json
"""
import os, sys, json
import pandas as pd
import requests

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pipe_common import de, load_gse83453

OUT = r"C:\Users\98757\Desktop\生信分析\论文_选题B"
META = "https://maayanlab.cloud/sigcom-lincs/metadata-api/"
DATA_API = "https://maayanlab.cloud/sigcom-lincs/data-api/api/v1/"

G1, gr1 = load_gse83453()
r1 = de(G1, gr1, "Calcified", "Normal")
sig = r1[(r1["P.Value"] < 0.05) & (r1["log2FC"].abs() > 0.585)]
up = sig[sig["log2FC"] > 0].index.tolist()
dn = sig[sig["log2FC"] < 0].index.tolist()
print("上调 DEG:", len(up), "下调 DEG:", len(dn))

# 1) symbol → 实体 UUID
r = requests.post(META + "entities/find", json={
    "filter": {"where": {"meta.symbol": {"inq": up + dn}}, "fields": ["id", "meta.symbol"]},
}, timeout=120)
r.raise_for_status()
entities = r.json()
sym2id = {e["meta"]["symbol"]: e["id"] for e in entities}
up_e = [sym2id[g] for g in up if g in sym2id]
dn_e = [sym2id[g] for g in dn if g in sym2id]
print("解析到实体: up", len(up_e), "/", len(up), "| down", len(dn_e), "/", len(dn))

# 2) ranktwosided 富集（双向：一次返回 mimickers + reversers）
query = {"up_entities": up_e, "down_entities": dn_e, "limit": 200,
         "database": "l1000_cp"}
r = requests.post(DATA_API + "enrich/ranktwosided", json=query, timeout=300)
r.raise_for_status()
dec = json.JSONDecoder()
objs, idx, text = [], 0, r.text.strip()
while idx < len(text):
    obj, end = dec.raw_decode(text, idx)
    objs.append(obj)
    idx = end
    while idx < len(text) and text[idx] in " \r\n\t":
        idx += 1
res = next((o for o in objs if "results" in o), objs[-1])
json.dump(res, open(os.path.join(OUT, "sigcom_raw_ranktwosided.json"), "w", encoding="utf-8"), indent=1)
results = res["results"]
print("返回签名数:", len(results), "| 字段:", list(results[0].keys()) if results else None)

# 3) 签名 UUID → 元数据（perturbagen 名称、细胞系、浓度）
sigids = {x["uuid"]: x for x in results}
r = requests.post(META + "signatures/find", json={
    "filter": {"where": {"id": {"inq": list(sigids.keys())}}}}, timeout=180)
r.raise_for_status()
sigs = r.json()
print("元数据条数:", len(sigs))

rows = []
for s in sigs:
    uid = s["id"]
    sc = {k: v for k, v in sigids[uid].items() if k != "uuid"}
    m = s.get("meta", {})
    zsum = sc.get("z-sum", 0) or 0
    rows.append(dict(drug=m.get("pert_name", ""), cell=m.get("cell_line", ""),
                     dose=m.get("pert_dose", ""), time=m.get("pert_time", ""),
                     zsum=zsum,
                     direction="reverser (candidate)" if zsum < 0 else "mimicker (disease-like)",
                     scores=json.dumps(sc)))
df = pd.DataFrame(rows).sort_values("zsum")
df.to_csv(os.path.join(OUT, "TableS_cMAP_LINCS.csv"), index=False, encoding="utf-8-sig")

rev = df[df["zsum"] < 0]
mim = df[df["zsum"] > 0]
print("\nTop reversers（逆转候选药）:")
print(rev.head(10)[["drug", "cell", "dose", "zsum"]].to_string(index=False))
print("\nTop mimickers:")
print(mim.tail(5)[["drug", "cell", "zsum"]].to_string(index=False))

json.dump({"n_up": len(up), "n_down": len(dn), "n_entities_up": len(up_e), "n_entities_down": len(dn_e),
           "n_signatures_returned": len(results),
           "top_reversers": rev.head(10).drop(columns="scores").to_dict("records"),
           "top_mimickers": mim.tail(10).drop(columns="scores").to_dict("records")},
          open(os.path.join(OUT, "results_upgrade9.json"), "w", encoding="utf-8"), indent=1, default=str)
print("\nsaved TableS_cMAP_LINCS.csv + results_upgrade9.json")
