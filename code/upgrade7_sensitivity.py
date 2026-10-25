# -*- coding: utf-8 -*-
"""选题B 升级7：计算敏感性分析（回应审稿意见，无需湿实验）
1) 固定标准化 + 固定阈值：训练集 scaler/Youden 阈值直接用于外部队列（不再队列内重标准化）
2) Bootstrap 稳定性选择：200 次重抽样，报告 28 候选基因被三算法选中的频率
3) FDR 口径候选筛选敏感性：discovery 用 FDR<0.05 & |log2FC|>0.585 重筛
4) 免疫相关性的组内 Spearman 与疾病状态校正后的偏相关（rank 残差法）
"""
import os, sys, json
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegressionCV, LogisticRegression
from sklearn.svm import SVC
from sklearn.feature_selection import RFE
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pipe_common import bh_fdr, de, load_gse83453, load_gse153555, immune_scores, read_series_matrix, collapse_probes, DATA

OUT = r"C:\Users\98757\Desktop\生信分析\论文_选题B"
HUB = ["CD74", "MMP9"]
RES = {}

G1, gr1 = load_gse83453()
G2, gr2 = load_gse153555()
ytr = (gr1 == "Calcified").astype(int).values
trX = G1.T[HUB].values

# GSE12644（同 upgrade4 口径）
expr, titles = read_series_matrix(os.path.join(DATA, "GSE12644_series_matrix.txt.gz"))
gr3 = pd.Series(["Case" if "Case" in t else "Ctrl" for t in titles], index=expr.columns)
p2s = {k: v for k, v in json.load(open(os.path.join(DATA, "gpl570_probe2symbol.json"))).items() if v}
G3 = collapse_probes(expr, p2s)
if G3.values.max() > 100:
    G3 = np.log2(G3 + 1)
yte3 = (gr3 == "Case").astype(int).values

# ---------- 1) 固定标准化 + 固定阈值 ----------
sc = StandardScaler().fit(trX)
clf = LogisticRegression(penalty="l2", C=1.0, max_iter=5000).fit(sc.transform(trX), ytr)
cv_split = StratifiedKFold(5, shuffle=True, random_state=42)
p_cv = cross_val_predict(LogisticRegression(penalty="l2", C=1.0, max_iter=5000),
                         sc.transform(trX), ytr, cv=cv_split, method="predict_proba")[:, 1]
# 训练集 CV 预测上的 Youden 阈值
ths = np.linspace(p_cv.min(), p_cv.max(), 501)
youden_t = ths[np.argmax([((p_cv >= t) & (ytr == 1)).sum() / (ytr == 1).sum() +
                          ((p_cv < t) & (ytr == 0)).sum() / (ytr == 0).sum() - 1 for t in ths])]

def fixed_eval(Gx, yte, label):
    teX = Gx.loc[HUB].T.values if isinstance(Gx, pd.DataFrame) else None
    # 固定标准化：用训练集均值/方差直接变换外部原始值
    z = (teX - sc.mean_) / sc.scale_
    p = 1 / (1 + np.exp(-(clf.intercept_[0] + z @ clf.coef_[0])))
    auc = roc_auc_score(yte, p)
    pred = (p >= youden_t).astype(int)
    se = ((pred == 1) & (yte == 1)).sum() / max(1, (yte == 1).sum())
    sp = ((pred == 0) & (yte == 0)).sum() / max(1, (yte == 0).sum())
    RES[f"fixed_{label}"] = dict(auc=float(auc), sens=float(se), spec=float(sp))
    print(f"固定标准化+固定阈值 {label}: AUC={auc:.3f} sens={se:.3f} spec={sp:.3f}")

fixed_eval(G2, (gr2 == "AS").astype(int).values, "GSE153555")
fixed_eval(G3, yte3, "GSE12644")
RES["youden_threshold_train_cv"] = float(youden_t)

# ---------- 2) Bootstrap 稳定性选择 ----------
lists = json.load(open(os.path.join(DATA, "gene_lists.json"), encoding="utf-8"))
IMMUNE = set(lists["immune_kegg"]); SECRETED = set(lists["secreted_uniprot"])
r1 = de(G1, gr1, "Calcified", "Normal")
r2 = de(G2, gr2, "AS", "NC")
deg1 = set(r1[(r1["P.Value"] < 0.05) & (r1["log2FC"].abs() > 0.585)].index)
cand = deg1 & IMMUNE & SECRETED & set(r2.index)
both = r1.loc[list(cand), ["log2FC", "P.Value"]].join(r2.loc[list(cand), ["log2FC", "P.Value"]], lsuffix="_1", rsuffix="_2")
cand2 = both[(np.sign(both["log2FC_1"]) == np.sign(both["log2FC_2"])) & (both["P.Value_2"] < 0.05)]
genes = sorted(cand2.index.tolist())
X = G1.loc[genes].T.values

B = 200
rng = np.random.RandomState(42)
cnt = {a: pd.Series(0, index=genes) for a in ("lasso", "svm", "rf")}
cnt["consensus"] = pd.Series(0, index=genes)
idx_pos = np.where(ytr == 1)[0]; idx_neg = np.where(ytr == 0)[0]
for b in range(B):
    idx = np.concatenate([rng.choice(idx_pos, len(idx_pos), replace=True),
                          rng.choice(idx_neg, len(idx_neg), replace=True)])
    Xb = StandardScaler().fit_transform(X[idx]); yb = ytr[idx]
    lr = LogisticRegressionCV(Cs=20, penalty="l1", solver="liblinear", cv=StratifiedKFold(5),
                              max_iter=5000, random_state=42).fit(Xb, yb)
    la = set(g for g, c in zip(genes, lr.coef_[0]) if abs(c) > 1e-6)
    sv = set(g for g, s in zip(genes, RFE(SVC(kernel="linear", random_state=42),
               n_features_to_select=5).fit(Xb, yb).support_) if s)
    rf = RandomForestClassifier(n_estimators=500, random_state=42).fit(Xb, yb)
    top5 = set(pd.Series(rf.feature_importances_, index=genes).nlargest(5).index)
    for g in la: cnt["lasso"][g] += 1
    for g in sv: cnt["svm"][g] += 1
    for g in top5: cnt["rf"][g] += 1
    for g in la & sv & top5: cnt["consensus"][g] += 1
stab = pd.DataFrame({k: (v / B).round(3) for k, v in cnt.items()}).sort_values("consensus", ascending=False)
stab.to_csv(os.path.join(OUT, "TableS_stability_selection.csv"), encoding="utf-8-sig")
RES["stability"] = {g: stab.loc[g].to_dict() for g in HUB}
print("稳定性选择 CD74:", stab.loc["CD74"].to_dict(), "MMP9:", stab.loc["MMP9"].to_dict())

# ---------- 3) FDR 口径候选筛选敏感性 ----------
deg1_fdr = set(r1[(r1["FDR"] < 0.05) & (r1["log2FC"].abs() > 0.585)].index)
cand_fdr = deg1_fdr & IMMUNE & SECRETED & set(r2.index)
both_f = r1.loc[list(cand_fdr), ["log2FC", "P.Value"]].join(r2.loc[list(cand_fdr), ["log2FC", "P.Value"]], lsuffix="_1", rsuffix="_2")
cand2_fdr = both_f[(np.sign(both_f["log2FC_1"]) == np.sign(both_f["log2FC_2"])) & (both_f["P.Value_2"] < 0.05)]
RES["fdr_candidates"] = dict(n_deg_fdr=len(deg1_fdr), n_cand=len(cand_fdr), n_rep=len(cand2_fdr),
                             genes=sorted(cand2_fdr.index.tolist()),
                             hub_retained=all(g in cand2_fdr.index for g in HUB),
                             overlap_with_28=len(set(cand2_fdr.index) & set(genes)))
print("FDR 口径: 候选", len(cand_fdr), "复制", len(cand2_fdr), "| hub 保留:", RES["fdr_candidates"]["hub_retained"])

# ---------- 4) 组内相关 + 偏相关 ----------
sc1 = immune_scores(G1)
m_calc = gr1 == "Calcified"
rows = []
for g in HUB:
    ge = G1.loc[g].values
    for cell in sc1.columns:
        s = sc1[cell].values
        r_all = stats.spearmanr(ge, s).statistic
        r_in = stats.spearmanr(ge[m_calc.values], s[m_calc.values]).statistic  # 钙化组内 n=19
        # 偏相关（rank 残差法，校正疾病状态）
        rx = pd.Series(stats.rankdata(ge)); ry = pd.Series(stats.rankdata(s))
        dz = pd.Series(ytr.astype(float))
        sx = np.linalg.lstsq(np.c_[np.ones(len(dz)), dz], rx, rcond=None)[0]
        sy = np.linalg.lstsq(np.c_[np.ones(len(dz)), dz], ry, rcond=None)[0]
        resx = rx - (sx[0] + sx[1] * dz); resy = ry - (sy[0] + sy[1] * dz)
        r_par = stats.spearmanr(resx, resy).statistic
        p_par = stats.spearmanr(resx, resy).pvalue
        rows.append(dict(gene=g, cell=cell, r_pooled=round(float(r_all), 3),
                         r_within_calcified=round(float(r_in), 3),
                         r_partial_adj_state=round(float(r_par), 3), P_partial=float(p_par)))
pc = pd.DataFrame(rows)
pc["FDR_partial"] = bh_fdr(pc["P_partial"].values)
pc.round(4).to_csv(os.path.join(OUT, "TableS_partial_corr.csv"), index=False, encoding="utf-8-sig")
RES["partial_corr_top"] = pc.sort_values("P_partial").head(6).to_dict("records")
RES["partial_corr_n_sig"] = int((pc["FDR_partial"] < 0.05).sum())
print("偏相关 FDR<0.05 的 基因×细胞 对数:", RES["partial_corr_n_sig"], "/ 26")

json.dump(RES, open(os.path.join(OUT, "results_upgrade7.json"), "w", encoding="utf-8"), indent=1, ensure_ascii=False, default=str)
print("saved results_upgrade7.json + TableS_stability_selection.csv + TableS_partial_corr.csv")
