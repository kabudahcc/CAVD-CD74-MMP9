# -*- coding: utf-8 -*-
"""选题B 审稿回应：补充统计量
1) train/CV/external/GSE12644 AUC 的 bootstrap 95% CI
2) Youden 阈值下的敏感性/特异性（CV 与 external）
3) 外部校准斜率与截距
4) 嵌套 5 折 CV（特征选择内嵌）作为乐观偏倚的敏感性分析
"""
import os, sys, json
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score, roc_curve, brier_score_loss
from sklearn.svm import SVC
from sklearn.feature_selection import RFE
from sklearn.ensemble import RandomForestClassifier

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pipe_common import load_gse83453, load_gse153555, read_series_matrix, collapse_probes, DATA

OUT = r"C:\Users\98757\Desktop\生信分析\论文_选题B"
HUB = ["CD74", "MMP9"]
rng = np.random.default_rng(42)

CAND = ['C6','CCL19','CCL21','CCL5','CD14','CD209','CD48','CD74','CFI','CNTFR','CXCL12',
        'CXCL14','CXCL16','F5','IL17D','IL32','ISG15','LEPR','MASP1','MMP9','PLAU','PLAUR',
        'SERPINA1','SERPINE1','SPP1','TNFRSF11B','TNFSF10','VCAM1']


def boot_auc_ci(y, p, n=2000):
    y = np.asarray(y); p = np.asarray(p)
    aucs = []
    for _ in range(n):
        idx = rng.integers(0, len(y), len(y))
        if len(np.unique(y[idx])) < 2:
            continue
        aucs.append(roc_auc_score(y[idx], p[idx]))
    return float(np.percentile(aucs, 2.5)), float(np.percentile(aucs, 97.5))


def youden_sens_spec(y, p):
    fpr, tpr, thr = roc_curve(y, p)
    j = tpr - fpr
    i = int(np.argmax(j))
    return float(tpr[i]), float(1 - fpr[i]), float(thr[i])


# ---------- 数据 ----------
G1, gr1 = load_gse83453()
ytr = (gr1 == "Calcified").astype(int).values
Xtr_all = G1.T[CAND].values
trX = G1.T[HUB].values

G2, gr2 = load_gse153555()
m = gr2.isin(["AS", "NC", "AI"])
yte = (gr2[m] == "AS").astype(int).values
teX = G2.loc[HUB, m].T.values

expr, titles = read_series_matrix(os.path.join(DATA, "GSE12644_series_matrix.txt.gz"))
gr3 = pd.Series(["Case" if "Case" in t else "Ctrl" for t in titles], index=expr.columns)
p2s = {k: v for k, v in json.load(open(os.path.join(DATA, "gpl570_probe2symbol.json"))).items() if v}
G3 = collapse_probes(expr, p2s)
if G3.values.max() > 100:
    G3 = np.log2(G3 + 1)
y3 = (gr3 == "Case").astype(int).values
X3 = G3.loc[HUB].T.values

# ---------- 主模型 ----------
sc = StandardScaler().fit(trX)
clf = LogisticRegression(penalty="l2", C=1.0, max_iter=5000).fit(sc.transform(trX), ytr)
p_train = clf.predict_proba(sc.transform(trX))[:, 1]
cv = StratifiedKFold(5, shuffle=True, random_state=42)
p_cv = cross_val_predict(LogisticRegression(penalty="l2", C=1.0, max_iter=5000),
                         sc.transform(trX), ytr, cv=cv, method="predict_proba")[:, 1]
sc2 = StandardScaler().fit(teX)
logit_ext = clf.intercept_[0] + (sc2.transform(teX) * clf.coef_[0]).sum(axis=1)
p_ext = 1 / (1 + np.exp(-logit_ext))
sc3 = StandardScaler().fit(X3)
logit_3 = clf.intercept_[0] + (sc3.transform(X3) * clf.coef_[0]).sum(axis=1)
p_3 = 1 / (1 + np.exp(-logit_3))

res = {}
for name, (yy, pp) in dict(train=(ytr, p_train), cv=(ytr, p_cv),
                           ext_GSE153555=(yte, p_ext), untouched_GSE12644=(y3, p_3)).items():
    auc = roc_auc_score(yy, pp)
    lo, hi = boot_auc_ci(yy, pp)
    se, sp, thr = youden_sens_spec(yy, pp)
    res[name] = dict(auc=float(auc), ci_lo=lo, ci_hi=hi,
                     sens=se, spec=sp, youden_thr=thr)
    print(f"{name}: AUC={auc:.3f} (95% CI {lo:.3f}-{hi:.3f}) sens={se:.3f} spec={sp:.3f}")

# ---------- 外部校准斜率/截距 ----------
def cal_slope_intercept(y, logit):
    lr = LogisticRegression(penalty=None, max_iter=5000).fit(logit.reshape(-1, 1), y)
    return float(lr.intercept_[0]), float(lr.coef_[0][0])

logit_cv = np.log(np.clip(p_cv, 1e-9, 1 - 1e-9) / (1 - np.clip(p_cv, 1e-9, 1 - 1e-9)))
ci_cv = cal_slope_intercept(ytr, logit_cv)
ci_ext = cal_slope_intercept(yte, logit_ext)
res["calibration"] = dict(cv_intercept=ci_cv[0], cv_slope=ci_cv[1],
                          ext_intercept=ci_ext[0], ext_slope=ci_ext[1])
print(f"calibration CV: intercept={ci_cv[0]:.3f} slope={ci_cv[1]:.3f}; "
      f"external: intercept={ci_ext[0]:.3f} slope={ci_ext[1]:.3f}")

# ---------- 嵌套 CV（特征选择内嵌于每一折） ----------
def select_genes(X, y, genes):
    Xz = StandardScaler().fit_transform(X)
    l1 = LogisticRegressionCV(penalty="l1", solver="liblinear", Cs=20, cv=5,
                              max_iter=5000, random_state=0).fit(Xz, y)
    lasso = set(np.array(genes)[np.abs(l1.coef_[0]) > 1e-8])
    rfe = RFE(SVC(kernel="linear", C=1.0), n_features_to_select=5).fit(Xz, y)
    svm = set(np.array(genes)[rfe.support_])
    rf = RandomForestClassifier(n_estimators=500, random_state=0).fit(Xz, y)
    top5 = set(np.array(genes)[np.argsort(rf.feature_importances_)[::-1][:5]])
    hub = lasso & svm & top5
    if not hub:
        votes = {g: (g in lasso) + (g in svm) + (g in top5) for g in genes}
        hub = {g for g, v in votes.items() if v >= 2}
    if not hub:
        hub = set(svm)
    return sorted(hub)

outer = StratifiedKFold(5, shuffle=True, random_state=7)
p_nested = np.zeros(len(ytr))
pick_log = []
for tr, te in outer.split(Xtr_all, ytr):
    sel = select_genes(Xtr_all[tr], ytr[tr], CAND)
    pick_log.append(sel)
    idx = [CAND.index(g) for g in sel]
    s = StandardScaler().fit(Xtr_all[tr][:, idx])
    m2 = LogisticRegression(penalty="l2", C=1.0, max_iter=5000).fit(s.transform(Xtr_all[tr][:, idx]), ytr[tr])
    p_nested[te] = m2.predict_proba(s.transform(Xtr_all[te][:, idx]))[:, 1]
    print("outer fold selected:", sel)

auc_n = roc_auc_score(ytr, p_nested)
lo, hi = boot_auc_ci(ytr, p_nested)
res["nested_cv"] = dict(auc=float(auc_n), ci_lo=lo, ci_hi=hi,
                        folds=pick_log,
                        cd74_picked=sum("CD74" in s for s in pick_log),
                        mmp9_picked=sum("MMP9" in s for s in pick_log))
print(f"nested CV AUC={auc_n:.3f} (95% CI {lo:.3f}-{hi:.3f}); "
      f"CD74 picked {res['nested_cv']['cd74_picked']}/5, MMP9 picked {res['nested_cv']['mmp9_picked']}/5")

json.dump(res, open(os.path.join(OUT, "results_upgrade6.json"), "w"), indent=1)
print("saved results_upgrade6.json")
