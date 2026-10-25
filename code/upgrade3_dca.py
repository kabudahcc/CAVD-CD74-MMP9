# -*- coding: utf-8 -*-
"""选题B 升级3：校准曲线 + DCA + 列线图（CD74+MMP9 双基因模型）
训练 GSE83453，外部验证 GSE153555（AS vs NC）
"""
import os, sys, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score, brier_score_loss

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pipe_common import load_gse83453, load_gse153555

plt.rcParams["font.family"] = "Arial"
OUT = r"C:\Users\98757\Desktop\生信分析\论文_选题B"
HUB = ["CD74", "MMP9"]

G1, gr1 = load_gse83453()
G2, gr2 = load_gse153555()
ytr = (gr1 == "Calcified").astype(int).values
trX = G1.T[HUB].values
te_mask = gr2.isin(["AS", "NC", "AI"])
G2s = G2.loc[:, te_mask]
gr2s = gr2[te_mask]
yte = (gr2s == "AS").astype(int).values  # 与手稿口径一致：AS vs NC+AI
teX = G2s.loc[HUB].T.values

sc = StandardScaler().fit(trX)
clf = LogisticRegression(penalty="l2", C=1.0, max_iter=5000).fit(sc.transform(trX), ytr)
print("coef:", dict(zip(HUB, clf.coef_[0])), "intercept:", clf.intercept_)

p_train = clf.predict_proba(sc.transform(trX))[:, 1]
cv = StratifiedKFold(5, shuffle=True, random_state=42)
p_cv = cross_val_predict(LogisticRegression(penalty="l2", C=1.0, max_iter=5000),
                         sc.transform(trX), ytr, cv=cv, method="predict_proba")[:, 1]
# 外部验证：批次差异下按外部队列内部 z-score 标准化（与原管线一致），系数沿用训练集
sc2 = StandardScaler().fit(teX)
logit_ext = clf.intercept_[0] + (sc2.transform(teX) * clf.coef_[0]).sum(axis=1)
p_ext = 1 / (1 + np.exp(-logit_ext))
print("AUC train/cv/ext:",
      roc_auc_score(ytr, p_train), roc_auc_score(ytr, p_cv), roc_auc_score(yte, p_ext))
print("Brier cv/ext:", brier_score_loss(ytr, p_cv), brier_score_loss(yte, p_ext))


def calibration_bins(y, p, n_bins=8):
    df = pd.DataFrame({"y": y, "p": p}).sort_values("p")
    df["bin"] = pd.qcut(df["p"].rank(method="first"), n_bins, labels=False)
    g = df.groupby("bin").agg(pred=("p", "mean"), obs=("y", "mean"), n=("y", "size"))
    return g


# ---------- Fig8 校准曲线 ----------
fig, ax = plt.subplots(figsize=(5.6, 5.4), dpi=300)
ax.plot([0, 1], [0, 1], ls="--", c="grey", lw=1, label="Ideal")
for (yy, pp, lab, c, m) in [
        (ytr, p_train, f"GSE83453 apparent (AUC={roc_auc_score(ytr, p_train):.3f})", "#C0504D", "o"),
        (ytr, p_cv, f"GSE83453 5-fold CV (AUC={roc_auc_score(ytr, p_cv):.3f})", "#E8862E", "s"),
        (yte, p_ext, f"GSE153555 replication (AUC={roc_auc_score(yte, p_ext):.3f})", "#4F81BD", "^")]:
    g = calibration_bins(yy, pp, 8 if len(yy) > 12 else 5)
    ax.plot(g["pred"], g["obs"], marker=m, c=c, lw=1.8, ms=6, label=lab)
ax.set_xlabel("Predicted probability of CAVD")
ax.set_ylabel("Observed proportion")
ax.set_title("Calibration of the CD74+MMP9 diagnostic model")
ax.legend(fontsize=8, loc="upper left")
ax.text(0.55, 0.12, f"Brier (CV) = {brier_score_loss(ytr, p_cv):.3f}\nBrier (replication) = {brier_score_loss(yte, p_ext):.3f}",
        fontsize=8.5, transform=ax.transAxes)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "Fig8_calibration.png"), bbox_inches="tight")
plt.close(fig)

# ---------- Fig9 DCA ----------
def net_benefit(y, p, thresholds):
    n = len(y)
    out = []
    for t in thresholds:
        pred = p >= t
        tp = ((pred == 1) & (y == 1)).sum()
        fp = ((pred == 1) & (y == 0)).sum()
        out.append(tp / n - fp / n * (t / (1 - t)))
    return np.array(out)

ths = np.linspace(0.01, 0.80, 79)
fig, ax = plt.subplots(figsize=(7.5, 5.4), dpi=300)
prev_tr = ytr.mean(); prev_te = yte.mean()
ax.plot(ths, net_benefit(ytr, p_cv, ths), c="#E8862E", lw=2,
        label=f"Model, GSE83453 (5-fold CV)")
ax.plot(ths, net_benefit(yte, p_ext, ths), c="#4F81BD", lw=2,
        label=f"Model, GSE153555 (replication)")
ax.plot(ths, prev_tr - (1 - prev_tr) * ths / (1 - ths), c="grey", ls=":", lw=1.4, label="Treat all (train prevalence)")
ax.plot(ths, prev_te - (1 - prev_te) * ths / (1 - ths), c="lightgrey", ls=":", lw=1.4, label="Treat all (replication prevalence)")
ax.axhline(0, c="black", lw=1.4, label="Treat none")
ax.set_ylim(-0.1, 0.8)
ax.set_xlim(0, 0.8)
ax.set_xlabel("Threshold probability")
ax.set_ylabel("Net benefit")
ax.set_title("Decision curve analysis: CD74+MMP9 model")
ax.legend(fontsize=8, bbox_to_anchor=(1.02, 1), loc="upper left")
plt.tight_layout()
plt.savefig(os.path.join(OUT, "Fig9_DCA.png"), bbox_inches="tight")
plt.close(fig)

# ---------- Fig10 列线图（静态绘制） ----------
# 以标准化表达值计分：每个基因的系数映射到 0-100 分
Zs = sc.transform(trX)
zmin, zmax = Zs.min(axis=0), Zs.max(axis=0)
coef = clf.coef_[0]
# 每基因在其观测范围内贡献的分值范围
ranges = coef * (zmax - zmin)
scale = 100 / np.abs(ranges).sum()
pts_per_gene = coef * scale  # 每单位 z 的分值

fig, ax = plt.subplots(figsize=(8.4, 4.8), dpi=300)
ax.axis("off")
y_rows = {"CD74": 0.86, "MMP9": 0.66, "Total points": 0.36, "Predicted probability": 0.14}

def draw_axis(y, vmin, vmax, to_pts, label, n_ticks=6, labels_below=True):
    ax.text(0.02, y, label, fontsize=10, fontweight="bold", va="center")
    ticks = np.linspace(vmin, vmax, n_ticks)
    pts = to_pts(ticks)
    for t, pt in zip(ticks, pts):
        x = 0.16 + 0.74 * (pt - 0) / 100.0
        ax.plot([x, x], [y - 0.012, y + 0.012], c="black", lw=0.8)
        tl = f"{t:.1f}".replace("-0.0", "0.0")
        if labels_below:
            ax.text(x, y - 0.045, tl, ha="center", fontsize=7.5)
        else:
            ax.text(x, y + 0.022, tl, ha="center", fontsize=7.5)
    ax.plot([0.16, 0.90], [y, y], c="black", lw=1)

# 分值轴：各基因 z 值 -> 分值（相对其最小值平移使最小为 0）
for i, g in enumerate(HUB):
    vmin, vmax = zmin[i], zmax[i]
    to_pts = lambda z, i=i: (z - zmin[i]) * pts_per_gene[i]
    draw_axis(y_rows[g], vmin, vmax, to_pts, f"{g} (z-score)")

# Points 顶部参考轴
y = 0.97
ax.text(0.02, y, "Points", fontsize=10, fontweight="bold", va="center")
for pt in range(0, 101, 10):
    x = 0.16 + 0.74 * pt / 100.0
    ax.plot([x, x], [y - 0.012, y + 0.012], c="black", lw=0.8)
    ax.text(x, y + 0.022, str(pt), ha="center", fontsize=7.5)
ax.plot([0.16, 0.90], [y, y], c="black", lw=1)

# Total points -> probability
# total = sum_i (z_i - zmin_i)*pts_i ; logit = intercept + sum coef*z = intercept + sum coef*(zmin + pts/pts_per)
inter = clf.intercept_[0]
base = inter + (coef * zmin).sum()
def prob_of_total(tp):
    logit = base + tp / scale
    return 1 / (1 + np.exp(-logit))
totals = np.linspace(0, 100, 6)
y = y_rows["Total points"]
ax.text(0.02, y, "Total points", fontsize=10, fontweight="bold", va="center")
for tp in totals:
    x = 0.16 + 0.74 * tp / 100.0
    ax.plot([x, x], [y - 0.012, y + 0.012], c="black", lw=0.8)
    ax.text(x, y + 0.022, f"{tp:.0f}", ha="center", fontsize=7.5)
ax.plot([0.16, 0.90], [y, y], c="black", lw=1)

y = y_rows["Predicted probability"]
ax.text(0.02, y, "P(CAVD)", fontsize=10, fontweight="bold", va="center")
probs = [0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95]
for pr in probs:
    logit = np.log(pr / (1 - pr))
    tp = (logit - base) * scale
    if 0 <= tp <= 100:
        x = 0.16 + 0.74 * tp / 100.0
        ax.plot([x, x], [y - 0.012, y + 0.012], c="#C00000", lw=0.9)
        ax.text(x, y - 0.045, f"{pr:.2f}", ha="center", fontsize=7.5, color="#C00000")
ax.plot([0.16, 0.90], [y, y], c="black", lw=1)
ax.set_xlim(0, 1); ax.set_ylim(0, 1.05)
ax.set_title("Nomogram: CD74 + MMP9 diagnostic model for CAVD", fontsize=11)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "Fig10_nomogram.png"), bbox_inches="tight")
plt.close(fig)

res = dict(auc_train=float(roc_auc_score(ytr, p_train)),
           auc_cv=float(roc_auc_score(ytr, p_cv)),
           auc_ext=float(roc_auc_score(yte, p_ext)),
           brier_cv=float(brier_score_loss(ytr, p_cv)),
           brier_ext=float(brier_score_loss(yte, p_ext)),
           coef={g: float(c) for g, c in zip(HUB, clf.coef_[0])},
           intercept=float(clf.intercept_[0]))
json.dump(res, open(os.path.join(OUT, "results_upgrade3.json"), "w"), indent=1)
print(json.dumps(res, indent=1))
print("saved Fig8/Fig9/Fig10")
