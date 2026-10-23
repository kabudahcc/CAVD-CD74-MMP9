# -*- coding: utf-8 -*-
"""选题B 升级1：CD74/MMP9 -> CAVS(FinnGen R12) 两样本 MR
暴露：eQTL Catalogue r7 / GTEx v8 cis-eQTL 窗口
结局：FinnGen R12 I9_CAVS_OPERATED
流程：工具变量筛选(距离剪枝) -> rsid 谐和(等位基因对齐, 回文剔除) -> Wald ratio / IVW
"""
import os, csv, math, json
import pandas as pd
import numpy as np
from scipy import stats

MR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "mr")

GENES = {
    "CD74": dict(chrom="5", center=150_407_000),
    "MMP9": dict(chrom="20", center=46_012_000),
}

# ---------- 载入 ----------
def load_exposure(path, source):
    df = pd.read_csv(path)
    df["source"] = source
    for c in ("pos", "beta", "se", "p"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["maf"] = pd.to_numeric(df.get("maf"), errors="coerce")
    df = df.dropna(subset=["pos", "beta", "se", "p"])
    df = df[df["rsid"].astype(str).str.startswith("rs")]
    return df


def load_gtex_cache(path, source):
    rows = []
    for ln in open(path):
        f = ln.rstrip("\n").split("\t")
        if len(f) != 18:
            continue
        try:
            rows.append(dict(rsid=f[17], pos=int(f[13]), ref=f[14], alt=f[15],
                             beta=float(f[8]), se=float(f[9]), p=float(f[2]),
                             maf=float(f[5]) if f[5] not in ("NA", "") else np.nan,
                             source=source))
        except ValueError:
            continue
    return pd.DataFrame(rows)


def clump_by_distance(df, thresh, win=250_000):
    """简化 LD 剪枝：按 p 升序，保留与已选 SNP 距离 > win 的"""
    df = df[df["p"] < thresh].sort_values("p")
    picked = []
    for _, r in df.iterrows():
        if all(abs(r["pos"] - p["pos"]) > win for p in picked):
            picked.append(r)
    return pd.DataFrame(picked) if picked else pd.DataFrame()


def f_stat(beta, se):
    return (beta / se) ** 2


# ---------- 谐和 ----------
def harmonize(exp, out):
    out = out.rename(columns={"beta": "beta_out", "se": "se_out", "p": "p_out",
                              "ref": "ref_out", "alt": "alt_out"})
    m = exp.merge(out, on="rsid", suffixes=("_exp", ""))
    keep = []
    for _, r in m.iterrows():
        e_ref, e_alt = r["ref"], r["alt"]
        o_ref, o_alt = r["ref_out"], r["alt_out"]
        alleles = {e_ref, e_alt}
        if {o_ref, o_alt} != alleles or len(alleles) != 2:
            continue  # 等位基因不匹配（含indel差异）
        flip = (o_ref == e_alt)  # 结局效应等位基因是 alt_out；若 alt_out==e_ref 则需翻转
        # 结局 beta 是对 alt_out 的；暴露 beta 是对 alt(exp) 的
        if o_alt == e_alt:
            b_out = r["beta_out"]
            af_out = r["af_alt"]
        elif o_alt == e_ref:
            b_out = -r["beta_out"]
            af_out = 1 - r["af_alt"]
        else:
            continue
        # 回文 SNP 检查
        pal = {tuple(sorted((e_ref, e_alt)))} & {("A", "T"), ("C", "G")}
        maf_e = r.get("maf_exp", np.nan)
        if pal:
            mafs = [x for x in (maf_e, min(af_out, 1 - af_out)) if pd.notna(x)]
            if mafs and all(0.4 < x < 0.6 for x in mafs):
                continue  # 模糊回文剔除
        keep.append(dict(rsid=r["rsid"], pos=r["pos"], beta_exp=r["beta"], se_exp=r["se"],
                         p_exp=r["p"], beta_out=b_out, se_out=r["se_out"],
                         p_out=r["p_out"], ea=e_alt, source=r["source"]))
    return pd.DataFrame(keep)


# ---------- MR ----------
def mr_estimate(h):
    n = len(h)
    if n == 0:
        return None
    wald = h["beta_out"] / h["beta_exp"]
    se_wald = h["se_out"] / h["beta_exp"].abs()
    if n == 1:
        b, s = wald.iloc[0], se_wald.iloc[0]
        method = "Wald ratio"
        p = 2 * stats.norm.sf(abs(b / s))
        q = np.nan
    else:
        w = 1 / se_wald**2
        b = float((w * wald).sum() / w.sum())
        # 乘性随机效应：残差乘数
        q = float((w * (wald - b) ** 2).sum())
        phi = max(1.0, q / (n - 1))
        s = math.sqrt(phi / w.sum())
        p = 2 * stats.norm.sf(abs(b / s))
        method = "IVW (multiplicative random effects)"
    return dict(n_snp=n, method=method, beta=b, se=s, p=p,
                OR=math.exp(b), OR_lo=math.exp(b - 1.96 * s), OR_hi=math.exp(b + 1.96 * s),
                Q=q if n > 1 else np.nan)


def main():
    results = []
    harm_store = {}
    # 暴露数据集选择：每个基因按 minP 选主数据集 + 复制数据集
    # 仅报告主要暴露数据集（每个基因/组织组合的最强来源）
    plan = {
        "MMP9": [("Quach_monocyte", "eqtl_Quach_monocyte_MMP9.csv", 5e-8)],
        "CD74": [("OneK1K_NK", "eqtl_OneK1K_NK_CD74.csv", 1e-4),
                 ("GTEx_WholeBlood", "gtex_Whole_blood_CD74.txt", 1e-4)],
    }
    for gene, dsets in plan.items():
        # 唯一结局：FinnGen R12 I9_CAVS_OPERATED（钙化性主动脉瓣狭窄，手术表型）
        for outcome_label, outcome_fn in [("CAVS_operated", f"finngen_{gene}.csv")]:
            out = pd.read_csv(os.path.join(MR, outcome_fn))
            for source, fn, thresh in dsets:
                path = os.path.join(MR, fn)
                if not os.path.exists(path):
                    continue
                exp = load_gtex_cache(path, source) if fn.endswith(".txt") else load_exposure(path, source)
                ins = clump_by_distance(exp, thresh)
                if ins.empty:
                    continue
                ins["F"] = f_stat(ins["beta"], ins["se"])
                h = harmonize(ins, out)
                if h.empty:
                    continue
                est = mr_estimate(h)
                est.update(gene=gene, exposure=source, outcome=outcome_label, thresh=thresh,
                           n_instruments=len(ins), n_harmonised=len(h),
                           median_F=float(ins["F"].median()), min_F=float(ins["F"].min()),
                           min_p_exp=float(ins["p"].min()))
                results.append(est)
                harm_store[f"{gene}_{source}_{outcome_label}"] = h
                print(f"{gene}/{source}/{outcome_label}: IVs={len(ins)} harmonised={len(h)} "
                      f"medianF={est['median_F']:.1f} OR={est['OR']:.3f} "
                      f"({est['OR_lo']:.3f}-{est['OR_hi']:.3f}) p={est['p']:.3g}")
    res = pd.DataFrame(results)
    res.to_csv(os.path.join(MR, "Table_MR_results.csv"), index=False)
    for k, h in harm_store.items():
        h.to_csv(os.path.join(MR, f"harmonised_{k}.csv"), index=False)
    print("\nsaved Table_MR_results.csv")
    return res, harm_store


if __name__ == "__main__":
    main()
