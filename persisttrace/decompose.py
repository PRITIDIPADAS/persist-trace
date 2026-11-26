"""Persistence-specific signal decomposition.

For each cell line and gene g:
    P_g = lfc(DTP vs control),  A_g = lfc(acute vs control)
    P_g = beta * A_g + S_g
beta is a genome-wide robust slope (how much of the acute stress response is still carried
by persisters). S_g is the part of the persister change that the acute response does not explain.
"""
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


def robust_slope(x, y, c=1.345, iters=50):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    b = np.dot(x, y) / max(np.dot(x, x), 1e-12)
    for _ in range(iters):
        r = y - b * x
        s = 1.4826 * np.median(np.abs(r - np.median(r))) or 1.0
        u = np.abs(r / s)
        w = np.where(u <= c, 1.0, c / u)
        nb = np.dot(w * x, y) / max(np.dot(w * x, x), 1e-12)
        if abs(nb - b) < 1e-8:
            break
        b = nb
    return float(b)


def decompose_line(res):
    P, A = res["dtp_vs_control"], res["acute_vs_control"]
    idx = P.index.intersection(A.index)
    P, A = P.loc[idx], A.loc[idx]
    informative = (P.base_mean > 20) & (A.se < 1.5) & (P.se < 1.5)
    beta = robust_slope(A.lfc[informative], P.lfc[informative])
    S = P.lfc - beta * A.lfc
    se = np.sqrt(P.se ** 2 + (beta * A.se) ** 2)
    out = pd.DataFrame({"P": P.lfc, "A": A.lfc, "S": S, "se_S": se, "z_S": S / se,
                        "se_P": P.se, "z_P": P.lfc / P.se})
    return beta, out


def random_effects(est, se):
    """DerSimonian-Laird random-effects meta-analysis across cell lines, row-wise."""
    est, se = np.asarray(est, float), np.asarray(se, float)
    ok = np.isfinite(est) & np.isfinite(se) & (se > 0)
    w = np.where(ok, 1 / np.where(ok, se, 1) ** 2, 0.0)
    e = np.where(ok, est, 0.0)
    k = ok.sum(axis=1)
    sw = w.sum(axis=1)
    fixed = np.divide((w * e).sum(axis=1), sw, out=np.full(len(sw), np.nan), where=sw > 0)
    q = (w * (e - fixed[:, None]) ** 2 * ok).sum(axis=1)
    c = sw - np.divide((w ** 2).sum(axis=1), sw, out=np.zeros(len(sw)), where=sw > 0)
    tau2 = np.maximum(0, np.divide(q - (k - 1), c, out=np.zeros(len(c)), where=c > 0))
    wr = np.where(ok, 1 / (np.where(ok, se, 1) ** 2 + tau2[:, None]), 0.0)
    swr = wr.sum(axis=1)
    mu = np.divide((wr * e).sum(axis=1), swr, out=np.full(len(swr), np.nan), where=swr > 0)
    mu_se = np.sqrt(np.divide(1, swr, out=np.full(len(swr), np.nan), where=swr > 0))
    i2 = np.where(q > 0, np.maximum(0, (q - (k - 1)) / np.where(q > 0, q, 1)), 0.0)
    return mu, mu_se, tau2, i2, k


def meta(per_line, value, se_col, lines=None):
    lines = lines or list(per_line)
    genes = sorted(set.intersection(*[set(per_line[l].index) for l in lines]))
    est = np.column_stack([per_line[l].loc[genes, value] for l in lines])
    se = np.column_stack([per_line[l].loc[genes, se_col] for l in lines])
    mu, mu_se, tau2, i2, k = random_effects(est, se)
    z = mu / mu_se
    p = 2 * stats.norm.sf(np.abs(z))
    p = np.where(np.isfinite(p), p, 1.0)
    sig = np.abs(est / se) > 2
    agree_up = ((est > 0) & sig).sum(axis=1)
    agree_dn = ((est < 0) & sig).sum(axis=1)
    return pd.DataFrame({"mu": mu, "se": mu_se, "z": z, "p": p,
                         "fdr": multipletests(p, method="fdr_bh")[1],
                         "tau2": tau2, "I2": i2, "k": k,
                         "lines_up": agree_up, "lines_down": agree_dn}, index=genes)


def run(results, lines=None):
    lines = lines or list(results)
    betas, per_line = {}, {}
    for line in lines:
        betas[line], per_line[line] = decompose_line(results[line])
    m_S = meta(per_line, "S", "se_S", lines)
    m_P = meta(per_line, "P", "se_P", lines)
    return betas, per_line, m_S, m_P


def signature(meta_df, n, min_lines=3, anchor=None):
    """Top up and down genes. With anchor (the treated-vs-control meta table), a gene only counts
    if the persister change itself is significant in the same direction, so genes that merely
    lost the acute response are not called persister-specific."""
    df = meta_df[np.isfinite(meta_df.z)]
    if anchor is not None:
        a = anchor.reindex(df.index)
        df = df[(np.sign(a.z) == np.sign(df.z)) & (a.z.abs() > 2)]
    up = df[(df.lines_up >= min_lines)].sort_values("z", ascending=False).head(n // 2)
    dn = df[(df.lines_down >= min_lines)].sort_values("z").head(n // 2)
    if len(up) < n // 4:
        up = df.sort_values("z", ascending=False).head(n // 2)
    if len(dn) < n // 4:
        dn = df.sort_values("z").head(n // 2)
    return list(up.index), list(dn.index)
