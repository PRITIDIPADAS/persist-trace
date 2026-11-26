"""Population-level reversibility classes and transcriptional memory.

A bulk sample returning toward baseline is compatible with reversal, but it does not show that
the same cells reverted. Every label here is a population-level statement.
"""
import numpy as np
import pandas as pd

CLASSES = [
    "acute_only",
    "dtp_reversed_by_24h",
    "dtp_persists_24h_reversed_late",
    "dtp_persists_after_regrowth",
    "dtp_overshoot",
    "dtp_unresolved",
    "not_changed",
]


def _sig(r, alpha, lfc):
    return (r.padj < alpha) & (r.lfc.abs() >= lfc)


def classify_line(res, alpha=0.05, lfc=1.0):
    g = res["dtp_vs_control"].index
    get = {k: v.reindex(g) for k, v in res.items()}
    a, d = get["acute_vs_control"], get["dtp_vs_control"]
    ws, wl = get["wo_short_vs_control"], get["wo_long_vs_control"]
    wld = get["wo_long_vs_dtp"]
    sig_a, sig_d = _sig(a, alpha, lfc), _sig(d, alpha, lfc)
    sign = np.sign(d.lfc)
    still_s = _sig(ws, alpha, lfc / 2) & (np.sign(ws.lfc) == sign)
    still_l = _sig(wl, alpha, lfc / 2) & (np.sign(wl.lfc) == sign)
    over = _sig(wl, alpha, lfc / 2) & (np.sign(wl.lfc) == -sign)
    moved_back = (wld.padj < alpha) & (np.sign(wld.lfc) == -sign)
    cls = pd.Series("not_changed", index=g, dtype=object)
    cls[sig_a & ~sig_d] = "acute_only"
    dtp = sig_d
    cls[dtp] = "dtp_unresolved"
    cls[dtp & ~still_s & ~still_l] = "dtp_reversed_by_24h"
    cls[dtp & still_s & ~still_l & moved_back] = "dtp_persists_24h_reversed_late"
    cls[dtp & still_l] = "dtp_persists_after_regrowth"
    cls[dtp & over] = "dtp_overshoot"
    return cls


def consensus(per_line_cls, min_lines=3):
    df = pd.DataFrame(per_line_cls)
    def pick(row):
        vals = row.dropna()
        top = vals.value_counts()
        if len(top) and top.iloc[0] >= min_lines:
            return top.index[0], "reproducible"
        changed = vals[vals != "not_changed"]
        if len(changed) >= 2:
            return "context_dependent", "context_dependent"
        if len(changed) == 1:
            return "single_line", "single_line"
        return "not_changed", "reproducible"
    picked = df.apply(pick, axis=1, result_type="expand")
    picked.columns = ["class", "status"]
    return df.join(picked)


def memory(res):
    """Fraction of the persister change still present after washout: m = lfc(wo vs ctrl) / lfc(dtp vs ctrl)."""
    d = res["dtp_vs_control"]
    out = pd.DataFrame(index=d.index)
    for k in ["wo_short", "wo_long"]:
        w = res[f"{k}_vs_control"].reindex(d.index)
        out[f"m_{k}"] = w.lfc / d.lfc
    out["dtp_lfc"] = d.lfc
    out["dtp_padj"] = d.padj
    return out


def fit_decay(m_short, m_long, t_short, t_long):
    """Least-squares fit of m(t) = exp(-k t) through m(0)=1 and the two washout points.

    Returns k (per day) and half-life, or NaN where the data are not monotone decay."""
    ts, tl = np.asarray(t_short, float), np.asarray(t_long, float)
    ms, ml = np.asarray(m_short, float), np.asarray(m_long, float)
    ok = (ms > 0.02) & (ml > 0.02) & (ms <= 1.5) & (ml <= 1.5)
    ys = np.log(np.clip(ms, 1e-6, None))
    yl = np.log(np.clip(ml, 1e-6, None))
    k = -(ts * ys + tl * yl) / (ts ** 2 + tl ** 2)
    k = np.where(ok & (k > 0), k, np.nan)
    return k, np.log(2) / k


def summarize_memory(mem_by_line, alpha=0.05, lfc=1.0):
    rows = []
    for line, m in mem_by_line.items():
        sel = m[(m.dtp_padj < alpha) & (m.dtp_lfc.abs() >= lfc)]
        for k in ["m_wo_short", "m_wo_long"]:
            v = sel[k].replace([np.inf, -np.inf], np.nan).dropna().clip(-1, 2)
            rows.append({"cell_line": line, "timepoint": k[2:], "n_genes": int(len(v)),
                         "median_memory": float(v.median()) if len(v) else None,
                         "q25": float(v.quantile(.25)) if len(v) else None,
                         "q75": float(v.quantile(.75)) if len(v) else None,
                         "hist": np.histogram(v, bins=30, range=(-1, 2))[0].tolist()})
    return rows
