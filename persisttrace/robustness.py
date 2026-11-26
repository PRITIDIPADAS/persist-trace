import numpy as np
import pandas as pd

from . import decompose, reversibility


def _jaccard(a, b):
    a, b = set(a), set(b)
    return len(a & b) / max(len(a | b), 1)


def threshold_grid(results, cfg, ref_cls):
    rows = []
    ref = ref_cls["class"]
    altered = ref[ref.str.startswith("dtp_")].index
    for a in cfg["robustness"]["alpha_grid"]:
        for l in cfg["robustness"]["lfc_grid"]:
            per = {line: reversibility.classify_line(res, a, l) for line, res in results.items()}
            cons = reversibility.consensus(per, cfg["thresholds"]["reproducible_lines"])["class"]
            same = (cons.reindex(altered) == ref.reindex(altered)).mean() if len(altered) else np.nan
            rows.append({"alpha": a, "lfc": l, "agreement_with_reference": float(same),
                         "n_dtp_altered": int(cons.str.startswith("dtp_").sum())})
    return pd.DataFrame(rows)


def leave_one_line_out(per_line, n):
    lines = list(per_line)
    full = decompose.meta(per_line, "S", "se_S", lines)
    up, dn = decompose.signature(full, n, anchor=decompose.meta(per_line, "P", "se_P", lines))
    ref = set(up + dn)
    rows = []
    for drop in lines:
        keep = [l for l in lines if l != drop]
        m = decompose.meta(per_line, "S", "se_S", keep)
        u, d = decompose.signature(m, n, min_lines=2, anchor=decompose.meta(per_line, "P", "se_P", keep))
        rows.append({"dropped": drop, "jaccard_top": _jaccard(ref, u + d)})
    return pd.DataFrame(rows)


def beta_sensitivity(results):
    rows = []
    for line, res in results.items():
        P, A = res["dtp_vs_control"], res["acute_vs_control"]
        idx = P.index.intersection(A.index)
        x, y = A.lfc.loc[idx].values, P.lfc.loc[idx].values
        ols = float(np.dot(x, y) / np.dot(x, x))
        rob = decompose.robust_slope(x, y)
        S_ols, S_rob = y - ols * x, y - rob * x
        r = float(np.corrcoef(S_ols, S_rob)[0, 1])
        rows.append({"cell_line": line, "beta_ols": ols, "beta_robust": rob, "residual_correlation": r})
    return pd.DataFrame(rows)


def label(value, stable=0.8, moderate=0.6):
    if value is None or not np.isfinite(value):
        return "not evaluated"
    return "stable" if value >= stable else ("moderately sensitive" if value >= moderate else "highly sensitive")
