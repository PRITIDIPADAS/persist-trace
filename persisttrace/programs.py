import io
import json
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

ENRICHR = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=MSigDB_Hallmark_2020"


def _cache(path, loader):
    path = Path(path)
    if path.exists():
        return pd.read_csv(path), "cache"
    df, origin = loader()
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    return df, origin


def load_hallmark(cache):
    def loader():
        try:
            import decoupler as dc
            df = dc.op.hallmark(organism="human")[["source", "target"]]
            return df, "decoupler/omnipath MSigDB hallmark"
        except Exception:
            req = urllib.request.Request(ENRICHR, headers={"User-Agent": "PersistTrace"})
            text = urllib.request.urlopen(req, timeout=120).read().decode()
            rows = []
            for line in text.splitlines():
                parts = [p.split(",")[0] for p in line.split("\t") if p]
                rows += [(parts[0], g) for g in parts[1:]]
            return pd.DataFrame(rows, columns=["source", "target"]), "Enrichr MSigDB_Hallmark_2020"
    return _cache(cache, loader)


def load_collectri(cache):
    def loader():
        import decoupler as dc
        df = dc.op.collectri(organism="human")
        return df[["source", "target", "weight"]], "decoupler/omnipath CollecTRI"
    return _cache(cache, loader)


def ulm(values, net, min_n=5):
    """Univariate linear model: t-value of the slope of gene statistic on regulon weight."""
    v = values.dropna()
    net = net[net.target.isin(v.index)]
    rows = []
    y = v.values - v.values.mean()
    for src, sub in net.groupby("source"):
        if len(sub) < min_n:
            continue
        x = pd.Series(0.0, index=v.index)
        w = sub.groupby("target")["weight"].mean() if "weight" in sub else pd.Series(1.0, index=sub.target.unique())
        x.loc[w.index] = w.values
        xc = x.values - x.values.mean()
        sxx = np.dot(xc, xc)
        b = np.dot(xc, y) / sxx
        resid = y - b * xc
        df = len(y) - 2
        se = np.sqrt(np.dot(resid, resid) / df / sxx)
        t = b / se
        rows.append((src, len(w), t, 2 * stats.t.sf(abs(t), df)))
    out = pd.DataFrame(rows, columns=["source", "n_targets", "score", "p"]).set_index("source")
    if len(out):
        out["fdr"] = multipletests(out.p, method="fdr_bh")[1]
    return out


def program_scores(stat_by_symbol, sets):
    s = sets.assign(weight=1.0)
    return ulm(stat_by_symbol, s, min_n=10)


def sample_program_matrix(logcpm_sym, sets, programs):
    z = logcpm_sym.sub(logcpm_sym.mean(axis=1), axis=0).div(logcpm_sym.std(axis=1) + 1e-9, axis=0)
    out = {}
    for p in programs:
        genes = sets.loc[sets.source == p, "target"]
        genes = [g for g in genes if g in z.index]
        if len(genes) >= 10:
            out[p] = z.loc[genes].mean(axis=0)
    return pd.DataFrame(out).T
