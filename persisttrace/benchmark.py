"""Leave-one-cell-line-out benchmark.

Signatures, feature selection and model fitting only ever see the training cell lines.
The held-out line is centred on its own sample mean, which uses no labels."""
import warnings

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, f1_score, roc_auc_score
from sklearn.neighbors import NearestCentroid
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from . import decompose


def _models(seed):
    m = {
        "nearest_centroid": lambda: NearestCentroid(),
        "logistic_l2": lambda: make_pipeline(StandardScaler(), LogisticRegression(C=0.05, max_iter=5000)),
        "elastic_net": lambda: make_pipeline(StandardScaler(), LogisticRegression(
            penalty="elasticnet", solver="saga", l1_ratio=0.5, C=0.1, max_iter=5000)),
        "random_forest": lambda: RandomForestClassifier(n_estimators=500, random_state=seed),
    }
    try:
        from xgboost import XGBClassifier
        m["xgboost"] = lambda: XGBClassifier(n_estimators=300, max_depth=3, learning_rate=0.05,
                                             subsample=0.8, random_state=seed, verbosity=0)
    except ImportError:
        pass
    return m


def centre_by_line(logcpm, samples):
    out = logcpm.copy()
    for _, sub in samples.groupby("cell_line"):
        cols = list(sub.gsm)
        out[cols] = out[cols].sub(out[cols].mean(axis=1), axis=0)
    return out


def score(z, up, dn):
    up = [g for g in up if g in z.index]
    dn = [g for g in dn if g in z.index]
    return z.loc[up].mean(axis=0) - z.loc[dn].mean(axis=0)


def _auc(y, s):
    try:
        return float(roc_auc_score(y, s))
    except ValueError:
        return None


def run(logcpm, samples, per_line, cfg):
    seed = cfg["seed"]
    n = cfg["benchmark"]["signature_size"]
    lines = sorted(samples.cell_line.unique())
    genes = sorted(set.intersection(*[set(per_line[l].index) for l in lines]) & set(logcpm.index))
    z = centre_by_line(logcpm.loc[genes], samples)
    classes = cfg["design"]["conditions"]
    y_all = samples.set_index("gsm").condition
    models = _models(seed)
    folds, preds = [], []
    for held in lines:
        train = [l for l in lines if l != held]
        mS, mP = _meta(per_line, train)
        sigs = {"naive_treated_vs_control": decompose.signature(mP, n),
                "persister_specific": decompose.signature(mS, n, anchor=mP)}
        test = samples[samples.cell_line == held]
        tr = samples[samples.cell_line.isin(train)]
        fold = {"held_out": held}
        for name, (up, dn) in sigs.items():
            sc = score(z[list(test.gsm)], up, dn)
            cond = test.set_index("gsm").condition.loc[sc.index]
            for a, b in [("dtp", "acute"), ("dtp", "control"), ("dtp", "wo_long")]:
                m = cond.isin([a, b])
                fold[f"{name}|auroc_{a}_vs_{b}"] = _auc((cond[m] == a).astype(int), sc[m])
            held_line = per_line[held]
            col = "z_S" if name == "persister_specific" else "z_P"
            top = up + dn
            sign = pd.Series([1] * len(up) + [-1] * len(dn), index=top)
            hz = held_line.reindex(top)[col]
            fold[f"{name}|replicated_fraction"] = float(((np.sign(hz) == sign) & (hz.abs() > 2)).mean())
            train_meta = (mS if name == "persister_specific" else mP).z
            common = train_meta.index.intersection(held_line.index)
            rho = stats.spearmanr(train_meta.loc[common], held_line.loc[common, col], nan_policy="omit")[0]
            fold[f"{name}|spearman_genomewide"] = float(rho)
        var_genes = z[list(tr.gsm)].var(axis=1).sort_values(ascending=False).head(cfg["benchmark"]["top_variable"]).index
        feats = {"top_variable": list(var_genes),
                 "naive_signature": sigs["naive_treated_vs_control"][0] + sigs["naive_treated_vs_control"][1],
                 "persister_signature": sigs["persister_specific"][0] + sigs["persister_specific"][1]}
        ytr = y_all.loc[tr.gsm].map(classes.index).values
        yte = y_all.loc[test.gsm].map(classes.index).values
        for fname, g in feats.items():
            Xtr, Xte = z.loc[g, list(tr.gsm)].T.values, z.loc[g, list(test.gsm)].T.values
            for mname, make in models.items():
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    clf = make().fit(Xtr, ytr)
                    pr = clf.predict(Xte)
                for gsm, t, p in zip(test.gsm, yte, pr):
                    preds.append({"held_out": held, "features": fname, "model": mname, "gsm": gsm,
                                  "true": classes[t], "pred": classes[int(p)]})
        folds.append(fold)
    preds = pd.DataFrame(preds)
    return pd.DataFrame(folds), preds, summarize(preds, classes, seed)


def _meta(per_line, lines):
    m_S = decompose.meta(per_line, "S", "se_S", lines)
    m_P = decompose.meta(per_line, "P", "se_P", lines)
    return m_S, m_P


def summarize(preds, classes, seed, boots=1000):
    rng = np.random.default_rng(seed)
    rows = []
    for (f, m), sub in preds.groupby(["features", "model"]):
        t, p = sub.true.values, sub.pred.values
        ba = balanced_accuracy_score(t, p)
        f1 = f1_score(t, p, average="macro", labels=classes, zero_division=0)
        bs = []
        for _ in range(boots):
            i = rng.integers(0, len(t), len(t))
            bs.append(balanced_accuracy_score(t[i], p[i]))
        per_fold = sub.groupby("held_out").apply(lambda d: balanced_accuracy_score(d.true, d.pred)).to_dict()
        dtp_acute = sub[sub.true.isin(["dtp", "acute"])]
        confusion = float(((dtp_acute.true == "dtp") & (dtp_acute.pred == "acute")).sum() +
                          ((dtp_acute.true == "acute") & (dtp_acute.pred == "dtp")).sum()) / max(len(dtp_acute), 1)
        rows.append({"features": f, "model": m, "balanced_accuracy": ba, "macro_f1": f1,
                     "ci_low": float(np.percentile(bs, 2.5)), "ci_high": float(np.percentile(bs, 97.5)),
                     "dtp_acute_confusion": confusion, "per_fold": per_fold, "n": len(t)})
    return pd.DataFrame(rows).sort_values("balanced_accuracy", ascending=False)
