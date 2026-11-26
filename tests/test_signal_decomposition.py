import numpy as np
import pandas as pd

from persisttrace import decompose, synthetic


def test_robust_slope_recovers_beta():
    rng = np.random.default_rng(0)
    x = rng.normal(0, 1, 5000)
    y = 0.7 * x + rng.normal(0, 0.1, 5000)
    y[:200] += 5
    assert abs(decompose.robust_slope(x, y) - 0.7) < 0.03


def test_random_effects_matches_fixed_when_homogeneous():
    est = np.array([[1.0, 1.0, 1.0, 1.0]])
    se = np.array([[0.5, 0.5, 0.5, 0.5]])
    mu, mu_se, tau2, i2, k = decompose.random_effects(est, se)
    assert np.isclose(mu[0], 1.0) and np.isclose(mu_se[0], 0.25) and tau2[0] == 0 and k[0] == 4


def test_stress_genes_shrink_more_than_persister_genes(synthetic_run):
    work = synthetic_run["work"]
    mS = pd.read_csv(work / "tables" / "meta_persister_specific.csv.gz", index_col=0)
    mP = pd.read_csv(work / "tables" / "meta_treated_vs_control.csv.gz", index_col=0)
    truth = synthetic.make()[3]
    truth.index = truth.index.astype(int)
    df = pd.DataFrame({"S": mS.mu, "P": mP.mu}).join(truth)
    g = df.groupby("truth").mean()
    shrink = 1 - g.S / g.P
    assert shrink["stress_shared"] > 0.3
    assert abs(shrink["dtp_persistent"]) < 0.15


def test_signature_anchor_excludes_acute_only(synthetic_run):
    work = synthetic_run["work"]
    mS = pd.read_csv(work / "tables" / "meta_persister_specific.csv.gz", index_col=0)
    mP = pd.read_csv(work / "tables" / "meta_treated_vs_control.csv.gz", index_col=0)
    truth = synthetic.make()[3]
    truth.index = truth.index.astype(int)
    up, dn = decompose.signature(mS, 100, anchor=mP)
    assert (truth.reindex(up + dn) == "acute_only").sum() <= 2
