import numpy as np
import pandas as pd

from persisttrace import reversibility, synthetic


def _res(lfcs, padj=1e-6):
    idx = ["g"]
    return {k: pd.DataFrame({"lfc": [v], "padj": [padj]}, index=idx) for k, v in lfcs.items()}


def test_missing_long_washout_is_not_called_reversed():
    r = _res({"acute_vs_control": 0, "dtp_vs_control": 3, "wo_short_vs_control": 3,
              "wo_long_vs_control": np.nan, "wo_long_vs_dtp": np.nan})
    r["wo_long_vs_control"]["padj"] = np.nan
    r["wo_long_vs_dtp"]["padj"] = np.nan
    cls = reversibility.classify_line(r)
    assert cls["g"] != "dtp_persists_24h_reversed_late"


def test_persistent_gene_classified():
    r = _res({"acute_vs_control": 0, "dtp_vs_control": 3, "wo_short_vs_control": 3,
              "wo_long_vs_control": 2.5, "wo_long_vs_dtp": -0.5})
    assert reversibility.classify_line(r)["g"] == "dtp_persists_after_regrowth"


def test_decay_fit_recovers_rate():
    k_true = 0.4
    k, hl = reversibility.fit_decay(np.array([np.exp(-k_true)]), np.array([np.exp(-k_true * 5)]), 1, 5)
    assert np.isclose(k[0], k_true) and np.isclose(hl[0], np.log(2) / k_true)


def test_synthetic_classes_track_truth(synthetic_run):
    cls = pd.read_csv(synthetic_run["work"] / "tables" / "reversibility_classes.csv.gz", index_col=0)
    truth = synthetic.make()[3]
    truth.index = truth.index.astype(int)
    df = cls.join(truth)
    assert (df[df.truth == "acute_only"]["class"] == "acute_only").mean() > 0.9
    assert (df[df.truth == "null"]["class"] == "not_changed").mean() > 0.95
