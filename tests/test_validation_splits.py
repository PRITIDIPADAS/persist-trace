import pandas as pd


def test_each_line_held_out_once(synthetic_run):
    preds = pd.read_csv(synthetic_run["work"] / "tables" / "benchmark_predictions.csv")
    for (f, m), sub in preds.groupby(["features", "model"]):
        assert sub.groupby("held_out").gsm.nunique().eq(15).all()
        assert sub.gsm.is_unique


def test_test_samples_belong_to_held_out_line(synthetic_run):
    preds = pd.read_csv(synthetic_run["work"] / "tables" / "benchmark_predictions.csv")
    manifest = pd.read_csv(synthetic_run["work"] / "manifests" / "sample_manifest.csv").set_index("gsm")
    assert (manifest.loc[preds.gsm, "cell_line"].values == preds.held_out.values).all()


def test_signatures_built_without_held_out_line():
    import numpy as np
    from persisttrace import benchmark
    lines = {}
    rng = np.random.default_rng(0)
    for l in ["A", "B", "C", "D"]:
        lines[l] = pd.DataFrame({"S": rng.normal(size=50), "se_S": 1.0, "P": rng.normal(size=50), "se_P": 1.0})
    lines["D"]["S"] = 1e6
    mS, mP = benchmark._meta(lines, ["A", "B", "C"])
    assert mS.mu.abs().max() < 100
