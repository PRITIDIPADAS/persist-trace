import pandas as pd
import pytest

from persisttrace import samples


@pytest.mark.parametrize("text,line,cond", [
    ("PC9_DMSO_rep1", "PC9", "control"),
    ("NCI-H1975 osimertinib 24h", "H1975", "acute"),
    ("HCC827 DTP 21d", "HCC827", "dtp"),
    ("PC9 washout 24h", "PC9", "wo_short"),
    ("HCC2935 osimertinib DTP washout 7 days", "HCC2935", "wo_long"),
])
def test_labels(text, line, cond):
    assert samples.cell_line(text) == line
    assert samples.condition(text) == cond


def test_unbalanced_design_fails():
    design = {"cell_lines": ["PC9"], "conditions": ["control", "dtp"], "expected_replicates": 3}
    df = pd.DataFrame({"gsm": ["a", "b", "c", "d"], "cell_line": "PC9",
                       "condition": ["control"] * 3 + ["dtp"]})
    with pytest.raises(samples.MetadataError):
        samples.validate(df, design)


def test_sample_counts_not_cell_counts(synthetic_run):
    meta = synthetic_run["res"]["meta"]
    assert meta["n_samples"] == 60
    assert len(synthetic_run["res"]["samples"]) == 60
