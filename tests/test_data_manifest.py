import hashlib

from persisttrace import geo


def test_manifest_records_checksum(tmp_path):
    f = tmp_path / "x.txt"
    f.write_bytes(b"persist")
    m = geo.Manifest(tmp_path / "m.csv")
    m.add("GSE0", "test", "https://example.org/x", f)
    m.write()
    assert m.rows[0]["sha256"] == hashlib.sha256(b"persist").hexdigest()
    assert (tmp_path / "m.sha256").read_text().startswith(m.rows[0]["sha256"])


def test_series_matrix_parser(tmp_path):
    text = "\n".join([
        '!Sample_title\t"PC9_DMSO_1"\t"PC9_DTP_1"',
        '!Sample_geo_accession\t"GSM1"\t"GSM2"',
        '!Sample_characteristics_ch1\t"cell line: PC9"\t"cell line: PC9"',
        '!Sample_characteristics_ch1\t"treatment: DMSO"\t"treatment: osimertinib 21 days"',
    ])
    p = tmp_path / "m.txt"
    p.write_text(text)
    df = geo.parse_series_matrix(p)
    assert list(df.geo_accession) == ["GSM1", "GSM2"]
    assert df.loc[1, "ch_treatment"] == "osimertinib 21 days"


def test_synthetic_run_has_no_manifest(synthetic_run):
    assert synthetic_run["res"]["meta"]["synthetic"] is True
    assert synthetic_run["res"]["gates"]["A_data_integrity"] == "synthetic run"


def test_author_counts_matched_by_title(tmp_path):
    import pandas as pd
    samples = pd.DataFrame({"gsm": ["GSM1", "GSM2"], "title": ["PC9_DMSO_1 [RNA-seq]", "PC9_osi_DTP_1 [RNA-seq]"]})
    t = pd.DataFrame({"gene": ["ENSG01.3", "ENSG02.1"], "PC9_DMSO_1": [10.4, 0], "PC9_osi_DTP_1": [3, 7.6], "extra": [1, 1]})
    p = tmp_path / "counts.txt.gz"
    t.to_csv(p, sep="\t", index=False)
    counts, report, name = geo.read_author_counts([p], samples)
    assert name == "counts.txt.gz"
    assert list(counts.columns) == ["GSM1", "GSM2"]
    assert counts.loc["ENSG01.3", "GSM1"] == 10
