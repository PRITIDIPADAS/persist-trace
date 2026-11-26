import gzip
import io
import tarfile

import numpy as np
import pandas as pd
from scipy import io as sio
from scipy import sparse

from persisttrace import external, geo


def _mtx(m, path):
    buf = io.BytesIO()
    sio.mmwrite(buf, sparse.coo_matrix(m))
    with gzip.open(path, "wb") as f:
        f.write(buf.getvalue())


def _fake_cells(rng, genes, n, shift):
    base = rng.gamma(2, 2, len(genes))
    m = rng.poisson(base[:, None] * np.exp(shift[:, None]) * np.ones((1, n)))
    return m


def test_external_modules_run_on_fake_files(tmp_path, monkeypatch):
    rng = np.random.default_rng(0)
    genes = [f"G{i}" for i in range(600)] + external.MAPK_OUTPUT + ["EGFR"] + [f"MT-{i}" for i in range(3)]
    sets = {"HALLMARK_E2F_TARGETS": set(genes[:40]), "HALLMARK_G2M_CHECKPOINT": set(genes[40:80])}
    sigs = {"naive_treated_vs_control": (genes[100:130], genes[130:160]),
            "persister_specific": (genes[200:230], genes[230:260])}
    files = {}
    n = len(genes)
    sh_a = np.zeros(n); sh_b = np.zeros(n)
    sh_b[[genes.index(g) for g in external.MAPK_OUTPUT]] = -2
    sh_b[200:230] = 1.5
    m = np.hstack([_fake_cells(rng, genes, 300, sh_a), _fake_cells(rng, genes, 300, sh_b)])
    _mtx(m, tmp_path / "m.mtx.gz"); files["GSE275199_matrix.mtx.gz"] = tmp_path / "m.mtx.gz"
    pd.DataFrame({0: genes, 1: genes, 2: "Gene Expression"}).to_csv(tmp_path / "f.tsv.gz", sep="\t", header=False, index=False)
    files["GSE275199_features.tsv.gz"] = tmp_path / "f.tsv.gz"; files["GSE275199_barcodes.tsv.gz"] = tmp_path / "f.tsv.gz"
    hg = ["hg38_" + g for g in genes] + ["mm10_Actb", "mm10_Gapdh"]
    with tarfile.open(tmp_path / "raw.tar", "w") as t:
        for i, gsm in enumerate(["GSM7790998", "GSM7790999", "GSM7791000", "GSM7791001"]):
            sh = sh_b if i % 2 else sh_a
            mm = np.vstack([_fake_cells(rng, genes, 200, sh), np.zeros((2, 200), int)])
            _mtx(mm, tmp_path / "x.mtx.gz")
            t.add(tmp_path / "x.mtx.gz", arcname=f"{gsm}_1_matrix.mtx.gz")
            pd.DataFrame({0: hg, 1: hg, 2: "Gene Expression"}).to_csv(tmp_path / "y.tsv.gz", sep="\t", header=False, index=False)
            t.add(tmp_path / "y.tsv.gz", arcname=f"{gsm}_1_features.tsv.gz")
    files["GSE243562_RAW.tar"] = tmp_path / "raw.tar"
    cols = ["Day 0 G1 clone 1", "Day 0 S clone 1", "Day 0 G2 clone 1", "Day 28 G1 clone 1", "Day 28 S clone 1", "Day 28 G2 clone 1"]
    cc = pd.DataFrame(rng.normal(6, 1, (len(genes), 6)), columns=cols)
    cc.insert(0, "Feature", genes)
    cc.to_csv(tmp_path / "cc.txt.gz", sep="\t", index=False)
    files["GSE335846_RNA-seq_for_GEO.txt.gz"] = tmp_path / "cc.txt.gz"
    monkeypatch.setattr(external, "_fetch", lambda acc, name, raw, man: files[name])
    out = external.run(tmp_path, geo.Manifest(tmp_path / "m.csv"), sigs, sets, 0)
    assert out["errors"] == {}, out["errors"]
    sc = out["pc9_single_cell"]
    assert sc["exposed_cells"] > 0 and len(sc["points"]) > 0
    assert "persister_specific" in sc["auroc_exposed_vs_unexposed"]
    assert set(out["pdx_single_cell"]["auroc_residual_vs_untreated"]["persister_specific"]) == {"YU-003", "YU-006"}
    assert len(out["cell_cycle_sorted"]["summary"]) == 6
