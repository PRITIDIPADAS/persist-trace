import json

import pytest

from persisttrace import pipeline


@pytest.fixture(scope="session")
def synthetic_run(tmp_path_factory):
    work = tmp_path_factory.mktemp("work")
    site = tmp_path_factory.mktemp("site")
    pipeline.main(["--source", "synthetic", "--work", str(work), "--site", str(site),
                   "--skip-networks", "--skip-media"])
    res = json.loads((site / "data" / "results_synthetic.json").read_text())
    genes = json.loads((site / "data" / "genes_synthetic.json").read_text())
    return {"work": work, "site": site, "res": res, "genes": genes}
