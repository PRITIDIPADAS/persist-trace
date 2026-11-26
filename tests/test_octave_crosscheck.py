import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from persisttrace import model, reversibility

OCT = Path(os.environ.get("PT_OCTAVE_DIR", "build/octave"))


@pytest.mark.skipif(not (OCT / "model_octave.csv").exists(), reason="Octave output not present")
def test_octave_model_matches_python():
    o = pd.read_csv(OCT / "model_octave.csv").set_index("quantity").value
    t, s, p = model.simulate()
    assert np.isclose(o["lambda_on"], model.growth_rate(model.DEFAULT["on"]), rtol=1e-6)
    assert np.isclose(o["lambda_off"], model.growth_rate(model.DEFAULT["off"]), rtol=1e-6)
    assert np.isclose(o["total_day42"], s[-1] + p[-1], rtol=1e-6)


@pytest.mark.skipif(not (OCT / "memory_octave.csv").exists(), reason="Octave memory fit not present")
def test_octave_memory_fit_matches_python():
    inp = pd.read_csv(OCT / "memory_input.csv")
    out = pd.read_csv(OCT / "memory_octave.csv")
    k_py, _ = reversibility.fit_decay(inp.m_short.values, inp.m_long.values, 1.0, inp.t_long.values)
    both = np.isfinite(k_py) & np.isfinite(out.k.values)
    assert (np.isfinite(k_py) == np.isfinite(out.k.values)).mean() > 0.999
    assert np.allclose(k_py[both], out.k.values[both], rtol=1e-6)
