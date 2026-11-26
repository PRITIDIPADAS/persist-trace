import numpy as np

from persisttrace import model


def test_drug_on_shrinks_and_withdrawal_regrows():
    t, s, p = model.simulate()
    total = s + p
    i21 = np.searchsorted(t, 21)
    assert total[i21] < total[0]
    assert total[-1] > total[i21]
    assert p[i21] / total[i21] > p[0] / total[0]


def test_growth_rate_is_dominant_eigenvalue():
    for phase in ("on", "off"):
        m = model.matrix(model.DEFAULT[phase])
        assert np.isclose(model.growth_rate(model.DEFAULT[phase]), max(np.linalg.eigvals(m).real))
