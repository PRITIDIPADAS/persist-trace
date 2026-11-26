"""Two-state switching model of sensitive (S) and persister (P) cells.

    d/dt [S, P] = M(D) [S, P],   M(D) = [[r_S(D) - a(D),  b(D)],
                                         [a(D),           r_P(D) - b(D)]]

D = 1 while the drug is on, 0 after withdrawal. This is a theory panel with illustrative
parameters. It is not fitted to cell counts, because GSE193258 has no cell counts."""
import numpy as np
from scipy.linalg import expm

DEFAULT = {
    "on": {"r_S": -0.8, "r_P": 0.0, "a": 0.05, "b": 0.0},
    "off": {"r_S": 0.6, "r_P": 0.1, "a": 0.001, "b": 0.3},
}


def matrix(p):
    return np.array([[p["r_S"] - p["a"], p["b"]], [p["a"], p["r_P"] - p["b"]]])


def growth_rate(p):
    m = matrix(p)
    tr, det = np.trace(m), np.linalg.det(m)
    return float((tr + np.sqrt(tr ** 2 - 4 * det + 0j)).real / 2)


def simulate(params=DEFAULT, on_days=21, off_days=21, dt=0.25, s0=0.999, p0=0.001):
    t, s, p = [0.0], [s0], [p0]
    x = np.array([s0, p0])
    for phase, days in (("on", on_days), ("off", off_days)):
        step = expm(matrix(params[phase]) * dt)
        for _ in range(int(round(days / dt))):
            x = step @ x
            t.append(t[-1] + dt)
            s.append(x[0])
            p.append(x[1])
    return np.array(t), np.array(s), np.array(p)
