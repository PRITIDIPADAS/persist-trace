"""Synthetic counts with known ground truth. Used only by tests and local dry runs.
Nothing produced here is ever exported to the website."""
import numpy as np
import pandas as pd

LIT = ["TACSTD2", "BTG1", "ING4", "PNRC1", "CALCOCO1", "AURKB", "BRD4", "TEAD1", "CD74", "APOBEC3A", "ETV4", "AXL"]

PATTERNS = {
    "stress_shared": dict(acute=2.0, dtp=2.0, wo_short=1.0, wo_long=0.0),
    "acute_only": dict(acute=2.0, dtp=0.0, wo_short=0.0, wo_long=0.0),
    "dtp_reversible": dict(acute=0.0, dtp=2.5, wo_short=2.0, wo_long=0.0),
    "dtp_persistent": dict(acute=0.0, dtp=2.5, wo_short=2.5, wo_long=2.3),
    "dtp_down": dict(acute=-0.3, dtp=-2.0, wo_short=-1.5, wo_long=-0.2),
}


def make(n_genes=1500, seed=7, lines=("PC9", "H1975", "HCC827", "HCC2935")):
    rng = np.random.default_rng(seed)
    conds = ["control", "acute", "dtp", "wo_short", "wo_long"]
    genes = [str(100000 + i) for i in range(n_genes)]
    symbols = [f"SYN{i:04d}" for i in range(n_genes)]
    symbols[: len(LIT)] = LIT
    truth = np.array(["null"] * n_genes, dtype=object)
    k = 60
    start = len(LIT)
    for name in PATTERNS:
        truth[start:start + k] = name
        start += k
    truth[: len(LIT)] = "dtp_persistent"
    base = rng.lognormal(5, 1.2, n_genes)
    samples, cols = [], []
    for li, line in enumerate(lines):
        line_eff = rng.normal(0, 0.4, n_genes)
        for cond in conds:
            for rep in range(1, 4):
                lfc = np.zeros(n_genes)
                if cond != "control":
                    for name, pat in PATTERNS.items():
                        lfc[truth == name] = pat[cond] * (0.4 if (name == "dtp_persistent" and li == 3) else 1.0)
                mu = base * np.exp2(line_eff + lfc) * rng.uniform(0.8, 1.2)
                disp = 0.05
                counts = rng.negative_binomial(1 / disp, 1 / (1 + mu * disp))
                gsm = f"GSM9{li}{conds.index(cond)}{rep:02d}"
                cols.append(pd.Series(counts, index=genes, name=gsm))
                samples.append({"gsm": gsm, "title": f"{line}_{cond}_rep{rep}", "cell_line": line,
                                "condition": cond, "replicate": rep,
                                "washout_days": {"wo_short": 1.0, "wo_long": 3.0 + 2 * li}.get(cond)})
    counts = pd.concat(cols, axis=1)
    annot = pd.DataFrame({"Symbol": symbols}, index=genes)
    return counts, annot, pd.DataFrame(samples), pd.Series(truth, index=genes, name="truth")
