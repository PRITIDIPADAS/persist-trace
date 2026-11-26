import logging
import warnings

import numpy as np
import pandas as pd

CONTRASTS = {
    "acute_vs_control": ("acute", "control"),
    "dtp_vs_control": ("dtp", "control"),
    "dtp_vs_acute": ("dtp", "acute"),
    "wo_short_vs_control": ("wo_short", "control"),
    "wo_long_vs_control": ("wo_long", "control"),
    "wo_short_vs_dtp": ("wo_short", "dtp"),
    "wo_long_vs_dtp": ("wo_long", "dtp"),
}


def filter_genes(counts, min_count, min_samples):
    return counts[(counts >= min_count).sum(axis=1) >= min_samples]


def log_cpm(counts, prior=1.0):
    lib = counts.sum(axis=0)
    return np.log2((counts + prior).div(lib + 2 * prior, axis=1) * 1e6)


def deseq_line(counts, samples, seed=0):
    """DESeq2 (pydeseq2) Wald tests for every contrast within one cell line."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats

    meta = samples.set_index("gsm")[["condition"]].copy()
    mat = counts[meta.index].T
    inference = DefaultInference(n_cpus=1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        logging.getLogger("pydeseq2").setLevel(logging.ERROR)
        dds = DeseqDataSet(counts=mat, metadata=meta, design="~condition", refit_cooks=True,
                           inference=inference, quiet=True)
        dds.deseq2()
        out = {}
        for name, (num, den) in CONTRASTS.items():
            st = DeseqStats(dds, contrast=["condition", num, den], inference=inference, quiet=True)
            st.summary()
            r = st.results_df[["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]].copy()
            r.columns = ["base_mean", "lfc", "se", "stat", "p", "padj"]
            r["padj"] = r["padj"].fillna(1.0)
            r["p"] = r["p"].fillna(1.0)
            out[name] = r
    return out


def run_all(counts, samples, cfg):
    f = cfg["filter"]
    results, kept = {}, {}
    for line, sub in samples.groupby("cell_line"):
        c = filter_genes(counts[list(sub.gsm)], f["min_count"], f["min_samples"])
        kept[line] = c.index
        results[line] = deseq_line(c, sub)
    return results, kept


def stack(results, contrast, field):
    return pd.DataFrame({line: res[contrast][field] for line, res in results.items()})
