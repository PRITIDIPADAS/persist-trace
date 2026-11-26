"""Validation outside GSE193258: cell-cycle sorted persisters, single PC9 cells and PDX tumours.

Signatures come only from the four-line bulk analysis and are frozen before these data are read."""
import gzip
import io
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import io as sio
from scipy import sparse, stats
from sklearn.metrics import roc_auc_score
from sklearn.mixture import GaussianMixture

from . import geo

FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{acc}/suppl/{name}"
MAPK_OUTPUT = ["DUSP6", "SPRY4", "SPRY2", "ETV4", "ETV5", "FOSL1", "EGR1", "PHLDA1", "CCND1", "MYC", "DUSP4", "IER3"]


def _url(acc, name):
    return FTP.format(stub=geo.series_stub(acc), acc=acc, name=name)


def _fetch(acc, name, raw, manifest):
    url = _url(acc, name)
    path = geo.download(url, Path(raw) / name)
    manifest.add(acc, "external_validation", url, path)
    return path


def score_matrix(z, up, dn):
    up = [g for g in up if g in z.index]
    dn = [g for g in dn if g in z.index]
    s = np.zeros(z.shape[1])
    if up:
        s = s + np.asarray(z.loc[up].mean(axis=0)).ravel()
    if dn:
        s = s - np.asarray(z.loc[dn].mean(axis=0)).ravel()
    return s, len(up), len(dn)


def _auc(y, s):
    try:
        return float(roc_auc_score(y, s))
    except ValueError:
        return None


def cell_cycle_sorted(raw, manifest, sigs):
    """GSE335846: PC9 sorted into G1, S and G2 before and after 28 days of osimertinib."""
    p = _fetch("GSE335846", "GSE335846_RNA-seq_for_GEO.txt.gz", raw, manifest)
    df = pd.read_csv(p, sep="\t")
    samples = [c for c in df.columns if c.startswith("Day ")]
    expr = df.groupby("Feature")[samples].mean()
    expr = expr[(expr.max(axis=1) > 2)]
    z = expr.sub(expr.mean(axis=1), axis=0).div(expr.std(axis=1) + 1e-9, axis=0)
    rows = []
    for name, (up, dn) in sigs.items():
        s, nu, nd = score_matrix(z, up, dn)
        for c, v in zip(samples, s):
            day, phase = c.split()[1], c.split()[2]
            rows.append({"signature": name, "sample": c, "day": int(day), "phase": phase, "score": float(v),
                         "genes_up": nu, "genes_down": nd})
    t = pd.DataFrame(rows)
    summary = []
    for name, sub in t.groupby("signature"):
        for ph, g in sub.groupby("phase"):
            d0, d28 = g[g.day == 0].score.mean(), g[g.day == 28].score.mean()
            summary.append({"signature": name, "phase": ph, "day0": d0, "day28": d28, "shift": d28 - d0})
    return {"accession": "GSE335846", "samples": t.to_dict("records"), "summary": summary,
            "note": "Values are the depositors' normalised log expression, averaged over transcripts per gene."}


def _read_mtx(path):
    with gzip.open(path, "rb") as f:
        return sparse.csr_matrix(sio.mmread(f))


def _lognorm(counts):
    lib = np.asarray(counts.sum(axis=0)).ravel()
    x = counts.multiply(1e4 / np.maximum(lib, 1)).tocsr()
    x.data = np.log1p(x.data)
    return x


def _gene_z(x, genes, wanted):
    idx = {g: i for i, g in enumerate(genes)}
    keep = [g for g in dict.fromkeys(wanted) if g in idx]
    if not keep:
        return pd.DataFrame()
    sub = x[[idx[g] for g in keep]].toarray()
    sub = (sub - sub.mean(axis=1, keepdims=True)) / (sub.std(axis=1, keepdims=True) + 1e-9)
    return pd.DataFrame(sub, index=keep)


def _embed(x, genes, seed, n_hvg=2000, n_pc=30):
    import scanpy as sc
    import anndata as ad
    a = ad.AnnData(x.T.tocsr().astype(np.float32))
    a.var_names = pd.Index(genes).astype(str)
    a.var_names_make_unique()
    sc.pp.highly_variable_genes(a, n_top_genes=min(n_hvg, a.n_vars - 1))
    a = a[:, a.var.highly_variable].copy()
    sc.pp.scale(a, max_value=10)
    sc.tl.pca(a, n_comps=n_pc, random_state=seed)
    sc.pp.neighbors(a, random_state=seed)
    sc.tl.umap(a, random_state=seed)
    clusters = {}
    for res in (0.3, 0.6, 1.0):
        sc.tl.leiden(a, resolution=res, random_state=seed, key_added=f"leiden_{res}", flavor="igraph", n_iterations=2)
        clusters[res] = a.obs[f"leiden_{res}"].astype(int).values
    return a.obsm["X_umap"], a.obsm["X_pca"], clusters


def _cycle(x, genes, sets):
    out = {}
    for key, name in (("s", "HALLMARK_E2F_TARGETS"), ("g2m", "HALLMARK_G2M_CHECKPOINT")):
        g = sets.get(name, [])
        z = _gene_z(x, genes, g)
        out[key] = z.mean(axis=0).values if len(z) else np.zeros(x.shape[1])
    cyc = np.maximum(out["s"], out["g2m"])
    gm = GaussianMixture(2, random_state=0).fit(cyc.reshape(-1, 1))
    hi = int(np.argmax(gm.means_.ravel()))
    return out["s"], out["g2m"], gm.predict(cyc.reshape(-1, 1)) == hi


def _bimodality(v, seed=0):
    v = np.asarray(v).reshape(-1, 1)
    b1 = GaussianMixture(1, random_state=seed).fit(v).bic(v)
    b2 = GaussianMixture(2, random_state=seed).fit(v).bic(v)
    return float(b1 - b2)


def _sample_points(n, k, seed):
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(n, size=min(n, k), replace=False))


def pc9_single_cell(raw, manifest, sigs, sets, seed):
    """GSE275199: PC9 and H358 pooled in one 10x run. Hashtag counts are not public, so cell line and
    drug exposure are inferred, and every inference is reported."""
    acc = "GSE275199"
    feats = pd.read_csv(_fetch(acc, "GSE275199_features.tsv.gz", raw, manifest), sep="\t", header=None)
    _fetch(acc, "GSE275199_barcodes.tsv.gz", raw, manifest)
    counts = _read_mtx(_fetch(acc, "GSE275199_matrix.mtx.gz", raw, manifest))
    genes = feats[1].astype(str).values
    lib = np.asarray(counts.sum(axis=0)).ravel()
    ngenes = np.asarray((counts > 0).sum(axis=0)).ravel()
    mt = np.array([g.startswith("MT-") for g in genes])
    mtfrac = np.asarray(counts[mt].sum(axis=0)).ravel() / np.maximum(lib, 1)
    ok = (ngenes >= 500) & (mtfrac < 0.25)
    counts = counts[:, ok]
    x = _lognorm(counts)
    umap, pcs, clusters = _embed(x, genes, seed)
    from sklearn.cluster import KMeans
    km = KMeans(2, n_init=10, random_state=seed).fit(pcs[:, :10]).labels_
    egfr = _gene_z(x, genes, ["EGFR"]).values.ravel() if "EGFR" in set(genes) else np.zeros(x.shape[1])
    pc9_label = int(np.argmax([egfr[km == k].mean() for k in (0, 1)]))
    is_pc9 = km == pc9_label
    mapk = _gene_z(x, genes, MAPK_OUTPUT).mean(axis=0).values
    gm = GaussianMixture(2, random_state=seed).fit(mapk[is_pc9].reshape(-1, 1))
    lo = int(np.argmin(gm.means_.ravel()))
    exposed = np.zeros(len(mapk), bool)
    exposed[is_pc9] = gm.predict(mapk[is_pc9].reshape(-1, 1)) == lo
    s_sc, g2m_sc, cycling = _cycle(x, genes, sets)
    scores, n_used, auc = {}, {}, {}
    for name, (up, dn) in sigs.items():
        up2 = [g for g in up if g not in MAPK_OUTPUT]
        dn2 = [g for g in dn if g not in MAPK_OUTPUT]
        z = _gene_z(x, genes, up2 + dn2)
        s, nu, nd = score_matrix(z, up2, dn2)
        scores[name], n_used[name] = s, (nu, nd)
        auc[name] = _auc(exposed[is_pc9].astype(int), s[is_pc9])
    ps = scores["persister_specific"]
    m_exp = exposed & is_pc9
    sub_cl = clusters[0.6]
    states = []
    for c in np.unique(sub_cl[m_exp]):
        sel = m_exp & (sub_cl == c)
        if sel.sum() < 30:
            continue
        states.append({"cluster": int(c), "cells": int(sel.sum()), "cycling_fraction": float(cycling[sel].mean()),
                       "persister_score": float(ps[sel].mean()), "s_score": float(s_sc[sel].mean()),
                       "g2m_score": float(g2m_sc[sel].mean())})
    from sklearn.metrics import adjusted_rand_score
    stability = float(adjusted_rand_score(clusters[0.3][m_exp], clusters[1.0][m_exp])) if m_exp.sum() > 50 else None
    rho = stats.spearmanr(ps[m_exp], np.maximum(s_sc, g2m_sc)[m_exp])[0] if m_exp.sum() > 50 else None
    keep = _sample_points(len(ps), 4000, seed)
    pts = [[round(float(umap[i, 0]), 3), round(float(umap[i, 1]), 3), int(is_pc9[i]), int(exposed[i]), int(cycling[i]),
            round(float(ps[i]), 3), int(sub_cl[i])] for i in keep]
    return {
        "accession": acc, "cells_after_qc": int(ok.sum()), "cells_total": int(len(ok)),
        "pc9_cells": int(is_pc9.sum()), "other_line_cells": int((~is_pc9).sum()),
        "pc9_assignment": f"k-means on 10 PCs, cluster with higher EGFR expression (mean z {egfr[is_pc9].mean():.2f} vs {egfr[~is_pc9].mean():.2f})",
        "exposed_cells": int(m_exp.sum()), "exposed_fraction_of_pc9": float(m_exp.sum() / max(is_pc9.sum(), 1)),
        "exposure_assignment": "two-component Gaussian mixture on MAPK output genes; low component = drug-exposed",
        "auroc_exposed_vs_unexposed": auc, "signature_genes_used": n_used,
        "cycling_fraction_exposed": float(cycling[m_exp].mean()) if m_exp.any() else None,
        "cycling_fraction_unexposed": float(cycling[is_pc9 & ~exposed].mean()) if (is_pc9 & ~exposed).any() else None,
        "persister_vs_cycling_spearman": None if rho is None else float(rho),
        "bimodality_delta_bic": _bimodality(ps[m_exp]) if m_exp.sum() > 50 else None,
        "states": states, "cluster_stability_ari": stability,
        "points": pts, "point_fields": ["umap1", "umap2", "is_pc9", "drug_exposed", "cycling", "persister_score", "cluster"],
        "caveat": "Hashtag counts were not deposited, so day 4 and day 11 cells cannot be separated and drug exposure is inferred from MAPK output.",
    }


def pdx_single_cell(raw, manifest, sigs, sets, seed):
    """GSE243562: two EGFR-mutant PDXs, untreated and residual after 6 weeks of osimertinib."""
    acc = "GSE243562"
    tar = _fetch(acc, "GSE243562_RAW.tar", raw, manifest)
    meta = {"GSM7790998": ("YU-003", "untreated"), "GSM7790999": ("YU-003", "residual"),
            "GSM7791000": ("YU-006", "untreated"), "GSM7791001": ("YU-006", "residual")}
    mats, feats = {}, {}
    with tarfile.open(tar) as t:
        for m in t.getmembers():
            gsm = m.name.split("_")[0]
            data = t.extractfile(m).read()
            if m.name.endswith("matrix.mtx.gz"):
                mats[gsm] = sparse.csr_matrix(sio.mmread(io.BytesIO(gzip.decompress(data))))
            elif m.name.endswith("features.tsv.gz"):
                feats[gsm] = pd.read_csv(io.BytesIO(gzip.decompress(data)), sep="\t", header=None)
    rows, per_sample = [], []
    blocks, labels = [], []
    genes_h = None
    for gsm, (pdx, cond) in meta.items():
        f, c = feats[gsm], mats[gsm]
        names = f[1].astype(str)
        human = names.str.startswith("hg38_").values
        hc = np.asarray(c[human].sum(axis=0)).ravel()
        tot = np.asarray(c.sum(axis=0)).ravel()
        frac = hc / np.maximum(tot, 1)
        g = names[human].str.replace("hg38_", "", regex=False).values
        ch = c[human]
        ng = np.asarray((ch > 0).sum(axis=0)).ravel()
        ok = (frac > 0.9) & (ng >= 300)
        per_sample.append({"gsm": gsm, "pdx": pdx, "condition": cond, "barcodes": int(len(tot)),
                           "human_tumour_cells": int(ok.sum()), "median_human_fraction": float(np.median(frac))})
        if genes_h is None:
            genes_h = g
        blocks.append(ch[:, ok])
        labels += [(pdx, cond)] * int(ok.sum())
    x = _lognorm(sparse.hstack(blocks).tocsr())
    lab = pd.DataFrame(labels, columns=["pdx", "condition"])
    s_sc, g2m_sc, cycling = _cycle(x, genes_h, sets)
    out = {"accession": acc, "samples": per_sample, "auroc_residual_vs_untreated": {}, "cycling_fraction": {}}
    for name, (up, dn) in sigs.items():
        z = _gene_z(x, genes_h, up + dn)
        s, nu, nd = score_matrix(z, up, dn)
        for pdx, sub in lab.groupby("pdx"):
            i = sub.index.values
            out["auroc_residual_vs_untreated"].setdefault(name, {})[pdx] = _auc((sub.condition == "residual").astype(int).values, s[i])
        if name == "persister_specific":
            lab["persister_score"] = s
    for (pdx, cond), sub in lab.groupby(["pdx", "condition"]):
        out["cycling_fraction"][f"{pdx} {cond}"] = float(cycling[sub.index.values].mean())
    umap, _, clusters = _embed(x, genes_h, seed)
    keep = _sample_points(len(lab), 4000, seed)
    out["points"] = [[round(float(umap[i, 0]), 3), round(float(umap[i, 1]), 3), lab.pdx[i], lab.condition[i],
                      int(cycling[i]), round(float(lab.persister_score[i]), 3)] for i in keep]
    out["point_fields"] = ["umap1", "umap2", "pdx", "condition", "cycling", "persister_score"]
    return out


def run(work, manifest, sigs, sets, seed):
    raw = Path(work) / "raw" / "external"
    raw.mkdir(parents=True, exist_ok=True)
    out, errors = {}, {}
    for key, fn in (("cell_cycle_sorted", lambda: cell_cycle_sorted(raw, manifest, sigs)),
                    ("pc9_single_cell", lambda: pc9_single_cell(raw, manifest, sigs, sets, seed)),
                    ("pdx_single_cell", lambda: pdx_single_cell(raw, manifest, sigs, sets, seed))):
        try:
            out[key] = fn()
        except Exception as err:
            import traceback
            errors[key] = f"{type(err).__name__}: {err}"
            traceback.print_exc()
    out["errors"] = errors
    return out
