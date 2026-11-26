import argparse
import datetime as dt
import json
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from . import (CONDITION_LABELS, __version__, benchmark, contrasts, decompose, evidence, geo, graph,
               programs, reversibility, robustness, samples as smp, synthetic)

ROOT = Path(__file__).resolve().parents[1]


def log(msg):
    print(f"[{dt.datetime.now():%H:%M:%S}] {msg}", flush=True)


def load(cfg, source, work):
    acc = cfg["primary_dataset"]
    if source == "synthetic":
        counts, annot, sample_df, _ = synthetic.make()
        return counts, annot, sample_df, pd.DataFrame(), {"source": "synthetic"}
    man = geo.Manifest(work / "manifests" / "file_manifest.csv")
    got = geo.fetch_primary(acc, cfg, work / "raw", man)
    man.write()
    matrix = geo.parse_series_matrix(got["matrix"])
    matrix.to_csv(work / "manifests" / "series_matrix_samples.csv", index=False)
    try:
        sample_df = smp.build_sample_table(matrix, cfg["design"])
    except smp.MetadataError as err:
        dbg = pd.DataFrame([{"gsm": r.geo_accession, "text": smp._text(r)} for _, r in matrix.iterrows()])
        dbg.to_csv(work / "manifests" / "metadata_debug.csv", index=False)
        man.write()
        raise SystemExit(f"metadata parsing failed: {err}")
    if "counts" in got:
        counts, annot = geo.read_ncbi_counts(got["counts"], got.get("annot"))
        got["counts_file"] = Path(got["counts"]).name
    else:
        counts, report, fname = geo.read_author_counts(got.get("suppl_files", []), sample_df)
        (work / "manifests" / "suppl_report.txt").write_text(
            "supplementary files: " + ", ".join(got.get("suppl_names", [])) + "\n" + "\n".join(report))
        if counts is None:
            man.write()
            raise SystemExit("no supplementary table matched the 60 sample titles; see suppl_report.txt")
        got["counts_file"] = fname
        annot = geo.symbol_table(counts.index, work / "raw", man)
    man.write()
    counts = smp.align_counts(counts, sample_df)
    log(f"counts: {counts.shape[0]} genes x {counts.shape[1]} samples from {got['counts_file']}")
    return counts, annot, sample_df, pd.DataFrame(man.rows), got


def qc(counts, sample_df):
    lib = counts.sum(axis=0)
    detected = (counts > 0).sum(axis=0)
    lc = contrasts.log_cpm(counts)
    rows = []
    for line, sub in sample_df.groupby("cell_line"):
        cols = list(sub.gsm)
        corr = lc[cols].corr(method="spearman")
        mean_r = (corr.sum() - 1) / (len(cols) - 1)
        zr = (mean_r - mean_r.mean()) / (mean_r.std() + 1e-9)
        for g in cols:
            rows.append({"gsm": g, "library_size": int(lib[g]), "genes_detected": int(detected[g]),
                         "mean_within_line_spearman": float(mean_r[g]), "outlier_flag": bool(zr[g] < -3)})
    return pd.DataFrame(rows)


def pca(z, n=3):
    x = z.T.values
    x = x - x.mean(axis=0)
    u, s, vt = np.linalg.svd(x, full_matrices=False)
    var = s ** 2 / (s ** 2).sum()
    return u[:, :n] * s[:n], var[:n]


def condition_means(lc, sample_df):
    out = {}
    for line, sub in sample_df.groupby("cell_line"):
        cols = list(sub.gsm)
        z = lc[cols].sub(lc[cols].loc[:, sub[sub.condition == "control"].gsm].mean(axis=1), axis=0)
        out[line] = pd.DataFrame({c: z[list(sub[sub.condition == c].gsm)].mean(axis=1)
                                  for c in ["control", "acute", "dtp", "wo_short", "wo_long"]})
    return out


def rank_percentile(meta_df, gid):
    if gid not in meta_df.index:
        return None
    r = meta_df.z.abs().rank(ascending=False, pct=True)
    return float(r.loc[gid])


def main(argv=None):
    ap = argparse.ArgumentParser(description="PersistTrace analysis pipeline")
    ap.add_argument("--source", choices=["geo", "synthetic"], default="geo")
    ap.add_argument("--work", default="build")
    ap.add_argument("--site", default="site")
    ap.add_argument("--skip-networks", action="store_true")
    ap.add_argument("--skip-media", action="store_true")
    args = ap.parse_args(argv)

    cfg = yaml.safe_load((ROOT / "config" / "analysis.yaml").read_text())
    np.random.seed(cfg["seed"])
    work = Path(args.work)
    for d in ["raw", "manifests", "tables", "figures"]:
        (work / d).mkdir(parents=True, exist_ok=True)
    started = dt.datetime.now(dt.timezone.utc)

    log(f"loading {args.source} data")
    counts, annot, sample_df, manifest, got = load(cfg, args.source, work)
    sample_df.drop(columns=["raw_text"], errors="ignore").to_csv(work / "manifests" / "sample_manifest.csv", index=False)
    sym = annot["Symbol"].reindex(counts.index).fillna(pd.Series(counts.index, index=counts.index))
    symbol_to_id = {s: g for g, s in sym.items()}

    log("quality control")
    qc_df = qc(counts, sample_df)
    qc_df.to_csv(work / "tables" / "qc.csv", index=False)

    log("DESeq2 contrasts per cell line")
    results, kept = contrasts.run_all(counts, sample_df, cfg)
    for line, res in results.items():
        for name, r in res.items():
            r.assign(symbol=sym.reindex(r.index).values).to_csv(work / "tables" / f"de_{line}_{name}.csv.gz")

    log("signal decomposition")
    lines = sorted(results)
    betas, per_line, mS, mP = decompose.run(results, lines)
    mS.assign(symbol=sym.reindex(mS.index).values).to_csv(work / "tables" / "meta_persister_specific.csv.gz")
    mP.assign(symbol=sym.reindex(mP.index).values).to_csv(work / "tables" / "meta_treated_vs_control.csv.gz")

    log("reversibility and memory")
    th = cfg["thresholds"]
    per_cls = {l: reversibility.classify_line(results[l], th["alpha"], th["lfc"]) for l in lines}
    cls = reversibility.consensus(per_cls, th["reproducible_lines"])
    cls.assign(symbol=sym.reindex(cls.index).values).to_csv(work / "tables" / "reversibility_classes.csv.gz")
    mem = {l: reversibility.memory(results[l]) for l in lines}
    mem_summary = reversibility.summarize_memory(mem, th["alpha"], th["lfc"])
    wo_days = sample_df[sample_df.condition == "wo_long"].groupby("cell_line").washout_days.median().to_dict()
    decay = {}
    for l in lines:
        d = wo_days.get(l)
        sel = mem[l][(mem[l].dtp_padj < th["alpha"]) & (mem[l].dtp_lfc.abs() >= th["lfc"])]
        if d is None or not np.isfinite(d) or d <= 1:
            decay[l] = {"long_washout_days": None, "status": "washout duration not in GEO metadata"}
            continue
        k, hl = reversibility.fit_decay(sel.m_wo_short, sel.m_wo_long, 1.0, d)
        hl = pd.Series(hl, index=sel.index).dropna()
        decay[l] = {"long_washout_days": float(d), "status": "fitted", "n_genes": int(len(hl)),
                    "median_half_life_days": float(hl.median()) if len(hl) else None,
                    "iqr": [float(hl.quantile(.25)), float(hl.quantile(.75))] if len(hl) else None}
        sel.assign(half_life_days=hl).to_csv(work / "tables" / f"memory_{l}.csv.gz")
    for l in lines:
        mem[l].to_csv(work / "tables" / f"memory_all_{l}.csv.gz")
    oct_rows = []
    for l in lines:
        d = decay[l].get("long_washout_days")
        sel = mem[l][(mem[l].dtp_padj < th["alpha"]) & (mem[l].dtp_lfc.abs() >= th["lfc"])]
        if d:
            oct_rows.append(pd.DataFrame({"line": l, "gene_id": sel.index, "m_short": sel.m_wo_short,
                                          "m_long": sel.m_wo_long, "t_long": d}))
    (work / "octave").mkdir(exist_ok=True)
    if oct_rows:
        pd.concat(oct_rows).replace([np.inf, -np.inf], np.nan).dropna().to_csv(work / "octave" / "memory_input.csv", index=False)

    log("programs and regulators")
    sets, set_origin, tf_scores, regulons, net_origin = None, "not loaded", pd.DataFrame(), None, "not loaded"
    prog_S = prog_P = pd.DataFrame()
    if not args.skip_networks:
        try:
            sets, set_origin = programs.load_hallmark(work / "raw" / "hallmark.csv")
        except Exception as err:
            set_origin = f"failed: {err}"
        try:
            regulons, net_origin = programs.load_collectri(work / "raw" / "collectri.csv")
        except Exception as err:
            net_origin = f"failed: {err}"
    zS_sym = pd.Series(mS.z.values, index=sym.reindex(mS.index).values).groupby(level=0).mean()
    zP_sym = pd.Series(mP.z.values, index=sym.reindex(mP.index).values).groupby(level=0).mean()
    if sets is not None:
        prog_S = programs.program_scores(zS_sym, sets)
        prog_P = programs.program_scores(zP_sym, sets)
    if regulons is not None:
        tf_scores = programs.ulm(zS_sym, regulons, min_n=10)
        tf_scores.to_csv(work / "tables" / "tf_activity_persister_specific.csv")
    lc = contrasts.log_cpm(counts)
    lc_sym = lc.groupby(sym.reindex(lc.index).values).mean()
    prog_time = {}
    if sets is not None:
        focus = list(prog_S.sort_values("p").head(12).index)
        for p in ["HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_DNA_REPAIR",
                  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                  "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
                  "HALLMARK_P53_PATHWAY", "HALLMARK_MYC_TARGETS_V1"]:
            if p not in focus:
                focus.append(p)
        for line, sub in sample_df.groupby("cell_line"):
            m = programs.sample_program_matrix(lc_sym[list(sub.gsm)], sets, focus)
            prog_time[line] = {p: {c: float(m.loc[p, list(sub[sub.condition == c].gsm)].mean())
                                   for c in cfg["design"]["conditions"]} for p in m.index}

    log("evidence passports")
    sources, lit_recs, lit_programs = evidence.load_literature(ROOT / "evidence" / "literature.yaml")
    up_ids, dn_ids = decompose.signature(mS, 60, anchor=mP)
    cand_ids = up_ids + dn_ids
    cand = [sym[g] for g in cand_ids]
    cand += [g for g in lit_recs.gene.unique() if g in symbol_to_id and g not in cand]
    if len(tf_scores):
        cand += [t for t in tf_scores[tf_scores.fdr < 0.05].sort_values("p").head(10).index
                 if t in symbol_to_id and t not in cand]
    passports = evidence.build(cand, mP, mS, betas, per_line, cls, tf_scores if len(tf_scores) else pd.DataFrame(columns=["score", "fdr", "n_targets"]),
                               regulons, lit_recs, sources, sample_df, cfg, symbol_to_id)
    set_map = {s: set(g.target) for s, g in sets.groupby("source")} if sets is not None else {}
    top_prog = list(prog_S.sort_values("p").head(8).index) if len(prog_S) else []
    egraph = graph.build(passports, set_map, top_prog, sources)

    log("leave-one-cell-line-out benchmark")
    folds, preds, bench = benchmark.run(lc, sample_df, per_line, cfg)
    folds.to_csv(work / "tables" / "benchmark_folds.csv", index=False)
    preds.to_csv(work / "tables" / "benchmark_predictions.csv", index=False)
    bench.to_csv(work / "tables" / "benchmark_summary.csv", index=False)

    log("robustness")
    grid = robustness.threshold_grid(results, cfg, cls)
    lolo = robustness.leave_one_line_out(per_line, cfg["benchmark"]["signature_size"])
    beta_sens = robustness.beta_sensitivity(results)

    log("embedding")
    genes_common = sorted(set.intersection(*[set(per_line[l].index) for l in lines]))
    zc = benchmark.centre_by_line(lc.loc[genes_common], sample_df)
    top = zc.var(axis=1).sort_values(ascending=False).head(2000).index
    coords, var = pca(zc.loc[top])
    raw_coords, raw_var = pca(lc.loc[genes_common].loc[lc.loc[genes_common].var(axis=1).sort_values(ascending=False).head(2000).index])

    log("hypotheses and gates")
    fold_mean = folds.mean(numeric_only=True)
    h = []
    a_n, a_d = fold_mean.get("naive_treated_vs_control|auroc_dtp_vs_acute"), fold_mean.get("persister_specific|auroc_dtp_vs_acute")
    h.append({"id": "H1", "text": "Removing the acute stress component makes persister signatures separate persisters from acute treatment better in an unseen cell line.",
              "metric": "mean held-out AUROC, DTP vs acute", "naive": a_n, "decomposed": a_d,
              "result": _verdict(a_d, a_n)})
    r_n, r_d = fold_mean.get("naive_treated_vs_control|replicated_fraction"), fold_mean.get("persister_specific|replicated_fraction")
    h.append({"id": "H1b", "text": "Persister-specific signature genes replicate in the held-out cell line at least as often as naive treated-vs-control genes.",
              "metric": "fraction of signature genes significant with the same sign in the held-out line", "naive": r_n, "decomposed": r_d,
              "result": _verdict(r_d, r_n)})
    dist = cls["class"][cls.status == "reproducible"].value_counts().to_dict()
    n_rev = sum(v for k, v in dist.items() if k in ("dtp_reversed_by_24h", "dtp_persists_24h_reversed_late"))
    n_per = dist.get("dtp_persists_after_regrowth", 0)
    h.append({"id": "H2", "text": "Persister programmes differ in how they behave after drug withdrawal.",
              "metric": "reproducible genes per withdrawal class", "classes": dist,
              "result": "supported" if n_rev > 0 and n_per > 0 else "not supported"})
    tf_sig = int((tf_scores.fdr < 0.05).sum()) if len(tf_scores) else None
    h.append({"id": "H3", "text": "Regulatory inference adds information beyond expression alone.",
              "metric": "TFs with significant persister-specific activity (FDR < 0.05)", "value": tf_sig,
              "result": "not evaluated" if tf_sig is None else ("partially supported" if tf_sig > 0 else "not supported"),
              "note": "Activity is inferred from target expression, so this is supporting, not independent, evidence."})
    marker_neg = [p["symbol"] for p in passports if p["claim"]["key"] in ("negative", "contested")]
    h.append({"id": "H4", "text": "Some strong persister markers lack evidence that they are needed for persister formation.",
              "metric": "passports labelled marker-without-support or contested", "genes": marker_neg,
              "result": "supported" if marker_neg else "not supported in this run"})
    h.append({"id": "H5", "text": "Different persister subpopulations coexist.",
              "result": "not testable with bulk RNA-seq",
              "note": "Needs single-cell data (GSE243562, GSE275199). Planned, not run."})
    lit_ids = {g: symbol_to_id[g] for g in lit_recs.gene.unique() if g in symbol_to_id}
    ranks = [{"gene": g, "naive_percentile": rank_percentile(mP, i), "persister_specific_percentile": rank_percentile(mS, i)}
             for g, i in lit_ids.items()]
    h.append({"id": "H6", "text": "Evidence-aware ranking places literature-supported genes higher than expression-only ranking.",
              "metric": "rank percentile (lower is better) of genes with published perturbation records", "ranks": ranks,
              "result": "exploratory", "note": "Tiny, literature-biased reference set. Not a validation."})

    best = bench.iloc[0].to_dict() if len(bench) else {}
    gates = {
        "A_data_integrity": "passed" if args.source == "geo" and len(manifest) >= 2 else ("synthetic run" if args.source == "synthetic" else "failed"),
        "B_analytical_correctness": "passed" if sample_df.groupby(["cell_line", "condition"]).size().eq(3).all() else "failed",
        "C_scientific_interpretation": "passed" if all(e.get("method") for e in egraph["edges"]) else "failed",
        "D_ml_validity": "passed" if len(bench) and set(preds.held_out) == set(lines) else "failed",
        "E_robustness": "passed" if len(grid) and len(lolo) else "not yet evaluated",
        "F_software_reproducibility": "evaluated in CI (see tests)",
        "G_judge_facing": "evaluated by rehearsal, not by code",
        "H_release_integrity": "passed" if args.source == "geo" else "synthetic run",
    }

    log("export")
    cond_means = condition_means(lc, sample_df)
    gene_rows = {}
    for gid in genes_common:
        s = sym.get(gid, gid)
        if s in gene_rows:
            continue
        gene_rows[s] = {
            "id": gid,
            "S": _r(mS.z.get(gid)), "P": _r(mP.z.get(gid)),
            "lfc": _r(mP.mu.get(gid)), "slfc": _r(mS.mu.get(gid)),
            "c": cls["class"].get(gid), "st": cls["status"].get(gid),
            "t": {l: [_r(v, 2) for v in cond_means[l].loc[gid].values] for l in lines},
        }
    sample_points = []
    for i, (_, r) in enumerate(sample_df.iterrows()):
        q = qc_df.set_index("gsm").loc[r.gsm]
        sample_points.append({"gsm": r.gsm, "title": r.title, "line": r.cell_line, "condition": r.condition,
                              "rep": int(r.replicate), "pc": [_r(v) for v in coords[i]],
                              "pc_raw": [_r(v) for v in raw_coords[i]],
                              "lib": int(q.library_size), "outlier": bool(q.outlier_flag)})
    out = {
        "meta": {"project": "PersistTrace", "version": __version__, "source": args.source,
                 "synthetic": args.source == "synthetic", "dataset": cfg["primary_dataset"],
                 "generated_utc": started.isoformat(timespec="seconds"),
                 "python": platform.python_version(), "n_samples": int(len(sample_df)),
                 "n_genes_tested": int(len(genes_common)), "cell_lines": lines,
                 "conditions": CONDITION_LABELS, "pca_variance": [_r(v) for v in var],
                 "pca_raw_variance": [_r(v) for v in raw_var],
                 "networks": {"hallmark": set_origin, "collectri": net_origin},
                 "counts_source": got.get("source"), "counts_file": got.get("counts_file")},
        "manifest": manifest.to_dict("records") if len(manifest) else [],
        "samples": sample_points,
        "betas": {l: _r(b) for l, b in betas.items()},
        "beta_sensitivity": beta_sens.to_dict("records"),
        "decomposition_scatter": _scatter(per_line, sym),
        "reversibility": {"distribution": dist,
                          "status_counts": cls.status.value_counts().to_dict(),
                          "per_line": {l: per_cls[l].value_counts().to_dict() for l in lines}},
        "memory": {"summary": mem_summary, "decay": decay},
        "programs": {"persister_specific": _top_table(prog_S, 20), "treated_vs_control": _top_table(prog_P, 20),
                     "time_course": prog_time, "literature": lit_programs},
        "tf_activity": _top_table(tf_scores, 25) if len(tf_scores) else [],
        "passports": passports,
        "graph": egraph,
        "benchmark": {"summary": bench.to_dict("records"), "folds": folds.to_dict("records")},
        "robustness": {"threshold_grid": grid.to_dict("records"), "leave_one_line_out": lolo.to_dict("records"),
                       "labels": {"threshold_agreement": robustness.label(grid.agreement_with_reference.min()),
                                  "signature_lolo": robustness.label(lolo.jaccard_top.min(), 0.5, 0.3)}},
        "hypotheses": h, "gates": gates, "qc": qc_df.to_dict("records"),
        "sources": sources,
    }
    site = Path(args.site) / "data"
    site.mkdir(parents=True, exist_ok=True)
    tag = "" if args.source == "geo" else "_synthetic"
    (site / f"results{tag}.json").write_text(json.dumps(_clean(out), separators=(",", ":")))
    (site / f"genes{tag}.json").write_text(json.dumps(_clean(gene_rows), separators=(",", ":")))
    (work / "results_summary.json").write_text(json.dumps(_clean({k: out[k] for k in ["meta", "betas", "hypotheses", "gates", "reversibility", "memory"]}), indent=1))
    if not args.skip_media:
        log("figures and video")
        from . import figures
        figures.render_all(out, gene_rows, work / "figures", Path(args.site) / "media", synthetic=args.source == "synthetic")
    log("done")
    return 0


def _verdict(new, old, margin=0.02):
    if new is None or old is None or not np.isfinite(new) or not np.isfinite(old):
        return "not evaluated"
    if new > old + margin:
        return "supported"
    if new < old - margin:
        return "contradicted"
    return "no clear difference"


def _scatter(per_line, sym, n=1500):
    out = {}
    rng = np.random.default_rng(0)
    for l, df in per_line.items():
        d = df[np.isfinite(df.P) & np.isfinite(df.A)]
        idx = rng.choice(len(d), size=min(n, len(d)), replace=False)
        sub = d.iloc[idx]
        out[l] = [[_r(a, 2), _r(p, 2), _r(s, 2), sym.get(g, g)] for g, a, p, s in zip(sub.index, sub.A, sub.P, sub.S)]
    return out


def _top_table(df, n):
    if df is None or not len(df):
        return []
    t = df.sort_values("p").head(n).reset_index()
    return [{k: (_r(v, 4) if isinstance(v, (float, np.floating)) else v) for k, v in r.items()} for r in t.to_dict("records")]


def _r(x, n=3):
    try:
        x = float(x)
    except (TypeError, ValueError):
        return None
    return round(x, n) if np.isfinite(x) else None


def _clean(o):
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else round(float(o), 5)
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, (pd.Timestamp,)):
        return o.isoformat()
    if o is pd.NA or (isinstance(o, float) and np.isnan(o)):
        return None
    return o


if __name__ == "__main__":
    sys.exit(main())
