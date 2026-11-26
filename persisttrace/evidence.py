"""Biological Evidence Passports.

Six dimensions are scored separately and never merged into one number. The final label is the
strongest claim the evidence permits, chosen by explicit rules (see CLAIM_RULES)."""
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

DTP_OUTCOMES = ("DTP formation", "DTP maintenance", "DTP vulnerability")

CLAIM_RULES = [
    ("contested", "Contested mechanism", "Published functional studies disagree; the context of each is shown."),
    ("dependency", "Perturbation-supported dependency", "Independent perturbation changed a persister outcome in the cited model."),
    ("negative", "Marker without functional support", "Enriched in persisters, but a perturbation test found no effect on DTP formation."),
    ("regulator", "Candidate regulator", "Inferred regulatory activity changes in persisters; no perturbation evidence in this context."),
    ("marker", "DTP-associated marker", "Reproducibly associated with the persister state beyond acute stress. Association only."),
    ("stress", "Shared stress response", "Change is largely explained by the acute drug response."),
    ("weak", "Insufficient evidence", "The data do not support a persister-specific interpretation."),
]


def load_literature(path):
    lit = yaml.safe_load(Path(path).read_text())
    recs = pd.DataFrame(lit["records"])
    return lit["sources"], recs, lit.get("programs", [])


def _status(ok, partial=False):
    return "supported" if ok else ("partial" if partial else "not supported")


def perturbation_status(recs):
    if recs is None or recs.empty:
        return "none", []
    func = recs[recs.kind.isin(["genetic_perturbation", "functional_and_pharmacological"])]
    if func.empty:
        return "none", []
    dtp = func[func.outcome.str.contains("DTP", na=False)]
    pos = dtp[dtp.effect.isin(["required", "promotes"])]
    neg = dtp[dtp.effect == "none"]
    if len(pos) and len(neg):
        return "contested", list(func.index)
    if len(pos):
        return "supported", list(func.index)
    if len(neg):
        return "negative", list(func.index)
    return "other_context", list(func.index)


def build(symbols, mP, mS, betas, per_line, cls, tf_scores, regulons, lit_recs, sources, samples, cfg,
          symbol_to_id):
    alpha = cfg["thresholds"]["alpha"]
    need = cfg["thresholds"]["reproducible_lines"]
    beta_mean = float(np.mean(list(betas.values())))
    passports = []
    for sym in symbols:
        gid = symbol_to_id.get(sym)
        if gid is None or gid not in mP.index:
            continue
        p, s = mP.loc[gid], mS.loc[gid]
        direction = "up" if p.mu > 0 else "down"
        n_lines = int(p.lines_up if p.mu > 0 else p.lines_down)
        enrich = _status(p.fdr < alpha and n_lines >= need, p.fdr < alpha and n_lines >= 2)
        same = np.sign(s.mu) == np.sign(p.mu)
        a_part = abs(p.mu - s.mu)
        acute_spec = "supported" if (s.fdr < alpha and same and abs(s.mu) >= a_part) else (
            "acute-driven" if a_part > abs(s.mu) and enrich != "not supported" else "uncertain")
        line_cls = {l: cls.loc[gid, l] if gid in cls.index else None for l in per_line}
        rev = cls.loc[gid, "class"] if gid in cls.index else "unresolved"
        rev_status = cls.loc[gid, "status"] if gid in cls.index else "unresolved"
        regs = []
        if regulons is not None and len(regulons):
            for _, e in regulons[regulons.target == sym].iterrows():
                if e.source in tf_scores.index:
                    t = tf_scores.loc[e.source]
                    consistent = np.sign(t.score) * np.sign(e.weight) == np.sign(p.mu)
                    regs.append({"tf": e.source, "weight": float(e.weight), "tf_score": float(t.score),
                                 "tf_fdr": float(t.fdr), "consistent": bool(consistent),
                                 "kind": "inferred (CollecTRI prior, activity from this data)"})
        regs = sorted(regs, key=lambda r: r["tf_fdr"])[:6]
        own_tf = None
        if sym in tf_scores.index:
            t = tf_scores.loc[sym]
            own_tf = {"score": float(t.score), "fdr": float(t.fdr), "n_targets": int(t.n_targets)}
        recs = lit_recs[lit_recs.gene == sym] if len(lit_recs) else lit_recs
        pstat, _ = perturbation_status(recs)
        lit_reg = recs[recs.kind.isin(["chromatin", "regulatory"])] if len(recs) else recs
        reg_ok = any(r["tf_fdr"] < alpha and r["consistent"] for r in regs) or len(lit_reg) > 0 or (
            own_tf is not None and own_tf["fdr"] < alpha)
        contradictions = []
        if p.lines_up and p.lines_down:
            contradictions.append(f"Direction differs across cell lines ({int(p.lines_up)} up, {int(p.lines_down)} down).")
        if p.I2 > 0.75:
            contradictions.append(f"High between-line heterogeneity (I2 = {p.I2:.2f}).")
        if pstat == "contested":
            contradictions.append("Published perturbation studies disagree.")
        if acute_spec == "acute-driven":
            contradictions.append("Most of the persister change is already present after 24 h of drug.")
        if rev_status == "context_dependent":
            contradictions.append("Withdrawal behaviour differs between cell lines.")

        if pstat == "contested":
            key = "contested"
        elif pstat == "supported" and enrich != "not supported":
            key = "dependency"
        elif pstat == "negative" and enrich != "not supported":
            key = "negative"
        elif own_tf is not None and own_tf["fdr"] < alpha and enrich != "not supported":
            key = "regulator"
        elif enrich == "supported" and acute_spec == "supported":
            key = "marker"
        elif enrich != "not supported" and acute_spec == "acute-driven":
            key = "stress"
        else:
            key = "weak"
        label, meaning = next((l, m) for k, l, m in CLAIM_RULES if k == key)

        lines = {}
        for l, df in per_line.items():
            if gid in df.index:
                r = df.loc[gid]
                lines[l] = {"A": _r(r.A), "P": _r(r.P), "S": _r(r.S), "z_S": _r(r.z_S), "class": line_cls.get(l)}
        gsms = samples[samples.condition.isin(["control", "dtp", "acute"])].gsm.tolist()
        passports.append({
            "symbol": sym, "gene_id": gid, "direction": direction,
            "dimensions": {
                "persistence_enrichment": {"status": enrich, "meta_lfc": _r(p.mu), "fdr": _r(p.fdr, 6),
                                           "lines_agreeing": n_lines},
                "acute_stress_specificity": {"status": acute_spec, "persister_specific_lfc": _r(s.mu),
                                             "fdr": _r(s.fdr, 6), "acute_explained_lfc": _r(p.mu - s.mu)},
                "reversibility": {"status": rev_status, "class": rev, "per_line": line_cls},
                "regulatory_support": {"status": "supported" if reg_ok else "not found", "regulators": regs,
                                       "own_tf_activity": own_tf,
                                       "literature": lit_reg.to_dict("records") if len(lit_reg) else []},
                "perturbation_support": {"status": pstat,
                                         "records": recs[~recs.kind.isin(["chromatin", "regulatory"])].to_dict("records") if len(recs) else []},
                "reproducibility": {"status": _status(n_lines >= need, n_lines >= 2), "lines": n_lines,
                                    "of": len(per_line), "I2": _r(p.I2)},
            },
            "per_line": lines,
            "claim": {"key": key, "label": label, "meaning": meaning},
            "contradictions": contradictions,
            "trace": {"contrasts": ["dtp_vs_control", "acute_vs_control", "persister-specific residual S = P - beta*A"],
                      "beta_mean": _r(beta_mean), "test": "DESeq2 Wald per line, DerSimonian-Laird random effects across lines, BH FDR",
                      "samples": gsms, "dataset": cfg["primary_dataset"], "publication": sources.get("criscione2022", {})},
        })
    return passports


def _r(x, n=3):
    try:
        x = float(x)
    except (TypeError, ValueError):
        return None
    return round(x, n) if np.isfinite(x) else None
