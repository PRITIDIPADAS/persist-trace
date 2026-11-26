"""Molecular Escape Graph. Every edge carries a type, a method and a source."""

EDGE_TYPES = {
    "observed_expression": "Measured in GSE193258 (DESeq2, this pipeline)",
    "inferred_regulatory": "CollecTRI prior network with activity inferred from this data",
    "pathway_membership": "MSigDB hallmark gene set membership",
    "published_perturbation": "Perturbation experiment reported in the cited paper",
    "published_regulatory": "Regulatory relationship reported in the cited paper",
    "published_chromatin": "Chromatin result reported in the cited paper",
}


def build(passports, program_sets, top_programs, sources):
    nodes, edges = {}, []

    def node(i, kind, **kw):
        nodes.setdefault(i, {"id": i, "kind": kind, **kw})

    node("state:dtp", "state", label="Persister state")
    for p in passports:
        g = f"gene:{p['symbol']}"
        node(g, "gene", label=p["symbol"], claim=p["claim"]["key"], direction=p["direction"])
        d = p["dimensions"]["persistence_enrichment"]
        edges.append({"source": "state:dtp", "target": g, "type": "observed_expression",
                      "sign": 1 if p["direction"] == "up" else -1, "value": d["meta_lfc"], "fdr": d["fdr"],
                      "evidence": "GSE193258", "method": EDGE_TYPES["observed_expression"]})
        for r in p["dimensions"]["regulatory_support"]["regulators"][:3]:
            if r["tf_fdr"] < 0.05:
                t = f"tf:{r['tf']}"
                node(t, "tf", label=r["tf"], score=r["tf_score"])
                edges.append({"source": t, "target": g, "type": "inferred_regulatory",
                              "sign": 1 if r["weight"] > 0 else -1, "value": r["tf_score"], "fdr": r["tf_fdr"],
                              "evidence": "CollecTRI", "method": EDGE_TYPES["inferred_regulatory"]})
        for rec in p["dimensions"]["perturbation_support"]["records"]:
            s = f"src:{rec['source']}"
            node(s, "source", label=_cite(sources.get(rec["source"], {})), doi=sources.get(rec["source"], {}).get("doi"))
            edges.append({"source": s, "target": g, "type": "published_perturbation", "effect": rec.get("effect"),
                          "outcome": rec.get("outcome"), "evidence": rec["source"],
                          "method": EDGE_TYPES["published_perturbation"]})
        for rec in p["dimensions"]["regulatory_support"]["literature"]:
            kind = "published_chromatin" if rec["kind"] == "chromatin" else "published_regulatory"
            s = f"src:{rec['source']}"
            node(s, "source", label=_cite(sources.get(rec["source"], {})), doi=sources.get(rec["source"], {}).get("doi"))
            edges.append({"source": s, "target": g, "type": kind, "effect": rec.get("effect"),
                          "evidence": rec["source"], "method": EDGE_TYPES[kind]})
        for prog in top_programs:
            members = program_sets.get(prog, set())
            if p["symbol"] in members:
                node(f"prog:{prog}", "program", label=prog.replace("HALLMARK_", "").replace("_", " ").title())
                edges.append({"source": g, "target": f"prog:{prog}", "type": "pathway_membership",
                              "evidence": "MSigDB", "method": EDGE_TYPES["pathway_membership"]})
    return {"nodes": list(nodes.values()), "edges": edges, "edge_types": EDGE_TYPES}


def _cite(src):
    if not src:
        return "unknown source"
    return f"{src.get('journal', '')} {src.get('year', '')}".strip()
