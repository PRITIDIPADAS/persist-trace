def test_every_edge_has_source_and_method(synthetic_run):
    g = synthetic_run["res"]["graph"]
    ids = {n["id"] for n in g["nodes"]}
    for e in g["edges"]:
        assert e["source"] in ids and e["target"] in ids
        assert e.get("method") and e.get("evidence")
        assert e["type"] in g["edge_types"]


def test_inferred_edges_never_labelled_published(synthetic_run):
    for e in synthetic_run["res"]["graph"]["edges"]:
        if e["type"] == "inferred_regulatory":
            assert e["evidence"] == "CollecTRI"
