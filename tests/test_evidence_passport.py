import pandas as pd

from persisttrace import evidence


def test_contested_when_studies_disagree():
    recs = pd.DataFrame([
        {"gene": "X", "kind": "genetic_perturbation", "outcome": "DTP formation", "effect": "none"},
        {"gene": "X", "kind": "functional_and_pharmacological", "outcome": "DTP maintenance", "effect": "required"},
    ])
    assert evidence.perturbation_status(recs)[0] == "contested"


def test_other_context_is_not_promoted():
    recs = pd.DataFrame([{"gene": "X", "kind": "genetic_perturbation",
                          "outcome": "apoptosis under acute treatment", "effect": "sensitizes"}])
    assert evidence.perturbation_status(recs)[0] == "other_context"


def test_no_dependency_without_perturbation(synthetic_run):
    for p in synthetic_run["res"]["passports"]:
        if p["claim"]["key"] == "dependency":
            assert p["dimensions"]["perturbation_support"]["status"] == "supported"
            assert p["dimensions"]["perturbation_support"]["records"]


def test_six_dimensions_always_present(synthetic_run):
    want = {"persistence_enrichment", "acute_stress_specificity", "reversibility",
            "regulatory_support", "perturbation_support", "reproducibility"}
    for p in synthetic_run["res"]["passports"]:
        assert set(p["dimensions"]) == want
        assert p["trace"]["samples"] and p["trace"]["dataset"]


def test_trop2_is_contested(synthetic_run):
    p = next(p for p in synthetic_run["res"]["passports"] if p["symbol"] == "TACSTD2")
    assert p["claim"]["key"] == "contested"


def test_literature_sources_resolve():
    sources, recs, progs = evidence.load_literature("evidence/literature.yaml")
    for s in list(recs.source) + [p["source"] for p in progs]:
        assert s in sources and sources[s]["doi"] and sources[s]["pmid"]
