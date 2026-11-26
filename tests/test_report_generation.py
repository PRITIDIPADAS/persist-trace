import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BANNED = [r"\bcure[sd]?\b", r"clinically validated", r"first AI to", r"proves? that", r"diagnos(e|is) patients"]


def test_results_json_is_finite(synthetic_run):
    text = json.dumps(synthetic_run["res"])
    assert "NaN" not in text and "Infinity" not in text


def test_hypotheses_have_verdicts(synthetic_run):
    for h in synthetic_run["res"]["hypotheses"]:
        assert h["result"]


def test_no_unsupported_claims_in_public_text():
    files = list((ROOT / "site").glob("*.html")) + list((ROOT / "site" / "js").glob("*.js")) + [ROOT / "README.md"]
    allowed = re.compile(r"(not|never|no|does not|doesn't|cannot)[^.]{0,80}$", re.I)
    for f in files:
        if not f.exists():
            continue
        for line in f.read_text().splitlines():
            for pat in BANNED:
                m = re.search(pat, line, re.I)
                if m and not allowed.search(line[: m.start()]):
                    raise AssertionError(f"{f.name}: unsupported claim '{m.group(0)}' in: {line.strip()[:120]}")


def test_no_em_dashes_in_public_text():
    for f in list((ROOT / "site").rglob("*.html")) + list((ROOT / "docs").glob("*.md")) + [ROOT / "README.md"]:
        if f.exists():
            assert "\u2014" not in f.read_text(), f"em dash in {f}"
