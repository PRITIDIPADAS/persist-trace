import re

import pandas as pd

LINES = {"pc9": "PC9", "h1975": "H1975", "ncih1975": "H1975", "hcc827": "HCC827", "hcc2935": "HCC2935"}


class MetadataError(ValueError):
    pass


def _text(row):
    parts = [str(v) for k, v in row.items() if k in ("title", "source_name_ch1") or k.startswith("ch_")]
    return " | ".join(parts)


def cell_line(text):
    flat = re.sub(r"[^a-z0-9]", "", text.lower())
    hits = {v for k, v in LINES.items() if k in flat}
    return hits.pop() if len(hits) == 1 else None


def _days(text):
    t = text.lower()
    m = re.search(r"(\d+(?:\.\d+)?)\s*-?\s*(d\b|day|days)", t)
    if m:
        return float(m.group(1))
    m = re.search(r"(\d+(?:\.\d+)?)\s*-?\s*(h\b|hr|hrs|hour|hours)", t)
    if m:
        return float(m.group(1)) / 24
    return None


def condition(text):
    t = re.sub(r"[_\-/]", " ", text.lower())
    wash = re.search(r"wash\s*out|washed|\bwo\b|\bwo\d|drug holiday|release|recover|regrow", t)
    if wash:
        if re.search(r"short|(?<!\d)24\s*h|(?<!\d)1\s*d(ay)?\b|\bst\b", t):
            return "wo_short"
        if re.search(r"long|\blt\b|late", t):
            return "wo_long"
        d = _days(t[wash.start():]) or _days(t)
        if d is not None:
            return "wo_short" if d <= 1.5 else "wo_long"
        return None
    if re.search(r"\bdtp\b|persist|21\s*d|21\s*days|chronic|long term treat", t):
        return "dtp"
    if re.search(r"dmso|vehicle|control|untreated|parental|naive", t):
        return "control"
    if re.search(r"acute|24\s*h|osimertinib|osi\b|treated", t):
        return "acute"
    return None


def build_sample_table(matrix_df, design):
    rows = []
    for _, r in matrix_df.iterrows():
        text = _text(r)
        cond = condition(text)
        rows.append({
            "gsm": r["geo_accession"],
            "title": r.get("title", ""),
            "cell_line": cell_line(text),
            "condition": cond,
            "washout_days": _days(text) if cond in ("wo_short", "wo_long") else None,
            "raw_text": text,
        })
    df = pd.DataFrame(rows)
    df["replicate"] = df.groupby(["cell_line", "condition"]).cumcount() + 1
    validate(df, design)
    return df


def validate(df, design):
    problems = []
    bad = df[df.cell_line.isna() | df.condition.isna()]
    if len(bad):
        problems.append(f"{len(bad)} samples could not be labelled: {list(bad.gsm)[:8]}")
    counts = df.dropna().groupby(["cell_line", "condition"]).size()
    exp = design["expected_replicates"]
    for line in design["cell_lines"]:
        for cond in design["conditions"]:
            n = int(counts.get((line, cond), 0))
            if n != exp:
                problems.append(f"{line}/{cond}: {n} samples, expected {exp}")
    if problems:
        raise MetadataError("; ".join(problems))


def align_counts(counts, samples):
    missing = [g for g in samples.gsm if g not in counts.columns]
    if missing:
        raise MetadataError(f"count matrix lacks samples: {missing[:8]}")
    return counts[list(samples.gsm)]
