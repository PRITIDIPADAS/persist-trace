import csv
import datetime as dt
import gzip
import hashlib
import io
import re
import time
import urllib.request
from pathlib import Path

import pandas as pd

UA = {"User-Agent": "PersistTrace/1.0 (science fair research; contact via GitHub PRITIDIPADAS)"}


def series_stub(acc):
    return acc[:-3] + "nnn"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def download(url, dest, tries=4):
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    last = None
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=180) as r, open(dest, "wb") as out:
                out.write(r.read())
            if dest.stat().st_size == 0:
                raise IOError("empty response")
            return dest
        except Exception as err:
            last = err
            time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"download failed after {tries} tries: {url}: {last}")


class Manifest:
    fields = ["accession", "role", "url", "path", "bytes", "sha256", "retrieved_utc"]

    def __init__(self, path):
        self.path = Path(path)
        self.rows = []

    def add(self, accession, role, url, path):
        path = Path(path)
        self.rows.append({
            "accession": accession,
            "role": role,
            "url": url,
            "path": str(path),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
            "retrieved_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        })

    def write(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with open(self.path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=self.fields)
            w.writeheader()
            w.writerows(self.rows)
        with open(self.path.with_suffix(".sha256"), "w") as f:
            for r in self.rows:
                f.write(f"{r['sha256']}  {Path(r['path']).name}\n")


def read_text(path):
    raw = Path(path).read_bytes()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    return raw.decode("utf-8", errors="replace")


def parse_series_matrix(path):
    """Return one row per GSM with title, source and all characteristics joined."""
    cols = {}
    char_lines = []
    for line in read_text(path).splitlines():
        if not line.startswith("!Sample_"):
            continue
        key, *vals = next(csv.reader([line], delimiter="\t"))
        vals = [v.strip().strip('"') for v in vals]
        key = key[len("!Sample_"):]
        if key == "characteristics_ch1":
            char_lines.append(vals)
        elif key not in cols:
            cols[key] = vals
    n = len(cols.get("geo_accession", []))
    df = pd.DataFrame({k: v for k, v in cols.items() if len(v) == n})
    for i, vals in enumerate(char_lines):
        name = vals[0].split(":")[0].strip().lower().replace(" ", "_") if vals and ":" in vals[0] else f"char_{i}"
        df[f"ch_{name}"] = [v.split(":", 1)[1].strip() if ":" in v else v for v in vals]
    return df


def read_ncbi_counts(counts_path, annot_path=None):
    counts = pd.read_csv(counts_path, sep="\t", index_col=0)
    counts.index = counts.index.astype(str)
    symbols = None
    if annot_path is not None:
        annot = pd.read_csv(annot_path, sep="\t", index_col=0, low_memory=False)
        annot.index = annot.index.astype(str)
        keep = [c for c in ["Symbol", "Description", "GeneType"] if c in annot.columns]
        symbols = annot[keep]
    return counts.astype(int), symbols


def list_suppl(acc, url_tmpl, dest):
    url = url_tmpl.format(stub=series_stub(acc), acc=acc)
    download(url, dest)
    names = re.findall(r'href="([^"?/][^"]*)"', read_text(dest))
    return url, [n for n in names if acc in n]


def fetch_primary(acc, cfg, raw_dir, manifest):
    raw_dir = Path(raw_dir)
    g = cfg["geo"]
    out = {}
    m_url = g["matrix_url"].format(stub=series_stub(acc), acc=acc)
    out["matrix"] = download(m_url, raw_dir / f"{acc}_series_matrix.txt.gz")
    manifest.add(acc, "series_matrix", m_url, out["matrix"])
    c_url = g["counts_url"].format(acc=acc)
    try:
        out["counts"] = download(c_url, raw_dir / f"{acc}_raw_counts_GRCh38.p13_NCBI.tsv.gz")
        if out["counts"].read_bytes()[:2] != b"\x1f\x8b":
            out.pop("counts").unlink()
            raise RuntimeError("NCBI returned a web page instead of a counts file")
        manifest.add(acc, "ncbi_raw_counts", c_url, out["counts"])
        out["source"] = "NCBI GEO uniform counts (GRCh38.p13)"
        a_url = g["annot_url"]
        try:
            a = download(a_url, raw_dir / "Human.GRCh38.p13.annot.tsv.gz")
            if a.read_bytes()[:2] == b"\x1f\x8b":
                out["annot"] = a
                manifest.add("GRCh38.p13", "ncbi_gene_annotation", a_url, a)
            else:
                a.unlink()
        except RuntimeError:
            pass
    except RuntimeError as err:
        url, names = list_suppl(acc, g["suppl_url"], raw_dir / "suppl_index.html")
        out["suppl_names"] = names
        out["source"] = f"author supplementary files ({err})"
        ranked = sorted(names, key=lambda n: (0 if "count" in n.lower() else 1 if re.search(r"\.(txt|tsv|csv)(\.gz)?$|\.xlsx?$", n, re.I) else 2, n))
        for name in ranked:
            if not re.search(r"\.(txt|tsv|csv)(\.gz)?$|\.xlsx?$", name, re.I):
                continue
            f_url = url + name
            path = download(f_url, raw_dir / name)
            manifest.add(acc, "author_supplementary", f_url, path)
            out.setdefault("suppl_files", []).append(path)
    return out


def read_table_any(path):
    path = Path(path)
    name = path.name.lower()
    if name.endswith((".xlsx", ".xls")):
        return pd.read_excel(path, index_col=0)
    return pd.read_csv(path, sep=None, engine="python", index_col=0, compression="infer")


def _norm(s):
    return re.sub(r"[^a-z0-9]", "", str(s).lower().replace("[rnaseq]", ""))


def read_author_counts(paths, samples):
    """Find the supplementary table whose columns match the sample titles; return genes x GSM counts."""
    want = {_norm(t.split("[")[0]): g for t, g in zip(samples.title, samples.gsm)}
    want.update({_norm(g): g for g in samples.gsm})
    report = []
    for p in paths:
        try:
            df = read_table_any(p)
        except Exception as err:
            report.append(f"{Path(p).name}: unreadable ({err})")
            continue
        report.append(f"{Path(p).name}: shape {df.shape}; columns {list(df.columns[:8])}; index {list(df.index[:5])}")
        mapping = {c: want[_norm(c)] for c in df.columns if _norm(c) in want}
        if len(set(mapping.values())) >= len(samples):
            num = df[list(mapping)].apply(pd.to_numeric, errors="coerce")
            num.columns = [mapping[c] for c in num.columns]
            num = num.loc[:, ~num.columns.duplicated()].dropna(how="all").fillna(0)
            num.index = num.index.astype(str)
            return num.round().clip(lower=0).astype(int), report, Path(p).name
    return None, report, None


HGNC_URLS = [
    "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
    "https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt",
]


def symbol_table(index, raw_dir, manifest):
    """Map gene IDs to symbols. Ensembl IDs are mapped through HGNC; anything else is kept as is."""
    ids = pd.Index(index).astype(str)
    stripped = ids.str.replace(r"\.\d+$", "", regex=True)
    key = "ensembl_gene_id" if stripped.str.startswith("ENSG").mean() > 0.5 else (
        "entrez_id" if stripped.str.fullmatch(r"\d+").mean() > 0.5 else None)
    if key is None:
        return pd.DataFrame({"Symbol": ids.str.split("|").str[-1]}, index=ids)
    if key == "entrez_id":
        url = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
        try:
            path = download(url, Path(raw_dir) / "Homo_sapiens.gene_info.gz")
            manifest.add("NCBI Gene", "gene_symbols", url, path)
            h = pd.read_csv(path, sep="\t", usecols=["GeneID", "Symbol"], dtype=str)
            m = dict(zip(h.GeneID, h.Symbol))
            return pd.DataFrame({"Symbol": [m.get(s, s) for s in stripped]}, index=ids)
        except Exception:
            pass
    for url in HGNC_URLS:
        try:
            path = download(url, Path(raw_dir) / "hgnc_complete_set.txt")
            manifest.add("HGNC", "gene_symbols", url, path)
            h = pd.read_csv(path, sep="\t", usecols=["symbol", key], dtype=str, low_memory=False).dropna()
            h[key] = h[key].str.replace(r"\.0$", "", regex=True)
            m = dict(zip(h[key], h.symbol))
            return pd.DataFrame({"Symbol": [m.get(s, s) for s in stripped]}, index=ids)
        except Exception:
            continue
    return pd.DataFrame({"Symbol": stripped}, index=ids)
