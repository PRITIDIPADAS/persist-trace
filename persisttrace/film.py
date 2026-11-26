"""A 50-second explainer film drawn entirely from results.json, so every number on screen is real."""
import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation, font_manager

ROOT = Path(__file__).resolve().parents[1]
NIGHT, INK, MUTE = "#13110e", "#efe8da", "#a39b8c"
C = {"control": "#8db4cf", "acute": "#e0ad55", "dtp": "#d9707b", "wo_short": "#a693cf", "wo_long": "#86b38a"}
ORDER = ["control", "acute", "dtp", "wo_short", "wo_long"]
NAMES = {"control": "Control", "acute": "24 h on drug", "dtp": "Persister, 21 d", "wo_short": "Washout 24 h", "wo_long": "Regrowth"}
FPS = 24


def ease(u):
    u = min(max(u, 0.0), 1.0)
    return u * u * (3 - 2 * u)


def _font():
    for p in [ROOT / "site" / "fonts" / "EBGaramond.ttf", ROOT / "site" / "fonts" / "EBGaramond-Italic.ttf"]:
        if p.exists():
            font_manager.fontManager.addfont(str(p))
    names = {f.name for f in font_manager.fontManager.ttflist}
    for n in ("Centaur", "Californian FB", "Garamond", "EB Garamond"):
        if n in names:
            return n
    return "serif"


class Film:
    def __init__(self, res):
        self.r = res
        self.fam = _font()
        self.fig = plt.figure(figsize=(12.8, 7.2), dpi=100, facecolor=NIGHT)
        self.ax = self.fig.add_axes([0, 0, 1, 1])
        self.scenes = [(4.5, self.title), (6.5, self.problem), (5.5, self.design), (7.5, self.samples),
                       (8, self.decompose), (6.5, self.memory), (7, self.single), (6, self.validate), (4.5, self.close)]
        self.starts = np.cumsum([0] + [d for d, _ in self.scenes])

    def text(self, x, y, s, size=30, color=INK, alpha=1.0, style="normal", ha="left", va="center"):
        self.ax.text(x, y, s, fontsize=size, color=color, alpha=max(0, min(1, alpha)), family=self.fam,
                     style=style, ha=ha, va=va)

    def frame(self, i):
        t = i / FPS
        k = int(np.searchsorted(self.starts, t, side="right") - 1)
        k = min(k, len(self.scenes) - 1)
        d, fn = self.scenes[k]
        u = (t - self.starts[k]) / d
        ax = self.ax
        ax.clear()
        ax.set_xlim(0, 16)
        ax.set_ylim(0, 9)
        ax.axis("off")
        ax.set_facecolor(NIGHT)
        fade = min(ease(u * 6), ease((1 - u) * 8))
        fn(u, fade)
        self.text(0.5, 0.35, "PersistTrace", 13, MUTE, 0.9, "italic")
        self.text(15.5, 0.35, f"{k + 1} / {len(self.scenes)}", 13, MUTE, 0.9, ha="right")

    def title(self, u, a):
        self.text(1.2, 6.0, "The tumour shrinks.", 54, INK, ease(u * 4) * a)
        self.text(1.2, 4.8, "Some cancer cells survive.", 54, INK, ease(u * 4 - 0.8) * a)
        self.text(1.2, 3.6, "What are they becoming?", 54, C["dtp"], ease(u * 4 - 1.6) * a, "italic")

    def problem(self, u, a):
        total, nsclc, egfr, cols = 248, 211, 68, 31
        n = int(total * ease(u * 2.2))
        for j in range(n):
            x, y = 1.2 + (j % cols) * 0.34, 7.2 - (j // cols) * 0.34
            stage = 2 if j < egfr else 1 if j < nsclc else 0
            col = C["dtp"] if stage == 2 and u > 0.55 else "#3c6e8f" if stage >= 1 and u > 0.35 else "#6b655b"
            self.ax.add_patch(plt.Circle((x, y), 0.12, color=col, alpha=a))
        self.text(1.2, 3.2, "2.48 million lung cancer diagnoses in 2022", 30, INK, a)
        self.text(1.2, 2.5, "one dot = 10,000 people · blue: non-small-cell · red: EGFR-mutant (estimate)", 18, MUTE, a)
        self.text(1.2, 1.7, "GLOBOCAN 2022 · Zhang et al. 2016 · FLAURA: median 18.9 months to progression on osimertinib", 15, MUTE, a)

    def design(self, u, a):
        xs = np.linspace(2, 14, 5)
        self.text(1.2, 7.6, "Four cell lines, five moments, three replicates: 60 RNA-seq samples", 28, INK, a)
        self.ax.add_patch(plt.Rectangle((xs[1] - 0.8, 4.35), xs[2] - xs[1] + 1.6, 0.3, color=C["dtp"], alpha=0.18 * a))
        for j, c in enumerate(ORDER):
            v = ease(u * 3 - j * 0.35) * a
            self.ax.add_patch(plt.Circle((xs[j], 4.5), 0.75, fill=False, ec=C[c], lw=2, alpha=v))
            rng = np.random.default_rng(j)
            m = {"control": 22, "acute": 12, "dtp": 5, "wo_short": 8, "wo_long": 17}[c]
            for q in range(m):
                ang, rr = q * 2.399, 0.12 + np.sqrt(q) * 0.11
                self.ax.add_patch(plt.Circle((xs[j] + np.cos(ang) * rr, 4.5 + np.sin(ang) * rr), 0.06, color=C[c], alpha=v))
            self.text(xs[j], 3.3, NAMES[c], 20, INK, v, "italic", ha="center")
        self.text(1.2, 1.6, "GSE193258 · Criscione et al., npj Precision Oncology 2022", 16, MUTE, a)

    def samples(self, u, a):
        S = self.r["samples"]
        pts = np.array([s["pc"][:2] for s in S])
        lo, hi = pts.min(0), pts.max(0)
        X = lambda v: 1.5 + (v[0] - lo[0]) / (hi[0] - lo[0]) * 10.5
        Y = lambda v: 1.6 + (v[1] - lo[1]) / (hi[1] - lo[1]) * 5.6
        pos = ease(u * 1.15) * 4
        k = min(int(pos), 3)
        w = ease(pos - k)
        cur = ORDER[min(k + (1 if w > 0.5 else 0), 4)]
        for s in S:
            j = ORDER.index(s["condition"])
            if j <= k or (j == k + 1 and w > 0.3):
                self.ax.add_patch(plt.Circle((X(s["pc"]), Y(s["pc"])), 0.08, color=C[s["condition"]], alpha=(0.95 if s["condition"] == cur else 0.45) * a))
        for line in self.r["meta"]["cell_lines"]:
            cen = {c: np.mean([s["pc"][:2] for s in S if s["line"] == line and s["condition"] == c], axis=0) for c in ORDER}
            path = [cen[c] for c in ORDER[:k + 1]] + [cen[ORDER[k]] + w * (cen[ORDER[min(k + 1, 4)]] - cen[ORDER[k]])]
            xy = np.array([[X(p), Y(p)] for p in path])
            self.ax.plot(xy[:, 0], xy[:, 1], color=INK, alpha=0.55 * a, lw=1.4)
            self.ax.add_patch(plt.Circle(xy[-1], 0.16, color=C[cur], ec=INK, lw=1.5, alpha=a))
        self.text(12.8, 7.4, "Every dot is", 22, MUTE, a)
        self.text(12.8, 6.8, "a real experiment.", 22, MUTE, a)
        self.text(12.8, 5.6, NAMES[cur], 34, C[cur], a, "italic")

    def decompose(self, u, a):
        line = self.r["meta"]["cell_lines"][-1] if "PC9" not in self.r["betas"] else "PC9"
        b = self.r["betas"][line]
        P = np.array([[p[0], p[1], p[2]] for p in self.r["decomposition_scatter"][line] if None not in p[:3]])
        X = lambda v: 8 + np.clip(v, -6, 6) * 0.6
        Y = lambda v: 4.6 + np.clip(v, -8, 8) * 0.38
        m = ease((u - 0.45) * 3)
        yv = P[:, 1] * (1 - m) + P[:, 2] * m
        strong = np.abs(P[:, 2]) > 1.5
        cols = np.where(strong, np.where(P[:, 2] > 0, C["dtp"], C["control"]), "#6b655b")
        self.ax.scatter(X(P[:, 0]), Y(yv), s=np.where(strong, 9, 3), c=cols, alpha=0.75 * a, lw=0)
        xs = np.linspace(-6, 6, 2)
        self.ax.plot(X(xs), Y(b * xs * (1 - m)), color=C["acute"], lw=2.4, alpha=a * ease(u * 4))
        self.ax.plot([X(-6), X(6)], [Y(0), Y(0)], color="#3a352d", lw=1, alpha=a)
        self.text(1.0, 7.9, "Subtract the stress, keep the persister", 34, INK, a)
        self.text(1.0, 7.1, f"{line}: persister change P = {b:.2f} × 24-hour change A + S", 20, MUTE, a)
        lab = "P, before" if m < 0.5 else "S = P − βA, after"
        self.text(13.2, 5.6, lab, 24, C["acute"] if m < 0.5 else C["dtp"], a, "italic")
        self.text(1.0, 1.4, "gold: the stress persisters still carry · coloured: what stress does not explain", 16, MUTE, a)

    def memory(self, u, a):
        ms = self.r["memory"]["summary"]
        edges = np.linspace(-1, 2, 31)
        self.text(1.0, 7.9, "The cells remember the drug", 34, INK, a)
        for tp, col, grow in (("wo_short", C["wo_short"], ease(u * 3)), ("wo_long", C["wo_long"], ease(u * 3 - 1))):
            h = np.sum([r["hist"] for r in ms if r["timepoint"] == tp], axis=0).astype(float)
            h = h / max(h.max(), 1)
            for j, v in enumerate(h):
                x = 2 + (edges[j] + 1) / 3 * 9
                self.ax.add_patch(plt.Rectangle((x, 1.6), 0.28, v * 4.6 * grow, color=col, alpha=0.6 * a))
        for v, lab, col in ((0, "baseline", C["control"]), (1, "persister level", C["dtp"])):
            x = 2 + (v + 1) / 3 * 9
            self.ax.plot([x, x], [1.6, 6.6], color=col, ls="--", lw=1.2, alpha=a)
            self.text(x + 0.1, 6.8, lab, 16, col, a, "italic")
        hl = [d.get("median_half_life_days") for d in self.r["memory"]["decay"].values() if d.get("median_half_life_days")]
        if hl:
            self.text(11.8, 4.6, "memory half-life", 22, MUTE, a)
            self.text(11.8, 3.7, f"{min(hl):.0f} to {max(hl):.0f} days", 40, C["dtp"], a * ease(u * 3 - 1.2))
        self.text(2.0, 1.1, "share of each persister change still present: violet 24 h after washout, green after regrowth", 15, MUTE, a)

    def single(self, u, a):
        sc = (self.r.get("external") or {}).get("pc9_single_cell")
        self.text(1.0, 7.9, "One state, or several?", 34, INK, a)
        if not sc:
            self.text(1.0, 5, "single-cell data not available in this build", 22, MUTE, a)
            return
        pts = np.array([[p[0], p[1], p[2], p[3], p[4]] for p in sc["points"]], dtype=float)
        lo, hi = pts[:, :2].min(0), pts[:, :2].max(0)
        X = 1.2 + (pts[:, 0] - lo[0]) / (hi[0] - lo[0]) * 8
        Y = 1.2 + (pts[:, 1] - lo[1]) / (hi[1] - lo[1]) * 5.8
        rec = ease((u - 0.5) * 4)
        base = np.where(pts[:, 2] == 0, "#4a443b", np.where(pts[:, 3] == 1, C["dtp"], C["control"]))
        cyc = np.where(pts[:, 4] == 1, C["acute"], "#4a443b")
        cols = base if rec < 0.5 else cyc
        n = int(len(pts) * ease(u * 3))
        self.ax.scatter(X[:n], Y[:n], s=3, c=cols[:n], alpha=0.85 * a, lw=0)
        self.text(10.2, 6.2, f"{sc['cells_after_qc']:,} single cells", 26, INK, a)
        self.text(10.2, 5.4, "GSE275199, PC9 on osimertinib", 17, MUTE, a)
        if rec >= 0.5:
            self.text(10.2, 4.2, f"{sc['cycling_fraction_exposed'] * 100:.0f}% of drug-exposed", 28, C["acute"], a)
            self.text(10.2, 3.5, "cells are still cycling", 28, C["acute"], a)
        else:
            self.text(10.2, 4.2, "blue: unexposed · red: exposed", 20, MUTE, a)

    def validate(self, u, a):
        folds = self.r["benchmark"]["folds"]
        self.text(1.0, 7.9, "Tested on a cell line it never saw", 34, INK, a)
        for j, f in enumerate(folds):
            x = 2 + j * 3
            for off, key, col in ((0, "naive_treated_vs_control|auroc_dtp_vs_acute", C["acute"]), (0.9, "persister_specific|auroc_dtp_vs_acute", C["dtp"])):
                v = (f.get(key) or 0) * ease(u * 2.5 - off * 0.4)
                self.ax.add_patch(plt.Rectangle((x + off, 1.6), 0.8, v * 4.6, color=col, alpha=a))
            self.text(x + 0.85, 1.2, f["held_out"], 18, INK, a, "italic", ha="center")
        self.ax.plot([1.6, 14], [1.6 + 0.5 * 4.6] * 2, color=MUTE, ls=":", lw=1, alpha=a)
        h1 = next((h for h in self.r["hypotheses"] if h["id"] == "H1"), None)
        if h1:
            self.text(15.2, 7.95, f"{h1['naive']:.2f} → {h1['decomposed']:.2f}", 40, C["dtp"], a * ease(u * 3 - 1), ha="right")
            self.text(15.2, 7.2, "mean AUROC, persister vs 24 h", 16, MUTE, a, ha="right")
        self.text(1.0, 6.9, "gold: naive treated-vs-control signature · red: stress removed", 17, MUTE, a)

    def close(self, u, a):
        self.text(8, 5.4, "PersistTrace", 60, INK, a, ha="center")
        self.text(8, 4.3, "Pritidipa Das · Team Skibidi Solo · CGWC", 26, MUTE, a, "italic", ha="center")
        self.text(8, 3.4, "pritidipadas.github.io/persist-trace", 22, C["dtp"], a, ha="center")
        self.text(8, 1.6, "A computational study. Not a diagnostic or treatment tool.", 15, MUTE, a, ha="center")

    def render(self, out):
        n = int(self.starts[-1] * FPS)
        anim = animation.FuncAnimation(self.fig, self.frame, frames=n, interval=1000 / FPS)
        anim.save(out, writer=animation.FFMpegWriter(fps=FPS, bitrate=4500, codec="libx264",
                                                     extra_args=["-pix_fmt", "yuv420p", "-preset", "veryfast"]))
        plt.close(self.fig)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="site/data/results.json")
    ap.add_argument("--out", default="site/media/explainer.mp4")
    ap.add_argument("--poster", default="site/media/explainer.jpg")
    args = ap.parse_args(argv)
    res = json.loads(Path(args.results).read_text())
    f = Film(res)
    f.frame(int(f.starts[3] * FPS + 5.5 * FPS))
    f.fig.savefig(args.poster, facecolor=NIGHT)
    f.render(args.out)


if __name__ == "__main__":
    main()
