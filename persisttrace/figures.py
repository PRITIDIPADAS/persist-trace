from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation

BG, FG, MUTED = "#07070f", "#ecebf5", "#8d8aa6"
COND = {"control": "#7dd3fc", "acute": "#fbbf24", "dtp": "#e879f9", "wo_short": "#a78bfa", "wo_long": "#34d399"}
ORDER = ["control", "acute", "dtp", "wo_short", "wo_long"]
LINE_MARK = {"PC9": "o", "H1975": "s", "HCC827": "^", "HCC2935": "D"}


def _style():
    plt.rcParams.update({
        "figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG,
        "axes.edgecolor": "#2a2840", "axes.labelcolor": FG, "xtick.color": MUTED, "ytick.color": MUTED,
        "text.color": FG, "font.family": "DejaVu Sans", "font.size": 10, "axes.grid": False,
        "svg.fonttype": "none",
    })


def _stamp(fig, synthetic):
    if synthetic:
        fig.text(0.5, 0.5, "SYNTHETIC TEST DATA", ha="center", va="center", fontsize=34, color="#ff4d6d",
                 alpha=0.35, rotation=20, weight="bold")


def pca_plot(res, path, synthetic):
    fig, ax = plt.subplots(figsize=(7, 5.2))
    for s in res["samples"]:
        ax.scatter(s["pc"][0], s["pc"][1], c=COND[s["condition"]], marker=LINE_MARK.get(s["line"], "o"),
                   s=46, edgecolors="none", alpha=0.9)
    v = res["meta"]["pca_variance"]
    ax.set_xlabel(f"PC1 ({v[0] * 100:.0f}% variance)")
    ax.set_ylabel(f"PC2 ({v[1] * 100:.0f}% variance)")
    for c in ORDER:
        ax.scatter([], [], c=COND[c], label=res["meta"]["conditions"][c], s=40)
    ax.legend(frameon=False, fontsize=8, loc="best")
    ax.set_title("60 RNA-seq samples, centred within each cell line", loc="left", fontsize=11)
    _stamp(fig, synthetic)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def decomposition_plot(res, path, synthetic):
    lines = list(res["decomposition_scatter"])
    fig, axes = plt.subplots(1, len(lines), figsize=(3.2 * len(lines), 3.4), sharex=True, sharey=True)
    for ax, l in zip(np.atleast_1d(axes), lines):
        pts = np.array([[p[0], p[1]] for p in res["decomposition_scatter"][l] if p[0] is not None and p[1] is not None])
        ax.scatter(pts[:, 0], pts[:, 1], s=3, c="#a78bfa", alpha=0.35, edgecolors="none", rasterized=True)
        b = res["betas"][l]
        xs = np.linspace(-6, 6, 10)
        ax.plot(xs, b * xs, c="#fbbf24", lw=1.4)
        ax.set_title(f"{l}  beta = {b:.2f}", fontsize=10)
        ax.set_xlim(-6, 6)
        ax.set_ylim(-8, 8)
        ax.set_xlabel("acute log2FC")
    np.atleast_1d(axes)[0].set_ylabel("persister log2FC")
    _stamp(fig, synthetic)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def memory_plot(res, path, synthetic):
    rows = res["memory"]["summary"]
    lines = sorted({r["cell_line"] for r in rows})
    fig, axes = plt.subplots(1, len(lines), figsize=(3.2 * len(lines), 3), sharey=True)
    edges = np.linspace(-1, 2, 31)
    for ax, l in zip(np.atleast_1d(axes), lines):
        for r in rows:
            if r["cell_line"] != l:
                continue
            col = COND["wo_short"] if r["timepoint"] == "wo_short" else COND["wo_long"]
            ax.stairs(r["hist"], edges, fill=True, alpha=0.45, color=col)
        ax.axvline(1, c=COND["dtp"], lw=1, ls="--")
        ax.axvline(0, c=COND["control"], lw=1, ls="--")
        ax.set_title(l, fontsize=10)
        ax.set_xlabel("fraction of persister change remaining")
    _stamp(fig, synthetic)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def benchmark_plot(res, path, synthetic):
    rows = res["benchmark"]["summary"]
    feats = ["top_variable", "naive_signature", "persister_signature"]
    models = sorted({r["model"] for r in rows})
    fig, ax = plt.subplots(figsize=(7.5, 3.6))
    w = 0.26
    cols = ["#7dd3fc", "#fbbf24", "#e879f9"]
    for i, f in enumerate(feats):
        vals = [next((r["balanced_accuracy"] for r in rows if r["model"] == m and r["features"] == f), np.nan) for m in models]
        lo = [next((r["ci_low"] for r in rows if r["model"] == m and r["features"] == f), np.nan) for m in models]
        hi = [next((r["ci_high"] for r in rows if r["model"] == m and r["features"] == f), np.nan) for m in models]
        x = np.arange(len(models)) + (i - 1) * w
        ax.bar(x, vals, w, color=cols[i], label=f.replace("_", " "))
        ax.errorbar(x, vals, yerr=[np.array(vals) - lo, np.array(hi) - vals], fmt="none", ecolor=FG, lw=0.8)
    ax.axhline(0.2, c=MUTED, ls=":", lw=1)
    ax.set_xticks(np.arange(len(models)), [m.replace("_", " ") for m in models], rotation=15)
    ax.set_ylabel("balanced accuracy (held-out line)")
    ax.set_ylim(0, 1.05)
    ax.legend(frameon=False, fontsize=8)
    _stamp(fig, synthetic)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def _centroids(res):
    out = {}
    for s in res["samples"]:
        out.setdefault(s["line"], {}).setdefault(s["condition"], []).append(s["pc"][:2])
    return {l: {c: np.mean(v, axis=0) for c, v in d.items()} for l, d in out.items()}


def video(res, mp4, gif, synthetic, fps=30, seconds_per_step=2.0):
    cents = _centroids(res)
    pts = np.array([s["pc"][:2] for s in res["samples"]])
    pad = 0.15 * (pts.max(axis=0) - pts.min(axis=0))
    lo, hi = pts.min(axis=0) - pad, pts.max(axis=0) + pad
    fig, ax = plt.subplots(figsize=(8, 4.5), dpi=120)
    fig.subplots_adjust(0.02, 0.02, 0.98, 0.9)
    steps = len(ORDER) - 1
    frames = int(fps * seconds_per_step * steps + fps * 1.5)

    def draw(f):
        ax.clear()
        ax.set_xlim(lo[0], hi[0])
        ax.set_ylim(lo[1], hi[1])
        ax.set_xticks([])
        ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
        t = min(f / (fps * seconds_per_step), steps)
        k = int(np.floor(t))
        u = t - k
        u = u * u * (3 - 2 * u)
        shown = ORDER[: min(k + 2, len(ORDER))]
        for s in res["samples"]:
            if s["condition"] in shown[: k + 1] or (s["condition"] == shown[-1] and u > 0.6):
                ax.scatter(*s["pc"][:2], s=18, c=COND[s["condition"]], alpha=0.35, edgecolors="none",
                           marker=LINE_MARK.get(s["line"], "o"))
        for l, c in cents.items():
            path = [c[o] for o in ORDER[: k + 1]]
            if k < steps:
                path.append(c[ORDER[k]] + u * (c[ORDER[k + 1]] - c[ORDER[k]]))
            path = np.array(path)
            ax.plot(path[:, 0], path[:, 1], c=FG, lw=1.2, alpha=0.7)
            ax.scatter(*path[-1], s=90, c=COND[ORDER[min(k + (1 if u > 0.5 else 0), steps)]], edgecolors=FG, linewidths=1.2, zorder=5)
            ax.text(path[-1][0], path[-1][1], f"  {l}", fontsize=8, color=FG, va="center")
        label = res["meta"]["conditions"][ORDER[min(k + (1 if u > 0.5 else 0), steps)]]
        ax.set_title(f"Real GSE193258 samples  |  {label}", loc="left", fontsize=11, color=FG)
        _stamp(fig, synthetic)

    anim = animation.FuncAnimation(fig, draw, frames=frames, interval=1000 / fps)
    anim.save(mp4, writer=animation.FFMpegWriter(fps=fps, bitrate=1800, codec="libx264",
                                                 extra_args=["-pix_fmt", "yuv420p", "-preset", "veryfast"]))
    plt.close(fig)
    import subprocess
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-i", str(mp4), "-vf",
                    "fps=15,scale=640:-1:flags=lanczos,split[a][b];[a]palettegen=max_colors=64[p];[b][p]paletteuse",
                    str(gif)], check=False)


def render_all(res, gene_rows, fig_dir, media_dir, synthetic=False):
    _style()
    fig_dir, media_dir = Path(fig_dir), Path(media_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)
    media_dir.mkdir(parents=True, exist_ok=True)
    tag = "_synthetic" if synthetic else ""
    pca_plot(res, fig_dir / f"fig1_embedding{tag}.svg", synthetic)
    decomposition_plot(res, fig_dir / f"fig2_decomposition{tag}.svg", synthetic)
    memory_plot(res, fig_dir / f"fig3_memory{tag}.svg", synthetic)
    benchmark_plot(res, fig_dir / f"fig4_benchmark{tag}.svg", synthetic)
    for f in fig_dir.glob(f"*{tag}.svg"):
        (media_dir / f.name).write_bytes(f.read_bytes())
    video(res, media_dir / f"trajectory{tag}.mp4", media_dir / f"trajectory{tag}.gif", synthetic)
