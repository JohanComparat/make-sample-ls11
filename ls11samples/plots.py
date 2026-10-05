"""Diagnostic figures (matplotlib, PNG)."""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

import numpy as np

INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"


def _style(ax):
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)


def vlim_planes(data: Mapping, samples: list[dict], curves: Mapping, path: Path, primary: str = "") -> Path:
    """Mr-z and M*-z planes (galaxy density), completeness limits and volume-limited sample boxes."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    z = np.asarray(data["BEST_Z"], float)
    panels = [("Mr", "ABSMAG_R", "M$_r$ (rest-frame DECam r)", (-24.5, -16), True),
              ("Mstar", "LOGMSTAR", f"log M* / M$_\\odot$ ({primary})", (8.5, 12.2), False)]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, (kind, col, label, ylim, invert) in zip(axes, panels):
        if col not in data:
            ax.set_visible(False)
            continue
        y = np.asarray(data[col], float)
        ok = np.isfinite(y) & np.isfinite(z)
        ax.hexbin(z[ok], y[ok], gridsize=70, extent=(0, 0.6, *ylim), cmap="Blues", bins="log", mincnt=1, linewidths=0)
        if kind in curves:
            zc, ml = curves[kind]
            ax.plot(zc, ml, color=INK, lw=1.6, label="completeness limit")
        for s in samples:
            if s["kind"] != kind:
                continue
            lo, hi = (s["hi"], max(s["lo"], ylim[0])) if kind == "Mr" else (s["lo"], min(s["hi"], ylim[1]))
            ax.add_patch(Rectangle((s["z0"], min(lo, hi)), s["z1"] - s["z0"], abs(hi - lo), fill=False,
                                   ec=INK2, lw=0.8))
        ax.set_xlim(0, 0.6)
        ax.set_ylim(ylim[::-1] if invert else ylim)
        ax.set_xlabel("BEST_Z", fontsize=9, color=INK2)
        ax.set_ylabel(label, fontsize=9, color=INK2)
        ax.legend(frameon=False, fontsize=8, loc="lower right" if not invert else "upper right")
        _style(ax)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path
