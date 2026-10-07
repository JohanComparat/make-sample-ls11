"""Validation of a run: completeness and row alignment of the per-sweep products, missing values,
distributions of the stellar masses and absolute magnitudes, and the agreement between codes.

Galaxies are the selected objects with STAR_FLAG & ``vlim.exclude_star_flag`` == 0 and BEST_Z
inside the redshift range of the SED grids (``kcorr.z_range``).
"""

from __future__ import annotations

import logging
from collections.abc import Sequence
from pathlib import Path

import numpy as np

from . import io
from .catalog import CODE_COLUMNS

log = logging.getLogger(__name__)


def nmad(x) -> float:
    x = np.asarray(x, np.float64)
    x = x[np.isfinite(x)]
    return float(1.4826 * np.median(np.abs(x - np.median(x)))) if x.size else float("nan")


def _pct(x, q=(5, 50, 95)) -> list[float]:
    x = np.asarray(x, np.float64)
    x = x[np.isfinite(x)]
    return [round(float(v), 4) for v in np.percentile(x, q)] if x.size else [float("nan")] * len(q)


def validate(paths, cfg: dict, sweeps: Sequence[Path], codes: Sequence[str] | None = None,
             frac: float = 0.02, ref: str | None = None, seed: int = 0) -> dict:
    """Summary dict of the run made with ``cfg`` over ``sweeps`` (see the module docstring)."""
    tag = cfg["tag"]
    codes = list(codes or cfg["sed"]["codes"])
    ref = ref or cfg["sed"].get("primary", codes[0])
    star_bits = int(cfg["vlim"].get("exclude_star_flag", 3))
    zmax = float(cfg["kcorr"].get("z_range", [0.0, 1.0])[1])
    rng = np.random.default_rng(seed)
    missing, misaligned = [], []
    n_all = n_gal = n_empty = 0
    nan = {c: [0, 0] for c in codes}
    keep = {c: {k: [] for k in CODE_COLUMNS} for c in codes}
    for s in sweeps:
        fsel = paths.product(s, tag)
        if not fsel.exists():
            missing.append([tag, s.name])
            continue
        nsel = int(io.read_header(fsel, 1)["NSEL"])
        n_all += nsel
        if nsel == 0:
            n_empty += 1
            continue
        sel = io.read_table(fsel, ["SWEEP_ROW", "STAR_FLAG", "BEST_Z"])
        gal = ((np.asarray(sel["STAR_FLAG"]).astype(int) & star_bits) == 0) & (sel["BEST_Z"] < zmax)
        n_gal += int(gal.sum())
        pick = gal & (rng.random(nsel) < frac)
        got = {}
        for c in codes:
            f = paths.product(s, c)
            if not f.exists():
                missing.append([c, s.name])
                continue
            d = io.read_table(f)
            if len(d["SWEEP_ROW"]) != nsel or not np.array_equal(d["SWEEP_ROW"], sel["SWEEP_ROW"]):
                misaligned.append([c, s.name])
                continue
            bad = ~np.isfinite(d["LOGMSTAR"])
            nan[c][0] += int(bad.sum())
            nan[c][1] += int((bad & gal).sum())
            got[c] = d
        if len(got) == len(codes):           # subsample only where every code exists (aligned)
            for c in codes:
                for k in CODE_COLUMNS:
                    keep[c][k].append(got[c][k][pick])
    out = {"tag": tag, "n_sweeps": len(sweeps), "n_empty_sweeps": n_empty, "n_objects": n_all,
           "n_galaxies": n_gal, "missing": missing, "misaligned": misaligned,
           "nan_fraction": {c: {"all": nan[c][0] / max(n_all, 1), "galaxies": nan[c][1] / max(n_gal, 1)}
                            for c in codes},
           "percentiles_5_50_95": {}, "vs_ref": {}, "sizes_gb": {}}
    v = {c: {k: np.concatenate(keep[c][k]) if keep[c][k] else np.zeros(0) for k in CODE_COLUMNS} for c in codes}
    for c in codes:
        out["percentiles_5_50_95"][c] = {k: _pct(v[c][k]) for k in CODE_COLUMNS}
    if ref in codes:
        for c in codes:
            if c == ref:
                continue
            dm = v[c]["LOGMSTAR"] - v[ref]["LOGMSTAR"]
            dr = v[c]["MABS_R"] - v[ref]["MABS_R"]
            out["vs_ref"][f"{c}-{ref}"] = {"dlogm_median": float(np.nanmedian(dm)) if dm.size else None,
                                          "dlogm_nmad": nmad(dm), "dmr_median": float(np.nanmedian(dr)) if dr.size else None,
                                          "dmr_nmad": nmad(dr)}
    for name in (tag, *codes):
        d = paths.product_dir(name)
        out["sizes_gb"][name] = round(sum(f.stat().st_size for f in d.glob("*.fits")) / 1e9, 3) if d.exists() else 0.0
    rf = paths.rand_file(tag)
    if rf.exists():
        h = io.read_header(rf, 1)
        out["randoms"] = {"n": int(h["NRAND"]), "area_deg2": float(h["AREA"]), "n_sweeps": int(h.get("NSWEEPS", -1)),
                          "size_gb": round(rf.stat().st_size / 1e9, 3)}
    return out


def _fmt(x, spec: str) -> str:
    """``x`` formatted with ``spec``; 'n/a' when there is no value (a code without outputs)."""
    return "n/a" if x is None or not np.isfinite(x) else format(x, spec)


def report(out: dict) -> str:
    lines = [f"run {out['tag']}: {out['n_sweeps']} sweeps ({out['n_empty_sweeps']} empty), "
             f"{out['n_objects']} objects, {out['n_galaxies']} galaxies",
             f"missing files: {len(out['missing'])} {out['missing'][:5]}",
             f"misaligned files: {len(out['misaligned'])} {out['misaligned'][:5]}"]
    for c, f in out["nan_fraction"].items():
        lines.append(f"{c:9s} NaN fraction: all {f['all']:.4f}, galaxies {f['galaxies']:.4f}")
    lines.append("galaxies, subsample, percentiles 5/50/95:")
    for c, p in out["percentiles_5_50_95"].items():
        lines.append(f"  {c:9s} " + "  ".join(f"{k} {v}" for k, v in p.items()))
    for k, v in out["vs_ref"].items():
        lines.append(f"  {k}: dlogM median {_fmt(v['dlogm_median'], '+.3f')} NMAD {_fmt(v['dlogm_nmad'], '.3f')}; "
                     f"dMr median {_fmt(v['dmr_median'], '+.3f')} NMAD {_fmt(v['dmr_nmad'], '.3f')}")
    if "randoms" in out:
        r = out["randoms"]
        lines.append(f"randoms: {r['n']} over {r['area_deg2']:.2f} deg2 ({r['n_sweeps']} sweeps), {r['size_gb']} GB")
    lines.append("sizes (GB): " + ", ".join(f"{k} {v}" for k, v in out["sizes_gb"].items()))
    return "\n".join(lines)
