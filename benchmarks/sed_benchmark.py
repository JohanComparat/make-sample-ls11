#!/usr/bin/env python
"""Benchmark of the stellar-mass / K-correction codes at fixed redshift on LS DR11 BGS-like galaxies.

All codes see the same galaxies, photometry (dereddened grizW1W2 + error floor), filter curves and
redshift (BEST_Z). Outputs in $LS11_OUT/benchmark/ (or --out):

  <code>.fits                  per-galaxy results (LOGMSTAR, LO, HI, LOGSFR, MABS_R, CHI2) + timing
  <code>_zlo.fits, _zhi.fits   refits of a photo-z subset at BEST_Z -/+ BEST_Z_ERR
  report.md, *.png             metrics and figures

    python benchmarks/sed_benchmark.py --codes kcorrect,lephare,cigale,eazy,dsps --n 10000
"""

from __future__ import annotations

import argparse
import logging
import os
import time
from pathlib import Path

import numpy as np

from ls11samples import catalog, io
from ls11samples.config import load_config
from ls11samples.env import get_paths, nproc
from ls11samples.sed import OUTPUT, get_backend
from ls11samples.sed.common import BANDS

log = logging.getLogger("benchmark")
PHOT = ["RA", "DEC"] + [f"{p}_{b}" for p in ("FLUX", "FLUX_IVAR", "MW_TRANSMISSION") for b in BANDS]
N_FULL = 20e6       # ~770 deg^-2 x ~26 000 deg^2 (DR11 south footprint, before masking)


def nmad(x):
    x = np.asarray(x)[np.isfinite(x)]
    return float(1.4826 * np.median(np.abs(x - np.median(x)))) if x.size else np.nan


def select_sample(data: dict, n: int, rng) -> np.ndarray:
    """Half spectroscopic redshifts (or all of them if fewer), the rest random photo-z galaxies."""
    spec = np.flatnonzero(data["Z_SOURCE"] == 1)
    phot = np.flatnonzero(data["Z_SOURCE"] > 1)
    ns = min(spec.size, n // 2)
    idx = np.concatenate([rng.choice(spec, ns, replace=False), rng.choice(phot, min(phot.size, n - ns), replace=False)])
    return np.sort(idx)


def run_code(code, cfg, t, outdir, overwrite, tag=""):
    path = outdir / f"{code}{tag}.fits"
    if path.exists() and not overwrite:
        return io.read_table(path), dict(io.read_header(path, 1))
    t0 = time.time()
    backend = get_backend(code)(cfg)
    t1 = time.time()
    res = backend.fit(t)
    t2 = time.time()
    threads = {"kcorrect": 1, "lephare": int(os.environ.get("OMP_NUM_THREADS", os.cpu_count())),
               "dsps": os.cpu_count()}.get(code, nproc())      # threads the code used (upper bound)
    hdr = {"CODE": code, "NGAL": len(t["BEST_Z"]), "T_INIT": t1 - t0, "T_FIT": t2 - t1, "THREADS": threads}
    io.write_table(path, {"LS_ID_DR11": t["LS_ID_DR11"], **{k: np.asarray(v, np.float32) for k, v in res.items()}},
                   header=hdr, extname="RESULTS")
    log.info("%s%s: init %.1fs, fit %.1fs for %d galaxies", code, tag, t1 - t0, t2 - t1, len(t["BEST_Z"]))
    return res, hdr


def match_dr10(data: dict, path: Path, radius_arcsec: float = 1.0) -> dict | None:
    """DR10 LePhare masses (LPH_MASS_BEST) and redshifts matched by position, or None."""
    if not path.exists():
        return None
    from scipy.spatial import cKDTree

    pos = io.read_table(path, ["RA", "DEC"])
    ra0, ra1 = np.min(data["RA"]) - 0.1, np.max(data["RA"]) + 0.1
    d0, d1 = np.min(data["DEC"]) - 0.1, np.max(data["DEC"]) + 0.1
    rows = np.flatnonzero((pos["RA"] > ra0) & (pos["RA"] < ra1) & (pos["DEC"] > d0) & (pos["DEC"] < d1))
    if rows.size == 0:
        return None
    dr10 = io.read_table(path, ["RA", "DEC", "LPH_MASS_BEST", "LPH_MAG_ABS1", "Z_SPEC", "Z_PHOT_MEAN", "Z_PHOT_MEAN_I"], rows=rows)

    def xyz(ra, dec):
        ra, dec = np.radians(ra), np.radians(dec)
        return np.c_[np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)]

    dist, j = cKDTree(xyz(dr10["RA"], dr10["DEC"])).query(xyz(data["RA"], data["DEC"]))
    ok = dist < np.radians(radius_arcsec / 3600)
    z10 = np.where(dr10["Z_PHOT_MEAN_I"] > -9, dr10["Z_PHOT_MEAN_I"], dr10["Z_PHOT_MEAN"])
    z10 = np.where(dr10["Z_SPEC"] > -9, dr10["Z_SPEC"], z10)
    out = {"MATCH": ok, "LOGMSTAR": np.where(ok, dr10["LPH_MASS_BEST"][j], np.nan),
           "MABS_R": np.where(ok, dr10["LPH_MAG_ABS1"][j], np.nan), "BEST_Z": np.where(ok, z10[j], np.nan)}
    out["LOGMSTAR"] = np.where(out["LOGMSTAR"] > 0, out["LOGMSTAR"], np.nan)
    return out


# --------------------------------------------------------------------------- figures
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]     # categorical slots 1-5
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"


def _style(ax):
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)


def figures(res: dict, ref: str, z: np.ndarray, outdir: Path, timing: dict) -> list[str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    names = []
    others = [c for c in res if c != ref]
    # 1. logM(code) vs logM(ref), one panel per code (density, single hue)
    fig, axes = plt.subplots(1, len(others), figsize=(3.2 * len(others), 3.3), sharey=True, squeeze=False)
    for ax, c in zip(axes[0], others):
        x, y = res[ref]["LOGMSTAR"], res[c]["LOGMSTAR"]
        ok = np.isfinite(x) & np.isfinite(y)
        ax.hexbin(x[ok], y[ok] - x[ok], gridsize=45, extent=(8.5, 12, -1, 1), cmap="Blues", mincnt=1, bins="log", linewidths=0)
        ax.axhline(0, color=INK2, lw=1)
        ax.set_title(c, fontsize=10, color=INK)
        ax.set_xlabel(f"log M* ({ref})", fontsize=9, color=INK2)
        _style(ax)
    axes[0, 0].set_ylabel("log M* (code) − log M* (ref)", fontsize=9, color=INK2)
    fig.tight_layout()
    fig.savefig(outdir / "dlogm_vs_logm.png", dpi=130)
    plt.close(fig)
    names.append("dlogm_vs_logm.png")
    # 2. median offset vs z, one line per code (<= 4 series: legend + direct labels)
    zb = np.arange(0.0, 0.62, 0.05)
    zc = 0.5 * (zb[1:] + zb[:-1])
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    for i, c in enumerate(others[:4]):
        d = res[c]["LOGMSTAR"] - res[ref]["LOGMSTAR"]
        med = [np.nanmedian(d[(z >= a) & (z < b)]) if np.sum((z >= a) & (z < b) & np.isfinite(d)) > 30 else np.nan
               for a, b in zip(zb[:-1], zb[1:])]
        ax.plot(zc, med, color=SERIES[i], lw=2, marker="o", ms=4, label=c)
        last = np.flatnonzero(np.isfinite(med))
        if last.size:
            ax.annotate(c, (zc[last[-1]], med[last[-1]]), xytext=(4, 0), textcoords="offset points",
                        fontsize=8, color=INK2, va="center")
    ax.axhline(0, color=INK2, lw=1)
    ax.set_xlabel("BEST_Z", fontsize=9, color=INK2)
    ax.set_ylabel(f"median log M* − {ref}", fontsize=9, color=INK2)
    ax.legend(frameon=False, fontsize=8, loc="best")
    _style(ax)
    fig.tight_layout()
    fig.savefig(outdir / "dlogm_vs_z.png", dpi=130)
    plt.close(fig)
    names.append("dlogm_vs_z.png")
    # 3. extrapolated cost for the full DR11 footprint (single series)
    codes = list(timing)
    cpuh = [timing[c]["cpu_h_full"] for c in codes]
    fig, ax = plt.subplots(figsize=(6.4, 0.5 + 0.45 * len(codes)))
    ax.barh(codes, cpuh, color=SERIES[0], height=0.55)
    for y, v in enumerate(cpuh):
        ax.annotate(f"{v:,.1f} CPU-h" if v < 10 else f"{v:,.0f} CPU-h", (v, y), xytext=(4, 0),
                    textcoords="offset points", va="center",
                    fontsize=8, color=INK)
    ax.set_xscale("log")
    ax.set_xlabel(f"CPU-hours for {N_FULL / 1e6:.0f}M galaxies (extrapolated)", fontsize=9, color=INK2)
    _style(ax)
    fig.tight_layout()
    fig.savefig(outdir / "cost.png", dpi=130)
    plt.close(fig)
    names.append("cost.png")
    return names


# --------------------------------------------------------------------------- report
def report(res, ref, data, sel, timing, zsens, dr10, outdir, figs) -> Path:
    z = data["BEST_Z"][sel]
    spec = data["Z_SOURCE"][sel] == 1
    gr = (data["MAG_G"] - data["MAG_R"])[sel]
    lines = [f"# Stellar-mass benchmark, LS DR11 BGS-like (local strip)", "",
             f"{len(z)} galaxies ({spec.sum()} with spectroscopic redshift), fixed z = BEST_Z, "
             f"grizW1W2 with the same filter curves and error floor for every code. Reference: **{ref}**.", "",
             "## Cost", "", "| code | init [s] | fit [s] | threads | ms / galaxy / thread | CPU-h for 20M |",
             "|---|---|---|---|---|---|"]
    for c, tm in timing.items():
        lines.append(f"| {c} | {tm['init']:.1f} | {tm['fit']:.1f} | {tm['threads']} | "
                     f"{1e3 * tm['fit'] * tm['threads'] / tm['n']:.2f} | "
                     + (f"{tm['cpu_h_full']:,.1f} |" if tm['cpu_h_full'] < 10 else f"{tm['cpu_h_full']:,.0f} |"))
    lines += ["", "## Agreement with the reference", "",
              "Δ = log M*(code) − log M*(ref); NMAD = 1.4826 median|Δ − median Δ|.", "",
              "| code | failed [%] | median χ² | median Δ | NMAD | median Δ (spec-z) | median Δ, g−r>0.8 | median Δ, g−r<0.6 | median ΔMr |",
              "|---|---|---|---|---|---|---|---|---|"]
    for c, r in res.items():
        d = r["LOGMSTAR"] - res[ref]["LOGMSTAR"]
        dm = r["MABS_R"] - res[ref]["MABS_R"]
        lines.append(f"| {c} | {100 * np.mean(~np.isfinite(r['LOGMSTAR'])):.2f} | {np.nanmedian(r['CHI2']):.2f} | "
                     f"{np.nanmedian(d):+.3f} | {nmad(d):.3f} | {np.nanmedian(d[spec]):+.3f} | "
                     f"{np.nanmedian(d[gr > 0.8]):+.3f} | {np.nanmedian(d[gr < 0.6]):+.3f} | "
                     + (f"{np.nanmedian(dm):+.3f} |" if np.isfinite(dm).any() else "n/a |"))
    lines += ["", "Median Δ in redshift bins:", "", "| code | " + " | ".join(
        f"{a:.1f}–{b:.1f}" for a, b in zip(np.arange(0, 0.6, 0.1), np.arange(0.1, 0.7, 0.1))) + " |",
              "|---|" + "---|" * 6]
    for c, r in res.items():
        d = r["LOGMSTAR"] - res[ref]["LOGMSTAR"]
        cells = [f"{np.nanmedian(d[(z >= a) & (z < a + 0.1)]):+.3f}" if np.sum((z >= a) & (z < a + 0.1)) > 30 else "–"
                 for a in np.arange(0, 0.6, 0.1)]
        lines.append(f"| {c} | " + " | ".join(cells) + " |")
    if zsens:
        lines += ["", "## Photo-z sensitivity", "",
                  "Refit of photo-z galaxies at BEST_Z ± BEST_Z_ERR: half the difference of the two "
                  "log M*, i.e. the mass error from the photo-z alone.", "",
                  "| code | N | median ½\\|Δ log M*\\| | 84th pct |", "|---|---|---|---|"]
        for c, v in zsens.items():
            lines.append(f"| {c} | {v['n']} | {v['med']:.3f} | {v['p84']:.3f} |")
    if dr10 is not None:
        lines += ["", "## Comparison with the DR10 LePhare masses", "",
                  f"{dr10['n']} galaxies matched within 1″ with |Δz| < 0.005 between DR10 and DR11 BEST_Z.", "",
                  "| code | median log M*(DR11 code) − DR10 LPH_MASS_BEST | NMAD |", "|---|---|---|"]
        for c, v in dr10["rows"].items():
            lines.append(f"| {c} | {v[0]:+.3f} | {v[1]:.3f} |")
    lines += ["", "## Figures", ""] + [f"![{f}]({f})" for f in figs] + [""]
    path = outdir / "report.md"
    path.write_text("\n".join(lines))
    return path


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--codes", default="kcorrect,lephare,cigale,eazy,dsps")
    p.add_argument("--ref", default="lephare")
    p.add_argument("--n", type=int, default=10000)
    p.add_argument("--n-zsens", type=int, default=1000)
    p.add_argument("--seed", type=int, default=2)
    p.add_argument("--out", default=None)
    p.add_argument("--dr10", default=os.path.expanduser("~/data/legacysurvey/dr10/sweep/MergeALL_BGSlike_LPH.fits"))
    p.add_argument("--overwrite", action="store_true")
    args = p.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(levelname)s: %(message)s")
    cfg = load_config()
    paths = get_paths()
    outdir = Path(args.out) if args.out else paths.run_dir(cfg["tag"], "benchmark")
    outdir.mkdir(parents=True, exist_ok=True)
    data = catalog.load(paths, paths.sweeps(), cfg["tag"], sweep_columns=PHOT, mag_bands=["G", "R", "Z"])
    data = io.take(data, (np.asarray(data["STAR_FLAG"]).astype(int) & 2) == 0)
    rng = np.random.default_rng(args.seed)
    sel = select_sample(data, args.n, rng)
    t = io.take(data, sel)
    codes = args.codes.split(",")
    ref = args.ref if args.ref in codes else codes[0]

    res, timing = {}, {}
    for c in codes:
        r, h = run_code(c, cfg, t, outdir, args.overwrite)
        res[c] = r
        thr = int(h["THREADS"])
        timing[c] = {"init": float(h["T_INIT"]), "fit": float(h["T_FIT"]), "threads": thr, "n": int(h["NGAL"]),
                     "cpu_h_full": float(h["T_FIT"]) * thr / int(h["NGAL"]) * N_FULL / 3600}

    zsens = {}
    phot = np.flatnonzero((t["Z_SOURCE"] > 1) & np.isfinite(t["BEST_Z_ERR"]) & (t["BEST_Z_ERR"] > 0))
    if args.n_zsens and phot.size:
        sub = io.take(t, np.sort(rng.choice(phot, min(args.n_zsens, phot.size), replace=False)))
        for c in codes:
            lo = dict(sub, BEST_Z=np.clip(sub["BEST_Z"] - sub["BEST_Z_ERR"], 0.002, None))
            hi = dict(sub, BEST_Z=sub["BEST_Z"] + sub["BEST_Z_ERR"])
            rlo, _ = run_code(c, cfg, lo, outdir, args.overwrite, "_zlo")
            rhi, _ = run_code(c, cfg, hi, outdir, args.overwrite, "_zhi")
            dd = 0.5 * np.abs(rhi["LOGMSTAR"] - rlo["LOGMSTAR"])
            zsens[c] = {"n": int(np.isfinite(dd).sum()), "med": float(np.nanmedian(dd)), "p84": float(np.nanpercentile(dd, 84))}

    dr10 = None
    m = match_dr10(t, Path(args.dr10))
    if m is not None:
        same = m["MATCH"] & (np.abs(m["BEST_Z"] - t["BEST_Z"]) < 0.005)
        dr10 = {"n": int(same.sum()), "rows": {}}
        for c in codes:
            d = (res[c]["LOGMSTAR"] - m["LOGMSTAR"])[same]
            dr10["rows"][c] = (float(np.nanmedian(d)), nmad(d))

    figs = figures(res, ref, t["BEST_Z"], outdir, timing)
    path = report(res, ref, data, sel, timing, zsens, dr10, outdir, figs)
    log.info("report: %s", path)


if __name__ == "__main__":
    main()
