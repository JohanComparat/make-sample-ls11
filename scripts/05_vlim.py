#!/usr/bin/env python
"""Step 5: Mr and M* volume-limited samples with their randoms ($LS11_OUT/<tag>/vlim/).

Reads the per-sweep selection and code files (and RA, DEC, EBV, g/r/z fluxes from the sweeps at
SWEEP_ROW). M* from ``sed.primary`` (or --code), Mr from ``vlim.mr_code``; objects with
STAR_FLAG & ``vlim.exclude_star_flag`` are left out.
"""

import logging

import numpy as np

from ls11samples import catalog, cli, export, io, kcorr, vlim
from ls11samples.config import config_hash
from ls11samples.cosmo import distmod
from ls11samples.sed.common import BANDS

log = logging.getLogger("ls11samples.step5")
PHOT = [f"{p}_{b}" for p in ("FLUX", "FLUX_IVAR", "MW_TRANSMISSION") for b in BANDS]


def mr_curve_from_complete_set(data, sweeps, cfg, r_lim):
    """Mr completeness limit from the r-band K_p(z) of a complete low-z SED set (kcorrect refit)."""
    fitter = kcorr.KcorrectV5(cfg, mc=0)
    idx = kcorr.complete_set_rows(data["BEST_Z"], data["MABS_R"], data["MAG_R"], r_lim, cfg,
                                  z_max=cfg["kcorr"].get("complete_set_zmax", 0.1))
    parts = []
    for i in np.unique(data["SWEEP_INDEX"][idx]):
        rows = idx[data["SWEEP_INDEX"][idx] == i]
        t = catalog.sweep_rows(sweeps[i], data["SWEEP_ROW"][rows], PHOT)
        t["BEST_Z"] = data["BEST_Z"][rows]
        parts.append(t)
    coeffs = fitter.fit(io.concat(parts))["KC_COEFFS"]
    coeffs = coeffs[coeffs.sum(axis=1) > 0]
    zg = np.arange(0.01, 0.61, 0.01)
    kp = fitter.k_curve(coeffs, zg, cfg["vlim"].get("completeness_percentile", 95))
    return vlim.mr_limit_from_kcurve(zg, kp, r_lim, cfg)


def main():
    p = cli.parser(__doc__)
    p.add_argument("--code", default=None,
                   help="stellar-mass code defining the M* samples (default: sed.primary; another code "
                        "writes to vlim_<code>/)")
    p.add_argument("--seed", type=int, default=1)
    args = p.parse_args()
    cfg, paths, _ = cli.setup(args)
    tag = cfg["tag"]
    vc = cfg["vlim"]
    primary = args.code or cfg["sed"]["primary"]
    mr_code = vc.get("mr_code", "kcorrect")
    sweeps = [s for s in paths.sweeps() if paths.product(s, tag).exists()]
    codes = [c for c in dict.fromkeys([primary, mr_code, *cfg["sed"]["codes"]])
             if all(paths.product(s, c).exists() for s in sweeps)]
    for c in (primary, mr_code):
        if c not in codes:
            raise SystemExit(f"{c} outputs missing for some sweeps: run step 4 first")
    data = catalog.load(paths, sweeps, tag, sweep_columns=["RA", "DEC", "EBV"], codes=codes,
                        mag_bands=["G", "R", "Z"])
    star = (np.asarray(data["STAR_FLAG"]).astype(int) & int(vc.get("exclude_star_flag", 2))) != 0
    if star.any():
        data = io.take(data, ~star)
    data["LOGMSTAR"] = data[f"LOGMSTAR_{primary.upper()}"]
    data["MABS_R"] = data[f"MABS_R_{mr_code.upper()}"]
    data["KCORR_R"] = (data["MAG_R"] - distmod(data["BEST_Z"], cfg) - data["MABS_R"]).astype(np.float32)
    rhdr = dict(io.read_header(paths.rand_file(tag), 1))
    area = float(rhdr["AREA"])
    r_lim = float(cfg["galaxy"]["r_range"][1])
    log.info("%d objects (%d stars left out), codes %s; M* from %s, Mr from %s", len(data["RA"]),
             int(star.sum()), codes, primary, mr_code)

    try:
        mr_curve = mr_curve_from_complete_set(data, sweeps, cfg, r_lim)
    except ImportError:
        log.warning("kcorrect not importable: Mr limit from the observed K-correction percentiles")
        mr_curve = None
    samples, curves = vlim.define_samples(data, cfg, r_lim, mr_curve=mr_curve)
    outdir = paths.run_dir(tag, "vlim" if primary == cfg["sed"]["primary"] else f"vlim_{primary}")
    rand = io.read_table(paths.rand_file(tag), ["RA", "DEC", "EBV"])
    rng = np.random.default_rng(args.seed)
    header = {"AREA": area, "RLIM": r_lim, "TAG": tag, "SEDCODE": primary, "SEDCODES": ",".join(codes),
              "MRCODE": mr_code, "COSMO": f"FlatLCDM H0={cfg['cosmology']['H0']} Om0={cfg['cosmology']['Om0']}",
              "CFGHASH": config_hash(cfg), "FPHASH": rhdr.get("FPHASH", ""), "PARENT": tag}
    summary = []
    for s in samples:
        if s["n"] == 0:
            continue
        r = vlim.make_randoms(rand, np.asarray(data["BEST_Z"])[s["sel"]], vc["n_rand_factor"], rng)
        export.write_sample(outdir, s, data, r, header)
        summary.append(vlim.summary_row(s, data, area, cfg))
        log.info("%s", s["name"])
    if summary:
        io.write_table(outdir / "vlim_boundaries.fits", {k: np.array([row[k] for row in summary]) for k in summary[0]},
                       header=header, extname="SAMPLES",
                       extra=[(f"LIMIT_{k.upper()}", {"Z": v[0], "LIMIT": v[1]}, None) for k, v in curves.items()])
    try:
        from ls11samples.plots import vlim_planes

        vlim_planes(data, samples, curves, outdir / "vlim_planes.png", primary)
    except ImportError:
        log.warning("matplotlib not available: no vlim_planes.png")
    export.write_manifest(outdir / "manifest.yaml", cfg, summary,
                          {"area_deg2": area, "tag": tag, "sed_code": primary, "sed_codes": codes,
                           "mr_code": mr_code, "rand_file": str(paths.rand_file(tag))})


if __name__ == "__main__":
    main()
