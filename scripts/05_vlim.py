#!/usr/bin/env python
"""Step 5: Mr and M* volume-limited samples with their randoms ($LS11_OUT/vlim/)."""

import logging

import numpy as np

from ls11samples import bgsl, cli, export, io, kcorr, vlim
from ls11samples.config import config_hash
from ls11samples.merge import join_by_id

log = logging.getLogger("ls11samples.step5")
DATA_COLS = ["RA", "DEC", "EBV", "BEST_Z", "BEST_Z_ERR", "Z_SOURCE", "LS_ID_DR11", "MAG_G", "MAG_R", "MAG_Z"]


def main():
    p = cli.parser(__doc__)
    p.add_argument("--code", default=None, help="stellar-mass code (default: sed.code in the config)")
    p.add_argument("--seed", type=int, default=1)
    args = p.parse_args()
    cfg, paths, _ = cli.setup(args)
    code = args.code or cfg["sed"]["code"]
    names = [s.name for s in paths.sweeps()]
    data = io.read_table(paths.data_file, DATA_COLS)
    dhdr = dict(io.read_header(paths.data_file, 1))
    area = float(dhdr["AREA"])
    r_lim = float(cfg["galaxy"]["r_range"][1])

    kc = join_by_id(data["LS_ID_DR11"], [paths.out / "kcorr" / kcorr.output_name(n) for n in names],
                    ["KCORR_R", "ABSMAG_R", "ABSMAG_R01", "KC_COEFFS"])
    data.update({k: kc[k] for k in ("KCORR_R", "ABSMAG_R", "ABSMAG_R01")})
    ms = join_by_id(data["LS_ID_DR11"], [paths.out / "mstar" / code / ("MS-" + n) for n in names], ["LOGMSTAR"])
    data["LOGMSTAR"] = ms["LOGMSTAR"]
    log.info("%d galaxies, %d with K-corrections, %d with %s masses", len(data["RA"]),
             kc["FOUND"].sum(), np.isfinite(ms["LOGMSTAR"]).sum(), code)

    mr_curve = None
    try:
        fitter = kcorr.KcorrectV5(cfg)
        sed_set = kcorr.complete_sed_set({**kc, "KC_COEFFS": kc["KC_COEFFS"]}, data["BEST_Z"], r_lim, cfg,
                                         z_max=cfg["kcorr"].get("complete_set_zmax", 0.1))
        zg = np.arange(0.01, 0.61, 0.01)
        kp = fitter.k_curve(sed_set, zg, cfg["vlim"].get("completeness_percentile", 95))
        mr_curve = vlim.mr_limit_from_kcurve(zg, kp, r_lim, cfg)
    except ImportError:
        log.warning("kcorrect not importable: Mr limit from the observed K-correction percentiles")

    samples, curves = vlim.define_samples(data, cfg, r_lim, mr_curve=mr_curve)
    outdir = paths.outdir("vlim")
    rand = io.read_table(paths.rand_file, ["RA", "DEC", "EBV"])
    rng = np.random.default_rng(args.seed)
    header = {"AREA": area, "RLIM": r_lim, "SEDCODE": code, "KCCODE": "kcorrect",
              "COSMO": f"FlatLCDM H0={cfg['cosmology']['H0']} Om0={cfg['cosmology']['Om0']}",
              "CFGHASH": config_hash(cfg), "FPHASH": dhdr.get("FPHASH", ""), "PARENT": paths.data_file.name}
    summary = []
    for s in samples:
        if s["n"] == 0:
            continue
        r = vlim.make_randoms(rand, np.asarray(data["BEST_Z"])[s["sel"]], cfg["vlim"]["n_rand_factor"], rng)
        export.write_sample(outdir, s, data, r, header)
        summary.append(vlim.summary_row(s, data, area, cfg))
        log.info("%s", s["name"])
    if summary:
        io.write_table(outdir / "vlim_boundaries.fits", {k: np.array([row[k] for row in summary]) for k in summary[0]},
                       header=header, extname="SAMPLES",
                       extra=[(f"LIMIT_{k.upper()}", {"Z": v[0], "LIMIT": v[1]}, None) for k, v in curves.items()])
    export.write_manifest(outdir / "manifest.yaml", cfg, summary,
                          {"area_deg2": area, "sed_code": code, "data_file": str(paths.data_file),
                           "rand_file": str(paths.rand_file)})


if __name__ == "__main__":
    main()
