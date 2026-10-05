import numpy as np
import pytest
from astropy.coordinates import SkyCoord

from ls11samples import bgsl, randoms, redshift, selection


def test_galactic_b_matches_astropy(rng):
    ra, dec = rng.uniform(0, 360, 1000), np.degrees(np.arcsin(rng.uniform(-1, 1, 1000)))
    b = SkyCoord(ra=ra, dec=dec, unit="deg").galactic.b.deg
    assert np.allclose(selection.galactic_b(ra, dec), b, atol=1e-5)


def test_is_south():
    ra = np.array([180.0, 0.0, 180.0, 10.0])
    dec = np.array([40.0, 40.0, 20.0, -10.0])
    assert selection.is_south(ra, dec).tolist() == [False, True, True, True]


def _footprint_table(n, rng):
    return {"RA": rng.uniform(0, 5, n), "DEC": rng.uniform(-5, 0, n), "EBV": rng.uniform(0, 0.15, n),
            "MASKBITS": rng.choice([0, 0, 0, 2, 2048, 4096, 8192, 1 << 16], n),
            "NOBS_G": rng.integers(0, 3, n), "NOBS_R": rng.integers(0, 3, n), "NOBS_Z": rng.integers(0, 3, n),
            "GALDEPTH_G": rng.uniform(-1, 10, n), "GALDEPTH_R": rng.uniform(0, 10, n),
            "GALDEPTH_Z": rng.uniform(0, 10, n)}


def test_footprint_same_for_galaxies_and_randoms(cfg, rng):
    t = _footprint_table(5000, rng)
    t["FITBITS"] = np.zeros(5000, np.int16)
    cuts = selection.footprint_cuts(t, cfg, fitbits=t["FITBITS"])
    flags = selection.sel_flags(cuts)
    assert np.array_equal(flags == 0, selection.footprint_mask(t, cfg))
    # SUB_BLOB (16) is not rejected, BRIGHT / MEDIUM / GALAXY / CLUSTER are
    mb = t["MASKBITS"]
    assert not cuts["maskbits"][np.isin(mb, [2, 2048, 4096, 8192])].any()
    assert cuts["maskbits"][np.isin(mb, [0, 1 << 16])].all()


def test_keep_sga_exempts_large_galaxies(cfg):
    t = {"RA": np.zeros(2), "DEC": np.zeros(2), "EBV": np.zeros(2), "MASKBITS": np.array([4096, 4096]),
         "NOBS_G": np.ones(2), "NOBS_R": np.ones(2), "NOBS_Z": np.ones(2),
         "GALDEPTH_G": np.ones(2), "GALDEPTH_R": np.ones(2), "GALDEPTH_Z": np.ones(2)}
    fit = np.array([0, 1 << 9])
    assert selection.footprint_cuts(t, cfg, fitbits=fit)["maskbits"].tolist() == [False, False]
    cfg["galaxy"]["keep_sga"] = True
    assert selection.footprint_cuts(t, cfg, fitbits=fit)["maskbits"].tolist() == [False, True]


def test_best_z_priority(cfg):
    pz = {"Z_SPEC": np.array([0.1, -99, -99, -99, 0.3]),
          "SURVEY": np.array(["DESI", "", "", "", "COSMOS2015"]),
          "Z_PHOT_MEAN_I": np.array([0.2, 0.25, -99, -99, 0.31]), "Z_PHOT_STD_I": np.full(5, 0.02),
          "Z_PHOT_MEAN": np.array([0.2, 0.26, 0.4, -99, 0.32]), "Z_PHOT_STD": np.full(5, 0.03)}
    z = redshift.best_z(pz, cfg)
    assert z["Z_SOURCE"].tolist() == [1, 2, 3, 0, 2]
    assert np.allclose(z["BEST_Z"][[0, 1, 2, 4]], [0.1, 0.25, 0.4, 0.31])
    assert np.isnan(z["BEST_Z"][3]) and z["BEST_Z_ERR"][0] == 0


def test_boxset_grid_equals_loop(rng):
    boxes = [(0, 5, -5, 0), (5, 10, -5, 0), (355, 360, 30, 35)]
    ra, dec = rng.uniform(-1, 361, 20000) % 360, rng.uniform(-10, 40, 20000)
    bs = randoms.BoxSet(boxes)
    assert bs.grid
    assert np.array_equal(bs.contains(ra, dec), randoms.io.in_boxes(ra, dec, boxes))


def test_cutflow_monotonic(cfg, first_sweep):
    t, flow, hdr = bgsl.select_sweep(first_sweep, cfg)
    assert np.all(np.diff(flow["N_PASS_CUMUL"]) <= 0)
    assert flow["N_PASS_CUMUL"][-1] == hdr["NSEL"] == len(t["SWEEP_ROW"])
    assert set(t) == {"LS_ID_DR11", "SWEEP_ROW", "BEST_Z", "BEST_Z_ERR", "Z_SOURCE", "STAR_FLAG"}
    assert np.all(np.diff(t["SWEEP_ROW"]) > 0) and np.all(t["LS_ID_DR11"] > 0)
    assert np.isfinite(t["BEST_Z"]).mean() > 0.99
    assert hdr["FPHASH"] and "BRIGHT(1)" in hdr["MASKREJ"]


def test_fracarea_uniform(rng):
    from ls11samples import maps
    n = 400_000
    ra = rng.uniform(0, 360, n)
    dec = np.degrees(np.arcsin(rng.uniform(-1, 1, n)))
    dens = n / 41252.96
    m = maps.fracarea_map(maps.pixel_index(ra, dec, 8), 8, dens)
    assert m.mean() == pytest.approx(1, abs=0.01)


def _galaxy_table(n, rng):
    return {"RA": np.zeros(n), "TYPE": rng.choice(["PSF", "REX", "DUP", "SER"], n),
            "FITBITS": rng.choice([0, 0, 2, 4096, 512], n).astype(np.int16),
            "FLUX_G": rng.uniform(-1, 30, n), "FLUX_R": rng.uniform(0.5, 60, n), "FLUX_Z": rng.uniform(-1, 90, n),
            "MW_TRANSMISSION_G": np.ones(n), "MW_TRANSMISSION_R": np.ones(n), "MW_TRANSMISSION_Z": np.ones(n),
            "FIBERFLUX_R": rng.uniform(0, 5, n), "FIBERTOTFLUX_R": rng.uniform(0, 5, n),
            "GAIA_PHOT_G_MEAN_MAG": rng.choice([0.0, 18.0, 20.5], n),
            "FLUX_IVAR_G": rng.uniform(-1, 1, n), "FLUX_IVAR_R": np.ones(n), "FLUX_IVAR_Z": np.ones(n),
            **{f"{q}_{b}": rng.uniform(0, 1, n) for q in ("FRACMASKED", "FRACIN", "FRACFLUX") for b in "GRZ"}}


def test_r21_config_only_bits_and_magnitude(rng):
    from ls11samples.config import load_config, DEFAULT

    cfg = load_config(DEFAULT.parent / "bgs_r21_dr10bits.yaml")
    assert cfg["tag"] == "bgsr21" and cfg["galaxy"]["r_range"] == [13.0, 21.0]
    assert cfg["footprint"]["maskbits_reject"] == ["NPRIMARY", "BRIGHT", "MEDIUM", "GALAXY", "CLUSTER"]
    assert cfg["cosmology"]["H0"] == 67.74                     # inherited from default.yaml
    t = _galaxy_table(4000, rng)
    cuts = selection.galaxy_cuts(t, cfg)
    for name in ("type", "flux_ivar", "gaia", "colour", "quality", "fiber", "rfibtot"):
        assert cuts[name].all(), name
    fit_rejected = np.isin(t["FITBITS"], [2, 4096])            # FIT_BACKGROUND, GAIA_POINTSOURCE
    assert np.array_equal(cuts["fitbits"], ~fit_rejected)      # LARGEGALAXY (512) is kept
    r = 22.5 - 2.5 * np.log10(t["FLUX_R"])
    assert np.array_equal(cuts["rmag"], (r > 13) & (r <= 21))


def test_star_flag():
    t = {"TYPE": np.array(["PSF", "PSF", "REX", "SER"]), "GAIA_PHOT_G_MEAN_MAG": np.array([0.0, 18.0, 18.0, 25.0]),
         "FLUX_R": np.full(4, 10 ** (-0.4 * (18.0 - 22.5)))}           # r_raw = 18
    assert selection.star_flag(t).tolist() == [1, 3, 2, 0]
