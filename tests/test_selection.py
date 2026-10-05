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
    assert flow["N_PASS_CUMUL"][-1] == hdr["NSEL"] == np.count_nonzero(t["SEL_FLAGS"] == 0)
    s = t["SEL_FLAGS"] == 0
    assert np.all(t["MAG_R"][s] <= cfg["galaxy"]["r_range"][1])
    assert np.all(t["LS_ID_DR11"] > 0) and np.isfinite(t["BEST_Z"][s]).mean() > 0.99
    assert hdr["FPHASH"]


def test_fracarea_uniform(rng):
    from ls11samples import maps
    n = 400_000
    ra = rng.uniform(0, 360, n)
    dec = np.degrees(np.arcsin(rng.uniform(-1, 1, n)))
    dens = n / 41252.96
    m = maps.fracarea_map(maps.pixel_index(ra, dec, 8), 8, dens)
    assert m.mean() == pytest.approx(1, abs=0.01)
