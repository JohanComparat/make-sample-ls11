import numpy as np
import pytest

from ls11samples import export, io, vlim
from ls11samples.cosmo import distmod, shell_volume


def _mock(cfg, rng, n=300_000, rlim=19.5):
    """Flux-limited mock: z ~ volume, Mr ~ Gaussian, K = k1 * z with k1 in [0, 2] (95th pct: 1.9)."""
    z = 0.6 * rng.uniform(0, 1, n) ** (1 / 3)
    M = rng.normal(-20.5, 1.2, n)
    k = rng.uniform(0, 2, n) * z
    r = M + distmod(z, cfg) + k
    logm = 10.5 - 0.4 * (M + 20.5) + rng.normal(0, 0.1, n)
    keep = r <= rlim
    return {"BEST_Z": z[keep], "MABS_R": M[keep], "KCORR_R": k[keep], "MAG_R": r[keep],
            "LOGMSTAR": logm[keep], "MAG_G": r[keep] + 0.7, "MAG_Z": r[keep] - 0.4,
            "RA": rng.uniform(0, 5, keep.sum()), "DEC": rng.uniform(-5, 0, keep.sum()),
            "EBV": np.zeros(keep.sum(), np.float32)}


def test_mr_limit_matches_analytic(cfg, rng):
    d = _mock(cfg, rng)
    zc, ml = vlim.mr_limit_curve(d["BEST_Z"], d["KCORR_R"], 19.5, cfg, pct=95)
    expected = 19.5 - distmod(zc, cfg) - 1.9 * zc
    good = (zc > 0.03) & (zc < 0.35)   # observed-K percentile biased beyond (see vlim doc)
    assert np.allclose(ml[good], expected[good], atol=0.05)
    # the volume-limited sample is complete: no galaxy brighter than the threshold is lost
    thr = -21.0
    z1 = vlim.zmax_brighter(zc, ml, thr)
    assert 0.1 < z1 < 0.6
    zt = np.linspace(0.01, z1, 50)
    assert np.all(19.5 - distmod(zt, cfg) - 2.0 * zt >= thr - 0.05)


def test_mstar_limit_monotonic(cfg, rng):
    d = _mock(cfg, rng)
    zc, ml = vlim.mstar_limit_curve(d["BEST_Z"], d["LOGMSTAR"], d["MAG_R"], 19.5)
    assert np.all(np.diff(ml) >= 0)
    z1 = vlim.zmax_heavier(zc, ml, 10.5)
    z2 = vlim.zmax_heavier(zc, ml, 11.0)
    assert z1 < z2


def test_define_and_write_samples(cfg, rng, tmp_path):
    d = _mock(cfg, rng, n=100_000)
    samples, curves = vlim.define_samples(d, cfg, 19.5)
    kinds = {s["kind"] for s in samples}
    assert kinds == {"Mr", "Mstar"} and set(curves) == {"Mr", "Mstar"}
    s = samples[0]
    assert s["name"].startswith("LS11_VLIM_ANY_") and s["name"].endswith(f"_N_{s['n']:07d}")
    assert np.all(d["BEST_Z"][s["sel"]] <= s["z1"])
    rand = vlim.make_randoms({"RA": rng.uniform(0, 5, 10**6), "DEC": rng.uniform(-5, 0, 10**6),
                              "EBV": np.zeros(10**6)}, d["BEST_Z"][s["sel"]], 5, rng)
    assert len(rand["RA"]) == 5 * s["n"]
    files = export.write_sample(tmp_path, s, d, rand, {"AREA": 25.0})
    dt = io.read_table(files[0])
    assert {"RA", "DEC", "BEST_Z", "LPH_MASS_BEST", "WEIGHT_COMP"} <= set(dt)
    assert dt["RA"].dtype == np.float64 and len(dt["RA"]) == s["n"]
    assert set(io.read_table(files[1])) == {"RA", "DEC", "EBV", "Z"}
    col = io.read_table(files[2])
    assert np.array_equal(col["REDSHIFT"], dt["BEST_Z"].astype(np.float32))
    row = vlim.summary_row(s, d, 25.0, cfg)
    assert row["VOLUME"] == pytest.approx(shell_volume(s["z0"], s["z1"], 25.0, cfg))
