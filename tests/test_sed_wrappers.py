"""SED backends without running the heavy codes: LePhare and CIGALE through stand-ins (a fake
``lephare`` module, a fake ``pcigale`` command line), the shared helpers and filters, and DSPS / eazy
when installed (``LS11_SLOW_TESTS=1`` for eazy).
"""

import os
import sys
import types
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table

from ls11samples.sed import OUTPUT, PRODUCT, common, get_backend
from ls11samples.sed.common import BANDS, NANOMAGGY_CGS, NANOMAGGY_MJY

PIVOT = {"G": 4810.0, "R": 6420.0, "I": 7810.0, "Z": 9170.0, "W1": 33680.0, "W2": 46180.0}


def galaxies(n=8, rng=None):
    """Input table of ``n`` galaxies: dereddened fluxes 10..80 nmgy, one row without redshift."""
    rng = rng or np.random.default_rng(3)
    t = {"BEST_Z": np.linspace(0.05, 0.6, n).astype(np.float32)}
    t["BEST_Z"][-1] = np.nan
    for b in BANDS:
        t[f"FLUX_{b}"] = rng.uniform(10, 80, n).astype(np.float32)
        t[f"FLUX_IVAR_{b}"] = np.full(n, 4.0, np.float32)
        t[f"MW_TRANSMISSION_{b}"] = np.ones(n, np.float32)
    t["FLUX_IVAR_W2"][0] = 0.0                                     # one missing band
    return t


# --------------------------------------------------------------------------- common / filters
def test_common_photometry(cfg):
    t = galaxies()
    z, f, e = common.photometry(t, cfg)
    assert f.shape == e.shape == (8, 6) and np.isinf(e[0, 5])
    floor = 0.4 * np.log(10) * 0.02 * f[1, 0]
    assert e[1, 0] == pytest.approx(np.hypot(0.5, floor))
    assert common.valid_z(z).tolist() == [True] * 7 + [False]
    assert set(common.empty_result(3)) == set(OUTPUT) and set(PRODUCT) <= set(OUTPUT)


def test_get_backend_names():
    for code in ("kcorrect", "lephare", "cigale", "eazy", "dsps"):
        assert get_backend(code).name == code


def test_filters(tmp_path):
    pytest.importorskip("kcorrect")
    from ls11samples.sed import filters

    for b, lo, hi in (("G", 4700, 4900), ("R", 6300, 6500), ("Z", 9000, 9300), ("W1", 33000, 34500)):
        wl, tr = filters.read_response(b)
        assert np.all(np.diff(wl) > 0) and tr.min() >= 0 and lo < filters.pivot_wavelength(b) < hi
    files = filters.write_lephare_filters(tmp_path)
    assert [f.name for f in files] == ["decam_g.pb", "decam_r.pb", "decam_i.pb", "decam_z.pb",
                                       "wise_w1.pb", "wise_w2.pb"]
    assert np.loadtxt(files[1]).shape[1] == 2


# --------------------------------------------------------------------------- LePhare
class _Keyword:
    def __init__(self, name, value):
        self.name, self.value = name, value


def fake_lephare(calls):
    """A stand-in for the ``lephare`` module recording the calls of the backend."""
    dr = types.SimpleNamespace(
        config_to_required_files=lambda config: ["opa/OPACITY.dat", "filt/ls11/decam_g.pb", "ext/SB_calzetti.dat"],
        download_registry_from_github=lambda outfile: calls.append(("registry", outfile)),
        make_retriever=lambda registry_file, data_path: ("retriever", data_path),
        download_all_files=lambda retriever, files, ignore_registry: calls.append(("download", files)))

    def process(config, tab):
        n = len(tab)
        calls.append(("process", n))
        mass = np.linspace(9.0, 11.0, n)
        res = Table({"MASS_MED": mass, "MASS_INF": mass - 0.1, "MASS_SUP": mass + 0.2,
                     "MASS_BEST": mass + 0.01, "SFR_MED": np.zeros(n), "CHI_BEST": np.ones(n)})
        res["MASS_MED"][0] = -99.0                                 # LePhare sentinel -> MASS_BEST
        mabs = np.full((n, 6), -21.0)
        mabs[1, 1] = -99.9                                         # no absolute magnitude
        res["MAG_ABS()"] = mabs
        res["EMAG_ABS()"] = np.full((n, 6), 40.0)                  # m - M in LePhare 1.0, unused
        return res, []

    return types.SimpleNamespace(
        LEPHAREDIR="/nonexistent", keyword=_Keyword, data_retrieval=dr, process=process,
        read_config=lambda path: {"GAL_LIB_OUT": _Keyword("GAL_LIB_OUT", "BC03_LS11")},
        prepare=lambda config: calls.append(("prepare", sorted(config))))


@pytest.fixture
def lephare_env(monkeypatch, tmp_path):
    from ls11samples.sed import lephare_fit

    calls = []
    monkeypatch.setitem(sys.modules, "lephare", fake_lephare(calls))
    monkeypatch.setenv("LEPHAREDIR", str(tmp_path / "data"))
    monkeypatch.setenv("LEPHAREWORK", str(tmp_path / "work"))
    monkeypatch.setattr(lephare_fit, "write_lephare_filters", lambda d: calls.append(("filters", str(d))))
    monkeypatch.setattr(lephare_fit, "pivot_wavelength", lambda b: PIVOT[b])
    (tmp_path / "data" / "opa").mkdir(parents=True)
    (tmp_path / "data" / "opa" / "OPACITY.dat").write_text("present")
    return calls, tmp_path


def test_lephare_setup(cfg, lephare_env):
    from ls11samples.sed.lephare_fit import LephareBackend

    calls, tmp = lephare_env
    b = LephareBackend(cfg)
    assert b.config["FILTER_REP"].value == str(tmp / "data" / "filt")
    assert b.config["PARA_OUT"].value.endswith("output_ls11.para")
    downloads = [c for c in calls if c[0] == "download"]
    assert downloads == [("download", ["ext/SB_calzetti.dat"])]  # present and ls11 filters skipped
    assert [c[0] for c in calls].count("prepare") == 1           # no library yet -> built
    doc = tmp / "work" / "lib_mag" / "BC03_LS11.doc"
    doc.parent.mkdir(parents=True)
    doc.write_text("library ready")
    calls.clear()
    LephareBackend(cfg)
    assert "prepare" not in [c[0] for c in calls]                # finished library: reused
    LephareBackend(cfg, rebuild=True)
    assert [c[0] for c in calls].count("prepare") == 1
    doc.write_text("")                                           # interrupted build: rebuilt
    calls.clear()
    LephareBackend(cfg, build=False)
    assert "prepare" not in [c[0] for c in calls]


def test_lephare_table_and_fit(cfg, lephare_env):
    from ls11samples.sed.lephare_fit import LephareBackend

    calls, _ = lephare_env
    b = LephareBackend(cfg, build=False)
    t = galaxies()
    tab, ok = b.table(t)
    assert tab.colnames[:3] == ["id", "f0", "e0"] and ok.sum() == 7
    assert tab["f0"][1] == pytest.approx(t["FLUX_G"][1] * NANOMAGGY_CGS)
    assert tab["e5"][0] == -99.0 and tab["context"][0] == 31 and tab["context"][1] == 63
    out = b.fit(t)
    assert ("process", 7) in calls
    assert out["LOGMSTAR"][0] == pytest.approx(9.01)                 # MASS_MED sentinel -> MASS_BEST
    assert out["LOGMSTAR_ERR"][2] == pytest.approx(0.15)
    assert np.isnan(out["MABS_R"][1]) and np.isnan(out["MABS_R_ERR"][1])
    assert np.isnan(out["LOGMSTAR"][-1])                             # no redshift: not fitted
    # MABS_R_ERR: error of the observed band nearest to rest-frame r: r at z=0.05, i at z~0.21
    z, f, e = common.photometry(t, cfg)
    for i, band in ((0, 1), (3, 2)):
        lam = np.array([PIVOT[x] for x in BANDS])
        assert np.argmin(np.abs(lam - PIVOT["R"] * (1 + z[i]))) == band
        assert out["MABS_R_ERR"][i] == pytest.approx(2.5 / np.log(10) * e[i, band] / f[i, band], rel=1e-5)
    assert b.fit({k: v[-1:] for k, v in t.items()})["LOGMSTAR"].tolist() == [pytest.approx(np.nan, nan_ok=True)]


# --------------------------------------------------------------------------- CIGALE
def fake_pcigale(calls, existing_filters=""):
    """Stand-in for ``_run``: the pcigale command line and the filter-database one-liners."""
    from configobj import ConfigObj

    def run(args, cwd=None):
        calls.append([str(a) for a in args[1:2]] if args[1] == "-c" else args[1])
        if args[1] == "-c":
            return existing_filters if "SimpleDatabase" in args[2] else ""
        cwd = Path(cwd)
        if args[1] == "init":
            ConfigObj({"data_file": "", "sed_modules": [], "analysis_method": "", "cores": 1}).write(
                open(cwd / "pcigale.ini", "wb"))
        elif args[1] == "genconf":
            c = ConfigObj(str(cwd / "pcigale.ini"), encoding="UTF8")
            c["additionalerror"] = 0.1                             # what genconf writes
            c["sed_modules_params"] = {m: {} for m in c["sed_modules"]}
            c["analysis_params"] = {"variables": [], "bands": [], "save_best_sed": True}
            c.write()
        elif args[1] == "run":
            inp = np.loadtxt(cwd / "input.txt")
            n = len(inp)
            m = 10 ** np.linspace(9.5, 11.0, n)
            lnu = 1e21 * np.ones(n)
            (cwd / "out").mkdir()
            Table({"id": inp[:, 0].astype(int), "bayes.stellar.m_star": m, "bayes.stellar.m_star_err": 0.1 * m,
                   "bayes.sfh.sfr": np.ones(n), "bayes.param.restframe_Lnu(ls11.decam_r)": lnu,
                   "bayes.param.restframe_Lnu(ls11.decam_r)_err": 0.05 * lnu,
                   "best.chi_square": np.full(n, 2.0)}).write(cwd / "out" / "results.fits")
        return ""

    return run


@pytest.fixture
def cigale_env(monkeypatch):
    pytest.importorskip("configobj")
    from ls11samples.sed import cigale_fit

    calls = []
    monkeypatch.setattr(cigale_fit, "_run", fake_pcigale(calls, existing_filters="wise.W1 ls11.decam_g"))
    monkeypatch.setattr(cigale_fit, "read_response", lambda b: (np.array([1.0, 2.0, 3.0]), np.array([0.0, 1.0, 0.0])))
    return calls


def test_cigale_filters_registered(cigale_env):
    from ls11samples.sed import cigale_fit

    cigale_fit.ensure_filters()
    assert len([c for c in cigale_env if c == ["-c"]]) == 2           # list + add the 5 missing ones


def test_cigale_fit(cfg, cigale_env, tmp_path):
    from configobj import ConfigObj
    from ls11samples.sed.cigale_fit import TO_LUMIN, CigaleBackend

    b = CigaleBackend(cfg, workdir=tmp_path, cores=2, keep_runs=True)
    t = galaxies()
    out = b.fit(t)
    conf = ConfigObj(str(b.last_run / "pcigale.ini"), encoding="UTF8")
    assert float(conf["additionalerror"]) == 0.0 and conf["cores"] == "2"
    assert conf["sed_modules_params"]["sfhdelayed"]["tau_main"] == ["250", "500", "1000", "2000", "4000", "8000"]
    assert conf["sed_modules_params"]["nebular"]["logU"] == "-2.0"
    assert "param.restframe_Lnu(ls11.decam_r)" in conf["analysis_params"]["variables"]
    inp = np.loadtxt(b.last_run / "input.txt")
    assert inp.shape == (7, 14) and inp[0, 2] == pytest.approx(t["FLUX_G"][0] * NANOMAGGY_MJY, rel=1e-5)
    assert np.isnan(inp[0, 13])                                        # missing W2 error
    assert out["LOGMSTAR"][:7] == pytest.approx(np.linspace(9.5, 11.0, 7))
    assert out["LOGMSTAR_ERR"][:7] == pytest.approx(0.1 / np.log(10))
    assert out["MABS_R"][0] == pytest.approx(-2.5 * np.log10(1e21 / TO_LUMIN / 3.631e6))
    assert out["MABS_R_ERR"][0] == pytest.approx(2.5 / np.log(10) * 0.05)
    assert np.isnan(out["LOGMSTAR"][-1])
    b2 = CigaleBackend(cfg, workdir=tmp_path / "w2")
    b2.fit(t)
    assert not list((tmp_path / "w2").glob("run_*"))                  # removed after reading
    assert np.isnan(b2.fit({k: v[-1:] for k, v in t.items()})["LOGMSTAR"]).all()


# --------------------------------------------------------------------------- DSPS / eazy
def test_dsps_helpers():
    from ls11samples.sed.dsps_fit import _weighted_percentiles, age_weights, surviving_fraction

    assert surviving_fraction(0.0) == pytest.approx(1.0) and surviving_fraction(10.0) < surviving_fraction(1.0)
    lg_age = np.linspace(-4, 1.3, 107)
    w, fsurv, sfr = age_weights(lg_age, tau=1.0, age=1.0)
    assert w.sum() == pytest.approx(1.0) and 0.6 < fsurv < 1
    assert sfr == pytest.approx(np.exp(-1) / (1 - 2 / np.e) / 1e9, rel=1e-3)   # delayed tau, analytic
    v = np.tile(np.arange(101.0), (2, 1))
    lo, med, hi = _weighted_percentiles(v, np.ones_like(v), (16, 50, 84))
    assert med.tolist() == [50.0, 50.0] and lo[0] == pytest.approx(16, abs=1) and hi[0] == pytest.approx(84, abs=1)
    assert np.isnan(_weighted_percentiles(v, np.zeros_like(v), (50,))[0]).all()


def test_dsps_backend(cfg):
    pytest.importorskip("dsps")
    from ls11samples.sed.dsps_fit import SSP_FILE, DspsBackend

    drn = Path(os.environ.get("DSPS_DRN", "~/.cache/dsps")).expanduser()
    if not (drn / SSP_FILE).exists():
        pytest.skip("no DSPS SSP file")
    b = DspsBackend(cfg, grid={"tau_gyr": [1.0], "age_gyr": [3.0, 9.0], "z_over_zsun": [1.0], "av": [0.0, 0.3]},
                    zgrid=np.arange(0.05, 0.7, 0.05))
    out = b.fit(galaxies())
    assert np.isfinite(out["LOGMSTAR"][:7]).all() and np.isnan(out["LOGMSTAR"][-1])
    assert np.all(out["LOGMSTAR_LO"][:7] <= out["LOGMSTAR_HI"][:7])


@pytest.mark.skipif(os.environ.get("LS11_SLOW_TESTS") != "1", reason="slow: set LS11_SLOW_TESTS=1")
def test_eazy_backend(cfg, tmp_path):
    pytest.importorskip("eazy")
    from ls11samples.sed.eazy_fit import EazyBackend

    out = EazyBackend(cfg, workdir=tmp_path, n_proc=2).fit(galaxies())
    assert np.isfinite(out["LOGMSTAR"][:7]).all() and np.isnan(out["LOGMSTAR"][-1])
