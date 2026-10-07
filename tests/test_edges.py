"""Edge cases and error branches (all fast, no real data)."""

import logging
import os
import sys

import healpy as hp
import numpy as np
import pytest

from ls11samples import bgsl, catalog, io, maps, randoms, selection, validate, vlim
from ls11samples.config import config_hash
from ls11samples.env import get_paths


# --------------------------------------------------------------------------- io
def test_sweep_box_rejects_other_names():
    with pytest.raises(ValueError, match="not a sweep"):
        io.sweep_box("catalogue.fits")


def test_strings_and_leftovers(tmp_path):
    path = tmp_path / "t.fits"
    (tmp_path / f"t.fits.{os.getpid()}.part").write_text("crashed writer")     # removed first
    io.write_table(path, {"S": np.array([b"PSF ", b"REX "]), "X": np.arange(2)})
    assert list(io.read_table(path)["S"]) == ["PSF", "REX"]
    assert [list(t["S"]) for _, t in io.iter_rows(path, ["S"], 1)] == [["PSF"], ["REX"]]
    assert io.concat([{}, {"a": np.zeros(0)}]) == {}


def test_read_pz_checks(tmp_path):
    sweep = tmp_path / "11.0" / "sweep-000m005-005p000.fits"
    ids = {"RELEASE": np.full(3, 11010, np.int16), "BRICKID": np.arange(3, dtype=np.int32),
           "OBJID": np.arange(3, dtype=np.int32)}
    io.write_table(sweep, ids)
    with pytest.raises(FileNotFoundError):
        io.read_pz(sweep, np.arange(3), *ids.values())
    pz = io.pz_path(sweep)
    io.write_table(pz, {**{k: v[:2] for k, v in ids.items()}, "LS_ID_DR11": io.ls_id(11010, ids["BRICKID"][:2], ids["OBJID"][:2])})
    with pytest.raises(ValueError, match="different numbers of rows"):
        io.read_pz(sweep, np.arange(3), *ids.values())
    io.write_table(pz, {**ids, "LS_ID_DR11": io.ls_id(11010, ids["BRICKID"], ids["OBJID"][::-1])})
    with pytest.raises(ValueError, match="not row-matched"):
        io.read_pz(sweep, np.arange(3), *ids.values())


def test_catalog_empty_rows(tmp_path):
    assert catalog.sweep_rows(tmp_path / "none.fits", np.zeros(0, int), ["RA"])["RA"].size == 0


# --------------------------------------------------------------------------- bits / selection / bgsl
def test_header_without_bits():
    from ls11samples import bits

    with pytest.raises(ValueError, match="no MBIT_n"):
        bits.check_header({"SIMPLE": True})


def test_cutflow_skips_absent_cuts():
    flow = selection.cutflow({"nobs": np.array([True, False]), "rmag": np.array([True, True])})
    assert flow["CUT"].tolist() == ["nobs", "rmag"] and flow["N_PASS_CUMUL"].tolist() == [1, 1]


def test_existing_selection_kept(cfg, tmp_path, monkeypatch):
    monkeypatch.setenv("LS11_DIR", str(tmp_path))
    p = get_paths()
    sweep = p.sweep_dir / "sweep-000m005-005p000.fits"
    out = p.product(sweep, cfg["tag"])
    io.write_table(out, {"SWEEP_ROW": np.arange(2)}, header={"CFGHASH": config_hash(cfg)})
    assert bgsl.check_existing(out, cfg, overwrite=False)
    assert bgsl.run_sweep(sweep, cfg, p) == out                    # kept, the sweep is never read


# --------------------------------------------------------------------------- randoms / maps
def test_boxset_off_grid_and_density_warning(cfg, tmp_path, caplog):
    b = randoms.BoxSet([(1.0, 2.5, -3.0, -1.0)])
    assert not b.grid and b.contains(np.array([2.0, 3.0]), np.array([-2.0, -2.0])).tolist() == [True, False]
    path = tmp_path / "r.fits"
    io.write_table(path, {"RA": np.zeros(1)}, header={"DENSITY": 100.0})
    with caplog.at_level(logging.WARNING):
        assert randoms.file_density(path, cfg) == 100.0
    assert "differs" in caplog.text


def test_subsample_more_than_available(caplog):
    rng = np.random.default_rng(0)
    with caplog.at_level(logging.WARNING):
        idx = randoms.subsample(5, 8, rng)
    assert sorted(idx.tolist()) == [0, 1, 2, 3, 4] and "only 5 available" in caplog.text
    z = randoms.shuffle_z(np.array([0.1, 0.2]), 50, rng, weights=[1, 0])
    assert np.all(z == 0.1)


def test_maps_missing_quantity_and_gaia(cfg, tmp_path, caplog):
    cfg["maps"]["nsides"], cfg["maps"]["quantities"] = [4], ["EBV", "NOT_A_COLUMN"]
    cfg["maps"]["gaia_templates"] = ["nstar_faint"]
    acc = maps.MapAccumulator(cfg)
    acc.add({"RA": np.array([10.0, 200.0]), "DEC": np.array([0.0, 30.0]), "EBV": np.array([0.01, 0.02])})
    with caplog.at_level(logging.WARNING):
        written = acc.write(tmp_path, 1.0, gaia_dir=tmp_path / "no_gaia")
    assert "no Gaia map" in caplog.text
    ebv = hp.read_map(str(tmp_path / "systematics" / "0004" / "LS11_EBV_NSIDE_0004.fits"))
    assert np.sort(ebv[ebv != hp.UNSEEN]).tolist() == pytest.approx([0.01, 0.02]) and len(written) == 3


# --------------------------------------------------------------------------- vlim / plots / kcorr
def test_vlim_curve_edges(cfg):
    zc, out = vlim._binned_percentile(np.linspace(0, 0.6, 600), np.ones(600), np.array([0, 0.3, 0.6]), 50,
                                      mask=np.arange(600) < 300)
    assert out[0] == 1 and np.isnan(out[1])                         # masked out: too few
    assert vlim.zmax_brighter(np.array([0.1, 0.2]), np.array([-21.0, -19.0]), -20.0) == 0.1  # first bin short
    assert np.isnan(vlim.zmax_brighter(np.array([0.1]), np.array([-23.0]), -22.0))     # never complete
    assert vlim.zmax_brighter(np.array([0.1]), np.array([-23.0]), -24.0) == 0.1   # always complete
    data = {"BEST_Z": np.full(50, 0.5), "MABS_R": np.full(50, -19.0), "LOGMSTAR": np.full(50, 9.0),
            "MAG_R": np.full(50, 19.0)}
    cfg["vlim"]["absmag_r"]["z_min"] = 0.6                          # above every z_max: no sample
    samples, _ = vlim.define_samples(data, cfg, 19.5, mr_curve=(np.array([0.1, 0.5]), np.array([-30.0, -30.0])))
    assert [s for s in samples if s["kind"] == "Mr"] == []


def test_plots_without_masses(cfg, tmp_path):
    pytest.importorskip("matplotlib")
    from ls11samples.plots import vlim_planes

    data = {"BEST_Z": np.linspace(0.05, 0.5, 100), "MABS_R": np.linspace(-23, -18, 100)}
    path = vlim_planes(data, [{"kind": "Mr", "lo": -25.0, "hi": -20.0, "z0": 0.05, "z1": 0.2}],
                       {"Mr": (np.array([0.1, 0.3]), np.array([-19.0, -21.0]))}, tmp_path / "p.png")
    assert path.stat().st_size > 1000


def test_kcorr_edges(cfg):
    pytest.importorskip("kcorrect")
    from ls11samples import kcorr
    from ls11samples.sed.common import BANDS

    k = kcorr.KcorrectV5(cfg, mc=0)
    t = {"BEST_Z": np.array([np.nan, 2.5])}
    for b in BANDS:
        t[f"FLUX_{b}"], t[f"FLUX_IVAR_{b}"], t[f"MW_TRANSMISSION_{b}"] = np.ones(2), np.ones(2), np.ones(2)
    assert np.isnan(k.fit(t)["KC_LOGMSTAR"]).all()                 # no usable redshift
    z = np.linspace(0.02, 0.1, 500)
    idx = kcorr.complete_set_rows(z, np.full(500, -25.0), np.full(500, 14.0), 19.5, cfg, n_max=50)
    assert idx.size == 50 and np.all(np.diff(idx) > 0)
    from ls11samples.sed.kcorrect_fit import KcorrectBackend

    assert KcorrectBackend(cfg).fit({"BEST_Z": np.zeros(0)})["LOGMSTAR"].size == 0


# --------------------------------------------------------------------------- validate / SED helpers
def test_validate_missing_and_empty(cfg, tmp_path, monkeypatch):
    monkeypatch.setenv("LS11_DIR", str(tmp_path))
    p = get_paths()
    s1, s2 = p.sweep_dir / "sweep-000m005-005p000.fits", p.sweep_dir / "sweep-005m005-010p000.fits"
    io.write_table(p.product(s2, cfg["tag"]), {"SWEEP_ROW": np.zeros(0, np.int32), "STAR_FLAG": np.zeros(0, np.uint8),
                                               "BEST_Z": np.zeros(0, np.float32)}, header={"NSEL": 0})
    out = validate.validate(p, cfg, [s1, s2], codes=["kcorrect"])
    assert out["missing"] == [[cfg["tag"], s1.name]] and out["n_empty_sweeps"] == 1
    assert "n/a" in validate.report({**out, "vs_ref": {"a-b": {"dlogm_median": None, "dlogm_nmad": float("nan"),
                                                                "dmr_median": None, "dmr_nmad": 0.1}}})


def test_cigale_run_wrapper(monkeypatch, tmp_path):
    from ls11samples.sed import cigale_fit

    assert cigale_fit._run([sys.executable, "-c", "print('ok')"]).strip() == "ok"
    with pytest.raises(RuntimeError, match="failed"):
        cigale_fit._run([sys.executable, "-c", "import sys; sys.exit(3)"])
    names = " ".join(cigale_fit.FILTER_NAMES.values())
    monkeypatch.setattr(cigale_fit, "_run", lambda args, cwd=None: names)
    cigale_fit.ensure_filters()                                     # all registered: nothing to do
    monkeypatch.setenv("TMPDIR", str(tmp_path))
    b = cigale_fit.CigaleBackend({"sed": {}})
    assert b.workroot == tmp_path / "ls11_cigale" and b.workroot.is_dir()


def test_filters_need_kcorrect(monkeypatch):
    import importlib.util

    from ls11samples.sed import filters

    monkeypatch.setattr(importlib.util, "find_spec", lambda name: None)
    with pytest.raises(ImportError, match="kcorrect"):
        filters.responses_dir()


def test_git_commit_outside_a_checkout(monkeypatch):
    import subprocess

    from ls11samples import export

    def fail(*a, **k):
        raise subprocess.CalledProcessError(128, "git")

    monkeypatch.setattr(subprocess, "run", fail)
    assert export.git_commit() == "unknown"
