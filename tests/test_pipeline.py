"""End-to-end run of the pipeline scripts on the synthetic DR11 tree of ``fake_dr11``.

The scripts run in-process (runpy) so that coverage sees them; the tests run in file order and
share one tree (module-scoped fixture): select -> randoms + maps -> kcorrect fits -> volume-limited
samples -> export checks -> validation, plus the skip / overwrite / stale branches.
"""

import runpy
import sys
from pathlib import Path

import numpy as np
import pytest
import yaml

import fake_dr11
from ls11samples import io
from ls11samples.config import DEFAULT, load_config
from ls11samples.env import get_paths

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"


def run_script(name: str, *args: str) -> int:
    """Run scripts/<name> as __main__ with ``args``; returns its exit status."""
    argv = sys.argv
    sys.argv = [name, *args]
    try:
        runpy.run_path(str(SCRIPTS / name), run_name="__main__")
        return 0
    except SystemExit as e:
        return e.code if isinstance(e.code, int) else (0 if e.code is None else 1)
    finally:
        sys.argv = argv


@pytest.fixture(scope="module")
def tree(tmp_path_factory):
    root = tmp_path_factory.mktemp("ls11")
    dr11 = fake_dr11.build(root)
    cfg_file = root / "test.yaml"
    cfg_file.write_text(fake_dr11.CONFIG.format(default=DEFAULT, density=fake_dr11.DENSITY))
    mp = pytest.MonkeyPatch()
    for k in ("LS11_SWEEPS", "LS11_SWEEP_LIST", "LS11_OUT", "LS11_SWEEP_OUT", "LS11_RANDOMS"):
        mp.delenv(k, raising=False)
    mp.setenv("LS11_DIR", str(dr11))
    mp.setenv("LS11_CONFIG", str(cfg_file))
    mp.setenv("LS11_GAIA_MAPS", str(root / "gaia"))
    mp.setenv("LS11_NPROC", "1")
    yield {"root": root, "dr11": dr11, "cfg_file": cfg_file}
    mp.undo()


def test_select(tree):
    assert run_script("01_select.py", "--nproc", "1") == 0
    p, cfg = get_paths(), load_config()
    sweeps = p.sweeps()
    assert [s.name for s in sweeps] == list(fake_dr11.SWEEPS[:2])        # third: no photo-z
    for s in sweeps:
        f = p.product(s, cfg["tag"])
        sel = io.read_table(f)
        h = io.read_header(f, 1)
        flow = io.read_table(f, ext="CUTFLOW")
        assert h["NSEL"] == len(sel["SWEEP_ROW"]) == flow["N_PASS_CUMUL"][-1] > 1000
        assert set(sel) == {"LS_ID_DR11", "SWEEP_ROW", "BEST_Z", "BEST_Z_ERR", "Z_SOURCE", "STAR_FLAG"}
        assert (sel["STAR_FLAG"] > 0).any() and (sel["Z_SOURCE"] == 1).any()
        assert "MEDIUM(11)" in h["MASKREJ"]
    # existing outputs are kept (also through the process pool), another configuration is refused
    assert run_script("01_select.py", "--nproc", "2") == 0
    other = tree["root"] / "other.yaml"
    other.write_text(tree["cfg_file"].read_text().replace("[13.0, 21.0]", "[13.0, 20.0]"))
    with pytest.raises(RuntimeError, match="another configuration"):
        run_script("01_select.py", "--config", str(other))


def test_select_needs_sweeps(tree, monkeypatch):
    monkeypatch.setenv("LS11_SWEEPS", "sweep-nothing*.fits")
    assert run_script("01_select.py") != 0


def test_randoms_and_maps(tree, monkeypatch):
    assert run_script("02_randoms.py") == 0
    p, cfg = get_paths(), load_config()
    rf = p.rand_file(cfg["tag"])
    h = io.read_header(rf, 1)
    r = io.read_table(rf)
    assert set(r) == {"RA", "DEC", "EBV"} and h["NRAND"] == len(r["RA"])
    assert 0.6 * 2 * 24.5 < h["AREA"] < 2 * 25.0                       # two 5x5 deg boxes, masked
    assert h["NSWEEPS"] == 2
    run = p.run_dir(cfg["tag"])
    for n in (8, 16):
        assert (run / "footprint" / f"LS11_FRACAREA_NSIDE_{n:04d}.fits").exists()
        names = {f.name for f in (run / "systematics" / f"{n:04d}").glob("*.fits")}
        assert f"LS11_GALDEPTH_R_NSIDE_{n:04d}.fits" in names and f"GAIA_nstar_faint_NSIDE_{n:05d}.fits" in names
    mtime = rf.stat().st_mtime
    assert run_script("02_randoms.py", "--no-maps") == 0                 # up to date: kept
    assert rf.stat().st_mtime == mtime
    # another sweep list -> the random file is stale and redone
    lst = tree["root"] / "one_sweep.txt"
    lst.write_text(fake_dr11.SWEEPS[0] + "\n")
    monkeypatch.setenv("LS11_SWEEP_LIST", str(lst))
    assert run_script("02_randoms.py", "--no-maps") == 0
    assert io.read_header(rf, 1)["NSWEEPS"] == 1
    monkeypatch.delenv("LS11_SWEEP_LIST")
    assert run_script("02_randoms.py", "--no-maps") == 0
    assert io.read_header(rf, 1)["NSWEEPS"] == 2


def test_randoms_need_selection(tree, monkeypatch):
    monkeypatch.setenv("LS11_RANDOMS", "randoms-none-*.fits")
    assert run_script("02_randoms.py", "--overwrite") != 0                # no random file


def test_fit_kcorrect(tree):
    pytest.importorskip("kcorrect")
    assert run_script("04_stellar_mass.py", "--code", "kcorrect", "--fetch") == 0
    assert run_script("04_stellar_mass.py", "--code", "kcorrect", "--chunk", "1000") == 0
    p, cfg = get_paths(), load_config()
    for s in p.sweeps():
        sel = io.read_table(p.product(s, cfg["tag"]))
        c = io.read_table(p.product(s, "kcorrect"))
        assert np.array_equal(c["SWEEP_ROW"], sel["SWEEP_ROW"])
        assert set(c) == {"SWEEP_ROW", "LOGMSTAR", "LOGMSTAR_ERR", "MABS_R", "MABS_R_ERR"}
        gal = (sel["STAR_FLAG"] == 0) & (sel["BEST_Z"] < 1)
        assert np.isfinite(c["LOGMSTAR"][gal]).mean() > 0.95
        assert np.nanmedian(c["MABS_R"][gal]) < -18
        assert io.read_header(p.product(s, "kcorrect"), 1)["SELHASH"] == io.read_header(p.product(s, cfg["tag"]), 1)["CFGHASH"]
    mtime = p.product(p.sweeps()[0], "kcorrect").stat().st_mtime
    assert run_script("04_stellar_mass.py", "--code", "kcorrect") == 0  # up to date: kept
    assert p.product(p.sweeps()[0], "kcorrect").stat().st_mtime == mtime


def test_vlim_and_export(tree):
    pytest.importorskip("kcorrect")
    assert run_script("05_vlim.py") == 0
    p, cfg = get_paths(), load_config()
    vdir = p.run_dir(cfg["tag"], "vlim")
    manifest = yaml.safe_load((vdir / "manifest.yaml").read_text())
    kinds = {s["KIND"] for s in manifest["samples"]}
    assert kinds == {"Mr", "Mstar"}
    name = manifest["samples"][0]["NAME"]
    d = io.read_table(vdir / f"{name}_DATA.fits")
    r = io.read_table(vdir / f"{name}_RAND.fits")
    assert len(r["RA"]) == 2 * len(d["RA"]) and "LPH_MASS_BEST" in d and "LOGMSTAR_KCORRECT" in d
    assert (vdir / "vlim_boundaries.fits").exists() and (vdir / "vlim_planes.png").exists()
    assert run_script("06_export.py") == 0
    manifest = yaml.safe_load((vdir / "manifest.yaml").read_text())
    assert manifest["checks"]["errors"] == [] and 8 in manifest["sys_mapping"]["template_dirs"]
    # a broken sample is reported and the step fails
    (vdir / f"{name}_COLOUR.fits").unlink()
    assert run_script("06_export.py") == 1


def test_vlim_needs_masses(tree):
    assert run_script("05_vlim.py", "--code", "lephare") != 0            # no lephare outputs


def test_validate(tree):
    pytest.importorskip("kcorrect")
    assert run_script("validate_run.py", "--frac", "0.5") == 0
    p, cfg = get_paths(), load_config()
    out = yaml.safe_load((p.run_dir(cfg["tag"]) / "validation.yaml").read_text())
    assert out["n_sweeps"] == 2 and out["missing"] == [] and out["misaligned"] == []
    assert out["randoms"]["n_sweeps"] == 2 and set(out["sizes_gb"]) == {"testsel", "kcorrect"}
    assert run_script("validate_run.py", "--codes", "kcorrect,cigale") == 1  # cigale never run


def test_failure_branches(tree, monkeypatch):
    """Each step refuses inputs that do not belong together; step 5 survives without kcorrect and
    matplotlib. Runs last: it alters the tree."""
    pytest.importorskip("kcorrect")
    import shutil

    import fitsio

    p, cfg = get_paths(), load_config()
    s0, _, s2 = (p.sweep_dir / n for n in fake_dr11.SWEEPS)
    # step 2: a sweep (now with its photo-z) has no selection file
    shutil.copy(io.pz_path(s0), io.pz_path(s2))
    assert run_script("02_randoms.py") != 0
    io.pz_path(s2).unlink()
    # step 2: the selection was made with another footprint
    other = tree["root"] / "footprint.yaml"
    other.write_text(tree["cfg_file"].read_text() + "footprint:\n  ebv_max: 0.08\n")
    assert run_script("02_randoms.py", "--config", str(other)) != 0
    # step 5: randoms made from another set of sweeps
    lst = tree["root"] / "first.txt"
    lst.write_text(fake_dr11.SWEEPS[0] + "\n")
    monkeypatch.setenv("LS11_SWEEP_LIST", str(lst))
    assert run_script("05_vlim.py") != 0
    monkeypatch.delenv("LS11_SWEEP_LIST")
    # step 5 without kcorrect (observed K percentiles) and without matplotlib (no figure)
    monkeypatch.setitem(sys.modules, "kcorrect.kcorrect", None)
    monkeypatch.setitem(sys.modules, "matplotlib", None)
    vdir = p.run_dir(cfg["tag"], "vlim")
    (vdir / "vlim_planes.png").unlink()
    assert run_script("05_vlim.py") == 0
    assert not (vdir / "vlim_planes.png").exists()
    monkeypatch.undo()
    # step 4: a code file of another selection
    with fitsio.FITS(str(p.product(s0, "kcorrect")), "rw") as h:
        h[1].write_key("SELHASH", "0123456789ab")
    with pytest.raises(RuntimeError, match="another selection"):
        run_script("04_stellar_mass.py", "--code", "kcorrect")
    # step 6: wrong type, non-finite redshift, COLOUR not aligned with DATA
    name = yaml.safe_load((vdir / "manifest.yaml").read_text())["samples"][0]["NAME"]
    d = io.read_table(vdir / f"{name}_DATA.fits")
    d["RA"] = d["RA"].astype(np.float32)
    d["BEST_Z"][0] = np.nan
    io.write_table(vdir / f"{name}_DATA.fits", d)
    c = io.read_table(vdir / f"{name}_COLOUR.fits")
    io.write_table(vdir / f"{name}_COLOUR.fits", io.take(c, slice(1, None)))
    assert run_script("06_export.py") == 1
    errors = " ".join(yaml.safe_load((vdir / "manifest.yaml").read_text())["checks"]["errors"])
    assert "not float64" in errors and "non-finite" in errors and "not row-aligned" in errors
