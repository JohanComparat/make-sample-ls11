import numpy as np
import pytest

from ls11samples import catalog, io
from ls11samples.config import deep_merge
from ls11samples.env import get_paths


def test_deep_merge_replaces_leaves_and_merges_dicts():
    base = {"a": {"x": 1, "y": [1, 2]}, "b": 2}
    out = deep_merge(base, {"a": {"y": None}, "c": 3})
    assert out == {"a": {"x": 1, "y": None}, "b": 2, "c": 3} and base["a"]["y"] == [1, 2]


def test_product_paths(monkeypatch, tmp_path):
    monkeypatch.setenv("LS11_DIR", str(tmp_path))
    monkeypatch.delenv("LS11_SWEEP_OUT", raising=False)
    monkeypatch.delenv("LS11_OUT", raising=False)
    p = get_paths()
    s = p.sweep_dir / "sweep-000m005-005p000.fits"
    assert p.product(s, "lephare") == tmp_path / "south/sweep/11.0-lephare/sweep-000m005-005p000-lephare.fits"
    assert p.rand_file("bgsr21") == tmp_path / "south/bgsr21/LS11_bgsr21_RAND.fits"


def test_catalog_joins_and_checks_alignment(monkeypatch, tmp_path):
    monkeypatch.setenv("LS11_DIR", str(tmp_path))
    p = get_paths()
    sweep = p.sweep_dir / "sweep-000m005-005p000.fits"
    io.write_table(sweep, {"RA": np.arange(10.0), "DEC": -np.arange(10.0)})
    rows = np.array([1, 4, 7], np.int32)
    io.write_table(p.product(sweep, "sel"), {"LS_ID_DR11": rows.astype(np.int64) + 100, "SWEEP_ROW": rows,
                                             "BEST_Z": np.full(3, 0.1, np.float32)})
    io.write_table(p.product(sweep, "kcorrect"), {"SWEEP_ROW": rows, **{k: np.arange(3, dtype=np.float32)
                                                                        for k in catalog.CODE_COLUMNS}})
    t = catalog.load(p, [sweep], "sel", sweep_columns=["RA"], codes=["kcorrect"])
    assert t["RA"].tolist() == [1.0, 4.0, 7.0] and t["MABS_R_KCORRECT"].tolist() == [0, 1, 2]
    assert t["SWEEP_INDEX"].tolist() == [0, 0, 0]
    io.write_table(p.product(sweep, "cigale"), {"SWEEP_ROW": rows[::-1], **{k: np.zeros(3, np.float32)
                                                                            for k in catalog.CODE_COLUMNS}})
    with pytest.raises(ValueError, match="row-aligned"):
        catalog.load(p, [sweep], "sel", codes=["cigale"])
