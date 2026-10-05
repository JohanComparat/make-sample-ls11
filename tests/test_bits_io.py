import numpy as np
import pytest

from ls11samples import bits, io


def _header(mb=bits.MASKBITS, fb=bits.FITBITS):
    h = {f"MBIT_{k}": v for k, v in mb.items()}
    h.update({f"FBIT_{k}": v for k, v in fb.items()})
    return h


def test_check_header_ok_and_mismatch():
    bits.check_header(_header())
    wrong = dict(bits.MASKBITS)
    wrong[12], wrong[13] = "CLUSTER", "GALAXY"
    with pytest.raises(ValueError, match="MBIT"):
        bits.check_header(_header(mb=wrong))


def test_mask_value():
    assert bits.mask_value(["BRIGHT", "GALAXY", "CLUSTER"], bits.MASKBIT) == 2 + 4096 + 8192
    with pytest.raises(KeyError):
        bits.mask_value(["NOPE"], bits.MASKBIT)


def test_ls_id_and_names(tmp_path):
    assert io.ls_id(11010, 5, 7) == (11010 << 42) | (5 << 22) | 7
    s = tmp_path / "11.0" / "sweep-000m005-005p000.fits"
    assert io.pz_path(s) == tmp_path / "11.0-photo-z" / "sweep-000m005-005p000-pz.fits"
    assert io.sweep_box(s) == (0.0, 5.0, -5.0, 0.0)
    assert io.box_area([(0, 360, -90, 90)]) == pytest.approx(41252.96, rel=1e-5)


def test_write_read_roundtrip(tmp_path):
    t = {"A": np.arange(5, dtype=np.int64), "B": np.linspace(0, 1, 5).astype(np.float32),
         "S": np.array(["PSF", "REX", "DUP", "SER", "EXP"]), "V": np.ones((5, 3))}
    f = io.write_table(tmp_path / "x.fits", t, header={"KEY": 3}, extname="X",
                       extra=[("Y", {"C": np.arange(2)}, None)])
    r = io.read_table(f)
    assert np.array_equal(r["A"], t["A"]) and list(r["S"]) == list(t["S"]) and r["V"].shape == (5, 3)
    assert io.read_header(f, "X")["KEY"] == 3
    assert io.read_table(f, ext="Y")["C"].tolist() == [0, 1]
    chunks = list(io.iter_rows(f, ["A"], chunk=2))
    assert [c[0] for c in chunks] == [0, 2, 4]
