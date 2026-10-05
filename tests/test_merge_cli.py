import numpy as np
import pytest

from ls11samples import cli, io
from ls11samples.merge import join_by_id


def test_join_by_id(tmp_path):
    io.write_table(tmp_path / "a.fits", {"LS_ID_DR11": np.array([5, 1, 9]), "X": np.array([50.0, 10.0, 90.0]),
                                         "V": np.arange(6.0).reshape(3, 2)})
    io.write_table(tmp_path / "b.fits", {"LS_ID_DR11": np.array([7]), "X": np.array([70.0]),
                                         "V": np.ones((1, 2))})
    out = join_by_id(np.array([1, 7, 3, 9]), [tmp_path / "a.fits", tmp_path / "b.fits", tmp_path / "missing.fits"])
    assert out["FOUND"].tolist() == [True, True, False, True]
    assert np.allclose(out["X"][[0, 1, 3]], [10, 70, 90]) and np.isnan(out["X"][2])
    assert out["V"].shape == (4, 2) and np.all(np.isnan(out["V"][2]))


def test_join_rejects_duplicates(tmp_path):
    io.write_table(tmp_path / "a.fits", {"LS_ID_DR11": np.array([1, 1]), "X": np.zeros(2)})
    with pytest.raises(ValueError, match="duplicated"):
        join_by_id(np.array([1]), [tmp_path / "a.fits"])


def test_my_part_covers_everything():
    items = list(range(23))
    parts = [cli.my_part(items, i, 5) for i in range(5)]
    assert sorted(sum(parts, [])) == items
