"""Reading sweeps, their photo-z files and the randoms; writing FITS tables.

Tables are plain ``dict[str, np.ndarray]`` in native byte order (fast to slice and to concatenate);
:func:`write_table` turns them into a FITS binary table with a header.
"""

from __future__ import annotations

import re
from collections.abc import Iterator, Mapping, Sequence
from pathlib import Path

import fitsio
import numpy as np

Table = dict

_SWEEP = re.compile(r"(\d{3})([mp])(\d{3})-(\d{3})([mp])(\d{3})")


# --------------------------------------------------------------------------- identifiers / names
def ls_id(release, brickid, objid) -> np.ndarray:
    """LS_ID_DR11 = RELEASE << 42 | BRICKID << 22 | OBJID (as in the DR11 photo-z files)."""
    return ((np.asarray(release, np.int64) << 42) | (np.asarray(brickid, np.int64) << 22)
            | np.asarray(objid, np.int64))


def pz_path(sweep: str | Path) -> Path:
    """Row-matched photo-z file of a sweep: ``<ver>/X.fits`` -> ``<ver>-photo-z/X-pz.fits``."""
    sweep = Path(sweep)
    return sweep.parent.parent / f"{sweep.parent.name}-photo-z" / sweep.name.replace(".fits", "-pz.fits")


def sweep_box(name: str | Path) -> tuple[float, float, float, float]:
    """(ra_min, ra_max, dec_min, dec_max) of a sweep file name, e.g. sweep-000m005-005p000."""
    m = _SWEEP.search(Path(name).name)
    if m is None:
        raise ValueError(f"not a sweep file name: {name}")
    ra0, s0, d0, ra1, s1, d1 = m.groups()
    sign = {"m": -1.0, "p": 1.0}
    return float(ra0), float(ra1), sign[s0] * float(d0), sign[s1] * float(d1)


def in_boxes(ra: np.ndarray, dec: np.ndarray, boxes: Sequence[tuple]) -> np.ndarray:
    """True for points inside any box; boxes are half-open [min, max) like the sweep files."""
    ok = np.zeros(np.shape(ra), bool)
    for ra0, ra1, d0, d1 in boxes:
        ok |= (ra >= ra0) & (ra < ra1) & (dec >= d0) & (dec < d1)
    return ok


def box_area(boxes: Sequence[tuple]) -> float:
    """Sky area (deg^2) of non-overlapping RA/Dec boxes."""
    return float(sum((ra1 - ra0) * (180 / np.pi) * (np.sin(np.radians(d1)) - np.sin(np.radians(d0)))
                     for ra0, ra1, d0, d1 in boxes))


# --------------------------------------------------------------------------- reading
def _native(a: np.ndarray) -> np.ndarray:
    if a.dtype.kind in "SUO":
        return a
    return a.astype(a.dtype.newbyteorder("="), copy=False)


def read_table(path: str | Path, columns: Sequence[str] | None = None, rows=None,
               ext: int | str = 1) -> Table:
    """Columns of a FITS binary table as a dict of native-endian arrays."""
    with fitsio.FITS(str(path)) as f:
        if columns is None:
            columns = f[ext].get_colnames()
        data = f[ext].read(columns=list(columns), rows=rows)
    out = {}
    for c in columns:
        a = data[c]
        if a.dtype.kind == "S":
            a = np.char.strip(a.astype("U"))
        out[c] = _native(np.ascontiguousarray(a))
    return out


def read_header(path: str | Path, ext: int | str = 0) -> fitsio.FITSHDR:
    return fitsio.read_header(str(path), ext)


def nrows(path: str | Path, ext: int | str = 1) -> int:
    with fitsio.FITS(str(path)) as f:
        return f[ext].get_nrows()


def read_pz(sweep: str | Path, rows: np.ndarray, release, brickid, objid,
            columns: Sequence[str] | None = None) -> Table:
    """Photo-z rows of ``sweep`` (row-matched file), checked object by object against the sweep IDs."""
    path = pz_path(sweep)
    if not path.exists():
        raise FileNotFoundError(f"no photo-z file for {sweep} (expected {path})")
    if nrows(path) != nrows(sweep):
        raise ValueError(f"{path} and {sweep} have different numbers of rows")
    cols = ["LS_ID_DR11", "RELEASE", "BRICKID", "OBJID"] + [
        c for c in (columns or []) if c not in ("LS_ID_DR11", "RELEASE", "BRICKID", "OBJID")]
    pz = read_table(path, cols if columns is not None else None, rows=rows)
    ids = ls_id(release, brickid, objid)
    same = (ls_id(pz["RELEASE"], pz["BRICKID"], pz["OBJID"]) == ids) & (pz["LS_ID_DR11"] == ids)
    if not np.all(same):
        raise ValueError(f"{path} is not row-matched to {sweep}: {np.count_nonzero(~same)} rows differ")
    return pz


def iter_rows(path: str | Path, columns: Sequence[str], chunk: int, ext: int | str = 1
              ) -> Iterator[tuple[int, Table]]:
    """(first row, table) chunks of a large FITS table, reading only ``columns``."""
    n = nrows(path, ext)
    with fitsio.FITS(str(path)) as f:
        hdu = f[ext]
        for start in range(0, n, chunk):
            stop = min(n, start + chunk)
            data = hdu.read(columns=list(columns), rows=np.arange(start, stop))
            t = {}
            for c in columns:
                a = data[c]
                if a.dtype.kind == "S":
                    a = np.char.strip(a.astype("U"))
                t[c] = _native(np.ascontiguousarray(a))
            yield start, t


# --------------------------------------------------------------------------- writing
def to_records(table: Mapping[str, np.ndarray]) -> np.ndarray:
    names = list(table)
    n = len(table[names[0]]) if names else 0
    dtype = []
    for k in names:
        a = np.asarray(table[k])
        if a.dtype.kind == "U":
            a = a.astype("S")
        dtype.append((k, a.dtype, a.shape[1:]) if a.ndim > 1 else (k, a.dtype))
    rec = np.empty(n, dtype=dtype)
    for k in names:
        a = np.asarray(table[k])
        rec[k] = a.astype("S") if a.dtype.kind == "U" else a
    return rec


def write_table(path: str | Path, table: Mapping[str, np.ndarray], header: Mapping | None = None,
                extname: str = "DATA", extra: Sequence[tuple[str, Mapping, Mapping | None]] = ()) -> Path:
    """Write ``table`` (and optional ``extra`` (extname, table, header) HDUs) to ``path``.

    Written to a temporary name first and renamed, so a crashed job never leaves a truncated file.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".part")
    if tmp.exists():
        tmp.unlink()
    hdr = [{"name": k, "value": v} for k, v in (header or {}).items()]
    with fitsio.FITS(str(tmp), "rw", clobber=True) as f:
        f.write(to_records(table), extname=extname, header=hdr)
        for name, t, h in extra:
            f.write(to_records(t), extname=name,
                    header=[{"name": k, "value": v} for k, v in (h or {}).items()])
    tmp.replace(path)
    return path


def concat(tables: Sequence[Mapping[str, np.ndarray]]) -> Table:
    tables = [t for t in tables if t and len(next(iter(t.values()))) > 0]
    if not tables:
        return {}
    return {k: np.concatenate([t[k] for t in tables]) for k in tables[0]}


def take(table: Mapping[str, np.ndarray], idx) -> Table:
    return {k: v[idx] for k, v in table.items()}
