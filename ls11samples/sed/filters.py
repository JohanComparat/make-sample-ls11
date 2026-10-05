"""DECam grizY / WISE W1-W2 response curves shared by every SED code.

The curves are the kcorrect v5 ones (``kcorrect/data/responses/{decam,wise}_*.dat``: DECam total
system throughput with atmosphere, per photon; WISE relative response per photon), written in the
format each code reads, so all codes see identical filters.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

KCORRECT_NAMES = {"G": "decam_g", "R": "decam_r", "I": "decam_i", "Z": "decam_z",
                  "W1": "wise_w1", "W2": "wise_w2"}


def responses_dir() -> Path:
    import importlib.util

    spec = importlib.util.find_spec("kcorrect")
    if spec is None or spec.origin is None:
        raise ImportError("kcorrect (v5) is needed for the filter curves")
    return Path(spec.origin).parent / "data" / "responses"


def read_response(band: str) -> tuple[np.ndarray, np.ndarray]:
    """(wavelength [Angstrom], response per photon) of a band (G R I Z W1 W2)."""
    from astropy.io import ascii

    path = responses_dir() / f"{KCORRECT_NAMES[band]}.dat"
    dat = ascii.read(str(path), format="fixed_width")      # the kcorrect reader's format
    order = np.argsort(dat["lambda"])
    wl = np.asarray(dat["lambda"], np.float64)[order]
    tr = np.clip(np.asarray(dat["pass"], np.float64)[order], 0, None)
    return wl, tr


def write_two_column(path: Path, band: str, header: str = "") -> Path:
    wl, tr = read_response(band)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(path, np.c_[wl, tr], fmt="%.2f %.6e", header=header or f"{KCORRECT_NAMES[band]} "
               "(kcorrect v5 response, per photon); lambda[Angstrom] transmission")
    return path


def write_lephare_filters(lephare_dir: str | Path, bands=("G", "R", "I", "Z", "W1", "W2")) -> list[Path]:
    """$LEPHAREDIR/filt/ls11/<kcorrect name>.pb for each band (FILTER_LIST ls11/<name>.pb)."""
    out = Path(lephare_dir) / "filt" / "ls11"
    return [write_two_column(out / f"{KCORRECT_NAMES[b]}.pb", b) for b in bands]


def pivot_wavelength(band: str) -> float:
    wl, tr = read_response(band)
    return float(np.sqrt(np.trapezoid(tr * wl, wl) / np.trapezoid(tr / wl, wl)))
