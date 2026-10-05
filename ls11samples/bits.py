"""DR11 MASKBITS / FITBITS definitions, checked against the sweep headers.

The sweep primary header documents each bit (MBIT_n = name, FBIT_n = name). Cuts are configured by
bit NAME; :func:`check_header` makes sure the numbers below are the ones of the files being read.
The randoms carry the same MASKBITS (their header has no bit list, so the sweep check covers them).
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping

import numpy as np

MASKBITS = {
    0: "NPRIMARY", 1: "BRIGHT", 2: "SATUR_G", 3: "SATUR_R", 4: "SATUR_Z", 5: "ALLMASK_G",
    6: "ALLMASK_R", 7: "ALLMASK_Z", 8: "WISEM1", 9: "WISEM2", 10: "BAILOUT", 11: "MEDIUM",
    12: "GALAXY", 13: "CLUSTER", 14: "SATUR_I", 15: "ALLMASK_I", 16: "SUB_BLOB", 17: "RESOLVED",
    18: "MCLOUDS", 19: "WISE_GAIA",
}
FITBITS = {
    0: "FORCED_POINTSOURCE", 1: "FIT_BACKGROUND", 2: "HIT_RADIUS_LIMIT", 3: "HIT_SERSIC_LIMIT",
    4: "FROZEN", 5: "BRIGHT", 6: "MEDIUM", 7: "GAIA", 8: "TYCHO2", 9: "LARGEGALAXY", 10: "WALKER",
    11: "RUNNER", 12: "GAIA_POINTSOURCE", 13: "ITERATIVE",
}
MASKBIT = {v: k for k, v in MASKBITS.items()}
FITBIT = {v: k for k, v in FITBITS.items()}


def header_bits(header: Mapping, prefix: str) -> dict[int, str]:
    """{bit: name} from header keys ``<prefix>_<n>`` (e.g. MBIT_12 = 'GALAXY')."""
    out = {}
    for key in header.keys():
        if key.startswith(prefix + "_") and key[len(prefix) + 1:].isdigit():
            out[int(key[len(prefix) + 1:])] = str(header[key]).strip()
    return out


def check_header(header: Mapping) -> None:
    """Raise if the sweep header bit definitions differ from :data:`MASKBITS` / :data:`FITBITS`."""
    for prefix, ref in (("MBIT", MASKBITS), ("FBIT", FITBITS)):
        got = header_bits(header, prefix)
        if not got:
            raise ValueError(f"no {prefix}_n keywords in the sweep header")
        if got != ref:
            diff = {b: (got.get(b), ref.get(b)) for b in set(got) | set(ref) if got.get(b) != ref.get(b)}
            raise ValueError(f"{prefix} definitions differ from DR11 (file, expected): {diff}")


def mask_value(names: Iterable[str], table: Mapping[str, int]) -> np.int64:
    """OR of the bits named in ``names``; unknown names raise."""
    value = 0
    for n in names:
        if n not in table:
            raise KeyError(f"unknown bit name {n!r}; known: {sorted(table)}")
        value |= 1 << table[n]
    return np.int64(value)


def any_set(values: np.ndarray, mask: np.int64) -> np.ndarray:
    return (np.asarray(values).astype(np.int64) & mask) != 0
