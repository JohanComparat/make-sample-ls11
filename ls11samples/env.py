"""Paths, all from environment variables (laptop defaults; set them on CC-IN2P3).

LS11_DIR        Legacy Surveys DR11 root                 [~/data/legacysurvey/dr11]
LS11_REGION     south | north                             [south]
LS11_SWEEP_VER  sweep version directory                  [11.0]
LS11_SWEEPS     optional glob restricting the sweeps      [sweep-*.fits]
LS11_SWEEP_OUT  root of the per-sweep products           [$LS11_DIR/$LS11_REGION/sweep]
                (<root>/<ver>-<name>/<sweep>-<name>.fits, next to <ver>/ and <ver>-photo-z/)
LS11_RANDOMS    glob of random files (relative to $LS11_DIR/$LS11_REGION/randoms, or absolute)
                                                          [randoms-$LS11_REGION-1-0.fits]
LS11_OUT        root of the per-run products (randoms, maps, samples, logs): <root>/<tag>/
                                                          [$LS11_DIR/$LS11_REGION]
LS11_CONFIG     selection configuration                   [config/default.yaml of this repo]
LS11_GAIA_MAPS  directory of full-sky Gaia star-density maps (<nside:04d>/GAIA_*.fits)
                                                          [~/data/legacysurvey/dr10/systematics]
LS11_NPROC      worker processes                          [min(8, cpu count)]
"""

from __future__ import annotations

import glob
import os
from dataclasses import dataclass
from pathlib import Path


def _env_path(name: str, default: str) -> Path:
    return Path(os.path.expandvars(os.path.expanduser(os.environ.get(name, default))))


def nproc() -> int:
    return int(os.environ.get("LS11_NPROC", min(8, os.cpu_count() or 1)))


@dataclass(frozen=True)
class Paths:
    ls_dir: Path
    region: str
    sweep_ver: str
    sweep_glob: str
    randoms_glob: str
    out: Path
    sweep_out: Path
    gaia_maps: Path

    @property
    def region_dir(self) -> Path:
        return self.ls_dir / self.region

    @property
    def sweep_dir(self) -> Path:
        return self.region_dir / "sweep" / self.sweep_ver

    @property
    def pz_dir(self) -> Path:
        return self.region_dir / "sweep" / f"{self.sweep_ver}-photo-z"

    def sweeps(self) -> list[Path]:
        return sorted(self.sweep_dir.glob(self.sweep_glob))

    def random_files(self) -> list[Path]:
        pattern = self.randoms_glob
        if not os.path.isabs(pattern):
            pattern = str(self.region_dir / "randoms" / pattern)
        return sorted(Path(p) for p in glob.glob(pattern))

    # per-sweep products, next to the sweeps: <sweep_out>/<ver>-<name>/<sweep stem>-<name>.fits
    def product_dir(self, name: str) -> Path:
        return self.sweep_out / f"{self.sweep_ver}-{name}"

    def product(self, sweep: str | Path, name: str) -> Path:
        return self.product_dir(name) / f"{Path(sweep).stem}-{name}.fits"

    # per-run products: <out>/<tag>/...
    def run_dir(self, tag: str, *parts: str) -> Path:
        d = self.out.joinpath(tag, *parts)
        d.mkdir(parents=True, exist_ok=True)
        return d

    def rand_file(self, tag: str) -> Path:
        return self.out / tag / f"LS11_{tag}_RAND.fits"


def get_paths() -> Paths:
    ls_dir = _env_path("LS11_DIR", "~/data/legacysurvey/dr11")
    region = os.environ.get("LS11_REGION", "south")
    return Paths(
        ls_dir=ls_dir,
        region=region,
        sweep_ver=os.environ.get("LS11_SWEEP_VER", "11.0"),
        sweep_glob=os.environ.get("LS11_SWEEPS", "sweep-*.fits"),
        randoms_glob=os.path.expandvars(os.path.expanduser(
            os.environ.get("LS11_RANDOMS", f"randoms-{region}-1-0.fits"))),
        out=_env_path("LS11_OUT", str(ls_dir / region)),
        sweep_out=_env_path("LS11_SWEEP_OUT", str(ls_dir / region / "sweep")),
        gaia_maps=_env_path("LS11_GAIA_MAPS", "~/data/legacysurvey/dr10/systematics"),
    )
