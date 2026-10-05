"""Shared command-line helpers: logging, config/paths options, job-array slicing."""

from __future__ import annotations

import argparse
import logging
import os
from collections.abc import Sequence

from .config import load_config
from .env import get_paths, nproc


def parser(description: str) -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=description)
    p.add_argument("--config", default=None, help="selection YAML (default: $LS11_CONFIG or config/default.yaml)")
    p.add_argument("--nproc", type=int, default=None, help="worker processes (default: $LS11_NPROC)")
    p.add_argument("--overwrite", action="store_true", help="redo existing outputs")
    p.add_argument("-v", "--verbose", action="store_true")
    return p


def add_slicing(p: argparse.ArgumentParser) -> None:
    """--part i --nparts n: process items i, i+n, i+2n, ... (defaults: SLURM array variables)."""
    p.add_argument("--part", type=int, default=int(os.environ.get("SLURM_ARRAY_TASK_ID", 0)))
    p.add_argument("--nparts", type=int, default=int(os.environ.get("SLURM_ARRAY_TASK_COUNT", 1)))


def setup(args):
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s %(name)s %(levelname)s: %(message)s")
    cfg = load_config(args.config)
    paths = get_paths()
    n = args.nproc or nproc()
    logging.getLogger("ls11samples").info("config %s, output %s, %d processes", cfg["_path"], paths.out, n)
    return cfg, paths, n


def my_part(items: Sequence, part: int, nparts: int) -> list:
    return list(items)[part::nparts]
