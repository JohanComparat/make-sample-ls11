"""Selection configuration: a nested dict read from YAML, plus a short hash recorded in outputs."""

from __future__ import annotations

import hashlib
import os
from pathlib import Path

import yaml

DEFAULT = Path(__file__).resolve().parents[1] / "config" / "default.yaml"


def load_config(path: str | Path | None = None) -> dict:
    path = Path(path or os.environ.get("LS11_CONFIG", DEFAULT))
    with open(path) as f:
        cfg = yaml.safe_load(f)
    cfg["_path"] = str(path)
    return cfg


def config_hash(cfg: dict, section: str | None = None) -> str:
    """Short hash of the configuration (or one section of it), written as CFGHASH in headers."""
    c = {k: v for k, v in cfg.items() if not k.startswith("_")}
    if section is not None:
        c = c[section]
    text = yaml.safe_dump(c, sort_keys=True)
    return hashlib.sha1(text.encode()).hexdigest()[:12]
