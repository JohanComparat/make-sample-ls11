"""Selection configuration: a nested dict read from YAML, plus a short hash recorded in outputs.

A configuration may start with ``base: <file>`` (relative to its own directory): the base is read
first and this file's keys are deep-merged on top (mappings merged, everything else replaced, so
``null`` switches a cut off).
"""

from __future__ import annotations

import copy
import hashlib
import os
from pathlib import Path

import yaml

DEFAULT = Path(__file__).resolve().parents[1] / "config" / "default.yaml"


def deep_merge(base: dict, over: dict) -> dict:
    out = copy.deepcopy(base)
    for k, v in over.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = deep_merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


def _read(path: Path) -> dict:
    with open(path) as f:
        cfg = yaml.safe_load(f) or {}
    base = cfg.pop("base", None)
    if base:
        cfg = deep_merge(_read((path.parent / base).resolve()), cfg)
    return cfg


def load_config(path: str | Path | None = None) -> dict:
    path = Path(path or os.environ.get("LS11_CONFIG", DEFAULT)).resolve()
    cfg = _read(path)
    cfg.setdefault("tag", "bgsl")
    cfg["_path"] = str(path)
    return cfg


def config_hash(cfg: dict, section: str | None = None) -> str:
    """Short hash of the configuration (or one section of it), written as CFGHASH in headers."""
    c = {k: v for k, v in cfg.items() if not k.startswith("_")}
    if section is not None:
        c = c[section]
    text = yaml.safe_dump(c, sort_keys=True)
    return hashlib.sha1(text.encode()).hexdigest()[:12]
