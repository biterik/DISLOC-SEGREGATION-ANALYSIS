"""Tidy long results table (docs/DEFINITIONS.md sections 9, 10): one row per (snapshot, quantity, region, basis), with
the run, definition set and tool version on every row, and a JSON side file with the full provenance."""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path

COLUMNS = ("run_id", "snapshot", "timestep", "quantity", "region", "basis", "value", "definition_set", "tool_version")


def tool_version(versions: dict) -> str:
    v = versions.get("disegreg", "?")
    h = versions.get("git_hash")
    if h:
        v += f"+g{h[:10]}" + (".dirty" if versions.get("git_dirty") else "")
    return v


def write_table(path, rows, provenance: dict) -> Path:
    """Write rows (dicts with COLUMNS) as CSV plus <path>.provenance.json; both atomically."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".tmp{os.getpid()}")
    with open(tmp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS, extrasaction="raise")
        w.writeheader()
        for r in rows:
            w.writerow({k: (repr(v) if isinstance(v, float) else v) for k, v in r.items()})
    os.replace(tmp, path)
    side = path.with_name(path.name + ".provenance.json")
    tmp = side.with_name(side.name + f".tmp{os.getpid()}")
    tmp.write_text(json.dumps(provenance, indent=1, sort_keys=True, default=str))
    os.replace(tmp, side)
    return path


def read_table(path) -> list[dict]:
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        r["value"] = float(r["value"])
        r["timestep"] = int(r["timestep"])
    return rows
