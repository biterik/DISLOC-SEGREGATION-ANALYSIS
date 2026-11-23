"""Tidy long results table (docs/DEFINITIONS.md sections 9, 10): one row per (snapshot, quantity, region, basis), with
the run, definition set and tool version on every row, and a JSON side file with the full provenance."""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path

COLUMNS = ("run_id", "snapshot", "timestep", "quantity", "region", "basis", "value", "definition_set", "tool_version")
STAGE3_COLUMNS = ("run_id", "quantity", "region", "basis", "mean", "sem", "n_samples", "sem_plateau", "change", "change_se",
                  "flat", "flat_change_se", "t_first", "t_last", "definition_set", "tool_version")


def tool_version(versions: dict) -> str:
    v = versions.get("disegreg", "?")
    h = versions.get("git_hash")
    if h:
        v += f"+g{h[:10]}" + (".dirty" if versions.get("git_dirty") else "")
    return v


def write_table(path, rows, provenance: dict, columns=COLUMNS) -> Path:
    """Write rows (dicts with ``columns``) as CSV plus <path>.provenance.json; both atomically."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".tmp{os.getpid()}")
    with open(tmp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=columns, extrasaction="raise")
        w.writeheader()
        for r in rows:
            w.writerow({k: (repr(v) if isinstance(v, float) else v) for k, v in r.items()})
    os.replace(tmp, path)
    side = path.with_name(path.name + ".provenance.json")
    tmp = side.with_name(side.name + f".tmp{os.getpid()}")
    tmp.write_text(json.dumps(provenance, indent=1, sort_keys=True, default=str))
    os.replace(tmp, side)
    return path


_FLOAT = {"value", "mean", "sem", "change", "change_se", "t_first", "t_last"}
_INT = {"timestep", "n_samples"}
_BOOL = {"sem_plateau", "flat", "flat_change_se"}


def read_table(path) -> list[dict]:
    """Stage-2 or Stage-3 CSV back into dicts with typed values."""
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        for k in r:
            if k in _FLOAT:
                r[k] = float(r[k])
            elif k in _INT:
                r[k] = int(r[k])
            elif k in _BOOL:
                r[k] = r[k] == "True"
    return rows
