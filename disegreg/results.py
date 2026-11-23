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


# --------------------------------------------------------------------------------------------- across runs (Stage 4 input)
RESULTS_COLUMNS = ("project", "run_id", "potential", "host", "solute", "character", "temperature", "delta", "start_type",
                   "quantity", "region", "basis", "value", "error", "n_samples", "flat", "run_converged",
                   "window_fraction", "x_ref", "x_ref_source", "definition_set", "tool_version")
RUN_FIELDS = {"project": "project", "potential": "potential", "host": "host_element", "solute": "solute_element",
              "character": "character", "temperature": "temperature", "delta": "delta", "start_type": "start_type"}


def collect_runs(stage3_files) -> tuple[list[dict], dict]:
    """Merge Stage-3 tables of many runs into one tidy long table (DEFINITIONS 9): value = window mean, error =
    block-averaging SEM. Run metadata come from the manifest copy in each table's provenance file. Refuses to mix
    definition sets."""
    rows, sources, labels = [], [], set()
    for f in stage3_files:
        f = Path(f)
        prov = json.loads(f.with_name(f.name + ".provenance.json").read_text())
        man = prov.get("manifest_data", {})
        s3 = prov.get("stage3", {})
        labels.add(prov.get("definition_set"))
        meta = {k: man.get(v, "") for k, v in RUN_FIELDS.items()}
        for r in read_table(f):
            rows.append({**meta, "run_id": r["run_id"], "quantity": r["quantity"], "region": r["region"],
                         "basis": r["basis"], "value": r["mean"], "error": r["sem"], "n_samples": r["n_samples"],
                         "flat": r["flat"], "run_converged": s3.get("converged", ""),
                         "window_fraction": s3.get("window_fraction", ""), "x_ref": prov.get("x_ref", ""),
                         "x_ref_source": prov.get("x_ref_source", ""), "definition_set": r["definition_set"],
                         "tool_version": r["tool_version"]})
        sources.append({"file": str(f), "run_id": man.get("run_id"), "tool_version": prov.get("tool_version"),
                        "manifest": prov.get("manifest")})
    if len(labels) > 1:
        raise ValueError(f"tables from different definition sets {sorted(map(str, labels))} must not be merged")
    return rows, {"definition_set": labels.pop() if labels else None, "inputs": sources}
