"""Tidy long results table (docs/DEFINITIONS.md sections 9, 10): one row per (snapshot, quantity, region, basis, ref),
with the run, definition set and tool version on every row, and a JSON side file with the full provenance.

    basis  distance the quantity uses: r (the distance actually used, DEFINITIONS 3.4) or r_dxa; empty if none
    ref    reference composition of E, Gamma, alpha_xref, c_B_over_x_ref, x_ref: c_B (in-cell far field of the same
           snapshot), x_iso (isotherm, manifest), x_bc (bulk control, given on the command line) or given; empty if none
    role   (Stage 3 and across runs) for Gamma under r: headline = Gamma(20 A), check = Gamma(40 A), curve = other R
"""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path

COLUMNS = ("run_id", "snapshot", "timestep", "quantity", "region", "basis", "ref", "value", "definition_set",
           "tool_version")
STAGE3_COLUMNS = ("run_id", "quantity", "region", "basis", "ref", "mean", "sem", "n_samples", "sem_plateau", "change",
                  "change_se", "flat", "flat_rev3_sem", "t_first", "t_last", "definition_set", "tool_version")
GAMMA_CURVE_COLUMNS = ("run_id", "R", "basis", "ref", "role", "mean", "sem", "n_samples", "flat", "definition_set",
                       "tool_version")
PAPER_HOLD = "no paper number before bulk control"

# v0.2 (docs/PLAN_v0.2.md, rule 3): quantities that depend on an undecided [PROPOSAL v0.2] carry status "proposal" in the
# results table, and the table refuses them for role "headline". Modules that compute such quantities register them here.
STATUS_DEFINITION, STATUS_PROPOSAL = "definition", "proposal"
PROPOSAL_QUANTITIES: set[str] = set()


def quantity_status(quantity: str) -> str:
    return STATUS_PROPOSAL if quantity in PROPOSAL_QUANTITIES else STATUS_DEFINITION


class ProposalAsHeadline(ValueError):
    pass


def check_headline_status(rows) -> None:
    """Refuse a number that depends on an undecided proposal as a headline."""
    bad = [r for r in rows if r.get("role") == "headline" and r.get("status") == STATUS_PROPOSAL]
    if bad:
        raise ProposalAsHeadline(f"{len(bad)} row(s) with status proposal and role headline, e.g. "
                                 f"{bad[0].get('quantity')} {bad[0].get('region')} of {bad[0].get('run_id') or bad[0].get('runs')}")


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
            # float(v): numpy scalars are floats too, but their repr is "np.float64(...)" (numpy >= 2)
            w.writerow({k: (repr(float(v)) if isinstance(v, float) else v) for k, v in r.items()})
    os.replace(tmp, path)
    side = path.with_name(path.name + ".provenance.json")
    tmp = side.with_name(side.name + f".tmp{os.getpid()}")
    tmp.write_text(json.dumps(provenance, indent=1, sort_keys=True, default=str))
    os.replace(tmp, side)
    return path


_FLOAT = {"value", "mean", "sem", "change", "change_se", "t_first", "t_last"}
_INT = {"timestep", "n_samples"}
_BOOL = {"sem_plateau", "flat", "flat_rev3_sem", "check_required"}


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
        r.setdefault("ref", "")
    return rows


# ------------------------------------------------------------------------------------------ headline, check, paper hold
def gamma_role(quantity, region, basis, headline_radius: float, check_radius: float) -> str:
    """headline = Gamma(20 A), check = Gamma(40 A), curve = any other R; only for Gamma under r (DEFINITIONS 6)."""
    if quantity != "Gamma" or basis != "r":
        return ""
    if region == f"{float(headline_radius):g}":
        return "headline"
    if region == f"{float(check_radius):g}":
        return "check"
    return "curve"


def paper_hold(manifest_data: dict, rules) -> str:
    """PAPER_HOLD if the run has no bulk control and matches a rule of reference.paper_hold_without_bulk_control
    (D1 rev. 4: ACE Ni-Cu), else ''."""
    if manifest_data.get("bulk_control"):
        return ""
    pot = str(manifest_data.get("potential", "")).lower()
    for rule in rules or []:
        if (str(rule.get("potential_contains", "")).lower() in pot
                and rule.get("host", manifest_data.get("host_element")) == manifest_data.get("host_element")
                and rule.get("solute", manifest_data.get("solute_element")) == manifest_data.get("solute_element")):
            return PAPER_HOLD
    return ""


def gamma_curve(stage3_rows, headline_radius: float, check_radius: float) -> list[dict]:
    """Gamma(R) per run from the Stage-3 rows: one row per (basis, ref, R), R ascending."""
    out = []
    for r in stage3_rows:
        if r["quantity"] != "Gamma":
            continue
        out.append({"run_id": r["run_id"], "R": float(r["region"]), "basis": r["basis"], "ref": r["ref"],
                    "role": gamma_role("Gamma", r["region"], r["basis"], headline_radius, check_radius),
                    "mean": r["mean"], "sem": r["sem"], "n_samples": r["n_samples"], "flat": r["flat"],
                    "definition_set": r["definition_set"], "tool_version": r["tool_version"]})
    return sorted(out, key=lambda r: (r["basis"], r["ref"], r["R"]))


# --------------------------------------------------------------------------------------------- across runs (Stage 4 input)
RESULTS_COLUMNS = ("project", "run_id", "potential", "host", "solute", "character", "temperature", "delta", "start_type",
                   "quantity", "region", "basis", "ref", "role", "value", "error", "n_samples", "flat", "flat_rev3_sem",
                   "run_converged", "window_fraction", "drift_criterion", "check_required", "paper_flag", "bulk_control",
                   "x_ref_source", "definition_set", "tool_version", "status")
RUN_FIELDS = {"project": "project", "potential": "potential", "host": "host_element", "solute": "solute_element",
              "character": "character", "temperature": "temperature", "delta": "delta", "start_type": "start_type"}


def collect_runs(stage3_files) -> tuple[list[dict], dict]:
    """Merge Stage-3 tables of many runs into one tidy long table (DEFINITIONS 9): value = window mean, error =
    block-averaging SEM. Run metadata come from the manifest copy in each table's provenance file. Refuses to mix
    definition sets.

    role: headline / check / curve for Gamma under r; check_required: the fallback axis was active in the averaging
    window, so Gamma(40 A) must be quoted next to Gamma(20 A); paper_flag: PAPER_HOLD for runs without a bulk control
    that the definition set holds back (D1 rev. 4: ACE Ni-Cu)."""
    rows, sources, labels = [], [], set()
    for f in stage3_files:
        f = Path(f)
        prov = json.loads(f.with_name(f.name + ".provenance.json").read_text())
        man = prov.get("manifest_data", {})
        s3 = prov.get("stage3", {})
        labels.add(prov.get("definition_set"))
        meta = {k: man.get(v, "") for k, v in RUN_FIELDS.items()}
        hr, cr = s3.get("headline_radius", 20.0), s3.get("check_radius", 40.0)
        for r in read_table(f):
            rows.append({**meta, "run_id": r["run_id"], "quantity": r["quantity"], "region": r["region"],
                         "basis": r["basis"], "ref": r["ref"],
                         "role": gamma_role(r["quantity"], r["region"], r["basis"], hr, cr),
                         "value": r["mean"], "error": r["sem"], "n_samples": r["n_samples"],
                         "flat": r["flat"], "flat_rev3_sem": r.get("flat_rev3_sem", ""),
                         "run_converged": s3.get("converged", ""), "window_fraction": s3.get("window_fraction", ""),
                         "drift_criterion": s3.get("criterion", ""), "check_required": s3.get("check_required", ""),
                         "paper_flag": s3.get("paper_flag", ""), "bulk_control": man.get("bulk_control") or "",
                         "x_ref_source": man.get("x_ref_source", ""), "definition_set": r["definition_set"],
                         "tool_version": r["tool_version"], "status": quantity_status(r["quantity"])})
        sources.append({"file": str(f), "run_id": man.get("run_id"), "tool_version": prov.get("tool_version"),
                        "manifest": prov.get("manifest")})
    if len(labels) > 1:
        raise ValueError(f"tables from different definition sets {sorted(map(str, labels))} must not be merged")
    check_headline_status(rows)
    return rows, {"definition_set": labels.pop() if labels else None, "inputs": sources}
