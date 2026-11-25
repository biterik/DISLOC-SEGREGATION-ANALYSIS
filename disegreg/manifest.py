"""Stage 0: the run manifest (docs/DESIGN.md, Stages). One YAML file per run, written by the PROJECT instance in its own
folder; disegreg only reads and validates it.

Every field is listed in FIELDS with its meaning; ``disegreg manifest template`` prints a commented template and
``disegreg manifest check <file>`` validates a manifest. Paths may be absolute, start with ~, or be relative to the
manifest file's folder.
"""

from __future__ import annotations

import glob as _glob
import itertools
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import yaml

SCHEMA_VERSION = 1
K_B = 8.617333262e-5  # Boltzmann constant, eV/K (CODATA 2018)
DELTA_TOLERANCE = 0.02  # |Delta - (dmu - dmu_coex)/kT| allowed (rounding of dmu in file names / inputs)

# name: (required, type(s), meaning)
FIELDS = {
    "schema_version": (True, int, "manifest schema version, currently 1"),
    "project": (True, str, "project name (folder under ~/Desktop/SIMULATIONS)"),
    "thread": (False, str, "thread folder of the project"),
    "run_id": (True, str, "unique run name (normally the run folder name)"),
    "potential": (True, str, "interatomic potential label, e.g. LeeShim2004-MEAM"),
    "host_element": (True, str, "element of atom type 1"),
    "solute_element": (True, str, "element of atom type 2"),
    "character": (True, str, "dislocation character: d0 (screw), d30, d60, d90 (edge)"),
    "frame": (True, dict, "crystal directions (integer Miller indices) of simulation X (line), Y, Z (glide-plane normal)"),
    "frame_source": (True, str, "cell-builder file the frame was read from"),
    "a_cell": (True, (int, float), "lattice parameter the cell was built / scaled with (A); used for all shell cut-offs"),
    "a_cell_source": (True, str, "where a_cell comes from"),
    "temperature": (True, (int, float), "temperature T of the MD/MC run (K)"),
    "dmu": (True, (int, float), "chemical-potential difference mu_solute - mu_host (eV), in THIS convention"),
    "dmu_lammps": (True, (int, float), "the value as written in the LAMMPS input (eV)"),
    "dmu_lammps_sign": (True, int, "+1 or -1: dmu = dmu_lammps_sign * dmu_lammps"),
    "dmu_coex": (True, (int, float), "dmu at coexistence at this T (eV, same convention)"),
    "coex_line": (True, str, "which coexistence line: incoherent (bulk phase) or coherent (clustering line)"),
    "dmu_coex_source": (True, str, "where dmu_coex comes from"),
    "delta": (True, (int, float), "offset Delta = (dmu - dmu_coex) / kT; negative = below coexistence"),
    "dumps": (True, str, "glob of the QUENCHED snapshot dumps of this run"),
    "quenched": (True, bool, "must be true: only energy-minimised snapshots are analysed"),
    "md_timestep_ps": (False, (int, float), "MD timestep (ps), to convert dump step numbers to time"),
    "z_boundary": (True, str, "Z faces: fixed (boundary layers held fixed) or free (free surfaces)"),
    "z_gp": (True, (int, float), "glide-plane height in dump coordinates (A), from the cell-builder metadata"),
    "z_gp_source": (True, str, "how z_gp was obtained (builder file, scaling applied)"),
    "bulk_control": (False, (str, type(None)), "run id of the bulk control cell (no dislocation, same size, dmu, length)"),
    "bulk_control_reason": (False, str, "required if bulk_control is empty: why there is none"),
    "x_ref_source": (True, str, "bulk_control | isotherm (DEFINITIONS 5)"),
    "x_iso": (False, (int, float), "solute fraction of the potential's bulk SGC isotherm at (dmu, T)"),
    "x_iso_source": (False, str, "where x_iso comes from (required if x_iso is given)"),
    "start_type": (True, str, "bottom-up | top-down | seeded"),
    "seeded_from": (False, str, "run id the start configuration was taken from (start_type seeded)"),
    "mc_stream": (False, str, "MC stream file (total solute count vs step), for the drift test"),
    "mc_stream_column": (False, str, "solute-count column of the stream (default: n*_interior or n*_mcgroup)"),
    "notes": (False, str, "free text"),
    # v0.2 (optional, schema 1 unchanged: every v0.1 manifest stays valid)
    "cell_type": (False, str, "dislocation (default) | bulk_control: a cell without dislocation that gives x_bc to the runs "
                  "naming it in bulk_control (DEFINITIONS 5)"),
    "n_atoms": (False, int, "number of atoms in the cell; required for a bulk control (D2 section 15: its size is recorded)"),
}
# not required in a bulk-control manifest (no line, no glide plane, no reference of its own)
NOT_REQUIRED_FOR_BULK = ("character", "frame", "frame_source", "z_gp", "z_gp_source", "x_ref_source", "start_type")
CHOICES = {
    "character": ("d0", "d30", "d60", "d90"),
    "coex_line": ("incoherent", "coherent"),
    "z_boundary": ("fixed", "free"),
    "x_ref_source": ("bulk_control", "isotherm"),
    "start_type": ("bottom-up", "top-down", "seeded"),
    "cell_type": ("dislocation", "bulk_control"),
}


class ManifestError(ValueError):
    pass


def _check_frame(frame) -> list[str]:
    errs = []
    try:
        X, Y, Z = (np.asarray(frame[k], dtype=float) for k in ("X", "Y", "Z"))
    except Exception:  # noqa: BLE001 -- any malformed frame is reported the same way
        return ["frame: needs X, Y, Z, each three integers"]
    for k, v in zip("XYZ", (X, Y, Z)):
        if v.shape != (3,) or not np.allclose(v, np.round(v)) or not v.any():
            errs.append(f"frame.{k}: must be three integers, not all zero")
    if errs:
        return errs
    for (p, u), (q, w) in itertools.combinations((("X", X), ("Y", Y), ("Z", Z)), 2):
        if abs(u @ w) > 1e-9:
            errs.append(f"frame: {p} and {q} are not orthogonal")
    if not errs and np.cross(X, Y) @ Z <= 0:
        errs.append("frame: X, Y, Z is not right-handed")
    if not np.allclose(np.sort(np.abs(Z)) / np.abs(Z).max(), 1.0):
        errs.append("frame.Z: must be a <111> direction (glide-plane normal)")
    return errs


def validate(data: dict) -> list[str]:
    """List of problems (empty if the manifest is valid)."""
    errs = []
    if not isinstance(data, dict):
        return ["manifest is not a YAML mapping"]
    for k in data:
        if k not in FIELDS:
            errs.append(f"unknown field '{k}'")
    bulk = data.get("cell_type") == "bulk_control"
    if "cell_type" in data and data["cell_type"] not in CHOICES["cell_type"]:
        errs.append(f"cell_type: {data['cell_type']!r} not one of {CHOICES['cell_type']}")
    for k, (req, typ, _) in FIELDS.items():
        if k not in data:
            if req and not (bulk and k in NOT_REQUIRED_FOR_BULK):
                errs.append(f"missing field '{k}'")
            continue
        v = data[k]
        if isinstance(v, bool) and typ in (int, (int, float)):
            errs.append(f"{k}: expected a number, got {v!r}")
        elif not isinstance(v, typ):
            errs.append(f"{k}: wrong type {type(v).__name__}")
    if errs:
        return errs
    if data["schema_version"] != SCHEMA_VERSION:
        errs.append(f"schema_version {data['schema_version']} != {SCHEMA_VERSION}")
    for k, allowed in CHOICES.items():
        if k in data and data[k] not in allowed and k != "cell_type":
            errs.append(f"{k}: {data[k]!r} not one of {allowed}")
    if "frame" in data:
        errs += _check_frame(data["frame"])
    if data["quenched"] is not True:
        errs.append("quenched: must be true (DEFINITIONS 1: quenched snapshots only)")
    if data["dmu_lammps_sign"] not in (1, -1):
        errs.append("dmu_lammps_sign: must be +1 or -1")
    elif abs(data["dmu_lammps_sign"] * data["dmu_lammps"] - data["dmu"]) > 1e-9:
        errs.append("dmu != dmu_lammps_sign * dmu_lammps")
    if data["temperature"] <= 0 or data["a_cell"] <= 0:
        errs.append("temperature and a_cell must be positive")
    else:
        d = (data["dmu"] - data["dmu_coex"]) / (K_B * data["temperature"])
        if abs(d - data["delta"]) > DELTA_TOLERANCE:
            errs.append(f"delta {data['delta']} != (dmu - dmu_coex)/kT = {d:.4f}")
    if "n_atoms" in data and data["n_atoms"] <= 0:
        errs.append("n_atoms must be positive")
    if bulk:
        if data.get("bulk_control"):
            errs.append("a bulk-control manifest cannot name a bulk_control of its own")
        if "n_atoms" not in data:
            errs.append("a bulk-control manifest needs n_atoms (its cell size, D2 section 15)")
        return errs + (["x_iso given without x_iso_source"] if "x_iso" in data and not data.get("x_iso_source") else [])
    if not data.get("bulk_control") and not data.get("bulk_control_reason"):
        errs.append("bulk_control is empty: bulk_control_reason is required")
    if data["x_ref_source"] == "bulk_control" and not data.get("bulk_control"):
        errs.append("x_ref_source bulk_control needs a bulk_control run id")
    if data["x_ref_source"] == "isotherm" and "x_iso" not in data:
        errs.append("x_ref_source isotherm needs x_iso")
    if "x_iso" in data and not data.get("x_iso_source"):
        errs.append("x_iso given without x_iso_source")
    if data["start_type"] == "seeded" and not data.get("seeded_from"):
        errs.append("start_type seeded needs seeded_from")
    return errs


def _step_key(path: str):
    """Sort dumps by the last number in the file name (the step), then by name."""
    nums = re.findall(r"\d+", Path(path).name)
    return (int(nums[-1]) if nums else -1, path)


@dataclass(frozen=True)
class Manifest:
    path: Path
    data: dict

    def __getitem__(self, key):
        return self.data[key]

    def get(self, key, default=None):
        return self.data.get(key, default)

    @property
    def is_bulk_control(self) -> bool:
        return self.data.get("cell_type", "dislocation") == "bulk_control"

    @property
    def kT(self) -> float:
        return K_B * float(self.data["temperature"])

    def resolve(self, p: str) -> Path:
        q = Path(p).expanduser()
        return q if q.is_absolute() else (self.path.parent / q)

    def dump_files(self) -> list[Path]:
        """Dump files matching the glob, sorted by step number."""
        pattern = str(self.resolve(self.data["dumps"]))
        return [Path(p) for p in sorted(_glob.glob(pattern), key=_step_key)]

    def frame_matrix(self) -> np.ndarray:
        """Rows = unit crystal vectors of simulation X, Y, Z (rotation crystal -> simulation frame)."""
        m = np.array([self.data["frame"][k] for k in ("X", "Y", "Z")], dtype=float)
        return m / np.linalg.norm(m, axis=1)[:, None]


def load_manifest(path, strict: bool = True) -> Manifest:
    path = Path(path).expanduser().resolve()
    with open(path) as f:
        data = yaml.safe_load(f)
    errs = validate(data)
    if errs and strict:
        raise ManifestError(f"{path}:\n  " + "\n  ".join(errs))
    return Manifest(path=path, data=data)


def template() -> str:
    """Commented YAML template with every field."""
    lines = [f"# disegreg run manifest (schema {SCHEMA_VERSION}); required fields first. Paths: absolute, ~, or",
             "# relative to this file's folder. Check with:  disegreg manifest check <this file>"]
    for req in (True, False):
        lines.append("# ---- required" if req else "# ---- optional")
        for k, (r, typ, meaning) in FIELDS.items():
            if r != req:
                continue
            if k == "frame":
                lines += [f"frame:                 # {meaning}", "  X: [0, -1, 1]", "  Y: [2, -1, -1]", "  Z: [1, 1, 1]"]
                continue
            val = {"schema_version": str(SCHEMA_VERSION), "quenched": "true", "dmu_lammps_sign": "1"}.get(k, "")
            if k in CHOICES:
                val = CHOICES[k][0]
                meaning += f"  [{' | '.join(CHOICES[k])}]"
            lines.append(f"{k + ':':22} {val:10} # {meaning}")
    return "\n".join(lines) + "\n"
