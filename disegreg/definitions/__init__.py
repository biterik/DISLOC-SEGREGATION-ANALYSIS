"""Loading of the versioned definition sets (``D1.yaml``, ``D2.yaml``, ...).

A definition set is the complete list of numerical parameters of docs/DEFINITIONS.md. Every output of disegreg records
the label and the SHA-256 of the YAML file it used, so a number can always be traced to its definitions.
"""

from __future__ import annotations

import copy
import hashlib
from dataclasses import dataclass, field
from pathlib import Path

import yaml

DEFAULT_LABEL = "D2"   # v0.2: D2 = D1 rev. 4 unchanged + sections 12-16 (Erik 25 Nov 2026)
_DIR = Path(__file__).resolve().parent


@dataclass(frozen=True)
class Definitions:
    label: str
    status: str
    sha256: str
    path: str
    params: dict = field(repr=False)

    def __getitem__(self, key):
        # deep copy: callers cannot change the loaded set by accident
        return copy.deepcopy(self.params[key])


def available() -> list[str]:
    """Labels of all definition sets shipped with this version of disegreg."""
    return sorted(p.stem for p in _DIR.glob("D*.yaml"))


def load_definitions(label_or_path: str | Path = DEFAULT_LABEL) -> Definitions:
    """Load a shipped definition set by label ("D1") or any YAML file by path."""
    p = Path(label_or_path)
    if not p.suffix:
        p = _DIR / f"{label_or_path}.yaml"
        if not p.is_file():
            raise FileNotFoundError(f"definition set {label_or_path!r} not shipped; available: {available()}")
    raw = p.read_bytes()
    params = yaml.safe_load(raw)
    for key in ("label", "status"):
        if key not in params:
            raise ValueError(f"{p}: missing top-level key {key!r}")
    return Definitions(
        label=params["label"],
        status=params["status"],
        sha256=hashlib.sha256(raw).hexdigest(),
        path=str(p),
        params=params,
    )


PROPOSALS_LABEL = "proposals_v0.2"


def load_proposals(label_or_path: str | Path = PROPOSALS_LABEL) -> Definitions:
    """Parameters of the still open v0.2 proposals (``proposals_v0.2.yaml``; since 25 Nov only P5, ordered domains).
    NOT a definition set: every number computed with them carries status = proposal (docs/PLAN_v0.2.md, rule 3)."""
    p = Path(label_or_path)
    if p.suffix not in (".yaml", ".yml"):
        p = _DIR / f"{label_or_path}.yaml"
    return load_definitions(p)


def label_number(label: str) -> int | None:
    """2 for "D2"; None for a label that is not of the form D<n> (a YAML file given by path with its own label)."""
    s = str(label or "")
    return int(s[1:]) if s[:1] == "D" and s[1:].isdigit() else None
