"""Stage 1 cache: one compressed .npz per snapshot with everything OVITO computes (docs/DESIGN.md, Stages).

Stored: positions (float32, sorted by atom id), types, a-CNA and PTM labels, DXA lines (points, Burgers vectors, lengths,
families), the cell, and provenance (SHA-256 of the dump, Stage-1 parameters, versions). Distances and regions are NOT
cached: they are cheap and definitions may still change, so a definition change never forces OVITO to run again.

This module does not import OVITO; only computing a missing entry does (via disegreg.io_ovito).
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .geometry import Cell
from .provenance import sha256_file

CACHE_FORMAT = 1

# structure-type codes, identical in OVITO's a-CNA and PTM output
STRUCTURE_CODES = {"other": 0, "fcc": 1, "hcp": 2, "bcc": 3, "ico": 4}


@dataclass
class Snapshot:
    positions: np.ndarray          # (N, 3) float32, sorted by atom id
    types: np.ndarray              # (N,) int8 (1 host, 2 solute)
    cell_matrix: np.ndarray        # (3, 4) float64: columns a, b, c, origin
    pbc: np.ndarray                # (3,) bool
    acna: np.ndarray               # (N,) int8 a-CNA structure type (0 other, 1 fcc, 2 hcp, 3 bcc, 4 ico)
    ptm: np.ndarray                # (N,) int8 PTM structure type (same codes)
    seg_points: np.ndarray         # (P, 3) float64, all DXA polylines concatenated (continuous, may leave the cell)
    seg_offsets: np.ndarray        # (S + 1,) int64, polyline k = seg_points[seg_offsets[k]:seg_offsets[k+1]]
    seg_true_b: np.ndarray         # (S, 3) Burgers vector in the cluster's lattice frame (units of a)
    seg_spatial_b: np.ndarray      # (S, 3) Burgers vector in the simulation frame (A)
    seg_length: np.ndarray         # (S,)
    seg_family: np.ndarray         # (S,) str, e.g. "1/6<112>" (Shockley), "1/2<110>" (perfect), "other"
    seg_cluster: np.ndarray        # (S,) int32 DXA crystal cluster
    seg_infinite: np.ndarray       # (S,) bool, line through the periodic boundary
    seg_loop: np.ndarray           # (S,) bool
    ids: np.ndarray | None = None  # (N,) int64 atom ids; None if they are exactly 1..N
    meta: dict = field(default_factory=dict)

    @property
    def n_atoms(self) -> int:
        return len(self.types)

    @property
    def n_segments(self) -> int:
        return len(self.seg_offsets) - 1

    def cell(self) -> Cell:
        return Cell(self.cell_matrix, pbc=tuple(bool(x) for x in self.pbc))

    def polylines(self, families=None) -> list[np.ndarray]:
        """DXA polylines, optionally only those whose Burgers-vector family is in ``families``."""
        out = []
        for k in range(self.n_segments):
            if families is None or self.seg_family[k] in families:
                out.append(self.seg_points[self.seg_offsets[k]:self.seg_offsets[k + 1]])
        return out

    # ---------------------------------------------------------------- file I/O
    _ARRAYS = ("positions", "types", "cell_matrix", "pbc", "acna", "ptm", "seg_points", "seg_offsets", "seg_true_b",
               "seg_spatial_b", "seg_length", "seg_family", "seg_cluster", "seg_infinite", "seg_loop")

    def save(self, path) -> Path:
        """Write atomically (temporary file, then rename), so an interrupted job never leaves a broken entry."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        arrays = {k: getattr(self, k) for k in self._ARRAYS}
        if self.ids is not None:
            arrays["ids"] = self.ids
        meta = dict(self.meta, cache_format=CACHE_FORMAT)
        arrays["meta_json"] = np.array(json.dumps(meta, sort_keys=True))
        tmp = path.with_name(path.name + f".tmp{os.getpid()}")
        with open(tmp, "wb") as f:
            np.savez_compressed(f, **arrays)
        os.replace(tmp, path)
        return path

    @classmethod
    def load(cls, path) -> "Snapshot":
        with np.load(path, allow_pickle=False) as z:
            kw = {k: z[k] for k in cls._ARRAYS}
            kw["ids"] = z["ids"] if "ids" in z.files else None
            kw["meta"] = json.loads(str(z["meta_json"]))
        if kw["meta"].get("cache_format") != CACHE_FORMAT:
            raise ValueError(f"{path}: cache format {kw['meta'].get('cache_format')} != {CACHE_FORMAT}")
        return cls(**kw)


# -------------------------------------------------------------------- cache location and lookup
def default_cache_root() -> Path:
    """$DISEGREG_CACHE, else the agreed location: M5 ~/Desktop/SIMULATIONS/DISEGREG-CACHE, cmmg /u/$USER/PTMP/DISEGREG-CACHE."""
    env = os.environ.get("DISEGREG_CACHE")
    if env:
        return Path(env).expanduser()
    if platform.system() == "Darwin":
        return Path.home() / "Desktop" / "SIMULATIONS" / "DISEGREG-CACHE"
    return Path("/u") / os.environ.get("USER", "unknown") / "PTMP" / "DISEGREG-CACHE"


def stage1_params(definitions) -> dict:
    """The parameters that determine Stage-1 content (a change gives a new cache key)."""
    return {"structure": definitions["structure"], "definition_set": definitions.label}


def params_key(params: dict, ovito_version: str, pbc=None) -> str:
    key = {"params": params["structure"], "ovito": ovito_version, "format": CACHE_FORMAT}
    if pbc is not None:  # only then: the keys of all v0.1 cache entries stay as they were
        key["pbc_override"] = [bool(x) for x in pbc]
    blob = json.dumps(key, sort_keys=True)
    return hashlib.sha256(blob.encode()).hexdigest()[:8]


def cache_path(dump, key: str, cache_root=None) -> Path:
    """<cache root>/<project>/<thread>/<run>/<dump file name>.s1-<key>.npz -- the dump's last three folder levels, so
    runs of the same name in different projects never collide."""
    dump = Path(dump).resolve()
    root = Path(cache_root) if cache_root is not None else default_cache_root()
    name = dump.name[:-3] if dump.name.endswith(".gz") else dump.name
    levels = [p.name for p in list(dump.parents)[:3] if p.name][::-1]
    return root.joinpath(*levels) / f"{name}.s1-{key}.npz"


def get_snapshot(dump, definitions, cache_root=None, recompute: bool = False, verbose: bool = False,
                 pbc=None) -> Snapshot:
    """Stage-1 data for one dump: from the cache if present and the dump is unchanged (SHA-256), else computed with
    OVITO and written to the cache. ``pbc`` overrides the file's periodicity (LAMMPS data files, see io_ovito)."""
    from . import io_ovito  # imports OVITO

    params = stage1_params(definitions)
    path = cache_path(dump, params_key(params, io_ovito.ovito_version(), pbc), cache_root)
    sha = sha256_file(dump)
    if path.is_file() and not recompute:
        snap = Snapshot.load(path)
        if snap.meta.get("source_sha256") == sha:
            if verbose:
                print(f"cache hit   {path}")
            return snap
        if verbose:
            print(f"stale cache (dump changed) {path}")
    snap = io_ovito.compute_stage1(dump, definitions, source_sha256=sha, pbc=pbc)
    snap.save(path)
    if verbose:
        print(f"computed    {path}")
    return snap
