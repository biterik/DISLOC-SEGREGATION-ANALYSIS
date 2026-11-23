"""Exclusive region codes (docs/DEFINITIONS.md section 4) and the glide-plane height check (section 3.4). Pure numpy.

Every atom gets exactly one code, so region excesses add up:

    X excluded    within ``surface_exclusion`` of the lowest or highest atom in Z
    C core        r < core_radius
    F fault       a-CNA hcp, |z - z_gp| < fault_half_width, r >= core_radius
    S shell       core_radius <= r < shell_radius, not F
    T transition  shell_radius <= r < far_field_radius, not F
    B far field   r >= far_field_radius, a-CNA fcc
    O far, other  r >= far_field_radius, not fcc, not F

Precedence X > C > F > S, T > B, O. The table in DEFINITIONS 4 does not say whether an hcp atom near the glide plane at
r >= R_far is F or O; here it is F (as for S and T). In practice the fault ribbon ends well inside R_far.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .cache import STRUCTURE_CODES
from .geometry import Cell, densify_polyline

REGION_CODES = ("X", "C", "F", "S", "T", "B", "O")
REGION_NAMES = {"X": "excluded", "C": "core", "F": "fault", "S": "shell", "T": "transition", "B": "far field",
                "O": "far, other"}
INTERIOR_CODES = REGION_CODES[1:]
_IDX = {c: k for k, c in enumerate(REGION_CODES)}


def code_index(code: str) -> int:
    return _IDX[code]


def excluded_mask(positions, surface_exclusion: float):
    """Atoms within ``surface_exclusion`` of the lowest or highest atom in Z (region X). Returns (mask, z_lo, z_hi)."""
    z = np.asarray(positions, dtype=float).reshape(-1, 3)[:, 2]
    z_lo, z_hi = float(z.min()), float(z.max())
    return (z - z_lo < surface_exclusion) | (z_hi - z < surface_exclusion), z_lo, z_hi


def z_gp_from_lines(cell: Cell, shockley_polylines, axis_polylines=(), spacing: float = 0.25):
    """Glide-plane height from the structure: median Z of the (densified, i.e. length-weighted) Shockley line points,
    else of the fallback axes. Returns (z, source) with source "shockley" | "axis" | "none"."""
    for lines, source in ((shockley_polylines, "shockley"), (axis_polylines, "axis")):
        pts = [densify_polyline(p, spacing) for p in lines if len(p) >= 2]
        if pts:
            return float(np.median(np.vstack(pts)[:, 2])), source
    return float("nan"), "none"


@dataclass
class GlidePlane:
    z_gp: float            # value used for the regions
    source: str            # "manifest" | "lines" (no manifest value given)
    z_lines: float         # from the structure (NaN if neither Shockley lines nor axes)
    lines_source: str      # "shockley" | "axis" | "none"
    difference: float      # z_manifest - z_lines (NaN if one is missing)
    tolerance: float       # z_gp_tolerance_planes x a/sqrt(3)
    flag: bool             # |difference| > tolerance, or no value at all


def glide_plane(z_manifest, z_lines: float, lines_source: str, a_cell: float, tolerance_planes: float) -> GlidePlane:
    tol = tolerance_planes * a_cell / np.sqrt(3.0)
    if z_manifest is None:
        return GlidePlane(z_gp=z_lines, source="lines", z_lines=z_lines, lines_source=lines_source,
                          difference=float("nan"), tolerance=tol, flag=True)
    diff = float(z_manifest) - z_lines
    return GlidePlane(z_gp=float(z_manifest), source="manifest", z_lines=z_lines, lines_source=lines_source,
                      difference=diff, tolerance=tol, flag=bool(np.isnan(diff) or abs(diff) > tol))


def assign_regions(positions, acna, r, excluded, z_gp: float, params: dict) -> np.ndarray:
    """Region code index (into REGION_CODES) for every atom; ``params`` is the 'regions' section of a definition set."""
    z = np.asarray(positions, dtype=float).reshape(-1, 3)[:, 2]
    acna = np.asarray(acna)
    r = np.asarray(r, dtype=float)
    rc, rs, rf = float(params["core_radius"]), float(params["shell_radius"]), float(params["far_field_radius"])
    fcc = acna == STRUCTURE_CODES["fcc"]
    fault = (acna == STRUCTURE_CODES["hcp"]) & (np.abs(z - z_gp) < float(params["fault_half_width"])) & (r >= rc)
    code = np.full(len(r), _IDX["T"], dtype=np.int8)
    code[r < rs] = _IDX["S"]
    far = r >= rf
    code[far & fcc] = _IDX["B"]
    code[far & ~fcc] = _IDX["O"]
    code[fault] = _IDX["F"]
    code[r < rc] = _IDX["C"]
    code[np.asarray(excluded, dtype=bool)] = _IDX["X"]
    return code


def far_slab40_mask(positions, acna, excluded, z_gp: float, half_width: float) -> np.ndarray:
    """Optional column (DEFINITIONS 4): a-CNA fcc, |z - z_gp| > half_width, not X."""
    z = np.asarray(positions, dtype=float).reshape(-1, 3)[:, 2]
    return (np.asarray(acna) == STRUCTURE_CODES["fcc"]) & (np.abs(z - z_gp) > half_width) & ~np.asarray(excluded, bool)


def radial_bin_labels(edges) -> list[str]:
    e = [float(x) for x in edges]
    return [f"{lo:g}-{hi:g}" for lo, hi in zip(e[:-1], e[1:])] + [f">={e[-1]:g}"]


def radial_bin_index(r, edges) -> np.ndarray:
    """Bin index for [e0, e1), [e1, e2), ..., [e_last, inf); -1 below e0."""
    return np.searchsorted(np.asarray(edges, dtype=float), np.asarray(r, dtype=float), side="right") - 1
