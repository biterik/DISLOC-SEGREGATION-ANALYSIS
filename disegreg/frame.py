"""Crystal frame and Burgers vectors (pure numpy).

The crystal directions of the simulation axes X (line), Y, Z (glide-plane normal) come from the cell-builder metadata via
the manifest (never typed in). With F = rows of unit crystal vectors of X, Y, Z (``Manifest.frame_matrix()``), a vector
in the simulation frame v_sim has crystal components v_cr = F^T v_sim.

b and xi are taken from the SAME DXA segment: b = spatial Burgers vector (simulation frame, OVITO convention FS/RH),
xi = the segment's line direction (principal direction of its points, oriented from the first to the last point).
Character angle = angle between b and xi folded into [0, 90] degrees: 0 screw, 90 edge.
"""

from __future__ import annotations

from math import gcd

import numpy as np

# Burgers-vector families: sorted absolute components of DXA's true Burgers vector (units of the lattice parameter)
FAMILIES = {
    "1/6<112>": (1 / 6, 1 / 6, 1 / 3),  # Shockley partial
    "1/2<110>": (0.0, 0.5, 0.5),        # perfect
    "1/6<110>": (0.0, 1 / 6, 1 / 6),    # stair-rod
    "1/3<100>": (0.0, 0.0, 1 / 3),      # Hirth
    "1/3<111>": (1 / 3, 1 / 3, 1 / 3),  # Frank
}
SHOCKLEY = "1/6<112>"
PERFECT = "1/2<110>"


def burgers_family(b_true, tol: float = 1e-3) -> str:
    """Family name of a true Burgers vector (any sign / permutation), "other" if none matches."""
    s = np.sort(np.abs(np.asarray(b_true, dtype=float)))
    for name, ref in FAMILIES.items():
        if np.allclose(s, ref, atol=tol):
            return name
    return "other"


def to_crystal(v_sim, frame_matrix) -> np.ndarray:
    """Simulation-frame vector(s) -> crystal components (same length unit)."""
    return np.asarray(v_sim, dtype=float) @ np.asarray(frame_matrix, dtype=float)


def _fmt(ints) -> str:
    return "[" + " ".join(str(int(i)) for i in ints) + "]"


def burgers_label(b_cr_in_a, tol: float = 0.03) -> str:
    """Miller label of a Burgers vector in units of a, e.g. "1/6[1 -2 1]", "1/2[0 -1 1]"; "" if no denominator in
    1, 2, 3, 6 gives integers within ``tol`` (in units of a)."""
    v = np.asarray(b_cr_in_a, dtype=float)
    for d in (1, 2, 3, 6):
        w = v * d
        r = np.round(w)
        if np.all(np.abs(w - r) < tol * d) and r.any():
            g = gcd(gcd(int(abs(r[0])), int(abs(r[1]))), int(abs(r[2])))
            g = gcd(g, d)
            r, dd = r / g, d // g
            return _fmt(r) if dd == 1 else f"1/{dd}{_fmt(r)}"
    return ""


def direction_label(v, max_index: int = 6, tol: float = 0.03) -> str:
    """Smallest integer direction [u v w] (|components| <= max_index) parallel to v within ``tol`` (radians); "" if none."""
    v = np.asarray(v, dtype=float)
    u = v / np.linalg.norm(v)
    best = None
    rng = range(-max_index, max_index + 1)
    for cand in ((i, j, k) for i in rng for j in rng for k in rng):
        c = np.array(cand, float)
        nc = np.linalg.norm(c)
        if nc == 0 or gcd(gcd(abs(cand[0]), abs(cand[1])), abs(cand[2])) != 1:
            continue
        ang = np.arccos(np.clip(u @ c / nc, -1, 1))
        if ang < tol and (best is None or (np.abs(c).sum(), ang) < best[0]):
            best = ((np.abs(c).sum(), ang), cand)
    return _fmt(best[1]) if best else ""


def line_direction(points) -> np.ndarray:
    """Unit line direction of a polyline: principal axis of its points, oriented first -> last point."""
    p = np.asarray(points, dtype=float).reshape(-1, 3)
    if len(p) < 2:
        return np.full(3, np.nan)
    _, _, vt = np.linalg.svd(p - p.mean(axis=0), full_matrices=False)
    d = vt[0]
    if d @ (p[-1] - p[0]) < 0:
        d = -d
    return d


def character_angle(b, xi) -> float:
    """Angle between b and xi folded into [0, 90] degrees (0 screw, 90 edge)."""
    b = np.asarray(b, float)
    xi = np.asarray(xi, float)
    c = abs(b @ xi) / (np.linalg.norm(b) * np.linalg.norm(xi))
    return float(np.degrees(np.arccos(np.clip(c, 0.0, 1.0))))


def segment_table(snap, frame_matrix, a_cell: float) -> list[dict]:
    """One dict per DXA segment: family, length, b and xi in the crystal frame with Miller labels, character angle."""
    out = []
    for k in range(snap.n_segments):
        pts = snap.seg_points[snap.seg_offsets[k]:snap.seg_offsets[k + 1]]
        xi = line_direction(pts)
        b_sim = snap.seg_spatial_b[k]
        b_cr = to_crystal(b_sim, frame_matrix) / a_cell
        xi_cr = to_crystal(xi, frame_matrix)
        out.append({"segment": k, "family": str(snap.seg_family[k]), "length": float(snap.seg_length[k]),
                    "infinite": bool(snap.seg_infinite[k]), "b_sim": [float(x) for x in b_sim],
                    "b_crystal_a": [round(float(x), 5) for x in b_cr], "b_label": burgers_label(b_cr),
                    "xi_crystal": [round(float(x), 5) for x in xi_cr], "xi_label": direction_label(xi_cr),
                    "character_angle": character_angle(b_sim, xi) if np.all(np.isfinite(xi)) else float("nan")})
    return out


def total_burgers(table, x_axis=(1.0, 0.0, 0.0)) -> np.ndarray:
    """Sum of the crystal-frame b (units of a) of all infinite segments, each with its line sense along +X: the Burgers
    vector of the whole (dissociated) dislocation threading the cell."""
    tot = np.zeros(3)
    F_x = np.asarray(x_axis, float)
    for r in table:
        if r["infinite"]:
            sense = np.sign(np.asarray(r["xi_crystal"]) @ F_x) or 1.0
            tot += sense * np.asarray(r["b_crystal_a"])
    return tot
