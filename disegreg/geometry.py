"""Geometry on triclinic periodic cells -- pure numpy/scipy, no OVITO.

Implements docs/DEFINITIONS.md section 3 (distance to the line, line coverage, fallback core axis) and the periodic
neighbour search used by the order and cluster analysis (section 7).

Conventions
-----------
* A cell is given by its three edge vectors a, b, c (the COLUMNS of a 3x3 matrix, as in OVITO's ``SimulationCell``)
  and an origin. LAMMPS: a = (lx, 0, 0), b = (xy, ly, 0), c = (xz, yz, lz).
* The dislocation line runs along a (simulation X); a and b are periodic, c (glide-plane normal Z) is not.
* "Fractional" coordinates s are defined by r = origin + s_a a + s_b b + s_c c.
* Periodic images use the FULL cell vectors (including the xy tilt). scipy's periodic KD-tree only supports rectangular
  boxes, so images are built explicitly (line points) or as ghost atoms (neighbour search), then a non-periodic
  KD-tree is used.
* DXA polylines are continuous (OVITO: "the point sequence always forms a continuous line, which may lead outside the
  primary cell"). They are therefore densified in their continuous form FIRST and wrapped into the cell afterwards;
  wrapping first would create spurious chords across the cell.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree


# --------------------------------------------------------------------------------------------------------------------
# Cell
# --------------------------------------------------------------------------------------------------------------------
class Cell:
    """Triclinic simulation cell with per-axis periodicity."""

    def __init__(self, matrix, origin=None, pbc=(True, True, False)):
        m = np.asarray(matrix, dtype=float)
        if m.shape == (3, 4):  # OVITO SimulationCell matrix: columns a, b, c, origin
            if origin is None:
                origin = m[:, 3]
            m = m[:, :3]
        if m.shape != (3, 3):
            raise ValueError(f"cell matrix must be 3x3 or 3x4, got {m.shape}")
        self.H = m.copy()
        self.origin = np.zeros(3) if origin is None else np.asarray(origin, dtype=float).copy()
        self.pbc = tuple(bool(x) for x in pbc)
        if len(self.pbc) != 3:
            raise ValueError("pbc must have three entries")
        det = np.linalg.det(self.H)
        if abs(det) < 1e-12:
            raise ValueError("degenerate cell")
        self.Hinv = np.linalg.inv(self.H)

    @classmethod
    def lammps(cls, lx, ly, lz, xy=0.0, xz=0.0, yz=0.0, origin=(0.0, 0.0, 0.0), pbc=(True, True, False)):
        """Cell from LAMMPS box lengths and tilt factors."""
        H = np.array([[lx, xy, xz], [0.0, ly, yz], [0.0, 0.0, lz]], dtype=float)
        return cls(H, origin=origin, pbc=pbc)

    def __repr__(self):
        return f"Cell(a={self.a.tolist()}, b={self.b.tolist()}, c={self.c.tolist()}, pbc={self.pbc})"

    @property
    def a(self):
        return self.H[:, 0]

    @property
    def b(self):
        return self.H[:, 1]

    @property
    def c(self):
        return self.H[:, 2]

    @property
    def line_length(self) -> float:
        """L_x = length of the periodic cell vector along the line (a)."""
        return float(np.linalg.norm(self.a))

    @property
    def periodic_axes(self) -> list[int]:
        return [k for k in range(3) if self.pbc[k]]

    @property
    def volume(self) -> float:
        return float(abs(np.linalg.det(self.H)))

    def widths(self) -> np.ndarray:
        """Perpendicular widths of the cell (distance between opposite faces) for axes a, b, c."""
        a, b, c = self.a, self.b, self.c
        v = self.volume
        return np.array([v / np.linalg.norm(np.cross(b, c)),
                         v / np.linalg.norm(np.cross(c, a)),
                         v / np.linalg.norm(np.cross(a, b))])

    def to_fractional(self, r) -> np.ndarray:
        r = np.asarray(r, dtype=float)
        return (r - self.origin) @ self.Hinv.T

    def to_cartesian(self, s) -> np.ndarray:
        s = np.asarray(s, dtype=float)
        return s @ self.H.T + self.origin

    def wrap_fractional(self, s) -> np.ndarray:
        """Fractional coordinates wrapped into [0, 1) along the periodic axes (non-periodic axes untouched)."""
        s = np.array(s, dtype=float, copy=True)
        for k in self.periodic_axes:
            s[..., k] -= np.floor(s[..., k])
            s[..., k][s[..., k] >= 1.0] -= 1.0  # guards against s = -1e-17 -> 1.0
        return s

    def wrap(self, r) -> np.ndarray:
        """Cartesian positions wrapped into the primary cell along the periodic axes."""
        return self.to_cartesian(self.wrap_fractional(self.to_fractional(r)))

    def image_shifts(self, n: int = 1) -> np.ndarray:
        """Cartesian shift vectors of all periodic images with integer coefficients in [-n, n]; the zero shift first."""
        pax = self.periodic_axes
        coeffs = [np.zeros(3)]
        for combo in itertools.product(range(-n, n + 1), repeat=len(pax)):
            if any(combo):
                v = np.zeros(3)
                v[pax] = combo
                coeffs.append(v)
        return np.array(coeffs) @ self.H.T


# --------------------------------------------------------------------------------------------------------------------
# Polylines
# --------------------------------------------------------------------------------------------------------------------
def densify_polyline(points, spacing: float = 0.25) -> np.ndarray:
    """Insert points along a continuous polyline so that consecutive points are at most ``spacing`` apart.

    The original vertices are kept. With spacing h, the distance from a point at distance r from the line to the nearest
    dense point overestimates the true distance by at most h^2 / (8 r) (0.016 A for h = 0.25 A, r = 0.5 A).
    """
    p = np.asarray(points, dtype=float).reshape(-1, 3)
    if len(p) < 2:
        return p.copy()
    d = np.diff(p, axis=0)
    seg_len = np.linalg.norm(d, axis=1)
    n = np.maximum(1, np.ceil(seg_len / spacing).astype(np.int64))
    seg = np.repeat(np.arange(len(d)), n)
    t = (np.arange(n.sum()) - np.repeat(np.cumsum(n) - n, n)) / np.repeat(n, n)
    out = p[seg] + t[:, None] * d[seg]
    return np.vstack([out, p[-1:]])


def polyline_length(points) -> float:
    p = np.asarray(points, dtype=float).reshape(-1, 3)
    return float(np.linalg.norm(np.diff(p, axis=0), axis=1).sum()) if len(p) > 1 else 0.0


class LineDistance:
    """Shortest distance from atoms to a set of polylines over all periodic images (DEFINITIONS 3.1).

    Each polyline is densified in its continuous form, wrapped into the cell, replicated into the neighbouring periodic
    images by the full cell vectors, and put into one non-periodic KD-tree.
    """

    def __init__(self, cell: Cell, polylines, spacing: float = 0.25, n_images: int = 1):
        self.cell = cell
        self.spacing = float(spacing)
        dense, owner = [], []
        for k, pl in enumerate(polylines):
            d = densify_polyline(pl, spacing)
            if len(d):
                dense.append(d)
                owner.append(np.full(len(d), k, dtype=np.int64))
        self.n_lines = len(polylines)
        if not dense:
            self._tree = None
            return
        pts = cell.wrap(np.vstack(dense))
        owner = np.concatenate(owner)
        shifts = cell.image_shifts(n_images)
        self._points = (pts[None, :, :] + shifts[:, None, :]).reshape(-1, 3)
        self._owner = np.tile(owner, len(shifts))
        self._tree = cKDTree(self._points)

    @property
    def empty(self) -> bool:
        return self._tree is None

    def query(self, positions, workers: int = -1):
        """Return (r, line_index) per atom; r = inf and line_index = -1 when there is no line."""
        pos = np.asarray(positions, dtype=float).reshape(-1, 3)
        if self._tree is None:
            return np.full(len(pos), np.inf), np.full(len(pos), -1, dtype=np.int64)
        r, i = self._tree.query(self.cell.wrap(pos), k=1, workers=workers)
        return r, self._owner[i]


def line_coverage(cell: Cell, polylines, spacing: float = 0.25) -> float:
    """Projected line coverage f_cov (DEFINITIONS 3.2).

    Fraction of the cell length along the line vector a that is covered by at least one polyline: the union of the
    polylines' intervals in the fractional coordinate s_a (modulo 1). Overlapping lines (two partials) count once, so a
    fully constricted line and a fully dissociated line both give 1.
    """
    lo_all, hi_all = [], []
    for pl in polylines:
        d = densify_polyline(pl, spacing)
        if len(d) < 2:
            continue
        s = cell.to_fractional(d)[:, 0]
        lo = np.minimum(s[:-1], s[1:])
        hi = np.maximum(s[:-1], s[1:])
        shift = np.floor(lo)
        lo, hi = lo - shift, hi - shift
        over = hi > 1.0  # interval crosses the periodic boundary: split it
        lo_all += [lo, np.zeros(over.sum())]
        hi_all += [np.minimum(hi, 1.0), hi[over] - 1.0]
    if not lo_all:
        return 0.0
    lo = np.concatenate(lo_all)
    hi = np.concatenate(hi_all)
    order = np.argsort(lo, kind="stable")
    lo, hi = lo[order], hi[order]
    run_max = np.maximum.accumulate(hi)
    prev = np.concatenate([[-np.inf], run_max[:-1]])
    new = lo > prev
    starts = lo[new]
    last = np.concatenate([np.nonzero(new)[0][1:] - 1, [len(lo) - 1]])
    covered = float(np.sum(run_max[last] - starts))
    return min(covered, 1.0)


# --------------------------------------------------------------------------------------------------------------------
# Neighbour search with ghost atoms
# --------------------------------------------------------------------------------------------------------------------
def neighbor_pairs(cell: Cell, positions, cutoff: float):
    """All atom pairs closer than ``cutoff`` (minimum image, triclinic). Returns (i, j, d) with i < j, each bond once.

    Ghost atoms are made for atoms within ``cutoff`` of a periodic face (measured perpendicular to the face), then a
    non-periodic KD-tree is used. Requires cutoff < half the perpendicular width along every periodic axis.
    """
    pos = np.asarray(positions, dtype=float).reshape(-1, 3)
    n = len(pos)
    pax = cell.periodic_axes
    w = cell.widths()
    for k in pax:
        if cutoff >= 0.5 * w[k]:
            raise ValueError(f"cutoff {cutoff} >= half the cell width {w[k]:.3f} along axis {k}")
    s = cell.wrap_fractional(cell.to_fractional(pos))
    wrapped = cell.to_cartesian(s)
    pts, owner = [wrapped], [np.arange(n)]
    t = np.array([cutoff / w[k] for k in pax])
    for combo in itertools.product((-1, 0, 1), repeat=len(pax)):
        if not any(combo):
            continue
        shifted = s[:, pax] + np.array(combo)
        mask = np.all((shifted >= -t) & (shifted < 1.0 + t), axis=1)
        if mask.any():
            v = np.zeros(3)
            v[pax] = combo
            pts.append(wrapped[mask] + cell.H @ v)
            owner.append(np.nonzero(mask)[0])
    allpts = np.vstack(pts)
    owner = np.concatenate(owner)
    pairs = cKDTree(allpts).query_pairs(cutoff, output_type="ndarray")
    if len(pairs) == 0:
        return np.empty(0, np.int64), np.empty(0, np.int64), np.empty(0)
    pairs = pairs[(pairs[:, 0] < n) | (pairs[:, 1] < n)]  # drop ghost-ghost pairs
    d = np.linalg.norm(allpts[pairs[:, 0]] - allpts[pairs[:, 1]], axis=1)
    i, j = owner[pairs[:, 0]], owner[pairs[:, 1]]
    keep = i != j
    i, j, d = i[keep], j[keep], d[keep]
    lo, hi = np.minimum(i, j), np.maximum(i, j)
    _, first = np.unique(lo * n + hi, return_index=True)  # a boundary bond appears twice (i-ghost j, ghost i-j)
    return lo[first].astype(np.int64), hi[first].astype(np.int64), d[first]


# --------------------------------------------------------------------------------------------------------------------
# Fallback core axis from the structure (DEFINITIONS 3.3)
# --------------------------------------------------------------------------------------------------------------------
@dataclass
class CoreAxis:
    points: np.ndarray        # (n_slabs + 1, 3) continuous closed polyline; last point = first + a (+ winding * b)
    n_cluster_atoms: int      # atoms in the largest connected non-fcc cluster
    n_core_atoms: int         # of which core-type ("other") atoms, used for the centroids
    n_slabs: int
    n_filled_slabs: int       # slabs that contained core atoms (the others are interpolated)
    n_clusters: int           # number of connected non-fcc clusters found


def _periodic_centroid_unit(x, bin_index, n_bins):
    """Per-bin centroid of a fractional coordinate with period 1; NaN for empty bins.

    Two steps: the circular mean gives a reference; every value is then moved to its image nearest that reference and
    averaged arithmetically. This is the exact centroid of the minimum-image cloud (the circular mean alone is biased
    for broad, asymmetric distributions).
    """
    ang = 2.0 * np.pi * x
    cs = np.bincount(bin_index, np.cos(ang), minlength=n_bins)
    sn = np.bincount(bin_index, np.sin(ang), minlength=n_bins)
    cnt = np.bincount(bin_index, minlength=n_bins)
    ok = cnt > 0
    ref = np.zeros(n_bins)
    ref[ok] = np.arctan2(sn[ok], cs[ok]) / (2.0 * np.pi)
    xu = x - np.round(x - ref[bin_index])
    out = np.full(n_bins, np.nan)
    out[ok] = np.bincount(bin_index, xu, minlength=n_bins)[ok] / cnt[ok]
    return out, cnt


def fallback_core_axis(cell: Cell, positions, nonfcc_mask, core_mask, connectivity_cutoff: float,
                       slab_width: float = 5.0, smoothing_slabs: int = 3) -> CoreAxis | None:
    """Core axis from the structure, for snapshots where DXA finds too little of the line (DEFINITIONS 3.3).

    Parameters
    ----------
    nonfcc_mask : bool array, INTERIOR atoms that are not a-CNA fcc (excluded surface layers already removed).
    core_mask : bool array, atoms of core type (a-CNA "other"); only those inside the largest cluster are used.
    connectivity_cutoff : 1st-shell cut-off (0.854 a_cell in D1).

    Returns None when there is no non-fcc atom or the largest cluster contains no core-type atom.
    """
    pos = np.asarray(positions, dtype=float).reshape(-1, 3)
    nonfcc_mask = np.asarray(nonfcc_mask, dtype=bool)
    core_mask = np.asarray(core_mask, dtype=bool)
    idx = np.nonzero(nonfcc_mask)[0]
    if len(idx) == 0:
        return None
    i, j, _ = neighbor_pairs(cell, pos[idx], connectivity_cutoff)
    m = len(idx)
    graph = coo_matrix((np.ones(len(i), dtype=np.int8), (i, j)), shape=(m, m))
    n_comp, lab = connected_components(graph, directed=False)
    counts = np.bincount(lab)
    members = idx[lab == np.argmax(counts)]
    core_members = members[core_mask[members]]
    if len(core_members) == 0:
        return None

    s = cell.wrap_fractional(cell.to_fractional(pos[core_members]))
    n_slabs = max(1, int(round(cell.line_length / slab_width)))
    k = np.minimum((s[:, 0] * n_slabs).astype(np.int64), n_slabs - 1)
    sb, cnt = _periodic_centroid_unit(s[:, 1], k, n_slabs)
    sc = np.full(n_slabs, np.nan)
    ok = cnt > 0
    sc[ok] = np.bincount(k, s[:, 2], minlength=n_slabs)[ok] / cnt[ok]

    # make s_b continuous along the filled slabs; winding = net number of b periods over one a period (normally 0)
    f = np.nonzero(ok)[0]
    sb_f = np.unwrap(sb[f], period=1.0)
    sc_f = sc[f]
    # integer winding that makes the step from the last filled slab to the first one (one a-period later) smallest
    winding = float(np.round(sb_f[-1] - sb_f[0]))

    # periodic interpolation of empty slabs
    x_ext = np.concatenate([[f[-1] - n_slabs], f, [f[0] + n_slabs]])
    sb_ext = np.concatenate([[sb_f[-1] - winding], sb_f, [sb_f[0] + winding]])
    sc_ext = np.concatenate([[sc_f[-1]], sc_f, [sc_f[0]]])
    xs = np.arange(n_slabs)
    sb_full = np.interp(xs, x_ext, sb_ext)
    sc_full = np.interp(xs, x_ext, sc_ext)

    # running mean over smoothing_slabs, periodic in X
    h = max(0, int(smoothing_slabs) // 2)
    if h > 0 and n_slabs > 1:
        h = min(h, n_slabs - 1)
        ker = np.ones(2 * h + 1) / (2 * h + 1)
        sb_full = np.convolve(np.concatenate([sb_full[-h:] - winding, sb_full, sb_full[:h] + winding]), ker, "valid")
        sc_full = np.convolve(np.concatenate([sc_full[-h:], sc_full, sc_full[:h]]), ker, "valid")

    sa = (xs + 0.5) / n_slabs
    frac = np.column_stack([sa, sb_full, sc_full])
    frac = np.vstack([frac, frac[:1] + np.array([1.0, winding, 0.0])])
    return CoreAxis(points=cell.to_cartesian(frac), n_cluster_atoms=int(len(members)),
                    n_core_atoms=int(len(core_members)), n_slabs=n_slabs, n_filled_slabs=int(ok.sum()),
                    n_clusters=int(n_comp))
