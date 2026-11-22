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


def _pair_intervals(cell: Cell, dense: np.ndarray):
    """Intervals [lo, hi] in s_a (shifted so lo is in [0, 1)) of the consecutive point pairs of one continuous polyline."""
    s = cell.to_fractional(dense)[:, 0]
    lo = np.minimum(s[:-1], s[1:])
    hi = np.maximum(s[:-1], s[1:])
    shift = np.floor(lo)
    return lo - shift, hi - shift


def _merge_intervals(lo, hi):
    """Union of intervals [lo, hi] with lo in [0, 1), periodic with period 1 -> sorted disjoint (starts, ends) in [0, 1]."""
    lo = np.asarray(lo, dtype=float)
    hi = np.asarray(hi, dtype=float)
    if len(lo) == 0:
        return np.empty(0), np.empty(0)
    over = hi > 1.0  # crosses the periodic boundary: split
    lo = np.concatenate([lo, np.zeros(over.sum())])
    hi = np.concatenate([np.minimum(hi, 1.0), hi[over] - 1.0])
    order = np.argsort(lo, kind="stable")
    lo, hi = lo[order], hi[order]
    run_max = np.maximum.accumulate(hi)
    prev = np.concatenate([[-np.inf], run_max[:-1]])
    new = lo > prev
    last = np.concatenate([np.nonzero(new)[0][1:] - 1, [len(lo) - 1]])
    return lo[new], run_max[last]


def _intervals_of(cell, polylines, spacing):
    lo, hi = [], []
    for pl in polylines:
        d = densify_polyline(pl, spacing)
        if len(d) >= 2:
            a, b = _pair_intervals(cell, d)
            lo.append(a)
            hi.append(b)
    if not lo:
        return np.empty(0), np.empty(0)
    return np.concatenate(lo), np.concatenate(hi)


def line_coverage(cell: Cell, polylines, spacing: float = 0.25) -> float:
    """Projected line coverage f_cov (DEFINITIONS 3.2).

    Fraction of the cell length along the line vector a that is covered by at least one polyline: the union of the
    polylines' intervals in the fractional coordinate s_a (modulo 1). Overlapping lines (two partials) count once, so a
    fully constricted line and a fully dissociated line both give 1.
    """
    starts, ends = _merge_intervals(*_intervals_of(cell, polylines, spacing))
    return min(float(np.sum(ends - starts)), 1.0)


def slab_coverage_from_intervals(starts, ends, n_slabs: int) -> np.ndarray:
    """Covered fraction of each of n_slabs equal slabs in s_a, from disjoint intervals (starts, ends) in [0, 1]."""
    edges = np.arange(n_slabs + 1) / n_slabs
    cov = np.zeros(n_slabs)
    for s0, s1 in zip(starts, ends):
        cov += np.clip(np.minimum(edges[1:], s1) - np.maximum(edges[:-1], s0), 0.0, None)
    return np.minimum(cov * n_slabs, 1.0)


def slab_coverage(cell: Cell, polylines, n_slabs: int, spacing: float = 0.25) -> np.ndarray:
    """Per-slab projected coverage of polylines (slabs = n_slabs equal intervals of s_a)."""
    return slab_coverage_from_intervals(*_merge_intervals(*_intervals_of(cell, polylines, spacing)), n_slabs)


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
# Fallback core axes from the structure, and the distance actually used (DEFINITIONS 3.3, 3.4, rev. 3)
# --------------------------------------------------------------------------------------------------------------------
def n_slabs_for(cell: Cell, slab_width: float) -> int:
    """Number of X slabs: round(L_x / slab_width), at least 1 (actual width L_x / n)."""
    return max(1, int(round(cell.line_length / slab_width)))


@dataclass
class CoreAxis:
    """Axis of one line defect cluster: per-slab centroid (fractional s_b, s_c) at slab centres s_a = (k + 0.5) / n."""
    sb: np.ndarray            # (n_slabs,) continuous (unwrapped) fractional b coordinate
    sc: np.ndarray            # (n_slabs,) fractional c coordinate
    winding: float            # integer: net b periods per a period (normally 0)
    members: np.ndarray       # indices (into the positions passed in) of all atoms of the cluster
    n_core_atoms: int         # core-type atoms of the cluster, used for the centroids
    n_filled_slabs: int       # slabs that contained core atoms (the others are interpolated)

    @property
    def n_slabs(self) -> int:
        return len(self.sb)

    def fractional_at(self, sa) -> np.ndarray:
        """Axis position (s_a, s_b, s_c) at arbitrary (unwrapped) s_a, periodic linear interpolation."""
        sa = np.atleast_1d(np.asarray(sa, dtype=float))
        n = self.n_slabs
        u = sa * n - 0.5
        k0 = np.floor(u).astype(np.int64)
        t = u - k0

        def val(arr, k, w):
            return arr[np.mod(k, n)] + w * np.floor_divide(k, n)

        sb = (1 - t) * val(self.sb, k0, self.winding) + t * val(self.sb, k0 + 1, self.winding)
        sc = (1 - t) * val(self.sc, k0, 0.0) + t * val(self.sc, k0 + 1, 0.0)
        return np.column_stack([sa, sb, sc])

    def polyline(self, cell: Cell) -> np.ndarray:
        """Closed continuous polyline through the slab centres; last point = first + a (+ winding * b)."""
        n = self.n_slabs
        return cell.to_cartesian(self.fractional_at((np.arange(n + 1) + 0.5) / n))


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


def _axis_from_core_atoms(cell, core_pos, n_slabs, smoothing_slabs):
    s = cell.wrap_fractional(cell.to_fractional(core_pos))
    k = np.minimum((s[:, 0] * n_slabs).astype(np.int64), n_slabs - 1)
    sb, cnt = _periodic_centroid_unit(s[:, 1], k, n_slabs)
    ok = cnt > 0
    sc = np.full(n_slabs, np.nan)
    sc[ok] = np.bincount(k, s[:, 2], minlength=n_slabs)[ok] / cnt[ok]

    f = np.nonzero(ok)[0]
    sb_f = np.unwrap(sb[f], period=1.0)  # continuous along the filled slabs
    sc_f = sc[f]
    winding = float(np.round(sb_f[-1] - sb_f[0]))  # makes the step last -> first (one a-period later) smallest

    x_ext = np.concatenate([[f[-1] - n_slabs], f, [f[0] + n_slabs]])  # periodic interpolation of empty slabs
    sb_ext = np.concatenate([[sb_f[-1] - winding], sb_f, [sb_f[0] + winding]])
    sc_ext = np.concatenate([[sc_f[-1]], sc_f, [sc_f[0]]])
    xs = np.arange(n_slabs)
    sb_full = np.interp(xs, x_ext, sb_ext)
    sc_full = np.interp(xs, x_ext, sc_ext)

    h = max(0, int(smoothing_slabs) // 2)  # running mean, periodic in X
    if h > 0 and n_slabs > 1:
        h = min(h, n_slabs - 1)
        ker = np.ones(2 * h + 1) / (2 * h + 1)
        sb_full = np.convolve(np.concatenate([sb_full[-h:] - winding, sb_full, sb_full[:h] + winding]), ker, "valid")
        sc_full = np.convolve(np.concatenate([sc_full[-h:], sc_full, sc_full[:h]]), ker, "valid")
    return sb_full, sc_full, winding, int(ok.sum()), k


def line_defect_axes(cell: Cell, positions, nonfcc_mask, core_mask, connectivity_cutoff: float,
                     slab_width: float = 5.0, smoothing_slabs: int = 3, min_slab_fraction: float = 0.5):
    """Fallback core axes from the structure (DEFINITIONS 3.3, rev. 3).

    Parameters
    ----------
    nonfcc_mask : bool array, INTERIOR atoms that are not a-CNA fcc (excluded surface layers already removed).
    core_mask : bool array, atoms of core type (a-CNA "other").
    connectivity_cutoff : 1st-shell cut-off (0.854 a_cell in D1).

    Returns
    -------
    axes : list of CoreAxis, one per connected non-fcc cluster whose core atoms occupy >= min_slab_fraction of the
           slabs, sorted by cluster size (largest first).
    atom_axis : int array (len(positions)), index into ``axes`` for atoms of those clusters, -1 otherwise.
    n_clusters : number of connected non-fcc clusters found (all sizes).
    """
    pos = np.asarray(positions, dtype=float).reshape(-1, 3)
    nonfcc_mask = np.asarray(nonfcc_mask, dtype=bool)
    core_mask = np.asarray(core_mask, dtype=bool)
    atom_axis = np.full(len(pos), -1, dtype=np.int64)
    idx = np.nonzero(nonfcc_mask)[0]
    if len(idx) == 0:
        return [], atom_axis, 0
    i, j, _ = neighbor_pairs(cell, pos[idx], connectivity_cutoff)
    m = len(idx)
    graph = coo_matrix((np.ones(len(i), dtype=np.int8), (i, j)), shape=(m, m))
    n_comp, lab = connected_components(graph, directed=False)
    sizes = np.bincount(lab)
    n_slabs = n_slabs_for(cell, slab_width)

    # slab occupancy of the core atoms of every cluster, in one pass
    core_sel = core_mask[idx]
    s_core = cell.wrap_fractional(cell.to_fractional(pos[idx[core_sel]]))
    k_core = np.minimum((s_core[:, 0] * n_slabs).astype(np.int64), n_slabs - 1)
    occupied = np.zeros((n_comp, n_slabs), dtype=bool)
    occupied[lab[core_sel], k_core] = True
    qualifies = occupied.sum(axis=1) >= min_slab_fraction * n_slabs

    axes = []
    for c in sorted(np.nonzero(qualifies)[0], key=lambda c: -sizes[c]):
        members = idx[lab == c]
        core_members = members[core_mask[members]]
        sb, sc, winding, n_filled, _ = _axis_from_core_atoms(cell, pos[core_members], n_slabs, smoothing_slabs)
        atom_axis[members] = len(axes)
        axes.append(CoreAxis(sb=sb, sc=sc, winding=winding, members=members, n_core_atoms=int(len(core_members)),
                             n_filled_slabs=n_filled))
    return axes, atom_axis, int(n_comp)


def assign_points_to_axes(cell: Cell, points, positions, nonfcc_mask, atom_axis, cutoff: float, n_images: int = 1):
    """Axis index of the nearest interior non-fcc atom for every point, -1 if that atom is farther than ``cutoff``
    or belongs to no line defect cluster (DEFINITIONS 3.4, rev. 3)."""
    pts = np.asarray(points, dtype=float).reshape(-1, 3)
    out = np.full(len(pts), -1, dtype=np.int64)
    idx = np.nonzero(np.asarray(nonfcc_mask, dtype=bool))[0]
    if len(idx) == 0 or len(pts) == 0:
        return out
    shifts = cell.image_shifts(n_images)
    atoms = cell.wrap(np.asarray(positions, dtype=float)[idx])
    tree = cKDTree((atoms[None] + shifts[:, None]).reshape(-1, 3))
    d, k = tree.query(cell.wrap(pts), k=1, distance_upper_bound=cutoff)
    hit = np.isfinite(d)
    out[hit] = np.asarray(atom_axis)[idx[k[hit] % len(idx)]]
    return out


@dataclass
class LineDistances:
    r: np.ndarray                  # distance actually used (3.4)
    r_dxa: np.ndarray              # distance to the DXA polylines (3.1); inf without DXA lines
    r_axis: np.ndarray             # distance to all fallback axes, whole length (diagnostic); inf without axes
    mode: str                      # "polyline" | "mixed" | "axis" | "none"
    f_cov: float                   # projected coverage of all DXA lines (3.2)
    f_cov_axis: list               # per axis: projected coverage by the DXA points assigned to it
    active: list                   # per axis: bool array of active slabs
    n_slabs: int
    axes: list                     # CoreAxis objects
    n_clusters: int                # connected non-fcc clusters of any size


def line_distances(cell: Cell, positions, dxa_polylines, nonfcc_mask, core_mask, a_cell: float, *,
                   spacing: float = 0.25, n_images: int = 1, slab_width: float = 5.0, smoothing_slabs: int = 3,
                   connectivity_cutoff_a: float = 0.854, min_slab_fraction: float = 0.5,
                   dxa_assign_cutoff_a: float = 0.854, active_below_dxa_coverage: float = 0.5,
                   workers: int = -1) -> LineDistances:
    """All distances of DEFINITIONS 3.1-3.4 (rev. 3) for one snapshot.

    nonfcc_mask must contain INTERIOR non-fcc atoms only; core_mask marks a-CNA "other" atoms.
    """
    pos = np.asarray(positions, dtype=float).reshape(-1, 3)
    dxa_polylines = [np.asarray(p, dtype=float).reshape(-1, 3) for p in dxa_polylines]
    ld = LineDistance(cell, dxa_polylines, spacing=spacing, n_images=n_images)
    r_dxa, _ = ld.query(pos, workers=workers)
    f_cov = line_coverage(cell, dxa_polylines, spacing)
    n_slabs = n_slabs_for(cell, slab_width)

    axes, atom_axis, n_clusters = line_defect_axes(cell, pos, nonfcc_mask, core_mask, connectivity_cutoff_a * a_cell,
                                                   slab_width, smoothing_slabs, min_slab_fraction)
    f_cov_axis, active = [], []
    if axes:
        # assign every dense DXA point to an axis; a point pair counts for axis m when both points belong to m
        dense = [densify_polyline(p, spacing) for p in dxa_polylines if len(p) >= 2]
        lab = [assign_points_to_axes(cell, d, pos, nonfcc_mask, atom_axis, dxa_assign_cutoff_a * a_cell, n_images)
               for d in dense]
        for m, ax in enumerate(axes):
            lo, hi = [], []
            for d, l in zip(dense, lab):
                both = (l[:-1] == m) & (l[1:] == m)
                if both.any():
                    a_, b_ = _pair_intervals(cell, d)
                    lo.append(a_[both])
                    hi.append(b_[both])
            starts, ends = _merge_intervals(np.concatenate(lo) if lo else [], np.concatenate(hi) if hi else [])
            f_cov_axis.append(min(float(np.sum(ends - starts)), 1.0))
            active.append(slab_coverage_from_intervals(starts, ends, n_slabs) < active_below_dxa_coverage)
    # r_axis: distance to all axes (diagnostic). r: an atom whose X lies in an active slab of axis m also counts its
    # distance to axis m; everywhere else r = r_DXA, so the axis never reaches into slabs where DXA found the line.
    r_axis = np.full(len(pos), np.inf)
    r = r_dxa.copy()
    k_atom = np.minimum((cell.wrap_fractional(cell.to_fractional(pos))[:, 0] * n_slabs).astype(np.int64), n_slabs - 1)
    for ax, act in zip(axes, active):
        r_m, _ = LineDistance(cell, [ax.polyline(cell)], spacing=spacing, n_images=n_images).query(pos, workers=workers)
        r_axis = np.minimum(r_axis, r_m)
        use = act[k_atom]
        r[use] = np.minimum(r[use], r_m[use])

    if ld.empty:
        mode = "axis" if axes else "none"
    else:
        mode = "mixed" if any(a.any() for a in active) else "polyline"
    return LineDistances(r=r, r_dxa=r_dxa, r_axis=r_axis, mode=mode, f_cov=f_cov, f_cov_axis=f_cov_axis,
                         active=active, n_slabs=n_slabs, axes=axes, n_clusters=n_clusters)
