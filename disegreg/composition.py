"""Composition quantities per snapshot (docs/DEFINITIONS.md section 6). Pure numpy.

Everything is first reduced to COUNTS (atoms and solute atoms per region, per radial bin, inside radius R). The reference
composition x_ref is usually a time average from another run (bulk control or isotherm), known only later, so the
counts are what a snapshot stores; c_k, E_k and Gamma follow from them for any x_ref.

    c_k      = N_sol,k / N_k
    E_k      = c_k / x_ref
    Gamma(R) = [N_sol(r < R) - x_ref N(r < R)] / L_x     (solute atoms per A of dislocation line, interior atoms only)
    Gamma_k  = [N_sol,k - x_ref N_k] / L_x               (the same restricted to region k)
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .geometry import Cell, densify_polyline, n_slabs_for, _periodic_centroid_unit
from .regions import REGION_CODES, radial_bin_index, radial_bin_labels


@dataclass
class Counts:
    """Atom and solute counts of one snapshot (interior atoms only, except the region X entry)."""
    L_x: float
    regions: tuple             # region codes
    n_region: np.ndarray       # (7,) atoms per region code (X included, for completeness)
    nsol_region: np.ndarray    # (7,) solute atoms per region code
    bin_edges: np.ndarray      # radial bin edges; last bin open
    bin_labels: list
    n_bin: dict                # {"r": (nbins,), "r_dxa": (nbins,)} interior atoms per radial bin
    nsol_bin: dict
    radii: np.ndarray          # R for the cumulative counts
    n_inside: dict             # {"r": (nR,), "r_dxa": (nR,)} interior atoms with distance < R
    nsol_inside: dict
    n_far_slab40: int
    nsol_far_slab40: int


def count(solute, codes, r, r_dxa, edges, radii, L_x: float, far_slab40=None) -> Counts:
    solute = np.asarray(solute, dtype=bool)
    codes = np.asarray(codes)
    interior = codes != REGION_CODES.index("X")
    n_reg = np.bincount(codes, minlength=len(REGION_CODES)).astype(np.int64)
    ns_reg = np.bincount(codes[solute], minlength=len(REGION_CODES)).astype(np.int64)
    edges = np.asarray(edges, dtype=float)
    radii = np.asarray(sorted(set(float(x) for x in radii)), dtype=float)
    n_bin, ns_bin, n_in, ns_in = {}, {}, {}, {}
    for basis, dist in (("r", r), ("r_dxa", r_dxa)):
        d = np.asarray(dist, dtype=float)[interior]
        s = solute[interior]
        b = radial_bin_index(d, edges)
        ok = b >= 0
        n_bin[basis] = np.bincount(b[ok], minlength=len(edges)).astype(np.int64)
        ns_bin[basis] = np.bincount(b[ok & s], minlength=len(edges)).astype(np.int64)
        ds = np.sort(d)
        dss = np.sort(d[s])
        n_in[basis] = np.searchsorted(ds, radii, side="left").astype(np.int64)
        ns_in[basis] = np.searchsorted(dss, radii, side="left").astype(np.int64)
    fs = np.zeros(len(solute), bool) if far_slab40 is None else np.asarray(far_slab40, dtype=bool)
    return Counts(L_x=float(L_x), regions=REGION_CODES, n_region=n_reg, nsol_region=ns_reg, bin_edges=edges,
                  bin_labels=radial_bin_labels(edges), n_bin=n_bin, nsol_bin=ns_bin, radii=radii, n_inside=n_in,
                  nsol_inside=ns_in, n_far_slab40=int(fs.sum()), nsol_far_slab40=int((fs & solute).sum()))


def concentration(nsol, n):
    """N_sol / N, NaN where N = 0."""
    nsol = np.asarray(nsol, dtype=float)
    n = np.asarray(n, dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(n > 0, nsol / np.where(n > 0, n, 1.0), np.nan)


def excess(nsol, n, x_ref: float, L_x: float):
    """[N_sol - x_ref N] / L_x: Gamma(R) from inside-R counts, Gamma_k from region counts."""
    return (np.asarray(nsol, dtype=float) - float(x_ref) * np.asarray(n, dtype=float)) / float(L_x)


# ---------------------------------------------------------------------------------------------------------------- partials
@dataclass
class PartialSeparation:
    n_partials: int            # Shockley polylines found
    n_bins: int
    n_valid: int               # X bins in which exactly two partials (groups of Shockley polylines) are present
    mean: float                # A, NaN if no valid bin
    std: float
    per_bin: np.ndarray        # (n_bins,) Y distance, NaN where not exactly two


def partial_separation(cell: Cell, shockley_polylines, bin_width: float = 5.0, spacing: float = 0.25,
                       merge_distance: float = 2.0):
    """Y distance between two Shockley partials in X bins (DEFINITIONS 6).

    In each bin the (periodic) mean position of every Shockley polyline present is taken. Polylines closer than
    ``merge_distance`` in the Y-Z plane are pieces of the same partial (a line split into several DXA segments) and are
    grouped. The bin is valid if there are exactly two groups; the separation is their minimum-image distance along Y
    (the glide-plane direction normal to the line).
    """
    n = n_slabs_for(cell, bin_width)
    lines = [densify_polyline(p, spacing) for p in shockley_polylines if len(p) >= 2]
    sb = np.full((len(lines), n), np.nan)
    sc = np.full((len(lines), n), np.nan)
    for m, pts in enumerate(lines):
        s = cell.wrap_fractional(cell.to_fractional(pts))
        k = np.minimum((s[:, 0] * n).astype(np.int64), n - 1)
        sb[m], cnt = _periodic_centroid_unit(s[:, 1], k, n)
        ok = cnt > 0
        sc[m, ok] = np.bincount(k, s[:, 2], minlength=n)[ok] / cnt[ok]
    per_bin = np.full(n, np.nan)
    H = cell.H
    for k in range(n):
        present = np.nonzero(np.isfinite(sb[:, k]))[0]
        if len(present) < 2:
            continue
        # pieces of the same partial (DXA segments meeting in this bin) lie on top of each other: group them
        ref = sb[present[0], k]
        u = np.column_stack([sb[present, k] - np.round(sb[present, k] - ref), sc[present, k]])
        groups = []
        for m in range(len(present)):
            for g in groups:
                db = u[m, 0] - u[g[0], 0]
                dv = H[:, 1] * (db - np.round(db)) + H[:, 2] * (u[m, 1] - u[g[0], 1])
                if np.linalg.norm(dv[1:]) < merge_distance:
                    g.append(m)
                    break
            else:
                groups.append([m])
        if len(groups) == 2:
            for g in groups:  # members of a group onto the image nearest its first member
                u[g, 0] -= np.round(u[g, 0] - u[g[0], 0])
            (gb1, gc1), (gb2, gc2) = (u[g].mean(axis=0) for g in groups)
            dsb = gb1 - gb2
            dsb -= np.round(dsb)
            per_bin[k] = abs(dsb * H[1, 1] + (gc1 - gc2) * H[1, 2])
    valid = np.isfinite(per_bin)
    return PartialSeparation(n_partials=len(lines), n_bins=n, n_valid=int(valid.sum()),
                             mean=float(per_bin[valid].mean()) if valid.any() else float("nan"),
                             std=float(per_bin[valid].std()) if valid.any() else float("nan"), per_bin=per_bin)
