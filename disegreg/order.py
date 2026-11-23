"""Chemical order per snapshot, per region of the central solute atom (docs/DEFINITIONS.md section 7). Pure numpy/scipy.

    shells        1st: r < 0.854 a (12 neighbours in fcc), 2nd: 0.854 a <= r < 1.112 a (6 neighbours)
    p_k           solute fraction among the k-th-shell neighbours of the solute atoms of a region (pooled:
                  all solute neighbours / all neighbours of those solute atoms)
    alpha_k       Warren-Cowley parameter (p_k - c_ref) / (1 - c_ref), with c_ref = x_ref and with c_ref = c of the region
    R2            <solute 2nd neighbours per solute> / (6 c_region); 1 random, 0 complete avoidance, L1_2: 4, D0_22: 2.7
    clusters      solute atoms connected by 1st-shell bonds; every solute atom carries the size of its cluster, aggregated
                  by the region of that atom

Neighbours are all atoms (also those of region X); the central solute atoms of region X are not reported.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from .geometry import Cell, neighbor_pairs
from .regions import REGION_CODES

CLUSTER_SIZE_BINS = (1, 2, 3, 4, 5, 10, 20, 50)  # histogram of cluster sizes per solute atom: [1], [2], ..., [50, inf)


@dataclass
class Shells:
    i: np.ndarray       # pair index arrays (i < j), each bond once
    j: np.ndarray
    shell: np.ndarray   # 1 or 2


def shell_pairs(cell: Cell, positions, a_cell: float, shell1_cutoff_a: float = 0.854,
                shell2_cutoff_a: float = 1.112) -> Shells:
    c1, c2 = shell1_cutoff_a * a_cell, shell2_cutoff_a * a_cell
    i, j, d = neighbor_pairs(cell, positions, c2)
    return Shells(i=i, j=j, shell=np.where(d < c1, 1, 2).astype(np.int8))


def neighbour_counts(n_atoms: int, solute, shells: Shells):
    """Per atom and shell: number of neighbours and number of solute neighbours. Returns (n1, s1, n2, s2)."""
    solute = np.asarray(solute, dtype=bool)
    out = []
    for k in (1, 2):
        m = shells.shell == k
        i, j = shells.i[m], shells.j[m]
        n = np.bincount(i, minlength=n_atoms) + np.bincount(j, minlength=n_atoms)
        s = np.bincount(i, weights=solute[j], minlength=n_atoms) + np.bincount(j, weights=solute[i], minlength=n_atoms)
        out += [n.astype(np.int64), s.astype(np.int64)]
    return tuple(out)


def warren_cowley(p, c_ref):
    """(p - c_ref) / (1 - c_ref); NaN where c_ref is not in (0, 1) (alpha undefined for a pure region)."""
    p = np.asarray(p, dtype=float)
    c_ref = np.broadcast_to(np.asarray(c_ref, dtype=float), p.shape)
    ok = (c_ref > 0) & (c_ref < 1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(ok, (p - c_ref) / np.where(ok, 1.0 - c_ref, 1.0), np.nan)


@dataclass
class OrderByRegion:
    regions: tuple
    n_solute: np.ndarray        # central solute atoms per region
    c_region: np.ndarray        # solute concentration of the region (all atoms of the region)
    z1: np.ndarray              # mean 1st-shell coordination of the central solute atoms
    z2: np.ndarray
    p1: np.ndarray
    p2: np.ndarray
    alpha1_local: np.ndarray    # c_ref = c_region
    alpha2_local: np.ndarray
    alpha1_xref: np.ndarray     # c_ref = x_ref (NaN if x_ref is not given)
    alpha2_xref: np.ndarray
    R2: np.ndarray


def order_by_region(solute, codes, counts, x_ref=None) -> OrderByRegion:
    """Order parameters per region code; ``counts`` = neighbour_counts(...)."""
    solute = np.asarray(solute, dtype=bool)
    codes = np.asarray(codes)
    n1, s1, n2, s2 = counts
    nreg = len(REGION_CODES)
    n_all = np.bincount(codes, minlength=nreg).astype(float)
    cs = codes[solute]

    def tot(w):
        return np.bincount(cs, weights=w[solute], minlength=nreg).astype(float)

    nsol = np.bincount(cs, minlength=nreg).astype(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        c_reg = np.where(n_all > 0, nsol / n_all, np.nan)
        p1 = tot(s1) / tot(n1)
        p2 = tot(s2) / tot(n2)
        z1 = tot(n1) / nsol
        z2 = tot(n2) / nsol
        R2 = (tot(s2) / nsol) / (6.0 * c_reg)
    xr = np.nan if x_ref is None else float(x_ref)
    return OrderByRegion(regions=REGION_CODES, n_solute=nsol.astype(np.int64), c_region=c_reg, z1=z1, z2=z2, p1=p1, p2=p2,
                         alpha1_local=warren_cowley(p1, c_reg), alpha2_local=warren_cowley(p2, c_reg),
                         alpha1_xref=warren_cowley(p1, xr), alpha2_xref=warren_cowley(p2, xr), R2=R2)


# ---------------------------------------------------------------------------------------------------------------- clusters
def solute_cluster_sizes(n_atoms: int, solute, shells: Shells) -> np.ndarray:
    """Size of the 1st-shell-connected solute cluster of every atom (0 for host atoms)."""
    solute = np.asarray(solute, dtype=bool)
    idx = np.nonzero(solute)[0]
    size = np.zeros(n_atoms, dtype=np.int64)
    if len(idx) == 0:
        return size
    m = (shells.shell == 1) & solute[shells.i] & solute[shells.j]
    local = np.full(n_atoms, -1, dtype=np.int64)
    local[idx] = np.arange(len(idx))
    g = coo_matrix((np.ones(int(m.sum()), np.int8), (local[shells.i[m]], local[shells.j[m]])), shape=(len(idx),) * 2)
    _, lab = connected_components(g, directed=False)
    size[idx] = np.bincount(lab)[lab]
    return size


@dataclass
class ClustersByRegion:
    regions: tuple
    n_solute: np.ndarray
    mean_size: np.ndarray        # mean cluster size per solute atom (solute-weighted)
    frac_isolated: np.ndarray    # fraction of solute atoms without solute 1st neighbours
    max_size: np.ndarray
    size_bins: tuple
    histogram: np.ndarray        # (7, len(size_bins)) solute atoms per region and size bin


def clusters_by_region(size, solute, codes) -> ClustersByRegion:
    solute = np.asarray(solute, dtype=bool)
    nreg = len(REGION_CODES)
    cs = np.asarray(codes)[solute]
    sz = np.asarray(size)[solute]
    n = np.bincount(cs, minlength=nreg)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.bincount(cs, weights=sz, minlength=nreg) / n
        iso = np.bincount(cs, weights=(sz == 1), minlength=nreg) / n
    mx = np.zeros(nreg, dtype=np.int64)
    np.maximum.at(mx, cs, sz)
    b = np.searchsorted(np.asarray(CLUSTER_SIZE_BINS), sz, side="right") - 1
    hist = np.zeros((nreg, len(CLUSTER_SIZE_BINS)), dtype=np.int64)
    np.add.at(hist, (cs, b), 1)
    return ClustersByRegion(regions=REGION_CODES, n_solute=n.astype(np.int64), mean_size=mean, frac_isolated=iso,
                            max_size=mx, size_bins=CLUSTER_SIZE_BINS, histogram=hist)
