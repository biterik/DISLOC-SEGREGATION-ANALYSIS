"""Ordered-domain sizes ([PROPOSAL v0.2] P5), for Ni-Al against the L1_2 domains of Howard et al., Acta Mater. 289 (2025)
120887 (Ni-13Al after 3 h at 1050 C: 2.7 +- 0.6 nm, densest at dislocations). Pure numpy/scipy.

    ordered solute   a solute atom with at most ``max_solute_1st`` (0) solute 1st neighbours and at least ``min_solute_2nd``
                     (4) of its 6 2nd neighbours solute. In L1_2 (Ni3Al: Al on the cube corners) every Al has 0 Al 1st and
                     6 Al 2nd neighbours; in D0_22 4 of 6; in a random solution at x = 0.13 only ~0.3 % of the solute atoms
                     qualify (binomial, 6 trials, >= 4 hits) and ~19 % have no solute 1st neighbour.
    domain           ordered solute atoms connected by 2nd-shell bonds (the solute sublattice of L1_2 is simple cubic with
                     spacing a, i.e. 2nd-shell neighbours). Domains with fewer than ``min_domain_solute`` atoms are counted
                     separately (``small``), not as domains.
    border           with ``include_border`` (option, recommended in the proposal): the solute atoms that are 2nd neighbours
                     of a domain's ordered atoms but not ordered themselves (the domain surface: fewer than 4 solute 2nd
                     neighbours) belong to that domain (to the largest one if they touch several). Without them a sphere of
                     L1_2 comes out too small: d / D = 0.80 at D = 1.6 nm, 0.84 at 2.0 nm, 0.94 at 2.7 nm; with them
                     0.99-1.05 (tests/test_domains.py).
    size             equivalent diameter d = (6 V / pi)^(1/3), V = n_solute x ``volume_per_solute_a3`` x a^3 (L1_2: one
                     solute atom per cubic cell, 1 a^3).
    where            every domain belongs to the region (D1 codes) of the majority of its atoms (ties: the innermost, C first);
                     per region: number of domains, number-mean and largest d, fraction of the region's solute that is ordered.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from .regions import REGION_CODES

QUANTITIES = ("domain_frac_ordered", "domain_n", "domain_d_mean", "domain_d_max", "domain_n_small", "domain_d_all_mean",
              "domain_n_all")


@dataclass
class Domains:
    ordered: np.ndarray        # (N,) bool
    label: np.ndarray          # (N,) domain index of every ordered atom, -1 otherwise
    size: np.ndarray           # (n_domains,) ordered solute atoms per domain (incl. small ones)
    diameter: np.ndarray       # (n_domains,) equivalent diameter (A)
    region: np.ndarray         # (n_domains,) region code index of the domain
    min_size: int
    n_border: np.ndarray | None = None  # (n_domains,) border solute atoms added (include_border)

    @property
    def big(self) -> np.ndarray:
        return self.size >= self.min_size


def ordered_solute(solute, counts, max_solute_1st: int = 0, min_solute_2nd: int = 4) -> np.ndarray:
    """``counts`` = order.neighbour_counts(...) = (n1, s1, n2, s2)."""
    n1, s1, n2, s2 = counts
    return np.asarray(solute, bool) & (s1 <= max_solute_1st) & (s2 >= min_solute_2nd)


def find_domains(solute, counts, shells, codes, a_cell: float, max_solute_1st: int = 0, min_solute_2nd: int = 4,
                 min_domain_solute: int = 4, volume_per_solute_a3: float = 1.0, include_border: bool = False) -> Domains:
    ordered = ordered_solute(solute, counts, max_solute_1st, min_solute_2nd)
    n = len(ordered)
    idx = np.nonzero(ordered)[0]
    label = np.full(n, -1, dtype=np.int64)
    if len(idx) == 0:
        e = np.empty(0)
        return Domains(ordered, label, e.astype(np.int64), e, e.astype(np.int64), int(min_domain_solute))
    m = (shells.shell == 2) & ordered[shells.i] & ordered[shells.j]
    local = np.full(n, -1, dtype=np.int64)
    local[idx] = np.arange(len(idx))
    g = coo_matrix((np.ones(int(m.sum()), np.int8), (local[shells.i[m]], local[shells.j[m]])), shape=(len(idx),) * 2)
    nd, lab = connected_components(g, directed=False)
    label[idx] = lab
    size = np.bincount(lab, minlength=nd)
    n_border = np.zeros(nd, dtype=np.int64)
    if include_border:
        sol = np.asarray(solute, bool)
        m2 = shells.shell == 2
        i, j = shells.i[m2], shells.j[m2]
        cand = []
        for a, b in ((i, j), (j, i)):       # a ordered (domain member), b solute, not ordered
            k = ordered[a] & sol[b] & ~ordered[b]
            cand.append(np.column_stack([b[k], label[a[k]]]))
        cand = np.vstack(cand) if cand else np.empty((0, 2), np.int64)
        if len(cand):
            # each border atom joins the largest domain it touches
            order_ = np.lexsort((-size[cand[:, 1]], cand[:, 0]))
            cand = cand[order_]
            first = np.concatenate([[True], cand[1:, 0] != cand[:-1, 0]])
            n_border = np.bincount(cand[first, 1], minlength=nd)
    diam = (6.0 * (size + n_border) * volume_per_solute_a3 * a_cell ** 3 / np.pi) ** (1.0 / 3.0)
    codes = np.asarray(codes)
    nreg = len(REGION_CODES)
    hist = np.zeros((nd, nreg), dtype=np.int64)
    np.add.at(hist, (lab, codes[idx]), 1)
    region = np.argmax(hist, axis=1)  # ties -> lowest index = innermost code after X (X only if all atoms are X)
    return Domains(ordered, label, size, diam, region, int(min_domain_solute), n_border)


def records(dom: Domains, solute, codes, add) -> None:
    """Stage-2 records per region (not X) and for the whole interior ("all")."""
    solute = np.asarray(solute, bool)
    codes = np.asarray(codes)
    big = dom.big
    for k, code in enumerate(REGION_CODES):
        if code == "X":
            continue
        ns = int((solute & (codes == k)).sum())
        add("domain_frac_ordered", (dom.ordered & (codes == k)).sum() / ns if ns else float("nan"), region=code)
        sel = big & (dom.region == k)
        add("domain_n", int(sel.sum()), region=code)
        add("domain_d_mean", float(dom.diameter[sel].mean()) if sel.any() else float("nan"), region=code)
        add("domain_d_max", float(dom.diameter[sel].max()) if sel.any() else float("nan"), region=code)
        add("domain_n_small", int((~big & (dom.region == k)).sum()), region=code)
    inner = dom.region != REGION_CODES.index("X")
    add("domain_n_all", int((big & inner).sum()), region="all")
    add("domain_d_all_mean", float(dom.diameter[big & inner].mean()) if (big & inner).any() else float("nan"),
        region="all")
