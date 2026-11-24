"""0 K site-energy tables for the McLean prediction (docs/DEFINITIONS.md section 9, [PROPOSAL v0.2] P1). Pure numpy.

A site table lists, for the atoms of one line period of a relaxed 0 K cell (pure host, one dislocation), the segregation
energy E_i of ONE solute atom put on site i (negative = attraction). The projects' tables are read-only input; one adapter
per table format turns them into (atom id, E_i). Formats:

    segE0K-by-site   Ni-Cu thread 09 and Nb-Au thread 06 (same writer): "# char potential id kind dy_centre(A) dz(A)
                     r_or_dy_from_partial(A) Eseg_relaxed(eV) Eseg_rigid(eV) n_fire"; E_i = Eseg_relaxed
    csv              any table with a header line and the columns id, E_seg (eV), optional character, potential, kind

Core sites [PROPOSAL v0.2 P1]: the site's atom (by id) in the relaxed 0 K cell is closer than core_site_radius (5 A, the
D1 core radius) to a DXA line of that cell (Stage 1 run on the 0 K cell, same DXA parameters, all periodic images).

McLean (independent sites, DEFINITIONS 9): theta_i / (1 - theta_i) = x / (1 - x) exp(-E_i / kT), predicted c_C = mean theta_i
over the core sites. Its local uptake exponent d ln c_C / d ln x = <theta (1 - theta)> / ((1 - x) <theta>) <= 1.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .geometry import LineDistance
from .provenance import sha256_file
from .regions import excluded_mask


@dataclass
class SiteTable:
    source: str
    format: str
    sha256: str
    potential: str
    character: str
    ids: np.ndarray          # atom ids in the 0 K cell
    energy: np.ndarray       # E_i (eV), negative = attraction
    kind: np.ndarray         # site-set label of the table ("core", "core-p1", ...), "" if none

    def __len__(self):
        return len(self.ids)

    def degenerate_sites(self, tol: float = 1e-4, min_group: int = 3) -> int:
        """Sites whose E_i equals that of at least min_group - 1 other sites (to tol eV). In the projects' tables these are
        sites where the partial glided to the solute during the relaxation (Ni-Cu thread 09, 2026-11-04): their E_i is
        the deepest state, not their own."""
        e = np.sort(self.energy)
        if len(e) == 0:
            return 0
        breaks = np.nonzero(np.diff(e) > tol)[0] + 1
        sizes = np.diff(np.concatenate([[0], breaks, [len(e)]]))
        return int(sizes[sizes >= min_group].sum())


def _read_segE0K(path, potential, character):
    ids, energy, kind = [], [], []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            if p[0] != character or p[1] != potential:
                continue
            ids.append(int(p[2]))
            kind.append(p[3])
            energy.append(float(p[7]))
    return ids, energy, kind


def _read_csv(path, potential, character):
    ids, energy, kind = [], [], []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            if r.get("potential") not in (None, "", potential) or r.get("character") not in (None, "", character):
                continue
            ids.append(int(r["id"]))
            energy.append(float(r["E_seg"]))
            kind.append(r.get("kind", "") or "")
    return ids, energy, kind


FORMATS = {"segE0K-by-site": _read_segE0K, "csv": _read_csv}


def read_site_table(path, fmt: str, potential: str, character: str) -> SiteTable:
    if fmt not in FORMATS:
        raise ValueError(f"unknown site-table format {fmt!r} (known: {sorted(FORMATS)})")
    path = Path(path).expanduser()
    ids, energy, kind = FORMATS[fmt](path, potential, character)
    if not ids:
        raise ValueError(f"{path}: no sites for potential {potential!r}, character {character!r}")
    ids = np.asarray(ids, dtype=np.int64)
    u, n = np.unique(ids, return_counts=True)
    if (n > 1).any():
        raise ValueError(f"{path}: atom ids listed more than once for {potential} {character}: {u[n > 1][:5].tolist()}")
    return SiteTable(source=str(path), format=fmt, sha256=sha256_file(path), potential=potential, character=character,
                     ids=ids, energy=np.asarray(energy, dtype=float), kind=np.asarray(kind, dtype=str))


def lattice_period(miller_x, a: float) -> float:
    """Shortest fcc lattice translation along the integer direction [h k l] (made coprime): a/2 |hkl| if h + k + l is
    even, else a |hkl|. [0 -1 1] -> a/sqrt(2), [1 -2 1] -> a sqrt(6)/2."""
    v = np.asarray(miller_x, dtype=np.int64)
    g = np.gcd.reduce(np.abs(v[v != 0]))
    v = v // g
    n = float(np.linalg.norm(v))
    return a * n / 2.0 if v.sum() % 2 == 0 else a * n


@dataclass
class CoreSites:
    r: np.ndarray                 # distance of every tabulated site to the 0 K DXA lines (A)
    is_core: np.ndarray           # r < radius
    radius: float
    n_periods: int                # line periods in the 0 K cell along X
    n_core_cell: int              # interior atoms of the 0 K cell with r < radius
    n_core_expected: float        # per line period: n_core_cell / n_periods
    x_extent: float               # extent of the tabulated sites along X (A), should be < one period
    duplicate: np.ndarray | None = None  # site is a periodic copy (one line period along X) of an earlier table entry
    info: dict = field(default_factory=dict)

    @property
    def n_core_tabulated(self) -> int:
        return int(self.is_core.sum())

    @property
    def coverage(self) -> float:
        """Tabulated core sites / core atoms per period of the 0 K cell (1 = complete table)."""
        return self.n_core_tabulated / self.n_core_expected if self.n_core_expected > 0 else float("nan")


def core_sites(table: SiteTable, snap0, radius: float, miller_x, a_cell: float, spacing: float = 0.25,
               n_images: int = 1, surface_exclusion: float = 15.0) -> CoreSites:
    """Which tabulated sites are core sites of the relaxed 0 K cell ``snap0`` (a Stage-1 Snapshot) [PROPOSAL v0.2 P1]."""
    cell = snap0.cell()
    ids0 = snap0.ids if snap0.ids is not None else np.arange(1, snap0.n_atoms + 1)
    k = np.searchsorted(ids0, table.ids)
    bad = (k >= len(ids0)) | (ids0[np.minimum(k, len(ids0) - 1)] != table.ids)
    if bad.any():
        raise ValueError(f"{table.source}: {int(bad.sum())} site ids not in the 0 K cell {snap0.meta.get('source')}")
    lines = snap0.polylines()
    if not lines:
        raise ValueError(f"no DXA line in the 0 K cell {snap0.meta.get('source')}")
    ld = LineDistance(cell, lines, spacing=spacing, n_images=n_images)
    pos = snap0.positions.astype(float)
    r_sites = ld.query(pos[k])[0]
    excl, _, _ = excluded_mask(pos, surface_exclusion)
    r_all = ld.query(pos[~excl])[0]
    L = cell.line_length
    n_per = max(1, int(round(L / lattice_period(miller_x, a_cell))))
    n_core = int((r_all < radius).sum())
    s = cell.wrap_fractional(cell.to_fractional(pos[k]))[:, 0]
    # extent along X on the circle (the sites may straddle the periodic boundary)
    ss = np.sort(s)
    gaps = np.diff(np.concatenate([ss, ss[:1] + 1.0]))
    extent = (1.0 - gaps.max()) * L if len(ss) > 1 else 0.0
    dup, dE = periodic_copies(cell, pos[k], table.energy, cell.a / n_per)
    return CoreSites(r=r_sites, is_core=(r_sites < radius) & ~dup, radius=float(radius), n_periods=n_per,
                     n_core_cell=n_core, n_core_expected=n_core / n_per, x_extent=float(extent), duplicate=dup,
                     info={"line_length": L, "period": L / n_per, "dxa_segments": len(lines),
                           "n_periodic_copies": int(dup.sum()), "max_dE_periodic_copies": dE})


def periodic_copies(cell, positions, energy, period_vector, tol: float = 0.5):
    """Sites that repeat an EARLIER table entry shifted by a whole number of line periods along the line (a table meant
    for one period that reaches slightly beyond it lists those sites twice). Returns (mask, max |E difference| of the
    copies, NaN if none). Code decision v0.2: a copy is the same site, it is counted once."""
    from scipy.spatial import cKDTree

    p = cell.wrap(np.asarray(positions, float))
    t = np.asarray(period_vector, float)
    n = len(p)
    dup = np.zeros(n, bool)
    dE = []
    if n < 2:
        return dup, float("nan")
    tree = cKDTree(p)
    span = int(np.ceil(cell.line_length / max(np.linalg.norm(t), 1e-9)))
    for j in [m for m in range(-min(span, 3), min(span, 3) + 1) if m]:
        for i, hits in enumerate(tree.query_ball_point(cell.wrap(p + j * t), tol)):
            for h in hits:
                if h < i and not dup[i]:
                    dup[i] = True
                    dE.append(abs(float(energy[i]) - float(energy[h])))
    return dup, (max(dE) if dE else float("nan"))


def mclean_theta(energy, x: float, kT: float) -> np.ndarray:
    """Independent-site occupancy: theta / (1 - theta) = x / (1 - x) exp(-E / kT), computed without overflow."""
    e = np.asarray(energy, dtype=float)
    if not 0.0 < x < 1.0:
        return np.full(e.shape, float(x >= 1.0))
    # theta = 1 / (1 + (1 - x) / x exp(E / kT)); log form avoids overflow for deep sites
    z = np.log1p(-x) - np.log(x) + e / kT
    return 0.5 * (1.0 - np.tanh(0.5 * z))


@dataclass
class McLean:
    c_C: float            # predicted core concentration
    dc_dx: float          # derivative for the error propagation from x
    n_up: float           # local uptake exponent d ln c_C / d ln x (<= 1)
    n_sites: int          # sites averaged over
    theta_max: float


def mclean_prediction(energy_core, x: float, kT: float, n_missing: float = 0.0) -> McLean:
    """Predicted c_C over the core sites; ``n_missing`` core sites without a table entry count with E = 0 (theta = x),
    the option 'bulk' of [PROPOSAL v0.2 P1] (0 for 'tabulated_only')."""
    th = mclean_theta(energy_core, x, kT)
    n = len(th) + n_missing
    if n <= 0:
        return McLean(c_C=float("nan"), dc_dx=float("nan"), n_up=float("nan"), n_sites=0, theta_max=float("nan"))
    c = (th.sum() + n_missing * x) / n
    # d theta / d x = theta (1 - theta) / (x (1 - x))
    var = (np.sum(th * (1 - th)) + n_missing * x * (1 - x)) / n
    dcdx = var / (x * (1 - x)) if 0 < x < 1 else float("nan")
    n_up = dcdx * x / c if c > 0 else float("nan")
    return McLean(c_C=float(c), dc_dx=float(dcdx), n_up=float(n_up), n_sites=int(round(n)),
                  theta_max=float(th.max()) if len(th) else float("nan"))
