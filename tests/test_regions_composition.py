"""Regions (DEFINITIONS 4), composition counts, Gamma and partial separation (section 6)."""

import numpy as np
import pytest

from disegreg import composition, regions
from disegreg.cache import STRUCTURE_CODES
from disegreg.definitions import load_definitions
from disegreg.geometry import Cell

FCC, HCP, OTHER = STRUCTURE_CODES["fcc"], STRUCTURE_CODES["hcp"], STRUCTURE_CODES["other"]
REG = load_definitions()["regions"]


def test_excluded_layers():
    z = np.array([0.0, 14.9, 15.0, 50.0, 85.0, 85.1, 100.0])
    pos = np.column_stack([np.zeros_like(z), np.zeros_like(z), z])
    m, lo, hi = regions.excluded_mask(pos, 15.0)
    assert (lo, hi) == (0.0, 100.0)
    assert m.tolist() == [True, True, False, False, False, True, True]


def test_assign_regions_every_rule():
    # (r, z - z_gp, a-CNA label, excluded) -> expected code
    cases = [
        (2.0, 0.0, OTHER, True, "X"),     # X beats everything
        (4.9, 0.0, HCP, False, "C"),      # core beats fault
        (5.0, 1.0, HCP, False, "F"),
        (8.0, 3.9, HCP, False, "F"),
        (8.0, 4.1, HCP, False, "S"),      # hcp but off the glide plane
        (8.0, 0.0, FCC, False, "S"),
        (9.99, 0.0, OTHER, False, "S"),
        (10.0, 0.0, FCC, False, "T"),
        (30.0, 2.0, HCP, False, "F"),
        (59.9, 0.0, OTHER, False, "T"),
        (60.0, 10.0, FCC, False, "B"),
        (60.0, 10.0, OTHER, False, "O"),
        (70.0, 10.0, HCP, False, "O"),
        (70.0, 1.0, HCP, False, "F"),     # precedence F > O (stated in regions.py)
    ]
    z_gp = 100.0
    r = np.array([c[0] for c in cases])
    pos = np.column_stack([np.zeros(len(r)), np.zeros(len(r)), z_gp + np.array([c[1] for c in cases])])
    acna = np.array([c[2] for c in cases])
    exc = np.array([c[3] for c in cases])
    codes = regions.assign_regions(pos, acna, r, exc, z_gp, REG)
    assert [regions.REGION_CODES[k] for k in codes] == [c[4] for c in cases]


def test_far_slab40():
    pos = np.array([[0, 0, 100 + 41.0], [0, 0, 100 - 39.0], [0, 0, 100 - 41.0], [0, 0, 100 + 45.0]])
    m = regions.far_slab40_mask(pos, np.array([FCC, FCC, HCP, FCC]), np.array([False, False, False, True]), 100.0, 40.0)
    assert m.tolist() == [True, False, False, False]


def test_glide_plane_check():
    a = 3.52
    gp = regions.glide_plane(100.0, 100.0 + 0.9 * a / np.sqrt(3), "shockley", a, 1.0)
    assert gp.source == "manifest" and not gp.flag
    gp = regions.glide_plane(100.0, 100.0 + 1.1 * a / np.sqrt(3), "shockley", a, 1.0)
    assert gp.flag
    gp = regions.glide_plane(None, 101.0, "axis", a, 1.0)
    assert gp.z_gp == 101.0 and gp.source == "lines" and gp.flag  # no manifest value is always flagged


def test_z_gp_from_lines_is_length_weighted_median():
    cell = Cell.lammps(20.0, 20.0, 50.0)
    short = np.array([[0, 5, 30.0], [1, 5, 30.0]])   # 1 A at z = 30
    long = np.array([[0, 10, 25.0], [20, 10, 25.0]])  # 20 A at z = 25
    z, src = regions.z_gp_from_lines(cell, [short, long])
    assert src == "shockley" and z == pytest.approx(25.0)
    z, src = regions.z_gp_from_lines(cell, [], [long])
    assert src == "axis" and z == pytest.approx(25.0)
    assert regions.z_gp_from_lines(cell, [], [])[1] == "none"


def test_radial_bins():
    edges = REG["radial_bins"]
    assert regions.radial_bin_labels(edges)[0] == "0-5" and regions.radial_bin_labels(edges)[-1] == ">=100"
    idx = regions.radial_bin_index([0.0, 4.99, 5.0, 99.0, 100.0, 1e9], edges)
    assert idx.tolist() == [0, 0, 1, 6, 7, 7]


def test_counts_and_gamma_vs_brute_force():
    rng = np.random.default_rng(1)
    n = 20000
    r = rng.uniform(0, 130, n)
    r_dxa = r + rng.uniform(0, 3, n)
    codes = rng.integers(0, 7, n).astype(np.int8)
    solute = rng.random(n) < np.where(r < 10, 0.3, 0.01)
    edges = REG["radial_bins"]
    radii = [20.0, 40.0, 5.0, 10.0]
    L = 100.0
    cnt = composition.count(solute, codes, r, r_dxa, edges, radii, L)
    interior = codes != 0
    assert cnt.n_region.sum() == n and cnt.nsol_region.sum() == solute.sum()
    assert cnt.radii.tolist() == [5.0, 10.0, 20.0, 40.0]
    x_ref = 0.01
    for basis, d in (("r", r), ("r_dxa", r_dxa)):
        for R, ni, nsi in zip(cnt.radii, cnt.n_inside[basis], cnt.nsol_inside[basis]):
            m = interior & (d < R)
            assert ni == m.sum() and nsi == (m & solute).sum()
            g = composition.excess(nsi, ni, x_ref, L)
            assert g == pytest.approx(((m & solute).sum() - x_ref * m.sum()) / L)
        assert cnt.n_bin[basis].sum() == interior.sum()
        m = interior & (d >= 20) & (d < 40)
        assert cnt.n_bin[basis][3] == m.sum() and cnt.nsol_bin[basis][3] == (m & solute).sum()


def test_region_excesses_add_up():
    rng = np.random.default_rng(2)
    n = 5000
    codes = rng.integers(1, 5, n).astype(np.int8)  # C, F, S, T only
    solute = rng.random(n) < 0.05
    r = rng.uniform(0, 59, n)
    cnt = composition.count(solute, codes, r, r, [0, 5, 10, 20, 40, 60], [60.0], 50.0)
    gk = composition.excess(cnt.nsol_region, cnt.n_region, 0.02, 50.0)
    g = composition.excess(cnt.nsol_inside["r"][0], cnt.n_inside["r"][0], 0.02, 50.0)
    assert gk.sum() == pytest.approx(g)


def test_concentration_nan_for_empty():
    c = composition.concentration([1, 0], [4, 0])
    assert c[0] == 0.25 and np.isnan(c[1])


@pytest.mark.parametrize("xy", [0.0, 1.1])
def test_partial_separation(xy):
    cell = Cell.lammps(100.0, 60.0, 80.0, xy=xy)
    a, b, c = cell.a, cell.b, cell.c

    def line(sb, sa0=0.0, sa1=1.0, n=40):
        sa = np.linspace(sa0, sa1, n)
        return cell.to_cartesian(np.column_stack([sa, np.full(n, sb), np.full(n, 0.5)]))

    d = 18.0
    sb1, sb2 = 0.95, 0.95 + d / 60.0  # the second partial crosses the Y boundary
    p1 = line(sb1)
    p2a, p2b = line(sb2, 0.0, 0.5), line(sb2, 0.5, 1.0)  # split into two DXA segments
    ps = composition.partial_separation(cell, [p1, p2a, p2b], 5.0)
    assert ps.n_bins == 20
    assert ps.n_valid == 20  # the two pieces of partial 2 meet in bins 0 and 10 and are grouped there
    assert ps.mean == pytest.approx(d, abs=1e-6) and ps.std < 1e-6
    ps = composition.partial_separation(cell, [p1], 5.0)
    assert ps.n_valid == 0 and np.isnan(ps.mean)
    ps = composition.partial_separation(cell, [p1, p2a, p2b, line(0.5)], 5.0)  # three partials: no valid bin
    assert ps.n_valid == 0
