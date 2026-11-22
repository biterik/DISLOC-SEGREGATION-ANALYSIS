"""Synthetic tests of disegreg.geometry (numpy/scipy only)."""

import numpy as np
import pytest

from disegreg.geometry import (Cell, LineDistance, densify_polyline, fallback_core_axis, line_coverage,
                               neighbor_pairs, polyline_length)

from helpers import brute_force_line_distance, brute_force_pairs, fcc_in_cell

RNG = np.random.default_rng(20261122)
TOL = 0.02  # DEFINITIONS 3.1: max error < 0.02 A at r >= 0.5 A


def random_positions(cell, n, rng=RNG):
    s = rng.random((n, 3))
    return cell.to_cartesian(s)


def infinite_line(cell, start, n=40, amp_y=0.0, amp_z=0.0):
    """Continuous polyline along a from start to start + a (an 'infinite' DXA line), optionally wavy."""
    t = np.linspace(0.0, 1.0, n + 1)
    pts = start[None, :] + t[:, None] * cell.a[None, :]
    pts[:, 1] += amp_y * np.sin(2 * np.pi * t)
    pts[:, 2] += amp_z * np.cos(2 * np.pi * t)
    return pts


# ---------------------------------------------------------------- cell
def test_fractional_roundtrip_and_wrap():
    cell = Cell.lammps(40.0, 50.0, 60.0, xy=12.0, xz=1.0, yz=-2.0, origin=(-3.0, 5.0, 1.0))
    r = RNG.normal(scale=80.0, size=(200, 3))
    assert np.allclose(cell.to_cartesian(cell.to_fractional(r)), r)
    s = cell.to_fractional(cell.wrap(r))
    assert np.all((s[:, :2] >= 0) & (s[:, :2] < 1))
    assert np.allclose(s[:, 2], cell.to_fractional(r)[:, 2])  # z is not periodic: untouched
    # wrapping moves by lattice vectors only
    dn = cell.to_fractional(cell.wrap(r)) - cell.to_fractional(r)
    assert np.allclose(dn, np.round(dn))


def test_ovito_matrix_and_widths():
    cell = Cell.lammps(40.0, 50.0, 60.0, xy=12.0)
    m34 = np.column_stack([cell.H, [1.0, 2.0, 3.0]])
    c2 = Cell(m34)
    assert np.allclose(c2.origin, [1, 2, 3]) and np.allclose(c2.H, cell.H)
    w = cell.widths()
    assert np.isclose(w[1], 50.0) and np.isclose(w[2], 60.0)
    assert np.isclose(w[0], 40.0 * 50.0 / np.hypot(50.0, 12.0))
    assert len(cell.image_shifts(1)) == 9 and np.allclose(cell.image_shifts(1)[0], 0)
    assert np.isclose(cell.line_length, 40.0)


# ---------------------------------------------------------------- densify
def test_densify_keeps_vertices_and_spacing():
    p = np.array([[0, 0, 0], [1.0, 0, 0], [1.0, 3.3, 0], [1.0, 3.3, 0.1]])
    d = densify_polyline(p, 0.25)
    steps = np.linalg.norm(np.diff(d, axis=0), axis=1)
    assert steps.max() <= 0.25 + 1e-12
    for v in p:
        assert np.min(np.linalg.norm(d - v, axis=1)) < 1e-12
    assert np.isclose(polyline_length(d), polyline_length(p))


# ---------------------------------------------------------------- distances
def test_straight_line_orthogonal_analytic():
    cell = Cell.lammps(30.0, 80.0, 70.0)
    y0, z0 = 3.0, 35.0
    line = infinite_line(cell, np.array([0.0, y0, z0]), n=1)
    pos = random_positions(cell, 3000)
    r, _ = LineDistance(cell, [line]).query(pos)
    dy = pos[:, 1] - y0
    dy -= 80.0 * np.round(dy / 80.0)
    ref = np.hypot(dy, pos[:, 2] - z0)
    far = ref >= 0.5
    assert np.all(r >= ref - 1e-9)  # dense points never underestimate
    assert np.max(np.abs(r - ref)[far]) < TOL


@pytest.mark.parametrize("xy", [1.1, 12.0, -15.0])
def test_triclinic_wavy_line_vs_brute_force(xy):
    """Solutes near the y boundary of a tilted cell get the right distance (exact images incl. the xy tilt)."""
    cell = Cell.lammps(40.0, 50.0, 60.0, xy=xy)
    line = infinite_line(cell, np.array([0.0, 1.5, 30.0]), n=30, amp_y=2.0, amp_z=1.0)
    pos = random_positions(cell, 1500)
    pos[:300, 1] = cell.b[1] - RNG.random(300) * 3.0  # concentrate atoms near the upper y face
    pos[:300, 0] += xy * (pos[:300, 1] / cell.b[1])  # keep them inside the tilted cell
    pos[:300, 2] = 30.0 + RNG.normal(scale=4.0, size=300)
    r, _ = LineDistance(cell, [line]).query(pos)
    ref = brute_force_line_distance(cell, pos, [line])
    far = ref >= 0.5
    assert np.all(r >= ref - 1e-9)
    assert np.max(np.abs(r - ref)[far]) < TOL


def test_orthogonal_box_approximation_fails_for_tilt():
    """Documents the defect of the old Nb-Au approach: ignoring the xy tilt gives wrong distances at the y boundary."""
    xy = 12.0
    cell = Cell.lammps(40.0, 50.0, 60.0, xy=xy)
    line = infinite_line(cell, np.array([0.0, 1.0, 30.0]), n=30, amp_y=3.0)
    pos = cell.to_cartesian(np.column_stack([RNG.random(400), 0.97 + 0.03 * RNG.random(400), np.full(400, 0.5)]))
    ref = brute_force_line_distance(cell, pos, [line])
    ortho = Cell.lammps(40.0, 50.0, 60.0, xy=0.0)
    wrong = brute_force_line_distance(ortho, pos, [line])
    assert np.max(np.abs(wrong - ref)) > 1.0
    r, _ = LineDistance(cell, [line]).query(pos)
    assert np.max(np.abs(r - ref)[ref >= 0.5]) < TOL


def test_line_outside_primary_cell_same_result():
    cell = Cell.lammps(40.0, 50.0, 60.0, xy=8.0)
    line = infinite_line(cell, np.array([0.0, 20.0, 30.0]), n=20, amp_y=4.0)
    shifted = line + 2 * cell.a - cell.b
    pos = random_positions(cell, 500)
    r1, _ = LineDistance(cell, [line]).query(pos)
    r2, _ = LineDistance(cell, [shifted]).query(pos)
    assert np.allclose(r1, r2, atol=1e-9)


def test_two_lines_owner_and_no_line():
    cell = Cell.lammps(40.0, 100.0, 60.0)
    l1 = infinite_line(cell, np.array([0.0, 20.0, 30.0]), n=1)
    l2 = infinite_line(cell, np.array([0.0, 60.0, 30.0]), n=1)
    pos = np.array([[5.0, 22.0, 30.0], [5.0, 58.0, 31.0], [5.0, 95.0, 30.0]])
    r, k = LineDistance(cell, [l1, l2]).query(pos)
    assert k.tolist() == [0, 1, 0]  # the last atom is 25 A from line 0 through the y boundary
    assert np.allclose(r, [2.0, np.hypot(2, 1), 25.0], atol=1e-6)
    r0, k0 = LineDistance(cell, []).query(pos)
    assert np.all(np.isinf(r0)) and np.all(k0 == -1)


# ---------------------------------------------------------------- coverage
def segment_along_a(cell, s_from, s_to, y=10.0, z=30.0):
    base = np.array([0.0, y, z])
    return np.array([base + s_from * cell.a, base + s_to * cell.a])


@pytest.mark.parametrize("xy", [0.0, 7.0])
def test_line_coverage(xy):
    cell = Cell.lammps(40.0, 50.0, 60.0, xy=xy)
    full = infinite_line(cell, np.array([0.0, 10.0, 30.0]), n=5)
    full2 = infinite_line(cell, np.array([0.0, 25.0, 30.0]), n=5)
    assert np.isclose(line_coverage(cell, [full]), 1.0)
    assert np.isclose(line_coverage(cell, [full, full2]), 1.0)  # two partials: projected, not 2
    assert np.isclose(line_coverage(cell, [segment_along_a(cell, 0.1, 0.5)]), 0.4)
    assert np.isclose(line_coverage(cell, [segment_along_a(cell, 0.8, 1.1)]), 0.3)  # crosses the boundary
    assert np.isclose(line_coverage(cell, [segment_along_a(cell, -1.2, -0.9)]), 0.3)  # outside the primary cell
    assert np.isclose(line_coverage(cell, [segment_along_a(cell, 0.0, 0.5), segment_along_a(cell, 0.3, 0.6)]), 0.6)
    assert np.isclose(line_coverage(cell, [segment_along_a(cell, 0.0, 0.2), segment_along_a(cell, 0.5, 0.6)]), 0.3)
    assert np.isclose(line_coverage(cell, [segment_along_a(cell, 0.6, 0.2)]), 0.4)  # direction does not matter
    assert line_coverage(cell, []) == 0.0
    perpendicular = np.array([[10.0, 0.0, 30.0], [10.0, 40.0, 30.0]])
    assert line_coverage(cell, [perpendicular]) == pytest.approx(40.0 * abs(xy) / 50.0 / 40.0, abs=1e-9)


# ---------------------------------------------------------------- neighbours
def test_neighbor_pairs_fcc_coordination():
    a = 3.52
    cell = Cell.lammps(6 * a, 6 * a, 6 * a, pbc=(True, True, True))
    pos = fcc_in_cell(cell, a)
    assert len(pos) == 4 * 6 ** 3
    i, j, d = neighbor_pairs(cell, pos, 0.854 * a)
    coord = np.bincount(np.concatenate([i, j]), minlength=len(pos))
    assert np.all(coord == 12)
    assert np.allclose(d, a / np.sqrt(2))
    i2, j2, d2 = neighbor_pairs(cell, pos, 1.112 * a)
    coord2 = np.bincount(np.concatenate([i2, j2]), minlength=len(pos))
    assert np.all(coord2 == 18)  # 12 first + 6 second neighbours


def test_neighbor_pairs_free_z_surface():
    a = 3.52
    cell = Cell.lammps(6 * a, 6 * a, 6 * a, pbc=(True, True, False))
    pos = fcc_in_cell(cell, a)
    i, j, _ = neighbor_pairs(cell, pos, 0.854 * a)
    coord = np.bincount(np.concatenate([i, j]), minlength=len(pos))
    z = pos[:, 2]
    surface = (z < 0.1) | (z > z.max() - 0.1)
    assert np.all(coord[~surface] == 12)
    assert np.all(coord[surface] == 8)  # (001) surface atoms


@pytest.mark.parametrize("pbc", [(True, True, False), (True, True, True)])
def test_neighbor_pairs_triclinic_vs_brute_force(pbc):
    cell = Cell.lammps(14.0, 16.0, 15.0, xy=5.0, xz=1.0, yz=-2.0, pbc=pbc)
    pos = random_positions(cell, 300) + cell.a * 3  # also: positions not wrapped
    i, j, d = neighbor_pairs(cell, pos, 3.0)
    ref, dmin = brute_force_pairs(cell, pos, 3.0)
    assert set(zip(i.tolist(), j.tolist())) == ref
    assert np.allclose(d, dmin[i, j])


def test_neighbor_pairs_cutoff_too_large():
    cell = Cell.lammps(10.0, 10.0, 10.0)
    with pytest.raises(ValueError):
        neighbor_pairs(cell, np.zeros((2, 3)), 6.0)


# ---------------------------------------------------------------- fallback axis
def _fcc_with_core(xy, amp_y=0.0, y0=1.0, drop_slabs=()):
    """fcc lattice in a (possibly tilted) cell; a wavy core cylinder labelled 'other', an hcp ribbon, and a smaller
    detached non-fcc cluster that must be ignored. Returns cell, positions, nonfcc, core, true axis polyline."""
    a = 3.52
    cell = Cell.lammps(10 * a, 12 * a, 10 * a, xy=xy * a)  # xy = integer * a keeps the cell vectors lattice vectors
    pos = fcc_in_cell(cell, a)
    z0 = 5 * a + 0.3
    axis = infinite_line(cell, np.array([0.0, y0, z0]), n=200, amp_y=amp_y)
    r = brute_force_line_distance(cell, pos, [axis], n_images=1)
    core = r < 3.5
    # hcp ribbon: same glide plane, 3.5 - 12 A away in y (on one side), not core
    s = cell.wrap_fractional(cell.to_fractional(pos))
    ribbon = (~core) & (np.abs(pos[:, 2] - z0) < 1.3) & (r < 12.0)
    # detached small cluster far from the line
    blob_c = cell.to_cartesian([0.5, 0.55, 0.25])
    blob = np.linalg.norm(pos - blob_c, axis=1) < 3.0
    nslab = int(round(cell.line_length / 5.0))
    for k in drop_slabs:  # slabs where DXA-like 'other' labels are missing (core atoms labelled hcp instead)
        in_slab = np.minimum((s[:, 0] * nslab).astype(int), nslab - 1) == k
        core = core & ~in_slab
        ribbon = ribbon | (in_slab & (r < 3.5))
    nonfcc = core | ribbon | blob
    core_label = core | blob  # the blob is 'other' too: it must be excluded because it is not in the largest cluster
    return cell, pos, nonfcc, core_label, axis, a


@pytest.mark.parametrize("xy", [0, 1])
def test_fallback_axis_straight_line_near_y_boundary(xy):
    cell, pos, nonfcc, core, axis, a = _fcc_with_core(xy, amp_y=0.0, y0=1.0)
    res = fallback_core_axis(cell, pos, nonfcc, core, 0.854 * a)
    assert res is not None and res.n_clusters >= 2
    assert res.n_filled_slabs == res.n_slabs == 7
    assert np.allclose(res.points[-1] - res.points[0], cell.a)  # closed through the periodic boundary, no winding
    # tolerance = lattice graininess: the 'other' atoms of a 5 A slab are discrete fcc sites inside a 3.5 A cylinder,
    # so their centroid is off the ideal axis by up to ~0.35 A (the exact centroid is tested separately below)
    err = brute_force_line_distance(cell, res.points[:-1], [axis], n_images=1)
    assert err.max() < 0.4
    # distances of the atoms to the recovered axis vs to the true axis
    r_true = brute_force_line_distance(cell, pos, [axis], n_images=1)
    r_axis, _ = LineDistance(cell, [res.points]).query(pos)
    assert np.max(np.abs(r_axis - r_true)) < 0.4


@pytest.mark.parametrize("xy", [0, 1])
def test_fallback_axis_equals_independent_slab_centroids(xy):
    """Without smoothing, every axis point is exactly the centroid of the core atoms of its slab (computed here
    independently: minimum-image y around the known line height, plain mean of z)."""
    cell, pos, nonfcc, core, axis, a = _fcc_with_core(xy, amp_y=0.0, y0=1.0)
    res = fallback_core_axis(cell, pos, nonfcc, core, 0.854 * a, smoothing_slabs=1)
    blob_c = cell.to_cartesian([0.5, 0.55, 0.25])
    sel = core & (np.linalg.norm(pos - blob_c, axis=1) >= 3.0)  # core atoms of the line only
    s = cell.wrap_fractional(cell.to_fractional(pos[sel]))
    k = np.minimum((s[:, 0] * res.n_slabs).astype(int), res.n_slabs - 1)
    sb = s[:, 1] - np.round(s[:, 1] - 1.0 / cell.b[1])  # unwrap around y0 = 1 A (s_b = y / ly)
    for slab in range(res.n_slabs):
        m = k == slab
        expect_sb, expect_sc = sb[m].mean(), s[m, 2].mean()
        got = cell.to_fractional(res.points[slab])
        assert np.isclose(got[0], (slab + 0.5) / res.n_slabs)
        assert np.isclose(got[1] - np.round(got[1] - expect_sb), expect_sb, atol=1e-9)
        assert np.isclose(got[2], expect_sc, atol=1e-9)


def test_fallback_axis_wavy_line_and_empty_slabs():
    cell, pos, nonfcc, core, axis, a = _fcc_with_core(1, amp_y=1.5, y0=1.0, drop_slabs=(2, 3))
    res = fallback_core_axis(cell, pos, nonfcc, core, 0.854 * a)
    assert res.n_filled_slabs == res.n_slabs - 2
    err = brute_force_line_distance(cell, res.points[:-1], [axis], n_images=1)
    assert err.max() < 1.0  # smoothing over 3 slabs flattens the 1.5 A wave by ~0.4 A, plus interpolation


def test_fallback_axis_none_cases():
    cell, pos, nonfcc, core, axis, a = _fcc_with_core(0)
    assert fallback_core_axis(cell, pos, np.zeros(len(pos), bool), core, 0.854 * a) is None
    ribbon_only = nonfcc & ~core
    blob_far = np.zeros(len(pos), bool)
    assert fallback_core_axis(cell, pos, ribbon_only, blob_far, 0.854 * a) is None
