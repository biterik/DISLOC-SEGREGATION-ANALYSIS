"""Synthetic tests of disegreg.geometry (numpy/scipy only)."""

import numpy as np
import pytest

from disegreg.geometry import (Cell, LineDistance, densify_polyline, line_coverage, line_defect_axes,
                               line_distances, neighbor_pairs, polyline_length, slab_coverage)

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
    base = cell.to_cartesian([0.0, y / cell.b[1], z / cell.c[2]])  # on the s_a = 0 face, also in a tilted cell
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
    sl = slab_coverage(cell, [segment_along_a(cell, 0.1, 0.5)], 10)
    assert np.allclose(sl, [0, 1, 1, 1, 1, 0, 0, 0, 0, 0])
    sl = slab_coverage(cell, [segment_along_a(cell, 0.85, 1.05)], 10)
    assert np.allclose(sl, [0.5, 0, 0, 0, 0, 0, 0, 0, 0.5, 1])
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


# ---------------------------------------------------------------- fallback axes and the distance used (3.3, 3.4)
A = 3.52
CUT = 0.854 * A


def build(xy=0, ly_cells=12, cores=((1.0, 0.0),), partial_sep=None, drop_slabs=(), blob=True):
    """fcc lattice in a (possibly tilted) cell with synthetic structure labels.

    cores: (y0, amp_y) per dislocation. Without partial_sep: one compact core (radius 3.5 A, 'other') at y0 with an hcp
    ribbon on the glide plane. With partial_sep: two partial cores at y0 and y0 + partial_sep joined by the hcp ribbon
    (the line is still dissociated), axis in the middle. Optionally a small detached non-fcc blob far from the lines and
    slabs where the core atoms are labelled hcp instead of 'other'.
    Returns dict with cell, pos, nonfcc, core, axes (true axes), partials (true partial lines).
    """
    cell = Cell.lammps(10 * A, ly_cells * A, 10 * A, xy=xy * A)  # xy = integer * a keeps the cell vectors lattice vectors
    pos = fcc_in_cell(cell, A)
    z0 = 5 * A + 0.3
    core = np.zeros(len(pos), bool)
    ribbon = np.zeros(len(pos), bool)
    axes, partials = [], []
    for y0, amp in cores:
        if partial_sep is None:
            ax = infinite_line(cell, np.array([0.0, y0, z0]), n=200, amp_y=amp)
            r = brute_force_line_distance(cell, pos, [ax], n_images=1)
            core |= r < 3.5
            ribbon |= (r >= 3.5) & (np.abs(pos[:, 2] - z0) < 1.3) & (r < 12.0)
            axes.append(ax)
        else:
            p1 = infinite_line(cell, np.array([0.0, y0, z0]), n=200, amp_y=amp)
            p2 = infinite_line(cell, np.array([0.0, y0 + partial_sep, z0]), n=200, amp_y=amp)
            mid = infinite_line(cell, np.array([0.0, y0 + partial_sep / 2, z0]), n=200, amp_y=amp)
            rp = brute_force_line_distance(cell, pos, [p1, p2], n_images=1)
            rm = brute_force_line_distance(cell, pos, [mid], n_images=1)
            core |= rp < 3.5
            ribbon |= (rp >= 3.5) & (np.abs(pos[:, 2] - z0) < 1.3) & (rm < partial_sep / 2 + 1.0)
            axes.append(mid)
            partials += [p1, p2]
    s = cell.wrap_fractional(cell.to_fractional(pos))
    nslab = int(round(cell.line_length / 5.0))
    for k in drop_slabs:
        in_slab = np.minimum((s[:, 0] * nslab).astype(int), nslab - 1) == k
        ribbon |= core & in_slab
        core &= ~in_slab
    blob_mask = np.zeros(len(pos), bool)
    if blob:
        blob_mask = np.linalg.norm(pos - cell.to_cartesian([0.5, 0.55, 0.25]), axis=1) < 3.0
    nonfcc = core | ribbon | blob_mask
    return dict(cell=cell, pos=pos, nonfcc=nonfcc, core=core | blob_mask, axes=axes, partials=partials,
                line_core=core, s=s, z0=z0)


@pytest.mark.parametrize("xy", [0, 1])
def test_axis_straight_line_near_y_boundary(xy):
    b = build(xy)
    axes, atom_axis, n_clusters = line_defect_axes(b["cell"], b["pos"], b["nonfcc"], b["core"], CUT)
    assert len(axes) == 1 and n_clusters >= 2  # the blob is a cluster, but does not run along the line
    ax = axes[0]
    assert ax.n_filled_slabs == ax.n_slabs == 7
    pl = ax.polyline(b["cell"])
    assert np.allclose(pl[-1] - pl[0], b["cell"].a)  # closed through the periodic boundary, no winding
    # tolerance = lattice graininess: the core atoms of a 5 A slab are discrete fcc sites inside a 3.5 A cylinder, so their
    # centroid is off the ideal axis by up to ~0.35 A (the exact centroid is tested separately below)
    err = brute_force_line_distance(b["cell"], pl[:-1], b["axes"], n_images=1)
    assert err.max() < 0.4
    r_true = brute_force_line_distance(b["cell"], b["pos"], b["axes"], n_images=1)
    r_axis, _ = LineDistance(b["cell"], [pl]).query(b["pos"])
    assert np.max(np.abs(r_axis - r_true)) < 0.4
    assert np.all(atom_axis[b["line_core"]] == 0) and np.all(atom_axis[~b["nonfcc"]] == -1)


@pytest.mark.parametrize("xy", [0, 1])
def test_axis_equals_independent_slab_centroids(xy):
    """Without smoothing, every axis point is exactly the centroid of the core atoms of its slab (computed here
    independently: minimum-image y around the known line height, plain mean of the fractional z)."""
    b = build(xy)
    cell = b["cell"]
    (ax,), _, _ = line_defect_axes(cell, b["pos"], b["nonfcc"], b["core"], CUT, smoothing_slabs=1)
    s = b["s"][b["line_core"]]
    k = np.minimum((s[:, 0] * ax.n_slabs).astype(int), ax.n_slabs - 1)
    sb = s[:, 1] - np.round(s[:, 1] - 1.0 / cell.b[1])  # unwrap around y0 = 1 A (s_b = y / ly)
    for slab in range(ax.n_slabs):
        m = k == slab
        assert np.isclose(ax.sb[slab] - np.round(ax.sb[slab] - sb[m].mean()), sb[m].mean(), atol=1e-9)
        assert np.isclose(ax.sc[slab], s[m, 2].mean(), atol=1e-9)


def test_axis_wavy_line_and_empty_slabs():
    b = build(1, cores=((1.0, 1.5),), drop_slabs=(2, 3))
    (ax,), _, _ = line_defect_axes(b["cell"], b["pos"], b["nonfcc"], b["core"], CUT)
    assert ax.n_filled_slabs == ax.n_slabs - 2
    err = brute_force_line_distance(b["cell"], ax.polyline(b["cell"])[:-1], b["axes"], n_images=1)
    assert err.max() < 1.0  # smoothing over 3 slabs flattens the 1.5 A wave by ~0.4 A, plus interpolation


def test_axis_none_cases():
    b = build(0)
    assert line_defect_axes(b["cell"], b["pos"], np.zeros(len(b["pos"]), bool), b["core"], CUT)[0] == []
    ribbon_only = b["nonfcc"] & ~b["core"]  # no core-type atoms at all -> no axis
    assert line_defect_axes(b["cell"], b["pos"], ribbon_only, b["core"], CUT)[0] == []


def test_dipole_gives_two_axes():
    b = build(0, ly_cells=16, cores=((5.0, 0.0), (33.0, 0.0)))
    axes, atom_axis, _ = line_defect_axes(b["cell"], b["pos"], b["nonfcc"], b["core"], CUT)
    assert len(axes) == 2
    heights = sorted(float(np.mod(np.mean(ax.polyline(b["cell"])[:, 1]), b["cell"].b[1])) for ax in axes)
    assert np.allclose(heights, [5.0, 33.0], atol=0.4)


def test_line_distances_modes():
    b = build(0)
    cell, pos = b["cell"], b["pos"]
    r_true = brute_force_line_distance(cell, pos, b["axes"], n_images=1)
    # DXA finds the whole line -> polyline mode, r = r_DXA
    full = line_distances(cell, pos, b["axes"], b["nonfcc"], b["core"], A)
    assert full.mode == "polyline" and np.array_equal(full.r, full.r_dxa) and not full.active[0].any()
    assert np.isclose(full.f_cov, 1.0) and np.isclose(full.f_cov_axis[0], 1.0)
    # no DXA line -> the axis everywhere
    none_dxa = line_distances(cell, pos, [], b["nonfcc"], b["core"], A)
    assert none_dxa.mode == "axis" and np.allclose(none_dxa.r, none_dxa.r_axis) and none_dxa.active[0].all()
    assert np.max(np.abs(none_dxa.r - r_true)) < 0.4
    # neither DXA nor a line defect cluster
    nothing = line_distances(cell, pos, [], np.zeros(len(pos), bool), b["core"], A)
    assert nothing.mode == "none" and np.all(np.isinf(nothing.r))


def test_line_distances_partial_dxa_uses_axis_only_where_dxa_misses():
    b = build(0)
    cell, pos = b["cell"], b["pos"]
    found = b["axes"][0][(b["axes"][0][:, 0] <= 0.5 * cell.line_length + 1e-9)]  # DXA finds s_a in [0, 0.5]
    res = line_distances(cell, pos, [found], b["nonfcc"], b["core"], A)
    assert res.mode == "mixed"
    assert res.active[0].tolist() == [False, False, False, False, True, True, True]  # slab 3 is exactly half covered
    r_true = brute_force_line_distance(cell, pos, b["axes"], n_images=1)
    sa = b["s"][:, 0]
    in_found = (sa > 0.1) & (sa < 0.4)
    in_missed = (sa > 0.62) & (sa < 0.95)
    assert np.array_equal(res.r[in_found & (r_true < 5)], res.r_dxa[in_found & (r_true < 5)])
    assert np.max(np.abs(res.r[in_missed] - r_true[in_missed])) < 0.5
    far_from_found = (sa > 0.66) & (sa < 0.84)  # > 5.6 A along X from both ends of the found piece
    assert np.all(res.r_dxa[far_from_found & (r_true < 3)] > 5)  # without the axis these core atoms would not be core


def test_dissociated_line_fault_atoms_stay_out_of_the_core():
    """Proposal (1) of rev. 3: where DXA finds the two partials, the axis (mid-fault) is not used, so the middle of
    the stacking fault is not turned into core."""
    b = build(0, ly_cells=14, cores=((10.0, 0.0),), partial_sep=15.0)
    cell, pos = b["cell"], b["pos"]
    found = [p[p[:, 0] <= 0.6 * cell.line_length + 1e-9] for p in b["partials"]]
    res = line_distances(cell, pos, found, b["nonfcc"], b["core"], A)
    (ax,) = res.axes
    assert np.allclose(np.mean(ax.polyline(cell)[:, 1]), 17.5, atol=0.5)  # axis in the middle of the fault
    mid_fault = (np.abs(pos[:, 1] - 17.5) < 2.0) & (np.abs(pos[:, 2] - b["z0"]) < 1.3)
    covered = (b["s"][:, 0] > 0.05) & (b["s"][:, 0] < 0.5)
    sel = mid_fault & covered  # includes atoms within 2-4 A (along X) of the slabs where DXA misses the line
    assert sel.sum() >= 8
    assert np.all(res.r_axis[sel] < 2.5)  # the old rule (axis everywhere) would make these core atoms
    assert np.all(res.r[sel] > 5.0)       # rev. 3: they keep their distance to the partials


def test_dipole_one_line_found_by_dxa():
    """Proposal (2) of rev. 3: global coverage is 1 because line A is found completely, yet line B is missing; the
    per-cluster assignment activates the axis of line B only."""
    b = build(0, ly_cells=16, cores=((5.0, 0.0), (33.0, 0.0)))
    cell, pos = b["cell"], b["pos"]
    res = line_distances(cell, pos, [b["axes"][0]], b["nonfcc"], b["core"], A)
    assert np.isclose(res.f_cov, 1.0) and res.mode == "mixed"
    yA = [float(np.mean(ax.polyline(cell)[:, 1])) for ax in res.axes]
    iA, iB = (0, 1) if abs(yA[0] - 5.0) < abs(yA[1] - 5.0) else (1, 0)
    assert not res.active[iA].any() and res.active[iB].all()
    assert np.isclose(res.f_cov_axis[iA], 1.0) and res.f_cov_axis[iB] == 0.0
    r_true = brute_force_line_distance(cell, pos, b["axes"], n_images=1)
    assert np.max(np.abs(res.r - r_true)) < 0.5  # lattice graininess of the core centroid of line B (0.3 A in y, 0.44 A in z)
