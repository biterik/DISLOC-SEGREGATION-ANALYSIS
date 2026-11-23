"""Chemical order (DEFINITIONS 7): ordered reference structures with exact alpha and R2, random solid solution, clusters."""

import numpy as np
import pytest

from disegreg import order
from disegreg.geometry import Cell
from disegreg.regions import REGION_CODES
from helpers import fcc_in_cell

A = 3.52
B = REGION_CODES.index("B")


def periodic_fcc(n=(4, 4, 4)):
    cell = Cell(np.diag(np.array(n, float) * A), pbc=(True, True, True))
    pos = fcc_in_cell(cell, A)
    return cell, pos, np.round(pos / A * 2).astype(int)  # positions in units of a/2 (integers)


def analyse(cell, pos, solute, x_ref=None):
    sh = order.shell_pairs(cell, pos, A)
    nc = order.neighbour_counts(len(pos), solute, sh)
    codes = np.full(len(pos), B, dtype=np.int8)
    return order.order_by_region(solute, codes, nc, x_ref), sh


def test_shell_coordination_12_and_6():
    cell, pos, _ = periodic_fcc()
    sh = order.shell_pairs(cell, pos, A)
    n1, _, n2, _ = order.neighbour_counts(len(pos), np.zeros(len(pos), bool), sh)
    assert np.all(n1 == 12) and np.all(n2 == 6)


def test_L12():
    cell, pos, h = periodic_fcc()
    solute = np.all(h % 2 == 0, axis=1)  # cube corners
    assert solute.mean() == pytest.approx(0.25)
    ob, _ = analyse(cell, pos, solute, x_ref=0.25)
    assert ob.p1[B] == 0.0 and ob.p2[B] == 1.0
    assert ob.alpha1_local[B] == pytest.approx(-1 / 3) and ob.alpha2_local[B] == pytest.approx(1.0)
    assert ob.alpha1_xref[B] == pytest.approx(-1 / 3)
    assert ob.R2[B] == pytest.approx(4.0)


def test_D022():
    cell, pos, h = periodic_fcc((4, 4, 4))  # 4 cubes in z = 2 tetragonal cells
    z_even = (h[:, 2] % 4) == 0
    z_odd = (h[:, 2] % 4) == 2
    corner_xy = (h[:, 0] % 2 == 0) & (h[:, 1] % 2 == 0)
    centre_xy = (h[:, 0] % 2 == 1) & (h[:, 1] % 2 == 1)
    solute = (z_even & corner_xy) | (z_odd & centre_xy)  # body-centred tetragonal solute lattice (a, a, 2a)
    assert solute.mean() == pytest.approx(0.25)
    ob, _ = analyse(cell, pos, solute)
    assert ob.p1[B] == 0.0
    assert ob.p2[B] == pytest.approx(4 / 6)
    assert ob.R2[B] == pytest.approx(8 / 3)


def test_random_solution_alpha_zero():
    cell, pos, _ = periodic_fcc((10, 10, 10))
    rng = np.random.default_rng(3)
    solute = rng.random(len(pos)) < 0.1
    ob, _ = analyse(cell, pos, solute, x_ref=0.1)
    # ~400 solute atoms x 12 neighbours: statistical error of alpha1 ~ 0.015
    assert abs(ob.alpha1_local[B]) < 0.05 and abs(ob.alpha2_local[B]) < 0.07
    assert ob.R2[B] == pytest.approx(1.0, abs=0.15)
    assert ob.n_solute[0] == 0 and np.isnan(ob.p1[0])  # empty region -> NaN


def test_clusters():
    cell, pos, h = periodic_fcc()

    def at(*hh):
        return int(np.nonzero(np.all(h == np.array(hh), axis=1))[0][0])

    solute = np.zeros(len(pos), bool)
    dimer = [at(0, 0, 0), at(1, 1, 0)]
    chain = [at(4, 4, 4), at(5, 5, 4), at(6, 6, 4)]    # 1st-neighbour chain of three
    second = [at(0, 0, 4), at(2, 0, 4)]                 # 2nd neighbours only: two monomers
    across = [at(7, 4, 1), at(0, 5, 1)]                 # 1st neighbours across the periodic x boundary
    for idx in dimer + chain + second + across:
        solute[idx] = True
    sh = order.shell_pairs(cell, pos, A)
    size = order.solute_cluster_sizes(len(pos), solute, sh)
    assert size[dimer].tolist() == [2, 2] and size[chain].tolist() == [3, 3, 3]
    assert size[second].tolist() == [1, 1] and size[across].tolist() == [2, 2]
    assert size[~solute].max() == 0
    codes = np.full(len(pos), B, dtype=np.int8)
    cb = order.clusters_by_region(size, solute, codes)
    assert cb.n_solute[B] == 9 and cb.max_size[B] == 3
    assert cb.mean_size[B] == pytest.approx((2 * 2 + 3 * 3 + 2 * 1 + 2 * 2) / 9)
    assert cb.frac_isolated[B] == pytest.approx(2 / 9)
    assert cb.histogram[B, :3].tolist() == [2, 4, 3]  # monomers, dimers, trimers (solute atoms)
