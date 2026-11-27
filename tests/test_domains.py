"""Ordered-domain sizes ([PROPOSAL v0.2] P5): an L1_2 Ni3Al particle of known size in fcc Ni, and a random solution."""

import numpy as np
import pytest

from disegreg import domains, order
from disegreg.geometry import Cell
from disegreg.regions import REGION_CODES
from helpers import fcc_in_cell

A = 3.52


def l12_particle(radius, n=14, seed=0, x_random=0.0):
    """Cubic fcc cell n^3 unit cells, fully periodic; inside a sphere of ``radius`` around the centre the cube-corner
    sites are Al (L1_2), everything else Ni (or random Al with fraction x_random outside the sphere)."""
    cell = Cell(np.eye(3) * n * A, pbc=(True, True, True))
    pos = fcc_in_cell(cell, A) + 1e-3
    frac = (pos - 1e-3) / A
    corner = np.all(np.abs(frac - np.round(frac)) < 1e-6, axis=1)      # (0, 0, 0) sublattice = L1_2 corner sites
    c = np.full(3, n * A / 2)
    inside = np.linalg.norm(pos - c, axis=1) < radius
    rng = np.random.default_rng(seed)
    sol = (inside & corner) | (~inside & (rng.random(len(pos)) < x_random))
    return cell, pos, sol, inside


def analyse(cell, pos, sol, codes=None, include_border=False):
    sh = order.shell_pairs(cell, pos, A)
    nc = order.neighbour_counts(len(pos), sol, sh)
    codes = np.full(len(pos), REGION_CODES.index("B")) if codes is None else codes
    return domains.find_domains(sol, nc, sh, codes, A, include_border=include_border), nc


@pytest.mark.parametrize("radius, ratio", [(8.0, 0.80), (10.0, 0.84), (13.5, 0.94), (20.0, 0.94)])
def test_l12_sphere_without_border_is_too_small(radius, ratio):
    """The bias of the plain definition (numbers quoted in docs/OPEN-DECISIONS_v0.2.md, P5)."""
    cell, pos, sol, _ = l12_particle(radius, n=16)
    dom, _ = analyse(cell, pos, sol)
    assert dom.diameter.max() / (2 * radius) == pytest.approx(ratio, abs=0.01)


@pytest.mark.parametrize("radius", [8.0, 10.0, 13.5, 18.0])
def test_l12_sphere_with_border_gives_its_diameter(radius):
    cell, pos, sol, _ = l12_particle(radius, n=16)
    dom, _ = analyse(cell, pos, sol, include_border=True)
    assert dom.big.sum() == 1
    assert dom.size.max() + dom.n_border.max() == sol.sum()     # every Al of the particle, nothing else
    assert dom.diameter.max() == pytest.approx(2 * radius, rel=0.06)
    cell, pos, sol, inside = l12_particle(radius, n=16, x_random=0.13, seed=9)   # in a random 13 % matrix
    dom, _ = analyse(cell, pos, sol, include_border=True)
    assert dom.diameter.max() == pytest.approx(2 * radius, rel=0.10)


@pytest.mark.parametrize("radius", [10.0, 14.0, 18.0])
def test_l12_sphere_gives_its_diameter(radius):
    cell, pos, sol, inside = l12_particle(radius)
    dom, nc = analyse(cell, pos, sol)
    assert dom.big.sum() == 1                                   # one particle -> one domain
    n_al = int(sol.sum())
    assert dom.size.max() <= n_al and dom.size.max() >= 0.6 * n_al  # surface Al with < 4 Al 2nd neighbours drop out
    # interior Al (all 6 2nd neighbours Al) are always ordered
    _, s1, _, s2 = nc
    assert dom.ordered[sol & (s2 == 6)].all() and (s1[sol] == 0).all()
    assert dom.diameter[dom.big][0] == pytest.approx(2 * radius, rel=0.20)   # plain definition: up to -16 % (see above)
    d_all = (6 * n_al * A ** 3 / np.pi) ** (1 / 3)              # all Al of the particle: the continuum value 2 R
    assert d_all == pytest.approx(2 * radius, rel=0.08)


def test_random_solution_has_no_domains():
    cell, pos, sol, _ = l12_particle(0.0, x_random=0.13, seed=4)
    dom, _ = analyse(cell, pos, sol)
    frac = dom.ordered.sum() / sol.sum()
    assert frac < 0.02                                          # binomial: ~0.3 % x P(no solute 1st neighbour)
    assert dom.big.sum() == 0


def test_particle_in_random_matrix_and_records():
    cell, pos, sol, inside = l12_particle(14.0, x_random=0.13, seed=5)
    codes = np.where(inside, REGION_CODES.index("C"), REGION_CODES.index("B"))
    dom, _ = analyse(cell, pos, sol, codes)
    assert dom.big.sum() >= 1
    big = np.argmax(dom.size)
    assert REGION_CODES[dom.region[big]] == "C"
    out = {}
    domains.records(dom, sol, codes, lambda q, v, region="": out.__setitem__((q, region), v))
    assert out[("domain_n", "C")] >= 1 and out[("domain_d_mean", "C")] == pytest.approx(28.0, rel=0.25)
    assert out[("domain_frac_ordered", "C")] > 0.6 and out[("domain_frac_ordered", "B")] < 0.05
    assert out[("domain_n_all", "all")] == dom.big.sum()


def test_no_solute():
    cell, pos, sol, _ = l12_particle(0.0)
    dom, _ = analyse(cell, pos, sol)
    assert len(dom.size) == 0 and not dom.ordered.any()
