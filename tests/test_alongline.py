"""Along-line core composition and spectrum ([PROPOSAL v0.2] P4): exact synthetic cases."""

import numpy as np
import pytest

from disegreg import alongline
from disegreg.definitions import load_definitions, load_proposals
from disegreg.geometry import Cell
from disegreg.results import PROPOSAL_QUANTITIES, quantity_status
from disegreg.stage2 import analyse_snapshot
from test_manifest_stage2 import A, synthetic_snapshot


def line_atoms(n_per_bin=40, n_bins=20, width=5.0, seed=0):
    cell = Cell(np.diag([n_bins * width, 50.0, 50.0]), pbc=(True, True, False))
    rng = np.random.default_rng(seed)
    x = (np.repeat(np.arange(n_bins), n_per_bin) + rng.random(n_bins * n_per_bin)) * width
    pos = np.column_stack([x, np.full_like(x, 25.0), np.full_like(x, 25.0)])
    return cell, pos


def test_periodic_particle_array_gives_its_wavelength():
    cell, pos = line_atoms()
    L = cell.line_length
    lam = 25.0                                        # particles every 25 A -> q = L / 25 = 4
    sol = (pos[:, 0] % lam) < 7.0                     # 7 A long solute-rich stretch per period
    al = alongline.profile(cell, pos, sol, np.ones(len(pos), bool), 5.0)
    assert al.n_bins == 20 and al.bin_width == pytest.approx(5.0)
    assert al.q_peak == 4 and al.lambda_peak == pytest.approx(lam)
    assert al.p_peak < 1e-6 and al.dispersion > 5
    assert al.c == pytest.approx(sol.mean())


def test_binomial_noise_has_flat_unit_spectrum():
    cell, pos = line_atoms(seed=1)
    rng = np.random.default_rng(2)
    S, D, P = [], [], []
    for _ in range(1500):
        sol = rng.random(len(pos)) < 0.2
        al = alongline.profile(cell, pos, sol, np.ones(len(pos), bool), 5.0)
        S.append(al.S)
        D.append(al.dispersion)
        P.append(al.p_peak)
    S = np.array(S)
    assert S.mean() == pytest.approx(1.0, abs=0.05)               # every q, binomial noise: mean 1
    assert np.abs(S.mean(axis=0) - 1).max() < 0.2
    assert np.mean(D) == pytest.approx(1.0, abs=0.03)
    assert np.mean(np.array(P) < 0.01) == pytest.approx(0.01, abs=0.006)   # the peak test: 1 % false alarms (n even)


def test_even_spread_and_degenerate_cores():
    cell, pos = line_atoms()
    k = (pos[:, 0] // 5.0).astype(int)
    sol = np.zeros(len(pos), bool)
    for b in range(20):                                             # exactly 8 of 40 solute in every bin
        sol[np.nonzero(k == b)[0][:8]] = True
    al = alongline.profile(cell, pos, sol, np.ones(len(pos), bool), 5.0)
    assert al.dispersion == pytest.approx(0.0, abs=1e-12) and np.allclose(al.S, 0.0)
    full = alongline.profile(cell, pos, np.ones(len(pos), bool), np.ones(len(pos), bool), 5.0)
    assert full.c == 1.0 and np.isnan(full.dispersion)              # all solute: nothing to analyse
    empty = alongline.profile(cell, pos, sol, np.zeros(len(pos), bool), 5.0)
    assert np.isnan(empty.c)
    assert alongline.binomial_peak_p(0.0, 10) == 1.0
    assert alongline.binomial_peak_p(np.log(10) + 10, 1) == pytest.approx(0.1 * np.exp(-10), rel=1e-6)
    from math import erfc
    assert alongline.binomial_peak_p(4.0, 0, 1) == pytest.approx(erfc(np.sqrt(2.0)))


def test_stage2_writes_alongline_with_proposal_status():
    snap, cell, pos, r_true, solute, z0 = synthetic_snapshot()
    d = load_definitions()
    res = analyse_snapshot(snap, d, A, z_gp=z0, proposals=load_proposals())
    q = {r["quantity"] for r in res.records}
    assert {"alongline_S", "alongline_dispersion", "alongline_N_bin"} <= q
    # the synthetic core is entirely solute: c = 1, no spectrum
    assert res.value("alongline_c", region="C") == 1.0
    assert all(quantity_status(x) == "proposal" for x in alongline.QUANTITIES)
    assert quantity_status("Gamma") == "definition" and "alongline_S" in PROPOSAL_QUANTITIES
    n_bins = int(res.value("alongline_n_bins", region="C"))
    assert sum(res.value("alongline_N_bin", region=f"x{k}") for k in range(n_bins)) == res.value("N", region="C")
    plain = analyse_snapshot(snap, d, A, z_gp=z0)                  # without proposals: D1 quantities only
    assert not any(r["quantity"].startswith("alongline") for r in plain.records)
