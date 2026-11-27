"""Along-line core composition and spectrum ([PROPOSAL v0.2] P4): exact synthetic cases."""

import numpy as np
import pytest

from disegreg import alongline
from disegreg.definitions import load_definitions
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


def test_stage2_writes_alongline_with_D2():
    snap, cell, pos, r_true, solute, z0 = synthetic_snapshot()
    d = load_definitions("D2")
    res = analyse_snapshot(snap, d, A, z_gp=z0)
    q = {r["quantity"] for r in res.records}
    assert {"alongline_S", "alongline_dispersion", "alongline_N_bin"} <= q
    # the synthetic core is entirely solute: c = 1, no spectrum
    assert res.value("alongline_c", region="C") == 1.0
    # D2 section 16: a definition in D2 tables, still a proposal in tables of D1 (the v0.2 demo of 24 Nov)
    assert all(quantity_status(x, "D2") == "definition" for x in alongline.QUANTITIES)
    assert all(quantity_status(x, "D1") == "proposal" for x in alongline.QUANTITIES)
    assert quantity_status("Gamma", "D1") == "definition" and "alongline_S" not in PROPOSAL_QUANTITIES
    assert quantity_status("domain_n", "D2") == "proposal"                       # P5 stays a proposal
    n_bins = int(res.value("alongline_n_bins", region="C"))
    assert sum(res.value("alongline_N_bin", region=f"x{k}") for k in range(n_bins)) == res.value("N", region="C")
    plain = analyse_snapshot(snap, load_definitions("D1"), A, z_gp=z0)   # D1: no along-line quantities
    assert not any(r["quantity"].startswith("alongline") for r in plain.records)


def test_window_test_on_the_window_averaged_spectrum():
    """D2 section 16: "periodic" from the p-value of the WINDOW-AVERAGED spectrum (review finding 2), not the mean of the
    per-snapshot p-values. A persistent moderate peak (S = 5 in every one of 10 snapshots) is periodic; per snapshot it is
    not significant."""
    assert alongline.window_peak_p(6.0, 1, 9, 1) == pytest.approx(alongline.binomial_peak_p(6.0, 9, 1), rel=1e-9)
    assert alongline.binomial_peak_p(5.0, 9, 1) > 0.05                    # one snapshot: not significant
    assert alongline.window_peak_p(5.0, 10, 9, 1) < 1e-5                  # window of 10: clearly periodic (2.7e-7)

    def s3(S, n_snap, n_bins=20, L=100.0):
        rows = [{"quantity": "alongline_S", "region": f"q{q}", "mean": v, "n_samples": n_snap, "t_first": 0, "t_last": 1}
                for q, v in zip(range(1, n_bins // 2 + 1), S)]
        return rows + [{"quantity": "alongline_n_bins", "region": "C", "mean": float(n_bins), "n_samples": n_snap},
                       {"quantity": "L_x", "region": "", "mean": L, "n_samples": n_snap}]

    S = np.ones(10)
    S[3] = 5.0                                                           # q = 4: lambda = 25 A
    out = {r["quantity"]: r["mean"] for r in alongline.window_test(s3(S, 10), 0.01, None)}
    assert out["alongline_window_q_peak"] == 4 and out["alongline_window_lambda_peak"] == pytest.approx(25.0)
    assert out["alongline_window_periodic"] == 1.0 and out["alongline_window_p_peak"] < 1e-5
    # binomial noise: false "periodic" in ~1 % of windows (snapshots independent)
    rng = np.random.default_rng(7)
    hits, trials = 0, 4000
    for _ in range(trials):
        n_snap = 5
        per = rng.exponential(1.0, size=(n_snap, 10))
        per[:, -1] = rng.chisquare(1, size=n_snap)                       # q = n/2 = 10: one degree of freedom
        res = {r["quantity"]: r["mean"] for r in alongline.window_test(s3(per.mean(axis=0), n_snap), 0.01, None)}
        hits += res["alongline_window_periodic"] == 1.0
    assert 0.005 < hits / trials < 0.016
    assert alongline.window_test([], 0.01, None) == []


def test_window_test_effective_number_of_snapshots():
    """D2 rev. 3 (Erik 26 Nov, N3 option b): the window test uses n_eff = n_s (SEM_naive / SEM_block)^2 from the block
    averaging of c_C over the same snapshots; without a plateau n_eff = 1, flagged."""
    from disegreg.timeseries import block_sem

    def s3(S, n_snap, n_bins=20, L=100.0):
        rows = [{"quantity": "alongline_S", "region": f"q{q}", "mean": v, "n_samples": n_snap, "t_first": 0, "t_last": 1}
                for q, v in zip(range(1, n_bins // 2 + 1), S)]
        return rows + [{"quantity": "alongline_n_bins", "region": "C", "mean": float(n_bins), "n_samples": n_snap},
                       {"quantity": "L_x", "region": "", "mean": L, "n_samples": n_snap}]

    S = np.ones(10)
    S[3] = 2.2
    rng = np.random.default_rng(3)
    # uncorrelated c_C (plateau at level 0): n_eff = n_s, the same p as for independent snapshots
    c_iid = 0.05 + 0.002 * rng.standard_normal(64)
    b = block_sem(c_iid)
    assert b.plateau and b.level == 0
    info = {}
    rows = {r["quantity"]: r for r in alongline.window_test(s3(S, 64), 0.01, c_iid, info=info)}
    assert rows["alongline_window_n_eff"]["mean"] == pytest.approx(64.0)
    assert rows["alongline_window_p_peak"]["mean"] == pytest.approx(alongline.window_peak_p(2.2, 64, 9, 1), rel=1e-12)
    assert rows["alongline_window_p_peak"]["sem_plateau"] is True and info["alongline_window"]["flag"] == ""
    # strongly correlated c_C (an AR(1) series, plateau reached): n_eff < n_s -> larger p, verdict can flip
    x = np.empty(256)
    x[0] = 0.0
    for i in range(1, 256):
        x[i] = 0.9 * x[i - 1] + rng.standard_normal()
    c_cor = 0.05 + 0.001 * x
    ne = alongline.window_n_eff(c_cor, 256)
    bb = block_sem(c_cor)
    assert bb.plateau and ne["n_eff"] == pytest.approx(256 * (bb.sem_naive / bb.sem) ** 2) and ne["n_eff"] < 256 / 4
    p_ind, p_eff = (alongline.window_peak_p(1.5, n, 9, 1) for n in (256, ne["n_eff"]))
    assert p_eff > 10 * p_ind
    # no plateau (2 or 3 snapshots always; here 3): n_eff = 1 and a flag
    info = {}
    rows = {r["quantity"]: r for r in alongline.window_test(s3(S, 3), 0.01, np.array([0.05, 0.051, 0.052]), info=info)}
    assert rows["alongline_window_n_eff"]["mean"] == 1.0 and rows["alongline_window_n_eff"]["sem_plateau"] is False
    assert rows["alongline_window_p_peak"]["mean"] == pytest.approx(alongline.binomial_peak_p(2.2, 9, 1), rel=1e-9)
    assert "without plateau" in info["alongline_window"]["flag"]
    # one snapshot: n_eff = 1 exactly, no flag
    one = alongline.window_n_eff(np.array([0.05]), 1)
    assert one["n_eff"] == 1.0 and one["flag"] == ""
    # n_eff stays within [1, n_s]
    for c in (c_iid, c_cor):
        assert 1.0 <= alongline.window_n_eff(c, len(c))["n_eff"] <= len(c)

