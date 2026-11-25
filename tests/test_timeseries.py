"""Stage 3 (DEFINITIONS 8): block averaging on correlated series with a known SEM, drift test, averaging window, MC stream."""

import numpy as np
import pytest

from disegreg import timeseries as ts
from disegreg.definitions import load_definitions

P = load_definitions()["timeseries"]


def ar1(n, phi, rng, sigma=1.0):
    e = rng.normal(0, sigma, n)
    x = np.empty(n)
    x[0] = e[0] / np.sqrt(1 - phi ** 2)
    for i in range(1, n):
        x[i] = phi * x[i - 1] + e[i]
    return x


@pytest.mark.parametrize("phi", [0.0, 0.5, 0.9])
def test_block_sem_matches_exact_ar1(phi):
    """AR(1): SEM^2 -> var / n x (1 + phi) / (1 - phi) for long series. Compared over 40 series (the SEM estimate of one
    series scatters by ~10-20 %)."""
    rng = np.random.default_rng(10)
    n = 2 ** 14
    ratios = []
    for _ in range(40):
        x = ar1(n, phi, rng)
        exact = np.sqrt(1 / (1 - phi ** 2) / n * (1 + phi) / (1 - phi))
        b = ts.block_sem(x)
        ratios.append(b.sem / exact)
        assert b.plateau
    assert np.mean(ratios) == pytest.approx(1.0, abs=0.12)
    if phi > 0:
        assert b.sem > 1.5 * b.sem_naive and b.n_eff < n / 2


def test_block_sem_short_and_constant():
    b = ts.block_sem([1.0, 2.0, 3.0])
    assert not b.plateau and b.sem == pytest.approx(np.std([1, 2, 3], ddof=1) / np.sqrt(3))
    assert ts.block_sem(np.ones(64)).sem == 0.0
    assert np.isnan(ts.block_sem([1.0]).sem)


def test_drift_test_pass_rates_on_flat_series():
    """The rev. 3 criterion |change| < 2 SEM rejects truly flat series about half of the time, because SE(change) =
    sqrt(12) SEM: expected pass rate P(|Z| < 2 / sqrt(12)) = 0.44. The rev. 4 criterion |change| < 2 SE(change): 0.95."""
    rng = np.random.default_rng(11)
    t = np.arange(2048.0)
    d1, prop, ratio = [], [], []
    for _ in range(300):
        d = ts.drift_test(t, ar1(2048, 0.8, rng) + 100)
        ratio.append(d.change_se / (np.sqrt(12) * d.sem))
        d1.append(d.flat_sem)
        prop.append(d.flat_change_se)
        assert d.passes() == d.flat_change_se  # rev. 4 default
    assert np.mean(ratio) == pytest.approx(1.0, abs=0.1)  # detrended SEM ~ SEM for a flat series
    assert np.mean(d1) == pytest.approx(0.44, abs=0.09)
    assert np.mean(prop) == pytest.approx(0.95, abs=0.04)


def test_drift_test_detects_drift():
    rng = np.random.default_rng(11)
    t = np.arange(3000.0)
    drifting = ar1(3000, 0.8, rng) + 100 + 0.002 * t  # change across the window 6, SE(change) ~ 0.3
    d = ts.drift_test(t, drifting)
    assert not d.flat_sem and not d.flat_change_se and d.change == pytest.approx(6.0, abs=1.0)
    with pytest.raises(ValueError):
        d.passes("nonsense")


def test_averaging_window_after_transient():
    rng = np.random.default_rng(12)
    t = np.arange(10000.0)
    y = 50 * np.exp(-t / 800.0) + ar1(10000, 0.5, rng, 0.5) + 10  # saturates after ~ 4000
    w = ts.averaging_window(t, y, P["drift_final_fraction"], P["min_window_fraction"], P["drift_factor"],
                            criterion=P["drift_criterion"])
    assert w.converged and 0.30 <= w.fraction <= 0.75
    assert abs(w.test.mean - 10) < 0.3
    y2 = y + 0.003 * t  # never flat
    for crit in ts.CRITERIA:
        w2 = ts.averaging_window(t, y2, criterion=crit)
        assert not w2.converged and w2.fraction == 0.30


def test_read_mc_stream_dedup_and_columns(tmp_path):
    f = tmp_path / "x_mc-stream.dat"
    f.write_text("# step time[ps] dmu[eV] nCu_all nCu_interior cCu_int\n"
                 "1000 1.0 0.88 10 10 0.1\n2000 2.0 0.88 11 11 0.1\n"
                 "# restart\n2000 2.0 0.88 12 12 0.1\n3000 3.0 0.88 nan nan nan\n3000 3.0 0.88 13 13 0.1\n"
                 "garbage line\n4000 4.0 0.88 14 14\n")
    s = ts.read_mc_stream(f)
    assert s.column("step").tolist() == [1000, 2000, 3000]
    assert s.column("nCu_interior").tolist() == [10, 12, 13]  # last row of a duplicated step; nan row and short row dropped
    assert ts.solute_count_column(s.columns) == "nCu_interior"
    assert ts.solute_count_column(["step", "time[ps]", "NZr_all", "NZr_mcgroup"]) == "NZr_mcgroup"
    assert ts.solute_count_column(["step", "nSol_all", "nSol_interior"]) == "nSol_interior"
    with pytest.raises(ValueError):
        ts.solute_count_column(["step", "foo"])
    assert ts.time_column(s.columns) == "time[ps]"


def test_run_stage3_window_from_stream(tmp_path):
    rng = np.random.default_rng(13)
    steps = np.arange(50000, 1000001, 1000)
    n = 300 + 200 * (1 - np.exp(-(steps - 50000) / 0.8e5)) + rng.normal(0, 3, len(steps))
    f = tmp_path / "s.dat"
    f.write_text("# step time[ps] nX_interior\n" + "\n".join(f"{s} {s / 1000:.1f} {v:.0f}" for s, v in zip(steps, n)))
    stream = ts.read_mc_stream(f)
    snap_steps = np.arange(100000, 1000001, 100000)
    rows = []
    for s in snap_steps:
        rows.append({"timestep": s, "quantity": "c", "region": "C", "basis": "", "value": 0.05 + rng.normal(0, 0.002)})
        rows.append({"timestep": s, "quantity": "N_sol_inside", "region": "20", "basis": "r", "value": 100.0})
    info, out = ts.run_stage3(rows, P, stream=stream, md_timestep_ps=0.001, criterion="change_se")
    assert info["converged"] and info["window_source"] == "mc_stream:nX_interior"
    assert info["window_t_end"] == pytest.approx(1000.0) and info["time_unit"] == "ps"
    c = next(r for r in out if r["quantity"] == "c")
    assert c["n_samples"] == sum(s * 0.001 >= info["window_t_start"] for s in snap_steps)
    assert abs(c["mean"] - 0.05) < 0.003
    assert {"flat", "flat_rev3_sem", "ref"} <= set(c)
    info2, out2 = ts.run_stage3(rows, P)  # no stream: key quantity, rev. 4 criterion, the rev. 3 one reported alongside
    assert info2["window_source"].startswith("stage2:") and info2["criterion"] == "change_se"
    assert "alt_sem_converged" in info2
    info3, out3 = ts.run_stage3(rows, P, criterion="sem")
    assert info3["criterion"] == "sem" and "alt_change_se_converged" in info3
    with pytest.raises(ValueError, match="unknown drift criterion"):
        ts.run_stage3(rows, P, criterion="d1")


def test_d1_rev4_drift_parameters():
    assert P["drift_criterion"] == "change_se" and P["drift_factor"] == 2.0
    assert P["drift_criterion"] == ts.DEFAULT_CRITERION


def test_stage3_counts_snapshots_after_the_stream():
    """v0.1.1: a stream that ends before the last snapshot leaves those snapshots out of every window -- they are counted
    (and `disegreg analyse` prints a WARNING); the window means of the snapshots used are unchanged."""
    from pathlib import Path

    from disegreg.timeseries import Stream, run_stage3

    rows = [{"quantity": "c", "region": "C", "basis": "", "ref": "", "timestep": t, "value": 0.1} for t in range(0, 200, 10)]
    st = np.arange(0, 100, 1.0)
    stream = Stream(path=Path("x"), columns=["step", "time", "nX_interior"],
                    data=np.column_stack([st, st, np.full(len(st), 50.0)]))
    params = {"drift_final_fraction": 0.3, "drift_factor": 2.0, "min_window_fraction": 0.3}
    info, out = run_stage3(rows, params, stream=stream)
    assert info["snapshots_total"] == 20 and info["snapshots_after_window_signal"] == 10
    full = Stream(path=Path("x"), columns=["step", "time", "nX_interior"],
                  data=np.column_stack([np.arange(0, 200, 1.0)] * 2 + [np.full(200, 50.0)]))
    assert run_stage3(rows, params, stream=full)[0]["snapshots_after_window_signal"] == 0
