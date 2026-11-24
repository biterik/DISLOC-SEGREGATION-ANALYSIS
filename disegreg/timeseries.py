"""Stage 3: time series of one run -- drift test, block averaging, averaging window (docs/DEFINITIONS.md section 8).

    block averaging   Flyvbjerg & Petersen, J. Chem. Phys. 91, 461 (1989): neighbouring samples are averaged pairwise again
                      and again; the standard error of the mean (SEM) estimated from the block means grows while the blocks
                      are still correlated and reaches a plateau once they are not. The plateau value is the SEM.
    drift test        over a final window: linear trend x window length (the change across the window) compared with
                      its own standard error, SE(change) = SEM_res x L x sqrt(n / S_tt) (= sqrt(12) SEM for evenly
                      spaced samples; SEM_res = block SEM of the residuals about the fitted line, so the trend does not
                      inflate it); "flat" if |change| < 2 SE(change) (criterion "change_se", D1 rev. 4, the default).
                      The rev. 3 rule of the old project scripts, |change| < 2 x the block SEM of the window mean
                      (criterion "sem"), passes a truly flat series only ~44 % of the time (tests/test_timeseries.py),
                      "change_se" ~95 %. Both verdicts are always reported.
    averaging window  the longest final stretch that passes the drift test, at least 30 % of the run (candidates in 5 %
                      steps of the run length). If even the last 30 % is not flat, the last 30 % is used and the run is
                      flagged not converged.

Primary drift signal: the interior solute count in the MC stream file (high time resolution); secondary: the Stage-2
quantities per snapshot, averaged over the snapshots inside the window found from the stream.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

MIN_BLOCKS = 4  # block levels with fewer block means are not used for the SEM
CRITERIA = ("change_se", "sem")  # drift criteria: D1 rev. 4 (default) and rev. 3 / old scripts
DEFAULT_CRITERION = "change_se"


# ---------------------------------------------------------------------------------------------------- block averaging
@dataclass
class BlockSEM:
    sem: float           # standard error of the mean
    level: int           # blocking level used (block length 2**level samples)
    plateau: bool        # False: no plateau reached -> sem is the largest estimate, probably still too small
    sem_naive: float     # sigma / sqrt(n), uncorrelated samples
    n: int

    @property
    def n_eff(self) -> float:
        """Effective number of independent samples."""
        return self.n * (self.sem_naive / self.sem) ** 2 if self.sem > 0 else float(self.n)


def block_sem(y) -> BlockSEM:
    """Flyvbjerg-Petersen standard error of the mean of a correlated series."""
    x = np.asarray(y, dtype=float)
    n = len(x)
    if n < 2:
        return BlockSEM(sem=float("nan"), level=0, plateau=False, sem_naive=float("nan"), n=n)
    sems, errs = [], []
    while len(x) >= MIN_BLOCKS:
        m = len(x)
        s = np.sqrt(x.var(ddof=1) / m)
        sems.append(s)
        errs.append(s / np.sqrt(2.0 * (m - 1)))
        x = 0.5 * (x[: 2 * (m // 2)][0::2] + x[: 2 * (m // 2)][1::2])
    naive = float(np.std(y, ddof=1) / np.sqrt(n))
    if not sems:  # 2 or 3 samples
        return BlockSEM(sem=naive, level=0, plateau=False, sem_naive=naive, n=n)
    if sems[0] == 0:
        return BlockSEM(sem=0.0, level=0, plateau=True, sem_naive=naive, n=n)
    # plateau: the first level from which the estimate no longer rises by more than its own error
    for k in range(len(sems) - 1):
        if sems[k + 1] - sems[k] < errs[k]:
            return BlockSEM(sem=float(sems[k]), level=k, plateau=True, sem_naive=naive, n=n)
    k = int(np.argmax(sems))
    return BlockSEM(sem=float(sems[k]), level=k, plateau=False, sem_naive=naive, n=n)


# --------------------------------------------------------------------------------------------------------- drift test
@dataclass
class DriftTest:
    t0: float
    t1: float
    n: int
    mean: float
    sem: float           # block-averaging SEM of the window mean
    plateau: bool
    slope: float         # per unit of t
    change: float        # slope x (t1 - t0)
    change_se: float     # standard error of the change, correlation-corrected with the block-averaging SEM
    flat_change_se: bool # D1 rev. 4: |change| < factor x change_se
    flat_sem: bool       # rev. 3 / old scripts: |change| < factor x sem

    def passes(self, criterion: str = DEFAULT_CRITERION) -> bool:
        if criterion == "change_se":
            return self.flat_change_se
        if criterion == "sem":
            return self.flat_sem
        raise ValueError(f"unknown drift criterion {criterion!r} (known: {CRITERIA})")


def drift_test(t, y, factor: float = 2.0) -> DriftTest:
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(t)
    if n < 2:
        return DriftTest(t0=float(t[0]) if n else np.nan, t1=float(t[-1]) if n else np.nan, n=n,
                         mean=float(y.mean()) if n else np.nan, sem=np.nan, plateau=False, slope=np.nan, change=np.nan,
                         change_se=np.nan, flat_change_se=False, flat_sem=False)
    b = block_sem(y)
    tc = t - t.mean()
    sxx = float(tc @ tc)
    slope = float(tc @ (y - y.mean()) / sxx) if sxx > 0 else 0.0
    span = float(t[-1] - t[0])
    change = slope * span
    # SE(change) from the block SEM of the DETRENDED series: a trend would otherwise inflate its own error bar
    b_res = block_sem(y - y.mean() - slope * tc)
    change_se = b_res.sem * span * np.sqrt(n / sxx) if sxx > 0 else np.nan

    def below(limit):
        return bool(np.isfinite(limit) and abs(change) < limit) or bool(limit == 0 and change == 0)

    return DriftTest(t0=float(t[0]), t1=float(t[-1]), n=n, mean=float(y.mean()), sem=b.sem, plateau=b.plateau,
                     slope=slope, change=change, change_se=float(change_se), flat_change_se=below(factor * change_se),
                     flat_sem=below(factor * b.sem))


@dataclass
class Window:
    t_start: float
    t_end: float
    fraction: float      # window length / run length
    converged: bool      # the last min_fraction of the run is flat
    test: DriftTest      # drift test of the chosen window


def averaging_window(t, y, final_fraction: float = 0.30, min_fraction: float = 0.30, factor: float = 2.0,
                     step: float = 0.05, criterion: str = DEFAULT_CRITERION) -> Window:
    """Longest final stretch (>= min_fraction of the run) that passes the drift test."""
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    t_run0, t_run1 = float(t[0]), float(t[-1])
    L = t_run1 - t_run0

    def test(frac):
        m = t >= t_run1 - frac * L - 1e-9 * max(L, 1.0)
        return drift_test(t[m], y[m], factor)

    last = test(final_fraction)
    if not last.passes(criterion):
        return Window(t_start=last.t0, t_end=last.t1, fraction=final_fraction, converged=False, test=last)
    fracs = np.arange(1.0, min_fraction - 1e-9, -step)
    for f in fracs:
        d = test(float(f))
        if d.passes(criterion):
            return Window(t_start=d.t0, t_end=d.t1, fraction=float(f), converged=True, test=d)
    return Window(t_start=last.t0, t_end=last.t1, fraction=final_fraction, converged=True, test=last)


# --------------------------------------------------------------------------------------------------------- MC stream
@dataclass
class Stream:
    path: Path
    columns: list
    data: np.ndarray       # rows sorted by step, duplicate steps (restart overlap) resolved to the LAST row

    def column(self, name) -> np.ndarray:
        return self.data[:, self.columns.index(name)]


def read_mc_stream(path) -> Stream:
    """MC stream file: '# step time[ps] ...' header, then numeric rows. Non-numeric or non-finite rows are skipped."""
    path = Path(path)
    columns, rows = None, {}
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                if columns is None:
                    names = line[1:].split()
                    if names and names[0] == "step":
                        columns = names
                continue
            p = line.split()
            if not p:
                continue
            try:
                v = [float(x) for x in p]
            except ValueError:
                continue
            if columns is not None and len(v) != len(columns):
                continue
            if not all(np.isfinite(v)):
                continue
            rows[int(v[0])] = v
    if columns is None:
        raise ValueError(f"{path}: no '# step ...' header line")
    steps = sorted(rows)
    return Stream(path=path, columns=columns, data=np.array([rows[s] for s in steps], dtype=float))


_COUNT_COLUMN = re.compile(r"^n\w*_(interior|mcgroup)$", re.IGNORECASE)


def solute_count_column(columns, override=None) -> str:
    """The interior solute count: nCu_interior (Ni-Cu), nSol_interior (Nb-Au), NZr_mcgroup (TURLO), unless given."""
    if override:
        if override not in columns:
            raise ValueError(f"column {override!r} not in the MC stream ({columns})")
        return override
    for c in columns:
        if _COUNT_COLUMN.match(c):
            return c
    raise ValueError(f"no interior solute-count column in {columns}; set mc_stream_column in the manifest")


def time_column(columns) -> str:
    for c in columns:
        if c.startswith("time"):
            return c
    return "step"


# ---------------------------------------------------------------------------------------------- Stage 3 of one run
def group_records(rows):
    """Stage-2 rows -> {(quantity, region, basis, ref): (timesteps, values)}, sorted by timestep."""
    out = {}
    for r in rows:
        key = (r["quantity"], r["region"], r["basis"], r.get("ref", ""))
        out.setdefault(key, []).append((int(r["timestep"]), float(r["value"])))
    return {k: tuple(np.array(x) for x in zip(*sorted(v))) for k, v in out.items()}


def run_stage3(stage2_rows, params: dict, stream: Stream | None = None, stream_column: str | None = None,
               md_timestep_ps: float | None = None, key_quantity=("N_sol_inside", "20", "r", ""),
               criterion: str | None = None):
    """Averaging window from the stream (or, without one, from the key Stage-2 quantity), then window means with
    block-averaging SEM and a drift test for every Stage-2 quantity.

    Snapshot times: timestep x md_timestep_ps (ps) if given, else the timestep itself; the stream is put on the same
    axis through its step column.
    criterion: drift criterion (None: the definition set's ``drift_criterion``); the other one is reported alongside.
    Returns (info dict, list of result rows).
    """
    ff, fac, mf = (float(params[k]) for k in ("drift_final_fraction", "drift_factor", "min_window_fraction"))
    criterion = criterion or params.get("drift_criterion", DEFAULT_CRITERION)
    if criterion not in CRITERIA:
        raise ValueError(f"unknown drift criterion {criterion!r} (known: {CRITERIA})")
    series = group_records(stage2_rows)
    scale = float(md_timestep_ps) if md_timestep_ps else 1.0
    info = {"time_unit": "ps" if md_timestep_ps else "timestep"}
    if stream is not None:
        col = solute_count_column(stream.columns, stream_column)
        ts = stream.column("step") * scale
        sig_t, sig_y = ts, stream.column(col)
        info.update(window_source=f"mc_stream:{col}", stream=str(stream.path), stream_samples=len(ts),
                    stream_t_first=float(ts[0]), stream_t_last=float(ts[-1]))
    else:
        if key_quantity not in series:
            raise ValueError(f"no MC stream and no key quantity {key_quantity} in the Stage-2 rows")
        steps, vals = series[key_quantity]
        sig_t, sig_y = steps * scale, vals
        info.update(window_source="stage2:" + "/".join(key_quantity))
    w = averaging_window(sig_t, sig_y, ff, mf, fac, criterion=criterion)
    alt = next(c for c in CRITERIA if c != criterion)
    w_alt = averaging_window(sig_t, sig_y, ff, mf, fac, criterion=alt)
    d = w.test
    info.update(criterion=criterion, window_t_start=w.t_start, window_t_end=w.t_end, window_fraction=w.fraction,
                converged=w.converged, signal_mean=d.mean, signal_sem=d.sem, signal_change=d.change,
                signal_change_se=d.change_se, signal_plateau=d.plateau,
                **{f"alt_{alt}_converged": w_alt.converged, f"alt_{alt}_window_fraction": w_alt.fraction})
    out = []
    for (q, reg, basis, ref), (steps, vals) in sorted(series.items()):
        t = steps * scale
        m = (t >= w.t_start - 1e-9) & (t <= w.t_end + 1e-9)
        if not m.any():
            # no snapshot inside the window (e.g. the stream runs past the last dump): use the last snapshot only
            m = t == t.max()
        tv, yv = t[m], vals[m]
        dt = drift_test(tv, yv, fac)
        out.append({"quantity": q, "region": reg, "basis": basis, "ref": ref, "mean": dt.mean, "sem": dt.sem,
                    "n_samples": int(m.sum()), "sem_plateau": dt.plateau, "change": dt.change, "change_se": dt.change_se,
                    "flat": dt.passes(criterion), "flat_rev3_sem": dt.flat_sem, "t_first": float(tv[0]),
                    "t_last": float(tv[-1])})
    return info, out
