"""Along-line core composition and its Fourier spectrum against binomial noise ([PROPOSAL v0.2] P4). Pure numpy.

Question (paper diagnostic for the "line defect phase"): does the core hold a CONTINUOUS line state (solute spread evenly
along the line) or a PERIODIC ARRAY of particles (solute-rich stretches separated by solute-poor ones)?

    bins       the core atoms (region C, r < 5 A) in n = round(L_x / 5 A) bins along X (actual width L_x / n)
    profile    N_k core atoms and s_k solute atoms per bin; c = sum s_k / sum N_k; fluctuation d_k = s_k - c N_k
    dispersion D = sum d_k^2 / ((1 - 1/n) sum N_k c (1 - c)): mean 1 for independent (binomial) placement (the factor
               1 - 1/n: one degree of freedom goes into c), >> 1 for particles (solute bunched in some bins), < 1 for a
               more even than random spread (an ordered, continuous line state)
    spectrum   S(q) = |sum_k d_k exp(-2 pi i q k / n)|^2 / (n Nbar c (1 - c)), q = 1 .. floor(n/2), Nbar = mean N_k;
               for binomial placement every S(q) has mean ~1 (exponentially distributed, independent between q).
               A periodic array of particles with spacing lambda gives a peak at q = L_x / lambda.
    peak test  S_max = max_q S(q); p = 1 - (1 - exp(-S_max))^m2 x (1 - erfc(sqrt(S_max / 2)))^m1: the chance that binomial
               noise alone gives a peak this high somewhere in the spectrum; m2 modes with two degrees of freedom (exp(1)),
               m1 = 1 for even n (the mode q = n/2 is real: S = Z^2, one degree of freedom), else 0.

Averaging S(q) over the snapshots of the window (Stage 3 does this for every Stage-2 quantity) sharpens the test: the mean
of n_s independent exp(1) values has standard deviation 1/sqrt(n_s). Snapshots of one window are correlated, so the test
(D2 rev. 3, section 16) uses the effective number of snapshots n_eff = n_s (SEM_naive / SEM_block)^2 of the block averaging
of the core concentration c_C over the same snapshots; without a plateau n_eff = 1 (flagged).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .geometry import Cell, n_slabs_for

# per-snapshot peak values: meaningful per snapshot only (their window mean is not the section-16 test) -> not in Stage 3
PER_SNAPSHOT_ONLY = ("alongline_S_max", "alongline_q_peak", "alongline_lambda_peak", "alongline_p_peak")
# the section-16 test on the window-averaged spectrum (Stage 3, window_test)
WINDOW_QUANTITIES = ("alongline_window_S_max", "alongline_window_q_peak", "alongline_window_lambda_peak",
                     "alongline_window_p_peak", "alongline_window_periodic", "alongline_window_n_eff")
QUANTITIES = ("alongline_n_bins", "alongline_bin_width", "alongline_c", "alongline_dispersion", "alongline_S",
              "alongline_S_max", "alongline_q_peak", "alongline_lambda_peak", "alongline_p_peak", "alongline_N_bin",
              "alongline_N_sol_bin") + WINDOW_QUANTITIES


@dataclass
class AlongLine:
    n_bins: int
    bin_width: float
    N: np.ndarray            # core atoms per bin
    s: np.ndarray            # core solute atoms per bin
    c: float                 # core solute fraction
    dispersion: float
    q: np.ndarray            # 1 .. floor(n/2)
    S: np.ndarray            # normalised spectrum
    S_max: float
    q_peak: int
    lambda_peak: float       # L_x / q_peak (A)
    p_peak: float


def binomial_peak_p(S_max: float, m: int, m1: int = 0) -> float:
    """P(max of m independent exp(1) variables and m1 independent chi^2_1 variables (mean 1) >= S_max)."""
    from math import erfc, sqrt

    if not np.isfinite(S_max) or m + m1 <= 0:
        return float("nan")
    if S_max <= 0:
        return 1.0
    log_below = m * np.log1p(-np.exp(-S_max)) + m1 * np.log1p(-erfc(sqrt(S_max / 2.0)))
    return float(-np.expm1(log_below))


def profile(cell: Cell, positions, solute, mask, bin_width: float = 5.0) -> AlongLine:
    """Along-line profile and spectrum of the atoms in ``mask`` (the core)."""
    pos = np.asarray(positions, float)[np.asarray(mask, bool)]
    sol = np.asarray(solute, bool)[np.asarray(mask, bool)]
    n = n_slabs_for(cell, bin_width)
    L = cell.line_length
    sx = cell.wrap_fractional(cell.to_fractional(pos))[:, 0] if len(pos) else np.empty(0)
    k = np.minimum((sx * n).astype(np.int64), n - 1)
    N = np.bincount(k, minlength=n).astype(float)
    s = np.bincount(k[sol], minlength=n).astype(float)
    nan = float("nan")
    q = np.arange(1, n // 2 + 1)
    if N.sum() == 0:
        return AlongLine(n, L / n, N, s, nan, nan, q, np.full(len(q), nan), nan, 0, nan, nan)
    c = s.sum() / N.sum()
    d = s - c * N
    var = N.sum() * c * (1 - c)
    if var <= 0:  # core entirely solute or entirely host: no fluctuation to analyse
        return AlongLine(n, L / n, N, s, float(c), nan, q, np.full(len(q), nan), nan, 0, nan, nan)
    disp = float((d * d).sum() / ((1.0 - 1.0 / n) * var)) if n > 1 else float("nan")
    F = np.fft.fft(d)
    S = (np.abs(F[q]) ** 2) / (n * (N.mean() * c * (1 - c)))
    j = int(np.argmax(S)) if len(S) else -1
    smax = float(S[j]) if j >= 0 else nan
    qp = int(q[j]) if j >= 0 else 0
    return AlongLine(n_bins=n, bin_width=L / n, N=N, s=s, c=float(c), dispersion=disp, q=q, S=S, S_max=smax,
                     q_peak=qp, lambda_peak=L / qp if qp else nan,
                     p_peak=binomial_peak_p(smax, len(q) - (n % 2 == 0), int(n % 2 == 0)))


def records(al: AlongLine, add) -> None:
    """Write the Stage-2 records through ``add(quantity, value, region=...)`` (stage2.Stage2Result.add)."""
    add("alongline_n_bins", al.n_bins, region="C")
    add("alongline_bin_width", al.bin_width, region="C")
    add("alongline_c", al.c, region="C")
    add("alongline_dispersion", al.dispersion, region="C")
    add("alongline_S_max", al.S_max, region="C")
    add("alongline_q_peak", al.q_peak, region="C")
    add("alongline_lambda_peak", al.lambda_peak, region="C")
    add("alongline_p_peak", al.p_peak, region="C")
    for qq, v in zip(al.q, al.S):
        add("alongline_S", v, region=f"q{int(qq)}")
    for kk, (nn, ss) in enumerate(zip(al.N, al.s)):
        add("alongline_N_bin", nn, region=f"x{kk}")
        add("alongline_N_sol_bin", ss, region=f"x{kk}")


def window_peak_p(S_max: float, n_snap: float, m: int, m1: int = 0) -> float:
    """P(max over the modes of the WINDOW-AVERAGED spectrum >= S_max) under binomial noise, for n_snap independent
    snapshots: the mean of n_snap exp(1) variables is Gamma(n_snap, 1/n_snap), the mean of n_snap chi^2_1 variables (the
    q = n/2 mode of an even bin number) Gamma(n_snap / 2, 2 / n_snap). n_snap = 1 gives binomial_peak_p. n_snap may be
    an effective (non-integer) number of snapshots, n_eff (window_n_eff)."""
    from scipy.stats import gamma

    if not np.isfinite(S_max) or m + m1 <= 0 or n_snap < 1:
        return float("nan")
    if S_max <= 0:
        return 1.0
    below = m * np.log(gamma.cdf(S_max, a=n_snap, scale=1.0 / n_snap)) if m else 0.0
    if m1:
        below += m1 * np.log(gamma.cdf(S_max, a=n_snap / 2.0, scale=2.0 / n_snap))
    return float(-np.expm1(below))


def window_n_eff(core_c, n_snap: int, without_plateau: float = 1.0) -> dict:
    """Effective number of independent snapshots for the window test (D2 rev. 3, section 16, Erik 26 Nov, N3 option b):
    n_eff = n_snap x (SEM_naive / SEM_block)^2 of the block averaging (Flyvbjerg-Petersen, timeseries.block_sem) of the
    core concentration c_C over the window snapshots ``core_c``, clipped to [1, n_snap]. No plateau of the block averaging
    (always so for 2-3 snapshots): n_eff = 1 and a flag. One snapshot: n_eff = 1 exactly, no flag.
    ``core_c`` None: snapshots treated as independent (n_eff = n_snap) -- synthetic tests only."""
    from .timeseries import block_sem

    out = {"n_snapshots": int(n_snap), "n_eff_source": "block averaging of c (region C)", "c_C_plateau": "", "flag": ""}
    if core_c is None:
        out.update(n_eff=float(n_snap), n_eff_source="snapshots treated as independent")
        return out
    if n_snap <= 1:
        out.update(n_eff=1.0)
        return out
    b = block_sem(core_c)
    out["c_C_plateau"] = bool(b.plateau)
    if not b.plateau:
        out.update(n_eff=float(without_plateau), flag=f"block averaging of c_C without plateau ({len(core_c)} snapshots): "
                   f"n_eff = {without_plateau:g} in the along-line window test (D2 section 16)")
        return out
    f = (b.sem_naive / b.sem) ** 2 if b.sem > 0 else 1.0   # c_C constant in the window (SEM 0): no correlation measured
    out["n_eff"] = float(min(max(n_snap * f, 1.0), float(n_snap)))
    return out


def window_test(stage3_rows, p_threshold: float, core_c, info: dict | None = None,
                n_eff_without_plateau: float = 1.0) -> list[dict]:
    """The section-16 test (D2 rev. 3) on the window-averaged spectrum: from the Stage-3 window means of alongline_S
    (regions q1 .. q<n/2>) and their sample count n_s, the highest mode, its wavelength L_x / q, its p-value under binomial
    noise for n_eff effective snapshots (window_peak_p, window_n_eff; ``core_c`` = c_C of the window snapshots) and the
    verdict 1 = periodic (p < p_threshold) / 0 / NaN. Returns Stage-3-style rows (mean = the value; no SEM, no drift test;
    sem_plateau = the plateau of the c_C block averaging); ``info`` (Stage-3 info) gets the n_eff record."""
    S = {int(r["region"][1:]): r for r in stage3_rows if r["quantity"] == "alongline_S"}
    nb = next((r for r in stage3_rows if r["quantity"] == "alongline_n_bins"), None)
    L = next((r for r in stage3_rows if r["quantity"] == "L_x"), None)
    if not S or nb is None:
        return []
    n_bins = int(round(float(nb["mean"])))
    q = sorted(S)
    vals = np.array([float(S[k]["mean"]) for k in q])
    n_snap = int(min(int(S[k]["n_samples"]) for k in q))
    proto = S[q[0]]
    ne = window_n_eff(core_c, n_snap, n_eff_without_plateau)
    if info is not None:
        info["alongline_window"] = ne
    if not np.isfinite(vals).any():
        smax, qp, p = float("nan"), 0, float("nan")
    else:
        j = int(np.nanargmax(vals))
        smax, qp = float(vals[j]), int(q[j])
        m1 = int(n_bins % 2 == 0 and n_bins // 2 in S)
        p = window_peak_p(smax, ne["n_eff"], len(q) - m1, m1)
    lam = float(L["mean"]) / qp if (L is not None and qp) else float("nan")
    verdict = float(p < p_threshold) if np.isfinite(p) else float("nan")
    out = []
    for quantity, value in (("alongline_window_S_max", smax), ("alongline_window_q_peak", float(qp) if qp else float("nan")),
                            ("alongline_window_lambda_peak", lam), ("alongline_window_p_peak", p),
                            ("alongline_window_periodic", verdict), ("alongline_window_n_eff", ne["n_eff"])):
        out.append({"quantity": quantity, "region": "C", "basis": "", "ref": "", "mean": value, "sem": float("nan"),
                    "n_samples": n_snap, "sem_plateau": ne["c_C_plateau"], "change": float("nan"),
                    "change_se": float("nan"), "flat": "", "flat_rev3_sem": "", "t_first": proto.get("t_first"),
                    "t_last": proto.get("t_last")})
    return out

