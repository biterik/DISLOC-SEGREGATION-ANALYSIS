#!/usr/bin/env python3
"""Figures of the regression report from regression/results/*.csv. Made at their final slide size, fonts >= 14 pt,
colour-blind-safe (Okabe-Ito, validated), marker shape as second encoding.   python make_figures.py OUTDIR"""

import csv
import sys
from pathlib import Path

import matplotlib
import matplotlib.ticker

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

RES = Path(__file__).parent / "results"
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else RES / "fig"
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.size": 14, "axes.labelsize": 14, "axes.titlesize": 14, "xtick.labelsize": 13,
                     "ytick.labelsize": 13, "legend.fontsize": 13, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": "#e6e6e6", "legend.frameon": False, "savefig.dpi": 200})
C = {"blue": "#0072B2", "orange": "#E69F00", "green": "#009E73", "red": "#D55E00", "ink": "#333333"}
POT = {"LeeShim2004-MEAM": (C["blue"], "o", "Lee-Shim"), "CuNi-ACE": (C["orange"], "s", "ACE"),
       "EAMX-DawChandross2023": (C["green"], "D", "Au"), "Mendelev2022-EAMFS": (C["red"], "^", "Nb")}


def rows(name):
    with open(RES / name) as f:
        return list(csv.DictReader(f))


def fl(r, k):
    try:
        return float(r[k])
    except (KeyError, ValueError, TypeError):
        return float("nan")


def one_to_one(ax, vals, log=True, pad=1.6):
    v = np.asarray([x for x in vals if np.isfinite(x) and (x > 0 or not log)])
    lo, hi = (v.min() / pad, v.max() * pad) if log else (v.min() - 0.05 * np.ptp(v), v.max() + 0.05 * np.ptp(v))
    ax.plot([lo, hi], [lo, hi], color="#999999", lw=1.2, zorder=0)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    if log:
        ax.set_xscale("log")
        ax.set_yscale("log")


def scatter_panel(ax, data, xk, yk, rk=None, log=True, title="", xlabel="old project value", ylabel="disegreg D1"):
    allv = []
    for pot, (col, mk, lab) in POT.items():
        sel = [r for r in data if r.get("pot") == pot and r.get("status", "ok") == "ok"]
        if not sel:
            continue
        x = np.array([fl(r, xk) for r in sel])
        y = np.array([fl(r, yk) for r in sel])
        ax.scatter(x, y, s=70, color=col, marker=mk, edgecolor="white", linewidth=1.5, label=lab, zorder=3)
        if rk:
            ax.scatter(x, [fl(r, rk) for r in sel], s=150, facecolor="none", edgecolor=C["ink"], marker=mk, lw=1.0,
                       zorder=2)
        allv += list(x) + list(y)
    one_to_one(ax, allv, log=log)
    ax.set_title(title, loc="left")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)


def fig_nicu():
    d = rows("nicu_comparison.csv")
    fig, ax = plt.subplots(1, 2, figsize=(7.4, 5.3), constrained_layout=True)  # = size on the slide
    scatter_panel(ax[0], d, "c0_5_old", "c0_5_D1", rk="c0_5_old_recipe", title="c(r < 5 Å)")
    scatter_panel(ax[1], d, "enrich_old", "enrich_C_vs_cB_D1", rk="enrich_old_recipe",
                  title="c(r < 5 Å) / c_far")
    ax[0].legend(loc="upper left")
    ax[1].text(0.97, 0.04, "open: old recipe\non disegreg data", transform=ax[1].transAxes,
               ha="right", va="bottom", fontsize=13, color=C["ink"])
    fig.savefig(OUT / "F1_nicu_old-vs-D1.png")
    plt.close(fig)


def fig_nbau():
    d = rows("nbau_comparison.csv")
    fig, ax = plt.subplots(2, 3, figsize=(8.0, 5.5), constrained_layout=True)
    panels = [("c_core", True, "c core"), ("c_SF", True, "c fault"), ("c_far", True, "c far"),
              ("sep", False, "separation (Å)"), ("a1_core", False, "α₁ core"),
              ("r2_core", False, "R2 core")]
    for a, (k, log, t) in zip(ax.flat, panels):
        sub = [r for r in d if np.isfinite(fl(r, f"{k}_old")) and (fl(r, f"{k}_old") > 0 or not log)]
        scatter_panel(a, sub, f"{k}_old", f"{k}_D1", rk=f"{k}_old_recipe", log=log, title=t, xlabel="old", ylabel="D1")
        if not log:
            a.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3))
            a.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
    ax[0, 0].legend(loc="lower right", fontsize=13, handletextpad=0.1, borderaxespad=0.1)
    fig.savefig(OUT / "F2_nbau_old-vs-D1.png")
    plt.close(fig)


def fig_turlo():
    d = {r["run"]: r for r in rows("turlo_comparison.csv")}
    fig, ax = plt.subplots(1, 2, figsize=(7.6, 5.0), constrained_layout=True)
    labels = ["old", "D1, r", "D1, r_DXA"]
    keys = ["NZr_excess_r<20_old", "NZr_excess_r<20_D1", "NZr_excess_r<20_D1_rdxa"]
    cols = [C["ink"], C["blue"], C["orange"]]
    hatches = ["", "", "//"]
    for a, run, title in ((ax[0], "R3", "R3 (mixed mode)"), (ax[1], "R4", "R4 (polyline)")):
        v = [fl(d[run], k) for k in keys]
        a.bar(range(3), v, color=cols, hatch=hatches, edgecolor="white", width=0.7)
        for i, x in enumerate(v):
            a.text(i, x, f"{x:.0f}" if x > 20 else f"{x:.1f}", ha="center", va="bottom", fontsize=14)
        a.set_xticks(range(3), labels)
        a.set_title(title, loc="left")
        a.set_ylabel("Zr excess, r < 20 Å (atoms)")
        a.grid(axis="x", visible=False)
        a.set_axisbelow(True)
    fig.savefig(OUT / "F3_turlo_R3-R4.png")
    plt.close(fig)


def fig_gamma_xref():
    """Gamma(R) of Ni-Cu d60 Lee-Shim -0.25 kT (last snapshot) for x_ref = x_iso and x_ref = c_B."""
    import glob

    from disegreg.cache import get_snapshot
    from disegreg.definitions import load_definitions
    from disegreg.stage2 import analyse_snapshot

    d = load_definitions()
    dump = glob.glob(str(Path("~/Desktop/SIMULATIONS").expanduser() / "Cu-AT-DISLOCS-Ni-EAM-MEAM-ACE"
                         / "07_MDMC-CU-SEGREGATION-T300-COEX-OFFSETS" / "mdmc300K_d60_LeeShim2004-MEAM_0.25kT*"
                         / "*_quenched_step830000.dump"))[0]
    snap = get_snapshot(dump, d)
    res = analyse_snapshot(snap, d, 3.53481, z_gp=125.5106, x_ref=6.8526e-3)
    L, cB = res.value("L_x"), res.value("c", "B")
    R = [5, 10, 20, 40, 60, 80, 100]
    n = np.array([res.value("N_inside", f"{r}", "r") for r in R])
    ns = np.array([res.value("N_sol_inside", f"{r}", "r") for r in R])
    fig, ax = plt.subplots(figsize=(8.2, 5.0), constrained_layout=True)
    for x, col, mk, lab in ((6.8526e-3, C["orange"], "s", f"x_ref = x_iso = {6.8526e-3:.5f}"),
                            (cB, C["blue"], "o", f"x_ref = c_B = {cB:.5f}")):
        ax.plot(R, (ns - x * n) / L, color=col, marker=mk, ms=9, lw=2, label=lab, mec="white", mew=1.5)
    ax.axhline(0, color="#999999", lw=1.2, zorder=0)
    ax.set_xlabel("R (Å)")
    ax.set_ylabel("Γ(R) (Cu atoms per Å of line)")
    ax.set_title("Ni-Cu d60, Lee-Shim, Δ = −0.25", loc="left")
    ax.legend(loc="lower left")
    fig.savefig(OUT / "F4_nicu_gamma-vs-xref.png")
    plt.close(fig)


if __name__ == "__main__":
    for f in (fig_nicu, fig_nbau, fig_turlo, fig_gamma_xref):
        try:
            f()
            print("ok", f.__name__)
        except FileNotFoundError as e:
            print("skip", f.__name__, e)
