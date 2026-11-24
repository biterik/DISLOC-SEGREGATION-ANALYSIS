"""Plots in Erik's house style (docs/PLAN_v0.2.md rule 6; CLAUDE.md): every font >= 13 pt at slide size, colour-blind-safe
palette, SYMMETRIC colour scale for every map of a signed quantity. matplotlib only (Agg backend, no screen needed).

Palette: Okabe & Ito (2008), order blue, vermillion, green, orange, pink, sky blue -- checked with a colour-vision-
deficiency simulator (worst neighbouring pair Delta E 9.6 for deuteranopia, >= 8 required). Three of the colours are
light on white (contrast < 3:1), so every series also carries its own marker shape and is labelled directly.
Colour follows the entity, never its rank: d0 is always blue, d30 vermillion, d60 green, d90 orange.
Signed maps: blue (negative) -- grey (zero) -- red (positive), limits +-max|value|.

    gamma_curve      Gamma(R) of one run, every reference, under r (solid) and r_DXA (dashed)
    timeseries       one Stage-2 quantity vs time with the averaging window, its mean +- SEM and the drift verdict
    vs_delta         a results-table quantity vs the offset Delta, one series per character (colour) and potential (marker)
    uptake           Stage-4 uptake exponents vs Delta, one panel per reference, with the McLean bound n = 1
    regions          concentration per region of one run (bar chart, log scale) with the reference compositions
    signed_map       any 2-D map of a signed quantity (e.g. solute excess in the Y-Z plane), symmetric colour scale
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

PALETTE = ("#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7", "#56B4E9")
CHARACTER_COLOUR = {"d0": PALETTE[0], "d30": PALETTE[1], "d60": PALETTE[2], "d90": PALETTE[3]}
CHARACTER_NAME = {"d0": "screw (0°)", "d30": "30°", "d60": "60°", "d90": "edge (90°)"}
MARKERS = ("o", "s", "^", "D", "v", "P")
DIVERGING = ("#2a78d6", "#f0efec", "#e34948")   # blue - neutral grey - red
TEXT, TEXT_2, GRID = "#0b0b0b", "#52514e", "#d9d8d4"
MIN_FONT = 13.0
SLIDE = (10.0, 5.625)                             # inches: 16:9 slide at full width
REGION_LABEL = {"C": "core\nC", "F": "fault\nF", "S": "shell\nS", "T": "transition\nT", "B": "far field\nB",
                "O": "far, other\nO", "interior": "bulk\ninterior"}

STYLE = {
    "font.size": 15, "axes.titlesize": 16, "axes.labelsize": 15, "xtick.labelsize": 14, "ytick.labelsize": 14,
    "legend.fontsize": 13, "legend.title_fontsize": 13, "figure.titlesize": 17, "font.family": "sans-serif",
    "axes.edgecolor": TEXT_2, "axes.labelcolor": TEXT, "xtick.color": TEXT_2, "ytick.color": TEXT_2,
    "text.color": TEXT, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8, "axes.axisbelow": True,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2.0, "lines.markersize": 8,
    "legend.frameon": False, "figure.dpi": 100, "savefig.dpi": 200, "savefig.bbox": "tight",
    "axes.prop_cycle": None,  # filled below
}


def _mpl():
    import matplotlib

    matplotlib.use("Agg", force=False)
    import matplotlib.pyplot as plt
    from cycler import cycler

    style = dict(STYLE, **{"axes.prop_cycle": cycler(color=list(PALETTE))})
    return matplotlib, plt, style


def diverging_cmap():
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list("disegreg_diverging", list(DIVERGING), N=256)


def symmetric_limits(values) -> tuple[float, float]:
    """(-m, m) with m = max |finite value| (1 if there is none): the colour scale of a signed map."""
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    m = float(np.abs(v).max()) if v.size and np.abs(v).max() > 0 else 1.0
    return -m, m


def check_fonts(fig, minimum: float = MIN_FONT) -> list[str]:
    """Texts of the figure smaller than ``minimum`` pt (empty = house style kept)."""
    small = []
    for t in fig.findobj(lambda a: hasattr(a, "get_fontsize") and hasattr(a, "get_text")):
        if t.get_text().strip() and t.get_visible() and t.get_fontsize() < minimum - 1e-6:
            small.append(f"{t.get_text()!r}: {t.get_fontsize():.1f} pt")
    return small


def save(fig, out, formats=("png",)) -> list[Path]:
    """Write the figure (after checking the font rule) as <out>.<fmt> for every format; returns the paths."""
    small = check_fonts(fig)
    if small:
        raise ValueError(f"fonts below {MIN_FONT} pt: {small[:5]}")
    out = Path(out)
    paths = []
    for fmt in formats:
        p = out.with_suffix(f".{fmt}")
        p.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(p)
        paths.append(p)
    _mpl()[1].close(fig)
    return paths


def _figure(nrows=1, ncols=1, size=SLIDE, **kw):
    matplotlib, plt, style = _mpl()
    with matplotlib.rc_context(style):
        fig, ax = plt.subplots(nrows, ncols, figsize=size, constrained_layout=True, **kw)
    return fig, ax, style


def _ctx():
    matplotlib, _, style = _mpl()
    return matplotlib.rc_context(style)


def _ref_label(ref: str) -> str:
    return {"c_B": "vs far field $c_B$", "x_iso": "vs isotherm $x_{iso}$", "x_bc": "vs bulk control $x_{bc}$",
            "given": "vs given $x_{ref}$"}.get(ref, ref)


# ----------------------------------------------------------------------------------------------------- Gamma(R)
def gamma_curve(curve_rows, title: str = "", headline_radius: float = 20.0, check_radius: float = 40.0):
    """Gamma(R) from a <run>_gammaR_<D>.csv (results.gamma_curve rows)."""
    with _ctx():
        fig, ax, _ = _figure()
        refs = sorted({r["ref"] for r in curve_rows})
        for k, ref in enumerate(refs):
            col, mk = PALETTE[k % len(PALETTE)], MARKERS[k % len(MARKERS)]
            for basis, ls in (("r", "-"), ("r_dxa", "--")):
                rs = sorted((r for r in curve_rows if r["ref"] == ref and r["basis"] == basis), key=lambda r: r["R"])
                if not rs:
                    continue
                R = np.array([r["R"] for r in rs])
                g = np.array([r["mean"] for r in rs])
                e = np.array([r["sem"] for r in rs])
                ax.errorbar(R, g, yerr=np.where(np.isfinite(e), e, 0), color=col, ls=ls, marker=mk if basis == "r" else None,
                            capsize=3, label=f"{_ref_label(ref)}, " + ("r (used)" if basis == "r" else "$r_{DXA}$"))
        ax.axhline(0, color=TEXT_2, lw=1)
        ax.axvline(headline_radius, color=TEXT_2, lw=1, ls=":")
        ax.annotate(f"headline Γ({headline_radius:g} Å)", (headline_radius, 1), xycoords=("data", "axes fraction"),
                    xytext=(4, -4), textcoords="offset points", va="top", fontsize=13, color=TEXT_2)
        ax.axvline(check_radius, color=GRID, lw=1, ls=":")
        ax.set_xlabel("radius R from the line (Å)")
        ax.set_ylabel("line excess Γ(R) (solute atoms / Å)")
        ax.set_title(title or "Line excess Γ(R)")
        ax.legend(loc="best")
    return fig


# --------------------------------------------------------------------------------------------------- time series
def timeseries(t, y, window=None, mean=None, sem=None, flat=None, ylabel="", title="", time_label="time (ns)",
               slope=None):
    """One quantity vs time; ``window`` = (t_start, t_end) shaded, its mean +- SEM as a band, drift verdict in the title."""
    t, y = np.asarray(t, float), np.asarray(y, float)
    with _ctx():
        fig, ax, _ = _figure()
        if window is not None:
            ax.axvspan(window[0], window[1], color=GRID, alpha=0.5, lw=0, label="averaging window")
        ax.plot(t, y, color=PALETTE[0], marker="o", ms=5 if len(t) > 50 else 8, lw=1.5 if len(t) > 50 else 2,
                label="snapshots" if len(t) <= 50 else "stream")
        if mean is not None and window is not None:
            ax.hlines(mean, window[0], window[1], color=PALETTE[1], lw=2.5, label="window mean")
            if sem is not None and np.isfinite(sem):
                ax.fill_between([window[0], window[1]], mean - sem, mean + sem, color=PALETTE[1], alpha=0.25, lw=0,
                                label="± SEM (block averaging)")
            if slope is not None and np.isfinite(slope):
                tc = 0.5 * (window[0] + window[1])
                ax.plot(window, [mean + slope * (window[0] - tc), mean + slope * (window[1] - tc)], color=PALETTE[2],
                        ls="--", lw=2, label="linear trend")
        verdict = "" if flat is None else ("  —  flat (no drift detected)" if flat else "  —  DRIFT")
        ax.set_title((title or ylabel) + verdict)
        ax.set_xlabel(time_label)
        ax.set_ylabel(ylabel)
        ax.legend(loc="best")
    return fig


# --------------------------------------------------------------------------------------- quantity vs Delta, uptake
def vs_delta(rows, ylabel: str, title: str = "", log: bool = False):
    """rows: dicts with character, potential, delta, value, error (e.g. results-table rows of one quantity)."""
    pots = sorted({r["potential"] for r in rows})
    with _ctx():
        fig, ax, _ = _figure()
        for ch in ("d0", "d30", "d60", "d90"):
            for k, pot in enumerate(pots):
                rs = sorted((r for r in rows if r["character"] == ch and r["potential"] == pot),
                            key=lambda r: float(r["delta"]))
                if not rs:
                    continue
                x = np.array([float(r["delta"]) for r in rs])
                v = np.array([float(r["value"]) for r in rs])
                e = np.array([float(r["error"]) if r["error"] != "" else np.nan for r in rs])
                ax.errorbar(x, v, yerr=np.where(np.isfinite(e), e, 0), color=CHARACTER_COLOUR[ch],
                            marker=MARKERS[k % len(MARKERS)], ls="-" if k == 0 else "--", capsize=3,
                            label=f"{CHARACTER_NAME[ch]}, {pot}" if len(pots) > 1 else CHARACTER_NAME[ch])
        if log:
            ax.set_yscale("log")
        ax.set_xlabel("offset Δ = (Δμ − Δμ$_{coex}$) / kT")
        ax.set_ylabel(ylabel)
        ax.set_title(title or ylabel)
        ax.legend(loc="best", ncol=2 if len(pots) > 1 else 1)
    return fig


def uptake(stage4_rows, title: str = "Core uptake exponent between neighbouring offsets"):
    """Stage-4 rows of analysis 'uptake' (n_up); x = midpoint of the two offsets; one panel per reference (no shared
    second axis), the McLean bound n_up = 1 as a line."""
    rows = [r for r in stage4_rows if r["analysis"] == "uptake"]
    refs = [x for x in ("c_B", "x_iso", "x_bc") if any(r["ref"] == x for r in rows)] or sorted({r["ref"] for r in rows})
    pots = sorted({r["potential"] for r in rows})
    with _ctx():
        fig, axes, _ = _figure(1, max(1, len(refs)), size=(SLIDE[0], SLIDE[1]), sharey=True, squeeze=False)
        for ax, ref in zip(axes[0], refs):
            for ch in ("d0", "d30", "d60", "d90"):
                for k, pot in enumerate(pots):
                    rs = sorted((r for r in rows if r["ref"] == ref and r["character"] == ch and r["potential"] == pot
                                 and r["family"] == "up"), key=lambda r: float(r["delta"]))
                    if not rs:
                        continue
                    x = np.array([0.5 * (float(r["delta"]) + float(r["delta_2"])) for r in rs])
                    v = np.array([float(r["value"]) for r in rs])
                    e = np.array([float(r["error"]) for r in rs])
                    ax.errorbar(x, v, yerr=np.where(np.isfinite(e), e, 0), color=CHARACTER_COLOUR[ch],
                                marker=MARKERS[k % len(MARKERS)], ls="-" if k == 0 else "--", capsize=3,
                                label=f"{CHARACTER_NAME[ch]}, {pot}" if len(pots) > 1 else CHARACTER_NAME[ch])
            ax.axhline(1.0, color=TEXT_2, lw=1.5, ls=":")
            ax.annotate("independent sites: n ≤ 1", (0, 1.0), xycoords=("axes fraction", "data"), xytext=(4, 4),
                        textcoords="offset points", fontsize=13, color=TEXT_2)
            ax.set_title(_ref_label(ref))
            ax.set_xlabel("mid offset Δ")
        axes[0][0].set_ylabel("uptake exponent $n_{up}$ = Δ ln $c_C$ / Δ ln $x_{ref}$")
        axes[0][-1].legend(loc="best", ncol=1)
        fig.suptitle(title)
    return fig


# -------------------------------------------------------------------------------------------------------- regions
def regions(rows, run_id: str, refs: dict | None = None, title: str = ""):
    """Concentration c per region of one run (results-table rows, quantity c), log scale; reference compositions as
    horizontal lines."""
    order = ["C", "F", "S", "T", "B", "O", "interior"]
    rs = {r["region"]: r for r in rows if r["run_id"] == run_id and r["quantity"] == "c" and r["region"] in order}
    regs = [g for g in order if g in rs and float(rs[g]["value"]) > 0]
    with _ctx():
        fig, ax, _ = _figure()
        x = np.arange(len(regs))
        v = np.array([float(rs[g]["value"]) for g in regs])
        e = np.array([float(rs[g]["error"]) if rs[g]["error"] != "" else np.nan for g in regs])
        ax.bar(x, v, width=0.7, color=PALETTE[0], edgecolor="white", linewidth=2)
        ax.errorbar(x, v, yerr=np.where(np.isfinite(e), e, 0), color=TEXT, ls="none", capsize=4)
        for k, (name, val) in enumerate((refs or {}).items()):
            ax.axhline(val, color=PALETTE[1 + k % 5], ls=("--", ":", "-.")[k % 3], lw=2, label=f"{name} = {val:.3g}")
        for xi, vi in zip(x, v):
            ax.annotate(f"{vi:.3g}", (xi, vi), xytext=(0, 4), textcoords="offset points", ha="center", fontsize=13)
        ax.set_xticks(x, [REGION_LABEL[g] for g in regs])
        ax.set_yscale("log")
        ax.set_ylabel("solute fraction c")
        ax.set_title(title or f"Solute fraction per region — {run_id}")
        if refs:
            ax.legend(loc="best")
    return fig


# ---------------------------------------------------------------------------------------------------- signed map
def signed_map(x_edges, y_edges, values, xlabel: str, ylabel: str, cbar_label: str, title: str = "", aspect="equal"):
    """2-D map of a signed quantity on bin edges; colour scale symmetric about zero (+-max|value|), grey at zero."""
    from matplotlib.colors import Normalize

    v = np.asarray(values, float)
    lo, hi = symmetric_limits(v)
    with _ctx():
        fig, ax, _ = _figure()
        m = ax.pcolormesh(x_edges, y_edges, v.T, cmap=diverging_cmap(), norm=Normalize(lo, hi), shading="flat")
        ax.set_aspect(aspect)
        ax.grid(False)
        cb = fig.colorbar(m, ax=ax)
        cb.set_label(cbar_label)
        cb.ax.tick_params(labelsize=14)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title or cbar_label)
    return fig
