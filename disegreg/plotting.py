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
    out = str(out)
    for fmt in formats:  # a name given WITH its extension: drop it (run ids contain dots, so no with_suffix)
        if out.endswith(f".{fmt}"):
            out = out[: -len(fmt) - 1]
    paths = []
    for fmt in formats:
        p = Path(f"{out}.{fmt}")
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


def _wrap(text: str, width: int = 62) -> str:
    """Break a long title (run ids are ~60 characters) at '_' / ' ' so it fits the slide width at >= 13 pt."""
    import textwrap

    out = []
    for line in text.split("\n"):
        if len(line) <= width:
            out.append(line)
            continue
        out += textwrap.wrap(line.replace("_", "_\u200b"), width=width, break_on_hyphens=False)
    return "\n".join(x.replace("\u200b", "") for x in out)


def _value_limits(ax, values, include=(), pad: float = 0.15):
    """y range from the VALUES (and ``include``), not from the error bars: one noisy pair must not flatten the others.
    Returns True if some error bar is clipped (then say so in the plot)."""
    v = np.asarray([x for x in list(values) + list(include) if np.isfinite(x)], float)
    if not v.size:
        return False
    lo, hi = float(v.min()), float(v.max())
    span = hi - lo if hi > lo else max(abs(hi), 1.0)
    ax.set_ylim(lo - pad * span, hi + pad * span)
    return True


QUANTITY_LABEL = {"c:C": "core solute fraction $c_C$", "c:B": "far-field solute fraction $c_B$",
                  "Gamma:20:r": "line excess Γ(20 Å) (atoms / Å)", "Gamma:40:r": "line excess Γ(40 Å) (atoms / Å)",
                  "E:C": "core enrichment $E_C$", "N_sol_inside:20:r": "solute atoms within 20 Å",
                  "partial_separation_mean": "partial separation (Å)", "alongline_dispersion:C": "along-line dispersion D"}


def quantity_label(key: str) -> str:
    """'quantity:region:basis:ref' -> readable axis label (the key itself if unknown)."""
    k = key.rstrip(":")
    parts = k.split(":")
    ref = parts[3] if len(parts) > 3 else ""
    base = QUANTITY_LABEL.get(":".join(parts[:3]).rstrip(":"), k)
    return base + (f", {_ref_label(ref)}" if ref else "")


def _legend_row(fig_axes_names, ncols=1, height=(5, 1.1), size=None):
    """Figure with the plot panels on top and a row reserved for the legend (it never covers data or axis labels)."""
    _, plt, _ = _mpl()
    names = list(fig_axes_names)
    fig, ax_d = plt.subplot_mosaic([names, ["legend"] * len(names)], figsize=size or (SLIDE[0], SLIDE[1] + 1.4),
                                   height_ratios=list(height), constrained_layout=True)
    ax_d["legend"].axis("off")
    return fig, [ax_d[n] for n in names], ax_d["legend"]


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
        ax.set_title(_wrap(title or "Line excess Γ(R)"), fontsize=14)
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
        ax.set_title(_wrap((title or ylabel) + verdict), fontsize=14)
        ax.set_xlabel(time_label)
        ax.set_ylabel(ylabel)
        ax.legend(loc="best")
    return fig


# --------------------------------------------------------------------------------------- quantity vs Delta, uptake
def vs_delta(rows, ylabel: str, title: str = "", log: bool = False):
    """rows: dicts with character, potential, delta, value, error (e.g. results-table rows of one quantity)."""
    pots = sorted({r["potential"] for r in rows})
    with _ctx():
        fig, (ax,), leg = _legend_row(["main"])
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
        ax.set_title(_wrap(title or ylabel), fontsize=14)
        h, lab = ax.get_legend_handles_labels()
        leg.legend(h, lab, loc="center", ncol=2 if len(lab) > 4 else max(1, len(lab)), frameon=False)
    return fig


def uptake(stage4_rows, title: str = "Core uptake exponent between neighbouring offsets"):
    """Stage-4 rows of analysis 'uptake' (n_up); x = midpoint of the two offsets; one panel per reference (no shared
    second axis), the McLean bound n_up = 1 as a line."""
    rows = [r for r in stage4_rows if r["analysis"] == "uptake"]
    refs = [x for x in ("c_B", "x_iso", "x_bc") if any(r["ref"] == x for r in rows)] or sorted({r["ref"] for r in rows})
    pots = sorted({r["potential"] for r in rows})
    with _ctx():
        _, plt, _ = _mpl()
        names = [f"p{k}" for k in range(max(1, len(refs)))]
        fig, ax_d = plt.subplot_mosaic([names, ["legend"] * len(names)], figsize=(SLIDE[0], SLIDE[1] + 1.4),
                                       height_ratios=[5, 1.1], constrained_layout=True)
        axes = [[ax_d[n] for n in names]]
        for ax in axes[0][1:]:
            ax.sharey(axes[0][0])
            ax.tick_params(labelleft=False)
        leg_ax = ax_d["legend"]
        leg_ax.axis("off")
        n_unknown = [0]
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
                    unk = ~np.isfinite(e)
                    if unk.any():  # error unknown (e.g. a single snapshot): open symbol, not a seemingly exact point
                        ax.plot(x[unk], v[unk], ls="none", marker=MARKERS[k % len(MARKERS)], mfc="white",
                                mec=CHARACTER_COLOUR[ch], mew=2, ms=9, zorder=5)
                        n_unknown[0] += int(unk.sum())
            ax.axhline(1.0, color=TEXT_2, lw=1.5, ls=":")
            ax.annotate("independent sites: n ≤ 1", (0, 1.0), xycoords=("axes fraction", "data"), xytext=(4, 4),
                        textcoords="offset points", fontsize=13, color=TEXT_2)
            ax.set_title(_ref_label(ref))
            ax.set_xlabel("mid offset Δ")
        # y range without the noise-dominated exponents (known error > max(1, |n| / 2)); those are
        # (e.g. a far field with a handful of solute atoms) and would flatten every other point
        up = [r for r in rows if r["family"] == "up" and np.isfinite(float(r["value"]))]
        # an UNKNOWN error (single snapshot) does not make a value noisy -- only a known large error does
        good = [float(r["value"]) for r in up
                if not np.isfinite(float(r["error"])) or float(r["error"]) < max(1.0, 0.5 * abs(float(r["value"])))]
        _value_limits(axes[0][0], good or [float(r["value"]) for r in up], include=(0.0, 1.0))
        lo, hi = axes[0][0].get_ylim()
        off = sum(1 for r in up if not lo <= float(r["value"]) <= hi)
        clipped = sum(1 for r in up if lo <= float(r["value"]) <= hi and np.isfinite(float(r["error"]))
                      and (float(r["value"]) - float(r["error"]) < lo or float(r["value"]) + float(r["error"]) > hi))
        notes = ["open symbols: error unknown (single snapshot)"] if n_unknown[0] else []
        note = "; ".join(x for x in notes + [f"{off} poorly determined value{'s' if off != 1 else ''} off the plot" if off else "",
                                     f"{clipped} error bar{'s' if clipped != 1 else ''} clipped" if clipped else ""] if x)
        if note:
            axes[0][0].annotate(_wrap(note, 34), (0.02, 0.02), xycoords="axes fraction", ha="left", va="bottom", fontsize=13,
                                color=TEXT_2, bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 2})
        from matplotlib.ticker import MaxNLocator
        for ax in axes[0]:
            ax.xaxis.set_major_locator(MaxNLocator(4))
        axes[0][0].set_ylabel("uptake exponent $n_{up}$ = Δ ln $c_C$ / Δ ln $x_{ref}$")
        handles, labels = axes[0][0].get_legend_handles_labels()
        for ax in axes[0][1:]:
            for h, lab in zip(*ax.get_legend_handles_labels()):
                if lab not in labels:
                    handles.append(h)
                    labels.append(lab)
        leg_ax.legend(handles, labels, loc="center", ncol=min(2, max(1, len(labels))), frameon=False)
        fig.suptitle(title)
    return fig


# -------------------------------------------------------------------------------------------------------- regions
def regions(rows, run_id: str, refs: dict | None = None, title: str = ""):
    """Solute fraction c per region of one run (results-table rows, quantity c) as dots with error bars on a log axis (bars
    would have no meaningful baseline on a log scale); reference compositions as horizontal lines."""
    from matplotlib.ticker import FuncFormatter

    order = ["C", "F", "S", "T", "B", "O", "interior"]
    rs = {r["region"]: r for r in rows if r["run_id"] == run_id and r["quantity"] == "c" and r["region"] in order}
    regs = [g for g in order if g in rs and float(rs[g]["value"]) > 0]
    with _ctx():
        fig, ax, _ = _figure()
        x = np.arange(len(regs))
        v = np.array([float(rs[g]["value"]) for g in regs])
        e = np.array([float(rs[g]["error"]) if rs[g]["error"] != "" else np.nan for g in regs])
        ax.errorbar(x, v, yerr=np.where(np.isfinite(e), e, 0), color=PALETTE[0], marker="o", ms=10, ls="none", capsize=5,
                    label="window mean ± SEM")
        for k, (name, val) in enumerate((refs or {}).items()):
            ax.axhline(val, color=PALETTE[1 + k % 5], ls=("--", ":", "-.")[k % 3], lw=2, label=f"{name} = {val:.3g}")
        for xi, vi in zip(x, v):
            ax.annotate(f"{vi:.3g}", (xi, vi), xytext=(10, 0), textcoords="offset points", va="center", fontsize=13)
        ax.set_xticks(x, [REGION_LABEL[g] for g in regs])
        ax.set_xlim(-0.5, len(regs) - 0.3)
        ax.set_yscale("log")
        fmt = FuncFormatter(lambda y, _: f"{y:.3g}")
        ax.yaxis.set_major_formatter(fmt)
        lo, hi = ax.get_ylim()
        # less than a decade on the axis: label the minor ticks too, else only one number would be shown
        ax.yaxis.set_minor_formatter(fmt if hi / lo < 10 else FuncFormatter(lambda y, _: ""))
        ax.tick_params(axis="y", which="minor", labelsize=14)
        ax.set_ylabel("solute fraction c")
        ax.set_title(_wrap(title or f"Solute fraction per region\n{run_id}"), fontsize=14)
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
        ax.set_title(_wrap(title or cbar_label), fontsize=14)
    return fig


# ---------------------------------------------------------------------------------- along-line spectrum (P4)
def alongline_spectrum(q, S_mean, S_sem, line_length: float, n_snapshots: int, title: str = "",
                       n_eff: float | None = None):
    """Mean normalised spectrum S(q) of the core composition along the line, with the binomial-noise level 1 and its
    95 % band for the mean of n_eff independent spectra (exp(1) values: sd 1/sqrt(n_eff)); n_eff = the effective number
    of snapshots of the window test (D2 section 16, from the Stage-3 provenance), else n_snapshots."""
    q = np.asarray(q, float)
    with _ctx():
        fig, ax, _ = _figure()
        ne = float(n_eff) if n_eff is not None and np.isfinite(n_eff) else float(n_snapshots)
        band = 1.96 / np.sqrt(max(ne, 1.0))
        lab = f"{n_snapshots} snapshot{'s' if n_snapshots != 1 else ''}"
        if n_eff is not None and n_snapshots > 1:
            lab += f", effective {ne:.3g}"
        ax.axhspan(max(0.0, 1 - band), 1 + band, color=GRID, alpha=0.6, lw=0, label=f"binomial noise ± 95 % ({lab})")
        ax.axhline(1.0, color=TEXT_2, lw=1.5, ls=":")
        ax.errorbar(q, S_mean, yerr=np.where(np.isfinite(S_sem), S_sem, 0), color=PALETTE[0], marker="o", capsize=3,
                    label="core, mean over the window")
        ax.set_xlabel("mode q (wavelength L$_x$ / q)")
        ax.set_ylabel("S(q) / binomial noise")
        sec = ax.secondary_xaxis("top", functions=(lambda x: line_length / np.maximum(x, 1e-9),
                                                   lambda w: line_length / np.maximum(w, 1e-9)))
        qs = [k for k in (1, 2, 3, 5, 10, 20) if len(q) and k <= q.max()]
        sec.set_xticks([line_length / k for k in qs], [f"{line_length / k:.0f}" for k in qs])
        sec.set_xlabel("wavelength (Å)")
        sec.tick_params(labelsize=14)
        ax.set_title(_wrap(title or "Along-line core composition: spectrum"), fontsize=14)
        ax.legend(loc="best")
    return fig
