"""Stage 2: everything per snapshot that follows cheaply from the Stage-1 cache (docs/DESIGN.md, Stages).

distances (3.1-3.4) -> glide-plane check (3.4) -> regions (4) -> counts, Gamma, partial separation (6) -> order, clusters (7).

The result is a flat list of records (quantity, region, basis, value), the input of the tidy results table. Counts are
always written, so Gamma can be recomputed for any x_ref later; Gamma, E and alpha vs x_ref only if x_ref is given.
"""

from __future__ import annotations

import time
from collections import Counter
from dataclasses import dataclass, field

import numpy as np

from . import composition, order, regions
from .cache import STRUCTURE_CODES, Snapshot
from .frame import SHOCKLEY
from .geometry import line_coverage, line_distances


@dataclass
class Stage2Result:
    records: list = field(default_factory=list)   # dicts: quantity, region, basis, value
    info: dict = field(default_factory=dict)      # scalars and flags (mode, coverage, z_gp check, timings)
    # per-atom arrays, kept for plots and tests (not written to the tables)
    r: np.ndarray | None = None
    r_dxa: np.ndarray | None = None
    codes: np.ndarray | None = None
    cluster_size: np.ndarray | None = None

    def add(self, quantity, value, region="", basis=""):
        v = float(value) if value is not None else float("nan")
        self.records.append({"quantity": quantity, "region": region, "basis": basis, "value": v})

    def value(self, quantity, region="", basis=""):
        for rec in self.records:
            if rec["quantity"] == quantity and rec["region"] == region and rec["basis"] == basis:
                return rec["value"]
        raise KeyError((quantity, region, basis))


def analyse_snapshot(snap: Snapshot, definitions, a_cell: float, z_gp=None, x_ref=None, workers: int = -1,
                     keep_arrays: bool = False) -> Stage2Result:
    """All Stage-2 quantities of one snapshot.

    a_cell : lattice parameter the cell was built with (manifest), for the neighbour-shell cut-offs.
    z_gp   : glide-plane height from the manifest; None -> taken from the lines (and flagged).
    x_ref  : reference solute fraction; None -> counts and local quantities only.
    """
    t0 = time.perf_counter()
    res = Stage2Result()
    cell = snap.cell()
    pos = snap.positions.astype(float)
    inp = definitions["input"]
    solute = snap.types == int(inp["type_solute"])
    acna = snap.acna
    dist, reg, comp, ordp = (definitions[k] for k in ("distance", "regions", "composition", "order"))
    fb = dist["fallback_axis"]

    # ---- distances (3.1 - 3.4)
    excluded, z_lo, z_hi = regions.excluded_mask(pos, float(reg["surface_exclusion"]))
    nonfcc = ~excluded & (acna != STRUCTURE_CODES["fcc"])
    core = acna == STRUCTURE_CODES["other"]
    lines = snap.polylines()
    ld = line_distances(cell, pos, lines, nonfcc, core, a_cell, spacing=float(dist["polyline_spacing"]),
                        n_images=int(dist["periodic_images"]), slab_width=float(fb["slab_width"]),
                        smoothing_slabs=int(fb["smoothing_slabs"]), connectivity_cutoff_a=float(fb["connectivity_cutoff_a"]),
                        min_slab_fraction=float(fb["min_slab_fraction"]), dxa_assign_cutoff_a=float(fb["dxa_assign_cutoff_a"]),
                        active_below_dxa_coverage=float(fb["active_below_dxa_coverage"]), workers=workers)
    shockley = snap.polylines(families={SHOCKLEY})
    f_shockley = line_coverage(cell, shockley, float(dist["polyline_spacing"]))
    fam = Counter(snap.seg_family.tolist())
    res.info.update(mode=ld.mode, f_cov=ld.f_cov, f_shockley=f_shockley, n_axes=len(ld.axes), n_clusters=ld.n_clusters,
                    n_slabs=ld.n_slabs, f_cov_axis=[float(f) for f in ld.f_cov_axis],
                    n_active_slabs=[int(a.sum()) for a in ld.active], segments=dict(fam),
                    non_shockley_segments=int(sum(v for k, v in fam.items() if k != SHOCKLEY)),
                    length_per_family={f: float(snap.seg_length[snap.seg_family == f].sum()) for f in fam},
                    z_lo=z_lo, z_hi=z_hi)
    for q in ("f_cov", "f_shockley", "n_axes", "n_clusters"):
        res.add(q, res.info[q])
    res.add("n_active_slabs", sum(res.info["n_active_slabs"]))
    for f, length in res.info["length_per_family"].items():
        res.add("dxa_length", length, region=f)

    # ---- glide plane (3.4) and regions (4)
    z_lines, src = regions.z_gp_from_lines(cell, shockley, [ax.polyline(cell) for ax in ld.axes],
                                           float(dist["polyline_spacing"]))
    gp = regions.glide_plane(z_gp, z_lines, src, a_cell, float(dist["z_gp_tolerance_planes"]))
    res.info.update(z_gp=gp.z_gp, z_gp_source=gp.source, z_gp_lines=gp.z_lines, z_gp_lines_source=gp.lines_source,
                    z_gp_difference=gp.difference, z_gp_flag=gp.flag)
    for q in ("z_gp", "z_gp_lines", "z_gp_difference"):
        res.add(q, res.info[q])
    res.add("z_gp_flag", gp.flag)
    if not np.isfinite(gp.z_gp):
        raise ValueError("no glide-plane height: no z_gp given and neither Shockley lines nor fallback axes found")
    codes = regions.assign_regions(pos, acna, ld.r, excluded, gp.z_gp, reg)
    slab40 = regions.far_slab40_mask(pos, acna, excluded, gp.z_gp, float(reg["far_slab40_half_width"]))

    # ---- composition (6)
    radii = list(comp["headline_radii"]) + list(reg["radial_bins"][1:])
    cnt = composition.count(solute, codes, ld.r, ld.r_dxa, reg["radial_bins"], radii, cell.line_length, slab40)
    L = cnt.L_x
    res.add("L_x", L)
    for k, code in enumerate(cnt.regions):
        res.add("N", cnt.n_region[k], region=code)
        res.add("N_sol", cnt.nsol_region[k], region=code)
        res.add("c", composition.concentration(cnt.nsol_region[k], cnt.n_region[k]), region=code)
    res.add("N", cnt.n_far_slab40, region="far_slab40")
    res.add("N_sol", cnt.nsol_far_slab40, region="far_slab40")
    res.add("c", composition.concentration(cnt.nsol_far_slab40, cnt.n_far_slab40), region="far_slab40")
    for basis in ("r", "r_dxa"):
        for lab, n, ns in zip(cnt.bin_labels, cnt.n_bin[basis], cnt.nsol_bin[basis]):
            res.add("N_bin", n, region=lab, basis=basis)
            res.add("N_sol_bin", ns, region=lab, basis=basis)
            res.add("c_bin", composition.concentration(ns, n), region=lab, basis=basis)
        for R, n, ns in zip(cnt.radii, cnt.n_inside[basis], cnt.nsol_inside[basis]):
            res.add("N_inside", n, region=f"{R:g}", basis=basis)
            res.add("N_sol_inside", ns, region=f"{R:g}", basis=basis)
            if x_ref is not None:
                res.add("Gamma", composition.excess(ns, n, x_ref, L), region=f"{R:g}", basis=basis)
    if x_ref is not None:
        res.add("x_ref", x_ref)
        for k, code in enumerate(cnt.regions):
            if code == "X":
                continue
            res.add("Gamma_region", composition.excess(cnt.nsol_region[k], cnt.n_region[k], x_ref, L), region=code)
            res.add("E", composition.concentration(cnt.nsol_region[k], cnt.n_region[k]) / x_ref, region=code)
    ps = composition.partial_separation(cell, shockley, float(comp["partial_separation_bin"]),
                                        float(dist["polyline_spacing"]), float(comp["partial_merge_distance"]))
    res.add("partial_separation_mean", ps.mean)
    res.add("partial_separation_std", ps.std)
    res.add("partial_separation_valid_bins", ps.n_valid)
    res.info.update(partial_separation_bins=ps.n_bins, n_shockley_polylines=ps.n_partials)

    # ---- order and clusters (7)
    sh = order.shell_pairs(cell, pos, a_cell, float(ordp["shell1_cutoff_a"]), float(ordp["shell2_cutoff_a"]))
    nc = order.neighbour_counts(len(pos), solute, sh)
    ob = order.order_by_region(solute, codes, nc, x_ref)
    size = order.solute_cluster_sizes(len(pos), solute, sh)
    cb = order.clusters_by_region(size, solute, codes)
    for k, code in enumerate(ob.regions):
        if code == "X":
            continue
        for q in ("z1", "z2", "p1", "p2", "alpha1_local", "alpha2_local", "R2"):
            res.add(q, getattr(ob, q)[k], region=code)
        if x_ref is not None:
            res.add("alpha1_xref", ob.alpha1_xref[k], region=code)
            res.add("alpha2_xref", ob.alpha2_xref[k], region=code)
        res.add("cluster_mean_size", cb.mean_size[k], region=code)
        res.add("cluster_frac_isolated", cb.frac_isolated[k], region=code)
        res.add("cluster_max_size", cb.max_size[k], region=code)
    res.info["timing_s"] = round(time.perf_counter() - t0, 3)
    if keep_arrays:
        res.r, res.r_dxa, res.codes, res.cluster_size = ld.r, ld.r_dxa, codes, size
    return res
