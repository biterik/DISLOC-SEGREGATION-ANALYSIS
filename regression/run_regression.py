#!/usr/bin/env python3
"""Regression of disegreg (definition set D1) against the existing project analyses (docs/DESIGN.md, Testing). M5 only.

For every reference number three values are written:
  old        the number in the project's table (read, never recomputed by the project script)
  old_recipe the OLD definition re-implemented on disegreg's Stage-1 data (same snapshots): reproduces 'old' if the
             inputs (dumps, DXA lines, labels) agree -- this separates input differences from definition changes
  D1         the disegreg value under definition set D1

    python run_regression.py nbau | nicu | turlo  [--out DIR]

Read-only on ~/Desktop/SIMULATIONS (except the Stage-1 cache in DISEGREG-CACHE). Output: CSV per project in --out.
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import re
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from disegreg import order
from disegreg.cache import STRUCTURE_CODES, get_snapshot
from disegreg.definitions import load_definitions
from disegreg.regions import REGION_CODES
from disegreg.stage2 import analyse_snapshot

SIM = Path("~/Desktop/SIMULATIONS").expanduser()
D1 = load_definitions("D1")
FRAMES = {"d0": ([0, -1, 1], [2, -1, -1], [1, 1, 1]), "d60": ([0, -1, 1], [2, -1, -1], [1, 1, 1]),
          "d30": ([1, -2, 1], [1, 0, -1], [1, 1, 1]), "d90": ([1, -2, 1], [1, 0, -1], [1, 1, 1])}
FCC, HCP = STRUCTURE_CODES["fcc"], STRUCTURE_CODES["hcp"]
IDX = {c: k for k, c in enumerate(REGION_CODES)}


def frame_matrix(ch):
    F = np.array(FRAMES[ch], float)
    return F / np.linalg.norm(F, axis=1)[:, None]


def rec(res, q, region="", basis=""):
    try:
        return res.value(q, region, basis)
    except KeyError:
        return float("nan")


def mean_sem(v):
    v = np.asarray([x for x in v if x is not None and np.isfinite(x)], float)
    if len(v) > 1:
        return float(v.mean()), float(v.std(ddof=1) / np.sqrt(len(v)))
    return (float(v[0]), float("nan")) if len(v) else (float("nan"), float("nan"))


def d1_snapshot(dump, a_cell, ch, x_ref, z_gp=None):
    snap = get_snapshot(dump, D1, verbose=True)
    res = analyse_snapshot(snap, D1, a_cell, z_gp=z_gp, x_ref=x_ref, frame_matrix=frame_matrix(ch) if ch else None,
                           keep_arrays=True)
    return snap, res


def write(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})
    print(f"-> {path} ({len(rows)} rows)")


# ------------------------------------------------------------------------------------------------------------- Nb-Au
NBAU = SIM / "Nb-Au-AT-DISLOCS-Ni-ACE-EAM" / "07_MDMC-SOLUTE-SEGREGATION-T300-WHOLE-INTERIOR"
NBAU_TABLE = NBAU / "analysis" / "harvest_2026-11-21" / "time-avg_07_2026-11-21.dat"
NBAU_XB = {("EAMX-DawChandross2023", "0.25"): 1.6204e-05, ("EAMX-DawChandross2023", "0.50"): 1.2650e-05,
           ("EAMX-DawChandross2023", "1.00"): 7.6509e-06, ("Mendelev2022-EAMFS", "0.25"): 6.8084e-02,
           ("Mendelev2022-EAMFS", "0.50"): 6.4393e-02, ("Mendelev2022-EAMFS", "1.00"): 5.7479e-02}  # analyse_07_harvest.py
NBAU_T0 = {"EAMX-DawChandross2023": 1.45, "Mendelev2022-EAMFS": 0.75}  # ns, start of the old averaging window


def nbau_a_cell():
    out = {}
    for line in open(NBAU / "mdmc-parameters.dat"):
        p = line.split()
        if p and not line.startswith("#"):
            out[(p[0], f"{float(p[2]):.2f}")] = float(p[4])  # potential solute offset_kT dmu_eV a_xT_A ...
    return out


def nbau_old_recipe(snap, a_cell):
    """region_order.py + dxa_snapshot.py of 2026-11-21: orthogonal-box periodic distance to DXA points (segments >= 5 A,
    0.8 A spacing), regions core / shell (incl. hcp) / SF / far (fcc, |z - zc| > 40), alpha vs c_far, 2nd shell 1.10 a."""
    pos = snap.positions.astype(float)
    cell = snap.cell()
    lx, ly = cell.H[0, 0], cell.H[1, 1]
    lines = [p for p, L in zip(snap.polylines(), snap.seg_length) if len(p) >= 2 and L >= 5]
    zc = float(np.median([p[:, 2].mean() for p in lines]))
    zlo, zhi = pos[:, 2].min() + 15, pos[:, 2].max() - 15
    inner = (pos[:, 2] > zlo) & (pos[:, 2] < zhi)
    P = []
    for q in lines:
        for A, B in zip(q[:-1], q[1:]):
            n = max(2, int(np.linalg.norm(B - A) / 0.8))
            P += list(A + (B - A) * np.linspace(0, 1, n)[:, None])
    P = np.asarray(P)
    P[:, 0] %= lx
    P[:, 1] %= ly
    w = pos.copy()
    w[:, 0] %= lx
    w[:, 1] %= ly
    tree = cKDTree(np.column_stack([P[:, 0], P[:, 1], P[:, 2] - zlo + 1]), boxsize=[lx, ly, 1e6])
    dist, _ = tree.query(np.column_stack([w[:, 0], w[:, 1], w[:, 2] - zlo + 1]))
    st = snap.acna
    reg = {"core": inner & (dist < 5), "shell": inner & (dist >= 5) & (dist < 10),
           "SF": inner & (st == HCP) & (np.abs(pos[:, 2] - zc) < 4) & (dist >= 5),
           "far": inner & (st == FCC) & (np.abs(pos[:, 2] - zc) > 40)}
    sol = snap.types == 2
    c_far = sol[reg["far"]].sum() / max(reg["far"].sum(), 1)
    sh = order.shell_pairs(cell, pos, a_cell, 0.854, 1.10)
    n1, s1, n2, s2 = order.neighbour_counts(len(pos), sol, sh)
    out = {"zc": zc}
    for name, m in reg.items():
        idx = m & sol
        out[f"c_{name}"] = float(idx.sum() / max(m.sum(), 1))
        p1 = s1[idx].sum() / max(n1[idx].sum(), 1)
        out[f"a1_{name}"] = float((p1 - c_far) / (1 - c_far)) if idx.any() else float("nan")
        out[f"r2_{name}"] = float(s2[idx].mean() / (6 * out[f"c_{name}"])) if idx.any() else float("nan")
    out["N5"] = int((reg["core"] & sol).sum())
    big = [p for p, b, L in zip(snap.polylines(), snap.seg_spatial_b, snap.seg_length)
           if len(p) >= 2 and L >= 5 and np.linalg.norm(b) > 0.5]
    if len(big) == 2:
        dy = abs(big[1][:, 1].mean() - big[0][:, 1].mean())
        out["sep"] = float(min(dy, ly - dy))
    return out


def nbau(out_dir: Path):
    old = {}
    with open(NBAU_TABLE) as f:
        hdr = None
        for line in f:
            if line.startswith("# pot "):
                hdr = line[2:].split()
            elif not line.startswith("#") and line.strip():
                v = dict(zip(hdr, line.split()))
                old[(v["pot"], v["ch"], v["off"])] = v
    acell = nbau_a_cell()
    rx = re.compile(r"_quenched_step(\d+)\.dump$")
    rows, per_snap = [], []
    for (pot, ch, off), ov in sorted(old.items()):
        runs = [d for d in glob.glob(str(NBAU / f"mdmc300K_{ch}_{pot}_{off}kT-below-coex_dmu*"))
                if Path(d).is_dir() and re.search(r"dmu[-\d.]+$", d)]
        if len(runs) != 1:
            print("skip", pot, ch, off, runs)
            continue
        dumps = sorted((int(rx.search(p).group(1)), p) for p in glob.glob(runs[0] + "/*_quenched_step*.dump"))
        dumps = [(s, p) for s, p in dumps if s / 1e6 >= NBAU_T0[pot] - 1e-9]
        a = acell[(pot, off)]
        xb = NBAU_XB[(pot, off)]
        S = []
        for step, p in dumps:
            snap, res = d1_snapshot(p, a, ch, xb)
            o = nbau_old_recipe(snap, a)
            c_far40 = rec(res, "c", "far_slab40")
            p1 = {k: rec(res, "p1", k) for k in "CFSB"}
            s = {"step": step, "N5_old_recipe": o["N5"], "c_core_old_recipe": o["c_core"], "c_shell_old_recipe": o["c_shell"],
                 "c_SF_old_recipe": o["c_SF"], "c_far_old_recipe": o["c_far"], "a1_core_old_recipe": o["a1_core"],
                 "a1_far_old_recipe": o["a1_far"], "r2_core_old_recipe": o["r2_core"], "r2_far_old_recipe": o["r2_far"],
                 "sep_old_recipe": o.get("sep", float("nan")),
                 "N5_D1": rec(res, "N_sol", "C"), "c_core_D1": rec(res, "c", "C"), "c_shell_D1": rec(res, "c", "S"),
                 "c_SF_D1": rec(res, "c", "F"), "c_far40_D1": c_far40, "c_B_D1": rec(res, "c", "B"),
                 "a1_core_vs_cfar40_D1": float((p1["C"] - c_far40) / (1 - c_far40)),
                 "a1_far_vs_cfar40_D1": float((p1["B"] - c_far40) / (1 - c_far40)),
                 "alpha1_local_C_D1": rec(res, "alpha1_local", "C"), "alpha1_local_B_D1": rec(res, "alpha1_local", "B"),
                 "r2_core_D1": rec(res, "R2", "C"), "r2_far_D1": rec(res, "R2", "B"),
                 "sep_D1": rec(res, "partial_separation_mean"), "Gamma20_D1": rec(res, "Gamma", "20", "r"),
                 "Gamma40_D1": rec(res, "Gamma", "40", "r"), "E_C_D1": rec(res, "E", "C"), "mode": res.info["mode"],
                 "f_cov": res.info["f_cov"], "n_segments": sum(res.info["segments"].values()),
                 "b_total": res.info.get("b_total_label", ""), "z_gp_lines": res.info["z_gp_lines"], "zc_old": o["zc"]}
            S.append(s)
            per_snap.append({"pot": pot, "ch": ch, "off": off, **s})
        r = {"pot": pot, "el": ov["el"], "ch": ch, "off": off, "n_snap_old": int(ov["n_snap"]), "n_snap": len(S)}
        pairs = [("N5", "N5", "N5_D1"), ("c_core", "c_core", "c_core_D1"), ("c_shell", "c_shell", "c_shell_D1"),
                 ("c_SF", "c_SF", "c_SF_D1"), ("c_far", "c_far", "c_far40_D1"), ("a1_core", "a1_core", "a1_core_vs_cfar40_D1"),
                 ("a1_far", "a1_far", "a1_far_vs_cfar40_D1"), ("r2_core", "r2_core", "r2_core_D1"),
                 ("r2_far", "r2_far", "r2_far_D1"), ("sep", "sep", "sep_D1")]
        for name, old_key, d1_key in pairs:
            r[f"{name}_old"] = float(ov[old_key])
            r[f"{name}_old_recipe"] = mean_sem([s[f"{name}_old_recipe"] for s in S])[0]
            r[f"{name}_D1"], r[f"{name}_D1_sem"] = mean_sem([s[d1_key] for s in S])
        for k in ("c_B_D1", "Gamma20_D1", "Gamma40_D1", "E_C_D1", "alpha1_local_C_D1", "alpha1_local_B_D1"):
            r[k], r[k + "_sem"] = mean_sem([s[k] for s in S])
        r["enr5_old"] = float(ov["enr5"])
        r["modes"] = "/".join(sorted({s["mode"] for s in S}))
        r["f_cov_min"] = min(s["f_cov"] for s in S) if S else float("nan")
        r["b_total"] = "/".join(sorted({s["b_total"] for s in S}))
        rows.append(r)
    write(out_dir / "nbau_comparison.csv", rows)
    write(out_dir / "nbau_per_snapshot.csv", per_snap)


# ------------------------------------------------------------------------------------------------------------- Ni-Cu
NICU = SIM / "Cu-AT-DISLOCS-Ni-EAM-MEAM-ACE" / "07_MDMC-CU-SEGREGATION-T300-COEX-OFFSETS"
NICU_TABLE = NICU / "analysis-2026-11-21_dxa-latest" / "dxa-cu-summary_2026-11-21.dat"


def nicu_params():
    out = {}
    for line in open(NICU / "mdmc-parameters.dat"):
        p = line.split()
        if p and not line.startswith("#"):
            out[(p[0], f"{float(p[1]):.2f}")] = (float(p[3]), float(p[4]))  # a_xT, xCu_overall
    return out


def nicu_old_recipe(snap):
    """dxa-cu-near-partials.py: distance to the MEAN (y, z) of each DXA line, periodic in y only; 15 A z margin;
    bins 0-5-10-20-40-80; far >= 80 A (all structures)."""
    pos = snap.positions.astype(float)
    ly = snap.cell().H[1, 1]
    zlo, zhi = pos[:, 2].min() + 15, pos[:, 2].max() - 15
    inner = (pos[:, 2] > zlo) & (pos[:, 2] < zhi)
    dist = np.full(len(pos), 1e9)
    for p in snap.polylines():
        dy = pos[:, 1] - p[:, 1].mean()
        dy -= ly * np.round(dy / ly)
        dist = np.minimum(dist, np.hypot(dy, pos[:, 2] - p[:, 2].mean()))
    cu = snap.types == 2
    out = {}
    for lo, hi in ((0, 5), (5, 10), (10, 20)):
        m = inner & (dist >= lo) & (dist < hi)
        out[f"c{lo}_{hi}"] = cu[m].sum() / max(m.sum(), 1)
        out[f"n{lo}_{hi}"] = (int(cu[m].sum()), int(m.sum()))
    far = inner & (dist >= 80)
    out["c_far"] = cu[far].sum() / max(far.sum(), 1)
    out["n_far"] = int(far.sum())
    ys = sorted(p[:, 1].mean() for p in snap.polylines())
    if len(ys) == 2:
        out["sep"] = min(ys[1] - ys[0], ly - (ys[1] - ys[0]))
    return out


def nicu(out_dir: Path):
    par = nicu_params()
    rows = []
    for line in open(NICU_TABLE):
        if line.startswith("#") or not line.strip():
            continue
        p = line.split()
        ch, pot, off, t_ps = p[0], p[1], p[2], int(p[3])
        c05, c510, c1020, cfar, enr, err, sep = (float(x) for x in (p[5], p[6], p[7], p[8], p[9], p[11], p[12]))
        dumps = glob.glob(str(NICU / f"mdmc300K_{ch}_{pot}_{off}kT-below-coex_dmu*" / f"*_quenched_step{t_ps * 1000}.dump"))
        if not dumps:
            print("not on the M5:", ch, pot, off, t_ps)
            rows.append({"ch": ch, "pot": pot, "off": off, "t_ps": t_ps, "status": "dump not on the M5"})
            continue
        a, x_iso = par.get((pot, off), (3.52, float("nan")))
        snap, res = d1_snapshot(dumps[0], a, ch, x_iso)
        o = nicu_old_recipe(snap)
        cB = rec(res, "c", "B")
        rows.append({
            "ch": ch, "pot": pot, "off": off, "t_ps": t_ps, "status": "ok",
            "c0_5_old": c05 / 100, "c0_5_old_recipe": o["c0_5"], "c0_5_D1": rec(res, "c_bin", "0-5", "r"),
            "c5_10_old": c510 / 100, "c5_10_old_recipe": o["c5_10"], "c5_10_D1": rec(res, "c_bin", "5-10", "r"),
            "c10_20_old": c1020 / 100, "c10_20_old_recipe": o["c10_20"], "c10_20_D1": rec(res, "c_bin", "10-20", "r"),
            "c_far_old": cfar / 100, "c_far_old_recipe": o["c_far"], "c_B_D1": cB, "x_iso": x_iso,
            "enrich_old": enr, "enrich_old_err": err, "enrich_old_recipe": o["c0_5"] / o["c_far"] if o["c_far"] else float("nan"),
            "enrich_C_vs_cB_D1": rec(res, "c", "C") / cB if cB else float("nan"), "E_C_vs_xiso_D1": rec(res, "E", "C"),
            "c_C_D1": rec(res, "c", "C"), "c_S_D1": rec(res, "c", "S"), "c_F_D1": rec(res, "c", "F"),
            "Gamma20_D1": rec(res, "Gamma", "20", "r"), "Gamma40_D1": rec(res, "Gamma", "40", "r"),
            "sep_old": sep, "sep_old_recipe": o.get("sep", float("nan")), "sep_D1": rec(res, "partial_separation_mean"),
            "mode": res.info["mode"], "f_cov": res.info["f_cov"], "b_total": res.info.get("b_total_label", ""),
            "segments": json.dumps(res.info["segments"]),
            # atoms per bin of the old recipe, for the one-solute-atom floor of the tolerance check (check_tolerance.py)
            "N0_5_old_recipe": o["n0_5"][1], "N5_10_old_recipe": o["n5_10"][1], "N10_20_old_recipe": o["n10_20"][1],
            "Nfar_old_recipe": o["n_far"]})
    write(out_dir / "nicu_comparison.csv", rows)


# ------------------------------------------------------------------------------------------------------------- TURLO
TURLO = SIM / "TURLO-REDO" / "50_TURLO-METHOD-VS-OURS-CU-ZR"
TURLO_TABLE = TURLO / "ANALYSIS-FINAL-2NS" / "t50_snapshot-summary.dat"


def turlo(out_dir: Path):
    old = {}
    for line in open(TURLO_TABLE):
        if line.startswith("# run "):
            hdr = line[2:].split()
        elif not line.startswith("#") and line.strip():
            v = line.split()
            old[v[0]] = dict(zip(hdr, v[: len(hdr)]))
    rows = []
    for tag in ("R4", "R3"):
        run = glob.glob(str(TURLO / f"t50_{tag}_d90-*_800K"))[0]
        dumps = sorted(glob.glob(run + "/*_quenched_step*.dump"), key=lambda p: int(re.search(r"step(\d+)", p).group(1)))
        dump = dumps[-1]
        lx = float(open(dump).readlines()[5].split()[1]) - float(open(dump).readlines()[5].split()[0])
        a = lx / (24 * np.sqrt(6) / 2)  # X = [1 -2 1]: 24 repeats of a sqrt(6)/2 (cell length / repeat, see report)
        snap = get_snapshot(dump, D1, verbose=True)
        res0 = analyse_snapshot(snap, D1, a, frame_matrix=frame_matrix("d90"), keep_arrays=True)  # x_ref unknown yet
        cB = rec(res0, "c", "B")
        res = analyse_snapshot(snap, D1, a, frame_matrix=frame_matrix("d90"), x_ref=cB, keep_arrays=True)
        L = rec(res, "L_x")
        ov = old[tag]
        inside = {R: (rec(res, "N_inside", f"{R}", "r"), rec(res, "N_sol_inside", f"{R}", "r")) for R in (10, 20)}
        inside_dxa = {R: (rec(res, "N_inside", f"{R}", "r_dxa"), rec(res, "N_sol_inside", f"{R}", "r_dxa")) for R in (10, 20)}
        rows.append({
            "run": tag, "dump": Path(dump).name, "a_cell_derived": a, "mode": res.info["mode"], "f_cov": res.info["f_cov"],
            "f_shockley": res.info["f_shockley"], "n_axes": res.info["n_axes"],
            "n_active_slabs": sum(res.info["n_active_slabs"]), "n_slabs": res.info["n_slabs"],
            "segments": json.dumps(res.info["segments"]), "b_total": res.info.get("b_total_label", ""),
            "c_far_old": float(ov["c_far"]), "c_B_D1": cB,
            "c_r<10_old": float(ov["c_r<10"]), "c_r<10_D1": inside[10][1] / inside[10][0],
            "c_r<10_D1_rdxa": inside_dxa[10][1] / max(inside_dxa[10][0], 1),
            "c_r<20_old": float(ov["c_r<20"]), "c_r<20_D1": inside[20][1] / inside[20][0],
            "NZr_excess_r<20_old": float(ov["NZr_excess_r<20"]),
            "NZr_excess_r<20_D1": rec(res, "Gamma", "20", "r") * L,
            "NZr_excess_r<20_D1_rdxa": rec(res, "Gamma", "20", "r_dxa") * L,
            "Gamma20_D1": rec(res, "Gamma", "20", "r"), "Gamma20_D1_rdxa": rec(res, "Gamma", "20", "r_dxa"),
            "Gamma40_D1": rec(res, "Gamma", "40", "r"), "Gamma40_D1_rdxa": rec(res, "Gamma", "40", "r_dxa"),
            "c_C_D1": rec(res, "c", "C"), "c_F_D1": rec(res, "c", "F"), "N_C_D1": rec(res, "N", "C"),
            "WC_alpha1_far_old": float(ov["WC_alpha1_far"]), "alpha1_local_B_D1": rec(res, "alpha1_local", "B"),
            "frac_hcp_old": float(ov["frac_hcp"]), "frac_other_old": float(ov["frac_other"]),
            "frac_hcp_acna": float(np.mean(snap.acna == HCP)), "frac_other_acna": float(np.mean(snap.acna == 0)),
            "frac_ico_ptm": float(np.mean(snap.ptm == STRUCTURE_CODES["ico"]))})
    write(out_dir / "turlo_comparison.csv", rows)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("which", choices=("nbau", "nicu", "turlo"))
    ap.add_argument("--out", default=str(Path(__file__).parent / "results"))
    args = ap.parse_args()
    {"nbau": nbau, "nicu": nicu, "turlo": turlo}[args.which](Path(args.out))
