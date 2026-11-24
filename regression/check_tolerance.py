#!/usr/bin/env python3
"""Regression with a tolerance: does the OLD recipe on disegreg's Stage-1 data reproduce the OLD project numbers?

DXA line positions are reproducible only to a few tenths of an A between machines (DEFINITIONS 2: cmmg vs M5, same
OVITO 3.17.2: one partial displaced by 0.09 A on average, 0.5 A at most). The old tables were computed on cmmg, the
old recipe here on the M5, so DXA-derived numbers are compared with a tolerance instead of exactly:

    positions (partial separation)      |old - old_recipe| <= 0.5 A
    concentrations, solute counts        <= max(2 % of the value, one solute atom)   [one atom: 1 / N_bin for a
                                                                                       concentration, 1 for a count]
    order parameters (alpha1, R2)        <= max(2 % of the value, 0.01)

plus half a unit of the last digit printed in the old table (rounding). Verdicts per value:
    exact      within the rounding of the old table
    tolerance  outside the rounding, inside the tolerance
    OUTSIDE    outside the tolerance -> look at it

    python check_tolerance.py [--results DIR]       (reads nicu_comparison.csv, nbau_comparison.csv; writes
                                                      tolerance_check.csv next to them)
These tolerances belong to the regression test only; they are not part of a definition set.
"""

from __future__ import annotations

import argparse
import csv
from decimal import Decimal
from pathlib import Path

POSITION_TOL = 0.5       # A
REL_TOL = 0.02           # 2 %
ORDER_ABS_TOL = 0.01

# project: (key columns, [(quantity, kind, column holding the atoms per bin for the one-atom floor or None)])
CHECKS = {
    "nicu": (("ch", "pot", "off", "t_ps"),
             [("c0_5", "conc", "N0_5_old_recipe"), ("c5_10", "conc", "N5_10_old_recipe"),
              ("c10_20", "conc", "N10_20_old_recipe"), ("c_far", "conc", "Nfar_old_recipe"), ("sep", "position", None)]),
    "nbau": (("pot", "ch", "off"),
             [("N5", "count", None), ("c_core", "conc", None), ("c_shell", "conc", None), ("c_SF", "conc", None),
              ("c_far", "conc", None), ("a1_core", "order", None), ("a1_far", "order", None), ("r2_core", "order", None),
              ("r2_far", "order", None), ("sep", "position", None)]),
}


def rounding(text: str) -> float:
    """Half a unit of the last printed digit, e.g. '0.1259' -> 5e-5, '19' -> 0.5, '1.2e-05' -> 5e-8."""
    exp = Decimal(text.strip()).as_tuple().exponent
    return 0.5 * 10.0 ** exp


def tolerance(kind: str, value: float, n_bin: float | None) -> float:
    if kind == "position":
        return POSITION_TOL
    if kind == "order":
        return max(REL_TOL * abs(value), ORDER_ABS_TOL)
    one_atom = 1.0 if kind == "count" else (1.0 / n_bin if n_bin else 0.0)
    return max(REL_TOL * abs(value), one_atom)


def check(results: Path) -> list[dict]:
    out = []
    for project, (keys, quantities) in CHECKS.items():
        f = results / f"{project}_comparison.csv"
        if not f.is_file():
            print(f"missing {f}")
            continue
        for row in csv.DictReader(open(f)):
            if row.get("status", "ok") != "ok":
                continue
            for q, kind, ncol in quantities:
                o_txt, r_txt = row.get(f"{q}_old", ""), row.get(f"{q}_old_recipe", "")
                if not o_txt or not r_txt or "nan" in (o_txt.lower(), r_txt.lower()):
                    continue
                o, r = float(o_txt), float(r_txt)
                n_bin = float(row[ncol]) if ncol and row.get(ncol) else None
                rnd = rounding(o_txt)
                tol = tolerance(kind, o, n_bin) + rnd
                d = abs(o - r)
                verdict = "exact" if d <= rnd * (1 + 1e-9) else ("tolerance" if d <= tol else "OUTSIDE")
                out.append({"project": project, "run": " ".join(row[k] for k in keys), "quantity": q, "kind": kind,
                            "old": o_txt, "old_recipe": f"{r:.6g}", "diff": f"{d:.3g}",
                            "rel_diff": f"{d / abs(o):.3g}" if o else "", "tolerance": f"{tol:.3g}", "verdict": verdict})
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", default=str(Path(__file__).parent / "results"))
    args = ap.parse_args()
    results = Path(args.results)
    rows = check(results)
    if not rows:
        return 1
    with open(results / "tolerance_check.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    for project in CHECKS:
        sel = [r for r in rows if r["project"] == project]
        runs = sorted({r["run"] for r in sel})
        worst = {run: min(("OUTSIDE", "tolerance", "exact").index(r["verdict"]) for r in sel if r["run"] == run)
                 for run in runs}
        n = {v: sum(1 for x in worst.values() if x == i) for i, v in enumerate(("OUTSIDE", "tolerance", "exact"))}
        print(f"{project}: {len(runs)} runs -- all values exact: {n['exact']}, within tolerance: {n['tolerance']}, "
              f"OUTSIDE: {n['OUTSIDE']}")
        tol_runs = sorted(run for run, v in worst.items() if v == 1)
        if tol_runs and len(tol_runs) <= 5:
            print(f"  within tolerance: {'; '.join(tol_runs)}")
        close = [r for r in sel if r["verdict"] != "exact"]
        if close:
            r = max(close, key=lambda r: float(r["diff"]) / float(r["tolerance"]))
            print(f"  closest call: {r['run']} {r['quantity']} old {r['old']} recipe {r['old_recipe']} diff {r['diff']}"
                  f" = {float(r['diff']) / float(r['tolerance']):.0%} of the tolerance {r['tolerance']}")
        for r in sel:
            if r["verdict"] == "OUTSIDE":
                print(f"  OUTSIDE {r['run']:40} {r['quantity']:8} old {r['old']:>10} recipe {r['old_recipe']:>10}"
                      f"  diff {r['diff']} (tol {r['tolerance']})")
    print(f"-> {results / 'tolerance_check.csv'}")
    return 1 if any(r["verdict"] == "OUTSIDE" for r in rows) else 0


if __name__ == "__main__":
    raise SystemExit(main())
