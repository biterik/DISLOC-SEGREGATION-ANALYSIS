#!/usr/bin/env python3
"""v0.2 demo: Stage 4 + summary tables + figures from the final demo analysis (M5).

    python regression/v02_demo_summary.py [--in regression/output/v0.2-D2] [--out regression/results/v0.2-demo] [--label D2]

Input: the Stage-2/3 tables written by (definition set D2, the default since 25 Nov; the 24 Nov demo used D1 + proposals,
regression/output/v0.2-final)
    disegreg analyse examples/stage4/manifests/{nbau,nicu,turlo}/*.yaml --out regression/output/v0.2-D2
Output (small, committed): stage4.csv (+ provenance), results.csv is NOT copied (large); summary_*.csv; fig/*.png.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

from disegreg.results import read_table

ROOT = Path(__file__).resolve().parents[1]
DIS = [sys.executable, "-m", "disegreg.cli"]


def run(*args):
    print("$ disegreg", " ".join(map(str, args)))
    subprocess.run(DIS + [str(a) for a in args], check=True)


def write(path, rows, keys=None):
    keys = keys or list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.4g}" if isinstance(v, float) else v) for k, v in r.items()})
    print(f"-> {path} ({len(rows)} rows)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", default=str(ROOT / "regression/output/v0.2-D2"))
    ap.add_argument("--out", default=str(ROOT / "regression/results/v0.2-demo"))
    ap.add_argument("--label", default="D2", help="definition set of the input tables (file names <run>_stage2_<label>.csv)")
    a = ap.parse_args()
    inp, out = Path(a.inp), Path(a.out)
    (out / "fig").mkdir(parents=True, exist_ok=True)
    big = inp / "_tables"
    big.mkdir(exist_ok=True)
    run("table", inp, "--out", big / "results.csv")
    run("across", inp, "--out", out / "stage4.csv", "--sites", ROOT / "examples/stage4/sites.yaml")
    s4 = read_table(out / "stage4.csv")
    for r in s4:
        for k in ("value", "error", "delta", "delta_2"):
            r[k] = float(r[k]) if r[k] != "" else float("nan")
    res = read_table(big / "results.csv")
    for r in res:
        r["value"] = float(r["value"]) if r["value"] != "" else float("nan")
        r["error"] = float(r["error"]) if r["error"] != "" else float("nan")
        r["delta"] = float(r["delta"])

    # ---- uptake exponents
    up = [r for r in s4 if r["analysis"] == "uptake"]
    write(out / "summary_uptake.csv", up, ["project", "potential", "character", "family", "delta", "delta_2", "ref", "value",
                                           "error", "flags", "paper_flag", "runs"])
    # ---- McLean
    mc = [r for r in s4 if r["analysis"] == "mclean" and r["quantity"] in ("c_C_mclean", "c_C_measured_over_mclean")]
    rows = {}
    for r in mc:
        k = (r["potential"], r["character"], r["delta"], r["ref"], r["runs"])
        rows.setdefault(k, {"potential": k[0], "character": k[1], "delta": k[2], "ref": k[3], "flags": r["flags"],
                            "paper_flag": r["paper_flag"]})
        rows[k][r["quantity"]] = r["value"]
        rows[k][r["quantity"] + "_err"] = r["error"]
    cC = {(r["run_id"]): r for r in res if r["quantity"] == "c" and r["region"] == "C"}
    for k, v in rows.items():
        v["c_C_measured"] = cC[k[4]]["value"] if k[4] in cC else float("nan")
    write(out / "summary_mclean.csv", sorted(rows.values(), key=lambda v: (v["potential"], v["character"], v["delta"], v["ref"])))
    # ---- hysteresis
    hy = [r for r in s4 if r["analysis"] == "hysteresis"]
    write(out / "summary_hysteresis.csv", hy, ["potential", "character", "delta", "quantity", "region", "ref", "value", "error",
                                               "status", "runs", "flags"])
    # ---- along line (P4) and domains (P5), per run, window means
    keep = {("alongline_dispersion", "C"), ("alongline_S_max", "C"), ("alongline_q_peak", "C"),
            ("alongline_lambda_peak", "C"), ("alongline_p_peak", "C"), ("alongline_c", "C"),
            ("domain_frac_ordered", "C"), ("domain_frac_ordered", "B"), ("domain_n_all", "all"), ("c", "C")}
    per = {}
    for r in res:
        if (r["quantity"], r["region"]) in keep and r["ref"] == "":
            d = per.setdefault(r["run_id"], {"run_id": r["run_id"], "potential": r["potential"], "character": r["character"],
                                             "delta": r["delta"], "start_type": r["start_type"], "n_samples": r["n_samples"]})
            d[f"{r['quantity']}:{r['region']}"] = r["value"]
    # window-averaged spectrum and its peak test (mean of n independent exp(1) -> Gamma(n, 1/n))
    from scipy.stats import gamma as gamma_dist
    spec = {}
    for r in res:
        if r["quantity"] == "alongline_S" and r["ref"] == "":
            spec.setdefault(r["run_id"], {})[int(r["region"][1:])] = (r["value"], int(r["n_samples"]))
    for rid, d in per.items():
        if rid not in spec:
            continue
        q = sorted(spec[rid])
        S = np.array([spec[rid][k][0] for k in q])
        n = max(1, spec[rid][q[0]][1])
        j = int(np.nanargmax(S)) if np.isfinite(S).any() else -1
        if j < 0:
            continue
        p1 = gamma_dist.sf(S[j], a=n, scale=1.0 / n)          # one mode
        d["S_mean_max"] = float(S[j])
        d["q_mean_max"] = q[j]
        d["p_mean_max"] = float(1 - (1 - p1) ** len(q))        # any of the modes
    write(out / "summary_alongline_domains.csv", sorted(per.values(), key=lambda d: (d["potential"], d["character"], d["delta"])))
    # ---- site tables
    prov = json.loads((out / "stage4.csv.provenance.json").read_text())
    write(out / "summary_site_tables.csv", prov.get("site_tables", []),
          ["potential", "character", "n_sites", "n_core_tabulated", "n_core_expected_per_period", "coverage", "n_periods",
           "x_extent_A", "x_extent_distinct_A", "n_periodic_copies", "period_A", "n_degenerate_sites", "E_min_core", "E_mean_core", "flags"])
    # ---- figures
    fig = out / "fig"
    run("plot", "uptake", out / "stage4.csv", "--out", fig / "uptake_all")
    for pot in ("LeeShim2004-MEAM", "Mendelev2022-EAMFS", "EAMX-DawChandross2023"):
        f = out / f"_stage4_{pot}.csv"
        rows4 = [r for r in read_table(out / "stage4.csv") if r["potential"] == pot]
        if rows4:
            with open(f, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=list(rows4[0]))
                w.writeheader()
                w.writerows(rows4)
            run("plot", "uptake", f, "--out", fig / f"uptake_{pot}")
            f.unlink()
    run("plot", "delta", big / "results.csv", "--quantity", "c:C", "--log", "--out", fig / "cC_vs_delta")
    run("plot", "delta", big / "results.csv", "--quantity", "Gamma:20:r:x_iso", "--out", fig / "Gamma20_xiso_vs_delta")
    for rid in ("mdmc300K_d60_LeeShim2004-MEAM_0.10kT-below-coex_dmu0.8856", "mdmc300K_d60_CuNi-ACE_0.10kT-below-coex_dmu1.1373",
                "mdmc300K_d0_Mendelev2022-EAMFS_0.50kT-below-coex_dmu-3.0234",
                "t51_U20_d90-SGC-dmu-3.5598-fromCu_Mendelev2009-EAMFS_800K",
                "t51_D20_d90-SGC-dmu-3.5598-fromR3_Mendelev2009-EAMFS_800K"):
        s2 = inp / f"{rid}_stage2_{a.label}.csv"
        if s2.is_file():
            run("plot", "alongline", s2, "--out", fig / f"alongline_{rid[:40]}")
            run("plot", "timeseries", s2, "--quantity", "c:C", "--out", fig / f"timeseries_cC_{rid[:40]}")
            run("plot", "gamma", inp / f"{rid}_gammaR_{a.label}.csv", "--out", fig / f"gammaR_{rid[:40]}")


if __name__ == "__main__":
    main()
