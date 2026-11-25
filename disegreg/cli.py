"""Command-line entry point ``disegreg``. Sub-commands are added stage by stage (see docs/DESIGN.md)."""

import argparse
import sys
from pathlib import Path

from . import __version__
from .definitions import DEFAULT_LABEL, load_definitions


def _cmd_version(_args) -> int:
    d = load_definitions()
    print(f"disegreg {__version__}  definition set {d.label} ({d.status})")
    return 0


def _cmd_stage1(args) -> int:
    from collections import Counter

    from .cache import default_cache_root, get_snapshot

    d = load_definitions(args.definitions)
    root = args.cache or default_cache_root()
    print(f"definition set {d.label} | cache {root}")
    for dump in args.dumps:
        snap = get_snapshot(dump, d, cache_root=root, recompute=args.recompute, verbose=True)
        fam = Counter(snap.seg_family.tolist())
        lengths = {f: round(float(snap.seg_length[snap.seg_family == f].sum()), 1) for f in fam}
        t = snap.meta.get("timings_s", {})
        print(f"  {snap.n_atoms} atoms | {snap.n_segments} DXA segments {dict(fam)} | length per family (A) {lengths}"
              f" | OVITO time {sum(t.values()):.1f} s")
    return 0


def _cmd_manifest(args) -> int:
    from .manifest import load_manifest, template, validate

    if args.action == "template":
        print(template(), end="")
        return 0
    bad = 0
    for f in args.files:
        m = load_manifest(f, strict=False)
        errs = validate(m.data)
        if errs:
            bad += 1
            print(f"INVALID {m.path}")
            for e in errs:
                print(f"  - {e}")
        else:
            n = len(m.dump_files())
            print(f"OK      {m.path}  ({m['run_id']}, Delta = {m['delta']:+.2f}, {n} dump file(s) found)")
    return 1 if bad else 0


def _md_step(dump, fallback) -> int:
    """MD/MC step of a snapshot: the number after 'step' in the file name (the state that was quenched); the dump's own
    TIMESTEP is larger by the minimisation steps."""
    import re
    from pathlib import Path

    m = re.search(r"step(\d+)", Path(dump).name)
    return int(m.group(1)) if m else int(fallback)


def find_x_bc(run_id: str, label: str, folders, min_solute: int = 100):
    """x_bc of a bulk control from its Stage-3 table <run_id>_stage3_<label>.csv (first folder that has it):
    window mean and block-averaging SEM of the interior solute fraction. None if no table is found."""
    import json
    from pathlib import Path

    from .results import read_table
    from .stage2 import BULK_REGION

    for folder in folders:
        f = Path(folder) / f"{run_id}_stage3_{label}.csv"
        if not f.is_file():
            continue
        rows = {(r["quantity"], r["region"]): r for r in read_table(f)}
        c, ns = rows.get(("c", BULK_REGION)), rows.get(("N_sol", BULK_REGION))
        if c is None:
            raise ValueError(f"{f}: no interior solute fraction -- is this a bulk-control table?")
        n_sol = ns["mean"] if ns else float("nan")
        side = f.with_name(f.name + ".provenance.json")
        s3 = (json.loads(side.read_text()).get("stage3") or {}) if side.is_file() else {}
        return {"source": str(f), "value": c["mean"], "sem": c["sem"], "n_samples": c["n_samples"], "flat": c["flat"],
                "n_solute": n_sol, "min_solute": min_solute, "enough_solute": bool(n_sol >= min_solute),
                "n_atoms": s3.get("n_atoms"), "last_md_step": s3.get("last_md_step")}
    return None


def bulk_control_check(x_bc_info: dict, n_atoms_run: int, last_md_step_run, definitions) -> dict:
    """Is the bulk control acceptable as x_bc for this run (DEFINITIONS 5; D2 section 15)?
    Same size or larger: yes. Smaller: D2 (smaller_cell_allowed) yes, but not for protocol diagnostic (iv); D1 no (x_bc is
    then reported as a diagnostic only). Unknown size (table without n_atoms): no. The run length (last MD step) is compared
    and a difference flagged; the per-atom swap rate is not in the data and stays the project's responsibility."""
    from .results import DIAG_IV_NOT_TESTED, DIAG_IV_NOT_USABLE, DIAG_IV_PENDING

    bc = (definitions.params if hasattr(definitions, "params") else definitions).get("bulk_control") or {}
    n_bc = x_bc_info.get("n_atoms")
    out = {"n_atoms_run": int(n_atoms_run), "flags": []}
    if n_bc is None:
        out.update(size_ratio=float("nan"), smaller=None, accepted=False, diagnostic_iv=DIAG_IV_NOT_TESTED)
        out["flags"].append("bulk-control size unknown (no n_atoms in its Stage-3 provenance): x_bc reported as a diagnostic")
    else:
        smaller = int(n_bc) < int(n_atoms_run)
        allowed = bool(bc.get("smaller_cell_allowed", False))
        out.update(size_ratio=int(n_bc) / int(n_atoms_run), smaller=smaller, accepted=(not smaller) or allowed)
        if not smaller:
            out["diagnostic_iv"] = DIAG_IV_PENDING
        else:
            iv = bool(bc.get("smaller_cell_for_diagnostic_iv", False))
            out["diagnostic_iv"] = DIAG_IV_PENDING if iv else DIAG_IV_NOT_USABLE
            if not allowed:
                out["flags"].append(f"bulk control smaller than the dislocation cell ({n_bc} < {n_atoms_run} atoms): "
                                    f"not a bulk control under {getattr(definitions, 'label', '?')} (same size), x_bc "
                                    "reported as a diagnostic")
    lb, lr = x_bc_info.get("last_md_step"), last_md_step_run
    out.update(last_md_step_bulk_control=lb, last_md_step_run=lr,
               same_run_length=(None if lb is None or lr is None else int(lb) == int(lr)))
    if out["same_run_length"] is False:
        out["flags"].append(f"run length differs: bulk control to step {lb}, run to step {lr}")
    return out


def _analyse_bulk(m, d, root, env, tv, args, props=None) -> int:
    """Stage 2 + 3 of a bulk-control run: x_bc per snapshot, then its window mean (DEFINITIONS 5, v0.2)."""
    import time
    from pathlib import Path

    from .cache import get_snapshot
    from .results import STAGE3_COLUMNS, write_table
    from .stage2 import BULK_REGION, analyse_bulk_snapshot
    from .timeseries import read_mc_stream, run_stage3

    dumps = m.dump_files()
    if args.last:
        dumps = dumps[-args.last:]
    if not dumps:
        print(f"{m['run_id']}: no dump files match {m['dumps']}")
        return 1
    rows, snaps = [], []
    for dump in dumps:
        t0 = time.perf_counter()
        snap = get_snapshot(dump, d, cache_root=root, verbose=args.verbose)
        res = analyse_bulk_snapshot(snap, d, float(m["a_cell"]), proposals=props)
        ts = _md_step(dump, snap.meta.get("timestep", -1))
        for r in res.records:
            rows.append({"run_id": m["run_id"], "snapshot": Path(dump).name, "timestep": ts, **r,
                         "definition_set": d.label, "tool_version": tv})
        snaps.append({"dump": str(dump), "source_sha256": snap.meta.get("source_sha256"), "md_step": ts,
                      "n_atoms": int(snap.n_atoms), "info": dict(res.info)})
        if m.get("n_atoms") is not None and int(snap.n_atoms) != int(m["n_atoms"]):
            print(f"  WARNING: {snap.n_atoms} atoms in the snapshot, manifest n_atoms = {m['n_atoms']}")
        print(f"{m['run_id']} (bulk control) step {ts}: x_bc {res.value('c', region=BULK_REGION):.5g} "
              f"({res.info['n_solute_interior']} of {res.info['n_interior']} interior atoms)  "
              f"({time.perf_counter() - t0:.1f} s)")
    out = Path(args.out).expanduser() / f"{m['run_id']}_stage2_{d.label}.csv"
    prov = {"environment": env, "tool_version": tv, "definition_set": d.label, "definitions_sha256": d.sha256,
            "definitions_status": d.status, "manifest": str(m.path), "manifest_data": m.data, "cell_type": "bulk_control",
            "snapshots": snaps}
    write_table(out, rows, prov)
    stream = None
    if m.get("mc_stream") and m.resolve(m["mc_stream"]).is_file():
        stream = read_mc_stream(m.resolve(m["mc_stream"]))
    try:
        info3, rows3 = run_stage3(rows, d["timeseries"], stream=stream, stream_column=m.get("mc_stream_column"),
                                  md_timestep_ps=m.get("md_timestep_ps"), criterion=args.drift_criterion,
                                  key_quantity=("N_sol", BULK_REGION, "", ""))
    except ValueError as e:
        print(f"  Stage 3 skipped: {e}")
        return 0
    for r in rows3:
        r.update(run_id=m["run_id"], definition_set=d.label, tool_version=tv)
    info3.update(n_atoms=snaps[-1]["n_atoms"], n_atoms_manifest=m.get("n_atoms"),
                 last_md_step=max(s["md_step"] for s in snaps))   # for the size / run-length check (D2 section 15)
    out3 = Path(args.out).expanduser() / f"{m['run_id']}_stage3_{d.label}.csv"
    write_table(out3, rows3, {**prov, "stage3": info3}, columns=STAGE3_COLUMNS)
    c = next(r for r in rows3 if r["quantity"] == "c" and r["region"] == BULK_REGION)
    print(f"  Stage 3: x_bc = {c['mean']:.5g} +- {c['sem']:.2g} ({c['n_samples']} snapshots, "
          f"{'converged' if info3['converged'] else 'NOT converged'})  -> {out3}")
    return 0


def _cmd_analyse(args) -> int:
    import time
    from pathlib import Path

    from .cache import default_cache_root, get_snapshot
    from .manifest import load_manifest
    from .provenance import environment
    from .results import (GAMMA_CURVE_COLUMNS, STAGE3_COLUMNS, diagnostic_iv, gamma_curve, paper_hold, tool_version,
                          write_table)
    from .stage2 import analyse_snapshot
    from .timeseries import read_mc_stream, run_stage3

    d = load_definitions(args.definitions)
    root = args.cache or default_cache_root()
    env = environment()
    tv = tool_version(env)
    status = 0
    ref_rule = d["reference"]
    comp = d["composition"]
    hr, cr = float(comp["headline_radius"]), float(comp["check_radius"])
    from .definitions import load_proposals

    props = None if args.no_proposals else load_proposals()  # open proposals (P5; status proposal)
    # bulk controls first, so the runs naming one find its x_bc (v0.2: automatic x_bc)
    manifests = sorted((load_manifest(mf) for mf in args.manifests), key=lambda m: not m.is_bulk_control)
    bulk_dirs = [Path(args.out).expanduser()] + [Path(b).expanduser() for b in (args.bulk_dir or [])]
    for m in manifests:
        if m.is_bulk_control:
            status |= _analyse_bulk(m, d, root, env, tv, args, props)
            continue
        # reference compositions (DEFINITIONS 5, rev. 4): c_B of the snapshot always; x_iso if the manifest has it;
        # x_bc from the bulk control's own Stage-3 table (automatic, v0.2) or a value from the command line
        refs = {}
        x_bc_info = None
        if args.x_ref is not None:
            refs["x_bc" if m["x_ref_source"] == "bulk_control" else "given"] = float(args.x_ref)
            if m["x_ref_source"] == "bulk_control":
                x_bc_info = {"source": "command line --x-ref", "value": float(args.x_ref), "sem": float("nan")}
        elif m.get("bulk_control"):
            x_bc_info = find_x_bc(m["bulk_control"], d.label, bulk_dirs,
                                  int(ref_rule.get("min_solute_in_bulk_control", 100)))
            if x_bc_info is None:
                print(f"{m['run_id']}: bulk control {m['bulk_control']} named, but no {m['bulk_control']}_stage3_"
                      f"{d.label}.csv in {[str(b) for b in bulk_dirs]} -- analyse it first (or give --x-ref)")
            else:
                refs["x_bc"] = x_bc_info["value"]
                print(f"{m['run_id']}: x_bc = {x_bc_info['value']:.5g} +- {x_bc_info['sem']:.2g} from {x_bc_info['source']}"
                      + ("" if x_bc_info["enough_solute"] else
                         f"  [only {x_bc_info['n_solute']:.0f} solute atoms < {x_bc_info['min_solute']}: x_iso is the "
                         "reference (DEFINITIONS 5, 2b), x_bc reported as a diagnostic]"))
        if m.get("x_iso") is not None:
            refs["x_iso"] = float(m["x_iso"])
        dumps = m.dump_files()
        if args.last:
            dumps = dumps[-args.last:]
        if not dumps:
            print(f"{m['run_id']}: no dump files match {m['dumps']}")
            status = 1
            continue
        rows, snaps = [], []
        for dump in dumps:
            t0 = time.perf_counter()
            snap = get_snapshot(dump, d, cache_root=root, verbose=args.verbose)
            res = analyse_snapshot(snap, d, float(m["a_cell"]), z_gp=float(m["z_gp"]), x_ref=refs or None,
                                   frame_matrix=m.frame_matrix(), c_B_ref=True, proposals=props)
            ts = _md_step(dump, snap.meta.get("timestep", -1))
            for r in res.records:
                rows.append({"run_id": m["run_id"], "snapshot": Path(dump).name, "timestep": ts, **r,
                             "definition_set": d.label, "tool_version": tv})
            info = dict(res.info)
            snaps.append({"dump": str(dump), "source_sha256": snap.meta.get("source_sha256"), "md_step": ts,
                          "dump_timestep": snap.meta.get("timestep"),
                          "ovito": snap.meta.get("provenance", {}).get("ovito"), "info": info})
            flag = "  z_gp FLAG" if res.info["z_gp_flag"] else ""
            flag += f"  b_total {res.info['b_total_label'] or res.info['b_total_crystal_a']}" \
                    f" ({res.info['b_total_character_angle']:.0f} deg)"
            print(f"{m['run_id']} step {ts}: mode {res.info['mode']}, f_cov {res.info['f_cov']:.3f}, "
                  f"c_C {res.value('c', region='C'):.4f}, c_B {res.value('c', region='B'):.5f}"
                  f"{flag}  ({time.perf_counter() - t0:.1f} s)")
        # is the bulk control acceptable (size, run length; D2 section 15)? Only for a bulk control found as a table
        if x_bc_info is not None and "enough_solute" in x_bc_info:
            x_bc_info.update(bulk_control_check(x_bc_info, snap.n_atoms, max(s["md_step"] for s in snaps), d))
            for fl in x_bc_info["flags"]:
                print(f"  bulk control: {fl}")
            if x_bc_info["smaller"] and x_bc_info["accepted"]:
                print(f"  bulk control: {x_bc_info['n_atoms']} atoms < {snap.n_atoms} in this cell -- x_bc accepted "
                      f"({d.label} section 15), diagnostic (iv) {x_bc_info['diagnostic_iv']}")
        has_bc = (bool(m.get("bulk_control")) and "x_bc" in refs
                  and (x_bc_info is None or (x_bc_info.get("enough_solute", True) and x_bc_info.get("accepted", True))))
        if not has_bc:
            missing = [r for r in ref_rule["without_bulk_control_report_against"] if r != "c_B" and r not in refs]
            if missing:
                print(f"{m['run_id']}: no bulk control and no {missing} in the manifest -- reported against c_B only")
        rule_text = ("x_bc (bulk control), c_B and x_iso as diagnostics" if has_bc else
                     "no bulk control: reported against both c_B and x_iso (DEFINITIONS 5, rev. 4)")
        hold = "" if has_bc else paper_hold({**m.data, "bulk_control": None},
                                            ref_rule.get("paper_hold_without_bulk_control"))
        diag_iv = diagnostic_iv(d, float(m["delta"]), x_bc_info)   # DEFINITIONS 5; D2 sections 15, 17
        out = Path(args.out).expanduser() / f"{m['run_id']}_stage2_{d.label}.csv"
        prov = {"environment": env, "tool_version": tv, "definition_set": d.label, "definitions_sha256": d.sha256,
                "definitions_status": d.status, "manifest": str(m.path), "manifest_data": m.data,
                "proposal_set": props.label if props else None, "proposals_sha256": props.sha256 if props else None,
                "x_ref": refs, "x_ref_source": m["x_ref_source"], "reference_rule": rule_text, "paper_flag": hold,
                "diagnostic_iv": diag_iv, "x_bc": x_bc_info, "snapshots": snaps}
        write_table(out, rows, prov)
        print(f"  -> {out} ({len(rows)} rows)")
        # ---- Stage 3: averaging window (from the MC stream if there is one), window means, drift tests
        stream = None
        if m.get("mc_stream") and m.resolve(m["mc_stream"]).is_file():
            stream = read_mc_stream(m.resolve(m["mc_stream"]))
        elif m.get("mc_stream"):
            print(f"  MC stream not found: {m.resolve(m['mc_stream'])} -> window from the snapshots")
        try:
            info3, rows3 = run_stage3(rows, d["timeseries"], stream=stream, stream_column=m.get("mc_stream_column"),
                                      md_timestep_ps=m.get("md_timestep_ps"), criterion=args.drift_criterion)
        except ValueError as e:
            print(f"  Stage 3 skipped: {e}")
            continue
        for r in rows3:
            r.update(run_id=m["run_id"], definition_set=d.label, tool_version=tv)
        fa = next((r["mean"] for r in rows3 if r["quantity"] == "fallback_active"), 0.0)
        info3.update(headline_radius=hr, check_radius=cr, fallback_active_fraction=fa, check_required=bool(fa > 0),
                     paper_flag=hold, diagnostic_iv=diag_iv)
        out3 = Path(args.out).expanduser() / f"{m['run_id']}_stage3_{d.label}.csv"
        write_table(out3, rows3, {**prov, "stage3": info3}, columns=STAGE3_COLUMNS)
        curve = gamma_curve(rows3, hr, cr)
        outg = Path(args.out).expanduser() / f"{m['run_id']}_gammaR_{d.label}.csv"
        write_table(outg, curve, {**prov, "stage3": info3}, columns=GAMMA_CURVE_COLUMNS)
        if info3.get("snapshots_after_window_signal"):
            print(f"  WARNING: {info3['snapshots_after_window_signal']} of {info3['snapshots_total']} snapshots are later "
                  f"than the end of the {info3['window_source']} signal and are NOT used (is the MC stream copy older than "
                  "the dumps?)")
        verdict = "converged" if info3["converged"] else "NOT converged"
        if not info3["signal_plateau"]:
            verdict += " (block averaging without plateau: error bar is a lower bound, verdict uncertain)"
        print(f"  Stage 3: window {info3['window_t_start']:.0f}-{info3['window_t_end']:.0f} {info3['time_unit']} "
              f"({info3['window_fraction']:.0%} of the run, {info3['window_source']}), {verdict}"
              f" [criterion {info3['criterion']}]  -> {out3}")
        for c in curve:
            if c["role"] == "headline" or (c["role"] == "check" and info3["check_required"]):
                print(f"  Gamma({c['R']:g} A) vs {c['ref']}: {c['mean']:.4g} +- {c['sem']:.2g} atoms/A"
                      + ("  [check: fallback axis active]" if c["role"] == "check" else "  [headline]"))
        print(f"  Gamma(R) curve -> {outg}")
        if hold:
            print(f"  FLAG: {hold} ({m['potential']} {m['host_element']}-{m['solute_element']}, DEFINITIONS 5)")
    return status


def _cmd_table(args) -> int:
    from pathlib import Path

    from .provenance import environment
    from .results import RESULTS_COLUMNS, collect_runs, tool_version, write_table

    files = sorted({str(p) for f in args.stage3 for p in (Path(f).glob("*_stage3_*.csv") if Path(f).is_dir() else [f])})
    if not files:
        print("no Stage-3 tables (*_stage3_*.csv) found")
        return 1
    rows, prov = collect_runs(files)
    env = environment()
    prov.update(environment=env, written_by=tool_version(env))
    write_table(args.out, rows, prov, columns=RESULTS_COLUMNS)
    runs = sorted({r["run_id"] for r in rows})
    print(f"{len(runs)} runs, {len(rows)} rows, definition set {prov['definition_set']} -> {args.out}")
    return 0


def load_site_cores(config_path, runs, definitions, cache_root=None, verbose=True):
    """Core sites + energies for every (project, potential, character) of the runs that the site configuration covers
    (core sites: D2 section 12, ``definitions["mclean"]``). Returns ({key: (E_core, n_missing, info)}, list of summaries)."""
    import numpy as np
    import yaml

    from .cache import get_snapshot
    from .sites import core_sites, read_site_table

    cfg = yaml.safe_load(Path(config_path).expanduser().read_text())
    base = Path(config_path).expanduser().resolve().parent
    mc = definitions["mclean"]
    min_cov = float(mc.get("min_coverage", 0.9))
    wanted = {}
    for r in runs:
        if not r.is_bulk_control:
            wanted.setdefault((r["project"], r["potential"], r["character"]), r)
    out, summary = {}, []
    for entry in cfg.get("site_tables", []):
        pots = entry.get("potentials") or [entry.get("potential")]
        chars = entry.get("characters") or [entry.get("character")]
        for pot in pots:
            for ch in chars:
                key = (entry["project"], pot, ch)
                if key not in wanted or key in out:
                    continue
                run = wanted[key]
                tpath = Path(str(entry["table"]).format(potential=pot, character=ch)).expanduser()
                cpath = Path(str(entry["cell_0K"]).format(potential=pot, character=ch)).expanduser()
                tpath, cpath = (q if q.is_absolute() else base / q for q in (tpath, cpath))
                table = read_site_table(tpath, entry.get("format", "segE0K-by-site"), pot, ch)
                snap0 = get_snapshot(cpath, definitions, cache_root=cache_root, verbose=verbose,
                                     pbc=entry.get("cell_0K_pbc"))
                cs = core_sites(table, snap0, float(mc["core_site_radius"]), run["frame"]["X"], float(run["a_cell"]),
                                spacing=float(definitions["distance"]["polyline_spacing"]),
                                surface_exclusion=float(definitions["regions"]["surface_exclusion"]))
                n_missing = max(0.0, cs.n_core_expected - cs.n_core_tabulated) if mc["missing_sites"] == "bulk" else 0.0
                e_core = table.energy[cs.is_core]
                flags = []
                if cs.coverage < min_cov:
                    flags.append(f"site table covers {cs.coverage:.0%} of the core sites")
                if cs.info["x_extent_distinct"] > cs.info["period"] + 0.5:
                    flags.append(f"distinct sites span {cs.info['x_extent_distinct']:.1f} A along X > one period "
                                 f"{cs.info['period']:.2f} A")
                n_deg = table.degenerate_sites()
                info = {"project": key[0], "potential": pot, "character": ch, "table": str(tpath), "table_sha256": table.sha256,
                        "format": table.format, "cell_0K": str(cpath), "cell_0K_sha256": snap0.meta.get("source_sha256"),
                        "n_sites": len(table), "n_core_tabulated": cs.n_core_tabulated,
                        "n_core_expected_per_period": round(cs.n_core_expected, 2), "coverage": round(cs.coverage, 3),
                        "n_periods": cs.n_periods, "period_A": round(cs.info["period"], 4), "x_extent_A": round(cs.x_extent, 3),
                        "x_extent_distinct_A": round(cs.info["x_extent_distinct"], 3),
                        "n_periodic_copies": cs.info["n_periodic_copies"],
                        "n_degenerate_sites": n_deg, "E_min_core": float(e_core.min()) if len(e_core) else float("nan"),
                        "E_mean_core": float(e_core.mean()) if len(e_core) else float("nan"),
                        "r_max_tabulated": float(cs.r.max()), "missing_sites": mc["missing_sites"], "n_missing": n_missing,
                        "flags": "; ".join(flags)}
                out[key] = (e_core, n_missing, info)
                summary.append(info)
                if verbose:
                    print(f"  sites {pot} {ch}: {len(table)} in table, {cs.n_core_tabulated} core (r < "
                          f"{mc['core_site_radius']} A of the 0 K DXA line) of {cs.n_core_expected:.1f} per period "
                          f"({cs.coverage:.0%}); E_min {info['E_min_core']:.3f} eV; {n_deg} degenerate"
                          + (f"  [{info['flags']}]" if flags else ""))
    return out, summary


def _cmd_across(args) -> int:
    from .across import STAGE4_COLUMNS, load_runs, stage4
    from .provenance import environment
    from .results import tool_version, write_table

    files = sorted({str(p) for f in args.stage3 for p in (Path(f).glob("*_stage3_*.csv") if Path(f).is_dir() else [f])})
    if not files:
        print("no Stage-3 tables (*_stage3_*.csv) found")
        return 1
    runs, label = load_runs(files)
    d = load_definitions(args.definitions or label or DEFAULT_LABEL)
    if label and label != d.label:
        print(f"Stage-3 tables are {label}, --definitions is {d.label}")
        return 1
    env = environment()
    tv = tool_version(env)
    site_cores, site_summary = ({}, [])
    if args.sites and "mclean" not in d.params:
        print(f"definition set {d.label} has no McLean core-site rule (D2 section 12) -- --sites ignored")
    elif args.sites:
        site_cores, site_summary = load_site_cores(args.sites, runs, d, cache_root=args.cache)
    if "hysteresis" not in d.params:
        print(f"definition set {d.label}: uptake exponents and hysteresis differences only (verdict, McLean, x_bc error: D2)")
    comp = d["composition"]
    rows = stage4(runs, d, site_cores, gamma_radii=(comp["headline_radius"], comp["check_radius"]))
    for r in rows:   # since D2 no Stage-4 row depends on an open proposal (P5 is not used here): proposal_set stays empty
        r.update(definition_set=d.label, proposal_set="", tool_version=tv)
    prov = {"environment": env, "tool_version": tv, "definition_set": d.label, "definitions_sha256": d.sha256,
            "proposal_set": None,
            "inputs": [{"file": r.file, "run_id": r.run_id, "tool_version": r.prov.get("tool_version")} for r in runs],
            "site_tables": site_summary, "sites_config": str(args.sites) if args.sites else None}
    write_table(args.out, rows, prov, columns=STAGE4_COLUMNS)
    n_sets = len({(r["project"], r["potential"], r["character"]) for r in rows})
    by = {}
    for r in rows:
        by[r["analysis"]] = by.get(r["analysis"], 0) + 1
    print(f"{len(runs)} runs, {n_sets} run sets, {len(rows)} rows {by} -> {args.out}")
    return 0


def _cmd_plot(args) -> int:
    """disegreg plot gamma | timeseries | delta | uptake | regions (house style, PNG/PDF)."""
    import json

    import numpy as np

    from . import plotting
    from .results import read_table

    fmts = tuple(args.format)
    kind = args.kind
    if kind == "gamma":
        f = Path(args.input)
        rows = read_table(f)
        refs = args.ref or sorted({r["ref"] for r in rows})
        rows = [r for r in rows if r["ref"] in refs]
        for r in rows:
            r["R"] = float(r["R"])
        run = rows[0]["run_id"] if rows else f.stem
        fig = plotting.gamma_curve(rows, title=f"Γ(R) — {run}")
        out = args.out or f.with_suffix("")
    elif kind == "timeseries":
        f = Path(args.input)  # <run>_stage2_<D>.csv; the Stage-3 table next to it gives the window
        q, region, basis, ref = (args.quantity.split(":") + ["", "", ""])[:4]
        rows = [r for r in read_table(f) if (r["quantity"], r["region"], r["basis"], r["ref"]) == (q, region, basis, ref)]
        if not rows:
            print(f"no rows for {args.quantity} in {f}")
            return 1
        rows.sort(key=lambda r: r["timestep"])
        f3 = f.with_name(f.name.replace("_stage2_", "_stage3_"))
        prov = json.loads(f.with_name(f.name + ".provenance.json").read_text())
        dt = prov.get("manifest_data", {}).get("md_timestep_ps")
        scale, unit = (dt * 1e-3, "time (ns)") if dt else (1.0, "MD step")
        t = np.array([r["timestep"] for r in rows]) * scale
        y = np.array([r["value"] for r in rows])
        window = mean = sem = flat = slope = None
        if f3.is_file():
            info3 = json.loads(f3.with_name(f3.name + ".provenance.json").read_text()).get("stage3", {})
            k = 1e-3 if info3.get("time_unit") == "ps" else 1.0
            window = (info3["window_t_start"] * k, info3["window_t_end"] * k) if dt else (info3["window_t_start"],
                                                                                          info3["window_t_end"])
            r3 = next((r for r in read_table(f3) if (r["quantity"], r["region"], r["basis"], r["ref"])
                       == (q, region, basis, ref)), None)
            if r3:
                mean, sem, flat = r3["mean"], r3["sem"], r3["flat"]
                span = (r3["t_last"] - r3["t_first"]) * (k if dt else 1.0)
                slope = r3["change"] / span if span > 0 else None
        fig = plotting.timeseries(t, y, window, mean, sem, flat, ylabel=plotting.quantity_label(args.quantity),
                                  title=f"{prov.get('manifest_data', {}).get('run_id', f.stem)}\n{args.quantity.rstrip(':')}",
                                  time_label=unit, slope=slope)
        out = args.out or f.with_name(f.stem + "_" + args.quantity.replace(":", "-").rstrip("-"))
    elif kind == "delta":
        f = Path(args.input)  # results table of `disegreg table`
        q, region, basis, ref = (args.quantity.split(":") + ["", "", ""])[:4]
        rows = [r for r in read_table(f) if (r["quantity"], r["region"], r["basis"], r["ref"]) == (q, region, basis, ref)]
        for r in rows:
            r["delta"] = float(r["delta"])
        fig = plotting.vs_delta(rows, ylabel=plotting.quantity_label(args.quantity), log=args.log)
        out = args.out or f.with_name(f.stem + "_vs-delta_" + args.quantity.replace(":", "-").rstrip("-"))
    elif kind == "uptake":
        f = Path(args.input)  # Stage-4 table of `disegreg across`
        rows = read_table(f)
        for r in rows:
            for k in ("value", "error", "delta", "delta_2"):
                r[k] = float(r[k]) if r[k] != "" else float("nan")
        if args.ref:
            rows = [r for r in rows if r["ref"] in args.ref]
        fig = plotting.uptake(rows)
        out = args.out or f.with_name(f.stem + "_uptake")
    elif kind == "regions":
        f = Path(args.input)  # results table
        rows = [r for r in read_table(f) if r["run_id"] == args.run and r["ref"] == ""]
        xr = {r["ref"]: r["value"] for r in read_table(f) if r["run_id"] == args.run and r["quantity"] == "x_ref"}
        fig = plotting.regions(rows, args.run, refs={f"x_ref ({k})": v for k, v in xr.items() if k != "c_B"})
        out = args.out or f.with_name(f"{args.run}_regions")
    elif kind == "alongline":
        f = Path(args.input)  # <run>_stage2_<D>.csv with the P4 quantities
        rows = read_table(f)
        prov = json.loads(f.with_name(f.name + ".provenance.json").read_text())
        run = prov.get("manifest_data", {}).get("run_id", f.stem)
        dt = prov.get("manifest_data", {}).get("md_timestep_ps")
        steps = sorted({r["timestep"] for r in rows if r["quantity"] == "alongline_N_bin"})
        if not steps:
            print(f"no along-line quantities in {f} (analysed with --no-proposals?)")
            return 1
        f3 = f.with_name(f.name.replace("_stage2_", "_stage3_"))
        if f3.is_file():  # only the snapshots of the averaging window
            info3 = json.loads(f3.with_name(f3.name + ".provenance.json").read_text()).get("stage3", {})
            k = (dt or 1.0)
            steps = [t for t in steps if info3["window_t_start"] - 1e-6 <= t * k <= info3["window_t_end"] + 1e-6] or steps
        nb = int(next(r["value"] for r in rows if r["quantity"] == "alongline_n_bins"))
        L = next(r["value"] for r in rows if r["quantity"] == "L_x")
        N = {(r["timestep"], int(r["region"][1:])): r["value"] for r in rows if r["quantity"] == "alongline_N_bin"}
        Ns = {(r["timestep"], int(r["region"][1:])): r["value"] for r in rows if r["quantity"] == "alongline_N_sol_bin"}
        d = np.zeros((nb, len(steps)))
        for j, t in enumerate(steps):
            n = np.array([N[(t, k)] for k in range(nb)])
            sn = np.array([Ns[(t, k)] for k in range(nb)])
            d[:, j] = sn - (sn.sum() / n.sum() if n.sum() else 0.0) * n
        S = {}
        for r in rows:
            if r["quantity"] == "alongline_S" and r["timestep"] in steps:
                S.setdefault(int(r["region"][1:]), []).append(r["value"])
        q = np.array(sorted(S))
        Sm = np.array([np.mean(S[k]) for k in q])
        Se = np.array([np.std(S[k], ddof=1) / np.sqrt(len(S[k])) if len(S[k]) > 1 else np.nan for k in q])
        tt = np.array(steps, float) * ((dt or 1.0) * (1e-3 if dt else 1.0))
        te = np.concatenate([[tt[0] - 0.5 * (tt[1] - tt[0]) if len(tt) > 1 else tt[0] - 0.5],
                             0.5 * (tt[1:] + tt[:-1]), [tt[-1] + 0.5 * (tt[-1] - tt[-2]) if len(tt) > 1 else tt[0] + 0.5]])
        base = args.out or f.with_name(f.stem + "_alongline")
        fig = plotting.signed_map(np.linspace(0, L, nb + 1), te, d, "position along the line X (Å)",
                                  "time (ns)" if dt else "MD step", "core solute excess per bin (atoms)",
                                  title=f"{run}: core composition along the line", aspect="auto")
        for p in plotting.save(fig, str(base) + "_map", fmts):
            print(f"-> {p}")
        fig = plotting.alongline_spectrum(q, Sm, Se, L, len(steps), title=f"{run}: along-line spectrum of the core")
        out = str(base) + "_spectrum"
    else:  # pragma: no cover -- argparse restricts the choices
        return 1
    for p in plotting.save(fig, out, fmts):
        print(f"-> {p}")
    return 0


def _cmd_slurm(args) -> int:
    import glob

    from .manifest import load_manifest
    from .slurm import CORE_H_LIMIT, stage1_array

    dumps = []
    for g in args.dumps or []:
        dumps += sorted(glob.glob(g))
    for mf in args.manifest or []:
        dumps += [str(p) for p in load_manifest(mf).dump_files()]
    dumps = list(dict.fromkeys(dumps))  # keep order, drop duplicates
    if not dumps:
        print("no dumps found")
        return 1
    script, lst, n_tasks, core_h = stage1_array(dumps, args.out, args.env, args.cache, chunk=args.chunk,
                                                partition=args.partition, cpus=args.cpus, mem=args.mem, time=args.time,
                                                max_parallel=args.max_parallel, definitions=args.definitions,
                                                version=__version__)
    print(f"{len(dumps)} dumps -> {n_tasks} array tasks; estimated {core_h:.1f} core-h")
    print(f"  script {script}\n  list   {lst}")
    if core_h > CORE_H_LIMIT:
        print(f"  WARNING: more than {CORE_H_LIMIT:.0f} core-h -- needs Erik's explicit OK before submitting")
    print(f"submit with:  sbatch {script}")
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="disegreg", description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("version", help="print tool version and default definition set").set_defaults(func=_cmd_version)
    p1 = sub.add_parser("stage1", help="read dumps with OVITO (a-CNA, PTM, DXA) and write the Stage-1 cache")
    p1.add_argument("dumps", nargs="+", help="LAMMPS dump files (quenched snapshots), plain or .gz")
    p1.add_argument("--cache", help="cache folder (default: $DISEGREG_CACHE or the agreed location for this machine)")
    p1.add_argument("--definitions", default=DEFAULT_LABEL, help=f"definition set label or YAML path (default {DEFAULT_LABEL})")
    p1.add_argument("--recompute", action="store_true", help="recompute even if a valid cache entry exists")
    p1.set_defaults(func=_cmd_stage1)
    pm = sub.add_parser("manifest", help="Stage 0: print a manifest template or check manifests")
    pm.add_argument("action", choices=("template", "check"))
    pm.add_argument("files", nargs="*", help="manifest YAML files (for check)")
    pm.set_defaults(func=_cmd_manifest)
    pa = sub.add_parser("analyse", help="Stage 2 for every dump of a run (Stage 1 from the cache, computed if missing)")
    pa.add_argument("manifests", nargs="+", help="run manifest YAML files")
    pa.add_argument("--out", default=".", help="folder for <run_id>_stage2_<D>.csv (default: current folder)")
    pa.add_argument("--x-ref", type=float, help="extra reference solute fraction, e.g. the bulk-control value (ref x_bc "
                    "if the manifest says bulk_control); c_B and the manifest's x_iso are always used")
    pa.add_argument("--last", type=int, help="only the last N dumps of each run")
    pa.add_argument("--cache", help="cache folder (default: $DISEGREG_CACHE or the agreed location for this machine)")
    pa.add_argument("--definitions", default=DEFAULT_LABEL, help=f"definition set label or YAML path (default {DEFAULT_LABEL})")
    pa.add_argument("--drift-criterion", choices=("change_se", "sem"), default=None,
                    help="change_se: |change| < 2 SE(change) (default, D1 rev. 4); sem: |change| < 2 SEM (rev. 3 / old "
                    "scripts); the other verdict is always stored as well")
    pa.add_argument("--bulk-dir", nargs="*", help="extra folders with bulk-control Stage-3 tables (default: --out only); "
                    "a run whose manifest names a bulk_control takes x_bc from <bulk_control>_stage3_<D>.csv")
    pa.add_argument("--no-proposals", action="store_true", help="quantities of the definition set only (no open-proposal quantities, i.e. no P5 ordered domains)")
    pa.add_argument("--verbose", action="store_true")
    pa.set_defaults(func=_cmd_analyse)
    ps = sub.add_parser("slurm", help="write a Slurm job array for Stage 1 (cmmg); prints the sbatch command")
    ps.add_argument("--dumps", nargs="*", help="glob(s) of dump files (quote them)")
    ps.add_argument("--manifest", nargs="*", help="manifest(s) whose dump glob to use")
    ps.add_argument("--out", required=True, help="sbatch script to write (the dump list goes next to it)")
    ps.add_argument("--env", default="/cmmc/ptmp/$USER/conda-envs/disegreg", help="conda env prefix with disegreg")
    ps.add_argument("--cache", default="/u/$USER/PTMP/DISEGREG-CACHE", help="Stage-1 cache folder")
    ps.add_argument("--chunk", type=int, default=20, help="dumps per array task (default 20)")
    ps.add_argument("--partition", default="s.cmmg")
    ps.add_argument("--cpus", type=int, default=4)
    ps.add_argument("--mem", default="16G")
    ps.add_argument("--time", default="02:00:00")
    ps.add_argument("--max-parallel", type=int, default=20, help="array tasks running at once")
    ps.add_argument("--definitions", default=DEFAULT_LABEL)
    ps.set_defaults(func=_cmd_slurm)
    pt = sub.add_parser("table", help="merge the Stage-3 tables of many runs into one tidy results table")
    pt.add_argument("stage3", nargs="+", help="*_stage3_<D>.csv files or folders containing them")
    pt.add_argument("--out", required=True, help="output CSV (a .provenance.json is written next to it)")
    pt.set_defaults(func=_cmd_table)
    px = sub.add_parser("across", help="Stage 4 across runs: uptake exponent, McLean, hysteresis (from Stage-3 tables)")
    px.add_argument("stage3", nargs="+", help="*_stage3_<D>.csv files or folders containing them")
    px.add_argument("--out", required=True, help="output CSV (a .provenance.json is written next to it)")
    px.add_argument("--sites", help="site-energy configuration YAML for the McLean prediction (examples/stage4/sites.yaml)")
    px.add_argument("--cache", help="Stage-1 cache folder for the 0 K cells (default: $DISEGREG_CACHE or the agreed location)")
    px.add_argument("--definitions", default=None, help="definition set label or YAML path (default: the tables' own)")
    px.set_defaults(func=_cmd_across)
    pp = sub.add_parser("plot", help="plots in the house style: gamma | timeseries | delta | uptake | regions | alongline")
    pp.add_argument("kind", choices=("gamma", "timeseries", "delta", "uptake", "regions", "alongline"))
    pp.add_argument("input", help="gamma: <run>_gammaR_<D>.csv; timeseries: <run>_stage2_<D>.csv; delta, regions: results "
                    "table of `disegreg table`; uptake: Stage-4 table of `disegreg across`")
    pp.add_argument("--quantity", default="N_sol_inside:20:r",
                    help="quantity:region:basis:ref (timeseries, delta), e.g. c:C or Gamma:20:r:x_iso")
    pp.add_argument("--ref", nargs="*", help="only these references (gamma, uptake)")
    pp.add_argument("--run", help="run id (regions)")
    pp.add_argument("--log", action="store_true", help="logarithmic y axis (delta)")
    pp.add_argument("--out", help="output file name without extension (default: next to the input)")
    pp.add_argument("--format", nargs="+", default=["png"], help="png and/or pdf (default png)")
    pp.set_defaults(func=_cmd_plot)
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
