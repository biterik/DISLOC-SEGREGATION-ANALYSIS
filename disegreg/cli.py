"""Command-line entry point ``disegreg``. Sub-commands are added stage by stage (see docs/DESIGN.md)."""

import argparse
import sys

from . import __version__
from .definitions import load_definitions


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


def _cmd_analyse(args) -> int:
    import time
    from pathlib import Path

    from .cache import default_cache_root, get_snapshot
    from .manifest import load_manifest
    from .provenance import environment
    from .results import GAMMA_CURVE_COLUMNS, STAGE3_COLUMNS, gamma_curve, paper_hold, tool_version, write_table
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
    for mf in args.manifests:
        m = load_manifest(mf)
        # reference compositions (DEFINITIONS 5, rev. 4): c_B of the snapshot always; x_iso if the manifest has it;
        # a value from the command line as x_bc (bulk control) or "given"
        refs = {}
        if args.x_ref is not None:
            refs["x_bc" if m["x_ref_source"] == "bulk_control" else "given"] = float(args.x_ref)
        if m.get("x_iso") is not None:
            refs["x_iso"] = float(m["x_iso"])
        has_bc = bool(m.get("bulk_control")) and "x_bc" in refs
        if not has_bc:
            missing = [r for r in ref_rule["without_bulk_control_report_against"] if r != "c_B" and r not in refs]
            if missing:
                print(f"{m['run_id']}: no bulk control and no {missing} in the manifest -- reported against c_B only")
        rule_text = ("x_bc (bulk control), c_B and x_iso as diagnostics" if has_bc else
                     "no bulk control: reported against both c_B and x_iso (DEFINITIONS 5, rev. 4)")
        # the hold is lifted only when x_bc is actually given, not when the manifest merely names a bulk control
        hold = "" if has_bc else paper_hold({**m.data, "bulk_control": None},
                                            ref_rule.get("paper_hold_without_bulk_control"))
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
                                   frame_matrix=m.frame_matrix(), c_B_ref=True)
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
        out = Path(args.out).expanduser() / f"{m['run_id']}_stage2_{d.label}.csv"
        prov = {"environment": env, "tool_version": tv, "definition_set": d.label, "definitions_sha256": d.sha256,
                "definitions_status": d.status, "manifest": str(m.path), "manifest_data": m.data,
                "x_ref": refs, "x_ref_source": m["x_ref_source"], "reference_rule": rule_text, "paper_flag": hold,
                "snapshots": snaps}
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
                     paper_flag=hold)
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
    p1.add_argument("--definitions", default="D1", help="definition set label or YAML path (default D1)")
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
    pa.add_argument("--definitions", default="D1", help="definition set label or YAML path (default D1)")
    pa.add_argument("--drift-criterion", choices=("change_se", "sem"), default=None,
                    help="change_se: |change| < 2 SE(change) (default, D1 rev. 4); sem: |change| < 2 SEM (rev. 3 / old "
                    "scripts); the other verdict is always stored as well")
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
    ps.add_argument("--definitions", default="D1")
    ps.set_defaults(func=_cmd_slurm)
    pt = sub.add_parser("table", help="merge the Stage-3 tables of many runs into one tidy results table")
    pt.add_argument("stage3", nargs="+", help="*_stage3_<D>.csv files or folders containing them")
    pt.add_argument("--out", required=True, help="output CSV (a .provenance.json is written next to it)")
    pt.set_defaults(func=_cmd_table)
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
