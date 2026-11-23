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


def _cmd_analyse(args) -> int:
    import time
    from pathlib import Path

    from .cache import default_cache_root, get_snapshot
    from .manifest import load_manifest
    from .provenance import environment
    from .results import tool_version, write_table
    from .stage2 import analyse_snapshot

    d = load_definitions(args.definitions)
    root = args.cache or default_cache_root()
    env = environment()
    tv = tool_version(env)
    status = 0
    for mf in args.manifests:
        m = load_manifest(mf)
        x_ref, x_src = args.x_ref, "command line"
        if x_ref is None and m["x_ref_source"] == "isotherm":
            x_ref, x_src = float(m["x_iso"]), "x_iso (manifest)"
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
            res = analyse_snapshot(snap, d, float(m["a_cell"]), z_gp=float(m["z_gp"]), x_ref=x_ref)
            ts = int(snap.meta.get("timestep", -1))
            for r in res.records:
                rows.append({"run_id": m["run_id"], "snapshot": Path(dump).name, "timestep": ts, **r,
                             "definition_set": d.label, "tool_version": tv})
            info = {k: v for k, v in res.info.items()}
            snaps.append({"dump": str(dump), "source_sha256": snap.meta.get("source_sha256"), "timestep": ts,
                          "ovito": snap.meta.get("provenance", {}).get("ovito"), "info": info})
            flag = "  z_gp FLAG" if res.info["z_gp_flag"] else ""
            print(f"{m['run_id']} step {ts}: mode {res.info['mode']}, f_cov {res.info['f_cov']:.3f}, "
                  f"c_C {res.value('c', region='C'):.4f}, c_B {res.value('c', region='B'):.5f}"
                  f"{flag}  ({time.perf_counter() - t0:.1f} s)")
        out = Path(args.out).expanduser() / f"{m['run_id']}_stage2_{d.label}.csv"
        prov = {"environment": env, "tool_version": tv, "definition_set": d.label, "definitions_sha256": d.sha256,
                "definitions_status": d.status, "manifest": str(m.path), "manifest_data": m.data,
                "x_ref": x_ref, "x_ref_source": x_src, "snapshots": snaps}
        write_table(out, rows, prov)
        print(f"  -> {out} ({len(rows)} rows)")
    return status


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
    pa.add_argument("--x-ref", type=float, help="reference solute fraction (default: x_iso if the manifest says isotherm)")
    pa.add_argument("--last", type=int, help="only the last N dumps of each run")
    pa.add_argument("--cache", help="cache folder (default: $DISEGREG_CACHE or the agreed location for this machine)")
    pa.add_argument("--definitions", default="D1", help="definition set label or YAML path (default D1)")
    pa.add_argument("--verbose", action="store_true")
    pa.set_defaults(func=_cmd_analyse)
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
