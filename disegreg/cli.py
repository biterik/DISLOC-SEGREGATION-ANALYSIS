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
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
