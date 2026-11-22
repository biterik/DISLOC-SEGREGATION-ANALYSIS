"""Command-line entry point ``disegreg``. Sub-commands are added stage by stage (see docs/DESIGN.md)."""

import argparse
import sys

from . import __version__
from .definitions import load_definitions


def _cmd_version(_args) -> int:
    d = load_definitions()
    print(f"disegreg {__version__}  definition set {d.label} ({d.status})")
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="disegreg", description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("version", help="print tool version and default definition set").set_defaults(func=_cmd_version)
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
