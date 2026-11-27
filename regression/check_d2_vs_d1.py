#!/usr/bin/env python3
"""D2 must give every D1 number unchanged: compare two analyses of the same runs, cell by cell.

    python regression/check_d2_vs_d1.py regression/output/v0.2-final regression/output/v0.2-D2

The first folder holds <run>_stage2|stage3|gammaR_D1.csv (D1, here with the v0.2 proposal quantities of 24 Nov), the second
the same runs analysed with D2. For every run and table: same rows (by key), same values in every column except
definition_set and tool_version -- for every D1 quantity (rows of the v0.2 additions, alongline_* and domain_*, are compared
separately and listed: since the review fix of 25 Nov, Stage 3 of D2 carries the section-16 window test instead of window
means of the per-snapshot peak values). Exit code 0 = every D1 quantity identical.
"""

from __future__ import annotations

import csv
import sys
from collections import Counter
from pathlib import Path

IGNORE = {"definition_set", "tool_version"}
V02 = ("alongline_", "domain_")   # quantities added in v0.2 (not D1)
KEYS = {"stage2": ("snapshot", "quantity", "region", "basis", "ref"), "stage3": ("quantity", "region", "basis", "ref"),
        "gammaR": ("R", "basis", "ref")}


def rows(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def main(a: str, b: str) -> int:
    a, b = Path(a), Path(b)
    bad, n_cells, n_tables = 0, 0, 0
    extra, v02a, v02b = Counter(), Counter(), Counter()
    for fa in sorted(a.glob("*_D1.csv")):
        kind = fa.stem.rsplit("_", 2)[-2]
        if kind not in KEYS:
            continue
        fb = b / fa.name.replace("_D1.csv", "_D2.csv")
        if not fb.is_file():
            print(f"MISSING in {b}: {fb.name}")
            bad += 1
            continue
        ra, rb = rows(fa), rows(fb)
        is_v02 = (lambda r: r.get("quantity", "").startswith(V02))
        for r in ra:
            if is_v02(r):
                v02a[(kind, r["quantity"])] += 1
        for r in rb:
            if is_v02(r):
                v02b[(kind, r["quantity"])] += 1
        ra, rb = [r for r in ra if not is_v02(r)], [r for r in rb if not is_v02(r)]
        ka = {tuple(r[k] for k in KEYS[kind]): r for r in ra}
        kb = {tuple(r[k] for k in KEYS[kind]): r for r in rb}
        missing = set(ka) - set(kb)
        for k in set(kb) - set(ka):
            extra[(kind, kb[k]["quantity"] if "quantity" in kb[k] else "gammaR")] += 1
        diff = [(k, c) for k in set(ka) & set(kb) for c in ka[k] if c not in IGNORE and ka[k][c] != kb[k].get(c)]
        n_cells += sum(len(r) - len(IGNORE) for r in ra)
        n_tables += 1
        if missing or diff:
            bad += 1
            print(f"DIFFERENT {fa.name}: {len(missing)} rows missing in D2, {len(diff)} cells differ, e.g. {diff[:3]}")
    labels = Counter((r["definition_set"]) for fb in b.glob("*_stage3_D2.csv") for r in rows(fb))
    print(f"{n_tables} tables, {n_cells} cells compared; {bad} table(s) different; D2 label in Stage 3: {dict(labels)}")
    if extra:
        print(f"rows only in D2 (D1 quantities): {dict(extra)}")
    changed = {k: (v02a.get(k, 0), v02b.get(k, 0)) for k in set(v02a) | set(v02b) if v02a.get(k, 0) != v02b.get(k, 0)}
    print(f"v0.2 quantities (alongline_*, domain_*) with a different row count (first, second): {changed or 'none'}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:3]))
