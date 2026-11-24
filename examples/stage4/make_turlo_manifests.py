#!/usr/bin/env python3
"""DEMO manifests for the hysteresis test (v0.2, proposal P2): the TURLO-REDO thread-51 Cu-Zr d90 runs at 800 K, bottom-up
(U: from pure Cu) and top-down (D: from the t50 R3 complexion) at Delta = 0, +1, +2 kT, and D at -0.5 kT.

Every number comes from TURLO-REDO's OWN generator, 52_PAPER-CONTRIBUTIONS-CU-ZR/DISEGREG/make-manifests.py (2026-11-24),
imported read-only: its module code only reads project files; its one write (into the run folder) is redirected here, and the
relative dump / stream paths are turned into ~/... paths to the run folders. The project owns the real manifests.

    python examples/stage4/make_turlo_manifests.py      (M5; writes examples/stage4/manifests/turlo/*.yaml)
"""

from __future__ import annotations

import builtins
import importlib.util
import os
from pathlib import Path

SIM = Path("~/Desktop/SIMULATIONS").expanduser()
GEN = SIM / "TURLO-REDO" / "52_PAPER-CONTRIBUTIONS-CU-ZR" / "DISEGREG" / "make-manifests.py"
OUT = Path(__file__).resolve().parent / "manifests" / "turlo"
TAGS = {"U00": "bottom-up", "U10": "bottom-up", "U20": "bottom-up", "D00": "top-down", "D10": "top-down",
        "D20": "top-down", "Dm05": "top-down"}


def main():
    spec = importlib.util.spec_from_file_location("turlo_make_manifests", GEN)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)          # reads the builder / coexistence / isotherm files, writes nothing
    for tag, start in TAGS.items():
        note = ("t50 R3 final thermal state (step 2050000, restart START_t50-R3-final-2ns_step2050000_restart.lmprst)"
                if start == "top-down" else "the relaxed 0 K pure-Cu d90 cell (33_RELAX-0K-D90-CELL)")
        mod.RUNS[tag] = (mod.T51, "t51", start, note)
    OUT.mkdir(parents=True, exist_ok=True)
    real_open = builtins.open

    def redirect(path, mode="r", *a, **kw):  # the generator's only write: <run folder>/<run_id>.disegreg.yaml
        if "w" in mode and str(path).endswith(".disegreg.yaml"):
            run_dir = os.path.dirname(path)
            return _Rewrite(OUT / os.path.basename(path).replace(".disegreg.yaml", ".yaml"), run_dir)
        return real_open(path, mode, *a, **kw)

    mod.open = redirect
    for tag in TAGS:
        mod.manifest(tag)


class _Rewrite:
    """Collects the manifest text; on close makes dumps / mc_stream point to the run folder (as ~/...) and marks it DEMO."""

    def __init__(self, out: Path, run_dir: str):
        self.out, self.run_dir, self.buf = out, run_dir.replace(str(Path.home()), "~", 1), []

    def write(self, text):  # the generator writes the whole manifest in one call and never closes the file
        self.buf.append(text)
        self._flush()

    def close(self):
        pass

    def _flush(self):
        lines = []
        for line in "".join(self.buf).splitlines():
            for key in ("dumps: ", "mc_stream: "):
                if line.startswith(key):
                    line = f"{key}{self.run_dir}/{line[len(key):]}"
            lines.append(line)
        text = ("# DEMO copy for disegreg v0.2 (examples/stage4/make_turlo_manifests.py), generated with TURLO-REDO's own "
                "make-manifests.py; the project owns the real manifest\n" + "\n".join(lines) + "\n")
        self.out.write_text(text)


if __name__ == "__main__":
    main()
