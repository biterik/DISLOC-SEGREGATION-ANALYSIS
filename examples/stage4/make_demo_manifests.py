#!/usr/bin/env python3
"""DEMO manifests for Stage 4 (v0.2): every Nb-Au thread-07 run (Mendelev 2022 Nb, EAM-X Au; 4 characters x 3 offsets)
and every Ni-Cu thread-07 run of Lee-Shim MEAM and ACE that has snapshots on the Mac.

Written by the analysis-tool instance from the project files (read-only), following the conventions of the projects'
own first manifests of 24 Nov (Nb-Au: ...0.25kT.../*.disegreg.yaml; Ni-Cu: ...d60_LeeShim2004-MEAM_0.10kT.../*.disegreg.yaml):
z_gp = builder glide gap 124.9847 A (at a = 3.52 A) x a_cell / 3.52; a_cell, x_iso from the thread's mdmc-parameters.dat;
Delta = (dmu - dmu_coex) / kT. The PROJECT instances own the real manifests; these are for the v0.2 demo only.

    python examples/stage4/make_demo_manifests.py      (M5; writes examples/stage4/manifests/{nbau,nicu}/*.yaml)
"""

from __future__ import annotations

import re
from pathlib import Path

import yaml

SIM = Path("~/Desktop/SIMULATIONS").expanduser()
OUT = Path(__file__).resolve().parent / "manifests"
KB = 8.617333262e-5
GAP_352 = 124.9847  # 00_CORRECTED-CELLS-A0-3.52/cells/gate-report_Ni_*: glide gap z at a = 3.52 A (all four characters)
FRAMES = {"d0": ([0, -1, 1], [2, -1, -1], [1, 1, 1]), "d60": ([0, -1, 1], [2, -1, -1], [1, 1, 1]),
          "d30": ([1, -2, 1], [1, 0, -1], [1, 1, 1]), "d90": ([1, -2, 1], [1, 0, -1], [1, 1, 1])}
DEMO = "DEMO manifest for disegreg v0.2 Stage 4 (examples/stage4/make_demo_manifests.py); the project owns the real one"


def params(path, solute_col: bool):
    """mdmc-parameters.dat -> {(potential, offset): (dmu, a, x, status text)}"""
    out = {}
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        p = line.split()
        if solute_col:  # Nb-Au: potential solute offset dmu a x status ...
            pot, off, dmu, a, x, rest = p[0], float(p[2]), float(p[3]), float(p[4]), float(p[5]), " ".join(p[6:])
        else:           # Ni-Cu: potential offset dmu a x status ...
            pot, off, dmu, a, x, rest = p[0], float(p[1]), float(p[2]), float(p[3]), float(p[4]), " ".join(p[5:])
        out[(pot, f"{off:.2f}")] = (dmu, a, x, rest)
    return out


def coex(rest: str) -> float:
    return float(re.search(r"dmu_coex\(300 K\) = (-?[\d.]+)", rest).group(1))


def tilde(p) -> str:
    """Paths as ~/... (the repository is public: no local user names in committed files)."""
    return str(p).replace(str(Path.home()), "~", 1)


def manifest(project, thread, run_dir: Path, pot, host, solute, ch, a, dmu, dmu_coex, x_iso, x_iso_source, start,
             seeded_from=None, stream_col=None, coex_source=""):
    run_id = run_dir.name
    d = {"schema_version": 1, "project": project, "thread": thread, "run_id": run_id, "potential": pot,
         "host_element": host, "solute_element": solute, "character": ch,
         "frame": dict(zip("XYZ", FRAMES[ch])), "frame_source": "Ni cells of 00_CORRECTED-CELLS-A0-3.52 (LAMMPS frame X = line)",
         "a_cell": a, "a_cell_source": "mdmc-parameters.dat a(x,T)", "temperature": 300, "dmu": dmu, "dmu_lammps": dmu,
         "dmu_lammps_sign": 1, "dmu_coex": dmu_coex, "coex_line": "incoherent",
         "dmu_coex_source": f"mdmc-parameters.dat {coex_source}".strip(),
         "delta": round((dmu - dmu_coex) / (KB * 300), 4), "dumps": tilde(run_dir / "*_quenched_step*.dump"),
         "quenched": True, "md_timestep_ps": 0.001, "z_boundary": "fixed", "z_gp": round(GAP_352 * a / 3.52, 4),
         "z_gp_source": f"builder glide gap {GAP_352} A at a = 3.52 A x a_cell / 3.52", "bulk_control": None,
         "bulk_control_reason": "no bulk control on the Mac yet (demo)", "x_ref_source": "isotherm", "x_iso": x_iso,
         "x_iso_source": x_iso_source, "start_type": start, "notes": DEMO}
    if seeded_from:
        d["seeded_from"] = seeded_from
    streams = sorted(run_dir.glob("*_mc-stream.dat"))
    if streams:
        d["mc_stream"] = tilde(streams[0])
        if stream_col:
            d["mc_stream_column"] = stream_col
    return d


def nbau():
    project, thread = "Nb-Au-AT-DISLOCS-Ni-ACE-EAM", "07_MDMC-SOLUTE-SEGREGATION-T300-WHOLE-INTERIOR"
    base = SIM / project / thread
    par = params(base / "mdmc-parameters.dat", solute_col=True)
    out = []
    for pot, solute in (("Mendelev2022-EAMFS", "Nb"), ("EAMX-DawChandross2023", "Au")):
        for ch in ("d0", "d30", "d60", "d90"):
            for off in ("0.25", "0.50", "1.00"):
                dmu, a, x, rest = par[(pot, off)]
                run = base / f"mdmc300K_{ch}_{pot}_{off}kT-below-coex_dmu{dmu:.4f}"
                if not run.is_dir() or not list(run.glob("*_quenched_step*.dump")):
                    continue
                out.append(manifest(project, thread, run, pot, "Ni", solute, ch, a, dmu, coex(rest), x,
                                    "mdmc-parameters.dat xSol_overall", "bottom-up", stream_col="nSol_interior",
                                    coex_source=rest.split(";")[0]))
    return out


def nicu():
    project, thread = "Cu-AT-DISLOCS-Ni-EAM-MEAM-ACE", "07_MDMC-CU-SEGREGATION-T300-COEX-OFFSETS"
    base = SIM / project / thread
    par = params(base / "mdmc-parameters.dat", solute_col=False)
    # 0.10 kT rows of mdmc-parameters.dat are "nominal, not used" (seeded runs): Lee-Shim x_iso as in the project's own
    # d60 0.10 kT manifest (isotherm up-branch, log-linear at dmu 0.8856); ACE: the nominal value, flagged
    x010 = {"LeeShim2004-MEAM": (7.8642e-3, "isotherm_LeeShim2004-MEAM.dat up-branch log-linear at dmu 0.8856 "
                                            "(as the project's d60 0.10 kT manifest)"),
            "CuNi-ACE": (3.5585e-4, "mdmc-parameters.dat 0.10 kT row, marked 'nominal, not used' -- DEMO ONLY")}
    seed = {"LeeShim2004-MEAM": "0.25", "CuNi-ACE": "0.50"}
    out = []
    for pot in ("LeeShim2004-MEAM", "CuNi-ACE"):
        for ch in ("d0", "d30", "d60", "d90"):
            for off in ("0.10", "0.25", "0.50", "1.00", "2.00"):
                if (pot, off) not in par:
                    continue
                dmu, a, x, rest = par[(pot, off)]
                run = base / f"mdmc300K_{ch}_{pot}_{off}kT-below-coex_dmu{dmu:.4f}"
                if not run.is_dir() or not list(run.glob("*_quenched_step*.dump")):
                    continue
                start, seeded, xs = "bottom-up", None, "mdmc-parameters.dat xCu_overall (thread-02 up-scan)"
                if off == "0.10":
                    x, xs = x010[pot]
                    dmu_s = par[(pot, seed[pot])][0]
                    start, seeded = "seeded", f"mdmc300K_{ch}_{pot}_{seed[pot]}kT-below-coex_dmu{dmu_s:.4f}"
                    a = par[(pot, seed[pot])][1]   # the seeded run continues the cell of its parent
                out.append(manifest(project, thread, run, pot, "Ni", "Cu", ch, a, dmu, coex(rest), x, xs, start, seeded,
                                    stream_col="nCu_interior", coex_source=rest.split(";")[0]))
    return out


if __name__ == "__main__":
    for name, ms in (("nbau", nbau()), ("nicu", nicu())):
        d = OUT / name
        d.mkdir(parents=True, exist_ok=True)
        for m in ms:
            (d / f"{m['run_id']}.yaml").write_text(f"# {DEMO}\n" + yaml.safe_dump(m, sort_keys=False, width=140))
        print(f"{name}: {len(ms)} manifests -> {d}")
