"""Automatic x_bc and Stage 4 through the command line, with OVITO on small synthetic dumps (skipped without OVITO)."""

import json

import numpy as np
import pytest
import yaml

pytest.importorskip("ovito")
pytestmark = pytest.mark.ovito

from helpers import fcc_in_cell, screw_dislocation_cell, write_lammps_dump  # noqa: E402

from disegreg.cache import get_snapshot  # noqa: E402
from disegreg.cli import main  # noqa: E402
from disegreg.definitions import load_definitions  # noqa: E402
from disegreg.geometry import Cell  # noqa: E402
from disegreg.results import read_table  # noqa: E402

A = 3.52


def base_manifest(run_id, dumps):
    return {"schema_version": 1, "project": "PROJ", "run_id": run_id, "potential": "none", "host_element": "Ni",
            "solute_element": "Cu", "a_cell": A, "a_cell_source": "test", "temperature": 300, "dmu": 0.5,
            "dmu_lammps": 0.5, "dmu_lammps_sign": 1, "dmu_coex": 0.5, "coex_line": "incoherent", "dmu_coex_source": "t",
            "delta": 0.0, "dumps": dumps, "quenched": True, "z_boundary": "free", "x_iso": 0.05, "x_iso_source": "t"}


def test_bulk_control_feeds_x_bc_automatically(tmp_path):
    rng = np.random.default_rng(11)
    # bulk control: perfect fcc, free z faces, ~15 % solute (> 100 solute atoms in the interior), new types per snapshot;
    # smaller than the dislocation cell (2304 < 3360 atoms) as allowed by D2 section 15
    bcell = Cell(np.diag([6 * A, 6 * A, 16 * A]), pbc=(True, True, False))
    bpos = fcc_in_cell(bcell, A) + 1e-3
    bdir = tmp_path / "PROJ" / "T" / "BULK"
    bdir.mkdir(parents=True)
    z = bpos[:, 2]
    interior = (z - z.min() >= 15) & (z.max() - z >= 15)
    x_steps = []
    for step in (1000, 2000, 3000, 4000):
        btypes = np.where(rng.random(len(bpos)) < 0.15, 2, 1)
        x_steps.append((btypes[interior] == 2).mean())
        write_lammps_dump(bdir / f"bulk_quenched_step{step}.dump", bcell, np.arange(1, len(bpos) + 1), btypes, bpos,
                          timestep=step, pbc=("pp", "pp", "ss"))
    mb = base_manifest("BULK", "PROJ/T/BULK/*_quenched_step*.dump")
    mb["cell_type"] = "bulk_control"
    mb["n_atoms"] = len(bpos)
    mb["swap_attempts_per_atom_per_step"] = 0.0105        # within 10 % of the run's 0.01 (D2 section 15)
    # dislocation run naming the bulk control
    cell, pos, (y0, z0) = screw_dislocation_cell()
    rdir = tmp_path / "PROJ" / "T" / "RUN"
    rdir.mkdir(parents=True)
    types = np.where(rng.random(len(pos)) < 0.05, 2, 1)
    for step in (1000, 2000):
        write_lammps_dump(rdir / f"screw_quenched_step{step}.dump", cell, np.arange(1, len(pos) + 1), types, pos,
                          timestep=step, pbc=("pp", "ss", "ss"))
    md = {**base_manifest("RUN", "PROJ/T/RUN/*_quenched_step*.dump"), "character": "d0",
          "frame": {"X": [1, -1, 0], "Y": [1, 1, -2], "Z": [1, 1, 1]}, "frame_source": "t", "z_gp": float(z0),
          "z_gp_source": "t", "bulk_control": "BULK", "x_ref_source": "bulk_control", "start_type": "bottom-up",
          "swap_attempts_per_atom_per_step": 0.01}
    fb, fd = tmp_path / "bulk.yaml", tmp_path / "run.yaml"
    fb.write_text(yaml.safe_dump(mb))
    fd.write_text(yaml.safe_dump(md))
    assert main(["manifest", "check", str(fb), str(fd)]) == 0
    out = tmp_path / "out"
    # the dislocation manifest first on the command line: the bulk control is still analysed first
    assert main(["analyse", str(fd), str(fb), "--out", str(out), "--cache", str(tmp_path / "cache")]) == 0
    s3b = {(r["quantity"], r["region"]): r for r in read_table(out / "BULK_stage3_D2.csv")}
    n_win = s3b[("c", "interior")]["n_samples"]
    x_bc = float(np.mean(x_steps[-n_win:]))   # window mean (Stage 3 picks the final stretch)
    assert s3b[("c", "interior")]["mean"] == pytest.approx(x_bc) and s3b[("N_sol", "interior")]["mean"] > 100
    prov = json.loads((out / "RUN_stage3_D2.csv.provenance.json").read_text())
    assert prov["x_bc"]["value"] == pytest.approx(x_bc) and prov["x_ref"]["x_bc"] == pytest.approx(x_bc)
    assert prov["x_bc"]["enough_solute"] and "x_bc (bulk control)" in prov["reference_rule"]
    # D2 section 15: the bulk control is smaller than the dislocation cell -> x_bc accepted, but not for diagnostic (iv)
    assert len(bpos) < len(pos)
    assert prov["x_bc"]["n_atoms"] == len(bpos) and prov["x_bc"]["n_atoms_run"] == len(pos)
    assert prov["x_bc"]["smaller"] and prov["x_bc"]["accepted"] and prov["x_bc"]["diagnostic_iv"].startswith("not usable")
    assert prov["x_bc"]["same_run_length"] is False          # bulk control to step 4000, run to step 2000: flagged
    assert any("run length differs" in f for f in prov["x_bc"]["flags"])
    prov3s = json.loads((out / "RUN_stage3_D2.csv.provenance.json").read_text())
    assert prov3s["stage3"]["diagnostic_iv"].startswith("not usable") and prov["diagnostic_iv"] == prov3s["stage3"]["diagnostic_iv"]
    rows = read_table(out / "RUN_stage2_D2.csv")
    assert {r["ref"] for r in rows if r["quantity"] == "Gamma"} == {"c_B", "x_iso", "x_bc"}
    # Stage 4 over both: the bulk control is not part of a run set; x_bc error rows for the dislocation run
    assert main(["across", str(out), "--out", str(tmp_path / "s4.csv")]) == 0
    s4 = read_table(tmp_path / "s4.csv")
    assert {r["analysis"] for r in s4} == {"x_bc_error"} and {r["runs"] for r in s4} == {"RUN"}
    # the table still works, with a status column
    assert main(["table", str(out), "--out", str(tmp_path / "all.csv")]) == 0
    allrows = read_table(tmp_path / "all.csv")
    assert ({r["diagnostic_iv"] for r in allrows if r["run_id"] == "RUN"}
            == {"not usable: bulk control smaller than the dislocation cell"})
    assert {r["diagnostic_iv"] for r in allrows if r["run_id"] == "BULK"} == {""}   # the bulk control itself
    from disegreg.results import PROPOSAL_QUANTITIES
    assert {r["status"] for r in allrows if r["quantity"] in PROPOSAL_QUANTITIES} == {"proposal"}
    assert {r["status"] for r in allrows if r["quantity"] not in PROPOSAL_QUANTITIES} == {"definition"}
    assert {r["quantity"] for r in allrows if r["status"] == "proposal"} >= {"domain_n"}
    assert {r["status"] for r in allrows if r["quantity"] == "alongline_S"} == {"definition"}      # D2 section 16
    # a bulk control with too few solute atoms (DEFINITIONS 5, 2b): x_bc only as a diagnostic, x_iso stays the reference
    defs = load_definitions()
    params = {**defs.params, "reference": {**defs.params["reference"], "min_solute_in_bulk_control": 10 ** 6}}
    strict = tmp_path / "strict.yaml"
    strict.write_text(yaml.safe_dump({**params, "label": "D2"}))
    out2 = tmp_path / "out2"
    assert main(["analyse", str(fd), "--out", str(out2), "--bulk-dir", str(out), "--cache", str(tmp_path / "cache"),
                 "--definitions", str(strict)]) == 0
    prov2 = json.loads((out2 / "RUN_stage3_D2.csv.provenance.json").read_text())
    assert prov2["x_bc"]["enough_solute"] is False and "no bulk control" in prov2["reference_rule"]
    assert prov2["stage3"]["diagnostic_iv"].startswith("not usable")   # a bulk control (too few solute), but smaller
    # swap rate 20 % off (D2 section 15): not accepted, x_iso (and c_B) as without a bulk control, reason written
    fd_fast = tmp_path / "run_fast.yaml"
    fd_fast.write_text(yaml.safe_dump({**md, "swap_attempts_per_atom_per_step": 0.0125}))
    out4 = tmp_path / "out4"
    assert main(["analyse", str(fd_fast), "--out", str(out4), "--bulk-dir", str(out), "--cache", str(tmp_path / "cache")]) == 0
    prov4 = json.loads((out4 / "RUN_stage3_D2.csv.provenance.json").read_text())
    assert prov4["x_bc"]["accepted"] is False and "swap rate differs" in prov4["x_bc"]["reasons"][0]
    assert "no bulk control" in prov4["reference_rule"] and "swap rate differs" in prov4["reference_rule"]
    assert prov4["x_ref"]["x_iso"] == 0.05 and "x_bc" not in prov4["x_ref"]
    # review finding 1: the rejected x_bc is only a diagnostic -- reference label x_bc_diagnostic, never headline, no Stage 4
    refs4 = {r["ref"] for r in read_table(out4 / "RUN_stage2_D2.csv") if r["quantity"] == "Gamma"}
    assert refs4 == {"c_B", "x_iso", "x_bc_diagnostic"}
    assert main(["table", str(out4), "--out", str(tmp_path / "all4.csv")]) == 0
    g20 = {r["ref"]: r["role"] for r in read_table(tmp_path / "all4.csv")
           if r["quantity"] == "Gamma" and r["region"] == "20" and r["basis"] == "r"}
    assert g20 == {"c_B": "headline", "x_iso": "headline", "x_bc_diagnostic": "diagnostic"}
    assert main(["across", str(out4), "--out", str(tmp_path / "s4b.csv")]) == 0
    assert not [r for r in read_table(tmp_path / "s4b.csv") if "x_bc" in r["ref"]]
    # Stage 3 carries the section-16 window test, not window means of per-snapshot peak values
    q3 = {r["quantity"] for r in read_table(out4 / "RUN_stage3_D2.csv")}
    assert "alongline_window_periodic" in q3 and "alongline_p_peak" not in q3
    # D1 (bulk control = same size): the smaller cell is not a bulk control, x_bc only as a diagnostic
    out3 = tmp_path / "out3"
    assert main(["analyse", str(fd), str(fb), "--out", str(out3), "--cache", str(tmp_path / "cache"),
                 "--definitions", "D1"]) == 0
    prov3 = json.loads((out3 / "RUN_stage3_D1.csv.provenance.json").read_text())
    assert prov3["x_bc"]["accepted"] is False and "no bulk control" in prov3["reference_rule"]
    assert not any(r["quantity"].startswith("alongline") for r in read_table(out3 / "RUN_stage2_D1.csv"))


def test_missing_bulk_control_table_is_reported(tmp_path, capsys):
    cell, pos, (y0, z0) = screw_dislocation_cell()
    rdir = tmp_path / "PROJ" / "T" / "RUN"
    rdir.mkdir(parents=True)
    write_lammps_dump(rdir / "screw_quenched_step1000.dump", cell, np.arange(1, len(pos) + 1), np.ones(len(pos), int),
                      pos, timestep=1000, pbc=("pp", "ss", "ss"))
    md = {**base_manifest("RUN", "PROJ/T/RUN/*_quenched_step*.dump"), "character": "d0", "potential": "CuNi-ACE",
          "frame": {"X": [1, -1, 0], "Y": [1, 1, -2], "Z": [1, 1, 1]}, "frame_source": "t", "z_gp": float(z0),
          "z_gp_source": "t", "bulk_control": "BULK", "x_ref_source": "bulk_control", "start_type": "bottom-up",
          "swap_attempts_per_atom_per_step": 0.01}
    fd = tmp_path / "run.yaml"
    fd.write_text(yaml.safe_dump(md))
    out = tmp_path / "out"
    assert main(["analyse", str(fd), "--out", str(out), "--cache", str(tmp_path / "cache")]) == 0
    assert "analyse it first" in capsys.readouterr().out
    prov = json.loads((out / "RUN_stage2_D2.csv.provenance.json").read_text())
    # ACE Ni-Cu that names a bulk control without having its x_bc yet stays held back (v0.2)
    assert prov["x_bc"] is None and prov["paper_flag"] == "no paper number before bulk control"


def test_pbc_override_for_lammps_data_file(tmp_path):
    cell = Cell(np.diag([6 * A, 6 * A, 8 * A]), pbc=(True, True, False))
    pos = fcc_in_cell(cell, A) + 1e-3
    f = tmp_path / "P" / "T" / "R" / "cell.lmps"
    f.parent.mkdir(parents=True)
    lines = ["LAMMPS data file", "", f"{len(pos)} atoms", "2 atom types", "", f"0 {6 * A} xlo xhi", f"0 {6 * A} ylo yhi",
             f"0 {8 * A} zlo zhi", "", "Masses", "", "1 58.69", "2 63.55", "", "Atoms # atomic", ""]
    lines += [f"{i + 1} 1 {x:.6f} {y:.6f} {z:.6f}" for i, (x, y, z) in enumerate(pos)]
    f.write_text("\n".join(lines) + "\n")
    d = load_definitions()
    s_default = get_snapshot(f, d, cache_root=tmp_path / "c")
    s_over = get_snapshot(f, d, cache_root=tmp_path / "c", pbc=(True, True, False))
    assert bool(s_default.pbc[2]) is True and bool(s_over.pbc[2]) is False
    assert s_over.meta["pbc_override"] == [True, True, False] and s_default.meta["pbc_override"] is None
    # periodic in z, the perfect slab is fcc everywhere; with free z faces the surface atoms are not fcc
    assert (s_default.acna == 1).all() and (s_over.acna != 1).any()
