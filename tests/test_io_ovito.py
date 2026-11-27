"""Stage 1 with OVITO on small synthetic dumps (skipped without OVITO)."""

import numpy as np
import pytest

pytest.importorskip("ovito")
pytestmark = pytest.mark.ovito

from helpers import fcc_in_cell, screw_dislocation_cell, write_lammps_dump  # noqa: E402

from disegreg import cache as cache_mod  # noqa: E402
from disegreg import io_ovito  # noqa: E402
from disegreg.cache import Snapshot, get_snapshot  # noqa: E402
from disegreg.definitions import load_definitions  # noqa: E402
from disegreg.geometry import Cell, line_distances  # noqa: E402

D1 = load_definitions("D1")
A = 3.52


@pytest.fixture
def tilted_fcc_dump(tmp_path):
    """Perfect fcc in a tilted cell (xy = a keeps the cell vectors lattice vectors), shuffled ids, two types, gzipped."""
    cell = Cell.lammps(8 * A, 8 * A, 8 * A, xy=A, origin=(-1.0, 2.0, 0.5), pbc=(True, True, False))
    pos = fcc_in_cell(cell, A)
    rng = np.random.default_rng(3)
    ids = rng.permutation(len(pos)) + 1
    types = np.where(rng.random(len(pos)) < 0.1, 2, 1)
    path = tmp_path / "PROJ" / "THREAD" / "RUN" / "snap_quenched_step1000.dump.gz"
    path.parent.mkdir(parents=True)
    write_lammps_dump(path, cell, ids, types, pos)
    return path, cell, ids, types, pos


def test_modifier_configuration():
    acna, ptm, dxa = io_ovito.make_modifiers(D1["structure"])
    assert ptm.rmsd_cutoff == pytest.approx(0.12)
    assert {s.name for s in ptm.structures if s.enabled} == {"Other", "FCC", "HCP", "BCC", "ICO"}
    assert dxa.trial_circuit_length == 14 and dxa.circuit_stretchability == 9 and dxa.line_coarsening_enabled


def test_stage1_tilted_perfect_fcc(tilted_fcc_dump):
    path, cell, ids, types, pos = tilted_fcc_dump
    s = io_ovito.compute_stage1(path, D1)
    order = np.argsort(ids)
    assert s.ids is None  # ids are 1..N after sorting
    assert np.allclose(s.positions, pos[order], atol=1e-4)
    assert np.array_equal(s.types, types[order])
    assert np.allclose(s.cell_matrix[:, :3], cell.H, atol=1e-6) and np.allclose(s.cell_matrix[:, 3], cell.origin)
    assert s.pbc.tolist() == [True, True, False]
    interior = (pos[order, 2] > cell.origin[2] + 3) & (pos[order, 2] < cell.origin[2] + cell.c[2] - 3)
    assert np.all(s.acna[interior] == 1) and np.all(s.ptm[interior] == 1)  # free z surfaces are not fcc
    assert s.n_segments == 0
    assert s.meta["timestep"] == 1000 and s.meta["provenance"]["ovito"] == "3.17.2"
    assert s.meta["stage1_params"]["definition_set"] == "D1"


def test_stage1_screw_dislocation_found_and_distances(tmp_path):
    cell, pos, (y0, z0) = screw_dislocation_cell()
    path = tmp_path / "screw.dump"
    write_lammps_dump(path, cell, np.arange(1, len(pos) + 1), np.ones(len(pos), int), pos, pbc=("pp", "ss", "ss"))
    s = io_ovito.compute_stage1(path, D1)
    assert s.n_segments == 1
    assert s.seg_family[0] == "1/2<110>" and s.seg_infinite[0]
    assert s.seg_length[0] == pytest.approx(cell.line_length, rel=0.01)
    assert np.allclose(np.abs(s.seg_spatial_b[0]), [A / np.sqrt(2), 0, 0], atol=0.02)
    (line,) = s.polylines()
    assert np.allclose(line[:, 1:].mean(axis=0), [y0, z0], atol=1.0)
    p = s.positions.astype(float)
    interior = (p[:, 1] > 5) & (p[:, 1] < cell.H[1, 1] - 25) & (p[:, 2] > 5) & (p[:, 2] < cell.H[2, 2] - 25)
    res = line_distances(s.cell(), p, s.polylines(), interior & (s.acna != 1), s.acna == 0, A)
    assert res.mode == "polyline" and res.f_cov == pytest.approx(1.0)
    r_ideal = np.hypot(p[:, 1] - y0, p[:, 2] - z0)
    assert np.max(np.abs(res.r - r_ideal)[r_ideal > 3]) < 1.2  # DXA line within ~1 A of the elastic centre


def test_get_snapshot_cache_hit_and_stale(tilted_fcc_dump, tmp_path, monkeypatch):
    path = tilted_fcc_dump[0]
    root = tmp_path / "CACHE"
    s1 = get_snapshot(path, D1, cache_root=root)
    files = list(root.rglob("*.npz"))
    assert len(files) == 1 and files[0].parts[-4:-1] == ("PROJ", "THREAD", "RUN")
    calls = []
    real = io_ovito.compute_stage1
    monkeypatch.setattr(io_ovito, "compute_stage1", lambda *a, **k: calls.append(1) or real(*a, **k))
    s2 = get_snapshot(path, D1, cache_root=root)
    assert calls == [] and np.array_equal(s1.positions, s2.positions)  # cache hit: OVITO not run
    import gzip
    data = gzip.open(path, "rt").read().replace("ITEM: TIMESTEP\n1000", "ITEM: TIMESTEP\n1001")  # changed dump
    with gzip.open(path, "wt") as f:
        f.write(data)
    s3 = get_snapshot(path, D1, cache_root=root)
    assert calls == [1] and s3.meta["timestep"] == 1001
    assert isinstance(Snapshot.load(files[0]), Snapshot)
    assert cache_mod.CACHE_FORMAT == 1


def test_cli_analyse_end_to_end(tmp_path):
    """manifest -> Stage 1 (OVITO) -> Stage 2 -> tidy CSV + provenance, through the command-line entry point."""
    import json

    import yaml

    from disegreg.cli import main
    from disegreg.results import read_table

    cell, pos, (y0, z0) = screw_dislocation_cell()
    run = tmp_path / "PROJ" / "THREAD" / "RUN"
    run.mkdir(parents=True)
    rng = np.random.default_rng(5)
    types = np.where(rng.random(len(pos)) < 0.05, 2, 1)
    for step in (2000, 1000):
        write_lammps_dump(run / f"screw_quenched_step{step}.dump", cell, np.arange(1, len(pos) + 1), types, pos,
                          timestep=step, pbc=("pp", "ss", "ss"))
    m = {"schema_version": 1, "project": "PROJ", "run_id": "RUN", "potential": "none", "host_element": "Ni",
         "solute_element": "Cu", "character": "d0", "frame": {"X": [1, -1, 0], "Y": [1, 1, -2], "Z": [1, 1, 1]},
         "frame_source": "test", "a_cell": A, "a_cell_source": "test", "temperature": 300, "dmu": 0.5,
         "dmu_lammps": 0.5, "dmu_lammps_sign": 1, "dmu_coex": 0.5, "coex_line": "incoherent", "dmu_coex_source": "test",
         "delta": 0.0, "dumps": "PROJ/THREAD/RUN/*_quenched_step*.dump", "quenched": True, "z_boundary": "free",
         "z_gp": float(z0), "z_gp_source": "test", "bulk_control": None, "bulk_control_reason": "test",
         "x_ref_source": "isotherm", "x_iso": 0.05, "x_iso_source": "test", "start_type": "bottom-up"}
    mf = tmp_path / "manifest.yaml"
    mf.write_text(yaml.safe_dump(m))
    out = tmp_path / "out"
    assert main(["manifest", "check", str(mf)]) == 0
    # the v0.1 path, explicitly with D1 (default since v0.2: D2, which gives the same D1 quantities)
    assert main(["analyse", str(mf), "--out", str(out), "--cache", str(tmp_path / "cache"), "--definitions", "D1"]) == 0
    rows = read_table(out / "RUN_stage2_D1.csv")
    assert [r["timestep"] for r in rows if r["quantity"] == "f_cov"] == [1000, 2000]  # sorted by step
    n = {r["region"]: r["value"] for r in rows if r["quantity"] == "N" and r["timestep"] == 1000}
    assert sum(n[c] for c in "XCFSTBO") == len(pos)
    assert all(r["definition_set"] == "D1" for r in rows)
    prov = json.loads((out / "RUN_stage2_D1.csv.provenance.json").read_text())
    assert prov["x_ref"] == {"x_iso": 0.05} and prov["x_ref_source"] == "isotherm" and len(prov["snapshots"]) == 2
    assert "both c_B and x_iso" in prov["reference_rule"] and prov["paper_flag"] == ""
    assert {r["ref"] for r in rows if r["quantity"] == "Gamma"} == {"c_B", "x_iso"}  # DEFINITIONS 5, rev. 4
    assert prov["snapshots"][0]["ovito"] == io_ovito.ovito_version()
    # Stage 3 was written next to Stage 2, and `table` merges it with the manifest metadata
    assert (out / "RUN_stage3_D1.csv").is_file()
    assert main(["table", str(out), "--out", str(tmp_path / "all.csv")]) == 0
    allrows = read_table(tmp_path / "all.csv")
    assert {r["run_id"] for r in allrows} == {"RUN"} and {r["character"] for r in allrows} == {"d0"}
    assert {r["diagnostic_iv"] for r in allrows} == {""}       # D1, Delta = 0, no bulk control (D1 section 5)
    head = [r for r in allrows if r["role"] == "headline"]
    assert {(r["quantity"], r["region"], r["basis"]) for r in head} == {("Gamma", "20", "r")}
    assert {r["ref"] for r in head} == {"c_B", "x_iso"} and {r["drift_criterion"] for r in head} == {"change_se"}
    assert {r["region"] for r in allrows if r["role"] == "check"} == {"40"}
    # the Gamma(R) curve of the run, 5-60 A (and the bin edges 80, 100), under r and r_dxa, for both references
    curve = read_table(out / "RUN_gammaR_D1.csv")
    radii = sorted({float(r["R"]) for r in curve})
    assert radii[:12] == [5.0 * k for k in range(1, 13)] and {r["ref"] for r in curve} == {"c_B", "x_iso"}
