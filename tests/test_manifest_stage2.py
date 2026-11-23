"""Stage 0 manifest validation, and Stage 2 end to end on a synthetic snapshot."""

import copy
from pathlib import Path

import numpy as np
import pytest
import yaml

from disegreg import manifest
from disegreg.cache import STRUCTURE_CODES, Snapshot
from disegreg.definitions import load_definitions
from disegreg.geometry import Cell
from disegreg.stage2 import analyse_snapshot
from helpers import fcc_in_cell

EXAMPLES = Path(__file__).resolve().parents[1] / "examples" / "manifests"


def demo():
    with open(next(EXAMPLES.glob("*.yaml"))) as f:
        return yaml.safe_load(f)


def test_example_manifests_are_valid():
    files = sorted(EXAMPLES.glob("*.yaml"))
    assert files
    for f in files:
        m = manifest.load_manifest(f)
        assert m.frame_matrix() @ m.frame_matrix().T == pytest.approx(np.eye(3))


@pytest.mark.parametrize("change, message", [
    ({"delta": -0.5}, "delta"),
    ({"frame": {"X": [1, 0, 0], "Y": [0, 1, 0], "Z": [1, 1, 1]}}, "orthogonal"),
    ({"frame": {"X": [0, 1, -1], "Y": [2, -1, -1], "Z": [1, 1, 1]}}, "right-handed"),
    ({"frame": {"X": [1, -1, 0], "Y": [0, 0, 1], "Z": [1, 1, 0]}}, "<111>"),
    ({"quenched": False}, "quenched"),
    ({"character": "d45"}, "character"),
    ({"bulk_control_reason": None}, "bulk_control_reason"),
    ({"dmu_lammps_sign": -1}, "dmu_lammps_sign"),
    ({"x_ref_source": "bulk_control"}, "bulk_control run id"),
    ({"start_type": "seeded"}, "seeded_from"),
    ({"temperature": True}, "temperature"),
    ({"surprise": 1}, "unknown field"),
])
def test_manifest_errors(change, message):
    d = demo()
    for k, v in change.items():
        if v is None:
            d.pop(k)
        else:
            d[k] = v
    errs = manifest.validate(d)
    assert any(message in e for e in errs), errs


def test_missing_and_template(tmp_path):
    d = demo()
    d.pop("z_gp")
    assert "missing field 'z_gp'" in manifest.validate(d)
    t = yaml.safe_load(manifest.template())
    assert set(t) == set(manifest.FIELDS)
    p = tmp_path / "m.yaml"
    p.write_text(yaml.safe_dump(d))
    with pytest.raises(manifest.ManifestError):
        manifest.load_manifest(p)
    assert manifest.load_manifest(p, strict=False).data["run_id"] == d["run_id"]


def test_dump_files_relative_and_sorted_by_step(tmp_path):
    d = demo()
    d["dumps"] = "run/*_quenched_step*.dump"
    (tmp_path / "run").mkdir()
    for s in (1000000, 99000, 830000):
        (tmp_path / "run" / f"x_quenched_step{s}.dump").write_text("")
    p = tmp_path / "m.yaml"
    p.write_text(yaml.safe_dump(d))
    m = manifest.load_manifest(p)
    assert [f.name for f in m.dump_files()] == [f"x_quenched_step{s}.dump" for s in (99000, 830000, 1000000)]
    assert m.kT == pytest.approx(manifest.K_B * 300)


# ------------------------------------------------------------------------------------------------- Stage 2, synthetic
A = 3.52


def synthetic_snapshot(seed=0):
    """Cubic fcc slab, X = [100] periodic line direction, Y periodic, Z free; one straight 'DXA' Shockley line at
    (y0, z0); a-CNA 'other' within 3 A of the line, hcp in a thin ribbon at z0, fcc elsewhere. All atoms with r < 5 A
    are solute; elsewhere solute with probability 0.01."""
    cell = Cell(np.diag([4 * A, 12 * A, 16 * A]), pbc=(True, True, False))
    pos = (fcc_in_cell(cell, A) + 1e-3).astype(np.float32).astype(float)  # as stored in the cache
    y0, z0 = 6 * A + 0.3, 8 * A + 0.2
    dy = pos[:, 1] - y0
    dy -= cell.H[1, 1] * np.round(dy / cell.H[1, 1])  # minimum image in Y
    r = np.hypot(dy, pos[:, 2] - z0)
    acna = np.full(len(pos), STRUCTURE_CODES["fcc"], np.int8)
    acna[(np.abs(pos[:, 2] - z0) < 1.0) & (r >= 3)] = STRUCTURE_CODES["hcp"]
    acna[r < 3] = STRUCTURE_CODES["other"]
    rng = np.random.default_rng(seed)
    types = np.where((r < 5) | (rng.random(len(pos)) < 0.01), 2, 1).astype(np.int8)
    line = np.array([[0.0, y0, z0], [cell.line_length, y0, z0]])
    snap = Snapshot(
        positions=pos.astype(np.float32), types=types, cell_matrix=np.column_stack([cell.H, cell.origin]),
        pbc=np.array(cell.pbc), acna=acna, ptm=acna.copy(), seg_points=line, seg_offsets=np.array([0, 2]),
        seg_true_b=np.array([[1 / 6, 1 / 6, -1 / 3]]), seg_spatial_b=np.zeros((1, 3)),
        seg_length=np.array([cell.line_length]), seg_family=np.array(["1/6<112>"]), seg_cluster=np.array([1]),
        seg_infinite=np.array([True]), seg_loop=np.array([False]))
    return snap, cell, pos, r, types == 2, z0


def test_stage2_synthetic_end_to_end():
    snap, cell, pos, r_true, solute, z0 = synthetic_snapshot()
    d = load_definitions()
    x_ref = 0.01
    res = analyse_snapshot(snap, d, A, z_gp=z0 + 0.5, x_ref=x_ref, keep_arrays=True)
    assert res.info["mode"] == "polyline" and res.info["f_cov"] == 1.0 and res.info["f_shockley"] == 1.0
    assert res.info["z_gp_lines"] == pytest.approx(z0) and not res.info["z_gp_flag"]
    far = r_true >= 0.5  # polyline densified to 0.25 A: error < 0.02 A at r >= 0.5 A (DEFINITIONS 3.1)
    assert res.r_dxa[far] == pytest.approx(r_true[far], abs=0.02)
    bound = np.sqrt(r_true ** 2 + 0.125 ** 2) - r_true + 1e-4  # point at most half a spacing away along the line
    assert np.all(np.abs(res.r_dxa - r_true) <= bound)

    z = pos[:, 2]
    interior = (z - z.min() >= 15) & (z.max() - z >= 15)
    L = cell.line_length
    assert res.value("L_x") == pytest.approx(L)
    # every interior atom has one code; X are the excluded layers
    n_total = sum(res.value("N", region=c) for c in "XCFSTBO")
    assert n_total == len(pos) and res.value("N", region="X") == (~interior).sum()
    core = interior & (r_true < 5)
    assert res.value("N", region="C") == core.sum() and res.value("c", region="C") == 1.0
    fault = interior & (snap.acna == STRUCTURE_CODES["hcp"]) & (np.abs(z - (z0 + 0.5)) < 4) & (r_true >= 5)
    assert res.value("N", region="F") == fault.sum() > 0
    for R in (5, 20, 40):
        m = interior & (r_true < R)
        g = ((m & solute).sum() - x_ref * m.sum()) / L
        assert res.value("Gamma", region=f"{R}", basis="r") == pytest.approx(g)
        assert res.value("Gamma", region=f"{R}", basis="r_dxa") == pytest.approx(g)
    assert res.value("E", region="C") == pytest.approx(100.0)
    assert res.value("partial_separation_valid_bins") == 0  # one Shockley line only
    # the core is one connected solute cluster along the whole periodic line
    assert res.value("cluster_max_size", region="C") == (r_true < 5).sum()
    assert np.isnan(res.value("alpha1_local", region="C"))  # c_C = 1: alpha undefined


def test_stage2_without_x_ref_and_without_z_gp():
    snap, *_ = synthetic_snapshot()
    res = analyse_snapshot(snap, load_definitions(), A)
    q = {r["quantity"] for r in res.records}
    assert "Gamma" not in q and "E" not in q and "alpha1_xref" not in q and "N_sol_inside" in q
    assert res.info["z_gp_source"] == "lines" and res.info["z_gp_flag"]


def test_stage2_flags_wrong_z_gp():
    snap, *_, z0 = synthetic_snapshot()
    res = analyse_snapshot(snap, load_definitions(), A, z_gp=z0 + 3.0)
    assert res.info["z_gp_flag"] and res.value("z_gp_difference") == pytest.approx(3.0)


def test_results_table_refuses_mixed_definition_sets(tmp_path):
    import json

    from disegreg.results import STAGE3_COLUMNS, collect_runs, write_table

    files = []
    for label in ("D1", "D2"):
        row = {"run_id": f"run{label}", "quantity": "c", "region": "C", "basis": "", "mean": 0.1, "sem": 0.01,
               "n_samples": 3, "sem_plateau": False, "change": 0.0, "change_se": 0.01, "flat": True,
               "flat_change_se": True, "t_first": 0.0, "t_last": 1.0, "definition_set": label, "tool_version": "x"}
        f = tmp_path / f"run{label}_stage3_{label}.csv"
        write_table(f, [row], {"definition_set": label, "manifest_data": demo(), "stage3": {"converged": True}},
                    columns=STAGE3_COLUMNS)
        files.append(f)
    rows, prov = collect_runs(files[:1])
    assert rows[0]["value"] == 0.1 and rows[0]["error"] == 0.01 and rows[0]["solute"] == "Cu"
    assert json.dumps(prov)
    with pytest.raises(ValueError, match="different definition sets"):
        collect_runs(files)
