"""Stage 4 (across runs), the McLean site tables and the bulk-control manifest -- synthetic, exact cases, no OVITO."""

import copy
import json

import numpy as np
import pytest
import yaml

from disegreg import across, sites
from disegreg.cache import STRUCTURE_CODES, Snapshot
from disegreg.cli import find_x_bc
from disegreg.definitions import load_definitions, load_proposals
from disegreg.geometry import Cell
from disegreg.manifest import validate
from disegreg.results import (PAPER_HOLD, STAGE3_COLUMNS, ProposalAsHeadline, check_headline_status, collect_runs,
                              read_table, write_table)
from disegreg.stage2 import BULK_REGION, analyse_bulk_snapshot
from helpers import fcc_in_cell
from test_manifest_stage2 import demo

D2 = load_definitions("D2")
PROPS = load_proposals()
KB = 8.617333262e-5


# ------------------------------------------------------------------------------------------------------ helpers
def s3row(q, region, mean, sem, basis="", ref=""):
    return {"run_id": "", "quantity": q, "region": region, "basis": basis, "ref": ref, "mean": mean, "sem": sem,
            "n_samples": 5, "sem_plateau": True, "change": 0.0, "change_se": sem, "flat": True, "flat_rev3_sem": True,
            "t_first": 0.0, "t_last": 1.0, "definition_set": "D1", "tool_version": "test"}


def write_run(tmp_path, run_id, delta, c_C, c_B, x_iso, sem=0.0, start="bottom-up", gamma20=None, converged=True,
              potential="LeeShim2004-MEAM", x_bc=None, paper_flag=""):
    man = demo()
    man.update(run_id=run_id, delta=delta, x_iso=x_iso, start_type=start, potential=potential)
    if start == "seeded":
        man["seeded_from"] = "x"
    refs = {"c_B": c_B, "x_iso": x_iso}
    if x_bc is not None:
        refs["x_bc"] = x_bc[0]
        man["bulk_control"] = "bulk"
    rows = [s3row("c", "C", c_C, sem * c_C), s3row("c", "B", c_B, sem * c_B), s3row("L_x", "", 100.0, 0.0),
            s3row("N_inside", "20", 18000.0, 10.0, "r"), s3row("N_inside", "40", 70000.0, 10.0, "r")]
    for ref, x in refs.items():
        rows.append(s3row("x_ref", "", x, 0.0, ref=ref))
        g = gamma20 if gamma20 is not None else (c_C - x) * 50
        rows.append(s3row("Gamma", "20", g, 0.1, "r", ref))
        rows.append(s3row("Gamma", "40", g + 1, 0.2, "r", ref))
    for r in rows:
        r["run_id"] = run_id
    f = tmp_path / f"{run_id}_stage3_D1.csv"
    prov = {"definition_set": "D1", "manifest_data": man,
            "stage3": {"converged": converged, "paper_flag": paper_flag},
            "x_bc": None if x_bc is None else {"value": x_bc[0], "sem": x_bc[1]}}
    write_table(f, rows, prov, columns=STAGE3_COLUMNS)
    return f


# --------------------------------------------------------------------------------------------------- uptake
def test_uptake_exponent_exact_power_law(tmp_path):
    """c_C = A x^n exactly -> n_up = n between every pair of offsets, for every reference proportional to x."""
    n_true = 3.5
    files = []
    for k, delta in enumerate((-1.0, -0.5, -0.25, -0.1)):
        x = 1e-3 * np.exp(delta)                         # dilute: x proportional to exp(Delta)
        files.append(write_run(tmp_path, f"r{k}", delta, c_C=2e4 * x ** n_true, c_B=0.9 * x, x_iso=x))
    runs, label = across.load_runs(files)
    assert label == "D1"
    rows = across.uptake_rows(next(iter(across.run_sets(runs))), across.run_sets(runs)[next(iter(across.run_sets(runs)))])
    assert len(rows) == 3 * 2                            # 3 neighbouring pairs x 2 references (c_B, x_iso)
    for r in rows:
        assert r["value"] == pytest.approx(n_true, rel=1e-12) and r["status"] == "definition"
        assert r["error"] == 0.0 and r["family"] == "up" and r["delta"] < r["delta_2"]


def test_uptake_error_propagation_matches_monte_carlo():
    c1, s1, c2, s2, x1, sx1, x2, sx2 = 0.10, 0.004, 0.35, 0.01, 1e-3, 2e-5, 2e-3, 3e-5
    n, e = across.uptake_exponent(c1, s1, c2, s2, x1, sx1, x2, sx2)
    assert n == pytest.approx(np.log(3.5) / np.log(2))
    rng = np.random.default_rng(1)
    k = 200_000
    ns = (np.log(rng.normal(c2, s2, k)) - np.log(rng.normal(c1, s1, k))) / (
        np.log(rng.normal(x2, sx2, k)) - np.log(rng.normal(x1, sx1, k)))
    assert e == pytest.approx(ns.std(), rel=0.03)
    assert np.isnan(across.uptake_exponent(0.0, 0, 0.1, 0, 1e-3, 0, 2e-3, 0)[0])   # c = 0: undefined
    assert np.isnan(across.uptake_exponent(0.1, 0, 0.2, 0, 1e-3, 0, 1e-3, 0)[0])   # same x: undefined
    n, e = across.uptake_exponent(0.1, np.nan, 0.2, 0.01, 1e-3, 0, 2e-3, 0)          # one run: single snapshot
    assert n == pytest.approx(1.0) and np.isnan(e)                                   # error unknown, not 0


def test_run_sets_families_flags_and_paper_flag(tmp_path):
    files = [write_run(tmp_path, "a", -0.5, 0.1, 1e-3, 1e-3),
             write_run(tmp_path, "b", -0.25, 0.3, 2e-3, 2e-3, start="seeded", converged=False),
             write_run(tmp_path, "c", -0.25, 0.35, 2e-3, 2e-3, start="top-down"),
             write_run(tmp_path, "d", -0.1, 0.6, 3e-3, 3e-3, start="top-down"),
             write_run(tmp_path, "e", -0.5, 0.1, 1e-3, 1e-3, potential="CuNi-ACE", paper_flag=PAPER_HOLD),
             write_run(tmp_path, "f", -0.25, 0.2, 2e-3, 2e-3, potential="CuNi-ACE", paper_flag=PAPER_HOLD)]
    runs, _ = across.load_runs(files)
    sets = across.run_sets(runs)
    assert len(sets) == 2
    rows = [r for k, rs in sets.items() for r in across.uptake_rows(k, rs)]
    up = [r for r in rows if r["potential"] == "LeeShim2004-MEAM" and r["family"] == "up"]
    assert {r["runs"] for r in up} == {"a|b"} and all("not converged: b" in r["flags"] for r in up)
    assert {r["runs"] for r in rows if r["family"] == "down"} == {"c|d"}
    ace = [r for r in rows if r["potential"] == "CuNi-ACE"]
    assert ace and all(r["paper_flag"] == PAPER_HOLD for r in ace)       # ACE Ni-Cu keeps its paper flag


# ----------------------------------------------------------------------------------------------- hysteresis
def test_hysteresis_pairs_and_verdict(tmp_path):
    files = [write_run(tmp_path, "up", -0.25, 0.30, 2e-3, 2e-3, gamma20=4.0),
             write_run(tmp_path, "down", -0.25, 0.30, 2e-3, 2e-3, start="top-down", gamma20=4.5),
             write_run(tmp_path, "up2", -0.5, 0.30, 2e-3, 2e-3, gamma20=4.0)]
    runs, _ = across.load_runs(files)
    (key, rs), = across.run_sets(runs).items()
    rows = across.hysteresis_rows(key, rs, D2["hysteresis"])
    assert {r["runs"] for r in rows} == {"up|down"}                       # same Delta only
    g = {(r["quantity"], r["ref"]): r for r in rows if r["region"] == "20"}
    d = g[("Gamma_up_minus_down", "x_iso")]
    assert d["value"] == pytest.approx(-0.5) and d["error"] == pytest.approx(0.1 * np.sqrt(2))
    v = g[("Gamma_hysteresis", "x_iso")]
    assert v["value"] == 1.0 and v["status"] == "definition" and v["error"] == pytest.approx(0.5 / (0.1 * np.sqrt(2)))
    c = [r for r in rows if r["quantity"] == "c_hysteresis"]
    assert c and c[0]["value"] == 0.0                                     # equal c_C: not different


# ---------------------------------------------------------------------------------- headline / status guard
def test_results_table_refuses_proposal_as_headline(tmp_path):
    rows = [{"role": "headline", "status": "definition"}, {"role": "curve", "status": "proposal"}]
    check_headline_status(rows)
    with pytest.raises(ProposalAsHeadline):
        check_headline_status(rows + [{"role": "headline", "status": "proposal", "quantity": "Gamma"}])
    f = write_run(tmp_path, "a", -0.5, 0.1, 1e-3, 1e-3)
    rows, _ = collect_runs([f])
    assert {r["status"] for r in rows} == {"definition"}
    from disegreg import results
    results.PROPOSAL_QUANTITIES.add("Gamma")                              # pretend Gamma were undecided
    try:
        with pytest.raises(ProposalAsHeadline):
            collect_runs([f])
    finally:
        results.PROPOSAL_QUANTITIES.discard("Gamma")


def test_stage4_driver_and_x_bc_error(tmp_path):
    files = [write_run(tmp_path, "a", -0.5, 0.1, 1e-3, 1e-3, x_bc=(1.1e-3, 2e-5)),
             write_run(tmp_path, "b", -0.25, 0.3, 2e-3, 2e-3, x_bc=(2.1e-3, 3e-5))]
    runs, _ = across.load_runs(files)
    rows = across.stage4(runs, D2)
    xr = {r["ref"] for r in rows if r["analysis"] == "uptake"}
    assert xr == {"c_B", "x_iso", "x_bc"}
    up = next(r for r in rows if r["analysis"] == "uptake" and r["ref"] == "x_bc")
    n, e = across.uptake_exponent(0.1, 0.0, 0.3, 0.0, 1.1e-3, 2e-5, 2.1e-3, 3e-5)
    assert up["value"] == pytest.approx(n) and up["error"] == pytest.approx(e) and e > 0
    xe = next(r for r in rows if r["analysis"] == "x_bc_error" and r["runs"] == "a" and r["region"] == "20")
    assert xe["error"] == pytest.approx(18000 / 100 * 2e-5)
    assert xe["value"] == pytest.approx(np.hypot(0.1, 18000 / 100 * 2e-5)) and xe["status"] == "definition"
    assert {r["status"] for r in rows} == {"definition"}                  # D2: no Stage-4 row depends on a proposal
    # tables of D1: only what D1 defines (uptake exponents, hysteresis differences); no x_bc error, no verdict
    d1_rows = across.stage4(runs, load_definitions("D1"))
    assert {r["analysis"] for r in d1_rows} == {"uptake"}
    assert [r["value"] for r in d1_rows] == [r["value"] for r in rows if r["analysis"] == "uptake"]


def test_hysteresis_with_D1_gives_differences_only(tmp_path):
    files = [write_run(tmp_path, "up", -0.25, 0.30, 2e-3, 2e-3, gamma20=4.0),
             write_run(tmp_path, "down", -0.25, 0.30, 2e-3, 2e-3, start="top-down", gamma20=4.5)]
    runs, _ = across.load_runs(files)
    rows = [r for r in across.stage4(runs, load_definitions("D1")) if r["analysis"] == "hysteresis"]
    assert rows and all(r["quantity"].endswith("_up_minus_down") for r in rows)
    rows2 = [r for r in across.stage4(runs, D2) if r["analysis"] == "hysteresis"]
    assert {r["quantity"] for r in rows2} == {"Gamma_up_minus_down", "Gamma_hysteresis", "c_up_minus_down", "c_hysteresis"}
    assert ([r["value"] for r in rows2 if r["quantity"].endswith("_up_minus_down")] == [r["value"] for r in rows])


# --------------------------------------------------------------------------------------------------- McLean
def test_mclean_theta_exact():
    kT = KB * 300
    assert sites.mclean_theta([0.0], 0.01, kT)[0] == pytest.approx(0.01)
    E = -0.1
    q = 0.01 / 0.99 * np.exp(-E / kT)
    assert sites.mclean_theta([E], 0.01, kT)[0] == pytest.approx(q / (1 + q))
    assert sites.mclean_theta([-50.0], 1e-6, kT)[0] == pytest.approx(1.0)      # no overflow
    assert sites.mclean_theta([+50.0], 1e-6, kT)[0] == pytest.approx(0.0)
    p = sites.mclean_prediction(np.array([-0.2, -0.1, 0.0, 0.05]), 1e-3, kT)
    th = sites.mclean_theta([-0.2, -0.1, 0.0, 0.05], 1e-3, kT)
    assert p.c_C == pytest.approx(th.mean()) and p.n_sites == 4 and p.n_up <= 1.0
    h = 1e-9                                                  # derivative and exponent against finite differences
    cp = sites.mclean_prediction(np.array([-0.2, -0.1, 0.0, 0.05]), 1e-3 + h, kT).c_C
    assert p.dc_dx == pytest.approx((cp - p.c_C) / h, rel=1e-4)
    assert sites.mclean_prediction(np.array([0.0, 0.0]), 0.02, kT).n_up == pytest.approx(1.0)  # all bulk-like: n = 1
    q = sites.mclean_prediction(np.array([-0.2]), 1e-3, kT, n_missing=3)
    assert q.c_C == pytest.approx((sites.mclean_theta([-0.2], 1e-3, kT)[0] + 3e-3) / 4) and q.n_sites == 4


def test_lattice_period():
    assert sites.lattice_period([0, -1, 1], 3.52) == pytest.approx(3.52 / np.sqrt(2))
    assert sites.lattice_period([1, -2, 1], 3.52) == pytest.approx(3.52 * np.sqrt(6) / 2)
    assert sites.lattice_period([2, -4, 2], 3.52) == pytest.approx(3.52 * np.sqrt(6) / 2)
    assert sites.lattice_period([1, 1, 1], 3.52) == pytest.approx(3.52 * np.sqrt(3))


def test_site_table_formats(tmp_path):
    f = tmp_path / "segE0K-by-site.dat"
    f.write_text("# char potential id kind dy_centre(A) dz(A) r_or_dy_from_partial(A) Eseg_relaxed(eV) Eseg_rigid(eV) n_fire\n"
                 "d60 P 10 core 1.0 1.0 1.4 -0.10 -0.09 100\n"
                 "d60 P 11 core-p1 2.0 1.0 0.3 -0.20 -0.19 100\n"
                 "d0 P 12 core 1.0 1.0 1.4 -0.30 -0.29 100\n"
                 "d60 Q 13 core 1.0 1.0 1.4 -0.40 -0.39 100\n")
    t = sites.read_site_table(f, "segE0K-by-site", "P", "d60")
    assert t.ids.tolist() == [10, 11] and t.energy.tolist() == [-0.10, -0.20] and t.kind.tolist() == ["core", "core-p1"]
    g = tmp_path / "s.csv"
    g.write_text("id,E_seg,character\n5,-0.1,d60\n6,-0.2,d0\n")
    assert sites.read_site_table(g, "csv", "P", "d60").ids.tolist() == [5]
    with pytest.raises(ValueError, match="unknown site-table format"):
        sites.read_site_table(g, "xlsx", "P", "d60")
    g.write_text("id,E_seg\n5,-0.1\n5,-0.2\n")
    with pytest.raises(ValueError, match="more than once"):
        sites.read_site_table(g, "csv", "P", "d60")
    with pytest.raises(ValueError, match="no sites"):
        sites.read_site_table(f, "segE0K-by-site", "P", "d90")
    t = sites.SiteTable("x", "csv", "", "P", "d0", np.arange(6), np.array([-.1, -.1, -.1, -.2, -.3, -.3]), np.array([""] * 6))
    assert t.degenerate_sites() == 3


def line_snapshot(a=3.52):
    """Perfect fcc slab, X = [100] periodic (line period a), a straight DXA line along X at the cell centre; Z free."""
    cell = Cell(np.diag([6 * a, 10 * a, 16 * a]), pbc=(True, True, False))
    pos = fcc_in_cell(cell, a) + 1e-3
    y0, z0 = 5 * a + 0.4, 8 * a + 0.3
    line = np.array([[0.0, y0, z0], [cell.line_length, y0, z0]])
    n = len(pos)
    snap = Snapshot(positions=pos.astype(np.float32), types=np.ones(n, np.int8),
                    cell_matrix=np.column_stack([cell.H, cell.origin]), pbc=np.array(cell.pbc),
                    acna=np.full(n, STRUCTURE_CODES["fcc"], np.int8), ptm=np.ones(n, np.int8), seg_points=line,
                    seg_offsets=np.array([0, 2]), seg_true_b=np.zeros((1, 3)), seg_spatial_b=np.zeros((1, 3)),
                    seg_length=np.array([cell.line_length]), seg_family=np.array(["1/6<112>"]),
                    seg_cluster=np.array([1]), seg_infinite=np.array([True]), seg_loop=np.array([False]))
    dy = pos[:, 1] - y0
    r = np.hypot(dy - cell.H[1, 1] * np.round(dy / cell.H[1, 1]), pos[:, 2] - z0)
    return snap, cell, pos.astype(np.float32).astype(float), r


def test_core_sites_on_a_straight_line():
    a = 3.52
    snap, cell, pos, r = line_snapshot(a)
    # the "table": every atom with r < 8 A in ONE line period (0 <= x < a), ids 1..N
    in_period = (pos[:, 0] < a - 1e-6) & (r < 8.0)
    ids = np.nonzero(in_period)[0] + 1
    E = np.where(r[in_period] < 5, -0.1, -0.01)
    t = sites.SiteTable("t", "csv", "", "P", "d0", ids, E, np.array([""] * len(ids)))
    cs = sites.core_sites(t, snap, 5.0, [1, 0, 0], a, surface_exclusion=15.0)
    assert cs.n_periods == 6 and cs.info["period"] == pytest.approx(a) and cs.info["a_0K"] == pytest.approx(a, rel=1e-5)
    # the period comes from the cell, not from the run's a_cell (an MD cell scaled by a few % must not change it)
    assert sites.core_sites(t, snap, 5.0, [1, 0, 0], 1.08 * a).n_periods == 6
    assert cs.r == pytest.approx(r[in_period], abs=0.02)
    assert cs.n_core_tabulated == int((r[in_period] < 5).sum())
    z = pos[:, 2]
    interior = (z - z.min() >= 15) & (z.max() - z >= 15)
    assert cs.n_core_cell == int((interior & (r < 5)).sum()) and cs.coverage == pytest.approx(1.0)
    assert cs.x_extent < a
    t2 = sites.SiteTable("t", "csv", "", "P", "d0", np.array([10 ** 9]), np.array([0.0]), np.array([""]))
    with pytest.raises(ValueError, match="not in the 0 K cell"):
        sites.core_sites(t2, snap, 5.0, [1, 0, 0], a)


def test_mclean_rows_against_measured(tmp_path):
    kT = KB * 300
    E = np.array([-0.15, -0.1, -0.05, 0.0])
    x = 2e-3
    c_pred = sites.mclean_theta(E, x, kT).mean()
    f = write_run(tmp_path, "a", -0.25, 2 * c_pred, 0.9 * x, x)
    runs, _ = across.load_runs([f])
    (key, rs), = across.run_sets(runs).items()
    rows = across.mclean_rows(key, rs, E, 0.0, {})
    r = {(q["quantity"], q["ref"]): q for q in rows}
    assert r[("c_C_mclean", "x_iso")]["value"] == pytest.approx(c_pred)
    assert r[("c_C_measured_over_mclean", "x_iso")]["value"] == pytest.approx(2.0)
    assert r[("c_C_mclean", "x_iso")]["error"] == 0.0                   # x_iso carries no error
    assert all(q["status"] == "definition" for q in rows)                # D2 section 12
    assert r[("n_up_mclean", "c_B")]["value"] <= 1.0


# ------------------------------------------------------------------------------------------- bulk control
def bulk_manifest():
    d = demo()
    for k in ("character", "frame", "frame_source", "z_gp", "z_gp_source", "x_ref_source", "start_type",
              "bulk_control", "bulk_control_reason", "mc_stream"):
        d.pop(k, None)
    d.update(cell_type="bulk_control", run_id="bulk_LeeShim_0.25kT", n_atoms=32000, swap_attempts_per_atom_per_step=0.01)
    return d


def test_bulk_control_manifest_validation():
    assert validate(bulk_manifest()) == []
    d = bulk_manifest()
    d["bulk_control"] = "other"
    assert any("cannot name a bulk_control" in e for e in validate(d))
    d = bulk_manifest()
    d.pop("n_atoms")                                  # D2 section 15: a bulk control records its size
    assert any("needs n_atoms" in e for e in validate(d))
    d = bulk_manifest()
    d.pop("swap_attempts_per_atom_per_step")          # ... and its swap rate
    assert any("swap_attempts_per_atom_per_step is required" in e for e in validate(d))
    d = {**demo(), "bulk_control": "bulk_LeeShim_0.25kT", "x_ref_source": "bulk_control"}
    assert any("swap_attempts_per_atom_per_step is required" in e for e in validate(d))   # a run naming one, too
    assert validate({**d, "swap_attempts_per_atom_per_step": 0.01}) == []
    d = bulk_manifest()
    d["cell_type"] = "precipitate"
    assert any("cell_type" in e for e in validate(d))
    d = demo()
    d.pop("z_gp")                                     # still required for a dislocation cell
    assert any("z_gp" in e for e in validate(d))
    assert validate({**demo(), "cell_type": "dislocation"}) == []


def test_bulk_snapshot_exact_and_find_x_bc(tmp_path):
    a = 3.52
    cell = Cell(np.diag([8 * a, 8 * a, 14 * a]), pbc=(True, True, False))
    pos = fcc_in_cell(cell, a) + 1e-3
    rng = np.random.default_rng(3)
    types = np.where(rng.random(len(pos)) < 0.03, 2, 1).astype(np.int8)
    n = len(pos)
    snap = Snapshot(positions=pos.astype(np.float32), types=types, cell_matrix=np.column_stack([cell.H, cell.origin]),
                    pbc=np.array(cell.pbc), acna=np.ones(n, np.int8), ptm=np.ones(n, np.int8),
                    seg_points=np.empty((0, 3)), seg_offsets=np.array([0]), seg_true_b=np.empty((0, 3)),
                    seg_spatial_b=np.empty((0, 3)), seg_length=np.empty(0), seg_family=np.array([], dtype="U10"),
                    seg_cluster=np.empty(0, np.int32), seg_infinite=np.empty(0, bool), seg_loop=np.empty(0, bool))
    d = load_definitions()
    res = analyse_bulk_snapshot(snap, d, a)
    z = pos[:, 2]
    interior = (z - z.min() >= 15) & (z.max() - z >= 15)
    assert res.value("N", BULK_REGION) == interior.sum()
    assert res.value("c", BULK_REGION) == pytest.approx((types[interior] == 2).mean())
    assert res.value("N", "X") + res.value("N", BULK_REGION) == n
    assert res.value("z1", BULK_REGION) == pytest.approx(12.0) and res.value("z2", BULK_REGION) == pytest.approx(6.0)
    # periodic in Z: no excluded layer
    snap.pbc = np.array([True, True, True])
    snap.cell_matrix = np.column_stack([np.diag([8 * a, 8 * a, 14 * a]), np.zeros(3)])
    assert analyse_bulk_snapshot(snap, d, a).value("N", BULK_REGION) == n
    # x_bc from the bulk control's Stage-3 table
    rows = [s3row("c", BULK_REGION, 0.0301, 0.0004), s3row("N_sol", BULK_REGION, 250.0, 3.0)]
    write_table(tmp_path / "bulk_stage3_D1.csv", rows, {"definition_set": "D1"}, columns=STAGE3_COLUMNS)
    info = find_x_bc("bulk", "D1", [tmp_path / "nothing", tmp_path])
    assert info["value"] == 0.0301 and info["sem"] == 0.0004 and info["enough_solute"]
    assert info["n_atoms"] is None                    # table without size: not accepted (bulk_control_check)
    write_table(tmp_path / "bulk_stage3_D2.csv", rows, {"definition_set": "D2", "stage3": {"n_atoms": 32000,
                                                                                         "last_md_step": 10 ** 6}},
                columns=STAGE3_COLUMNS)
    info2 = find_x_bc("bulk", "D2", [tmp_path])
    assert info2["n_atoms"] == 32000 and info2["last_md_step"] == 10 ** 6
    assert find_x_bc("bulk", "D1", [tmp_path], min_solute=1000)["enough_solute"] is False
    assert find_x_bc("missing", "D1", [tmp_path]) is None


def test_bulk_control_check_size_run_length_and_swap_rate():
    """DEFINITIONS 5 / D2 section 15 (Erik 25 Nov: E1 option a; fix: run length >= 0.95, swap rate within 10 %)."""
    from disegreg.cli import bulk_control_check

    d1 = load_definitions("D1")

    def bc(n, step, swap=1.0):
        return {"n_atoms": n, "last_md_step": step, "swap_attempts_per_atom_per_step": swap}

    same = bulk_control_check(bc(570000, 800000), 570000, 800000, D2, 1.0)
    assert same["accepted"] and not same["smaller"] and "P6" in same["diagnostic_iv"] and same["flags"] == []
    small = bulk_control_check(bc(32000, 800000), 570000, 800000, D2, 1.0)
    assert small["accepted"] and small["smaller"] and small["diagnostic_iv"].startswith("not usable")
    assert small["size_ratio"] == pytest.approx(32000 / 570000) and small["same_run_length"] is True
    assert small["reasons"] == [] and small["flags"] == []
    small_d1 = bulk_control_check(bc(32000, 800000), 570000, 800000, d1, 1.0)
    assert small_d1["accepted"] is False and "not a bulk control under D1" in small_d1["reasons"][0]
    assert bulk_control_check(bc(570000, 800000), 570000, 800000, d1, 1.0)["accepted"]
    unknown = bulk_control_check({"n_atoms": None}, 570000, 800000, D2, 1.0)
    assert unknown["accepted"] is False and unknown["same_run_length"] is None and unknown["diagnostic_iv"] == "not tested"
    # run length: >= 0.95 x the run's last MD step (exact boundary accepted)
    assert bulk_control_check(bc(32000, 760000), 570000, 800000, D2, 1.0)["accepted"]
    short = bulk_control_check(bc(32000, 759999), 570000, 800000, D2, 1.0)
    assert short["accepted"] is False and "too short" in short["reasons"][0] and "NOT accepted" in short["flags"][-1]
    # swap rate within 10 % of the run's
    assert bulk_control_check(bc(32000, 800000, 1.10), 570000, 800000, D2, 1.0)["accepted"]
    assert bulk_control_check(bc(32000, 800000, 0.90), 570000, 800000, D2, 1.0)["accepted"]
    fast = bulk_control_check(bc(32000, 800000, 1.2), 570000, 800000, D2, 1.0)
    assert fast["accepted"] is False and "swap rate differs" in fast["reasons"][0]
    missing = bulk_control_check(bc(32000, 800000, None), 570000, 800000, D2, 1.0)
    assert missing["accepted"] is False and "missing" in missing["reasons"][0]
    assert bulk_control_check(bc(32000, 800000), 570000, 800000, D2, None)["accepted"] is False
    # same dmu and T as the run, any size (DEFINITIONS 5; review finding 4)
    st = {"dmu": 0.88, "temperature": 300}
    ok = bulk_control_check({**bc(570000, 800000), "dmu": 0.88, "temperature": 300}, 570000, 800000, D2, 1.0, st)
    assert ok["accepted"]
    other = bulk_control_check({**bc(570000, 800000), "dmu": 0.85, "temperature": 300}, 570000, 800000, D2, 1.0, st)
    assert other["accepted"] is False and "dmu differs" in other["reasons"][0]
    hot = bulk_control_check({**bc(32000, 800000), "dmu": 0.88, "temperature": 600}, 570000, 800000, D2, 1.0, st)
    assert hot["accepted"] is False and "temperature differs" in hot["reasons"][0]
    assert bulk_control_check(bc(570000, 800000), 570000, 800000, D2, 1.0, st)["accepted"] is False   # dmu unknown
    # a same-size bulk control is not held to the 0.95 / 10 % rules (only flagged)
    same_short = bulk_control_check(bc(570000, 400000, 2.0), 570000, 800000, D2, 1.0)
    assert same_short["accepted"] and "run length differs" in same_short["flags"][0]


def test_proposals_file_holds_only_P5():
    d1 = load_definitions("D1")
    assert d1.status.startswith("rev. 4") and "mclean" not in d1.params
    assert PROPS.label == "P-v0.2" and set(PROPS.params) == {"label", "status", "ordered_domains"}
    assert D2["mclean"]["core_site_radius"] == 5.0 and D2["hysteresis"]["factor"] == 2.0
    from disegreg.definitions import available
    assert available() == ["D1", "D2"]                 # a proposal file is not a definition set


def test_core_sites_merge_periodic_copies():
    """A table that reaches beyond one line period lists some sites twice; the copies count once."""
    a = 3.52
    snap, cell, pos, r = line_snapshot(a)
    in_15 = (pos[:, 0] < 1.5 * a - 1e-6) & (r < 8.0)       # one and a half periods
    ids = np.nonzero(in_15)[0] + 1
    t = sites.SiteTable("t", "csv", "", "P", "d0", ids, np.full(len(ids), -0.1), np.array([""] * len(ids)))
    cs = sites.core_sites(t, snap, 5.0, [1, 0, 0], a)
    one = (pos[:, 0] < a - 1e-6) & (r < 8.0)
    assert cs.duplicate.sum() == in_15.sum() - one.sum() > 0
    assert cs.coverage == pytest.approx(1.0) and cs.info["max_dE_periodic_copies"] == 0.0
    assert cs.x_extent == pytest.approx(a, abs=1e-3) and cs.info["x_extent_distinct"] < a


def test_diagnostic_iv_column(tmp_path):
    """DEFINITIONS 5 (D1: "not tested" for Delta > 0 without bulk control); D2 section 17 (Erik 25 Nov): every run without a
    bulk control, whatever Delta; section 15: a smaller bulk control is "not usable"; same size: test = proposal P6."""
    from disegreg.results import DIAG_IV_NOT_TESTED, DIAG_IV_NOT_USABLE, DIAG_IV_PENDING, diagnostic_iv

    d1 = load_definitions("D1")
    for delta in (-2.0, -0.1, 0.0, 1.0):
        assert diagnostic_iv(D2, delta, None) == DIAG_IV_NOT_TESTED
        assert diagnostic_iv(D2, delta, {"source": "command line --x-ref", "value": 1e-3}) == DIAG_IV_NOT_TESTED
    assert diagnostic_iv(d1, -0.1, None) == "" and diagnostic_iv(d1, 0.0, None) == ""
    assert diagnostic_iv(d1, 1.0, None) == DIAG_IV_NOT_TESTED
    small = {"enough_solute": True, "smaller": True}
    assert diagnostic_iv(D2, -0.1, small) == DIAG_IV_NOT_USABLE
    # D1: a smaller cell is no bulk control -> the D1 rule for runs without one (review finding 3)
    assert diagnostic_iv(d1, 1.0, small) == DIAG_IV_NOT_TESTED and diagnostic_iv(d1, -0.5, small) == ""
    assert diagnostic_iv(D2, -0.1, {"enough_solute": True, "smaller": False}) == DIAG_IV_PENDING
    assert diagnostic_iv(D2, -0.1, {"enough_solute": True, "smaller": None}) == DIAG_IV_NOT_TESTED   # size unknown
    # the results table carries it from the Stage-3 provenance
    f = write_run(tmp_path, "a", -0.5, 0.1, 1e-3, 1e-3)
    rows, _ = collect_runs([f])
    assert "diagnostic_iv" in rows[0]
