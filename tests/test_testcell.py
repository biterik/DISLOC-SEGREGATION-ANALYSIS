"""The relaxed dissociated test cell (tests/data-generation/, relaxed on cmmg): Stage 1 + 2 end to end on a real
dissociated dislocation. Skipped until tests/data/testcell_relaxed.dump.gz exists."""

from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("ovito")
pytestmark = pytest.mark.ovito

DUMP = Path(__file__).parent / "data" / "testcell_relaxed.dump.gz"
FRAME = np.array([[0, -1, 1], [2, -1, -1], [1, 1, 1]], float)
FRAME /= np.linalg.norm(FRAME, axis=1)[:, None]
A = 3.52
Z0 = 64.016598  # glide plane of make_test_cell.py (default sizes)


@pytest.fixture(scope="module")
def result():
    if not DUMP.is_file():
        pytest.skip("relaxed test cell not yet produced (docs/CMMG.md, 'Test cell')")
    from disegreg import io_ovito
    from disegreg.definitions import load_definitions
    from disegreg.stage2 import analyse_snapshot

    d = load_definitions()
    snap = io_ovito.compute_stage1(DUMP, d)
    return snap, analyse_snapshot(snap, d, A, z_gp=Z0, x_ref=0.02, frame_matrix=FRAME, keep_arrays=True)


def test_dissociated_into_two_shockley_partials(result):
    snap, res = result
    assert res.info["mode"] == "polyline" and res.info["f_shockley"] == pytest.approx(1.0)
    assert res.info["segments"] == {"1/6<112>": 2}
    angles = sorted(round(r["character_angle"]) for r in res.info["segment_table"])
    assert abs(angles[0] - 30) <= 3 and abs(angles[1] - 90) <= 3
    # DXA b with xi along +X is the NEGATIVE of the builder's nominal b (make_test_cell.py: 1/2[1 -1 0]); the production
    # d60 cells show the same relation (nominal 1/2[-1 1 0] -> DXA 1/2[1 -1 0]). DEFINITIONS 1; not an analysis error.
    assert res.info["b_total_label"] == "1/2[-1 1 0]"
    assert res.info["b_total_character_angle"] == pytest.approx(60.0, abs=1.0)


def test_partials_glide_plane_and_regions(result):
    snap, res = result
    assert not res.info["z_gp_flag"]
    sep = res.value("partial_separation_mean")
    assert 5.0 < sep < 40.0 and res.value("partial_separation_valid_bins") == res.info["partial_separation_bins"]
    assert res.value("N", region="F") > 0  # the stacking-fault ribbon between the partials
    n = sum(res.value("N", region=c) for c in "XCFSTBO")
    assert n == snap.n_atoms
