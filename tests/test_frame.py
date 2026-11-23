"""Crystal frame (frame.py): rotation of b and xi, Miller labels, character angle."""

import numpy as np
import pytest

from disegreg import frame
from disegreg.cache import Snapshot

F_D60 = np.array([[0, -1, 1], [2, -1, -1], [1, 1, 1]], float)
F_D60 /= np.linalg.norm(F_D60, axis=1)[:, None]
A = 3.53


def test_labels():
    assert frame.burgers_label(np.array([1, -2, 1]) / 6) == "1/6[1 -2 1]"
    assert frame.burgers_label(np.array([0, -1, 1]) / 2) == "1/2[0 -1 1]"
    assert frame.burgers_label(np.array([1, 1, 1]) / 3) == "1/3[1 1 1]"
    assert frame.burgers_label(np.array([2, -1, -1]) / 6 + 0.01) == "1/6[2 -1 -1]"  # strained cell
    assert frame.burgers_label([0.123, 0.4, 0.0]) == ""
    assert frame.direction_label([0, -0.7, 0.7]) == "[0 -1 1]"
    assert frame.direction_label([1, -2, 1.02]) == "[1 -2 1]"
    assert frame.direction_label([0.13, 0.77, 0.31]) == "[1 5 2]"  # within 0.03 rad
    assert frame.direction_label([0.13, 0.77, 0.31], max_index=3) == ""


def test_rotation_and_character():
    b_cr = np.array([1, -1, 0]) / 2 * A           # perfect b of the d60 cells
    b_sim = F_D60 @ b_cr                            # crystal -> simulation
    assert frame.to_crystal(b_sim, F_D60) == pytest.approx(b_cr)
    assert frame.character_angle(b_sim, [1, 0, 0]) == pytest.approx(60.0)
    assert frame.character_angle(-b_sim, [1, 0, 0]) == pytest.approx(60.0)  # folded into [0, 90]
    assert frame.character_angle([0, 1, 0], [1, 0, 0]) == pytest.approx(90.0)


def test_line_direction_orientation():
    p = np.array([[0, 0, 0], [1, 0.01, 0], [2, -0.01, 0], [3, 0, 0]], float)
    assert frame.line_direction(p) == pytest.approx([1, 0, 0], abs=0.01)
    assert frame.line_direction(p[::-1]) == pytest.approx([-1, 0, 0], abs=0.01)


def test_segment_table_dissociated_d60():
    """Two infinite Shockley partials of the d60 dislocation: 1/6[2 -1 -1] (edge) + 1/6[1 -2 1] (30 deg) = 1/2[1 -1 0]."""
    L = 100.0
    lines = [np.array([[0, 40, 50], [L, 40, 50]], float), np.array([[0, 58, 50], [L, 58, 50]], float)]
    bs = [F_D60 @ (np.array([2, -1, -1]) / 6 * A), F_D60 @ (np.array([1, -2, 1]) / 6 * A)]
    snap = Snapshot(positions=np.zeros((1, 3), np.float32), types=np.ones(1, np.int8), cell_matrix=np.zeros((3, 4)),
                    pbc=np.array([True, True, False]), acna=np.ones(1, np.int8), ptm=np.ones(1, np.int8),
                    seg_points=np.vstack(lines), seg_offsets=np.array([0, 2, 4]), seg_true_b=np.zeros((2, 3)),
                    seg_spatial_b=np.array(bs), seg_length=np.array([L, L]), seg_family=np.array(["1/6<112>"] * 2),
                    seg_cluster=np.array([1, 1]), seg_infinite=np.array([True, True]), seg_loop=np.array([False] * 2))
    t = frame.segment_table(snap, F_D60, A)
    assert [r["b_label"] for r in t] == ["1/6[2 -1 -1]", "1/6[1 -2 1]"]
    assert [r["xi_label"] for r in t] == ["[0 -1 1]"] * 2
    assert [round(r["character_angle"]) for r in t] == [90, 30]
    tot = frame.total_burgers(t, x_axis=F_D60[0])
    assert frame.burgers_label(tot) == "1/2[1 -1 0]"
    assert frame.character_angle(tot, F_D60[0]) == pytest.approx(60.0)
