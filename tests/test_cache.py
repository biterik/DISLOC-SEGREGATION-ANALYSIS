"""Stage-1 cache without OVITO: file round trip, locations, Burgers families."""

import numpy as np
import pytest

from disegreg.cache import Snapshot, cache_path, default_cache_root
from disegreg.frame import burgers_family


def make_snapshot(n=50, ids=None):
    rng = np.random.default_rng(1)
    pts = [rng.random((4, 3)) * 10, rng.random((6, 3)) * 10]
    return Snapshot(
        positions=rng.random((n, 3)).astype(np.float32), types=rng.integers(1, 3, n).astype(np.int8),
        cell_matrix=np.column_stack([np.diag([10.0, 11.0, 12.0]), [0.0, 0.0, -1.0]]), pbc=np.array([True, True, False]),
        acna=np.ones(n, np.int8), ptm=np.ones(n, np.int8), seg_points=np.vstack(pts),
        seg_offsets=np.array([0, 4, 10]), seg_true_b=np.array([[1 / 6, 1 / 6, -1 / 3], [0.5, 0.5, 0.0]]),
        seg_spatial_b=np.zeros((2, 3)), seg_length=np.array([1.0, 2.0]),
        seg_family=np.array(["1/6<112>", "1/2<110>"], dtype="U10"), seg_cluster=np.array([1, 1], np.int32),
        seg_infinite=np.array([True, False]), seg_loop=np.array([True, False]), ids=ids,
        meta={"source_sha256": "abc", "note": "test"})


@pytest.mark.parametrize("ids", [None, np.arange(10, 60, dtype=np.int64)])
def test_snapshot_roundtrip(tmp_path, ids):
    s = make_snapshot(ids=ids)
    p = s.save(tmp_path / "x" / "snap.npz")
    assert p.is_file() and not list(p.parent.glob("*.tmp*"))
    t = Snapshot.load(p)
    for k in Snapshot._ARRAYS:
        assert np.array_equal(getattr(s, k), getattr(t, k)), k
    assert (t.ids is None) == (ids is None)
    assert t.meta["note"] == "test" and t.meta["cache_format"] == 1
    assert t.n_atoms == 50 and t.n_segments == 2
    assert len(t.polylines()) == 2 and len(t.polylines({"1/6<112>"})) == 1
    assert np.allclose(t.cell().H, np.diag([10.0, 11.0, 12.0])) and t.cell().pbc == (True, True, False)


def test_cache_path_uses_three_levels(tmp_path):
    dump = tmp_path / "PROJ" / "07_THREAD" / "run_d60" / "run_quenched_step100.dump.gz"
    p = cache_path(dump, "abcd1234", tmp_path / "CACHE")
    assert p == tmp_path / "CACHE" / "PROJ" / "07_THREAD" / "run_d60" / "run_quenched_step100.dump.s1-abcd1234.npz"


def test_default_cache_root_env(monkeypatch, tmp_path):
    monkeypatch.setenv("DISEGREG_CACHE", str(tmp_path / "c"))
    assert default_cache_root() == tmp_path / "c"
    monkeypatch.delenv("DISEGREG_CACHE")
    assert default_cache_root().name == "DISEGREG-CACHE"


@pytest.mark.parametrize("b,fam", [((1 / 6, -1 / 3, 1 / 6), "1/6<112>"), ((-0.5, 0, 0.5), "1/2<110>"),
                                   ((0, 1 / 6, -1 / 6), "1/6<110>"), ((1 / 3, 0, 0), "1/3<100>"),
                                   ((1 / 3, 1 / 3, -1 / 3), "1/3<111>"), ((0.2, 0.1, 0), "other")])
def test_burgers_family(b, fam):
    assert burgers_family(b) == fam
