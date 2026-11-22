"""Independent reference implementations for the tests (brute force, exact, slow)."""

import itertools

import numpy as np


def point_segment_distance(p, a, b):
    """Exact distance from points p (N,3) to segments a->b (M,3); returns (N, M)."""
    ab = b - a
    ap = p[:, None, :] - a[None, :, :]
    denom = np.einsum("mk,mk->m", ab, ab)
    t = np.clip(np.einsum("nmk,mk->nm", ap, ab) / np.where(denom > 0, denom, 1.0), 0.0, 1.0)
    closest = a[None] + t[..., None] * ab[None]
    return np.linalg.norm(p[:, None, :] - closest, axis=2)


def brute_force_line_distance(cell, positions, polylines, n_images=3):
    """Exact distance to continuous polylines, minimised over all periodic images with |n| <= n_images."""
    pos = np.asarray(positions, float)
    best = np.full(len(pos), np.inf)
    pax = cell.periodic_axes
    for pl in polylines:
        pl = np.asarray(pl, float)
        a, b = pl[:-1], pl[1:]
        for combo in itertools.product(range(-n_images, n_images + 1), repeat=len(pax)):
            v = np.zeros(3)
            v[pax] = combo
            shift = cell.H @ v
            best = np.minimum(best, point_segment_distance(pos, a + shift, b + shift).min(axis=1))
    return best


def brute_force_pairs(cell, positions, cutoff, n_images=2):
    """Set of (i, j), i < j, with minimum-image distance < cutoff (brute force over images)."""
    pos = np.asarray(positions, float)
    pax = cell.periodic_axes
    n = len(pos)
    dmin = np.full((n, n), np.inf)
    for combo in itertools.product(range(-n_images, n_images + 1), repeat=len(pax)):
        v = np.zeros(3)
        v[pax] = combo
        d = np.linalg.norm(pos[:, None, :] - (pos[None, :, :] + cell.H @ v), axis=2)
        dmin = np.minimum(dmin, d)
    i, j = np.nonzero(np.triu(dmin < cutoff, k=1))
    return set(zip(i.tolist(), j.tolist())), dmin


def fcc_in_cell(cell, a):
    """All fcc lattice sites (cubic orientation, lattice parameter a) inside the cell; cell vectors must be lattice vectors."""
    basis = np.array([[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]) * a
    lo = np.floor(np.minimum(0, cell.H.sum(axis=1).min()) / a) - 2
    span = np.abs(cell.H).sum(axis=1).max()
    hi = np.ceil(span / a) + 2
    g = np.arange(lo, hi + 1)
    grid = np.array(np.meshgrid(g, g, g, indexing="ij")).reshape(3, -1).T * a
    pts = (grid[:, None, :] + basis[None]).reshape(-1, 3) + cell.origin
    s = cell.to_fractional(pts)
    eps = 1e-9
    inside = np.all((s >= -eps) & (s < 1 - eps), axis=1)
    return pts[inside]
