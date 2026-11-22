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


def write_lammps_dump(path, cell, ids, types, positions, timestep=1000, pbc=("pp", "pp", "ss")):
    """Write a LAMMPS text dump (triclinic header if the cell is tilted), optionally gzipped (path ending in .gz)."""
    import gzip

    H, o = cell.H, cell.origin
    lx, ly, lz = H[0, 0], H[1, 1], H[2, 2]
    xy, xz, yz = H[0, 1], H[0, 2], H[1, 2]
    xlo, ylo, zlo = o
    tilted = any(abs(t) > 0 for t in (xy, xz, yz))
    lines = ["ITEM: TIMESTEP", str(timestep), "ITEM: NUMBER OF ATOMS", str(len(ids))]
    if tilted:
        lines.append("ITEM: BOX BOUNDS xy xz yz " + " ".join(pbc))
        lines.append(f"{xlo + min(0, xy, xz, xy + xz):.10f} {xlo + lx + max(0, xy, xz, xy + xz):.10f} {xy:.10f}")
        lines.append(f"{ylo + min(0, yz):.10f} {ylo + ly + max(0, yz):.10f} {xz:.10f}")
        lines.append(f"{zlo:.10f} {zlo + lz:.10f} {yz:.10f}")
    else:
        lines.append("ITEM: BOX BOUNDS " + " ".join(pbc))
        lines += [f"{xlo:.10f} {xlo + lx:.10f}", f"{ylo:.10f} {ylo + ly:.10f}", f"{zlo:.10f} {zlo + lz:.10f}"]
    lines.append("ITEM: ATOMS id type x y z c_PE_ATOM")
    body = "\n".join(f"{i} {t} {x:.8f} {y:.8f} {z:.8f} -4.45" for i, t, (x, y, z) in zip(ids, types, positions))
    text = "\n".join(lines) + "\n" + body + "\n"
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt") as f:
        f.write(text)


def screw_dislocation_cell(a=3.52, nx=8, ny=10, nz=7):
    """fcc crystal X = [1-10] (line, periodic), Y = [11-2], Z = [111], with the unrelaxed (isotropic) displacement field
    of a perfect screw dislocation b = a/2[1-10] along X; free surfaces in Y and Z (10 A vacuum).
    Returns (Cell, positions, (y0, z0))."""
    from disegreg.geometry import Cell

    R = np.array([[1, -1, 0], [1, 1, -2], [1, 1, 1]], float)
    R /= np.linalg.norm(R, axis=1)[:, None]
    lx, ly, lz = nx * a / np.sqrt(2), ny * a * np.sqrt(6) / 2, nz * a * np.sqrt(3)
    basis = np.array([[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
    n = int(np.ceil(max(lx, ly, lz) / a * 1.8)) + 2
    g = np.arange(-n, n + 1)
    grid = np.array(np.meshgrid(g, g, g, indexing="ij")).reshape(3, -1).T
    pts = ((grid[:, None, :] + basis[None]) * a).reshape(-1, 3) @ R.T + 1e-6
    keep = ((pts[:, 0] >= 0) & (pts[:, 0] < lx - 1e-4) & (pts[:, 1] >= 0) & (pts[:, 1] < ly)
            & (pts[:, 2] >= 0) & (pts[:, 2] < lz))
    pts = pts[keep]
    y0, z0 = ly / 2 + 0.7, lz / 2 + 0.4
    pts[:, 0] += a / np.sqrt(2) / (2 * np.pi) * np.arctan2(pts[:, 2] - z0, pts[:, 1] - y0)
    cell = Cell(np.diag([lx, ly + 20.0, lz + 20.0]), origin=(0.0, -10.0, -10.0), pbc=(True, False, False))
    return cell, pts, (y0, z0)
