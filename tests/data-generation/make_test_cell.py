#!/usr/bin/env python3
"""Small periodic d60 dislocation cell in Ni for the OVITO tests (docs/DESIGN.md, Testing). numpy only.

Same construction as the production cells (Cu-AT-DISLOCS-Ni-EAM-MEAM-ACE thread 00, gate reports), scaled down:
  frame       X = [0 -1 1] (line, periodic), Y = [2 -1 -1] (glide direction, periodic), Z = [1 1 1] (free)
  dislocation b = a/2 [1 -1 0] (60 deg), isotropic Volterra field (nu = 0.3) at (y0, z0), z0 midway between two (111)
              planes, cut along the glide plane on the LEFT (y < y0): there the half-crystals are displaced by +-b/2
  box         period vector along Y = (0, L0, 0) - b/2, i.e. LAMMPS ly = L0 - b_e/2 and tilt xy = -b_s/2, so the top
              half-crystal continues perfectly across the Y boundary and the bottom one with a lattice-vector shift;
              bottom atoms duplicated by the shorter period are removed
  solute      2 % of the atoms type 2 (Cu), random, seed fixed

Writes a LAMMPS data file (atom_style atomic, boundary p p s). The relaxation (in.relax_testcell, Fischer 2018 EAM)
runs on cmmg (docs/CMMG.md, section "Test cell"); the relaxed dump goes to tests/data/.

    python make_test_cell.py OUT.lmp [--nx 8 --ny 35 --nz 64 --a 3.52]
"""

import argparse

import numpy as np

FRAME = np.array([[0, -1, 1], [2, -1, -1], [1, 1, 1]], float)
B_CRYSTAL = np.array([1, -1, 0]) / 2.0  # units of a


def build(nx=8, ny=35, nz=64, a=3.52, nu=0.3, x_cu=0.02, seed=7):
    F = FRAME / np.linalg.norm(FRAME, axis=1)[:, None]
    lx = nx * a / np.sqrt(2)              # repeat along [0 -1 1]: a/sqrt(2)
    l0 = ny * a * np.sqrt(6) / 2          # repeat along [2 -1 -1]: a sqrt(6)/2 (1/2[2 -1 -1] is a lattice vector)
    d111 = a / np.sqrt(3)
    lz = nz * d111
    # perfect fcc, rotated into the frame, inside [0, lx) x [0, l0) x [0, lz)
    basis = np.array([[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
    n = int(np.ceil(max(lx, l0, lz) / a * 1.8)) + 2
    g = np.arange(-n, n + 1)
    grid = np.array(np.meshgrid(g, g, g, indexing="ij")).reshape(3, -1).T
    pts = ((grid[:, None, :] + basis[None]) * a).reshape(-1, 3) @ F.T
    eps = 1e-6
    keep = ((pts[:, 0] >= -eps) & (pts[:, 0] < lx - eps) & (pts[:, 1] >= -eps) & (pts[:, 1] < l0 - eps)
            & (pts[:, 2] >= -eps) & (pts[:, 2] < lz - d111 / 2))
    pts = pts[keep]
    pts[np.abs(pts) < eps] = 0.0

    b = F @ (B_CRYSTAL * a)               # simulation frame: (b_s, b_e, 0)
    b_s, b_e = b[0], b[1]
    y0 = 0.5 * l0 + 0.25 * a * np.sqrt(6) / 2
    z0 = (round(nz / 2) - 0.5) * d111     # midway between two (111) planes
    Y, Z = pts[:, 1] - y0, pts[:, 2] - z0
    r2 = Y ** 2 + Z ** 2
    th = np.arctan2(Z, Y)                 # cut at theta = +-pi: glide plane, y < y0
    u = np.zeros_like(pts)
    u[:, 0] = b_s / (2 * np.pi) * th
    u[:, 1] = b_e / (2 * np.pi) * (th + Y * Z / (2 * (1 - nu) * r2))
    u[:, 2] = -b_e / (2 * np.pi) * ((1 - 2 * nu) / (4 * (1 - nu)) * np.log(r2 / a ** 2) + (Y ** 2 - Z ** 2) / (4 * (1 - nu) * r2))
    pos = pts + u

    # box: a = (lx, 0, 0), b = (xy, ly, 0) = (0, l0, 0) - b/2, c = (0, 0, lz_box); wrap, then drop duplicates
    xy, ly = -b_s / 2, l0 - b_e / 2
    H = np.array([[lx, xy, 0], [0, ly, 0], [0, 0, 1.0]])
    s = np.linalg.solve(H, pos.T).T
    s[:, :2] -= np.floor(s[:, :2])
    pos = (H @ s.T).T
    from scipy.spatial import cKDTree

    shifts = [(i, j) for i in (-1, 0, 1) for j in (-1, 0, 1)]
    imgs = np.vstack([pos + i * H[:, 0] + j * H[:, 1] for i, j in shifts])
    owner = np.tile(np.arange(len(pos)), len(shifts))
    pairs = cKDTree(imgs).query_pairs(1.2, output_type="ndarray")
    pairs = np.unique(np.sort(np.column_stack([owner[pairs[:, 0]], owner[pairs[:, 1]]]), axis=1), axis=0)
    pairs = pairs[pairs[:, 0] != pairs[:, 1]]
    drop = np.unique(pairs[:, 1])
    pos = np.delete(pos, drop, axis=0)
    rng = np.random.default_rng(seed)
    types = np.where(rng.random(len(pos)) < x_cu, 2, 1)
    info = dict(lx=lx, ly=ly, xy=xy, l0=l0, b_s=b_s, b_e=b_e, y0=y0, z0=z0, n_atoms=len(pos), n_dropped=len(drop),
                zmin=pos[:, 2].min(), zmax=pos[:, 2].max())
    return pos, types, info


def write_data(path, pos, types, info):
    zlo, zhi = info["zmin"] - 0.5, info["zmax"] + 0.5
    with open(path, "w") as f:
        f.write(f"# disegreg test cell: Ni d60, b = a/2[1 -1 0], X=[0 -1 1] Y=[2 -1 -1] Z=[1 1 1]; "
                f"y0 = {info['y0']:.6f} z0 = {info['z0']:.6f} (unrelaxed, isotropic Volterra field)\n\n")
        f.write(f"{len(pos)} atoms\n2 atom types\n\n")
        f.write(f"0.0 {info['lx']:.10f} xlo xhi\n0.0 {info['ly']:.10f} ylo yhi\n{zlo:.10f} {zhi:.10f} zlo zhi\n")
        f.write(f"{info['xy']:.10f} 0.0 0.0 xy xz yz\n\nMasses\n\n1 58.6934\n2 63.546\n\nAtoms # atomic\n\n")
        for k, (p, t) in enumerate(zip(pos, types), start=1):
            f.write(f"{k} {t} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out")
    ap.add_argument("--nx", type=int, default=8)
    ap.add_argument("--ny", type=int, default=35)
    ap.add_argument("--nz", type=int, default=64)
    ap.add_argument("--a", type=float, default=3.52)
    args = ap.parse_args()
    pos, types, info = build(args.nx, args.ny, args.nz, args.a)
    write_data(args.out, pos, types, info)
    print({k: (round(v, 6) if isinstance(v, float) else v) for k, v in info.items()})
