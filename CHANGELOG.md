# Changelog

All notable changes to `disegreg` are listed here. Versions follow semantic versioning; every release also states
the definition set (D1, D2, ...) it implements (see docs/DEFINITIONS.md).

## [Unreleased]
### Added
- Design (docs/DESIGN.md), definition set D1 draft rev. 2 (docs/DEFINITIONS.md), TURLO-REDO feedback, LICENSE (BSD-3-Clause), CITATION.cff.
- Reproducible conda environment (environment.yml): Python 3.12, OVITO 3.17.2 pinned exactly (from conda.ovito.org), numpy 2.5,
  scipy 1.18; exact solved package list for macOS arm64 in env-lock/.
- Package skeleton (pyproject.toml, `disegreg version` command), definition set D1 as `disegreg/definitions/D1.yaml`
  with a loader that records label and SHA-256.
- `disegreg.geometry`: triclinic periodic cell, polyline densification, distance to DXA lines over all periodic images
  (DEFINITIONS 3.1), projected line coverage f_cov (3.2), fallback core axis from the structure (3.3), neighbour pairs
  with ghost atoms for triclinic cells.
- Definition set D1 draft rev. 3 (Erik, 22 Nov): fallback axis per line defect cluster (every cluster running along >= 50 % of
  the slabs, dipole cells get two axes), used only slab-wise where DXA misses the line (DXA points assigned to clusters by
  the nearest non-fcc atom); `geometry.line_distances` returns r, r_DXA, r_axis, mode and per-cluster coverage.
- Stage 1: `disegreg.io_ovito` (the only OVITO module: dump -> a-CNA, PTM fcc/hcp/bcc/ico RMSD 0.12, DXA, in separate
  pipelines), `disegreg.cache` (one compressed .npz per snapshot, keyed by Stage-1 parameters + OVITO version, validated by
  the dump's SHA-256, atomic writes), `disegreg stage1` command, `disegreg.frame` (Burgers-vector families),
  `disegreg.provenance`. DXA parameters recorded in D1.yaml.
- Verified on Ni-Cu d60: DXA ignores a preceding PTM / a-CNA; the old Ni-Cu "PTM" labels were DXA's internal ones.
- Stage 0: `disegreg.manifest` (schema 1: every field with its meaning; validation of frame orthogonality / handedness / <111>
  normal, Delta = (dmu - dmu_coex)/kT, dmu sign conversion, x_ref source, bulk-control reason), `disegreg manifest template|check`,
  demo manifest examples/manifests/Ni-Cu_d60_LeeShim2004-MEAM_0.25kT-below-coex.yaml.
- Stage 2: `disegreg.regions` (codes X C F S T B O, far_slab40 column, z_gp check against the Shockley lines / axes),
  `disegreg.composition` (counts per region / radial bin / inside R under r and r_DXA, c, E, Gamma(R), region excesses,
  partial separation with merging of split DXA pieces), `disegreg.order` (1st/2nd shells, Warren-Cowley vs x_ref and local c,
  R2, 1st-shell solute clusters), `disegreg.stage2` driver, `disegreg.results` (tidy CSV + provenance JSON),
  `disegreg analyse <manifest>`.
- D1 rev. 3 (23 Nov): Erik confirmed min_slab_fraction and active_below_dxa_coverage 0.5; new [CONFIRM] items: region precedence
  F before O, partial_merge_distance 2 A.
- Stage 3: `disegreg.timeseries` (Flyvbjerg-Petersen block averaging with plateau flag, drift test with the D1 criterion and
  the proposed SE(change) criterion side by side, averaging window, MC-stream reader with restart de-duplication); `analyse`
  writes `<run>_stage3_<D>.csv`; snapshot times from the MD/MC step in the file name.
- Crystal frame: b and xi of every DXA segment in crystal axes with Miller labels, character angle, total b of the line.
- `disegreg table` (tidy results table across runs, refuses mixed definition sets), `disegreg slurm` (Stage-1 job array
  for cmmg), env-lock/disegreg-linux-64.txt (solved on the M5 for glibc 2.38), docs/CMMG.md, docs/HOWTO_PROJECTS.md.
- Test cell: tests/data-generation/ (periodic Ni d60 cell, production construction scaled down; LAMMPS relaxation input and
  job script for cmmg); tests/test_testcell.py runs once the relaxed dump is in tests/data/.
- Regression against the project analyses: regression/run_regression.py (old value / old recipe on disegreg data / D1),
  regression/make_figures.py, results in regression/results/.
- Demo manifests for Nb-Au (d60 Mendelev) and TURLO (t50 R3).
- Synthetic tests against brute-force references (incl. the tilted-cell y-boundary case), exact L1_2 / D0_22 order parameters
  (alpha1 = -1/3, R2 = 4 and 8/3), cluster sizes across the periodic boundary, end-to-end Stage 2 and `disegreg analyse`, OVITO smoke test,
  GitHub Actions: synthetic tests on Python 3.11-3.13, OVITO tests on Linux, conda environment + all tests on macOS arm64.
