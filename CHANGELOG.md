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
- Synthetic tests against brute-force references (incl. the tilted-cell y-boundary case), OVITO smoke test,
  GitHub Actions: synthetic tests on Python 3.11-3.13, OVITO tests on Linux, conda environment + all tests on macOS arm64.
