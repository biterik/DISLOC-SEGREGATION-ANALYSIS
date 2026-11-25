# Changelog

All notable changes to `disegreg` are listed here. Versions follow semantic versioning; every release also states
the definition set (D1, D2, ...) it implements (see docs/DEFINITIONS.md).

## [Unreleased]
v0.2 development on branch v0.2-dev (docs/PLAN_v0.2.md). **New definition set D2** (25 Nov, Erik): D1 rev. 4 unchanged +
sections 12-17 of docs/DEFINITIONS.md (from proposals P1-P4 and questions E1, N1); D2 is the default of this branch, D1 stays
selectable (`--definitions D1`) and gives the v0.1 outputs. The ordered-domain sizes (P5) stay a proposal
(disegreg/definitions/proposals_v0.2.yaml), marked `status = proposal` in every table. Open decisions:
docs/OPEN-DECISIONS_v0.2.md. (v0.1.1 was released from main on 25 Nov, see [0.1.1] below; main merged into v0.2-dev the same day.)
### Added (D2, 25 Nov)
- `disegreg/definitions/D2.yaml`: every D1 parameter value for value (`extends: D1`) + sections `mclean` (12: core sites
  r < 5 A of the 0 K DXA line, tabulated sites only, coverage < 90 % flagged), `hysteresis` (13: different if
  |d| > 2 sqrt(SEM_up^2 + SEM_down^2), Gamma(20 A) and c_C), `bulk_control` (14: Gamma_error_total incl. SEM(x_bc);
  15: smaller bulk-control cells accepted for x_bc, not for diagnostic (iv)), `alongline` (16: 5 A bins, core atoms,
  periodic if p < 0.01). Default definition set = D2 (`DEFAULT_LABEL`); `disegreg across` takes the tables' own set.
- Bulk controls (D2 section 15): manifest field `n_atoms` (required for `cell_type: bulk_control`; schema stays 1);
  `disegreg analyse` checks the snapshots against it and records size and last MD step; a run taking x_bc gets
  `n_atoms_run`, `size_ratio`, `smaller`, `accepted`, `diagnostic_iv`, `same_run_length` and flags in its provenance.
  Under D1 a smaller cell is not accepted (x_bc reported as a diagnostic only, paper_flag kept).
- Tests: D2 = D1 + new sections (parameters), every D1 Stage-2 record bit-identical under D2 on a synthetic snapshot,
  Stage 4 with D1 tables = uptake + hysteresis differences only, bulk-control size / run-length check (150 tests).
- regression/check_d2_vs_d1.py: the 59 demo runs give identical tables under D1 (+ proposals) and D2 (2 318 883 cells; Stage 4
  identical except `status`); regression/results/v0.2-demo/ regenerated from D2.
- D2 rev. 2 (Erik 25 Nov, E1 fix): a smaller bulk control is accepted only if its last MD step >= 0.95 x the run's and its
  swap attempts per atom per MD step are within 10 % of the run's (new manifest field `swap_attempts_per_atom_per_step`,
  required for bulk controls and for runs naming one; D2.yaml `smaller_cell_min_run_length_fraction`,
  `smaller_cell_swap_rate_tolerance`); otherwise not accepted, x_iso (and c_B) used, reasons in `x_bc.reasons` and
  `reference_rule`.
- D2 rev. 2 (Erik 25 Nov, N1): section 17 / YAML section `diagnostic_iv`; results-table column `diagnostic_iv` (last
  column): "not tested" for every run without a bulk control (any Delta; with D1 the D1 rule: Delta > 0 only), "not usable"
  with a smaller bulk control, "test not yet defined" with a same-size one. The test itself: [PROPOSAL v0.2] P6 in
  docs/DEFINITIONS.md (bulk-control clusters vs a random alloy of the same composition + drift of x_bc), text only.
### Changed (D2, 25 Nov)
- Status: along-line quantities, McLean, the hysteresis verdict and Gamma_error_total are `definition` in D2 tables
  (along-line quantities in D1-labelled tables, e.g. the 24 Nov demo, stay `proposal`); only `domain_*` (P5) is `proposal`.
- `stage4()` takes the definition set instead of the proposals; `load_site_cores()` reads `definitions["mclean"]`
  (new parameter `min_coverage` 0.90, the value that was hard-coded); Stage-4 rows have an empty `proposal_set`.
- `disegreg analyse`: the along-line spectrum comes from the definition set (D2), no longer from the proposals file;
  `--no-proposals` now only leaves out P5.
- proposals_v0.2.yaml holds only `ordered_domains` (P5).
### Added (v0.2)
- Stage 4, `disegreg across` (`disegreg.across`): run sets ordered by Delta (up = bottom-up + seeded, down = top-down), uptake
  exponent n_up = Delta ln c_C / Delta ln x_ref between neighbouring offsets for every shared reference with propagated errors,
  hysteresis bottom-up vs top-down at the same Delta (differences; verdict = proposal P2), McLean prediction (proposal P1),
  error of Gamma vs x_bc incl. SEM(x_bc) (proposal P3). Output: tidy table with `status` per row.
- `disegreg.sites`: 0 K site-energy tables (adapters `segE0K-by-site` for Ni-Cu thread 09 / Nb-Au thread 06, and `csv`),
  core sites = r < 5 A of the 0 K DXA line, periodic copies merged, coverage and degenerate sites reported; McLean occupancy
  and its local exponent. Site configuration YAML (examples/stage4/sites.yaml).
- Bulk-control runs: manifest field `cell_type: bulk_control` (optional; schema 1 unchanged), `stage2.analyse_bulk_snapshot`,
  `disegreg analyse` analyses bulk controls first and feeds x_bc automatically into the runs naming them (`--bulk-dir`), with
  the D1 minimum of 100 solute atoms.
- Results table: last column `status` (definition | proposal); a proposal with role headline is refused.
- Stage 1: optional periodicity override (`pbc`) for LAMMPS data files; v0.1 cache keys unchanged.
- `disegreg plot gamma | timeseries | delta | uptake | regions` (`disegreg.plotting`, item B): house style (fonts >= 13 pt,
  Okabe-Ito palette checked for colour-vision deficiency, fixed colour per character, symmetric diverging scale for signed
  maps); optional dependency `plot` (matplotlib).
- Along-line core composition (`disegreg.alongline`, proposal P4, item C): core solute counts in ~5 A X bins, dispersion
  index, spectrum normalised to binomial noise, peak wavelength and p-value (exact for the q = n/2 mode); Stage-2 quantities
  `alongline_*` (status proposal); `disegreg plot alongline` (kymograph with symmetric scale + window-averaged spectrum).
- Ordered-domain sizes (`disegreg.domains`, proposal P5, item D): ordered solute (no solute 1st, >= 4 of 6 solute 2nd
  neighbours), 2nd-shell-connected domains with their surface atoms, equivalent diameter, per region and in bulk controls;
  Stage-2 quantities `domain_*` (status proposal). Synthetic L1_2 Ni3Al spheres 1.6-4 nm recovered to 1-5 %.
- `disegreg analyse` adds the proposal quantities by default (`--no-proposals`: D1 only).
- Stage 3 warns and counts snapshots later than the end of the window signal (MC stream), which v0.1.0 drops silently.
- examples/stage4/: demo manifests for 52 Nb-Au and Ni-Cu thread-07 runs (make_demo_manifests.py), site configuration.
- regression/v02_demo_summary.py and regression/results/v0.2-demo/: Stage 4 + summaries + figures of the 59 demo runs.
### Fixed (v0.2 branch; no v0.1.0 number affected)
- Uptake exponent: an unknown SEM of c_C (single snapshot) gives an unknown error (was treated as 0).
- McLean sites: line period from the 0 K cell's own lattice parameter (the run's a_cell gave 39 instead of 40 periods for
  Mendelev); periodic copies merged; the "more than one period" flag looks at the distinct sites.
- `plotting.save` keeps names with dots (run ids such as ..._0.10kT...); `with_suffix` cut them and files overwrote each other.
- `results.write_table` writes numpy numbers as plain numbers (numpy >= 2 repr was "np.float64(...)"); also in 0.1.1.
- `paper_flag` is dropped only when x_bc is actually used, not as soon as a bulk control is named; also in 0.1.1.
### Added (main, after v0.1.1)
- Zenodo DOI of v0.1.1 (10.5281/zenodo.24535067) in README.md and CITATION.cff.

## [0.1.1] - 2026-11-25
Definition set **D1, rev. 4** (unchanged). Bug-fix release: no number computed by v0.1.0 changes; only the column
`tool_version` and the provenance differ.
### Fixed
- Stage 3: snapshots later than the end of the window signal (the MC stream) lie outside every averaging window and were
  dropped without a message. They are still not used (D1 unchanged), but are now counted (`snapshots_after_window_signal`,
  `snapshots_total` in the Stage-3 provenance) and `disegreg analyse` prints a WARNING. Seen on Mac copies of MC streams
  that were rsynced before the dumps (Nb-Au thread 07: about half of the snapshots of every run unused).
- `results.write_table` writes numpy scalars as plain numbers (numpy >= 2 would have written "np.float64(...)"); v0.1.0
  always passed Python floats, so no v0.1.0 table differs.
- `paper_flag` ("no paper number before bulk control", ACE Ni-Cu) is lifted only when x_bc is actually given
  (`--x-ref` with `x_ref_source: bulk_control`), not already when the manifest names a bulk control. No v0.1.0 manifest
  names one, so no v0.1.0 output differs.
### Changed
- docs/HOWTO_PROJECTS.md: rsync the MC stream together with the dumps.
### Added
- Zenodo DOIs in README.md and CITATION.cff: concept 10.5281/zenodo.24513776, v0.1.0 10.5281/zenodo.24513777.
### Changed
- CITATION.cff abstract and README: state that v0.1 covers one run at a time plus the results table across runs;
  cross-run diagnostics (uptake exponent, McLean, bulk-control and hysteresis comparisons) are planned for v0.2.

## [0.1.0] - 2026-11-24
First release. Definition set **D1, rev. 4 (final)**: every number quoted with this version carries the label D1.
Stages 0-3: run manifest (Stage 0), OVITO a-CNA / PTM / DXA with a cache (Stage 1), distances to the lines with exact
triclinic periodic images and the structural fallback axis, regions, concentrations, line excess Gamma(R) with
Gamma(20 A) as headline, Warren-Cowley order, solute clusters (Stage 2), block-averaged time averages with the
SE(change) drift test (Stage 3), the tidy results table across runs, and full provenance. Reference compositions:
until a run has its bulk control, against both the in-cell far field c_B and the isotherm x_iso. Tested on macOS arm64
(M5) and Linux x86_64 (cmmg, GitHub Actions); regression against the Ni-Cu, Nb-Au and TURLO-REDO project analyses
in regression/. Requires OVITO 3.17.2 exactly (environment.yml). Planned for v0.2: Stage 4 (uptake exponent, McLean,
hysteresis), along-line Fourier analysis, ordered-domain sizes, plotting.

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
- D1 rev. 4 (Erik, 24 Nov; D1.yaml status "rev. 4 (2026-11-24)"): no open definition point left.
  - Drift test: `change_se` (|change| < 2 SE(change)) is the default (`timeseries.drift_criterion`, `drift_factor`); the
    rev. 3 rule is criterion `sem` and its verdict is stored in the extra column `flat_rev3_sem`. `--drift-criterion` takes
    `change_se | sem` (the old name `d1` is gone).
  - Reference composition: new column `ref` in every table (Stage 2, Stage 3, results). `analyse` reports E, Gamma,
    alpha vs x_ref and `c_B_over_x_ref` against the snapshot's own far field `c_B` and the manifest's `x_iso` (and against
    `--x-ref`, named `x_bc` for a bulk-control manifest). Rule `reference.paper_hold_without_bulk_control` (ACE Ni-Cu):
    column `paper_flag` = "no paper number before bulk control" until the run has its bulk control.
  - Gamma(20 A) is the headline: column `role` = headline / check (Gamma(40 A)) / curve for Gamma under r; `check_required`
    is true when the fallback axis was active in the window (per-snapshot quantity `fallback_active`). New per-run file
    `<run>_gammaR_<D>.csv` with the Gamma(R) curve, R = 5-60 A in 5 A steps (plus 80, 100), under r and r_DXA, per reference.
- tests/test_testcell.py expects the DXA Burgers vector 1/2[-1 1 0] (the builder's nominal b with opposite sign, DEFINITIONS 1).
- The relaxed test cell tests/data/testcell_relaxed.dump.gz (0.7 MB, relaxed on cmmg 24 Nov) is in the repository, so the
  test-cell tests run in GitHub Actions too.
- Regression with tolerances (regression/check_tolerance.py): positions 0.5 A, concentrations max(2 %, one solute atom),
  order parameters max(2 %, 0.01). Ni-Cu: 16 runs exact, the 3 ACE runs within tolerance; Nb-Au: 24 within tolerance.
  DXA is deterministic on one machine (four Stage-1 runs of the ACE d90 snapshot on the M5 bit-identical).
