# disegreg -- package design (state 22 Nov 2026, agreed with Erik in the claude.ai design session)

Repository: ~/DEVEL/DISLOC-SEGREGATION-ANALYSIS -> GitHub biterik/DISLOC-SEGREGATION-ANALYSIS (PUBLIC), licence BSD-3-Clause
(as Erik's other Python analysis tools hull0k, ovito-auto-viz, openbis-mcp-server). Python package `disegreg`, command-line tool `disegreg`
("dislocseg" was dropped: it reads as "dislocation segment"). Definitions: docs/DEFINITIONS.md (set D1, draft rev. 2, under Erik's review).

## Purpose
ONE analysis code for all simulation projects of the paper "Dislocation Defect Phases and how to find them" (Acta Materialia), so every number
uses the same definitions: Ni-Cu, Ni-Au/Ni-Nb, Ni-Al (new), TURLO-REDO Cu-Zr. Background and the paper's diagnostics (i)-(vi):
~/Desktop/PUBLISH/2026_Dislocation-Defect-Phases/CLAUDE-NOTES/PAPER-PLAN_2026-11-22.md and
~/Desktop/PUBLISH/2026_Dislocation-Defect-Phases/CLAUDE-NOTES/INSTRUCTIONS_ANALYSIS-TOOL_2026-11-22.md.

## Key principle
Only ONE module imports OVITO (io_ovito.py). Everything else is numpy/scipy, testable without OVITO.

## Layout
```
disegreg/
  manifest.py     Stage 0: run manifest (YAML), validated on load; `disegreg manifest draft <run dir>` pre-fills it from file names + project files
  io_ovito.py     the ONLY OVITO module: read dump (.gz, triclinic), a-CNA labels, PTM labels (fcc/hcp/bcc/ico, RMSD 0.12), DXA lines (b, xi, FS/RH)
  cache.py        Stage 1 cache read/write (compressed .npz per snapshot) + provenance
  frame.py        crystal frame from cell-builder metadata; b and xi into crystal axes (never typed in); Miller labels
  geometry.py     triclinic periodic images (full cell vectors incl. xy tilt), polylines densified to 0.25 A, non-periodic KD-tree;
                  ghost atoms for neighbour search (scipy's periodic KD-tree only handles rectangular boxes); line coverage; fallback core axis
  regions.py      exclusive region codes from definitions/D1.yaml (versioned separately from the code)
  composition.py  concentrations by region, radial bins, Gamma(R) (under r_DXA and the mixed distance), partial separation
  order.py        solute clusters, Warren-Cowley alpha1/alpha2 (vs x_ref and vs local c), 2nd-neighbour ratio R2
  alongline.py    (v0.2) core composition along the line + Fourier spectrum
  timeseries.py   Stage 3: drift test, block averaging (Flyvbjerg-Petersen), averaging window
  across.py       Stage 4: x_ref hierarchy, uptake exponent, McLean prediction, bottom-up vs top-down
  results.py      tidy long table (run, quantity, region, value, error, definition set, tool version) + provenance
  plotting.py     Erik's house style: fonts >= 13 pt at slide size, colour-blind-safe, SYMMETRIC colour scales for signed maps
  cli.py          disegreg manifest | stage1 | slurm | analyse | table | plot | version
  definitions/D1.yaml all numerical definition parameters (R_far = 60 A, core 5 A, shell 10 A, fault 4 A, surface 15 A, ...);
                  inside the package so a pip-installed release carries its definitions (versioned separately from the code)
docs/DEFINITIONS.md   (the paper's Methods section quotes it)   docs/DESIGN.md (this file)
tests/  regression/  examples/manifests/  CITATION.cff  LICENSE  CHANGELOG.md  pyproject.toml  .github/workflows/ci.yml
```

## Stages
- Stage 0 (manifest, per run, written by the PROJECT instance in its own folder): project, run id, potential label, host/solute element
  (type 1 host, type 2 solute), character d0/d30/d60/d90, crystal frame (X = line, Y, Z = (111) normal) from the cell-builder metadata,
  a_cell(T), T, dmu in the convention mu_solute - mu_host (with the conversion from the LAMMPS input sign), dmu_coex(T) and which line
  (incoherent / coherent), Delta = (dmu - dmu_coex)/kT, dump glob, timestep, z boundary layers, glide-plane height z_gp, bulk-control run id
  or the reason there is none, x_ref source, start type (bottom-up / top-down / seeded), MC-stream file.
- Stage 1 (expensive, once per snapshot, cached): positions (float32), types, a-CNA labels, PTM labels, DXA lines (points, segment types,
  b and xi in simulation AND crystal frame, lengths), cell, SHA-256 of the dump, versions. ~6-7 MB per snapshot (compressed .npz).
  Distances and regions are deliberately NOT cached: they are cheap (KD-tree, ~1 s) and definitions may still change, so a definition change
  never forces OVITO to run again. Cost ~1 core-minute per snapshot (~250 core-h for 100 runs x 150 snapshots).
- Stage 2 (cheap, per snapshot, from the cache): r_DXA, line coverage, fallback axis, r_axis, regions, compositions, Gamma(R), order, clusters.
- Stage 3 (per run): time series, drift test, averaging window, block-averaged means and errors.
- Stage 4 (across runs): x_ref, uptake exponent, McLean, hysteresis, the paper's results table.

## Where it runs (Erik's rule: computations on the M5 or on cmmg, never in Claude's cloud workspace)
- Stage 1 production: cmmg (x86_64 Linux, `pip install ovito` works), Slurm job arrays. The FULL time series exist only on cmmg.
- Development and regression: natively on the M5 macOS (OVITO wheels for macOS arm64). The Mac holds ~5 selected quenched snapshots per run
  (Ni-Cu: 366 dumps in 71 run dirs; Nb-Au: 359 in 51; TURLO: ~20 per run for t50/t51 after Erik's rsync of 22 Nov).
- NOT in the M5 Linux VM of the claude.ai desktop app (aarch64 Linux: no OVITO wheel) and NOT in Claude's cloud workspace.
- Cache: cmmg /u/$USER/PTMP/DISEGREG-CACHE; M5 ~/Desktop/SIMULATIONS/DISEGREG-CACHE. Created automatically on first use;
  overridable by environment variable DISEGREG_CACHE or --cache. If PTMP is purged the cache is simply rebuilt from the dumps.
- Pin the OVITO version exactly (DXA output can change between versions) and record it in every output.

## Data facts (checked 22 Nov)
- Ni cells: ~587 000-589 000 atoms, ~104-106 x 249-252 x 249-251 A, LAMMPS "ITEM: BOX BOUNDS xy xz yz pp pp ss" (triclinic, small xy tilt
  ~1.1 A in d0/d30/d60; d90 orthogonal), columns "id type x y z c_PE_ATOM", plain text ~25 MB per quenched dump, fixed z boundary layers.
- Frames (builder metadata, ~/Desktop/SIMULATIONS/Cu-AT-DISLOCS-Ni-EAM-MEAM-ACE/00_CORRECTED-CELLS-A0-3.52/):
  d0, d60: X = [0 -1 1], Y = [2 -1 -1], Z = [1 1 1];  d30, d90: X = [1 -2 1], Y = [1 0 -1], Z = [1 1 1].
- TURLO d90 Cu-Zr: X = [1 -2 1], Y = [1 0 -1], Z = [1 1 1], orthogonal box "pp pp ss", FREE (111) z faces, type 1 Cu, type 2 Zr, 800 K,
  dmu_coex(800 K) = -3.6977 eV (mu_Zr - mu_Cu), kT = 0.068936 eV.
- In every project a larger dmu (mu_solute - mu_host) means more solute (checked: Nb-Au Mendelev and EAM-X, TURLO).

## Starting material and regression references (READ ONLY -- never modify anything in the simulation folders)
- Nb-Au: ~/Desktop/SIMULATIONS/Nb-Au-AT-DISLOCS-Ni-ACE-EAM/07_MDMC-SOLUTE-SEGREGATION-T300-WHOLE-INTERIOR/analysis/harvest_2026-11-21/
  (region_order.py, analyse_07_harvest.py, stream_drift_07.py, to_npz.py, time-avg_07_2026-11-21.dat) and
  .../07_MDMC-SOLUTE-SEGREGATION-T300-WHOLE-INTERIOR/analysis/time-avg_EAMX-round1_2026-11-03/dxa_snapshot.py
- Ni-Cu: ~/Desktop/SIMULATIONS/Cu-AT-DISLOCS-Ni-EAM-MEAM-ACE/07_MDMC-CU-SEGREGATION-T300-COEX-OFFSETS/analysis-2026-11-21_dxa-latest/
  (dxa-cu-near-partials.py, mclean-vs-dxa.py), .../07_MDMC-CU-SEGREGATION-T300-COEX-OFFSETS/check-convergence-mdmc300K.py,
  .../10_ACE-BULK-CU-CLUSTERS-NEAR-COEX/harvest-thread10.py; v6.2 enrichments in REPORT_2026-11-21_status-deck-v6
- TURLO: ~/Desktop/SIMULATIONS/TURLO-REDO/50_TURLO-METHOD-VS-OURS-CU-ZR/analyse-thread50.py,
  ~/Desktop/SIMULATIONS/TURLO-REDO/51_COMPLEXION-ONSET-DMU-SCAN-CU-ZR/analyse-thread51.py; isotherm in TURLO-REDO/32_BULK-X-OF-DMU-AND-A-AT-T
- 0 K site energies for McLean: Ni-Cu thread 09, Nb-Au thread 06 (folders 09_SEGREGATION-ENERGY-0K-SINGLE-CU, 06_SEGREGATION-ENERGY-0K-SINGLE-SOLUTE).
Known defects of the old scripts (see DEFINITIONS.md section 11): Ni-Cu distance to the mean line; Nb-Au orthogonal box for tilted cells;
the Ni-Cu "PTM" labels are DXA's internal labels (read after DXA; verified 22 Nov on d60); old Nb-Au comment "R2 ~ 4 for D0_22" is wrong (L1_2: 4, D0_22: 2.7).

## Testing
- pytest. Synthetic tests without OVITO: perfect fcc cell; straight analytic line -> exact distances; triclinic image test (solute near the
  y boundary of a tilted cell); known solute placements -> exact c, Gamma, alpha, cluster sizes; coverage and fallback axis on constructed
  cases; drift test on synthetic flat/drifting series; block averaging on synthetic correlated series.
- OVITO tests: GitHub Actions (ubuntu x86_64, pip ovito) and natively on the M5. Needs one small dissociated-dislocation cell (~20k atoms),
  made ONCE with LAMMPS on cmmg -- prepare the input, Erik runs it.
- Regression (on M5/cmmg, not in CI -- data too big): reproduce Nb-Au time-avg_07_2026-11-21.dat and Ni-Cu v6.2 where definitions agree;
  document every number that changes and why. TURLO regression pair: t50 R4 final (polyline mode) vs t50 R3 final (mixed mode).

## Releases
- v0.1.0: Stages 0-3, regions, Gamma(R), radial profiles, Warren-Cowley, R2, clusters, drift test + block averaging, provenance, regression report;
  one demo run per project with manifests in examples/manifests/; a short "how to call disegreg" note per project instance; DEFINITIONS.md final.
- v0.2: Stage 4 (uptake exponent, McLean, hysteresis), along-line Fourier analysis, ordered-domain sizes (Ni-Al vs 2.7 nm), plotting module.
- Zenodo DOI for the release the paper cites (as for openbis-mcp-server).
