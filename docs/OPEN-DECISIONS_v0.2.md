# disegreg v0.2 -- open decisions (the one page to read; branch v0.2-dev)

State: 24 Nov 2026, written by the unattended Claude Code session (docs/PLAN_v0.2.md). Nothing here is merged into main,
nothing is tagged, D1 rev. 4 is unchanged (D1.yaml byte-identical). All proposals live in docs/DEFINITIONS.md, section
"[PROPOSAL v0.2]", their parameters in disegreg/definitions/proposals_v0.2.yaml (label P-v0.2).

Abbreviations: SEM = standard error of the mean (from block averaging); DXA = dislocation extraction algorithm (OVITO);
Delta = offset (dmu - dmu_coex) / kT, negative = below coexistence; c_C = solute fraction in the core (r < 5 A);
x_ref = reference solute fraction: c_B (in-cell far field), x_iso (isotherm), x_bc (bulk-control cell).

## 0. Bugs in v0.1.0 (may need a v0.1.1 on main, with Erik)
None found that changes a v0.1.0 number. Two robustness points, fixed on the branch only:
- `results.write_table` would have written a numpy number as the text "np.float64(0.25)" (numpy >= 2). v0.1.0 always passes
  plain Python numbers, so no v0.1.0 table is affected; the fix writes plain Python numbers identically.
- `paper_flag` (ACE Ni-Cu held back) was dropped as soon as a manifest NAMED a bulk control, even before its x_bc existed.
  No v0.1.0 manifest names a bulk control yet, so nothing is affected. v0.2 drops the flag only when x_bc is actually used.

## 1. What was done (short)
| item | state | commit |
|---|---|---|
| A. Stage 4 across runs | code + tests done; demo on 52 runs (Nb-Au, Ni-Cu) running | see git log |
| B. plotting | (next) | |
| C. along-line Fourier | (open) | |
| D. ordered domains Ni-Al | (open) | |
| E. project feedback | (open) | |

### A. Stage 4 (`disegreg across`)
- `disegreg across <Stage-3 tables or folders> --out stage4.csv [--sites sites.yaml]` -> one tidy table, every row with
  `status` = definition | proposal. Run sets = same project, potential, host, solute, character, T, coexistence line, ordered by
  Delta; "up" family = bottom-up + seeded starts, "down" family = top-down starts.
- Uptake exponent n_up = Delta ln c_C / Delta ln x_ref between neighbouring offsets, for every reference the two runs share,
  errors propagated (checked against a Monte-Carlo propagation to 3 %). ACE Ni-Cu rows keep `paper_flag`. Status: definition
  (D1 section 9).
- McLean prediction from the 0 K site-energy tables (adapter `segE0K-by-site` reads both Ni-Cu thread 09 and Nb-Au thread 06;
  a generic `csv` adapter for Ni-Al): status proposal (P1).
- Hysteresis bottom-up vs top-down at the same Delta: differences (definition) and the verdict (proposal P2). No top-down run is
  on the Mac yet -> tested on synthetic tables only.
- Automatic x_bc: a bulk control gets its own manifest (`cell_type: bulk_control`; character, frame, z_gp etc. not needed).
  `disegreg analyse` handles bulk controls first and writes x_bc per snapshot + its window mean; every run whose manifest names
  that bulk control takes x_bc from `<bulk_control>_stage3_D1.csv` (in --out or --bulk-dir). The D1 rule "fewer than 100 solute
  atoms in the bulk-control interior -> x_iso" is applied. No same-size bulk control is on the Mac yet (Ni-Cu thread 11 is
  prepared, not run) -> tested end to end with OVITO on a synthetic bulk cell.
- Code decisions taken (inside the agreed design, reported): schema stays 1 (only the optional field `cell_type` added);
  results table has a new last column `status`; Stage 1 accepts a periodicity override (`pbc`) for LAMMPS data files (the
  Nb-Au 0 K cells are data files, which OVITO would read as periodic in Z) -- the cache keys of all v0.1 entries are unchanged;
  a bulk control periodic in Z has no excluded layer, one with z faces excludes 15 A as in D1; two runs of one family at the same
  Delta are both used (every combination, flagged), never silently one.

## 2. Decisions needed

### P1. McLean: which 0 K sites count as core sites?  (DEFINITIONS [PROPOSAL v0.2] P1)
Proposal: a tabulated site is a core site if its atom lies within 5 A (the D1 core radius) of a DXA line of the relaxed 0 K
cell; table entries that repeat a site one line period further along X count once.
Demo, all 16 tables on the Mac (sites with r < 5 A of the 0 K DXA line, per line period):

| potential | d0 | d30 | d60 | d90 |
|---|---|---|---|---|
| Lee-Shim MEAM (Cu) | 30 of 32 (94 %) | 58 of 59 (98 %) | 33 of 34 (97 %) | 58 of 60 (98 %) |
| ACE (Cu) | 33 of 35 (94 %) | 57 of 60 (96 %) | 33 of 34 (97 %) | 58 of 61 (96 %) |
| Mendelev 2022 (Nb) | 33 of 35 (95 %) | 57 of 58 (99 %) | 33 of 34 (98 %) | 60 of 59 (101 %) |
| EAM-X (Au) | 33 of 35 (94 %) | 56 of 58 (97 %) | 33 of 33 (99 %) | **28 of 59 (47 %)** |

Findings: (a) the tables reach slightly beyond one line period (x extent 3.1 A for a 2.49 A period at d0), so 3-10 sites per table
are listed twice (identical E, |dE| = 0) -- merged; (b) EAM-X d90 misses half the core sites (the energy-based partial locator of
thread 06, misplaced by up to 8 A, see Ni-Cu thread 09 "FINDING 2") -> flagged, needs a project re-run or the option "bulk";
(c) "degenerate" sites (same E as >= 2 others: the partial glided to the solute during the relaxation) are 6-157 per table, i.e.
for most tables the majority -> c_C(McLean) is an upper bound of the dilute picture, as the old mclean-vs-dxa.py said.
Options:
1. **(recommended)** r < 5 A of the 0 K DXA line, missing sites left out ("tabulated_only"), coverage reported and < 90 % flagged.
2. Same, but missing core sites counted as bulk-like (E = 0, "bulk") -- lowers c_C(McLean) where coverage is poor.
3. The projects' own site sets (all tabulated sites, kind "core"/"core-p1/p2", up to r = 15-54 A) -- mixes core and shell sites,
   not comparable with the measured c_C (r < 5 A).
Demo numbers (measured / McLean per run): see section 3 when the demo run is done.

### P2. Hysteresis: when are bottom-up and top-down "different"?  (P2)
Proposal: |difference| > 2 x sqrt(SEM_up^2 + SEM_down^2), for Gamma(20 A) (every reference) and c_C. Recommended. Alternative:
3 x (fewer false alarms, ~0.3 %), or a verdict only from Gamma(20 A). No demo numbers: no top-down run on the Mac (Lee-Shim
top-down "A2" approved 22 Nov, not yet run).

### P3. Error bar of Gamma vs x_bc  (P3)
Proposal: total error = sqrt(SEM_Gamma^2 + (N(r < R) / L_x x SEM(x_bc))^2), row `Gamma_error_total`. Recommended (otherwise the
bulk-control uncertainty is silently dropped). For Gamma(20 A), N(r < 20 A) / L_x ~ 180 atoms per A, so SEM(x_bc) = 1e-5 adds
0.0018 atoms/A.

## 3. Demo numbers
(filled in when the demo run of 52 runs is finished)

## 4. Blocked / waiting
- Hysteresis demo: needs a top-down run (Lee-Shim A2) on the Mac.
- x_bc demo: needs a harvested same-size bulk control (Ni-Cu thread 11 prepared 24 Nov; ACE thread 10 cells are 32 000-atom
  isotherm cells with Cu-only dumps, not usable as x_bc input).
