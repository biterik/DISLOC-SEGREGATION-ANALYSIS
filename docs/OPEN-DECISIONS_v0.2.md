# disegreg v0.2 -- open decisions (the one page to read; branch v0.2-dev)

State: 24 Nov 2026, written by the unattended Claude Code session (docs/PLAN_v0.2.md). Nothing here is merged into main,
nothing is tagged, D1 rev. 4 is unchanged (D1.yaml byte-identical). All proposals live in docs/DEFINITIONS.md, section
"[PROPOSAL v0.2]", their parameters in disegreg/definitions/proposals_v0.2.yaml (label P-v0.2).
**Update 25 Nov:** Erik decided (section 5). Done: v0.1.1 released from main (0.1); P1-P4 and E1 are definition set D2
(docs/DEFINITIONS.md sections 12-16, disegreg/definitions/D2.yaml); P5 stays a proposal. What was done, the code decisions
taken on the way and the two new open points (N1, N2): section 6. Sections 0-4 below are the record of 24 Nov.

Abbreviations: SEM = standard error of the mean (from block averaging); DXA = dislocation extraction algorithm (OVITO);
Delta = offset (dmu - dmu_coex) / kT, negative = below coexistence; c_C = solute fraction in the core (r < 5 A);
x_ref = reference solute fraction: c_B (in-cell far field), x_iso (isotherm), x_bc (bulk-control cell).

## 0. Bugs in v0.1.0 (may need a v0.1.1 on main, with Erik)
**0.1 Snapshots after the end of the MC stream are silently dropped (v0.1.0 behaviour).** DECIDED 25 Nov: option (a),
released as v0.1.1 (section 6).
Stage 3 takes the averaging window from the MC stream; a snapshot later than the last stream line lies outside every window
and is not used -- without any message. Found on the Mac copies of the Nb-Au Mendelev runs: the stream ends at 1.53 ns
(rsynced earlier), the dumps go to 3.01 ns, so about half of the snapshots are ignored. No number is wrong for the snapshots
it does use, and on cmmg (stream and dumps from the same run folder) it does not happen. On the branch: a WARNING and the count
(`snapshots_after_window_signal` in the Stage-3 provenance); no number changed.
Options: (a) keep v0.1.0 behaviour + the warning (branch state; recommended for v0.1.1 = warning only, so no v0.1.0 number
changes); (b) when the stream ends before the last snapshot, take the window from the Stage-2 key quantity instead (changes
numbers for such runs); (c) refuse to run Stage 3 in that case. Recommendation: (a) now, and ask the projects to rsync the
stream together with the dumps (one line in docs/HOWTO_PROJECTS.md).

Two robustness points, fixed on the branch only (no v0.1.0 number affected):
- `results.write_table` would have written a numpy number as the text "np.float64(0.25)" (numpy >= 2). v0.1.0 always passes
  plain Python numbers, so no v0.1.0 table is affected; the fix writes plain Python numbers identically.
- `paper_flag` (ACE Ni-Cu held back) was dropped as soon as a manifest NAMED a bulk control, even before its x_bc existed.
  No v0.1.0 manifest names a bulk control yet, so nothing is affected. v0.2 drops the flag only when x_bc is actually used.

## 1. What was done (short)
| item | state | commit |
|---|---|---|
| A. Stage 4 across runs | code + tests done; demo on 59 runs (section 3) | see git log |
| B. plotting | `disegreg plot` gamma, timeseries, delta, uptake, regions, alongline; checked by eye on real Nb-Au data | see git log |
| C. along-line Fourier | code + tests done (proposal P4); demo numbers with the demo run | see git log |
| D. ordered domains Ni-Al | code + tests done (proposal P5); synthetic only -- no Ni-Al data on the Mac | see git log |
| E. project feedback | no FEEDBACK file after 24 Nov; the projects' first manifests reviewed (E1-E3 below) | see git log |

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

## 2. Decisions needed (24 Nov) -- P1-P4 and E1 DECIDED 25 Nov as recommended -> D2 sections 12-16; P5 still open

### P1. McLean: which 0 K sites count as core sites?  (DECIDED: option 1 -> D2 section 12)
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
Demo numbers (measured / McLean per run): section 3.

### P2. Hysteresis: when are bottom-up and top-down "different"?  (DECIDED: factor 2 -> D2 section 13)
Proposal: |difference| > 2 x sqrt(SEM_up^2 + SEM_down^2), for Gamma(20 A) (every reference) and c_C. Recommended. Alternative:
3 x (fewer false alarms, ~0.3 %), or a verdict only from Gamma(20 A). No demo numbers: no top-down run on the Mac (Lee-Shim
top-down "A2" approved 22 Nov, not yet run).

### P3. Error bar of Gamma vs x_bc  (DECIDED: as proposed -> D2 section 14)
Proposal: total error = sqrt(SEM_Gamma^2 + (N(r < R) / L_x x SEM(x_bc))^2), row `Gamma_error_total`. Recommended (otherwise the
bulk-control uncertainty is silently dropped). For Gamma(20 A), N(r < 20 A) / L_x ~ 180 atoms per A, so SEM(x_bc) = 1e-5 adds
0.0018 atoms/A.

### P4. Along-line core composition: continuous line state or periodic particles?  (DECIDED: 5 A, C only, p < 0.01 -> D2 section 16)
Proposal (DEFINITIONS P4): core atoms (r < 5 A) in X bins of ~5 A; dispersion index D (1 = random, >> 1 = particles,
< 1 = more even than random); spectrum S(q) normalised to binomial noise (mean 1), window-averaged; "periodic" if the highest
peak has p < 0.01 under binomial noise. Synthetic checks: a particle array with 25 A spacing in a 100 A cell gives q = 4,
lambda = 25.0 A, p < 1e-6, D > 5; random placement gives mean S = 1.00, mean D = 1.00, and 0.86 % false peaks at p < 0.01 (5.3 % at p < 0.05; 5000 trials;
the mode q = n/2 is treated with its one degree of freedom -- without that the false rate was 1.7 %).
Options: bins 5 A (recommended: ~40 core atoms per bin in the Ni cells, wavelengths 10-100 A) or 2.5 A (one (110) period at
d0/d60, noisier); atoms C only (recommended) or C + S; threshold p < 0.01 (recommended) or 0.05.
Demo numbers: section 3.

### P5. Ordered-domain sizes (Ni-Al vs 2.7 +- 0.6 nm)  (OPEN: stays a proposal until Ni-Al data exist)
Proposal (DEFINITIONS P5): ordered solute = no solute 1st neighbour and >= 4 of 6 solute 2nd neighbours; domain = 2nd-shell-
connected ordered solute atoms (>= 4) plus their solute 2nd neighbours (the surface); d = (6 n a^3 / pi)^(1/3).
Synthetic test, L1_2 Ni3Al sphere of diameter D in fcc Ni (a = 3.52 A), measured d / D:

| D (nm) | 1.6 | 2.0 | 2.7 | 3.2 | 4.0 |
|---|---|---|---|---|---|
| ordered atoms only (no surface) | 0.80 | 0.84 | 0.94 | 0.91 | 0.94 |
| with the surface atoms (recommended) | 1.05 | 0.99 | 1.02 | 1.00 | 0.99 |

The same sphere in a random Ni-13Al matrix: d / D within 10 %, no spurious domain in the matrix (random Ni-13Al: < 2 % of the Al
"ordered", no domain of >= 4). Options: k = 4 (recommended; counts D0_22-type order too) or k = 5 / 6 (L1_2 only, smaller
domains); with or without the surface atoms (recommended: with -- otherwise a 2 nm domain is reported as 1.7 nm).
Blocked: no Ni-Al MD/MC snapshot exists yet (the Ni-Al project folder on the Mac is empty).

### E. Project feedback (item E)
No docs/FEEDBACK_from-*.md and no feedback file in the project folders after 24 Nov. What the projects did on 24 Nov: each wrote
its first manifest (Ni-Cu d60 Lee-Shim -0.10 kT, Nb-Au d60 Mendelev -0.25 kT, TURLO t51 D20) -- all three pass
`disegreg manifest check` with the v0.2 code; TURLO-REDO has a manifest generator (52_.../DISEGREG/make-manifests.py) that
the v0.2 demo reuses read-only (our copy of D20 is identical to theirs). No bug report. Points that need you:
- **E1. Is a SMALL bulk control a bulk control?** (DECIDED: option (a) -> D2 section 15) Ni-Cu thread 11 (Lee-Shim, prepared 24 Nov) and thread 10 (ACE) are 32 000-atom
  periodic cells; D1 section 5 defines the bulk control as "same size, same dmu, same run length" as the ~570 000-atom dislocation
  cell. Its manifest notes "switch x_ref_source to bulk_control when thread 11 is harvested". Options: (a) accept a smaller cell
  for x_bc if it has the same per-atom swap rate and run length and >= 100 solute atoms in its interior (Lee-Shim at -0.10 kT:
  32 000 x 7.9e-3 = 250 Cu, enough), and record the size in the manifest -- recommended for x_bc; (b) keep "same size" strictly:
  thread 11 then counts as an isotherm point, not as x_bc. Note: diagnostic (iv) (does the dislocation-free crystal precipitate?)
  is weaker in a 32 000-atom cell (fewer nucleation sites, smaller clusters fit), so (a) should not be used for (iv).
- E2. TURLO-REDO's convergence item (thread 52, SI) used the rev. 3 drift rule (|trend| < 2 SEM); D1 rev. 4 uses 2 SE(change).
  Their own note says "recompute with disegreg when released" -- nothing for you, a reminder for the TURLO instance.
- E3. The Mac copies of the Nb-Au MC streams are older than the dumps (section 0.1): ask the Nb-Au instance to rsync streams and
  dumps together.

## 3. Demo numbers (24 Nov, M5; 59 runs; regression/results/v0.2-demo/, figures in fig/)
Input: DEMO manifests (examples/stage4/) for the Nb-Au thread-07 runs (Mendelev 2022 Nb, EAM-X Au; 4 characters x 3 offsets),
the Ni-Cu thread-07 runs with snapshots on the Mac (Lee-Shim MEAM, ACE), and the TURLO-REDO t51 bottom-up / top-down runs
(their own manifest generator, read-only). Snapshots: what is on the Mac (Ni-Cu: mostly 1-2 per run, so most Ni-Cu errors are
UNKNOWN; Nb-Au 15-18 per run, but see 0.1: in all 24 Nb-Au runs 7-8 snapshots lie after the end of the stream copy and are not
used; TURLO 20). 35 of 59 windows fail the drift test (mostly "last 30 % not flat" on the short Mac series) -- demo numbers,
not paper numbers. Reproduce: `python regression/v02_demo_summary.py` after the `disegreg analyse` call in its docstring.

**Uptake exponent n_up (vs x_iso; vs c_B in brackets), Lee-Shim MEAM Ni-Cu, bottom-up/seeded series** (fig/uptake_LeeShim2004-MEAM.png):

| offsets (kT) | screw d0 | 30 deg | 60 deg | edge d90 |
|---|---|---|---|---|
| -2.00 -> -1.00 | 1.21 (0.97) | 1.53 (1.16) | 1.60 +- 0.42 (1.16) | 2.41 +- 0.33 (1.75) |
| -1.00 -> -0.50 | 1.23 (1.58) | 0.89 (1.06) | 1.12 +- 0.25 (1.48) | 0.96 +- 0.31 (1.21) |
| -0.50 -> -0.25 | 0.18 (0.19) | 0.99 (1.16) | 1.17 +- 0.60 (1.17) | 1.06 +- 0.48 (1.15) |
| -0.25 -> -0.10 | 2.47 (2.36) | 3.22 (3.13) | **5.55 (4.17)** | **4.78 (4.09)** |

No error = one of the runs has a single snapshot on the Mac. The cooperative uptake near coexistence (60 deg and edge, "exponent 3-4"
in the paper plan) is reproduced with D1 definitions; between -2 and -0.25 kT n_up ~ 1 (independent sites). Mendelev Nb:
n_up 0.0-2.1, mostly < 1 (Nb near-saturated core, x ~ 0.07); EAM-X Au vs x_iso 0.3-1.4; EAM-X vs c_B is noise (far field
holds a handful of Au atoms: errors up to +-2900) -- shown, not interpreted. ACE (paper hold): 4.2-6.0 vs the NOMINAL x_iso of
the -0.10 kT row (marked "not used" by the project), 1.1-2.0 vs c_B -- an example of why D1 holds ACE until its bulk control.

**McLean (P1), measured c_C / predicted, vs x_iso** (summary_mclean.csv): Lee-Shim 0.2-0.9 far from coexistence, rising to
1.04-1.52 at -0.10 kT (d30 1.04, d0 1.32, d60 1.38, d90 1.52): the measured core overtakes the independent-site prediction
exactly where n_up > 1. Mendelev Nb 0.09-0.14 and EAM-X Au 0.06-0.32: McLean OVERpredicts (degenerate sites carry the deepest E;
Nb at x ~ 0.07 is not dilute) -- with the degenerate-site caveat P1 is an upper bound, so ratio < 1 is uninformative and ratio > 1
is the signal. ACE -0.10 kT 2.5-15 (paper hold, nominal x_iso).

**Hysteresis (P2), TURLO-REDO Cu-Zr d90, 800 K, bottom-up (from Cu) vs top-down (from the R3 complexion)** (summary_hysteresis.csv):

| Delta (kT) | Gamma(20 A) up - down (atoms/A), vs x_iso | |d| / s | c_C up - down | verdict (2 s) |
|---|---|---|---|---|
| 0.00 | -0.007 +- 0.009 | 0.8 | -0.0005 +- 0.0004 | not different |
| +1.00 | -0.046 +- 0.038 | 1.2 | -0.0022 +- 0.0018 | not different |
| +2.00 | **-1.40 +- 0.04** | 38 | **-0.124 +- 0.004** | **different** |

This is TURLO-REDO's own thread-51 finding (top-down: the complexion dissolves at 0 / +1 kT and persists at +2 kT), now with D1
definitions and an error criterion. The verdict does not depend on the factor 2 vs 3.

**Along-line spectrum (P4)** (summary_alongline_domains.csv, fig/alongline_*): ACE Ni-Cu: dispersion D = 4.5-9.4 and spectrum peaks
far above noise (p < 1e-5; d60 -0.10 kT: S = 36 at q = 2, i.e. two Cu-rich stretches 50 A apart) -> particles along the line
(single snapshots). Lee-Shim: D = 0.8-2.2, at -0.10 kT no significant peak (d0-d90: p = 0.09-0.58) -> an evenly filled
core, the continuous line state. Mendelev Nb: D = 0.24-0.57 -> MORE even than random (Nb-Nb avoidance along the core). EAM-X Au:
D ~ 1 (random). TURLO D20 complexion: D = 0.65, no peak -> continuous. Two TURLO runs with an empty core in some snapshots give
D = NaN (window mean over NaN); a note, not a bug.

**Ordered domains (P5)**: negative control only -- 0 domains in every run (Ni-Cu, Nb-Au, Cu-Zr), ordered-solute fraction 0 except
Mendelev Nb core 0.3-0.4 %. Ni-Al data needed.

**Site tables (P1)**: summary_site_tables.csv -- as in the P1 table above (now with the line period from the 0 K cell's own lattice
parameter: Mendelev cells were first given 39 instead of 40 periods because the MD a_cell is 1.4 % larger -- fixed).

## 4. Blocked / waiting
- Hysteresis demo: needs a top-down run (Lee-Shim A2) on the Mac.
- x_bc demo: needs a harvested same-size bulk control (Ni-Cu thread 11 prepared 24 Nov; ACE thread 10 cells are 32 000-atom
  isotherm cells with Cu-only dumps, not usable as x_bc input).

## 5. Erik's decisions (25 Nov 2026, 15:24, given in the claude.ai design session)
- 0.1: release v0.1.1 on main with option (a): WARNING + count only (`snapshots_after_window_signal`), no v0.1.0 number changes;
  plus one line in docs/HOWTO_PROJECTS.md: rsync the MC stream together with the dumps. Port only that fix (and the two
  robustness fixes if they do not change any v0.1.0 output) from v0.2-dev to main; tag v0.1.1, GitHub release (Zenodo makes the
  new version DOI automatically), CITATION.cff version/date/DOIs updated afterwards.
- P1: option 1 (r < 5 A of the 0 K DXA line, tabulated sites only, coverage reported, < 90 % flagged). EAM-X d90 flagged;
  ask the Nb-Au instance (via Erik) for a corrected 0 K site table.
- P2: |difference| > 2 x sqrt(SEM_up^2 + SEM_down^2), for Gamma(20 A) and c_C.
- P3: Gamma_error_total including the x_bc uncertainty, as proposed.
- P4: as recommended: 5 A bins, core atoms (C) only, "periodic" if p < 0.01.
- E1: option (a): a smaller bulk-control cell is accepted for x_bc (same per-atom swap rate and run length, >= 100 solute atoms in
  its interior, size recorded in the manifest), but NOT for diagnostic (iv).
- P5: NOT decided -- stays [PROPOSAL v0.2] until Ni-Al MD/MC data exist; its columns keep status = proposal.
- Label for the new definitions: still open (Erik asked in the same session: D2 recommended = D1 rev. 4 unchanged + sections
  P1-P4 / E1; every D1 quantity identical in D2).

## 6. Done on 25 Nov (Claude Code session, after Erik's decisions)
### 6.1 v0.1.1 (main; Erik authorised work on main for this release only)
- Ported from v0.2-dev: the Stage-3 count `snapshots_after_window_signal` / `snapshots_total` + the WARNING in `disegreg analyse`;
  `write_table` plain floats; `paper_flag` lifted only when x_bc is actually given (no v0.1.0 manifest names a bulk control,
  checked also in the three project manifests in the project folders). docs/HOWTO_PROJECTS.md: rsync the stream with the dumps.
- Checked: the three example manifests (Ni-Cu d60 Lee-Shim, Nb-Au d60 Mendelev, TURLO t50 R3) analysed with v0.1.0 (from the
  tag) and with v0.1.1 give identical Stage-2, Stage-3 and Gamma(R) tables, 141 470 cells, only `tool_version` differs. The
  Nb-Au run prints "WARNING: 8 of 17 snapshots are later than the end of the mc_stream:nSol_interior signal and are NOT used".
  104 tests pass; GitHub Actions green on main and on the tag.
- Tag v0.1.1, GitHub release https://github.com/biterik/DISLOC-SEGREGATION-ANALYSIS/releases/tag/v0.1.1, Zenodo DOI of v0.1.1:
  10.5281/zenodo.24535067 (concept DOI unchanged), now in CITATION.cff (version 0.1.1, 2026-11-25) and README.md.
### 6.2 D2 (branch v0.2-dev; not merged, not tagged)
- disegreg/definitions/D2.yaml = every D1 section value for value (tested) + sections `mclean`, `hysteresis`, `bulk_control`,
  `alongline`; docs/DEFINITIONS.md: header (D1 and D2), pointer lines [D2] in sections 5, 7, 9, new sections 12-16 and a D2
  revision log; the [PROPOSAL v0.2] part holds only P5. proposals_v0.2.yaml holds only `ordered_domains`.
- Status: along-line, McLean, hysteresis verdict, Gamma_error_total = `definition` in D2 tables; `domain_*` (P5) = `proposal`.
- Tests: 150 pass (D2 = D1 + new sections; every D1 Stage-2 record bit-identical under D2 on a synthetic snapshot; Stage 4
  with D1 tables; bulk-control size check; the OVITO end-to-end tests with D1 and with D2).
### 6.3 Code decisions taken (inside the agreed design, reported)
- D2 is the default definition set on the branch (`--definitions D1` still gives the v0.1 outputs; tested end to end with OVITO).
  `disegreg across` takes the definition set of the tables it reads. Stage 1 is the same for D1 and D2 (same cache entries).
- With D1 tables Stage 4 computes only what D1 defines: uptake exponents and hysteresis differences (Gamma(20 A), c_C), no
  verdict, no McLean (`--sites` ignored with a message), no x_bc error. Along-line rows in D1-labelled tables (the 24 Nov demo)
  keep `status = proposal`.
- E1 in code: a bulk-control manifest must state `n_atoms` (schema stays 1; `disegreg manifest check` says so); the snapshots
  are checked against it. A run taking x_bc records `n_atoms_run`, `size_ratio`, `smaller`, `accepted`, `diagnostic_iv`
  ("usable" / "not usable: bulk control smaller ..."), the last MD step of both runs and `same_run_length` (equal last MD step;
  a difference is flagged, not refused). A bulk-control table without size (made before 25 Nov) is not accepted. The per-atom
  swap rate is not in the data: the project states it in `notes` (docs/HOWTO_PROJECTS.md).
- `mclean.min_coverage: 0.90` is a D2 parameter (the value that was hard-coded in the code as the flag threshold).
### 6.4 New open points
- **N1. Diagnostic (iv) is not implemented** (since D1 rev. 2). DEFINITIONS section 5 says "for Delta > 0 without a bulk control
  the results table marks diagnostic (iv) as 'not tested'"; neither v0.1.0/v0.1.1 nor v0.2 has such a column (found while
  writing D2 section 15). No number is affected. Options: (a) v0.2 adds a results-table column `diagnostic_iv` = "not tested"
  (Delta > 0, no bulk control) / "not usable" (bulk control smaller, D2 section 15) / "tested" (same-size bulk control; the test
  itself = the bulk control's largest solute cluster vs the dislocation cell's far field, to be defined); (b) only the text in
  the provenance (state now). Recommendation: (a) with "not tested" / "not usable" now, the test itself defined when the first
  same-size bulk control exists.
- **N2. v0.2-dev is not merged with main** (your instruction). main now has v0.1.1 and the DOIs in README / CITATION.cff; the
  branch has the same code fixes. Before a v0.2 release main has to be merged into v0.2-dev (conflicts expected only in
  CHANGELOG.md, CITATION.cff, README.md, CLAUDE.md). Nothing to decide now, a reminder.
- Still open from before: P5 (until Ni-Al data), the corrected EAM-X d90 0 K site table (Nb-Au project, via you), the Ni-Cu
  bulk controls (thread 11) to be harvested.

