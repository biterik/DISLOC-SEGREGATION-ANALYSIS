# disegreg -- definitions (set D1, rev. 4 of 2026-11-24: final, released with disegreg v0.1.0)

Every number produced by `disegreg` carries the tool version (git hash) AND the definition-set label (D1, D2, ...).
A change to anything in this file creates a new label; the paper's Methods section quotes the label it used.
No [CONFIRM] or [PROPOSAL] items are left (24 Nov 2026).

Abbreviations: DXA = dislocation extraction algorithm (OVITO); a-CNA = adaptive common-neighbour analysis (OVITO);
CNA = common-neighbour analysis; PTM = polyhedral template matching; MD/MC = molecular dynamics / Monte Carlo;
SGC = semi-grand-canonical ensemble (fixed total atom number, fixed chemical-potential difference); SEM = standard error of the mean.

## 1. Input and conventions
- Input: QUENCHED (energy-minimised) snapshots only. Quenching moves atoms but swaps none, so compositions are those of the MD/MC state,
  while structure labels and DXA lines are free of thermal noise.
- Atom types: type 1 = host (Ni; Cu in TURLO Cu-Zr), type 2 = solute (Cu, Au, Nb, Al; Zr).
- Simulation frame: X = line direction (periodic), Y = in the glide plane, normal to the line (periodic), Z = glide-plane normal (non-periodic).
  The crystal directions of X, Y, Z are read from the cell-builder metadata via the run manifest, never typed in.
- Characters: d0 = screw, d30, d60, d90 = edge (angle between Burgers vector and line).
- Burgers vector b and line direction xi: taken from DXA in the simulation frame, rotated into the crystal frame; b and xi from the same segment;
  OVITO convention: left-hand start-finish (LH/SF, as stated in the OVITO manual), which gives the same b as the textbook
  right-hand finish-start (RH/FS, Hirth & Lothe) convention. b is reported with the line direction xi oriented along +X.
  The cell builders' nominal b (LEGO-TOOLS build_single_dislocation_cell.py; tests/data-generation/make_test_cell.py) is the
  NEGATIVE of the DXA b with xi = +X (checked 24 Nov: production d60 nominal 1/2[-1 1 0] -> DXA 1/2[1 -1 0]; test cell nominal
  1/2[1 -1 0] -> DXA 1/2[-1 1 0]). Concentrations do not depend on this sign; the paper states b and xi in the DXA convention.
- Lattice parameter a = the value the cell was built/scaled with, a_cell (manifest), used for all neighbour-shell cut-offs.
- Chemical-potential difference: dmu = mu_solute - mu_host (eV), stored in the manifest in THIS convention whatever sign the LAMMPS
  input uses (fix sgcmc in TURLO-REDO uses the opposite sign; the manifest records the conversion).
- Chemical-potential offset: Delta = (dmu - dmu_coex(T)) / kT. With dmu = mu_solute - mu_host a larger dmu means more solute, so
  Delta < 0: below coexistence (undersaturated), Delta > 0: supersaturated. The manifest states which coexistence line (incoherent bulk
  phase or coherent clustering line). File names such as "0.50kT-below-coex" correspond to Delta = -0.50.

## 2. Structure labels
- Region labels (fcc / hcp / other): a-CNA on the quenched snapshot. Parameter-free.
- Line detection: DXA with input structure fcc and OVITO's default parameters (recorded in D1.yaml). DXA does its own internal
  structure identification and does not use a preceding PTM or a-CNA result -- VERIFIED 2026-11-22 on Ni-Cu d60 (Lee-Shim MEAM,
  -0.25 kT, quenched step 830000, 568 280 atoms): the DXA lines are bit-identical with DXA alone, PTM -> DXA and a-CNA -> DXA, and the
  "Structure Type" read after PTM -> DXA equals DXA's own labels, not PTM's (it differs from a-CNA for 67 atoms, 0.012 %). The a-CNA
  labels are computed in a separate pipeline and are the only ones used for regions.
- Reproducibility (checked 24 Nov on one Ni-Cu ACE d90 snapshot, cmmg x86_64 vs M5 arm64, same OVITO 3.17.2 and parameters):
  positions, a-CNA and PTM labels bit-identical; DXA finds the same segments and Burgers vectors, one partial identical, the other
  displaced by 0.09 A on average (max 0.5 A, mean position 0.22 A, length 0.13 A). Line positions are therefore reproducible to a few
  tenths of an A; core concentrations of a single snapshot to ~1-2 %. Regression tests compare DXA-derived numbers with a tolerance.
  On ONE machine DXA is deterministic: four Stage-1 runs of the same snapshot on the M5 (23 and 24 Nov) gave bit-identical
  lines, so the cross-machine difference is floating-point arithmetic (x86_64 vs arm64), not thread order.
- Extra per-atom label, stored for every snapshot but NOT used for regions: PTM with the fcc, hcp, bcc and icosahedral templates,
  RMSD (root-mean-square deviation from the ideal template) cut-off 0.12. It distinguishes a core state with a definite local structure
  (e.g. the icosahedral Cu layer of the Cu-Zr complexion, which a-CNA calls "other") from a disordered core. [from TURLO-REDO feedback]

## 3. Distance to the line, r
### 3.1 From the DXA lines (default)
- r_DXA = shortest distance from an atom to the nearest DXA dislocation polyline, over ALL periodic images.
- Exact periodic treatment for triclinic cells: the polyline points are replicated by the full cell vectors (incl. the xy tilt) into the
  neighbouring images in X and Y, densified to 0.25 A spacing (max. error < 0.02 A at r >= 0.5 A), then queried with a non-periodic KD-tree
  (k-dimensional search tree).
- All DXA segments count (Shockley partials, perfect 1/2<110> segments, any other type); segment types and lengths are listed and
  non-Shockley segments flagged.

### 3.2 Line coverage [from TURLO-REDO feedback, modified]
- f_cov = fraction of the cell length L_x along X covered by at least one DXA segment (union of the segments' X intervals / L_x).
  Projected coverage, not total length / (2 L_x): a fully constricted line found as perfect segments has f_cov = 1, not 0.5.
- Also reported: f_Shockley = the same for Shockley segments only, and the total DXA length per type.

### 3.3 Fallback core axes from the structure [rev. 3]
- Line defect clusters: connected clusters of interior non-fcc atoms (a-CNA hcp + other + any other non-fcc label, 1st-shell
  connectivity, excluded layers of section 4 removed first -- otherwise free surfaces would be picked up).
  EVERY cluster that runs along the line counts: its core-type atoms ("other") occupy at least 50 % of the X slabs below (Erik 23 Nov).
  This keeps vacancies, isolated distorted solute sites and compact precipitates away from the line out, and gives a cell with two
  dislocations (dipole, e.g. TURLO t50 R1) two axes.
- One axis per such cluster: in X slabs of 5 A, the centroid (Y, Z) of the cluster's "other" atoms (core atoms, not the hcp ribbon);
  periodic in Y (minimum image around the circular mean, then the arithmetic mean); slabs without core atoms are interpolated;
  smoothed by a running mean over 3 slabs (periodic in X).
- r_axis = distance to all axes, same periodic-image treatment as 3.1 (diagnostic, always reported).
- Limitation: a coherent fcc precipitate is "fcc" and invisible to this fallback; DXA still sees the line in that case.

### 3.4 Distance actually used, and what is reported [rev. 3]
- Assignment of DXA to the axes: every DXA line point belongs to the cluster of its nearest interior non-fcc atom, if that atom is closer
  than the 1st-shell cut-off 0.854 a (else to none). Per cluster m and slab k: DXA coverage f_km = covered fraction of the slab (projected,
  as in 3.2, using only cluster m's DXA points).
- The axis of cluster m is ACTIVE only in the slabs where DXA covers less than half of the slab for that cluster, f_km < 0.5 (Erik 23 Nov).
  Where DXA has found the line, the axis is not used -- in a still-dissociated line the axis lies in the middle of the stacking fault, and
  using it there would turn fault atoms into core atoms.
- r, slab-wise: an atom whose X lies in an active slab of cluster m gets r = min(r_DXA, distance to axis m); every other atom r = r_DXA.
  (Slab-wise, not "distance to the active stretches": the end of a stretch would otherwise reach ~5 A into the neighbouring slab where
  DXA did find the partials.) Mode "polyline" if no slab is active, "mixed" otherwise ("axis" if there is no DXA line at all, "none" if
  there is neither a DXA line nor a line defect cluster).
- r_DXA and r_axis are BOTH computed for every snapshot, and Gamma(R) is reported under r_DXA and under r, so a time series that switches
  mode shows the size of the switch. Every snapshot records f_cov, f_Shockley, per cluster f_cov,m and the number of active slabs, and the mode.
- Glide-plane height z_gp: from the cell-builder metadata (manifest), checked against the median Z of the Shockley line points (or of
  the fallback axis); a difference larger than one (111) plane spacing is flagged.
- Regression pair: TURLO t50 R4 final (two complete Shockley partials, polyline mode) vs t50 R3 final (complexion, mixed mode), same cell
  and potential.

## 4. Regions (exclusive: every interior atom has exactly one code, so region excesses add up)
| code | name | definition |
|---|---|---|
| X | excluded | within 15 A of the lowest or highest atom in Z (fixed boundary layers or free surfaces) |
| C | core | r < 5 A |
| F | fault | a-CNA hcp, \|z - z_gp\| < 4 A, r >= 5 A (the stacking-fault ribbon; (111) spacing a/sqrt(3) ~ 2 A) |
| S | shell | 5 <= r < 10 A, not F |
| T | transition | 10 A <= r < R_far, not F |
| B | far field | r >= R_far, a-CNA fcc |
| O | far, other | r >= R_far, not fcc (reported; should be ~0) |
- R_far = 60 A (decided 22 Nov). Keeps ~75 % of the interior of the ~250 x 250 A cells (measured on Ni-Cu d60: 75.1 %).
- Precedence where two rows overlap: X > C > F > S, T > B, O. The table leaves open whether an hcp atom near the glide plane at
  r >= R_far is F or O; it is F (as for S and T). [confirmed by Erik 24 Nov; affects ~0 atoms, the fault ribbon ends well inside R_far]
- Optional extra column (not default): far_slab40 = a-CNA fcc, |z - z_gp| > 40 A, not X (the definition of the 21 Nov decks).
- Radial profile bins (A): 0-5, 5-10, 10-20, 20-40, 40-60, 60-80, 80-100, >= 100 -- the outer bins are the plateau check of the far field.

## 5. Reference composition x_ref
Order of preference, chosen per run in the manifest and always reported with its source:
1. x_bc: bulk control cell (no dislocation, same size, same dmu, same run length), interior atoms, time average with block-averaging error.
2. x_iso: the potential's own bulk SGC isotherm x(dmu, T) of the perfect crystal (project data), used when
   (a) there is no bulk control at this dmu (reason recorded in the manifest), or
   (b) the bulk control holds too few solute atoms (expected solute count in its interior < 100, e.g. EAM-X Au at x ~ 1e-5).
   Caveat: x_iso gives the reference composition but NOT protocol diagnostic (iv) (does the dislocation-free crystal precipitate at this
   dmu?). For Delta > 0 (supersaturated) without a bulk control the results table marks diagnostic (iv) as "not tested". [rule (a) from
   TURLO-REDO feedback]
3. c_B, the in-cell far field: diagnostic only. The ratio c_B / x_ref is reported for every run; a ratio clearly different from 1 is a result
   (e.g. long-range elastic enrichment by the dislocation's pressure field, or a cell-size effect).
Reporting rule (Erik, 24 Nov): until a run has its bulk control, enrichments and Gamma are reported against BOTH c_B and x_iso,
and no ACE Ni-Cu enrichment is quoted in the paper before its bulk control has decided between the two. With c_B as the reference,
each snapshot uses its OWN c_B (same cell, same snapshot); E and Gamma are then averaged over the window like every other
quantity (Erik, 24 Nov; Gamma is linear in x_ref, so this differs from using the window mean of c_B only by a small covariance
term, and the error bar includes the scatter of c_B). Reason (regression 23 Nov):
c_B / x_iso = 0.91-0.98 for Lee-Shim MEAM, but 1.6-2.6 for ACE at -0.10 kT (core enrichment 213-590 depending on the reference)
and 0.5-0.85 for ACE at -0.50 kT. If the bulk control also shows the ACE bulk clustering, c_B is the right reference.

## 6. Composition quantities (per snapshot)
- c_k = N_sol,k / N_k for each region k and each radial bin; enrichment E_k = c_k / x_ref.
- Line excess Gamma(R) = [ N_sol(r < R) - x_ref N(r < R) ] / L_x, atoms per A of dislocation, with L_x = length of the periodic cell vector along
  the line (the whole dissociated dislocation, both partials together); interior atoms only, all structure labels.
  Decomposition: Gamma = Gamma_C + Gamma_F + Gamma_S + Gamma_T (same formula restricted to each region).
  Headline number (Erik, 24 Nov): Gamma(20 A) -- the core, fault ribbon and near shell (~180 atoms per A of line in the Ni cells);
  less noise and ~3x less sensitive to an error in x_ref than Gamma(40 A). The full curve Gamma(R), R = 5-60 A, goes into the
  supplement. For snapshots in which the fallback axis is active (DXA misses part of the line, e.g. core complexions), Gamma(40 A) is
  quoted as the robust check: it hardly depends on the distance definition (TURLO t50 R3: Gamma(40) 4.73 vs 4.77, Gamma(20) 3.80 vs 3.42).
- Partial separation: with exactly two Shockley partials, the Y distance between them in X bins of 5 A; mean and standard deviation along the line.
  In each bin, Shockley polylines closer than 2 A to each other (Y-Z plane) are pieces of the same partial (a partial split into
  several DXA segments) and are merged; a bin counts if exactly two partials remain. The number of valid bins is reported. [2 A confirmed by Erik 24 Nov]
- Stored per snapshot: the COUNTS N_k, N_sol,k (regions, radial bins, inside R; under r and r_DXA). x_ref is usually a time average from
  another run, known only later; c_k, E_k and Gamma follow from the counts for any x_ref.

## 7. Chemical order (per snapshot, per region of the central solute atom)
- Neighbour shells: 1st r < 0.854 a (midway between a/sqrt(2) and a; 12 neighbours), 2nd 0.854 a <= r < 1.112 a (midway between a and
  sqrt(3/2) a; 6 neighbours).
- Warren-Cowley parameter alpha_k = (p_k - c_ref) / (1 - c_ref), p_k = solute fraction among the k-th-shell neighbours of solute atoms.
  Reported with c_ref = x_ref and with c_ref = the region's own c_k. alpha = 0 random, < 0 solute avoids solute, > 0 solute clusters.
- 2nd-neighbour ratio R2 = <solute 2nd neighbours per solute> / (6 c_k): 1 random, 0 complete avoidance.
  Ordered reference values (every solute without solute 1st neighbours, at c = 0.25):
  L1_2 (fcc with the solute on the cube corners and Ni on the face centres, as Ni3Al): all 6 2nd neighbours solute, R2 = 4.
  D0_22 (the same arrangement with every second cube layer shifted, as Ni3V): 4 of 6, R2 = 2.7.
- alpha is undefined (reported as NaN) for c_ref = 0 or 1, e.g. a core that is entirely solute.
- Central solute atoms: interior atoms (not X); their neighbours: all atoms, including those of region X.
- Clusters: solute atoms connected by 1st-shell bonds. Reported per solute atom as the size of the cluster it belongs to, aggregated by that
  atom's region (no cluster-to-region assignment needed).
- Planned for v0.2: ordered-domain sizes (for Ni-Al, against the 2.7 nm L1_2 domains of Howard et al. 2025); along-line core composition in
  5 A bins with its Fourier spectrum (periodic array of particles vs continuous line state).

## 8. Time series, convergence, averages (per run)
- Primary drift signal: total solute count from the MC stream file (high time resolution); secondary: the Stage-2 quantities per snapshot.
- Drift test (rev. 4, Erik 24 Nov): over the last 30 % of the run, fit a straight line; the change across the window
  (slope x window length) is compared with its own standard error SE(change), taken from the block-averaging SEM of the detrended
  series (SE(change) = sqrt(12) x SEM for evenly spaced samples). "Flat" if |change| < 2 SE(change): a truly flat series passes 95 %
  of the time. A passed test means "no drift detected", not proof of convergence; diagnostic (v) (bottom-up vs top-down) is the
  stronger check.
- Block averaging (Flyvbjerg-Petersen): consecutive samples are grouped into ever longer blocks until the block means are uncorrelated;
  the SEM is taken from the scatter of those block means.
- Averaging window: the longest final stretch that passes the drift test, at least 30 % of the run; window start/end are reported.
- History: rev. 3 used the rule of the old Ni-Cu / Nb-Au scripts, "|change| < 2 x SEM". That passes a truly flat series only ~44 %
  of the time (tests/test_timeseries.py), because the change has its own, ~3.5x larger standard error. The old verdict is still
  computed and stored next to the new one, for comparison with the 21 Nov decks.
- A block-averaging estimate without a plateau (correlation time comparable to the window) is flagged: its SEM is a lower bound.

## 9. Across runs
- Uptake exponent between neighbouring offsets: n_up = Delta ln c_C / Delta ln x_ref, errors propagated.
  Independent-site (McLean) segregation gives n_up <= 1; n_up > 1 indicates cooperative uptake.
- McLean prediction: each core site i with 0 K segregation energy E_i (project tables) is occupied with theta_i, where
  theta_i / (1 - theta_i) = x_ref / (1 - x_ref) exp(-E_i / kT); predicted c_C = mean theta_i, compared with the measured c_C.
- Hysteresis: bottom-up vs top-down runs at the same Delta, compared with block-averaging errors.
- Output: one tidy long table (run, quantity, region, value, error, definition set, tool version).

## 10. Provenance (in every output file)
Tool version (git hash), definition set, OVITO version (pinned), Python/numpy/scipy versions, SHA-256 of every input file, all parameters, host name.

## 11. What changes compared with the existing project scripts
| item | before | D1 |
|---|---|---|
| distance (Ni-Cu) | to the mean (y,z) of each line, periodic in y only | to the polyline, all images |
| distance (Nb-Au) | orthogonal periodic box (approximate for the tilted d0/d30/d60 cells) | exact triclinic images |
| far field (Ni-Cu) | r >= 80 A from the mean lines | r >= R_far from the polylines, fcc only |
| far field (Nb-Au) | fcc, > 40 A from the glide plane | r >= R_far (slab kept as optional column) |
| far field (TURLO) | > 60 A from every core | r >= R_far |
| structure (TURLO) | own conventional CNA, r_c = 3.14 A | a-CNA |
| structure (Ni-Cu) | "PTM" -- in fact DXA's internal labels (read after DXA; verified 22 Nov) | a-CNA |
| surface layers (TURLO) | 10 A | 15 A |
| shell (Nb-Au) | 5-10 A including hcp fault atoms | fault atoms only in F |
| core when DXA finds only part of the line | not handled (core atoms get a large r) | coverage f_cov + fallback axis (3.2-3.4) |
| z_gp | median Z of DXA points | cell-builder metadata, checked against DXA |
| 2nd-shell cut-off (Nb-Au) | 1.10 a | 1.112 a |
| reference | per-script (isotherm up-scan, in-cell far field) | x_ref hierarchy of section 5 |
Each change will be quantified on the regression runs (Nb-Au time-avg_07_2026-11-21.dat, Ni-Cu v6.2 enrichments).

## Revision log of this draft
- rev. 1 (22 Nov 18:20): first draft.
- rev. 2 (22 Nov): R_far = 60 A (Erik); dmu sign convention stated; TURLO-REDO feedback (docs/FEEDBACK_from-TURLO-REDO_2026-11-22.md):
  line coverage + fallback axis (coverage made projected, fallback restricted to core-type atoms and interior, both distances always
  reported, z_gp from builder metadata), x_iso also without a bulk control (with the diagnostic-(iv) caveat), PTM incl. icosahedral template as
  stored extra label. Proposed Gamma(20 A), Gamma(40 A) as headline numbers.
- rev. 3 (22 Nov, Erik accepted both proposals from the analysis-tool session): (1) the fallback axis is used only in the slabs where DXA
  misses the line (f_km < 0.5), not along the whole line -- otherwise mid-fault atoms of a partially found dissociated line become core;
  replaces the global switch "f_cov < 0.9 -> mixed"; (2) every line defect cluster that runs along the line (>= 50 % of the slabs) gets its
  own axis, not only the largest -- dipole cells; DXA points assigned to clusters by the nearest non-fcc atom. The axis is applied
  slab-wise (atoms in active slabs only), so it cannot reach into slabs where DXA found the line. Y centroid made exact
  (minimum image around the circular mean). Still label D1: no number has been produced with D1 yet.
  Same day: DXA parameters recorded (OVITO 3.17.2 defaults); the DXA-internal-labels question of section 2 verified (no change).
  23 Nov: Erik confirmed min_slab_fraction 0.5 and active_below_dxa_coverage 0.5 (values unchanged). Written down while implementing
  sections 4, 6, 7 (no value changed): region precedence (F before O) [CONFIRM]; merging of Shockley pieces for the partial
  separation, new parameter partial_merge_distance 2 A [CONFIRM]; counts stored per snapshot; alpha = NaN for a region with c = 0 or 1.
  23 Nov, section 8: drift-test criterion found statistically too strict -> PROPOSAL (not yet adopted; would be D2 or D1 rev. 4).
- rev. 4 (24 Nov, Erik): drift test |change| < 2 SE(change) adopted (section 8); x_ref reporting rule: both c_B and x_iso until the bulk
  control exists, no ACE Ni-Cu enrichment in the paper before it (section 5); headline Gamma(20 A), full Gamma(R) in the supplement,
  Gamma(40 A) as check where the fallback axis is active (section 6); region precedence F before O and partial_merge_distance 2 A
  confirmed (sections 4, 6). Documented: Burgers-vector sign relation between the cell builders and DXA, and the cross-machine DXA
  reproducibility (sections 1, 2). Label stays D1: no number for the paper has been produced with D1 yet.
  Same day: c_B as reference is taken per snapshot, not as a window mean (section 5); DXA deterministic on one machine
  (section 2). No value changed.
- final (24 Nov): D1 rev. 4 released with disegreg v0.1.0. Any further change creates D2.

---

# [PROPOSAL v0.2] -- open proposals, NOT part of D1 (written 24 Nov by the analysis-tool session, waiting for Erik)
Everything below is a proposal for disegreg v0.2. D1 rev. 4 above is unchanged. Each proposal gets a label (D1 extension or
D2) only when Erik decides it; until then every number that depends on it carries `status = proposal` in the output tables, and
the results table refuses such a number for role "headline". Parameter values: `disegreg/definitions/proposals_v0.2.yaml`
(label P-v0.2). Decisions, options and demo numbers: docs/OPEN-DECISIONS_v0.2.md.

## P1. McLean prediction: which 0 K sites are core sites (section 9)
- Input: a project's table of 0 K single-solute segregation energies E_i (one solute atom on site i of the relaxed 0 K pure-host
  cell with one dislocation; negative = attraction), sites identified by atom id, and that relaxed 0 K cell.
- Proposal: site i is a CORE site if its atom lies closer than 5 A (the D1 core radius) to a DXA line of the relaxed 0 K cell
  (Stage 1 with the D1 DXA parameters, all periodic images). The table is meant to cover one line period; sites that repeat an
  earlier entry shifted by a whole line period along X are the same site and count once.
- Predicted core concentration: c_C(McLean) = mean of theta_i over the core sites, with
  theta_i / (1 - theta_i) = x_ref / (1 - x_ref) exp(-E_i / kT), for every reference x_ref of the run (c_B, x_iso, x_bc).
  Reported: c_C(McLean), measured c_C / c_C(McLean), and the local exponent of the prediction d ln c_C / d ln x_ref (<= 1).
- Core sites missing from the table: option "tabulated_only" (mean over the tabulated core sites; default of the proposal) or
  "bulk" (missing core sites count with E = 0, theta = x_ref). The coverage (tabulated core sites / core atoms per line period of
  the 0 K cell) is always reported; below 90 % the rows are flagged.
- Caveat (Ni-Cu thread 09, 4 Nov): many tabulated sites relax to the same deepest state because the partial glides to the solute;
  their E_i is not their own. The number of such "degenerate" sites is reported. Then c_C(McLean) is an upper bound.

## P2. Hysteresis criterion (section 9)
- Pairs: a bottom-up and a top-down run of the same run set at the same Delta (|difference| < 0.001). Seeded runs are not compared.
- Compared: Gamma(20 A) under r (every reference both runs have) and c_C. Difference d = up - down, error s = sqrt(SEM_up^2 + SEM_down^2)
  (block-averaging SEMs of the two windows).
- Proposal: the two runs are DIFFERENT if |d| > 2 s (two independent normal errors: a false "different" in ~5 % of pairs).
  The difference and |d| / s are always reported, the verdict as 1 / 0.

## P3. Error of Gamma vs x_bc (section 5, code side)
- The Stage-3 SEM of Gamma vs x_bc holds only the time scatter of the dislocation run; x_bc is one number with its own SEM.
- Proposal: quote the total error sqrt(SEM_Gamma^2 + (N(r < R) / L_x x SEM(x_bc))^2) (Stage 4, row "Gamma_error_total").
