# disegreg -- definitions (set D1, DRAFT 2026-11-22 rev. 2, for Erik's review)

Every number produced by `disegreg` carries the tool version (git hash) AND the definition-set label (D1, D2, ...).
A change to anything in this file creates a new label; the paper's Methods section quotes the label it used.
Values marked **[CONFIRM]** still need Erik's OK.

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
  OVITO convention (Burgers circuit FS/RH = finish-to-start, right-handed).
- Lattice parameter a = the value the cell was built/scaled with, a_cell (manifest), used for all neighbour-shell cut-offs.
- Chemical-potential difference: dmu = mu_solute - mu_host (eV), stored in the manifest in THIS convention whatever sign the LAMMPS
  input uses (fix sgcmc in TURLO-REDO uses the opposite sign; the manifest records the conversion).
- Chemical-potential offset: Delta = (dmu - dmu_coex(T)) / kT. With dmu = mu_solute - mu_host a larger dmu means more solute, so
  Delta < 0: below coexistence (undersaturated), Delta > 0: supersaturated. The manifest states which coexistence line (incoherent bulk
  phase or coherent clustering line). File names such as "0.50kT-below-coex" correspond to Delta = -0.50.

## 2. Structure labels
- Region labels (fcc / hcp / other): a-CNA on the quenched snapshot. Parameter-free.
- Line detection: DXA with input structure fcc. As far as known, DXA does its own internal structure identification and does not use a
  preceding PTM or a-CNA result. [TO VERIFY on one real snapshot in v0.1.0.] The a-CNA labels are computed separately and are the only ones
  used for regions.
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

### 3.3 Fallback core axis from the structure, used when f_cov < 0.9
- Line defect cluster: the largest connected cluster of interior non-fcc atoms (a-CNA hcp + other, 1st-shell connectivity, excluded layers
  of section 4 removed first -- otherwise free surfaces would be picked up).
- Axis: in X slabs of 5 A, the centroid (Y, Z) of that cluster's "other" atoms (core atoms, not the hcp ribbon), with periodic-image-aware
  averaging in Y, smoothed by a running mean over 3 slabs (periodic in X).
- r_axis = distance to that axis, same periodic-image treatment as 3.1.
- Limitation: a coherent fcc precipitate is "fcc" and invisible to this fallback; DXA still sees the line in that case.

### 3.4 Distance actually used, and what is reported
- f_cov >= 0.9: r = r_DXA (mode "polyline"). f_cov < 0.9: r = min(r_DXA, r_axis) (mode "mixed"; with no DXA segment at all this is the axis).
- r_DXA and r_axis are BOTH computed for every snapshot, and Gamma(R) is reported under both, so a time series that switches mode shows
  the size of the switch. Every snapshot records f_cov, f_Shockley and the mode.
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
- R_far = 60 A (decided 22 Nov). Keeps ~75 % of the interior of the ~250 x 250 A cells.
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

## 6. Composition quantities (per snapshot)
- c_k = N_sol,k / N_k for each region k and each radial bin; enrichment E_k = c_k / x_ref.
- Line excess Gamma(R) = [ N_sol(r < R) - x_ref N(r < R) ] / L_x, atoms per A of dislocation, with L_x = length of the periodic cell vector along
  the line (the whole dissociated dislocation, both partials together); interior atoms only, all structure labels.
  Decomposition: Gamma = Gamma_C + Gamma_F + Gamma_S + Gamma_T (same formula restricted to each region).
  Proposed headline numbers: Gamma(20 A) and Gamma(40 A). They contain the whole core state even if the axis is off by a few A, so they
  are far less sensitive to the distance definition (3.4) than the core/shell split. [PROPOSAL]
- Partial separation: with exactly two Shockley partials, the Y distance between them in X bins of 5 A; mean and standard deviation along the line.

## 7. Chemical order (per snapshot, per region of the central solute atom)
- Neighbour shells: 1st r < 0.854 a (midway between a/sqrt(2) and a; 12 neighbours), 2nd 0.854 a <= r < 1.112 a (midway between a and
  sqrt(3/2) a; 6 neighbours).
- Warren-Cowley parameter alpha_k = (p_k - c_ref) / (1 - c_ref), p_k = solute fraction among the k-th-shell neighbours of solute atoms.
  Reported with c_ref = x_ref and with c_ref = the region's own c_k. alpha = 0 random, < 0 solute avoids solute, > 0 solute clusters.
- 2nd-neighbour ratio R2 = <solute 2nd neighbours per solute> / (6 c_k): 1 random, 0 complete avoidance.
  Ordered reference values (every solute without solute 1st neighbours, at c = 0.25):
  L1_2 (fcc with the solute on the cube corners and Ni on the face centres, as Ni3Al): all 6 2nd neighbours solute, R2 = 4.
  D0_22 (the same arrangement with every second cube layer shifted, as Ni3V): 4 of 6, R2 = 2.7.
- Clusters: solute atoms connected by 1st-shell bonds. Reported per solute atom as the size of the cluster it belongs to, aggregated by that
  atom's region (no cluster-to-region assignment needed).
- Planned for v0.2: ordered-domain sizes (for Ni-Al, against the 2.7 nm L1_2 domains of Howard et al. 2025); along-line core composition in
  5 A bins with its Fourier spectrum (periodic array of particles vs continuous line state).

## 8. Time series, convergence, averages (per run)
- Primary drift signal: total solute count from the MC stream file (high time resolution); secondary: the Stage-2 quantities per snapshot.
- Drift test (as in the Ni-Cu / Nb-Au scripts): over the last 30 % of the run, linear trend x window length compared with 2 x the
  block-averaging SEM; "flat" if smaller.
- Block averaging (Flyvbjerg-Petersen): consecutive samples are grouped into ever longer blocks until the block means are uncorrelated;
  the SEM is taken from the scatter of those block means.
- Averaging window: the longest final stretch that passes the drift test, at least 30 % of the run; window start/end are reported.

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
| structure (Ni-Cu) | "PTM" -- probably DXA's internal labels (read after DXA) | a-CNA [to verify] |
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
