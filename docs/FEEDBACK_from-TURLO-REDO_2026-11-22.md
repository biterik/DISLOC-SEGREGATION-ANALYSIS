# Feedback to disegreg from the TURLO-REDO instance (Cu-Zr, Mendelev 2009, 800 K) -- 22 Nov 2026

Written at Erik's request directly into this repository for the ANALYSIS-TOOL instance. Refers to docs/DEFINITIONS.md, set D1 (draft
2026-11-22). New file only; DEFINITIONS.md was not touched. Two problems that block the TURLO-REDO numbers for the paper, and one smaller wish.
They are not TURLO-specific: they will hit every Ni+X system as soon as a real line defect phase forms, i.e. exactly the case the paper is about.

## Problem 1 -- D1 has no core when DXA finds no partial lines (sections 3, 4, 6)
In D1 everything spatial hangs on the DXA polylines: r (section 3), the glide-plane height z_gp = median Z of the line points (section 3),
the regions C, F, S, T, B (section 4), and therefore Gamma(R) (section 6). In the defect-phase state the partials are gone and DXA finds
only part of the line, often not as Shockley partials:

| run (d90 cell, line length 108 A) | state | DXA line found |
|---|---|---|
| t50 R4 (coex - 0.5 kT) | intact stacking-fault ribbon | 2 x 108 A, two Shockley partials -- D1 works |
| t50 R3 (coex + 3.25 kT) | Zr-rich core complexion | 62 of 108 A, perfect 1/2[10-1] segments + one Shockley fragment |
| t51 D20 (coex + 2 kT, from R3) | complexion kept | 160 A in total, 110 A of it perfect 1/2[10-1] edge (constricted core) + Shockley fragments |
| t50 R1 (Turlo-Rupert dipole, 2 x 90 A) | complexion | 59 of 180 A, fragments |
| t50 R5 (Mendelev 2019, 4 ns) | precipitate on the line | 64 of 108 A, short Shockley fragments only |

Consequences in D1 as written:
- Where DXA finds nothing, the atoms of the defect phase get a large r: they fall into S/T/B, Gamma_C is underestimated and the far field
  is contaminated. The error is largest exactly in the states the paper wants to detect.
- z_gp from the median of fragments can be off the actual core layer.
- The "line length" in Gamma = .../L_x stays correct (L_x is the cell vector), but the region assignment does not.

Proposal (a new definition item, flagged per snapshot):
- Line coverage f_DXA = (DXA length assigned to this dislocation) / (expected length: 2 L_x for a dissociated line). Report it always.
- If f_DXA < a threshold (proposed 0.9) -> fallback core axis from the structure: per slab along X (e.g. 5 A), the centroid (Y, Z) of the
  largest connected cluster of non-fcc interior atoms (a-CNA), smoothed along X; r = distance to that axis with the same periodic-image
  treatment as for polylines. z_gp = median Z of that axis.
- Mixed option for partial coverage: r = min(distance to DXA polyline, distance to fallback axis).
- Output records which definition each snapshot used (polyline / fallback / mixed), so a time series that switches can be seen.
- Regression/test pair: t50 R4 final (D1 works) vs t50 R3 final (needs the fallback), same cell, same potential.

## Problem 2 -- x_ref when no bulk control exists at that dmu (section 5)
Section 5 allows x_iso only when the bulk control holds too few solute atoms. TURLO-REDO has exactly one bulk control, t50 R1b (VC-SGC,
c = 0.012, about coex + 3 kT, no precipitation in 2.2 ns). The scan in thread 51 (coex - 0.5, 0, +1, +2 kT, bottom-up and top-down) has none,
and running new controls was not part of the TURLO brief.
Proposal: allow x_iso also when no bulk control exists at that dmu (reason recorded in the manifest), i.e. rule 2 becomes
"no bulk control at this dmu, OR the control holds < 100 solute atoms". TURLO-REDO has the isotherm x(dmu, 800 K) of the perfect crystal for
Mendelev 2009 in TURLO-REDO/32_BULK-X-OF-DMU-AND-A-AT-T. c_B / x_ref is then the in-cell check, as in D1.

## Smaller wish -- icosahedral label (section 2)
The Cu-Zr core complexion contains an icosahedrally coordinated Cu layer between two Zr-rich layers. a-CNA labels these atoms "other",
so the region labels cannot distinguish a defect phase with a definite local structure from a disordered core. An optional extra per-atom
column (PTM with the icosahedral template enabled, RMSD cutoff 0.12 -- what TURLO-REDO used) would let the paper show this. Not needed for
the composition numbers.

## TURLO-REDO data you can use (on Erik's M5 once his rsync of 22 Nov has finished)
- d90 cells: X = [1-21] (line, periodic), Y = [10-1] (periodic), Z = [111] (free (111) faces, non-periodic); about 587 000 atoms; type 1 Cu,
  type 2 Zr; quenched snapshots every 100 ps, about 20 per run.
- Runs: t51 U00, U10, U20 (from the clean cell, coex + 0 / 1 / 2 kT), Dm05, D00, D10, D20 (from R3, coex - 0.5 / 0 / 1 / 2 kT);
  t50 R3 (coex + 3.25 kT), R4 (coex - 0.5 kT), R1b (bulk control, VC-SGC).
- Paths: ~/Desktop/SIMULATIONS/TURLO-REDO/51_COMPLEXION-ONSET-DMU-SCAN-CU-ZR/t51_*/ and
  ~/Desktop/SIMULATIONS/TURLO-REDO/50_TURLO-METHOD-VS-OURS-CU-ZR/t50_R{3,4,1b}_*/ (file names *_quenched_step*.dump).
- dmu values and coexistence: dmu_coex(800 K) = -3.6977 eV (TURLO sign convention: dmu_ours = mu_Zr - mu_Cu; LAMMPS fix sgcmc uses the
  opposite sign), kT = 0.068936 eV. TURLO-REDO will write the run manifests as soon as the manifest schema is released.
