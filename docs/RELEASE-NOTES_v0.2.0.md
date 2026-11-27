# disegreg v0.2.0 -- release notes (27 Nov 2026)

Shared analysis package for the paper *"Dislocation Defect Phases and how to find them"*.

**Definition set: D2, rev. 3** (26 Nov 2026) = **D1 rev. 4 unchanged** + sections 12-17 of
[docs/DEFINITIONS.md](https://github.com/biterik/DISLOC-SEGREGATION-ANALYSIS/blob/v0.2.0/docs/DEFINITIONS.md).
Every D1 quantity is defined and computed identically in D2: on the 59 demo runs (Ni-Cu, Nb-Au, Cu-Zr) all Stage-2, Stage-3
and Gamma(R) values of the D1 quantities are identical cell by cell under D1 and D2 (1 800 408 cells). D2 is the default; `--definitions D1`
gives the v0.1 outputs. Numbers quoted with this version carry the label D2.

## What is new
- **Stage 4, across runs** (`disegreg across`): run sets (same project, potential, host, solute, character, temperature and
  coexistence line) ordered by the offset Delta = (dmu - dmu_coex) / kT;
  - uptake exponent n_up = Delta ln c_C / Delta ln x_ref between neighbouring offsets, for every reference composition
    (in-cell far field c_B, isotherm x_iso, bulk control x_bc), with first-order error propagation (section 9);
  - McLean (independent-site) prediction of the core concentration from the projects' 0 K segregation-energy tables, core
    sites = r < 5 A of the 0 K dislocation line, coverage and "degenerate" sites reported (section 12);
  - hysteresis: bottom-up vs top-down runs at the same Delta, "different" if |difference| > 2 x combined standard error,
    for Gamma(20 A) and c_C (section 13);
  - total error of Gamma vs x_bc including the bulk control's own error (section 14).
- **Bulk-control runs** (`cell_type: bulk_control` in the manifest): analysed first, their x_bc passed on automatically to
  every run that names them. Every bulk control, whatever its size, is accepted only with the run's dmu and T, >= 100 solute
  atoms in its interior, a last MD step >= 0.95 x the run's and a swap rate within 10 % of the run's (section 15); a cell
  smaller than the dislocation cell is accepted under the same rules. A bulk control that is not accepted gives its x_bc only as a diagnostic (reference label
  `x_bc_diagnostic`, never a headline number, not used across runs); the run is then referred to x_iso and c_B, and the
  reason is recorded. New manifest fields `n_atoms`, `swap_attempts_per_atom_per_step` (manifest schema stays 1; every
  v0.1 manifest stays valid).
- **Along-line core composition** (section 16): core solute counts in 5 A bins along the line, dispersion index and spectrum
  normalised to binomial noise; "periodic" if the highest mode of the window-averaged spectrum has p < 0.01 under binomial
  noise: continuous line state vs periodic array of particles. Correlated snapshots are accounted for by the effective number
  of snapshots n_eff = n_s (SEM_naive / SEM_block)^2 from the block averaging of the core concentration over the same
  snapshots; without a plateau of that block averaging (always with fewer than 8 snapshots (2-7)) n_eff = 1 and the run is flagged.
- **Diagnostic (iv) column** (section 17): `diagnostic_iv` in the results table -- "not tested" for every run without a bulk
  control (any Delta), "not usable" with a smaller bulk control; with a same-size bulk control "test not yet defined" (P6).
- **Plots** (`disegreg plot gamma | timeseries | delta | uptake | regions | alongline`, optional dependency `plot` =
  matplotlib): fonts >= 13 pt, colour-blind-safe Okabe-Ito palette, symmetric colour scales for signed maps.
- Stage 1: periodicity override for LAMMPS data files (0 K cells); the cache keys of all v0.1 entries are unchanged.
- Results table: new columns `status` (definition | proposal) and `diagnostic_iv`.

## Still proposals (not part of D2; their numbers carry `status = proposal` and are refused as headline numbers)
- **P5, ordered-domain sizes** (L1_2 domains in Ni-Al, against 2.7 +- 0.6 nm of Howard et al. 2025): implemented, waiting for
  Ni-Al data.
- **P6, the diagnostic-(iv) test** itself (bulk-control clusters vs a random alloy of the same composition, and the drift of
  x_bc): written down only, not implemented.

## What is tested how
| feature | exact synthetic tests | on real project data |
|---|---|---|
| D2 = D1 for every D1 quantity | yes (bit-identical records) | yes: 59 demo runs, every cell identical |
| uptake exponent + errors | yes (exact values; errors vs Monte-Carlo propagation within 3 %) | yes: Lee-Shim MEAM Ni-Cu, Nb-Au (Mendelev, EAM-X) -- many Ni-Cu runs have only 1-2 snapshots on the Mac, so their errors are unknown |
| McLean prediction | yes (theta formula, site selection, periodic copies) | yes: all 16 0 K site tables of Ni-Cu and Nb-Au |
| hysteresis criterion | yes | yes: TURLO-REDO Cu-Zr bottom-up / top-down pairs at 0, +1, +2 kT |
| along-line spectrum and window test | yes (particle arrays; binomial noise: false "periodic" 0.86 % per snapshot and ~1 % per window of independent snapshots at p < 0.01; n_eff for uncorrelated, correlated (AR(1)) and too-short c_C series) | yes: 59 demo runs |
| plots | smoke tests | checked by eye on Ni-Cu, Nb-Au, Cu-Zr |
| **x_bc from a bulk control, acceptance rules (run length, swap rate, dmu / T, size), Gamma_error_total, diagnostic_iv** | yes, incl. OVITO end to end on synthetic cells | **no: synthetic only** -- no bulk-control run has been harvested yet |
| **ordered domains (P5, proposal)** | yes (L1_2 spheres 1.6-4 nm recovered to 1-5 %) | **synthetic only**; on real data only as a negative control (0 domains in Ni-Cu, Nb-Au, Cu-Zr) |

153 tests pass on macOS arm64 (M5) and Linux x86_64 (GitHub Actions). The demo numbers (docs/OPEN-DECISIONS_v0.2.md,
regression/results/v0.2-demo/) are illustrations from the snapshots on the Mac, not paper numbers.

## Unchanged
OVITO 3.17.2 exactly (environment.yml); v0.1 manifests, cache entries and the D1 numbers.

Full list of changes: [CHANGELOG.md](https://github.com/biterik/DISLOC-SEGREGATION-ANALYSIS/blob/v0.2.0/CHANGELOG.md).

Funded by the Deutsche Forschungsgemeinschaft (DFG) through Project A02 of the SFB 1394, Project ID 409476157.
