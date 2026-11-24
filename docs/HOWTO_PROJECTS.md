# How a project instance calls disegreg

For the Claude instances that own the simulation projects (Ni-Cu, Nb-Au, Ni-Al, TURLO-REDO). disegreg is used
read-only: never edit its code or definitions from a project; report problems to the analysis-tool instance.

## 1. Write one manifest per run (Stage 0)

```bash
disegreg manifest template > <run folder>/<run_id>.disegreg.yaml     # every field, with its meaning
disegreg manifest check <run folder>/<run_id>.disegreg.yaml           # must print OK
```
The manifest lives in the project folder (owned by the project). Rules that matter most:
- `dmu` in the convention mu_solute - mu_host; `dmu_lammps` / `dmu_lammps_sign` state how the LAMMPS input writes it.
- `frame`, `a_cell` and `z_gp` come from the cell-builder metadata, never typed from memory; say where in `*_source`.
- `x_ref_source`: `bulk_control` (with the run id) or `isotherm` (with `x_iso`); without a bulk control give the reason.
- Paths may be relative to the manifest file, so the same manifest works on the M5 and on cmmg if it sits in the run folder.
Demo manifests: `examples/manifests/` (Ni-Cu d60 Lee-Shim, Nb-Au d60 Mendelev, TURLO t50 R3).

## 2. Run the analysis

On the M5 (environment `disegreg`, see README):
```bash
conda run -n disegreg disegreg analyse <manifest(s)> --out <project analysis folder>
```
On cmmg (docs/CMMG.md): Stage 1 for many dumps as a job array (`disegreg slurm ...`, Erik submits), then
`/cmmc/ptmp/$USER/conda-envs/disegreg/bin/disegreg analyse <manifests> --out <folder>` (cheap: ~1 s per snapshot).

Output per run: `<run_id>_stage2_D1.csv` (one row per snapshot x quantity), `<run_id>_stage3_D1.csv` (window means with
block-averaging errors and drift tests), `<run_id>_gammaR_D1.csv` (the Gamma(R) curve, R = 5-60 A, for the supplement),
each with a `.provenance.json` (tool version, definition set, SHA-256 of the dumps, the manifest). Across runs:
`disegreg table <folders> --out results_D1.csv`.

Columns to know: `basis` = distance used (`r` the one actually used, `r_dxa` DXA only); `ref` = reference composition of
E, Gamma and alpha_xref (`c_B` the snapshot's own far field, `x_iso` the isotherm from the manifest, `x_bc` a bulk-control
value given with `--x-ref`). Without a bulk control every such number exists for both `c_B` and `x_iso` (DEFINITIONS 5).
`role` (results table): `headline` = Gamma(20 A) under r, `check` = Gamma(40 A), `curve` = other R.

## 3. Read the flags before quoting a number

- `mode` (Stage 2): `polyline` = DXA found the whole line; `mixed` = the fallback axis was used in some slabs; compare
  Gamma under `r` and `r_dxa` to see the size of the switch.
- `z_gp_flag`: manifest glide plane and the lines disagree by more than one (111) spacing -- check the manifest (or the
  dislocation moved out of its glide plane).
- Stage 3 `converged` and `signal_plateau`: without a block-averaging plateau the error bar is a lower bound. Drift
  criterion `change_se` (|change| < 2 SE(change)); the old rule's verdict is in `flat_rev3_sem`, for comparison only.
- `check_required`: the fallback axis was active in the window -> quote Gamma(40 A) (`role` check) next to Gamma(20 A).
- `paper_flag` = "no paper number before bulk control" (ACE Ni-Cu without a bulk control): do not put that run's
  enrichments or Gamma into the paper.
- Quote numbers with the definition set (D1) and the tool version from the provenance file.

## Project notes

| project | solute column of the MC stream | dmu sign in LAMMPS | remarks |
|---|---|---|---|
| Cu-AT-DISLOCS-Ni-EAM-MEAM-ACE (Ni-Cu) | `nCu_interior` | +1 (fix atom/swap) | builder metadata: 00_CORRECTED-CELLS-A0-3.52 gate reports (glide gap at a = 3.52 A: scale by a_cell / 3.52) |
| Nb-Au-AT-DISLOCS-Ni-ACE-EAM | `nSol_interior` | +1 (DMU = mu_Nb - mu_Ni) | EAM-X Au counts are tiny (x ~ 1e-5): x_ref from the isotherm (dilute law), not from a bulk control |
| TURLO-REDO (Cu-Zr) | `NZr_mcgroup` | -1 (fix sgcmc: DMU_SGC = -dmu) | free (111) z faces; complexion runs are `mixed`; PTM label 4 = icosahedral |

## v0.2 additions (branch v0.2-dev, not released -- for information; use v0.1.0 for numbers)
- **Bulk control runs**: write a manifest for the bulk-control cell itself with `cell_type: bulk_control` (character, frame,
  frame_source, z_gp, z_gp_source, x_ref_source, start_type are not needed) and name its run_id in the `bulk_control` field of
  every dislocation run at the same dmu. `disegreg analyse` takes bulk controls first and passes their x_bc on automatically
  (`--bulk-dir <folder>` if the bulk control was analysed into another folder). Whether a 32 000-atom cell counts as a bulk
  control is open (docs/OPEN-DECISIONS_v0.2.md, E1).
- **rsync the MC stream together with the dumps.** Stage 3 takes the averaging window from the stream; snapshots later than the
  end of the stream are not used (v0.2 prints a WARNING with their number).
- Across runs: `disegreg across <folder with *_stage3_D1.csv> --out stage4.csv [--sites sites.yaml]` (uptake exponent, McLean,
  bottom-up vs top-down). Plots: `disegreg plot gamma | timeseries | delta | uptake | regions | alongline ...`.
- Rows with `status = proposal` depend on an undecided v0.2 proposal: do not quote them in the paper.
