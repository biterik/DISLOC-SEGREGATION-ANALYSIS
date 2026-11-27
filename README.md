# disegreg

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24513776.svg)](https://doi.org/10.5281/zenodo.24513776)

Shared analysis of solute segregation at dislocations in atomistic MD/MC simulations (molecular dynamics / Monte
Carlo), written for the paper *"Dislocation Defect Phases and how to find them"*. One code and one versioned
definition set for all simulation projects, so that every number is computed the same way.

**Status: v0.2.0 (27 Nov 2026), Stages 0-4.** Definitions: [docs/DEFINITIONS.md](docs/DEFINITIONS.md)
(set D2, rev. 3, the default = D1 rev. 4 unchanged + sections 12-17; `--definitions D1` gives the v0.1 numbers).
Release notes: [docs/RELEASE-NOTES_v0.2.0.md](docs/RELEASE-NOTES_v0.2.0.md). Design: [docs/DESIGN.md](docs/DESIGN.md). Changes: [CHANGELOG.md](CHANGELOG.md).

## What it computes
From quenched LAMMPS snapshots with one dissociated dislocation: the distance of every atom to the dislocation line
(OVITO's dislocation extraction algorithm, DXA, with exact periodic images of triclinic cells, plus a structural fallback
when DXA finds only part of the line), exclusive regions (core, stacking fault, shell, transition, far field), solute
concentrations and the line excess Gamma(R), chemical short-range order (Warren-Cowley parameters), solute clusters,
time averages with block-averaging errors and drift tests, the along-line core composition and its spectrum, and one tidy
results table across runs. Across runs (Stage 4, v0.2): uptake exponent, McLean prediction from 0 K site energies,
bottom-up vs top-down hysteresis, reference composition from bulk-control cells; plots in the house style.

## Install
Reproducible conda environment (macOS arm64 and Linux x86_64):

    mamba env create -f environment.yml --channel-priority flexible
    conda activate disegreg
    pip install -e . --no-deps

Only `disegreg/io_ovito.py` needs OVITO (pinned to 3.17.2); everything else needs numpy, scipy and pyyaml.

## Use

    disegreg manifest template > my-run.yaml        # Stage 0: one manifest per run (written by the project)
    disegreg manifest check my-run.yaml
    disegreg stage1 /path/to/*_quenched_step*.dump   # Stage 1: OVITO once per snapshot, cached
    disegreg analyse my-run.yaml --out results/      # Stages 2+3: regions, Gamma(R), order; time averages -> CSV
    disegreg table results/ --out all_D2.csv         # one tidy table across runs
    disegreg across results/ --out stage4.csv        # Stage 4: uptake, McLean (--sites), hysteresis, x_bc error
    disegreg plot gamma results/<run>_gammaR_D2.csv  # plots: gamma | timeseries | delta | uptake | regions | alongline
    disegreg slurm --dumps "<glob>" --out job.slurm  # Stage-1 job array for cmmg (docs/CMMG.md)

Demo manifests: [examples/manifests/](examples/manifests/). Project instances: [docs/HOWTO_PROJECTS.md](docs/HOWTO_PROJECTS.md).
cmmg: [docs/CMMG.md](docs/CMMG.md).

## Tests

    pytest

## Licence and citation
BSD-3-Clause, see [LICENSE](LICENSE). How to cite: [CITATION.cff](CITATION.cff). Archived on Zenodo:
- all versions (concept DOI, resolves to the latest release): [10.5281/zenodo.24513776](https://doi.org/10.5281/zenodo.24513776)
- v0.1.0 (definition set D1 rev. 4): [10.5281/zenodo.24513777](https://doi.org/10.5281/zenodo.24513777)
- v0.1.1 (definition set D1 rev. 4; same numbers as v0.1.0, plus a warning for snapshots after the end of the MC stream):
  [10.5281/zenodo.24535067](https://doi.org/10.5281/zenodo.24535067)

Cite the DOI of the version you used, together with the definition set label (D1 or D2), so a number can be traced to its code
and definitions.
