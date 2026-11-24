# disegreg

Shared analysis of solute segregation at dislocations in atomistic MD/MC simulations (molecular dynamics / Monte
Carlo), written for the paper *"Dislocation Defect Phases and how to find them"*. One code and one versioned
definition set for all simulation projects, so that every number is computed the same way.

**Status: v0.1.0 (24 Nov 2026), Stages 0-3.** Definitions: [docs/DEFINITIONS.md](docs/DEFINITIONS.md)
(set D1, rev. 4, final). Design: [docs/DESIGN.md](docs/DESIGN.md). Changes: [CHANGELOG.md](CHANGELOG.md).

## What it computes
From quenched LAMMPS snapshots with one dissociated dislocation: the distance of every atom to the dislocation line
(OVITO's dislocation extraction algorithm, DXA, with exact periodic images of triclinic cells, plus a structural fallback
when DXA finds only part of the line), exclusive regions (core, stacking fault, shell, transition, far field), solute
concentrations and the line excess Gamma(R), chemical short-range order (Warren-Cowley parameters), solute clusters,
time averages with block-averaging errors and drift tests, and comparisons across runs.

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
    disegreg table results/ --out all_D1.csv         # one tidy table across runs
    disegreg slurm --dumps "<glob>" --out job.slurm  # Stage-1 job array for cmmg (docs/CMMG.md)

Demo manifests: [examples/manifests/](examples/manifests/). Project instances: [docs/HOWTO_PROJECTS.md](docs/HOWTO_PROJECTS.md).
cmmg: [docs/CMMG.md](docs/CMMG.md).

## Tests

    pytest

## Licence and citation
BSD-3-Clause, see [LICENSE](LICENSE). How to cite: [CITATION.cff](CITATION.cff).
