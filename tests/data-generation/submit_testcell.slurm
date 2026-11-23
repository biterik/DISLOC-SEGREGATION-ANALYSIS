#!/bin/bash
# disegreg test cell: 0 K relaxation of a 35 584-atom Ni d60 cell (minutes on 16 cores, < 1 core-h).
# Submit from the run folder (docs/CMMG.md, "Test cell"):  sbatch <repo>/tests/data-generation/submit_testcell.slurm
#SBATCH --job-name=disegreg-testcell
#SBATCH --partition=s.cmmg
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=%x-%j.out

set -euo pipefail
REPO=/cmmc/ptmp/$USER/DISLOC-SEGREGATION-ANALYSIS
POTDIR=/cmmc/ptmp/$USER/POTENTIALS
module purge
module load lammps/250722
export OMP_NUM_THREADS=1

gunzip -c "$REPO/tests/data-generation/testcell_unrelaxed.lmp.gz" > testcell_unrelaxed.lmp
cp "$REPO/tests/data-generation/in.relax_testcell" .
srun lmp -in in.relax_testcell -var POTDIR "$POTDIR" -log log.relax_testcell
gzip -f testcell_relaxed.dump
sha256sum testcell_relaxed.dump.gz
echo "done: $(pwd)/testcell_relaxed.dump.gz"
