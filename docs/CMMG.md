# disegreg on cmmg -- step by step (Erik runs every command; nothing here is run by an assistant)

Status: parts A and C done on 24 Nov 2026 (micromamba 2.5.0; `disegreg version` OK, OVITO 3.17.2, pytest 98 passed / 2 skipped;
test cell relaxed in 22 s, sha256 c7ab1858...f413). One Stage-1 job (part B) ran: ~19 s per 587k-atom snapshot on 4 cores.

Machine: cmmg (x86_64 Linux, SLES 15 SP7, glibc 2.38). App: Terminal on the M5, logged in to a cmmg login node
(`ssh <your user>@<cmmg login node>`; the real names are in the private `canon/local`).
Everything lives under PTMP (`/cmmc/ptmp/$USER`, the same place as `/u/$USER/PTMP`); nothing goes into `$HOME`.

## A. One-time setup (about 15 minutes, login node; no simulation runs here)

Copy each block into the cmmg terminal as it is.

```bash
# A1 -- folders and cache redirects (so nothing is written into $HOME)
export P=/cmmc/ptmp/$USER
mkdir -p $P/tmp $P/.conda-pkgs $P/.pip-cache $P/conda-envs $P/DISEGREG-JOBS
export TMPDIR=$P/tmp CONDA_PKGS_DIRS=$P/.conda-pkgs PIP_CACHE_DIR=$P/.pip-cache

# A2 -- the disegreg code (public GitHub repository)
cd $P && git clone https://github.com/biterik/DISLOC-SEGREGATION-ANALYSIS.git

# A3 -- micromamba (no module needed): the user's own micromamba, set up in June 2026 and loaded by the cluster .bashrc
#       (dotfiles repo); the other project envs (grace, HDNNP, ovito_env) were made with it. Check it is there:
/cmmc/u/$USER/.local/bin/micromamba --version

# A4 -- the environment, exactly the packages solved for cmmg (env-lock/disegreg-linux-64.txt, incl. OVITO 3.17.2)
/cmmc/u/$USER/.local/bin/micromamba create -y -p $P/conda-envs/disegreg -f $P/DISLOC-SEGREGATION-ANALYSIS/env-lock/disegreg-linux-64.txt

# A5 -- disegreg itself into that environment (editable: a later `git pull` updates it)
$P/conda-envs/disegreg/bin/pip install --no-cache-dir --no-deps -e $P/DISLOC-SEGREGATION-ANALYSIS

# A6 -- check
$P/conda-envs/disegreg/bin/disegreg version
$P/conda-envs/disegreg/bin/python -c "import ovito; print('OVITO', ovito.version_string)"
cd $P/DISLOC-SEGREGATION-ANALYSIS && $P/conda-envs/disegreg/bin/python -m pytest -q
```

What you should see:
- A3 prints a version (2.5.0 on 24 Nov 2026).
- A4 downloads ~295 packages (4 of them from conda.ovito.org) and ends without an error.
- A6: `disegreg 0.1.0  definition set D1 (...)`, then `OVITO 3.17.2`, then pytest ends with a line like
  `NN passed, 2 skipped in 20s` (the 2 skipped tests need the relaxed test cell, part C).

If A4 fails (e.g. "package not found"), use the general file instead and save the result as the new lock file:
```bash
/cmmc/u/$USER/.local/bin/micromamba create -y -p $P/conda-envs/disegreg -f $P/DISLOC-SEGREGATION-ANALYSIS/environment.yml
/cmmc/u/$USER/.local/bin/micromamba env export -p $P/conda-envs/disegreg --explicit > $P/disegreg-linux-64.txt
```
and send me the file `$P/disegreg-linux-64.txt` (part D shows how to copy it to the M5).

## B. Stage 1 for many runs (job array)

```bash
export P=/cmmc/ptmp/$USER
$P/conda-envs/disegreg/bin/disegreg slurm \
    --dumps "$P/<project>/<thread>/*/*_quenched_step*.dump" \
    --out $P/DISEGREG-JOBS/<project>_<thread>_stage1.slurm
```
It prints the number of dumps, array tasks and the estimated core-hours, and the `sbatch` line to use, e.g.
`sbatch /cmmc/ptmp/<user>/DISEGREG-JOBS/<project>_<thread>_stage1.slurm`. Above 20 000 core-h it prints a WARNING:
then do not submit before we have talked. Cached snapshots are skipped, so resubmitting the same script only computes
what is missing. Logs: `$P/DISEGREG-JOBS/slurm-logs/`. Cache: `/u/$USER/PTMP/DISEGREG-CACHE` (rebuilt from the dumps
if PTMP is ever purged).

Updating the code later: `cd /cmmc/ptmp/$USER/DISLOC-SEGREGATION-ANALYSIS && git pull`.

## C. Test cell (once; < 1 core-h)

A 35 584-atom Ni cell with one 60-degree dislocation and 2 % Cu (`tests/data-generation/make_test_cell.py`) is
relaxed at 0 K with the Fischer 2018 EAM potential (`Cu_Ni_Fischer_2018.eam.alloy` in `/cmmc/ptmp/$USER/POTENTIALS`),
so that the OVITO tests have a real dissociated dislocation.

```bash
mkdir -p /cmmc/ptmp/$USER/DISEGREG-TESTCELL && cd /cmmc/ptmp/$USER/DISEGREG-TESTCELL
sbatch /cmmc/ptmp/$USER/DISLOC-SEGREGATION-ANALYSIS/tests/data-generation/submit_testcell.slurm
```
After a few minutes (22 s wall on 24 Nov) `disegreg-testcell-<jobid>.out` ends with `RUN marker: ALL PHASES COMPLETE`, a sha256 line and
`done: /cmmc/ptmp/<user>/DISEGREG-TESTCELL/testcell_relaxed.dump.gz`.

## D. Copying a result to the M5

On the M5, in Terminal (not logged in to cmmg):
```bash
rsync -av <your user>@<cmmg login node>:/cmmc/ptmp/<your user>/DISEGREG-TESTCELL/testcell_relaxed.dump.gz /Users/<you>/DEVEL/DISLOC-SEGREGATION-ANALYSIS/tests/data/
```
