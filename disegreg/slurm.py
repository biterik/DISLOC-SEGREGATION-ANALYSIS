"""Slurm job array for Stage 1 on cmmg (docs/CMMG.md). disegreg only WRITES the script; Erik submits it.

One array task = one chunk of dumps, each read once by OVITO (a-CNA, PTM, DXA) and written to the Stage-1 cache.
Cached snapshots are skipped automatically (cache hit), so a resubmitted array only computes what is missing.
"""

from __future__ import annotations

import shlex
from pathlib import Path

SECONDS_PER_SNAPSHOT = 10.0     # estimate for a ~590k-atom dump on one cmmg node share (M5: 4.6 s); refine after the first run
CORE_H_LIMIT = 20000.0          # above this Erik's explicit OK is needed (CLAUDE.md)

TEMPLATE = """#!/bin/bash
# disegreg Stage 1 job array -- written by `disegreg slurm` ({version}, definition set {definitions})
# {n_dumps} dumps in {n_tasks} array tasks of <= {chunk}; list: {listfile}
# Estimated cost: {core_h:.1f} core-h ({sec:.0f} s per snapshot x {cpus} cores; refine after the first run)
#SBATCH --job-name={job_name}
#SBATCH --partition={partition}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}
#SBATCH --time={time}
#SBATCH --array=0-{last_task}%{max_parallel}
#SBATCH --output={logdir}/%x-%A_%a.out

set -euo pipefail
module purge                                  # nothing from modules is needed: the conda env is self-contained
export OVITO_THREAD_COUNT=$SLURM_CPUS_PER_TASK OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export DISEGREG_CACHE={cache}
export TMPDIR={tmpdir}
mkdir -p "$DISEGREG_CACHE" "$TMPDIR"

CHUNK={chunk}
START=$(( SLURM_ARRAY_TASK_ID * CHUNK ))
DUMPS=()
while IFS= read -r line; do DUMPS+=("$line"); done < <(sed -n "$(( START + 1 )),$(( START + CHUNK ))p" {listfile})
echo "task $SLURM_ARRAY_TASK_ID: ${{#DUMPS[@]}} dumps from line $START, host $(hostname), $(date -Is)"
{env}/bin/disegreg stage1 --definitions {definitions} "${{DUMPS[@]}}"
echo "done $(date -Is)"
"""


def stage1_array(dumps, out, env_prefix: str, cache: str, chunk: int = 20, partition: str = "s.cmmg", cpus: int = 4,
                 mem: str = "16G", time: str = "02:00:00", max_parallel: int = 20, definitions: str = "D1",
                 job_name: str = "disegreg-s1", tmpdir: str = "/cmmc/ptmp/$USER/tmp", version: str = "?"):
    """Write <out> (sbatch script) and <out>.dumps.txt (one dump per line). Returns (script, list, n_tasks, core_h)."""
    dumps = [str(Path(d)) for d in dumps]
    if not dumps:
        raise ValueError("no dumps")
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    listfile = out.with_name(out.name + ".dumps.txt")
    listfile.write_text("\n".join(dumps) + "\n")
    n_tasks = (len(dumps) + chunk - 1) // chunk
    core_h = len(dumps) * SECONDS_PER_SNAPSHOT * cpus / 3600.0
    logdir = out.parent / "slurm-logs"
    logdir.mkdir(exist_ok=True)
    out.write_text(TEMPLATE.format(version=version, definitions=definitions, n_dumps=len(dumps), n_tasks=n_tasks,
                                   chunk=chunk, listfile=shlex.quote(str(listfile.resolve())), core_h=core_h,
                                   sec=SECONDS_PER_SNAPSHOT, cpus=cpus, job_name=job_name, partition=partition, mem=mem,
                                   time=time, last_task=n_tasks - 1, max_parallel=max_parallel,
                                   logdir=str(logdir.resolve()), cache=cache, tmpdir=tmpdir, env=env_prefix))
    out.chmod(0o755)
    return out, listfile, n_tasks, core_h
