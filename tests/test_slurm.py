"""The Stage-1 Slurm array script: valid bash, right array size, chunking logic."""

import subprocess

from disegreg.cli import main


def test_slurm_script(tmp_path):
    dumps = []
    for k in range(45):
        f = tmp_path / "run" / f"x_quenched_step{k * 1000}.dump"
        f.parent.mkdir(exist_ok=True)
        f.write_text("")
        dumps.append(f)
    out = tmp_path / "jobs" / "s1.slurm"
    assert main(["slurm", "--dumps", str(tmp_path / "run" / "*.dump"), "--out", str(out), "--chunk", "20"]) == 0
    text = out.read_text()
    assert "#SBATCH --array=0-2%20" in text and "#SBATCH --partition=s.cmmg" in text
    assert "--definitions D1" in text and "DISEGREG_CACHE=/u/$USER/PTMP/DISEGREG-CACHE" in text
    lines = (tmp_path / "jobs" / "s1.slurm.dumps.txt").read_text().split()
    assert len(lines) == 45
    assert subprocess.run(["bash", "-n", str(out)]).returncode == 0
    # the chunk selection of the script, run for task 2: the last 5 dumps
    body = text.split("CHUNK=20\n", 1)[1].split("echo", 1)[0]  # the selection lines exactly as in the script
    snippet = "SLURM_ARRAY_TASK_ID=2; CHUNK=20\n" + body + 'echo ${#DUMPS[@]} ${DUMPS[0]}'
    r = subprocess.run(["bash", "-c", snippet], capture_output=True, text=True)
    n, first = r.stdout.split()
    assert n == "5" and first == lines[40]
