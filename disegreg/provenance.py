"""Provenance: tool version (incl. git hash), library versions, host -- recorded in every output (DEFINITIONS 10)."""

from __future__ import annotations

import hashlib
import platform
import socket
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from . import __version__

_REPO = Path(__file__).resolve().parent.parent


def git_info() -> dict:
    """Git hash and dirty flag of the disegreg source tree, if it is a git checkout (editable install)."""
    try:
        h = subprocess.run(["git", "-C", str(_REPO), "rev-parse", "HEAD"], capture_output=True, text=True, timeout=5)
        if h.returncode != 0:
            return {}
        st = subprocess.run(["git", "-C", str(_REPO), "status", "--porcelain", "--untracked-files=no"],
                            capture_output=True, text=True, timeout=5)
        return {"git_hash": h.stdout.strip(), "git_dirty": bool(st.stdout.strip())}
    except (OSError, subprocess.SubprocessError):
        return {}


def versions() -> dict:
    import numpy
    import scipy

    out = {"disegreg": __version__, **git_info(), "python": platform.python_version(), "numpy": numpy.__version__,
           "scipy": scipy.__version__}
    if "ovito" in sys.modules:
        out["ovito"] = sys.modules["ovito"].version_string
    return out


def environment() -> dict:
    return {**versions(), "host": socket.gethostname(), "platform": platform.platform(),
            "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds")}


def sha256_file(path, chunk: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while block := f.read(chunk):
            h.update(block)
    return h.hexdigest()
