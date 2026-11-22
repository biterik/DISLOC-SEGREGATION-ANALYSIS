"""Crystal frame and Burgers vectors (pure numpy). v0.1.0: Burgers-vector families; rotation of b and xi into the
crystal frame from the cell-builder metadata follows with the manifest (Stage 0)."""

from __future__ import annotations

import numpy as np

# Burgers-vector families: sorted absolute components of DXA's true Burgers vector (units of the lattice parameter)
FAMILIES = {
    "1/6<112>": (1 / 6, 1 / 6, 1 / 3),  # Shockley partial
    "1/2<110>": (0.0, 0.5, 0.5),        # perfect
    "1/6<110>": (0.0, 1 / 6, 1 / 6),    # stair-rod
    "1/3<100>": (0.0, 0.0, 1 / 3),      # Hirth
    "1/3<111>": (1 / 3, 1 / 3, 1 / 3),  # Frank
}
SHOCKLEY = "1/6<112>"
PERFECT = "1/2<110>"


def burgers_family(b_true, tol: float = 1e-3) -> str:
    """Family name of a true Burgers vector (any sign / permutation), "other" if none matches."""
    s = np.sort(np.abs(np.asarray(b_true, dtype=float)))
    for name, ref in FAMILIES.items():
        if np.allclose(s, ref, atol=tol):
            return name
    return "other"
