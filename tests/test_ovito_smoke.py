"""OVITO environment check (runs in CI on Linux and on the M5; skipped when OVITO is not installed)."""

import numpy as np
import pytest

ovito = pytest.importorskip("ovito")
pytestmark = pytest.mark.ovito

from helpers import fcc_in_cell  # noqa: E402

from disegreg.geometry import Cell  # noqa: E402


def test_ovito_version_pinned():
    assert ovito.version_string == "3.17.2"


def test_acna_and_dxa_on_perfect_fcc():
    from ovito.data import DataCollection, Particles, SimulationCell
    from ovito.modifiers import CommonNeighborAnalysisModifier, DislocationAnalysisModifier
    from ovito.pipeline import Pipeline, StaticSource

    a = 3.52
    cell = Cell.lammps(8 * a, 8 * a, 8 * a, pbc=(True, True, True))
    pos = fcc_in_cell(cell, a)
    dc = DataCollection()
    dc.objects.append(SimulationCell(matrix=np.column_stack([cell.H, cell.origin]), pbc=cell.pbc))
    p = Particles()
    p.create_property("Position", data=pos)
    dc.objects.append(p)
    pl = Pipeline(source=StaticSource(data=dc))
    pl.modifiers.append(CommonNeighborAnalysisModifier())
    pl.modifiers.append(DislocationAnalysisModifier(input_crystal_structure=DislocationAnalysisModifier.Lattice.FCC))
    d = pl.compute()
    assert d.attributes["CommonNeighborAnalysis.counts.FCC"] == len(pos)
    assert len(d.dislocations.segments) == 0
