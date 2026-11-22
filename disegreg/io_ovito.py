"""The ONLY module of disegreg that imports OVITO: Stage 1 of one snapshot (docs/DESIGN.md).

Reads a LAMMPS dump (plain or .gz, triclinic) once, then runs three independent pipelines on the same in-memory data:
a-CNA (region labels), PTM with the fcc/hcp/bcc/icosahedral templates (extra label), and DXA (dislocation lines).
Independent pipelines because every structure modifier writes the same "Structure Type" property: in a PTM -> DXA
pipeline the labels read afterwards are DXA's internal ones (verified 2026-11-22 on Ni-Cu d60; DXA lines are identical
with or without a preceding PTM / a-CNA, i.e. DXA does its own structure identification).
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import ovito
from ovito.io import import_file
from ovito.modifiers import CommonNeighborAnalysisModifier as _CNA
from ovito.modifiers import DislocationAnalysisModifier as _DXA
from ovito.modifiers import PolyhedralTemplateMatchingModifier as _PTM
from ovito.pipeline import Pipeline, StaticSource

from .cache import STRUCTURE_CODES, Snapshot, stage1_params
from .frame import burgers_family
from .provenance import environment, sha256_file

# the cache stores OVITO's integer codes; make sure they mean what cache.STRUCTURE_CODES says
for _name, _code in STRUCTURE_CODES.items():
    _key = {"other": "OTHER", "fcc": "FCC", "hcp": "HCP", "bcc": "BCC", "ico": "ICO"}[_name]
    assert int(getattr(_CNA.Type, _key)) == _code and int(getattr(_PTM.Type, _key)) == _code, _name

_PTM_TEMPLATES = {"fcc": _PTM.Type.FCC, "hcp": _PTM.Type.HCP, "bcc": _PTM.Type.BCC, "ico": _PTM.Type.ICO,
                  "sc": _PTM.Type.SC, "cubic_diamond": _PTM.Type.CUBIC_DIAMOND,
                  "hex_diamond": _PTM.Type.HEX_DIAMOND, "graphene": _PTM.Type.GRAPHENE}

def ovito_version() -> str:
    return ovito.version_string


def make_modifiers(structure: dict):
    """a-CNA, PTM and DXA modifiers configured from the 'structure' section of a definition set."""
    acna = _CNA(mode=_CNA.Mode.AdaptiveCutoff)
    p = structure["ptm_extra_label"]
    ptm = _PTM(rmsd_cutoff=float(p["rmsd_cutoff"]))
    wanted = {int(_PTM_TEMPLATES[t]) for t in p["templates"]}
    for st in ptm.structures:
        if st.id != int(_PTM.Type.OTHER):
            st.enabled = st.id in wanted
    d = structure["dxa"]
    if d["input_structure"] != "fcc":
        raise ValueError("only fcc input structure is supported")
    dxa = _DXA(input_crystal_structure=_DXA.Lattice.FCC,
               trial_circuit_length=int(d["trial_circuit_length"]),
               circuit_stretchability=int(d["circuit_stretchability"]),
               line_smoothing_enabled=int(d["line_smoothing_level"]) > 0,
               line_smoothing_level=int(d["line_smoothing_level"]),
               line_coarsening_enabled=bool(d["line_coarsening"]),
               defect_mesh_smoothing_level=int(d["defect_mesh_smoothing_level"]),
               only_perfect_dislocations=bool(d["only_perfect_dislocations"]))
    return acna, ptm, dxa


def _run(data, modifier):
    pl = Pipeline(source=StaticSource(data=data))
    pl.modifiers.append(modifier)
    return pl.compute()


def compute_stage1(dump, definitions, source_sha256: str | None = None) -> Snapshot:
    """Read one dump and compute everything Stage 1 caches."""
    dump = Path(dump).resolve()
    t0 = time.perf_counter()
    data = import_file(str(dump)).compute()
    t_read = time.perf_counter() - t0
    parts = data.particles
    ids = np.asarray(parts.identifiers, dtype=np.int64) if parts.identifiers is not None else np.arange(1, parts.count + 1)
    order = np.argsort(ids, kind="stable")
    ids = ids[order]
    types = np.asarray(parts.particle_types, dtype=np.int64)[order]
    if types.max() > 127:
        raise ValueError("atom type > 127 does not fit the int8 cache column")

    structure = definitions["structure"]
    acna_mod, ptm_mod, dxa_mod = make_modifiers(structure)
    timings = {"read": t_read}
    t = time.perf_counter()
    acna = np.asarray(_run(data, acna_mod).particles["Structure Type"], dtype=np.int8)[order]
    timings["acna"] = time.perf_counter() - t
    t = time.perf_counter()
    ptm = np.asarray(_run(data, ptm_mod).particles["Structure Type"], dtype=np.int8)[order]
    timings["ptm"] = time.perf_counter() - t
    t = time.perf_counter()
    dd = _run(data, dxa_mod)
    timings["dxa"] = time.perf_counter() - t

    segs = list(dd.dislocations.segments) if dd.dislocations is not None else []
    pts = [np.asarray(s.points, dtype=float).reshape(-1, 3) for s in segs]
    offsets = np.concatenate([[0], np.cumsum([len(p) for p in pts])]).astype(np.int64)
    true_b = np.array([s.true_burgers_vector for s in segs], dtype=float).reshape(-1, 3)
    families = np.array([burgers_family(b) for b in true_b], dtype="U10")

    meta = {
        "source": str(dump),
        "source_sha256": source_sha256 or sha256_file(dump),
        "source_bytes": dump.stat().st_size,
        "timestep": int(data.attributes.get("Timestep", -1)),
        "n_atoms": int(parts.count),
        "stage1_params": stage1_params(definitions),
        "definitions_sha256": definitions.sha256,
        "provenance": {**environment(), "ovito": ovito_version()},
        "timings_s": {k: round(v, 3) for k, v in timings.items()},
        "dxa_attributes": {k: float(v) for k, v in dd.attributes.items() if k.startswith("DislocationAnalysis.")},
    }
    return Snapshot(
        positions=np.asarray(parts.positions, dtype=np.float32)[order],
        types=types.astype(np.int8),
        cell_matrix=np.asarray(data.cell[...], dtype=float),
        pbc=np.array(data.cell.pbc, dtype=bool),
        acna=acna,
        ptm=ptm,
        seg_points=np.vstack(pts) if pts else np.empty((0, 3)),
        seg_offsets=offsets,
        seg_true_b=true_b,
        seg_spatial_b=np.array([s.spatial_burgers_vector for s in segs], dtype=float).reshape(-1, 3),
        seg_length=np.array([s.length for s in segs], dtype=float),
        seg_family=families,
        seg_cluster=np.array([s.cluster_id for s in segs], dtype=np.int32),
        seg_infinite=np.array([s.is_infinite_line for s in segs], dtype=bool),
        seg_loop=np.array([s.is_loop for s in segs], dtype=bool),
        ids=None if np.array_equal(ids, np.arange(1, len(ids) + 1)) else ids,
        meta=meta,
    )
