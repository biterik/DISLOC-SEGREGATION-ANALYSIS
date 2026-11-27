"""Stage 4: across runs (docs/DEFINITIONS.md section 9). Pure numpy; reads the Stage-3 tables of many runs.

    run sets      runs of the same project, potential, host, solute, character, temperature and coexistence line,
                  ordered by the offset Delta. Bulk controls are not part of a set. Within a set, the "up" family is
                  bottom-up + seeded starts, the "down" family top-down starts (code decision, v0.2).
    uptake        n_up = Delta ln c_C / Delta ln x_ref between neighbouring offsets of one family, for every reference
                  x_ref the runs share (c_B: window mean of the in-cell far field; x_iso: manifest, no error; x_bc: the
                  bulk control's window mean and SEM). Errors propagated to first order:
                  SE(n)^2 = [(s1/c1)^2 + (s2/c2)^2 + n^2 ((sx1/x1)^2 + (sx2/x2)^2)] / (ln x2 - ln x1)^2.
                  Independent-site (McLean) segregation gives n_up <= 1; n_up > 1 indicates cooperative uptake.
    hysteresis    bottom-up vs top-down runs at the same Delta: difference of Gamma(20 A) (every reference) and of c_C,
                  error sqrt(SEM_up^2 + SEM_down^2) (section 9); verdict "different" (D2 section 13): |difference| > 2 x error.
    McLean        predicted c_C from the 0 K site energies (disegreg.sites) for every reference, vs the measured c_C
                  (section 9; which sites are core sites: D2 section 12).
    x_bc error    total error of Gamma vs x_bc incl. SEM(x_bc) (D2 section 14).

Output: one tidy long table (STAGE4_COLUMNS), one row per (analysis, set, offset pair or run, quantity, ref). Every row
carries ``status`` (all Stage-4 rows are "definition" since D2; the column stays for the open proposal P5 and the table
still refuses a proposal row with role "headline", docs/PLAN_v0.2.md rule 3). With tables of definition set D1 only what
D1 defines is computed: uptake exponents and the hysteresis differences (no verdict, no McLean, no x_bc error).
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .results import (PAPER_HOLD, STATUS_DEFINITION, STATUS_PROPOSAL, ProposalAsHeadline,  # noqa: F401
                      check_headline_status, read_table)

SET_KEYS = ("project", "potential", "host_element", "solute_element", "character", "temperature", "coex_line")
DELTA_SAME = 1e-3  # two runs are "at the same Delta" if their offsets differ by less than this
STAGE4_COLUMNS = ("analysis", "project", "potential", "host", "solute", "character", "temperature", "family", "delta",
                  "delta_2", "quantity", "region", "basis", "ref", "role", "value", "error", "status", "runs", "flags",
                  "paper_flag", "definition_set", "proposal_set", "tool_version")


# ------------------------------------------------------------------------------------------------------------ runs
@dataclass
class Run:
    file: str
    run_id: str
    manifest: dict
    prov: dict
    rows: dict = field(repr=False)   # (quantity, region, basis, ref) -> Stage-3 row

    @property
    def stage3(self) -> dict:
        return self.prov.get("stage3", {})

    def __getitem__(self, key):
        return self.manifest.get(key)

    @property
    def delta(self) -> float:
        return float(self.manifest.get("delta", np.nan))

    @property
    def is_bulk_control(self) -> bool:
        return self.manifest.get("cell_type", "dislocation") == "bulk_control"

    @property
    def family(self) -> str:
        return "down" if self.manifest.get("start_type") == "top-down" else "up"

    @property
    def paper_flag(self) -> str:
        return self.stage3.get("paper_flag", "") or ""

    @property
    def converged(self) -> bool:
        return bool(self.stage3.get("converged", False))

    def value(self, quantity, region="", basis="", ref=""):
        """(window mean, block-averaging SEM); NaN if the run does not have it."""
        r = self.rows.get((quantity, region, basis, ref))
        return (float(r["mean"]), float(r["sem"])) if r else (float("nan"), float("nan"))

    def refs(self) -> list[str]:
        # a diagnostic-only reference (x_bc of a bulk control that was not accepted) is not a reference for Stage 4
        return sorted({k[3] for k in self.rows if k[0] == "x_ref" and k[3] and not k[3].endswith("_diagnostic")})

    def x_ref(self, ref: str):
        """(x_ref, its error) for one reference, as used by the Stage-3 Gamma of this run."""
        if ref == "c_B":
            return self.value("c", "B")
        if ref == "x_iso":
            x = self.manifest.get("x_iso")
            return (float(x), 0.0) if x is not None else (float("nan"), float("nan"))
        if ref == "x_bc":
            info = self.prov.get("x_bc") or {}
            if info.get("accepted") is False or info.get("enough_solute") is False:
                return float("nan"), float("nan")   # not accepted: diagnostic only (D2 section 15, DEFINITIONS 5 2b)
            if info.get("value") is not None:
                sem = info.get("sem")
                return float(info["value"]), float(sem) if sem is not None else float("nan")
        return self.value("x_ref", ref=ref)[0], float("nan")


def load_runs(stage3_files) -> tuple[list[Run], str]:
    """Stage-3 tables (with their provenance side files) -> runs. Refuses mixed definition sets."""
    runs, labels = [], set()
    for f in stage3_files:
        f = Path(f)
        prov = json.loads(f.with_name(f.name + ".provenance.json").read_text())
        labels.add(prov.get("definition_set"))
        rows = {(r["quantity"], r["region"], r["basis"], r["ref"]): r for r in read_table(f)}
        man = prov.get("manifest_data", {})
        runs.append(Run(file=str(f), run_id=man.get("run_id", f.stem), manifest=man, prov=prov, rows=rows))
    if len(labels) > 1:
        raise ValueError(f"tables from different definition sets {sorted(map(str, labels))} must not be combined")
    return runs, (labels.pop() if labels else "")


def run_sets(runs) -> dict:
    """{set key: runs sorted by Delta}; bulk controls are left out."""
    out = {}
    for r in runs:
        if r.is_bulk_control:
            continue
        out.setdefault(tuple(r[k] for k in SET_KEYS), []).append(r)
    return {k: sorted(v, key=lambda r: (r.delta, r.run_id)) for k, v in sorted(out.items(), key=lambda kv: str(kv[0]))}


def _base(set_key, **kw):
    meta = dict(zip(SET_KEYS, set_key))
    row = {c: "" for c in STAGE4_COLUMNS}
    row.update(project=meta["project"], potential=meta["potential"], host=meta["host_element"],
               solute=meta["solute_element"], character=meta["character"], temperature=meta["temperature"])
    row.update(kw)
    return row


def _flags(*runs) -> str:
    f = [f"not converged: {r.run_id}" for r in runs if not r.converged]
    return "; ".join(f)


def _paper(*runs) -> str:
    return PAPER_HOLD if any(r.paper_flag for r in runs) else ""


# ---------------------------------------------------------------------------------------------------------- uptake
def uptake_exponent(c1, s1, c2, s2, x1, sx1, x2, sx2):
    """n = (ln c2 - ln c1) / (ln x2 - ln x1) and its first-order propagated error. An unknown SEM of c_C (e.g. a run with
    a single snapshot) makes the error unknown (NaN); an unknown x_ref error is taken as 0 (the caller flags it)."""
    if not (c1 > 0 and c2 > 0 and x1 > 0 and x2 > 0) or x1 == x2:
        return float("nan"), float("nan")
    D = np.log(x2) - np.log(x1)
    n = (np.log(c2) - np.log(c1)) / D
    if not (np.isfinite(s1) and np.isfinite(s2)):
        return float(n), float("nan")
    rel = [(s / v) ** 2 for s, v in ((s1, c1), (s2, c2))]
    relx = [(s / v) ** 2 if np.isfinite(s) else 0.0 for s, v in ((sx1, x1), (sx2, x2))]
    return float(n), float(np.sqrt(sum(rel) + n * n * sum(relx)) / abs(D))


def uptake_rows(set_key, runs) -> list[dict]:
    rows = []
    for fam in ("up", "down"):
        series = [r for r in runs if r.family == fam]
        groups = []  # runs of this family grouped by offset
        for r in series:
            if groups and abs(r.delta - groups[-1][0].delta) < DELTA_SAME:
                groups[-1].append(r)
            else:
                groups.append([r])
        # several runs at one offset (e.g. a repeat): every combination is computed and flagged -- no silent choice
        pairs = [(a, b, len(ga) > 1 or len(gb) > 1) for ga, gb in zip(groups[:-1], groups[1:]) for a in ga for b in gb]
        for a, b, several in pairs:
            c1, s1 = a.value("c", "C")
            c2, s2 = b.value("c", "C")
            for ref in sorted(set(a.refs()) & set(b.refs())):
                x1, sx1 = a.x_ref(ref)
                x2, sx2 = b.x_ref(ref)
                n, e = uptake_exponent(c1, s1, c2, s2, x1, sx1, x2, sx2)
                flags = _flags(a, b)
                if not (np.isfinite(s1) and np.isfinite(s2)):
                    flags = "; ".join(x for x in (flags, "c_C error unknown (single snapshot?): error NaN") if x)
                if not np.isfinite(sx1) or not np.isfinite(sx2):
                    flags = "; ".join(x for x in (flags, "x_ref error unknown (taken as 0)") if x)
                if several:
                    flags = "; ".join(x for x in (flags, "several runs at one offset") if x)
                rows.append(_base(set_key, analysis="uptake", family=fam, delta=a.delta, delta_2=b.delta,
                                  quantity="n_up", region="C", ref=ref, value=n, error=e, status=STATUS_DEFINITION,
                                  runs=f"{a.run_id}|{b.run_id}", flags=flags, paper_flag=_paper(a, b)))
    return rows


# ------------------------------------------------------------------------------------------------------ hysteresis
def hysteresis_rows(set_key, runs, params: dict) -> list[dict]:
    """bottom-up vs top-down at the same Delta (seeded runs are not compared: their history is neither)."""
    up = [r for r in runs if r["start_type"] == "bottom-up"]
    down = [r for r in runs if r["start_type"] == "top-down"]
    factor = float(params["factor"]) if params.get("factor") is not None else None   # None: D1, differences only
    rows = []
    for a in up:
        for b in down:
            if abs(a.delta - b.delta) >= DELTA_SAME:
                continue
            for q, region, basis in params["quantities"]:
                refs = sorted(set(a.refs()) & set(b.refs())) if q == "Gamma" else [""]
                for ref in refs:
                    va, sa = a.value(q, region, basis, ref)
                    vb, sb = b.value(q, region, basis, ref)
                    diff, err = va - vb, float(np.hypot(sa, sb))
                    common = dict(analysis="hysteresis", delta=a.delta, delta_2=b.delta, region=region, basis=basis,
                                  ref=ref, runs=f"{a.run_id}|{b.run_id}", flags=_flags(a, b), paper_flag=_paper(a, b))
                    rows.append(_base(set_key, quantity=f"{q}_up_minus_down", value=diff, error=err,
                                      status=STATUS_DEFINITION, **common))
                    if factor is None:
                        continue
                    # 1 = different, 0 = not different, NaN = undecidable (a value or its error missing)
                    verdict = float(abs(diff) > factor * err) if np.isfinite(diff) and np.isfinite(err) else float("nan")
                    rows.append(_base(set_key, quantity=f"{q}_hysteresis", value=verdict,
                                      error=abs(diff) / err if err > 0 else float("nan"),
                                      status=STATUS_DEFINITION, **common))
    return rows


# ---------------------------------------------------------------------------------------------------------- McLean
def mclean_rows(set_key, runs, energy_core, n_missing: float, site_info: dict) -> list[dict]:
    """Predicted vs measured c_C per run and reference; error of the prediction from the x_ref error only.
    Rows: c_C_mclean, c_C_measured_over_mclean, n_up_mclean (local exponent of the prediction, <= 1)."""
    from .sites import mclean_prediction

    rows = []
    for r in runs:
        kT = 8.617333262e-5 * float(r["temperature"])
        c, sc = r.value("c", "C")
        for ref in r.refs():
            x, sx = r.x_ref(ref)
            p = mclean_prediction(energy_core, x, kT, n_missing)
            sp = abs(p.dc_dx) * sx if np.isfinite(sx) else float("nan")
            ratio = c / p.c_C if p.c_C > 0 else float("nan")
            sr = ratio * np.hypot(sc / c if c > 0 else np.nan, sp / p.c_C if np.isfinite(sp) else 0.0)
            common = dict(analysis="mclean", family=r.family, delta=r.delta, region="C", ref=ref, runs=r.run_id,
                          flags="; ".join(x for x in (_flags(r), site_info.get("flags", "")) if x),
                          paper_flag=_paper(r), status=STATUS_DEFINITION)
            rows.append(_base(set_key, quantity="c_C_mclean", value=p.c_C, error=sp, **common))
            rows.append(_base(set_key, quantity="c_C_measured_over_mclean", value=ratio, error=sr, **common))
            rows.append(_base(set_key, quantity="n_up_mclean", value=p.n_up, error=float("nan"), **common))
    return rows


# ---------------------------------------------------------------------------------------------- x_bc error (P3)
def xbc_error_rows(set_key, runs, radii) -> list[dict]:
    """Contribution of SEM(x_bc) to the error of Gamma(R) vs x_bc: N(r < R) / L_x x SEM(x_bc) (D2 section 14).
    The Stage-3 SEM of Gamma vs x_bc holds the time scatter only (x_bc is one number for the whole run)."""
    rows = []
    for r in runs:
        x, sx = r.x_ref("x_bc")
        if not np.isfinite(x):
            continue
        L = r.value("L_x")[0]
        for R in radii:
            n_in = r.value("N_inside", f"{float(R):g}", "r")[0]
            g, sg = r.value("Gamma", f"{float(R):g}", "r", "x_bc")
            e_x = n_in / L * sx if np.isfinite(sx) else float("nan")
            rows.append(_base(set_key, analysis="x_bc_error", family=r.family, delta=r.delta, quantity="Gamma_error_total",
                              region=f"{float(R):g}", basis="r", ref="x_bc", value=float(np.hypot(sg, e_x)),
                              error=e_x, status=STATUS_DEFINITION, runs=r.run_id, flags=_flags(r), paper_flag=_paper(r)))
    return rows


# -------------------------------------------------------------------------------------------------------- driver
def stage4(runs, definitions, site_cores: dict | None = None, gamma_radii=(20.0, 40.0)) -> list[dict]:
    """All Stage-4 rows. ``definitions``: the definition set of the tables (Definitions or its params dict); D2 sections
    hysteresis / bulk_control switch on the verdict and the x_bc error. ``site_cores``: {(project, potential, character):
    (E of core sites, n_missing, info)} (from D2 section 12; the caller loads them only with D2)."""
    p = getattr(definitions, "params", definitions)
    hy = p.get("hysteresis") or {"factor": None,                      # D1: differences of the headline and c_C only
                                 "quantities": [["Gamma", f"{float(gamma_radii[0]):g}", "r"], ["c", "C", ""]]}
    bc = p.get("bulk_control") or {}
    rows = []
    for key, rs in run_sets(runs).items():
        rows += uptake_rows(key, rs)
        rows += hysteresis_rows(key, rs, hy)
        meta = dict(zip(SET_KEYS, key))
        sc = (site_cores or {}).get((meta["project"], meta["potential"], meta["character"]))
        if sc is not None:
            rows += mclean_rows(key, rs, sc[0], sc[1], sc[2])
        if bc.get("propagate_x_ref_error", False):
            rows += xbc_error_rows(key, rs, gamma_radii)
    check_headline_status(rows)
    return rows
