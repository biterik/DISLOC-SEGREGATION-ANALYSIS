# CLAUDE.md -- disegreg (DISLOC-SEGREGATION-ANALYSIS)

You are the ANALYSIS-TOOL instance for Erik Bitzek's paper "Dislocation Defect Phases and how to find them". You own and develop ONE shared
analysis package, `disegreg`, used read-only by all simulation projects. Read, in this order, before doing anything:
1. docs/DESIGN.md            -- agreed design, data facts, where things run, starting material, testing, releases
2. docs/DEFINITIONS.md       -- definition set D1 (draft rev. 2), still under Erik's review
3. docs/FEEDBACK_from-TURLO-REDO_2026-11-22.md -- feedback already folded into DEFINITIONS rev. 2
4. ~/Desktop/PUBLISH/2026_Dislocation-Defect-Phases/CLAUDE-NOTES/PAPER-PLAN_2026-11-22.md
5. ~/Desktop/PUBLISH/2026_Dislocation-Defect-Phases/CLAUDE-NOTES/INSTRUCTIONS_ANALYSIS-TOOL_2026-11-22.md

## Hard rules
- NEVER modify anything under ~/Desktop/SIMULATIONS/ (project folders are read-only input), except the cache folder
  ~/Desktop/SIMULATIONS/DISEGREG-CACHE. Do not modify the paper folder except when Erik asks.
- Working mode: WORKSHOP with Erik. Propose and ask before deciding anything that changes a definition; every definition change gets a new
  label (D1 -> D2) and a line in the revision log of DEFINITIONS.md. Code decisions inside the agreed design you may take, and report.
- Computations only on the M5 (natively) or on cmmg -- never in a cloud workspace. Expensive runs (> ~20k core-h) only with Erik's explicit OK;
  you prepare, Erik submits on cmmg.
- When Erik must do something himself, give a complete step-by-step recipe: which machine (M5 / cmmg), which app (Terminal, browser),
  the absolute directory, the exact command, and what the output should look like. Always absolute paths. He does not know git/shell jargon.
- Spell out every abbreviation and explain every field-specific term on first use (outside atomistic metals/fracture assume no prior knowledge).
- Plots: fonts >= 13 pt at slide size, colour-blind-safe palettes, SYMMETRIC colour scales for every map of a signed quantity.
- Reports as SusMat-template pptx decks go to ~/Desktop/PUBLISH/2026_Dislocation-Defect-Phases/SIM-REPORTS/ (dated file name).
- Commits: small and descriptive; CHANGELOG.md kept up to date; tagged releases v0.1.0, ...; CITATION.cff funding text = SFB 1394 A02 only.

## Decisions taken by Erik (22 Nov 2026)
- Name: repository DISLOC-SEGREGATION-ANALYSIS (GitHub biterik/DISLOC-SEGREGATION-ANALYSIS, PUBLIC); package and CLI `disegreg`.
- Licence BSD-3-Clause (LICENSE written). CITATION.cff written (funding: SFB 1394 A02 only).
- Far field: r >= 60 A from every line incl. periodic images, a-CNA fcc, >= 15 A from z faces. The 40 A glide-plane slab only as an
  optional extra column, not default.
- Region labels: a-CNA. DXA uses its own internal structure identification (verify on one real snapshot).
- Offset Delta = (dmu - dmu_coex)/kT, negative = below coexistence; dmu = mu_solute - mu_host.
- Gamma(R) per unit cell length L_x along the line (both partials together), R from the nearest line; fault excess separately.
- Quenched snapshots only.
- Cache: cmmg /u/$USER/PTMP/DISEGREG-CACHE, M5 ~/Desktop/SIMULATIONS/DISEGREG-CACHE.
- Manifests owned by the project instances in their folders; the demo manifests live in examples/manifests/ of this repo.
- Testing: synthetic tests anywhere, OVITO tests in GitHub Actions and on the M5, regression on M5/cmmg run by Erik with your recipes.
- Scope: v0.1.0 and v0.2 as in docs/DESIGN.md.

## Open (waiting for Erik)
- Nothing for v0.1.0. All definition decisions were taken on 24 Nov (docs/DEFINITIONS.md rev. 4, revision log).

## State (22 Nov 2026, evening)
- Done: steps 1 and 2. M5 conda env `disegreg` (Python 3.12, OVITO 3.17.2 from conda.ovito.org, created with
  `mamba env create -f environment.yml --channel-priority flexible`; a-CNA + DXA smoke test passed). Git repo on GitHub
  biterik/DISLOC-SEGREGATION-ANALYSIS (public), branch main.
- cmmg facts: SLES 15 SP7, glibc 2.38 (OVITO 3.17.2 wheels/conda builds OK). Existing cmmg env /cmmc/ptmp/$USER/conda-envs/ovito_env has
  OVITO 3.14.1 -- leave it alone; disegreg gets its own env under PTMP from environment.yml (recipe still to write), plus
  env-lock/disegreg-linux-64.txt.
- OVITO pin 3.17.2 agreed; if regression shows DXA differences vs the old scripts, rerun one snapshot with the old OVITO version.

- Step 3 done: package skeleton, definitions/D1.yaml (inside the package: disegreg/definitions/), geometry.py + tests, CI.
- D1 rev. 3 (Erik accepted both proposals): per-cluster fallback axes, used slab-wise only where DXA misses. min_slab_fraction 0.5 and
  active_below_dxa_coverage 0.5 confirmed by Erik 23 Nov. Open [CONFIRM]s: region precedence F before O; partial_merge_distance 2 A.
- Step 4 (first half) done: io_ovito.py, cache.py, `disegreg stage1`; DXA-internal-labels question VERIFIED on Ni-Cu d60 (see
  DEFINITIONS section 2). First real cache entry: ~/Desktop/SIMULATIONS/DISEGREG-CACHE/Cu-AT-DISLOCS-Ni-EAM-MEAM-ACE/... (5.6 MB, 4.6 s).

- Step 4b done (23 Nov): regions, composition, order, manifest (schema 1), stage2 driver, results CSV, `disegreg manifest|analyse`;
  83 tests. Ni-Cu d60 demo (examples/manifests/): z_gp builder 125.51 vs lines 125.47 A; partial separation 18.8 +- 1.0 A;
  Stage 2 1.1 s/snapshot; c_B / x_iso = 0.93 (so Gamma(R) falls beyond ~20 A with x_ref = x_iso -- point for the regression report).

- 23 Nov (Erik: "please do 1-7"): Stage 3 (timeseries.py), crystal frame, `disegreg table`, `disegreg slurm`, linux-64 lock
  (solved on the M5), docs/CMMG.md, test cell (tests/data-generation/, waits for Erik's cmmg run), regression (regression/,
  deck SIM-REPORTS/2026-11-23_disegreg-D1-regression.pptx), demo manifests Ni-Cu / Nb-Au / TURLO, docs/HOWTO_PROJECTS.md.
  Regression: old recipe reproduces Ni-Cu 16/19 exactly (3 ACE runs: one partial 0.1-0.2 A off; NOT the OVITO version --
  3.14.1 = 3.17.2 on the M5; check one snapshot on cmmg), Nb-Au within 0.4 %; TURLO R3 mixed mode, axis recovers 10 %.
  Open for Erik: drift criterion proposal (DEFINITIONS 8), x_ref choice (Ni-Cu c_B/x_iso = 0.93), headline Gamma(20) vs
  Gamma(40), [CONFIRM]s (F before O, merge distance 2 A), DEFINITIONS rev. 3 review before tagging v0.1.0.

- 24 Nov (in the claude.ai design session, with Erik):
  - cmmg parts A and C done: env via Erik's own micromamba (/cmmc/u/$USER/.local/bin/micromamba, 2.5.0; NO conda module needed;
    the newest cluster conda module would be anaconda/3/2023.03, not used), OVITO 3.17.2, pytest on cmmg 98 passed / 2 skipped.
    docs/CMMG.md rewritten to the micromamba route (uncommitted).
  - Test cell relaxed on cmmg (22 s); copied to tests/data/testcell_relaxed.dump.gz (sha256 c7ab1858...f413, gitignored).
    On the M5: 99 passed, 1 FAILED: tests/test_testcell.py expects b_total_label "1/2[1 -1 0]", DXA gives "1/2[-1 1 0]".
    Cause: the cell builders' nominal b is the NEGATIVE of the DXA b with xi = +X (OVITO: LH/SF, same b as RH/FS) -- the
    production d60 cells show the same relation (nominal 1/2[-1 1 0] -> DXA 1/2[1 -1 0]). Not an analysis error.
  - Ni-Cu ACE d90 -0.10 kT snapshot (step 790000) through Stage 1 on cmmg (job array via `disegreg slurm`; ~19 s on 4 cores),
    compared with the M5 cache entry (same key fc0d2769): positions, a-CNA, PTM bit-identical; DXA: same segments and b, one partial
    identical, the other off by 0.09 A mean / 0.5 A max (mean y 0.22 A, length 0.13 A). So the 3 ACE regression mismatches are DXA's
    cross-machine non-reproducibility, not disegreg and not the OVITO version.
  - Erik's decisions (now in DEFINITIONS rev. 4): drift test |change| < 2 SE(change); x_ref: report vs both c_B and x_iso until a bulk
    control exists, no ACE Ni-Cu enrichment in the paper before it; headline Gamma(20 A), Gamma(R) 5-60 A in the supplement, Gamma(40 A)
    as check where the fallback axis is active; F before O and partial_merge_distance 2 A confirmed.
  - Files changed by that session: docs/DEFINITIONS.md, docs/CMMG.md, CLAUDE.md (committed 24 Nov as 8063d67).

- 24 Nov (Claude Code session on the M5): next steps 1-4 of the 24 Nov list done, committed locally, NOT pushed, NOT tagged.
  Code matches DEFINITIONS rev. 4 (drift criterion change_se default + column flat_rev3_sem; column `ref` with c_B and x_iso;
  paper_flag for ACE Ni-Cu without bulk control; role headline/check/curve, check_required, <run>_gammaR_D1.csv). 103 tests pass
  on the M5 (incl. the test cell). DXA deterministic on one machine (4 Stage-1 runs of ACE d90 790000 bit-identical).
  regression/check_tolerance.py: Ni-Cu 16 exact + 3 ACE within tolerance; Nb-Au 24 within tolerance; closest calls ~90 % of the
  tolerance (Ni-Cu d60 ACE -0.50 core c; Nb-Au EAM-X d0 -0.25 separation).

- 24 Nov: pushed; GitHub Actions green incl. the test cell on Linux x86_64. Released v0.1.0 (tag, GitHub release; Zenodo webhook
  active on the repo, so Zenodo mints the DOI from the release).
- Zenodo DOIs: concept 10.5281/zenodo.24513776 (all versions), v0.1.0 10.5281/zenodo.24513777 (cite this one in the paper,
  with D1); both in README.md and CITATION.cff.

## State v0.2 (branch v0.2-dev, worktree ~/DEVEL/DISLOC-SEGREGATION-ANALYSIS-v0.2, env disegreg-dev; unattended session)
- Rules: docs/PLAN_v0.2.md. Decisions for Erik: docs/OPEN-DECISIONS_v0.2.md (the one page). D1 unchanged (D1.yaml byte-identical);
  D2 (25 Nov) = D1 + sections 12-16 in docs/DEFINITIONS.md and disegreg/definitions/D2.yaml; the one open proposal (P5) in
  "[PROPOSAL v0.2]" + proposals_v0.2.yaml, its outputs carry status = proposal.
- 24 Nov: item A (Stage 4: across.py, sites.py, bulk-control manifests + automatic x_bc) done with tests; plotting module
  (item B) drafted. disegreg-dev now has the editable install of the worktree (it had none).
- 24 Nov: items C (alongline.py, P4) and D (domains.py, P5) done with synthetic tests; item B plots checked by eye on the Nb-Au
  demo (titles wrapped, legend rows, uptake y range from well-determined values, regions as dots on log axis).
- Demo outputs: Stage 2/3 (not in git) regression/output/v0.2-final/ (59 runs); summaries + figures (in git)
  regression/results/v0.2-demo/ (regression/v02_demo_summary.py). Findings in docs/OPEN-DECISIONS_v0.2.md section 3:
  Lee-Shim n_up 5.6 (60 deg) / 4.8 (edge) at -0.25 -> -0.10 kT; TURLO hysteresis different only at +2 kT; ACE core = particles.
- 25 Nov (Erik's decisions, OPEN-DECISIONS section 5; Erik authorised work on main for v0.1.1 only): v0.1.1 released from
  main (Stage-3 WARNING + count, two robustness fixes, HOWTO rsync line; D1 unchanged, outputs cell-by-cell identical to v0.1.0;
  Zenodo DOI v0.1.1 10.5281/zenodo.24535067). On the branch: definition set D2 = D1 rev. 4 unchanged + DEFINITIONS sections
  12-16 (P1-P4, E1); D2.yaml, default D2, P5 stays in proposals_v0.2.yaml. v0.2-dev NOT merged with main, v0.2 NOT tagged.

## Next steps
1. Hand each project instance its "how to call disegreg" note (docs/HOWTO_PROJECTS.md); a Ni-Al demo manifest is still missing
   (DESIGN: one demo run per project).
2. v0.2 as in docs/DESIGN.md (Stage 4, Fourier analysis, ordered domains, plotting); update setup-micromamba in CI (Node.js 20 warning).
