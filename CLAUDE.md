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
- Review of docs/DEFINITIONS.md rev. 2 (incl. the TURLO-REDO additions: line coverage, fallback core axis, x_iso without bulk control,
  PTM icosahedral label; proposal Gamma(20 A), Gamma(40 A) as headline numbers).

## State (22 Nov 2026, evening)
- Done: steps 1 and 2. M5 conda env `disegreg` (Python 3.12, OVITO 3.17.2 from conda.ovito.org, created with
  `mamba env create -f environment.yml --channel-priority flexible`; a-CNA + DXA smoke test passed). Git repo on GitHub
  biterik/DISLOC-SEGREGATION-ANALYSIS (public), branch main.
- cmmg facts: SLES 15 SP7, glibc 2.38 (OVITO 3.17.2 wheels/conda builds OK). Existing cmmg env /cmmc/ptmp/$USER/conda-envs/ovito_env has
  OVITO 3.14.1 -- leave it alone; disegreg gets its own env under PTMP from environment.yml (recipe still to write), plus
  env-lock/disegreg-linux-64.txt.
- OVITO pin 3.17.2 agreed; if regression shows DXA differences vs the old scripts, rerun one snapshot with the old OVITO version.

- Step 3 done: package skeleton, definitions/D1.yaml (inside the package: disegreg/definitions/), geometry.py + tests, CI.

## Next steps
4. io_ovito.py + cache.py; verify the DXA-internal-labels question on one real snapshot; then regions/composition/order; regression report.
