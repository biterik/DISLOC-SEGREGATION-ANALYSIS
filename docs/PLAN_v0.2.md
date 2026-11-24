# Plan for disegreg v0.2 -- brief for unattended work by the Claude Code session (Erik, 24 Nov 2026)

Erik wants v0.2 developed in parallel and largely without him. These rules make that safe; they come on top of CLAUDE.md.

## Working rules
1. SEPARATE FOLDER. Never switch branches in /Users/e.bitzek/DEVEL/DISLOC-SEGREGATION-ANALYSIS: the simulation instances read its
   docs, and the M5 conda env `disegreg` runs its code (editable install). It stays on `main` (v0.1.0 + docs-only commits).
   Develop in a git worktree: branch `v0.2-dev`, folder /Users/e.bitzek/DEVEL/DISLOC-SEGREGATION-ANALYSIS-v0.2, with its own conda env
   `disegreg-dev` (from environment.yml, editable install of the worktree). Never install the worktree into the env `disegreg`.
2. NO RELEASE. Commit freely on `v0.2-dev` and push that branch (GitHub Actions runs on it). Never push to main, never tag, never
   touch Zenodo. Merging and the v0.2.0 tag happen only with Erik.
3. NO DEFINITION DECISIONS. D1 rev. 4 stays unchanged. Every new quantity gets a written proposal in docs/DEFINITIONS.md (new
   section, marked [PROPOSAL v0.2], label to be decided: D1 extension or D2) and is listed in docs/OPEN-DECISIONS_v0.2.md with
   options, numbers from the demo runs, and a recommendation. Code for undecided quantities is built and tested, but its output
   columns carry `status = proposal` and the results table refuses them for role `headline`. Then continue with the next item;
   do not stop and wait.
4. Computations only on the M5 (natively) -- no cluster access, no cloud workspace. Use the snapshot copies on the Mac and the
   M5 cache /Users/e.bitzek/Desktop/SIMULATIONS/DISEGREG-CACHE (read; new entries allowed). Never write into project folders.
5. Every item: tests (synthetic exact cases first), CHANGELOG [Unreleased] on the branch, CLAUDE.md state on the branch.
   Keep the full test suite green; if a v0.1.0 result would change, stop that item and write it into OPEN-DECISIONS.
6. Erik's style rules still apply to everything he will read (abbreviations spelled out, plots >= 13 pt, colour-blind-safe,
   SYMMETRIC colour scales for signed maps; step-by-step recipes when he has to act).

## Scope, in this order
A. Stage 4 across runs (DEFINITIONS section 9 -- already defined, so no proposal needed except where marked):
   - run sets from manifests: same project/potential/character/T, ordered by Delta; bottom-up vs top-down pairs at the same Delta;
   - uptake exponent n_up = Delta ln c_C / Delta ln x_ref between neighbouring offsets, errors propagated, for each reference
     (c_B, x_iso, x_bc) as in the v0.1 tables; ACE Ni-Cu keeps paper_flag;
   - McLean prediction from the 0 K site-energy tables (Ni-Cu thread 09, Nb-Au thread 06; read-only; an adapter per table
     format), compared with the measured c_C; [PROPOSAL] which atoms count as core sites (r < 5 A of the relaxed 0 K line?);
   - hysteresis: Gamma(20 A) and c_C of bottom-up vs top-down at the same Delta; [PROPOSAL] the criterion for "different"
     (e.g. |difference| > 2 x combined block-averaging error);
   - automatic x_bc: a bulk-control run's own manifest (no dislocation) feeds x_bc into the runs that name it.
B. Plotting module (`disegreg plot`): Gamma(R) curves, time series with averaging window and drift verdict, c_C and Gamma(20 A)
   vs Delta per character, uptake-exponent plot, region bar charts; Erik's house style (see rule 6).
C. Along-line core composition: [PROPOSAL] solute count in the core in X bins (5 A?), its Fourier spectrum against binomial
   noise -- does the core hold a continuous line state or a periodic array of particles?
D. Ordered-domain sizes for Ni-Al: [PROPOSAL] definition of an L1_2-ordered solute (e.g. no solute 1st neighbours and >= k
   solute 2nd neighbours), domains = 2nd-shell-connected ordered atoms, size as equivalent diameter, in the core vs the far
   field; against the 2.7 +- 0.6 nm of Howard et al. 2025. Synthetic test: an L1_2 Ni3Al particle of known size in fcc Ni.
E. Anything the project instances report in docs/FEEDBACK_from-*.md files after 24 Nov: list it in OPEN-DECISIONS; fix bugs on
   the branch (a real bug in v0.1.0 -> write it at the top of OPEN-DECISIONS: it may need a v0.1.1 on main, with Erik).

## When Erik comes back
docs/OPEN-DECISIONS_v0.2.md is the single page he reads: what was done (short), what needs his decision (each with options and a
recommendation), what is blocked. Keep it current after every item.
