# Release pre-check for 0.9.3: jury over the commits since the 0.9.2 regression review

**Date:** 2026-10-04
**Purpose:** Overnight read-only check before tagging 0.9.3: of the 129 commits since `v0.9.2`, the first 89 were juried on 2026-09-23 (`dl20260923_regression_review_devel_vs_092.md`, no P0); this round covers the ~40 after `9c9428a` (61 files, +2999/−727 under `src/`, `tests/`, `examples/corc_periodic_bc`) and the release records for the whole range. Question: is there a blunder that should block or be fixed before the tag?
**Module:** `src/fem/kernel`, `src/fem/iwg`, `src/fem/maxwell`, `src/fem/thermal`, `src/mesh`, `src/math/graph`, `src/sparse`, `src/physics/materials`, `src/io`; `tests/` read; `CHANGELOG.md`, `CITATION.cff`, `CMakeLists.txt`
**AIs involved:** Claude Opus 5.5 (pre-registration, verification, reconciliation); Codex `gpt-6-astra`/high (kernel), Codex `gpt-5.6-terra`/high (physics, io+release); Grok `grok-4.7`/high (all three). All three Codex legs of the 23:45 dispatch failed on the ChatGPT usage limit (cause only in `_cross_review.log`); they were re-run with identical prompts after the 00:03 reset by a detached runner (`review_rel093_codex_rerun.sh`), so each exchange file carries one "Audit FAILED" Codex entry beside the real one
**Status:** COMPLETE — no P0. Three P1s (stale rho cache, create_field double free, release records), all three closed the same day (§Fixes) together with the ruled μ0 band and committed in one commit
**Verification:** reviewed, not verified, for the findings. Executed this session: `make check` **19/19** on HEAD `aaee5ad` in a private tree (release, SCLS `mkl`, Armadillo, `-DUSE_GASMODELS=ON`, shared libs) — the first run of the gastables/gasmodels suites on this state, which the nightly CI never configures (DR-157); the debug/MKL/Blaze private-tree run of 2026-10-03 (17/17 on `b01911f`, `src/` byte-identical to HEAD) stands. `check_doc_claims.py` 38/38, `check_wrapper_policy.py` clean, `check_practice_drift.py --strict` exit 0 (P2 notes only), `check_citations.py` exit 1 on the known "Amazon Linux 2023" false positive. Release build: two `-Wfloat-conversion` warnings, both pre-0.9.2 (`cl_GM_EoS_Hydrogen.cpp:182` `mNab` is a `real` count; `gasmodels/main.cpp:109`)
**Exchange:** `tmp/ai_exchange/review_rel093_kernel.md`, `review_rel093_physics.md`, `review_rel093_io.md` (frozen pre-registration, both auditor entries, verification, reconciliation); subjects `review_rel093_*.subject.md`

## Summary

Nothing found that makes 0.9.3 compute a wrong answer on a shipped deck, and no shipped
`examples/*/input.conf` newly aborts (all three voices checked the decks touched by the new
refusals). Every formulation change re-derived holds: the ghost-penalty SWIP operator
(ρ_w(∂h_m+∂h_s) = φ for continuous flux; its manufactured-flux test separates all four
operator forms at 1:1, 1:10, 10:1), the Modified-Kim field derivative, the coincident-knot
tangent limit, the Debye overflow-free identities, the QUAD16 facet fix. Three things should
be dealt with before the tag.

## Must address before tagging (P1)

1. **Stale rho cache defeats the Kohler fixes** (Codex + Grok independently; confirmed by
   source trace). `rho_database_is_current` (`fn_rho_database_is_current.hpp:47`) only asks
   whether the HDF5 file has an `RRR` dataset, a marker already shipped in 0.9.2; the cache
   name is label + integer RRR (`cl_Material_Metal.cpp:893`). A warm `Lead_RRR*.hdf5`,
   `WhiteTin…`, `Silver…`, or an alloy cache built from them (PbSn solder; same test at
   `cl_Material_Alloy.cpp:516`) is reloaded, so DR-191's switch fix and the silver low-field
   cubics silently do not apply: a cold and a warm run give different ρ(B). Fix options: a
   content-revision dataset bumped now (and on every future table-law change), or at minimum
   a CHANGELOG instruction to delete `*_RRR*.hdf5`. No test catches it
   (`test_MetalFieldDependence.cpp:211` builds without tables).
2. **`Kernel::create_field` refusal path double-frees hanging dofs in throw mode** (Codex
   astra, single raiser; confirmed by source trace). a3bdbb4 made the field owned by a
   `unique_ptr` until the push; `set_equation → collect_hanging_dofs` aliases `mDOFs` entries
   into `mHangingDOFs` (`cl_FEM_DofMgr_DofData.cpp:1957`) and can refuse inside
   `create_dofwise_t_matrices_master` (e.g. `:3333`) before the split (`:2014` serial,
   `:2336` parallel); `DofData::reset` (`:78-88`) then deletes both lists. 0.9.2 leaked here;
   HEAD double-frees. Debug/gtest/parallel-debug only (release aborts before unwinding), and
   the refusals in that window are rare — but the symptom (a double-free crash that buries the
   refusal text) is exactly what the ownership change set out to avoid. Cheap fix.
3. **Release records** (all three voices; still open from 2026-09-23). `CHANGELOG.md`
   `[Unreleased]` lists only the rutherford deck and the SCLS page; `CITATION.cff` is 0.9.2 /
   2026-09-18 while CMake says 0.9.3. Missing: the intended result changes (prism side-face
   increment 2×, thin-shell recovery rewrite, jc/n below 10 mT, ghost SWIP operator, pure-metal
   cp on a 1 K grid, silver/lead/tin Kohler, YBCO λ below 1 K, Kim and piecewise tangents,
   QUAD16 facet 1), the controller changes (floor stall suspension can accept a step 0.9.2
   rejected; the iteration line now prints the head residual — log parsers notice), the
   corc_periodic_bc deck (tolerance 1e-10, Anderson on), and the new refusals (typed sideset
   on a periodic plane; multi-source Dirichlet alias; RRR-token/key and builtin/curve
   conflicts; connector without its layer block; non-linear or μ≠μ0 conductor beside a thin
   shell; incomplete periodic node match; non-unit cochain coefficient; HDF5 rank ≠ 1).

## Rulings for Christian (design, not settled by vote)

- **Multi-source Dirichlet refusal** (`cl_FEM_SideSet.cpp:53`) also refuses the consistent
  case where every source sits on the same constant-Dirichlet sideset (Codex; my K-6). Grok
  found no refiner in `src/mesh` that produces such a dof today. Keep, or admit compatible
  sources?
- **μ = μ0 band is absolute** (`cl_MaxwellFactory.cpp:2924-2927`, `< BELFEM_EPSILON`). No
  shipped deck or plugin sets μ; 4π·1e-7 passes (Δ ≈ 7e-16), a hand-rounded 1.256637e-6
  aborts. Relative band?
- **Alloy cp still resampled on a 4 K grid** (`cl_Material_Alloy.cpp:315`) after the metals
  moved to 1 K; whether PbSn cp undershoots below 4 K needs a probe (Grok, single raiser).

## P2 (pre-existing or minor)

- `load_matrix_from_file` has no rank check (`hdf5_tools.hpp:838-849`) beside the new 1-D
  guards (Codex + Grok).
- Leak-named tests (`test_AdmissionRefusals` FreesTheFirst ×2, the material-factory pair)
  assert only the refusal text; valgrind is the real observer (both).
- Periodic matcher (`cl_Mesh_PeriodicityFactory.cpp:742-770`): my "writes past n" was refuted
  by Codex (flagged targets are skipped); a double hit plus a miss can still pass
  `tCount == n` with a non-bijective map. Still in no DR row.
- Iteration-line format change; connector refusal fires per holding rank (debug contract);
  `n()` constant fallback skips the n ≥ 1 clamp; `FacetTables` compares node sets, not order.
- `check_citations.py` red on a false positive; two pre-existing float-conversion warnings.

## Process notes

Pre-registration was written before dispatch. Of my own frozen findings, K-12 (overflow) was
refuted by Codex; K-13 (does the SWIP test pin the fix) I answered myself before the auditors
returned; my P-6 scale doubt was not taken up by either voice and stays low. I had missed the
two P1 defects, the stale cache and the double free — both sat one call away from the diff,
which is where the auditors earned the round. `git status` was clean after every Grok leg.

## Rulings (Christian, 2026-10-04, question form)

- rho cache: **content-revision marker** in the cache, bumped now and on every table-law change.
- create_field double free: **fix before the tag**, with a refusal regression under valgrind.
- multi-source Dirichlet refusal: **keep**; list it among the CHANGELOG refusals.
- μ = μ0 band: **relative band**.

## Fixes (same session, `todo/release093_precheck_fixes.md`)

Christian ruled the P1s and the μ band, approved plan v2 with source edits, and chose 1e-6 for the
band and unit tests for the double free (the full-path fixture filed as `test_gap_followups.md` R13).

- **Rho cache revision.** `gRhoDatabaseRevision = 1` in `fn_rho_database_is_current.hpp`; both
  savers write a `revision` dataset; the probe loads a cache only when `RRR` exists and `revision`
  is a one-element **integer** dataset equal to the constant, never throws, and closes every
  identifier. Every 0.9.2 cache rebuilds once with a warning naming the revision. Tests:
  `RhoDatabase.AnotherRevisionIsRebuilt` (different, two-entry, fractional, missing revision; the
  rebuild leaves a current file) and `RhoDatabase.AMetalTableCarriesTheRevision` (the Metal saver,
  which an alloy never reaches because it builds components with tables off).
- **DofData::reset deletes each dof once.** Unflag the non-null hanging dofs, flag the non-null
  `mDOFs`, delete only unflagged hanging dofs, then `mDOFs`. Tests: `DofDataReset.*` (shared,
  disjoint, null slots, stale flags) with a counting `Dof` subclass as the oracle. Positive
  control: on unfixed HEAD `SharedDofIsDeletedOnce` dies with SIGSEGV. The real refusal path is
  not reachable from node fixtures — `Mesh::finalize` flattens node source chains
  (`cl_Mesh.cpp:2012`), probed twice — so it needs a hanging-edge fixture (R13), which should also
  pin the pre-existing `tBefore` leak and the dangling mesh slots both auditors found on that path.
- **μ0 band.** `MaxwellFactory::mu_is_mu0`: relative 1e-6, a rounding allowance (seven-digit μ0
  passes, six-digit and |χ| ~ 1e-5 refuse). Test `Mu0Band.*`. Alloy components keep the absolute
  band (out of scope).
- **Records.** `CHANGELOG.md` `[Unreleased]` gained Changed / Refusals / Fixed for the whole
  v0.9.2..HEAD range, after a Codex language sweep (`gpt-5.6-luna`/medium); one sweep edit that
  moved the "not linear" refusal from the conductor to the shell was reverted.
  `CITATION.cff` is left at 0.9.2 until the tag date is known.

**Juries.** Plan round: Codex `gpt-6-astra`/high + Grok `grok-4.7`/high, no P0; both raised the
null-slot hazard in the flag pass (applied) and asked for a real refusal regression (probed,
routed, ruled). Code round: Codex astra/high on the pre-P3 diff found two P3s (a fractional
revision converted by `H5Dread` and passed; `std::vector` in a test); the first Grok leg died on
HTTP 402 "Grok Build usage balance exhausted", reported at once; Christian topped up and the Grok
leg re-ran on the post-P3 diff: one P3, a stale "roundoff band" comment, fixed. No P0/P1 in code.

**Gates (executed, Claude's private trees).** `make check` 17/17 (debug, SCLS `debug`, reference
BLAS, Armadillo) and 19/19 (release, SCLS `mkl`, Armadillo, gas models), zero warnings, after the
P3 fixes; valgrind 0 errors on `DofDataReset.*`, `Mu0Band.*`, `RhoDatabase.*` (debug tree).
`check_doc_claims` 38/38, `check_wrapper_policy` clean, `check_citations` the known false positive.

## Owed

- Christian: `CITATION.cff`
  version and date with the tag; his own `make check` in the shared trees if he wants one. The commit was made at his request from the follow-on session.
- `todo/test_gap_followups.md` R13 (hanging-edge refusal fixture); `/code-review ultra` as a
  breadth pass from the workstation terminal (it does not launch over Remote Control).
