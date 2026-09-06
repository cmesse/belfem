# Todo cleanup after three files left the published tree; DR-157 and DR-158 filed

**Date:** 2026-09-05
**Purpose:** Remove every dangling reference to three todo files Christian deleted, salvage the
still-open findings of one of them, and record two non-test findings in the debt register.
**Module:** todo, tests/doc, CLAUDE.md, src/homology/doc

## What happened

Christian removed `todo/doxygen_contradiction_sweep.md`, `todo/handoff_for_gregory_20260831.md`
and `todo/test_hardening_campaign.md` as internal working material that does not belong in the
published tree, and asked for a check of `todo/README.md`.

**The index was not clean.** All three files still had full entries plus mentions in the round-3
currentness summary, and one unrelated link (`parmetis_ptscotch_wiring.md`) still pointed at the
top level although the plan had moved to `closed/` on 2026-09-04. Beyond the index,
`src/homology/doc/README.md` carried a "Handoff" section linking to the deleted Gregory handoff,
and four closed records named the files in prose (`closed/doc_currentness_fixes.md`,
`closed/handoff_20260831_open_defects.md`, `closed/test_normals_and_pipette_coverage.md`, the
DR-109 row in `debt_register_closed.md`). All fixed on Christian's instruction; devlogs left as
written. Every markdown link under `todo/` now resolves.

Side effect worth knowing: the Normals/Pipette coverage gap was tracked only inside the removed
test-hardening plan. Its closed predecessor now says the gap is real and currently untracked.

## Salvage of the test-hardening plan (read from the HEAD copy)

Question asked: are any of D1-D16 still unresolved? Answer, high confidence from reading the
tree: every *mechanism* finding landed (D8 fan-out list, D9 non-skippable rank sentinel, D10
globals in `test_commmpi_main.cpp`, D13 variable reset in `Add_Test.cmake`, D14 three
`ShareReceive` tests, D16 hardening incl. `strtol` + `EXPECT_GE` floor, D7 probe run). Every
finding that was "fixed" by rewriting a step that was never implemented is still open as a
**missing test**: R3a-R3c (`PeriodicityFactory`, cap edges after the Maxwell edge rebuild,
both-hanging `entangle`), R4/R4b (magnetic free-dof np-invariance, `synch_source_field`), R5
(`reset_timestep` / restart Δt cap), R6 (`MaterialFactory` input-section route; the existing
test uses `create_material(label, RRR)` directly, the route the plan said cannot catch the
defect class), the SuperLU solver test (INC-409 residue), O5 (no prescribed-partition path in
`mesh::Partitioner`). Those are not re-filed anywhere; Christian excluded missing tests.

## Filed

- **DR-157 `[MIXED][P]` P2** — the nightly GitLab job never configures `-DUSE_GASMODELS=ON`.
  Christian ruled it should on 2026-08-30 and the plan recorded O2 as resolved, but all four
  matrix jobs in `.gitlab-ci.yml` still configure without it, so `test_gastables`/`test_gasmodels`
  never exist and the TARGET guard skips them silently. About 4,650 of ~4,975 physics test lines
  have never run in CI. One token in Christian's CI file, then observe the ctest log.
- **DR-158 `[CODE][W]` P3** — the Tier 2 `mpirun` is found via `find_program` (hints, then PATH)
  and never checked against the MPI the binary was linked with; `find_mpi.cmake` identifies the
  compiler's `mpi.h`, not the launcher. Crash-not-silent (the sentinel and verdict fold catch a
  singleton launch), so watch only. Was recorded only in the removed plan.

`check_doc_claims.py` 38/38 after the recount (`[P]` 8, `[W]` 3).

## Documentation fixed

- `CLAUDE.md` "Honest testing posture": "runs the **full** test suite" replaced by what the job
  does (a debug/release × Armadillo/Blaze matrix of `make check`, Tier 2 suites at 2 and 4 ranks,
  4-rank cases oversubscribed on the two-core runner), and the gas-model gap named as the fourth
  limit. This was R7's documentation half, never written.
- `tests/doc/tests_06_comm.md`: binary name corrected (`test_commmpi`, not `test_comm_mpi`),
  main file, ctest registration and the sentinel described, suite table refreshed from the source
  (35 tests in 8 suites; `CommReduce` and `ShareReceive` were missing, `CommSendRecv` had 9 of 16).

Reviewed, not verified: no build or test run in this session.
