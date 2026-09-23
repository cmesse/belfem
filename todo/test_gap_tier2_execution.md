# Test-Gap Tier 2 Execution: 48 Triage Rows Plus the Six Reopened Exclusions

**Date:** 2026-09-22
**Purpose:** Run Tier 2 of `tmp/devlog_sweep/triage.md` (§2, 48 rows) plus N-49 … N-53 reopened by
the exclusion re-read (`todo/test_gap_backlog.md` §3b; N-54 stays Tier 3) through the standing
loop, batched by suite: one plan per batch (a section per row), one plan jury per batch, Opus
implementation, gates in both trees, one code jury per batch diff, one commit per row. Consumer
inventory before any src fix that changes a value or refuses an input; rulings batched for
Christian at the end (the quiz).
**Module:** every suite under `tests/`; src fixes only where a row's defect is live and the jury confirms
**AIs involved:** Claude Fable then Opus 5 (orchestration, planning agents, juries, commits), Codex + Grok (juries; terra/high, astra/high on the critical list), Opus (implementation)
**Status:** IN PROGRESS on Opus 5 — 2026-09-22 18:30 the Fable credit window closed and five of
eight batch planners died mid-read; Christian moved the session to Opus and the work continues.
Three plans survive (`review_tgap2_b1_plan.md`, `_b2_plan.md`, `_b7_plan.md`); B3, B4, B5a, B5b
and B6 are being re-planned two at a time with named file lists (incident P1). Nothing has been
juried, implemented or committed for any Tier 2 row yet. B4 (mesh) waits until the peer's live
`src/mesh` edits settle — a planner citing `meshtools.cpp` line numbers mid-edit would cite
line numbers that move.

> **Standing approvals (carried from Tier 1, `memory: test-gap-overnight-run`):** tests and jury-confirmed
> local src fixes land without a per-row stop; Claude runs make/ctest; one commit per row, own files only,
> peer WIP never staged (`git diff --cached` empty before every add; shared CMake lines staged apart);
> a red not fixed in one attempt, a vendor quota failure, or an open ruling → register row and move on.
> One `make` per tree. No G4 mutation chain while a jury reads the tree (P5).

## 1. Batches

| batch | suite(s) | rows | plan slug |
|---|---|---|---|
| B1 | math, linalg, containers, core | T-05, T-14, T-19, T-28, T-53 (linalg half), T-54, N-02, N-05, N-23, N-52 | `review_tgap2_b1` |
| B2 | io, homology | T-04, T-06 (O9), T-07 (O1), T-43, T-56, T-77 | `review_tgap2_b2` |
| B3 | sparse | T-09, T-17, T-18, T-20, T-41, T-42, T-57, T-76 (MKL lane), N-22, N-28, N-29 | `review_tgap2_b3` |
| B4 | mesh | T-03, T-39, T-45, T-46 (mesh half), T-48, T-55, N-03, N-08 | `review_tgap2_b4` |
| B5 | fem, circuit, interpolation | T-15, T-21, T-31, T-35, T-37, T-46 (fem half), T-50, T-59, N-01, N-04, N-06, N-07, N-20, N-21, N-25, N-26, N-27, N-30, N-49, N-50, N-51, N-53 | `review_tgap2_b5a` / `_b5b` |
| B6 | physics, gasmodels | T-16, T-22, T-27, T-29, T-30 (O3), T-58, N-24 | `review_tgap2_b6` |
| B7 | comm | T-40, N-28 (comm half), T-53 (commmpi half) | `review_tgap2_b7` |

## 2. Ledger

| id | batch | re-read on `devel` | gate | commit | residue |
|---|---|---|---|---|---|
| (filled per row as it lands) | | | | | |

## 3. Rulings owed (the quiz)

- O1 (T-07) `Ascii` relative path: cwd or `$PWD`? — carried from the triage.
- O-B3-1 (T-20) `SolverParameters` has no `set_max_iterations`; add one (then the one-iteration cap
  is testable) or leave the clause untested.
- O-B3-2 (T-09) the deck arm bypasses the setter: validate in both places, or route the deck
  through the setter so there is one gate.
- O-B6-1 (T-30) the λ(300 K) reference values for nickel and iron, with DOIs.
- O-B6-2 (T-29) the C1 band at `mBmin` = 10 mT.
- O-B5a-1 (T-15) guard `gRhoMax > gRhoMin` at runtime, or document it as a caller precondition
  (the planner recommends documenting; `std::clamp` is UB on an inverted window).
- O-B5a-2 (N-06) strike the `unfold_row` clause: the symbol does not exist in `src/`.
- O3 (T-30) Ni/Fe λ(300 K) within 15 %: land `DISABLED_` until ruled.
- O5 (T-53) logical size vs capacity after shrink.
- O9 (T-06) `.bfm` foreign git hash: refuse or flag.

## 4. Incidents

- **P1 (credit exhaustion, 2026-09-22 ~18:20).** Eight batch planners were launched in parallel;
  three finished, five died on `rate_limit` / out-of-credits partway through their reads. The shape
  was too expensive: each planner re-reads a wide slice of the tree, and eight at once multiplied
  that. Cheaper shape for the resume: one planner at a time, a named file list in the brief instead
  of "find it", and the row's own suite file read first so a row that is ALREADY COVERED costs one
  read rather than ten.

## 4a. What B1 established (the one batch whose agent reported)

| row | label | note |
|---|---|---|
| T-05 | LIVE (diagnostic) | release branch of `Map`/`OrderedMap::operator()` is a bare "Key not found in map."; the debug branch names the key. Four sites, one `BELFEM_ERROR` each |
| T-14 | CONTRACT | `cl_Logger.hpp:122-126` fixed 2026-09-13, unpinned; new `tests/core/test_Logger.cpp` |
| T-19 | CONTRACT | `polyfit` tests stop at degree 2; the measured degree-8 fixture is 1.5e-7 on the Blaze route |
| T-28 | CONTRACT | the "Silver high" net's `Qy` is NOT monotone, so the `xi_by_y` half needs random monotone nets |
| T-53 (linalg) | LIVE (resize) + CONTRACT (stride) | both backends' `Vector::operator=( expression )` is a bare assignment; ruling O5 |
| T-54 | CONTRACT / behind a seam / BLOCKED | `residual_string` is a file-static needing a seam; the `throw_on_error` default is unobservable in a gtest binary; the Anderson third needs F2 |
| N-02 | ALREADY COVERED | `Cardano.CardanoDoubleRoot` and two siblings assert it |
| N-05 | CONTRACT | one debug-guarded length-mismatch case |
| N-23 | CONTRACT | both `gels` overloads suppress identically; tests cover values only |
| N-52 | **WRONG — strike the row** | `MatrixDebug.SubmatRowOutOfBoundsThrows` / `SubmatColOutOfBoundsThrows` (`test_Matrix.cpp:492-507`) already pin the bounds half; the exclusion re-read missed them |

## 4b. What B6 established (physics)

| row | label | note |
|---|---|---|
| T-16 | **CONTRACT** (planner said LIVE; Claude re-read and downgraded) | the six rho/lambda field-dependency flags are set inside `populate_rho_database` (`cl_Material_Alloy.cpp:509-515`), reached only under `mComputeTables`; the comment at `:504-508` says that is deliberate — the flags advertise that the database exists, and an alloy without it has no field-dependent path, unlike a metal whose Kohler flags are analytic. Today's T-51 test builds alloys with tables off and is therefore not pinning a defect |
| T-22 | **LIVE (physics half)** | the window guard is `cl_JcFunction_Database.hpp:124-130` and is a coverage test, not a width test; nothing constructs a `JcFunctionDatabase` today. The `bn_angle` identity needs a `Calculator` and moves to the fem batch |
| T-27 | **BLOCKED** | `tests/physics/CMakeLists.txt:60-64` gates the gas suites on `USE_GASMODELS`, OFF in both trees (DR-157) |
| T-29 | **LIVE, cheap** | `share/material/sp-ap.hdf5` ships and resolves through `material::data_file()` with `gBelfemDataPath` set by `Add_Test.cmake:41`; axes T 3→91 K, log10 B −2→1.6, angle 0→π in 1° steps, so `mBmin` really is 10 mT |
| T-30 | **LIVE, ruling-gated** | the two λ(300 K) magnitudes land as `DISABLED_` pins printing measured / reference / deviation; the Ni magnetization half is plannable now (`cl_Material_Nickel.cpp:58-67` cites Crangle & Goodman 1971, 10.1098/rspa.1971.0044) |
| T-58 | **WRONG as written; half covered, half blocked** | `≤ ρₙ` holds for Riva and Piecewise — the latter already covered by `HtsLaws.PiecewiseNormalRegimeAndAboveTcrit` — while `rho_powerlaw` is unbounded by design and clamped in `MaxwellData::compute_rho()`. The `share/fluid` half is gastables, blocked with T-27 |
| N-24 | **ALREADY COVERED**, one-line residue | `tests/io/test_filetools.cpp` pins the precedence ladder; its own comment names the residue — `tests/io` cannot include `physics/materials`, so the `"material"` forwarder is inspection-only. One physics case closes it |

B6 totals 16 tests across one new file and three existing ones; five rulings, two of them physics
VALUE questions (the λ(300 K) references with DOIs, and the C1 band at `mBmin`).

## 4c. What B3 established (sparse) — 42 tests, four rulings, two of which block three of them

| row | label | note |
|---|---|---|
| T-09 | **LIVE, two halves** | `set_relative_tolerance` (`cl_SolverParameters.cpp:190-195`) assigns with no validation while both siblings validate (`:199-206`, `:294-303`) — and the deck arm (`:64-70`) writes `mRelativeTolerance` directly, bypassing the setter, so validating the setter alone would not cover a deck. Both confirmed by direct read. Consumers: `strumpacktools.cpp:225`, `cl_SolverPETSC.cpp:36`, `:518`, `cl_FEM_Controller.cpp:3452` |
| T-17 | CONTRACT | every clause implemented; tests drive the public `set_strumpack_options`; the retired AUTO / explicit-only policy is not encoded |
| T-18 | CONTRACT | `warning_message` is public and decodes the bits including the ambiguity arm; the free phase copies all 80 INFO slots; the −9 ladder stays owned by `SolverPolicy.MumpsWorkspaceActionTable` |
| T-20 | CONTRACT, one clause blocked | the named defect class is CLOSED — a negative convergence reason is read, Allreduced and hard-errored (`cl_SolverPETSC.cpp:212-310`). No `set_max_iterations` exists, so the one-iteration cap is untestable → ruling |
| T-41 | CONTRACT (two clauses of three) | the enum round trip and `check_naupd` are sparse; the restart-budget floor lives in a private fem method and is not a sparse row |
| T-42 | **LIVE coverage gap** | `grep -c SUPERLU tests/sparse/test_Solver.cpp` = 0 with `USE_SUPERLU=ON` in both trees. Seven tests; the mutation needs an ASYMMETRIC matrix or it is invisible on a symmetric Laplacian |
| T-57 | three labels | `rcond` BLOCKED (`USE_SUITESPARSE=OFF` both trees); the duplicate-index throw is a CONTRACT through a debug-only `check_graph`; the `multiply` NaN-prefill is a LIVE CANDIDATE — the 2-argument overload fills `y` on the MKL path and says why, the 5-argument one does not |
| T-76 | CONTRACT | `USE_PARDISO=ON` in both trees; three restore sites, with a skip for the non-MKL arm |
| N-22 | CONTRACT | extend `CopyPreservesTuningMembers` (the iteration budget is unchecked and `metis nodendp` default-only, as the test's own comment confesses), plus a `synchronize()` round trip over all payload slots |
| N-28 | CONTRACT | the re-initialization guard exists and is deliberately first; the test pairs the throw with a pool-occupancy assertion so a guard moved below the create still fails |
| N-29 | CONTRACT | `init_defaults()` runs only in the constructor, so a pre-initialize arm survives; the `set_symmetric` half is fem and is named, not planned |

**Scope flag:** two rulings ask for new public setters (a max-iterations setter for T-20, and the
T-09 validation shape), so B3 is not a pure test batch — three of its tests are blocked behind
them. Its riskiest assumption is that `<StrumpackSparseSolver.hpp>` can be included from
`tests/sparse` at all; SuperLU's header is known to collide with both matrix backends
(memory `project_superlu_header_conflicts`), and if STRUMPACK's does too, the T-17 block falls
back to reading options off a constructed wrapper.

## 4d. What B5a established (fem contracts) — 13 rows, one live fix, one phantom

| row | label | note |
|---|---|---|
| T-15 | CONTRACT | `gRhoMin = 0`, `gRhoMax = 1e10` are set in one place and narrowed nowhere in `src/`, so the clamp is the identity on the row's window. Latent: `std::clamp` is UB on an inverted window and nothing enforces the ordering → ruling |
| T-21 | PARTLY COVERED | `test_GhostSwitch.cpp:147-165` already pins tag-differs-on-eta; the `…,11` / `…,13`, swapped-pair and other-`input curves` cases are open |
| T-31 | CONTRACT | the stride carries its own anti-aliasing rationale; both guards are debug-only, so the oracle is a sentinel round trip |
| T-35 | **WRONG as written** | a non-convergent Newton does not return false — `cl_IF_ElementMapper.cpp:369` is `BELFEM_ERROR( tCount < 100, … )`, so it aborts or throws (verified). The rest of the row is a contract |
| T-37 | CONTRACT, cheap | `IWG` has no pure virtuals, so the probe subclass is three lines |
| T-46 (fem half) | MOSTLY COVERED | `test_BdfTimestepMethod.cpp` covers the history levels and the warm-start round trip; open is `reset_fields` restoring the pre-shift step, which needs a field-linked IWG. The "upcoming Δt" clause is a Controller property, out of this batch |
| T-50 | CONTRACT ×2, one WRONG, **and the batch's only live fix** | `Solver`, `ElectricalCircuit` and `ElementMapper` already delete copy; **`Kernel` declares only a constructor and destructor while holding seven raw owning pointers** (verified: `cl_FEM_Kernel.hpp:104-112`), so it is implicitly copyable and a copy would double-free. Four-line header fix, compile-time, no consumer breaks since everything holds a `Kernel *`. `Controller::time()` already returns a const reference; the "Gram at h = 1e-4" clause is a red herring (that Gram is in-plane) |
| T-59 | PARTLY COVERED | `test_FacetIntegrationPoints.cpp:150-205` already does the mapped-point identity over every slave orientation for the 3-D types and both PENTA TS types. Open precisely: QUAD4TS / QUAD9TS run with orientations OFF, there is no 2-D volume entry so the edge-facet-versus-1-D case is untested, and `populate( 0 )` → `auto_integration_order`. The QUARTIC gap is in `mesh::interpolation_order()`, i.e. row T-03 |
| N-01 | CONTRACT | a layer is tagged `Buffer` iff its material lacks `rho`; `create_buffers` is private, so the test asserts the observable outcome |
| N-04 | CONTRACT, all four clauses hold | Buffer shares the Air φ branch; BufferSymmetry / AntiSymmetry register; the thermal IWG accepts Buffer |
| N-06 | CONTRACT + **phantom** | **`unfold_row` does not exist anywhere in `src/`** (verified: zero hits) → strike that clause. The self-containment invariant is real but unenforced: at `cl_FEM_Dof.cpp:133` `this` is not yet hanging, so a self-source passes the assert, and the assert is compiled out in release — the test must be a property walk over the built dof set |
| ghost dof count | CONTRACT | the IWG sums master and slave per entity class; the existing test pins the count through a different function |
| bn_angle identity | CONTRACT | public on the calculator; `acos(x) + acos(−x) = π` holds including the degenerate branch, where both sides return π/2 |

**Cross-cutting:** `tests/fem/test_fem_main.cpp:22` sets throw-on-error unconditionally, so
`BELFEM_ERROR` throws in that binary under release too — which makes the `#ifndef NDEBUG` guard
in `test_BdfTimestepMethod.cpp:296` stale. The correct guard for an assert-oracle test is
`#if BELFEM_ASSERTIONS_ACTIVE`, not `#ifndef NDEBUG`.

## 5. Audit trail

Exchange threads `tmp/ai_exchange/review_tgap2_<batch>.md`; plans `..._plan.md`; gate logs
`tmp/ai_exchange/tgap2_*.log`.
