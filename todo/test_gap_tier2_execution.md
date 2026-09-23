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
**Status:** IN PROGRESS on Opus 5 (orchestration back on Fable) — seven of eight batches planned;
the quiz is complete (17 rulings, §3a-3d); implementation starts with slice 1 (DR-185 + T-05/N-20a,
both src → both jury rounds) and the test-only rows of B1 (one jury round on the finished diff).
B4 (mesh) waits on the peer's live `src/mesh` edits.

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
- O-B5b-1 (DR-186) the commented-out recovery guards: restore, fix the predicate, or document the
  exposure. The comment records neither the deck nor the predicate that fired.
- O-B5b-2 (DR-187) the silent thermal deactivation: refuse by name, or log it where a user sees it.
- O-B5b-3 (N-53) is a fixture worth building for TET10 orientations 2 and 3, when the only evidence
  any nondefault orientation occurs is a single assertion of orientation 1?
- O3 (T-30) Ni/Fe λ(300 K) within 15 %: land `DISABLED_` until ruled.
- O5 (T-53) logical size vs capacity after shrink.
- O9 (T-06) `.bfm` foreign git hash: refuse or flag.

## 3a. Rulings taken (quiz, batch 1)

- **Jury rounds:** a Tier 2 row that only adds tests to an existing suite takes ONE round (plan and
  code together, on the finished diff); a row that touches `src/` keeps both. Standing rule.
- **DR-186:** investigate the Lawrencium trigger first — which deck, which predicate — then restore
  the guards (real misconfiguration) or fix the predicate (false positive). No change until known.
- **DR-187:** refuse by name — `BELFEM_ERROR` naming the group and its domain type replaces the
  silent `activate( false )`. Consumer inventory before it lands (a deck that hits it today
  solves with a missing contribution; afterwards it stops at setup).
- **T-09:** route the deck arm through `set_relative_tolerance` and validate there — one gate.

## 3b. Rulings taken (quiz, batch 2)

- **T-20:** add `SolverParameters::set_max_iterations` with the siblings' validation shape; no deck
  key, so the input contract does not change. Unblocks three sparse tests.
- **O5 (T-53):** logical size everywhere — the MPI routes deliver `n_rows × n_cols` after a shrink;
  capacity is recorded by the test, never asserted.
- **O9 (T-06):** a `.bfm` written by another git hash loads with a WARNING naming both hashes; not
  refused. Consistent with the failure-tolerant input rule.
- **O1 (T-07):** `Ascii` resolves a relative path against the process cwd (`getcwd()`), not
  `getenv( "PWD" )` — a src change with a consumer inventory (every caller passing a relative
  path; the tests already pass absolute ones).

## 3c. Rulings taken (quiz, batch 3)

- **T-15:** document `gRhoMax > gRhoMin` as a precondition at the definition plus one test asserting
  it; no runtime guard on the hot clamp.
- **N-06:** the `unfold_row` clause is STRUCK; the row keeps the self-containment invariant, pinned
  by a property walk over the built dof set.
- **N-53:** reachability first — one read of how `FaceFactory` assigns slave orientation on real
  meshes; if 2 and 3 cannot occur, the `switch` arms are dead code and that is the finding.
- **T-29:** the C1 band at 10 mT is measured on the first run, then set to ten times the
  measurement (the YBCO spline-band pattern).

## 3d. Rulings taken (quiz, batch 4) — the quiz is complete

- **T-30:** the `DISABLED_` pins use each material file's own cited source (DOI) for the 300 K
  value and print measured / reference / deviation; Christian corrects a number if the file's
  citation is not the one he trusts.
- **Proceed:** implementation starts with the unblocked cheapest slice — DR-185 (kernel copy
  deletions + `static_assert`) and the shared map-message fix (T-05 / N-20a), then the test-only
  rows of B1 and B6 — one batch in flight at a time; mesh waits on the peer.

## 3e. Handoff to the next session (written 2026-09-22 evening, session 8ab15769 closing)

Christian is reloading the session on Opus 5.5. State at handoff:

- **Plans:** seven of eight on disk under `tmp/ai_exchange/review_tgap2_{b1,b2,b3,b5a,b5b,b6,b7}_plan.md`
  (ephemeral: distil what you use into this file before the sweep). B4 (mesh) is unplanned and waits
  until the peer's live edits to `cl_Face`, `cl_FaceFactory`, `cl_Element_Factory`, `meshtools` and
  `cl_CurveFactory` are committed.
- **Rulings:** all seventeen in §3a–3d; none open. One is standing policy: a test-only row takes ONE
  jury round on the finished diff; a row touching `src/` keeps both.
- **In flight at close:**
  1. Slice-1 plan jury (DR-185 kernel copy deletions + T-05/N-20a map message), astra/high, slug
     `review_tgap2_s1`. It runs detached and appends to `tmp/ai_exchange/review_tgap2_s1.md`; read the
     two verdicts, verify their citations, reconcile, then implement on Opus and take the code round.
  2. Opus agent implementing B1's five test-only rows (T-14 `tests/core/test_Logger.cpp` new, T-19,
     T-28, N-05, N-23 in existing files). If its hand-back landed before close, its report is in §3f
     below; if not, `git status -- tests/` shows what it wrote and each file must be syntax-gated in
     both trees before anything else. Gate: `test_core`, `test_linalg`, `test_math` in both trees; then
     one jury round on the own-files diff; then one commit per row.
- **Order after that:** B6's three real rows (T-29 REBCO table, T-30 disabled pins, T-22 window
  guard), then B3 with its two new setters, then B5a/B5b, then B2/B7, then B4 when the peer is done.
  One batch in flight at a time.
- **Register:** DR-185, DR-186, DR-187 sit in the working copy of `todo/debt_register.md` beside the
  peer's uncommitted rows; the three shared index files (`todo/debt_register.md`, `todo/README.md`,
  `devlog/README.md`) carry this run's lines beside peer WIP and go in at Christian's session close,
  never staged from a fresh session.
- **Do not:** launch more than two planning agents at once (incident P1); run a G4 mutation chain
  while a jury reads the tree (Tier 1 incident P5); stage anything in the three shared files.

## 3f. B1 test-only rows landed; slice-1 plan jury reconciled (2026-09-22 late, same session on Opus 5.5)

- [x] **B1 test-only rows, one jury round, one commit each:** T-14 `d91d542`, T-19 `29fce0c`,
  T-28 `dcf0d27`, N-05 `143b864`, N-23 `e9e0d48`. The implementing agent was stopped mid-write; its
  five files were gated from scratch. Both trees are RELEASE today (Armadillo `build/`, Blaze
  `cmake-build-debug/`, MKL in both): core 58, linalg 285, math 361 pass in each. N-05 skips there
  and ran in a `-DDEBUG` probe instead.
- **What the gate changed:** the agent's T-19 bound (1e-5) failed on Armadillo at 1.2e-5. Measured
  worst coefficient: Armadillo/MKL 1.2e-5, Blaze/MKL 5.9e-6, reference LAPACK 6.8e-6, OpenBLAS
  3.6e-6; bound set to 1e-4. The comment's "normal equations miss by 7e+1" was measured at 7.8e+2.
- **Mutations:** thirteen, all on scratch copies shadowing the original (harness in the session
  scratchpad, the shared tree never mutated), all killed. The gels pair fails on its own overload's
  lines only, after the calls were wrapped in `EXPECT_NO_THROW`.
- **Code jury** (`review_tgap2_b1_code`, terra/high + grok-4.7/high): no blocking defect. Folded in:
  a strict-`<` level mutation survived (new boundary case kills it), an unopenable log read as empty,
  the `mt19937` cross-library claim, four comment numbers, the Silver net is the Q (medium) piece.
- [x] **Slice-1 plan jury reconciled** in `tmp/ai_exchange/review_tgap2_s1.md`: ten findings, none
  refuted, plan approved with them folded in, no second plan round. The implementation must: drop
  the `#ifndef NDEBUG` around `test_Map.cpp:351-385`; call both const `operator()`s; mutate `Solver`'s
  copy ctor to `= default` and assert all four traits incl. `! is_move_constructible_v`; register the
  ownership test with circuit include dirs (or put that assert in `tests/circuit`); fix the three docs
  stating the bare text (`container_usage_guide.md:246`, `examples/README.md:142`,
  `input_schema.yaml:2058`); say "unknown" for `key128_t` in the FaceFactory comment; pin `time_thermal()`.
- **DR-185 row corrected** in the working copy of `todo/debt_register.md` (not staged; goes in at
  Christian's close): the compiler does not supply the move, copy assignment is already deleted by
  the const members, and `mParams` is the eighth owned resource.
- [ ] **Next:** slice 1 implementation (Opus), both jury rounds because it touches `src/`; then B6.

## 3g. Overnight run 2026-09-22/23 (Opus 5.5, same approvals as Tier 1)

- [x] **Slice 1:** DR-185 `25dbc53` (Kernel copy/move deleted; `test_OwnershipContracts.cpp` pins
  four traits for Kernel, Solver, ElementMapper, ElectricalCircuit); T-05/N-20a `a32ad4e` (the
  release lookup names the key; container guide, examples README, schema, ThinShellFactory comment).
  Plan + code juries astra/grok. The FaceFactory comment (`cl_FaceFactory.cpp:565-569`, peer WIP)
  still narrates the old bare text; with a `key128_t` key the new text says "unknown" — mesh batch.
- [x] **B6:** T-22 + T-29 `85e32c8`, T-30 `505ddec`. T-29's measurement changed the plan: the
  low-field bridge engages at all angles only from 76 K up on sp-ap → **DR-189** (working copy of the
  register, not staged) + two `DISABLED_` pins. T-30 λ pins disabled: model −1.7 % (Ni vs 90.7) and
  −1.8 % (Fe vs 83.2).
- [x] **B3:** source half `0d6b726` (T-09 validation + deck through the setter; T-20
  `set_max_iterations( int )`; input reference + schema), tests `fe52394` (T-17, T-18, T-41, T-42,
  T-57, T-76, N-22, N-28, N-29 + the T-09/T-20 setter and PETSc-cap cases). The T-18 code jury found
  a live decoder defect (an unknown MUMPS warning bit dropped next to a known one) — fixed in
  `fe52394`, report-only callers. T-57's `multiply` NaN candidate is falsified on MKL (non-MKL
  lane unmeasured). R-11 dropped (ruled).
- [x] **B5a** `90f8901`: T-15 (precondition documented in `globals.hpp` + pinned), T-21, T-31,
  T-35 (the Newton abort pinned as today's behavior — policy question in the quiz), T-37, T-46 (fem),
  T-50 (rest), T-59 (pyramids/TET35 `DISABLED_`: no interpolation order), N-01, N-04, N-06, the ghost
  dof count and `bn_angle`. Jury astra + grok (grok rerun alone after max-turns), both folded in.
- [x] **B5b** `a6e923a`: N-07, N-20 (b)/(c), N-21 (b)/(c) reachable halves, N-25, N-49, N-50, N-51.
  Held: N-20 (a), N-26 (i) (peer-edited `cl_MaxwellFactory.cpp`), N-26 (ii)/N-27 (rulings), N-30,
  N-53. Refusal decks leak → DR-175 question in the quiz.
- [x] **B2 source:** T-07 `1b92533` (Ascii resolves relative paths against the process cwd, ruling
  O1; eight stale `$PWD` test comments rewritten), T-06 `9cb58e3` (foreign git stamp on a `.bfm`
  loads with a warning naming both, ruling O9). **T-04 HELD** (see the quiz item).
- [ ] **B2 tests** (T-43 a/c, T-56, T-77): written, green in both trees, in their jury round.
- **Held for the peer:** B4 (`src/mesh` WIP), B7 (`tests/comm` WIP).
- **Gates:** `make check` 17/17 in both trees (Armadillo `build/`, Blaze `cmake-build-debug/`, both
  release) after B3's first pass; post-jury B3 fixes re-run in the Blaze tree (sparse 132 + MPI
  14 at np 2/4).
- **For the morning quiz:**
  - DR-189: is the flat clamp below ~56 K on sp-ap acceptable (value error < 1 %, tangent kink),
    or should J0 below its measured range be extrapolated? The bridge comment's "measurements
    satisfy J0 >= jc( Bmin, theta )" is false there.
  - T-30: which λ(300 K) references — Ni 90.7, Fe 83.2 (the backlog row) or 80.2 (commonly
    quoted)? Neither was read from Touloukian / NBSIR 84-3007; the pins stay disabled until you say.
  - O-B3-2: add a `set_metis_nodendp` setter (then N-22 and the sync test cover it), or accept the
    gap?
  - O-B3-4: duplicate column indices within one CSR row — file as a Tier 3 row (planner's advice)?
  - DR-175 extension: `JcFunctionDatabase`'s constructor throws after allocating `mDatabase`, so
    each refused table leaks (now exercised twice per `make check`). Widen DR-175 or new row?
  - Register header: the `[W]` count says 14; the tree has 9 open `[W]` rows incl. DR-189
    (`check_doc_claims.py` fails on it; the mismatch predates this run).
  - (done) the input-reference sentence took its Codex prose pass.
  - **DR-187, the ruling needs narrowing (consumer inventory, 2026-09-23).** Both `default:` arms of
    the thermal dispatch (`cl_ThermalFactory.cpp`, blocks and sidesets) also carry BY-DESIGN
    exclusions: a thermal kernel over a deck with `Coil`, `Inactive`, `GeometryOnly`, `Ghost` or
    `Cut` blocks, and every interface / symmetry / periodic / terminal sideset, lands there. "Refuse
    by name" on the whole arm would abort decks that run correctly today. Proposal: list the
    by-design exclusions explicitly, keep them silent, and refuse only an unexpected type by name.
    Not implemented until you confirm.
  - **T-04 HELD (two plan rounds, both "request changes").** The HDF5 existence helpers are used
    as link probes for DATASETS: `group_exists` on `abstract`, `pinned`, `orphaned`
    (`cl_Mesh_BfmFile.cpp:411-435`) and `extra` (`cl_Spline.cpp:380`), plus 16 BfmFile probes on
    names written both as groups and datasets (periodic nodes/edges/faces/facets among them). A
    type-checked `group_exists` would silently drop that data on every `.bfm` reload. Options:
    (a) switch every probe to the kind its writer uses (or `link_exists`) first, proven by a `.bfm`
    round trip carrying abstract, pinned, orphaned and periodic data, then type-check the helpers
    [recommended]; (b) leave the helpers link-level and only fix their comments and the negative
    `H5Lexists` → true conversion. Evidence: `tmp/ai_exchange/review_tgap2_b2_t04.md`. T-43(b)
    waits with it.
  - **DR-175 widening:** three more refusal paths leak what was built before the throw, and the new
    tests now exercise them in every `make check`: `MaxwellBoundaryConditionFactory` and
    `ThermalBoundaryConditionFactory` (a condition and its `SourceFunction` per refused `expk`
    deck), `MaterialFactory`'s `critical temperature` refusal (one material), and
    `JcFunctionDatabase`'s window guard (its `Database`). Widen DR-175, or accept the bounded test
    leaks until the factories get exception safety?
  - **ElementMapper Newton:** a non-convergent inverse map aborts (`BELFEM_ERROR` at
    `cl_IF_ElementMapper.cpp:369`), while the coding philosophy wants expected algorithmic failure
    as a status return, and `cl_EdgeCutter.cpp:455` already tries another candidate on `false`.
    Keep the abort (the test pins it as today's behavior) or make it a status return?
  - **By-catch, not fixed:** `BfmFile::save` throws on a mesh with no facets (`save_facet_data`);
    a Maxwell `dirichlet` section refuses any dimensioned amplitude (its unit string stays empty).
    Register rows?
  - **B5b, held for the peer:** F-1/F-2 (named errors for an unresolved sideset id and an unknown
    layer material, O2) and the DR-186 recovery guards (O4) are in `cl_MaxwellFactory.cpp`, which
    carries a peer's uncommitted edits. The B5b test-only rows go ahead after B5a.

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

## 4e. What B5b established (fem admission, circuit, reopened) — six rows refuted, two defects filed

| row | label | note |
|---|---|---|
| N-07 | CONTRACT | DR-120 already fixed the code (validation before every cast); no test exists |
| N-20 (a) | **WRONG**, + a live quality defect | a missing sideset already aborts in both builds; but in release the message is the bare "Key not found in map." — the same defect as row T-05, so the two share one fix |
| N-20 (b) | COVERED for the circuit factory; CONTRACT for the untested Maxwell and thermal copies of the same guard | `grep expk examples/` finds nothing |
| N-20 (c) | **WRONG** | `"critical temperature"` is already refused three ways |
| N-21 (a) | PARTLY COVERED | the reader end is covered; the factory consumer half is Tier 3 |
| N-21 (b) | CONTRACT | the controller always writes the slot, so the IWG constructor default cannot survive |
| N-21 (c) | **WRONG** | conflates three slots: `eta` and `chi` allow zero (zero is the off state), while `k_reg` must be dimensioned and refuses zero |
| N-25 | CONTRACT | `set_currents` already leaves the surplus abstract dof free; the test pins the August ruling |
| N-26 (i) | **LIVE** (quality) | an unknown layer material label is stored unchecked and aborts later inside the material map, naming neither the stack nor the layer |
| N-26 (ii) | **WRONG — and worse than the row** | both guards are COMMENTED OUT (`cl_MaxwellFactory.cpp:2891-2917`, "commenting these errors because they thow on Lawrencium without cause"), so nothing rejects either configuration today → **DR-186** |
| N-27 | **WRONG — and worse than the row** | an unsupported thermal group is silently DEACTIVATED, not rejected (`cl_ThermalFactory.cpp:204-207`, `:224-227`) — the gantry null-jump shape → **DR-187** |
| N-30 | **WRONG** | per-field counters since DR-96; the starvation text in the header is past tense |
| N-49 | CONTRACT | `load_state` already restores the current; no test names the class |
| N-50 | CONTRACT | HEX8TS twin at two thicknesses |
| N-51 | CONTRACT | there is no sign swap in the tree; the battery is already flip-aware and its two callers simply never flip |
| N-53 | **Tier 3** | the two-TET10 fixture cannot feed `EF_TestVolume` (wrong executable, mesh-versus-fem element, and it owns its mesh) and reaches orientation 1 only; it names the fixture it would need |

**No new refusal fires on a shipped deck.** The two fixes offered replace an anonymous abort with a
named one: same decks abort, better message. The two changes that WOULD move an acceptance set —
re-enabling the commented guards, and turning the thermal silence into an error — are carried as
rulings, not fixes.

## 5. Audit trail

Exchange threads `tmp/ai_exchange/review_tgap2_<batch>.md`; plans `..._plan.md`; gate logs
`tmp/ai_exchange/tgap2_*.log`.
