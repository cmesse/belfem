# Currentness Sweep of `todo/` and the Debt Register (second pass)

**Date:** 2026-08-11
**Purpose:** Re-sweep `todo/` for currency and dispose of what is no longer live — plans to
`closed/` or `deferred/`, Status lines corrected, checkboxes ticked where reality had already
ticked them, register rows struck where the work is genuinely done.
**Module:** meta (task tracking)
**AIs involved:** Claude (sweep), Codex (prose pass), Christian (authorisation)
**Scope:** bookkeeping only — **no source edited.** Writing to `./todo/` and `./devlog/` is
allowed under the read-only default; source edits were explicitly out of scope.

---

## 1. The rule this pass ran under

The recurring failure in `todo/` is that register rows and Status lines are copies of *devlog
claims*, not checks against the *tree*. So the rule for this pass was: for every row and every
checkbox touched, the evidence is the source file, `git log`, or the example deck — never the
plan, never the devlog, never the handoff that commissioned the sweep.

Protocol §11 vocabulary applies throughout. Nothing is called **verified** unless an executable
gate actually ran — and where one did, it was run by an earlier session, not by this one (§8).
Static work is **reviewed**.

## 2. What moved

| file | to | the evidence for the disposition |
|---|---|---|
| `bfm_stale_cache_detection.md` | `closed/` | `hphirun` end to end after Christian's build: tag round-trips through HDF5, a 1.6 → 1.7 µm thickness edit rebuilds *naming the changed value*, a restored deck reuses silently, a genuinely pre-feature `.bfm` rebuilds exactly once (`2c0a7382`) |
| `example_deck_and_material_db_repair.md` | `closed/` | every §6 DoD item ticked against a run — `mpirun -np 2` builds the database in-run and marches timesteps, the genuine old-format `Copper_RRR30.hdf5` warns and regenerates, `examples/README.md`'s gmsh line executed verbatim, D4's replacement deck run from a clean directory to t = 300 ms with a warm restart |
| `coreduce_…performance_findings.md` | `deferred/` | its own status line has said "deferred" since 2026-07-03; DR-29 is deferred; the one live item is sequenced behind rectification tests that are only part done |
| `slepc_eigensolver_integration.md` | `deferred/` | deprioritized *by measurement* on 2026-08-10, before this sweep — the move records a decision, it does not make one |

Neither deferral is a judgement this pass invented. Both files already said what they were; the
directory did not.

**One residual was carried rather than dropped.** `example_deck_and_material_db_repair.md`'s O2 —
the `Tape_Quench/BuiltinMat` deck, dead since 2026-04-08 because its custom material is labelled
`buffer`, a label reserved by `8eae5f00` — is now `debt_register.md` **DR-65**. Closing a plan must
not be a way of losing its open questions.

## 3. The finding worth reading: DR-45(a) was worse than recorded

Every prior correction in this directory went one direction. DR-43 carried four open gasmodels
defects; all four were already fixed. DR-19 listed items already ticked in the plan it points at.
B6 sat done-but-unticked for four weeks. DR-52 was written up as an open P0 straight from a jury
devlog and found already fixed when someone checked the tree. Four for four: **the debt was already
paid and the register had not noticed.**

That quietly builds an expectation — a stale row is a *pessimistic* row. This sweep found the
counter-example.

Commit `5a2ddf81` fixed the worker-warning gate in `MUMPS::free` and in `MUMPS::solve`. But there
are **two** `solve` overloads — one takes a `Vector` RHS, one takes a `Matrix` — and only the first
was fixed. The Matrix-RHS overload still gated its warning branch on `this->rank() == 0`, so a
worker-raised warning was still discarded there. Including `INFO(1) = +1`, which is the entire
subject of DR-45's open half: MUMPS found indices outside the matrix and **dropped those entries**,
meaning it factorised a different matrix than the one assembled.

*(No line numbers, deliberately. The anchors this sweep read were the pre-`78ea534d` ones and
`78ea534d` rewrote that whole region — citing them now would hand a reader four numbers that no
longer mean what the sentence says. In the current tree the overloads are `cl_SolverMUMPS.cpp:243`
and `:401`, both routing positive `INFO(1)` through `MUMPS::check_warnings()` at `:1225`, and no
rank gate survives on a warning path. This is the file's own lesson applied to itself.)*

Grok raised this in the in-flight DR-45(b) jury round. This sweep did not take the vote as the
finding — it opened `cl_SolverMUMPS.cpp` and read both overloads. The claim held.

**Provenance matters here more than the defect does.** The row said "both `INFO` reads
(`MUMPS::free`, `MUMPS::solve`)" and had said so for a week. That phrasing is not wrong so much as
*singular where the code is plural*, and a singular noun is exactly the kind of thing a reader
confirms by nodding rather than by grepping. There are three `INFO` sites, not two.

The session that owns the row implemented both clauses and committed them (`78ea534d`) within the
hour, routing positive `INFO(1)` from both `solve` overloads through a new private
`MUMPS::check_warnings()` that raises `BELFEM_ERROR` on bit `+1` and reports the rest per rank;
`free()` stays report-only because it is teardown. **DR-45 is struck, on Christian's explicit instruction**, as a recorded exception to the register's
"fixed-but-still-gating rows stay live" rule — the same call he made for DR-42 and DR-49 earlier
the same day, and on the same grounds: the design work is finished and the residual is one run, not
a decision, a design question, or an unwritten test. **Struck ≠ verified applies with full force
here.** The fix compiles clean under `-DNDEBUG` and `-DDEBUG` with `-Werror` and has never
executed, and its regression risk is a *false abort on a healthy deck* — precisely the failure mode
of the `set_currents` guard reverted in `61bbfa11`, which had aborted every parallel run with a
current BC. That is a run gate, no amount of reading closes it, and the strike does not touch it:
it survives in the row's status column, which is why that column is never struck.

**A correction that is itself an instance of what this file documents.** The sentence above first
read "the `set_currents` guard that broke every parallel run on 2026-08-10 (DR-24)". Both halves
were wrong, and I inherited both from the DR-45 row rather than checking them. The revert is
`61bbfa11`, dated **2026-08-11**. And **DR-24 is not that incident** — it is the periodic free-cut
`mPrescribedCurrents` row, struck and closed by Christian on 2026-08-09; the parallel-run abort has
**no register row at all**. The association is not arbitrary, which is exactly why it survived: DR-24
had asked for a one-time assert on `set_currents`, and the reverted guard was a failed attempt at
precisely that, so the ID *reads* right at a glance. It survived a jury round, the row edit that
introduced it, and my own read. A plausible ID copied forward without being followed is the same
defect class as a Status line copied from a devlog — this sweep's whole subject, committed by the
sweep itself. Caught by the session that owns DR-45, not by me.

## 3b. Three instances of one defect class, in one day

Read separately, the findings above are three unrelated slips. Read together they are one:

| instance | the pointer | what it pointed at |
|---|---|---|
| stale line anchors | `cl_SolverMUMPS.cpp:238/:408/:563`, and the ~15 drifted sites in the rho/lambda doc inventory, and `fix_facet_masters` at `:1114` | the right thing, until the file moved under the citation |
| the DR-24 mislabel | a register ID | a *different, closed* row — the incident it meant has no row |
| "both `INFO` reads (`MUMPS::free`, `MUMPS::solve`)" | a prose enumeration | two of three sites; the noun was singular where the code is plural |

**In every case the pointer was correct when written and was never re-followed.** That is the
whole mechanism. None of the three is a reasoning error, none would have been caught by thinking
harder about the claim, and all three survived readers who were paying attention — the DR-24
mislabel got through a jury round, a row edit and two independent reads, because a plausible ID
*reads* right at a glance and confirming it costs a `grep` that nobody spends.

`CLAUDE.md` already carries the countermeasure for the first row, learned the expensive way from
`doc/input_file_reference.md`: **anchor by a searchable token, never a line number.** What today
suggests is that the rule was scoped too narrowly. It was written about `file:line` citations in the
input contract, but a register ID, a commit hash, a devlog filename and a prose enumeration are all
the same object — a reference a reader confirms by recognition rather than by lookup. The line
number is merely the instance that rots fastest.

Whether to widen the convention is Christian's call, and it is recorded here rather than acted on.
Observation contributed by the session that owns DR-45.

## 4. What was corrected without moving

- **`instruction_doc_currency.md`** — R4 was **done and unticked**. DR-63 was closed in `7432385a`
  by *amending the convention* rather than manufacturing 22 reference sections:
  `doc/documentation_guidelines.md:312` now names `doc/literature_references.md` as the single
  source of truth for citations, says module docs link to it, and records that an empty References
  heading is worse than its absence — so 1-of-23 is the expected state, not a gap. Also corrected:
  the plan said `scripts/check_doc_claims.py` guards 21 claims. It guards **32**, and re-running it
  in this sweep gives 32/32.
- **`iterate_refactor_plan.md`** — all three surviving bugs re-confirmed live and every anchor
  re-baselined. The ω-ordering asymmetry was read out of the source rather than carried over:
  `iterate_magnetic` sets then clamps (`cl_FEM_Controller.cpp:1569`/`:1570`), while both other
  paths clamp then set (`:875`/`:914` and `:1756`/`:1759`).
- **`rho_lambda_argument_convention.md`** — the doc-pass inventory was re-baselined and **grew**.
  Three of six files had drifted and the list was short by roughly fifteen sites, including a full
  `Metal::rho(real B, real beta, real T)` signature at
  `materials_contracts_and_invariants.md:314`. One **source** site joins them:
  `cl_Material.hpp:981`'s doc-comment still describes "three-parameter properties like
  `rho(B,angle,T)`" while the declaration 400 lines above it is `rho( const real T, const real B,
  const real beta )`. Comment text only — recorded, not fixed.
- **`maxwell_kernel_collapse_plan.md`** — both R12 sub-items still open; the architecture note is
  still unwritten — re-grepped: nothing under `src/fem/maxwell/doc/` or `src/fem/kernel/doc/`
  documents the *collapsed kernel set*. (`h_ghost` is named in three module docs, but only as the
  ghost-penalty stabilization, never as one of the five kernels that replaced 38 variants. The
  looser claim "no doc names any of the five" would have been false.). But DR-02's **Gate A instrumentation now exists**: the env-gated
  assembly dump is committed in `bc578b5e` at `cl_FEM_Controller.cpp:921` and `:1577`, writing the
  assembled system before the solve when `BELFEM_DUMP_SYSTEM` is set. It is no longer a
  commented-out call site. That does not revive Gate A — Christian cut it because the legacy
  kernels exist only in history — but a future equivalence check now has its mechanism.
- **Anchors** in `controller_picard_tapestack_regression.md` (`try_escalate_to_newton` `:542` →
  `:600`; the retracted-D10 residual site `:2286` → `:2788`, re-read and still holding) and
  `thin_shell_overhang_bug_analysis.md` (`fix_facet_masters` `:1114` → `:1243` — **+129 lines in
  two days**, which is the standing argument for locating by symbol).
- **The README index entry for the PID controller** was stale where the plan file was not: the
  index still said "uncommitted, not compiled" for work committed in `1d6ef305` on 2026-08-09. A
  reminder that the index drifts independently of the files it indexes.

## 5. What the tree confirmed unchanged

Six register rows were re-read at their cited sites and each was exactly as written: DR-07's
`real dbdT = 0.0` placeholder (`cl_FEM_Calculator.hpp:2879`), DR-13's `uint16_t` counters
(`cl_Vertex.hpp:62-63`), DR-26's `mPairVerdict` built at `cl_CutProcessor.cpp:785` and consumed
nowhere, DR-46's `DISABLED_GhostElementContract`, DR-53's two `false` fuse flags
(`cl_ThinShellFactory.hpp:347-348`), and DR-64's signed dot-product fix live at both producer
sites.

One claim did **not** survive its challenge intact, and the failure is instructive because the
conclusion was right. The FVM plans say "no `src/fvm` code change since 2026-07-09". `git log --
src/fvm` shows two later commits: `e4c02ac8` (2026-08-06) edits `src/fvm/CMakeLists.txt`, and
`4f2c11cd` (2026-08-07) adds `src/fvm/doc/README.md` plus four other `.md` files. No `.cpp` or
`.hpp` under `src/fvm` has changed since `1173454d`, so the plans' *verdict* — the module is
stalled — holds exactly. But "code change" is the wrong words for it, and a reader who runs that
`git log` concludes the plan is stale when it is only imprecise. My own first draft of this devlog
called both later commits documentation-only, which is also wrong: one is a build-system change.
Corrected in both places.

## 5b. Two test failures fixed, and one defect found underneath them

Not part of the sweep as commissioned — `make check` was run while it was in flight and
`test_linalg` came back with two failures, both `GesvdWorkBufferReuse`, both **real types only**
(`float`, `double`; the two complex instantiations passed).

**The test was wrong and the wrapper was right — but only about correctness.** The test helper
`gesvd_min_work()` encodes reference LAPACK's documented minimum workspace,
`max( 3·mn + max(m,n), 5·mn )`, and the test asserted that the size the wrapper *queried from
LAPACK* meets it. Measured with a standalone probe against the linked library: **MKL's `dgesvd`
query returns 7 for a 4×3 `'A','A'` problem where that formula gives 15**, and MKL then completes
with `lwork = 7` and `info = 0`. The reference minimum is *sufficient*, not *necessary*, and a
queried size may legitimately fall below it. Asserting one against the other is simply invalid.

Underneath that sat a genuine defect, now **DR-66**. The wrapper's acceptance test for a
caller-supplied buffer uses the same reference minimum. So the buffer it grows on its own first
call (7) is smaller than its own threshold (15): every later call re-enters the query branch,
issues an extra LAPACK call and two `set_size` round trips, and across changing shapes the shared
buffer *shrinks* — 4×3 then 3×2 takes it 7 → 5. "Reuse one buffer across a loop", the entire
reason the overload exists, never happens for a caller who sized from the query. Correctness is
untouched, which is exactly why only a length assertion caught it.

**The repair does not assert the defect away.** `GesvdWorkBufferReuse` now checks results and
usability instead of buffer lengths, with the reason stated inline; a new
`GesvdWorkBufferReuseAtReferenceMinimum` sizes its buffer so the reuse branch *is* taken and pins
it across three iterations. Green on all four LAPACK flavors, 257/257 in `test_linalg`. Whether to
lower the acceptance floor to the queried size or keep the conservative one and document that
query-sized buffers are not reusable is a design ruling, and DR-66 is open on it.

One method note, since this file is about checking rather than assuming: my first probe printed
`gesvd(...)` and `Work.length()` as two arguments of the same `printf` and reported the buffer
never growing at all. Evaluation order between the two is unsequenced, so the length was being
read before the call ran. The corrected probe — statements, not arguments — gave a completely
different and correct picture.

*(Also fixed in the same round, from a `make check` break in `tests/math`: `tDet` in
`fn_quaternion_from_rotation_matrix.hpp` was consumed only by a `BELFEM_ASSERT` and became an
unused variable under the `-DNDEBUG` the test tree compiles with, breaking `-Werror`. The
determinant now forms inside the assert, so it also stops being computed at all in release.)*

## 6. What this sweep could not settle

- **`make check-fast` has still never run.** `USE_TEST` is OFF in the shared tree, so the spline
  suite, `Hex8TbUnitCirculation` and the new `tests/math/test_GraphSpfa.cpp` have never executed.
  Every plan whose residual is that gate stays active with the gate named. Builds are Christian's.
- **Every "pending run" row is beyond a static pass by construction** — DR-02's Gate B thermal
  coverage, DR-52's greg3 A/B, DR-19's five R7 checks, DR-64's ferro A/B, and
  DR-15/18/22/34/38/40. A bookkeeping pass can confirm the code is in and the anchors are right.
  It cannot manufacture evidence that only a run produces.
- **DR-45's run gate** — helix, serial and on ≥2 ranks, completing without the new hard error
  firing.

## 7. Dispositions handed back to Christian

Unchanged from the 2026-08-09 sweep and deliberately not acted on: **FVM in or out of 1.0**
(re-checked — `add_subdirectory( fvm )` is still commented out at `src/CMakeLists.txt:12`, and the
last change to any `.cpp`/`.hpp` under `src/fvm` is still 2026-07-09, `1173454d`; stalled is not
deferred);
**`ngspice_parser_plan.md`**, a feature proposal with no consumer pressure that reads like a
`deferred/` candidate but is a product call; and the **`[seeded — confirm]` correction pass** on
`debt_register.md`, where roughly half the unstruck rows still carry severities and blocking-1.0
flags that are proposals from the 2026-08-05 seeding rather than Christian's judgement. New this
pass: **DR-65**, whether the `Tape_Quench/BuiltinMat` deck is still wanted.

## 8. Verification status

**This sweep built no BELFEM target and ran no solver deck.** Two exceptions, both narrow and both
in §5b: standalone probes and a standalone `test_linalg` binary were compiled and run from the
scratchpad against the prebuilt `.a`s, using the build tree's own flags, to settle the gesvd
question. That is executable evidence for §5b's claims specifically — the MKL query returning 7,
and 257/257 green — and for nothing else in this file. Its own claims are static source
traces at the `file:line` cited, plus `git log` / `git show` results — the middle of the evidence
ladder (protocol §11), not verified results. One executable check did run: `scripts/check_doc_claims.py`
(32/32), which checks documentation claims against the tree, not code behaviour.

The distinction that matters when reading §2: where this file says a plan closed *on run evidence*
— `hphirun` round-tripping the `.bfm` tag, `mpirun -np 2` building a rho database, a 12-timestep
corc run — that evidence was **recorded by the sessions that produced it**, and this sweep's
contribution was to confirm it exists and matches the plan's own definition of done. Relying on a
recorded run is a rung above a static read and a rung below running it again.

A collision worth recording: four of this sweep's `git mv` renames were staged at the moment they
ran and rode into a concurrent session's commit (`78ea534d`). The moves are intact and disclosed
in that commit's message; nothing was lost. `git mv` stages immediately, which is a trap for any
two sessions sharing a working tree.
