# Debt-register triage and the first run-gate batch

**Date:** 2026-08-23
**Purpose:** Full triage of the 49 unstruck debt-register rows; then, on Christian's
"let's start with the run only ones", execution of the cheap run gates. One new P1
defect found and registered (DR-100).

## Triage (three parallel read-only sweeps, all 49 then-unstruck rows)

Classification of what each row is *waiting on*, from the status columns as written —
no source re-verification this pass. Result recorded as a dated pass paragraph in
`todo/debt_register.md`:

- **2 closed-but-unstruck** → struck this session: DR-98 (HEX8 curl sign, closed with
  executed evidence), DR-99 (penalty input keys, closed with parse smoke). Live set 47.
- **14 run-only**, **5 ruling-only**, **14 mixed** (residue mostly runs/rulings),
  **14 genuine code** (none blocking-1.0; mostly P2/P3 hygiene).
- Of the blocking-1.0 rows (DR-02, 06, 18, 19, 22, 23, 78, 94), only DR-19 carries
  remaining code; DR-94 is pure prose + one ruling.

## Run gates executed (rebuilt tree at `0d10e6b0` + working-tree changes)

| gate | result |
|---|---|
| `make check-fast` (USE_TEST=ON, release tree) | **9/9 suites green**, 3 s, incremental build |
| DR-17 — PENTA6TS circulation tests | `test_fem --gtest_filter='*Penta6Ts*'` ran all 4 by name, **4/4 PASSED**; only Christian's sign-off remains |
| DR-23 — cohomology suite half | first-ever ctest execution of `test_homology`: 6 annulus + 2 periodic-bar fixtures green; T9 corc run + layer-tier repair fixture still open |
| DR-75 — circuit terminal-pair freeze | all three gates PASS: multi-group scratch deck aborts with the named error at `cl_ElectricalCircuitFactory.cpp:526` (label, group count, bracket hint all present); `examples/RLC_Circuit` and `examples/circuit` parse unchanged and converge |
| DR-79 (a) — keyless MUMPS stays exact | `examples/circuit` serial `-v 4`: `ICNTL(35) = 0` requested and effective, `CNTL(7) = 0.0`, zero BLR statistics blocks; every 'BLR' log hit is a parameter echo showing OFF. PASS; garber A/B still owed (needs a pre-fix binary) |
| DR-40 — costheta symmetry/current | deck ran end to end (4 ranks, full 12 s ramp, ~ -149 dB). No hidden current scaling: amplitude flows verbatim (`cl_MaxwellBoundaryConditionFactory.cpp:341` → `set_ramp`), no symmetry-conditioned factor exists in the BC path. No bisected cable: all 272 cable points strictly inside the quadrant; only Point 273 (origin bearing node) touches the axes. Awaiting strike |
| DR-84 — thin-shell layer normal B | `examples/2D_Tapestack` frame 5, layer block `tape_103_ybco` (7638 nodes): By nonzero (max 4.47e-2; identically 0 under the defect), edge-peaked (ends 3.14e-3 vs mid 1.20e-4, ~26×), consistent with nearby air (max 7.56e-2). **All three gate clauses PASS, verified by execution** |

## New defect: DR-100 (P1, blocking-1.0 candidate)

`examples/RLC_Circuit` **crashes under its own `./Allrun`** (derived np=4) with
MPI_ERR_IN_STATUS in `comm_check` (`commtools.cpp:49`) during "computing full
matrices" — root's Waitall inside `collect` serving
`SolverData::collect_matrices( true )` (`cl_FEM_DofMgr_SolverData.cpp:1850-1910`),
right after the first accepted BDF1 step.

Rank discrimination, measured: serial passes (ran 407 steps); **np=2 passes (ran 56+
steps)**; np=4 crashes deterministically. `examples/costheta` at np=4 passes the same
phase in 13 ms, and `examples/circuit` — same order-2 inductor, larger mesh — ran 40
steps clean at np=4, so the trigger is deck × rank-count and the inductor is cleared. The RLC deck's distinguishing features: 16 781
condensed hanging dofs, order-2 circuit inductor. Suspicions (unproven, recorded in
the row): the null-matrix scalar `send( tZero )` branch (`:1901-1903`, `int_t` scalar
vs root's `index_t` size recv), or a length mismatch whose guarding `BELFEM_ASSERT`
(`:1880`) is compiled out in release. Next step is a debug-build run or an MPI status
decode; any fix rides the standing plan+audit round.

## Debug-tree follow-up (same evening, after Christian reconfigured `build/`)

Christian reconfigured and rebuilt `build/` (USE_DEBUG=ON, USE_TEST=ON); debug `hphirun`
added on top. Two results:

- **DR-83 largely discharged:** the RLC deck ran under the debug binary (serial and
  np=4). Three of four revived `#ifdef DEBUG` sites executed — the SolverData
  `k < nnz` pair ran live through full-matrix assembly (armed, silent), `check_graph`
  got its first-ever firing via METIS at np=4 (no finding), `cl_SpMatrix.cpp:47`
  compiled in. Residue: the PARDISO checker — `examples/sidecoating` selects pardiso,
  but `build/` is `USE_PARDISO=OFF`, so it needs a reconfigure with `-DUSE_PARDISO=ON`.
- **DR-100 sharpened:** debug np=4 crashes identically with NO assert firing first —
  length tables consistent, failure genuinely at the MPI layer. Full demangled
  backtrace captured: `comm_check` ← sizes collect (`commtools.hpp:1099`) ←
  `collect<real>(Cell<Vector>&)` ← `collect_matrices(true)` ←
  `compute_full_matrices` ← **`Controller::save_IV`**. Two hypotheses refuted by
  source read (int_t/index_t width — both 32-bit here; tag asymmetry — `comm_tag` is
  min/max-symmetric). Next discriminator: a ~5-line per-request status decode probe
  in `comm_check` (probe policy). Blast radius: `examples/circuit` (same order-2
  inductor) ran 40 steps clean at np=4 — the inductor is cleared as the trigger.

## Second debug-tree round: PARDISO on → DR-83 fully discharged, DR-101 found

Christian reconfigured `build/` with `-DUSE_PARDISO=ON` and rebuilt. The
`examples/sidecoating` deck (the tree's only `library : pardiso` deck) ran under debug
`hphiTrun` in a scratch copy (plus the campaign's `libcustom.so` plugin): thermal
Picard 1 solved clean at −116.50 dB with the `pardisotools.f90:192` matrix checker
active — **DR-83's fourth and last site executed; the row is fully discharged** and
awaits only the strike.

Then the revival did exactly what DR-83 hoped: the revived
`SpMatrix::operator()` base assert fired on the very next assembly and exposed
**DR-101 (new P1 at the time, downgraded to P2 after the audit round below)**:
`PARDISO::initialize` flips the matrix to Fortran indexing (`cl_SolverPARDISO.cpp:240`,
`:309`) and never restores it, so `solve()`'s backup/restore never triggers and the
matrix stays 1-based forever. Blast radius source-verified: MUMPS restores correctly,
the other solvers set Cpp, and Christian's sidecoatings campaign runs petsc thermal —
**no existing result is impugned**; the exposure is the shipped example itself. The
scratch deck is the ready-made red/green gate for the fix round.

## DR-101 fixed the same night (full plan+audit round)

Christian ordered the fix. Plan pre-registered
(`tmp/ai_exchange/dr101_pardiso_indexing_base.md`), Codex + Grok audited blind in
parallel: **both land-as-written, zero blocking findings** — and both convergently
corrected the severity story, which the reconciliation adopted after re-reading the
source: `position()` (`cl_SpMatrix.hpp:836-885`) is base-aware and `allocate_values`
keeps an nnz+1 dump slot, so release `operator()` is accidentally base-safe and the
original "silent Jacobian corruption / OOB" mechanism was wrong twice. Real severity:
debug-fatal for the shipped deck, plus loud release failures in 0-based walkers
(`fn_create_graph_from_matrix.cpp`).

The landed edit (`cl_SolverPARDISO.cpp`, `PARDISO::initialize`): backup the entry
base after `Wrapper::initialize()`, restore-if-flipped after the
symbolic-factorization status check — body untouched (the `mParameters(1)` read must
keep seeing Fortran), restore-to-entry mirroring `solve()`. Gates: red captured
pre-fix (`:904` assert after thermal Picard 1); green post-fix — same debug deck runs
4+ BDF timesteps, thermal Picard −116..−118 dB, no assert; release `make check-fast`
9/9.

**Post-landing code-audit round, same session:** Codex and Grok both independently
re-read the landed diff against the pre-registered plan — **both "clean to keep",
high confidence, zero file:line defects.** Both verified line-for-line: backup
placement, `mParameters(1)` timing untouched, restore inside `#ifdef BELFEM_PARDISO`
after the status check, both `solve()` overloads unchanged, no shadowing. Grok's one
real finding was prose, not code: the register row's bold problem statement and
blocking-1.0 column still led with the refuted "silently corrupts the Jacobian / OOB
write" mechanism after the status cell had already been corrected — fixed in the same
edit (severity P1→P2, blocking flag candidate→no). Row left **unstruck** per the
DR-52/92 convention: fixed and gate-green does not strike a row that is still
uncommitted and un-blessed by Christian. Audit residue recorded in the register row
(inverted `cl_SolverPARDISO.hpp:30-31` comment; pre-existing Vector-solve double
symbolic factorization; optional lock-in test). Four voices total across two
sub-rounds (plan: Codex+Grok; code: Codex+Grok), no disagreement survived
reconciliation.

## Also observed (small, unregistered)

- `examples/scripts/Allrun:409` — `PIPESTATUS[1]: unbound variable` under `set -u`
  after the solver exits; appeared after both the RLC and costheta runs. Cosmetic
  (the run itself is complete by then), but it masks the solver's exit code.
- DR-06 is **not runnable as written**: the per-iteration min/max-T probe is gone from
  the Controller (probe-removal policy) and `cmake-build-debug/tape_hphiTrun` no longer
  exists — the deck evolved into the tapestack3d campaign. Needs Christian's call:
  re-add the probe temporarily (probes skip the audit round) or accept the tapestack3d
  campaign's hundreds of healthy coupled steps at the new rtol as the discharge.
- DR-83's original blocker (`build/`'s pinned cmake no longer existing) was resolved
  by Christian's reconfigure the same evening — see the debug-tree section above.

## Register maintenance

Same-session updates per protocol §11: DR-98/DR-99 struck; triage pass paragraph
added; status cells amended for DR-17, DR-23, DR-40, DR-75, DR-79; DR-100 appended.

## Status

DR-84's gate landed (table above). Remaining run-only rows need either campaign-scale
machinery (DR-76, 78, 87, 89, 90, 92 — memdumps, multi-hour replays in Christian's
`cmake-build-debug` campaign dirs) or Christian's build-tree reconfigure (DR-83).
Example dirs `RLC_Circuit`, `costheta`, `2D_Tapestack`, and `circuit` (mesh only)
carry tonight's run outputs; not cleaned, in case Christian wants to look.
