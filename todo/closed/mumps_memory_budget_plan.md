# MUMPS Memory Budget: Probe the Machine, Hand MUMPS ICNTL(23), Keep the Ladder Behind It

> **CLOSED 2026-09-03** (todo/ currentness sweep, round 3): code, tests and input contract landed 2026-09-01; the remnant is the R12 cold-start gate (`tape_quench` under MUMPS) and a Darwin compile. Status lines and checkboxes below are as they stood at closure and are not maintained.

**Date:** 2026-09-01
**Purpose:** When MUMPS runs out of factorization workspace (`INFOG(1) = -9`), BELFEM today walks
a blind percentage ladder (`ICNTL(14)`: 30 → 60 → 120 → 240) and, when that is exhausted, hands
the step to the timestep controller. The ladder's ceiling is a compile-time constant, the machine's
actual memory never enters the decision, and the log no longer says how far short MUMPS was. This
plan adds a cgroup-aware probe for available memory, reduces it to a rank-uniform per-process
budget, and hands that budget to MUMPS as `ICNTL(23)` — the strict bound MUMPS documents for the
case where raising `ICNTL(14)` stops helping — **with the `ICNTL(14)` ladder kept behind it** for the
residual `-9`/`-8` the guide says can still occur, and give-up reserved for `-19`.
**Module:** `src/sparse` (MUMPS wrapper + Fortran shim), `src/core` (probe), `src/comm` (MIN
reduction), `doc/input_*` (deck override)
**AIs involved:** Claude (exploration + plan + code), Codex + Grok (jury on the plan 2026-09-01,
jury on the code)
**Status:** IN PROGRESS — **code landed 2026-09-01 (R1–R11), code jury done, fixes applied; R12 gate
owed (Christian builds).** Plan jury 2026-09-01 (Codex `gpt-5.6-terra`/high, Grok `grok-4.6`/high):
P0 + eleven P1 folded in. Code jury same day, same depth: no P0; four P1 and seven P2 confirmed and
fixed, one refuted; one design change (the probe moved to `initialize()`, closing O4). Supersedes
`todo/closed/mumps_workspace_cap_raise.md`. Exchange: `tmp/ai_exchange/review_mumps_memory_budget.md`.

> **Scope guards:**
> - STRUMPACK, PETSc, PARDISO are OUT of scope. The probe is generic; only the MUMPS wrapper
>   consumes it in this plan.
> - Out-of-core factorization (`ICNTL(22)`) is OUT of scope. In-core only.
> - The `soft fail → controller cuts Δt` contract is KEPT. This plan changes what happens
>   *before* the soft fail, not the soft fail itself. `-19` under the controller is a soft fail
>   (expected algorithmic failure, never an abort).
> - No change to the STRUMPACK-first solver guidance.
> - Darwin code is written to the Mach documentation. **It is not the tree's `__APPLE__` idiom**
>   (that is `sysctlbyname`, `cl_SolverWrapper.cpp:372-379`) and **cannot be compiled or run in
>   this session** (Linux host). Owed in §7.

---

## 1. Current Behaviour and How It Fails

**The ladder.** `cl_SolverMUMPS.cpp:33-34` fixes `gDefaultMemoryRelaxation = 30` and
`gMaxMemoryRelaxation = 240`. `MUMPS::escalate_workspace()` (`:534-613`) doubles
`mIParameters( MemoryRelaxation )` on every rank after a rank-uniform `INFOG(1) = -9`, clamps
at the cap, and re-enters the collective with `JOB = 5`. It returns false for any other code
(`:541-544`), so `-8` (integer workarray IS too small — same remedy per the guide) is never
retried. The escalated value persists on the instance (`:557-566`); `free()` restores the default
(`:444`). There is no deck key for any of this.

**After the ladder.** A `-9` that survives the cap reaches the soft-fail arm (`:815-848`): the
wrapper flags failure, drops the analysis (`mMatrix = nullptr`) and returns; the controller cuts
the timestep. That contract was written for `INFOG(1) = -10` (`cl_FEM_Controller.cpp:113-121`).

**What the machine never learns.** Nothing in `src/` queries physical memory (grep for
`sysinfo`, `meminfo`, `host_statistics`, `getrusage`: no hits). The CPU-count probe in `Wrapper`
(`cl_SolverWrapper.cpp:284-385`) is the structural precedent: platform-guarded, affinity-mask
aware *because* a SLURM cgroup narrows the mask (`:288-295`), `ifstream` allowed on the setup path
(`:302-304`), unknown platform → 0 ("no budget, therefore no warning and no advice"). Two things
about it must **not** be copied: it sits under `#ifdef OMP` (`:284`), and its Darwin branch is
`sysctlbyname`.

**What the log no longer says.** The Aug 30 binary printed the raw DMUMPS lines with the
`INFOG(2)` shortfall; the Sep 1 binary silences the MUMPS streams below `-v 4`
(`mumpstools.f90:360-401`) and `print_soft_fail()` deliberately omits `INFOG(2)`
(`cl_SolverMUMPS.cpp:1755-1760`). The soft path — the only path the controller ever takes —
reports the exhausted percentage and nothing quantitative.

| Failure | Mechanism | Evidence |
|---|---|---|
| Cold start could not launch | `-9` at 240 on the first factorization of step 1, then every Δt cut (0.5 → 0.031 ms) failed identically | **transcript only**: the 12:55 `tape_quench` run, read in the 2026-09-01 session; its `out.txt` was overwritten by the 13:13 relaunch (jury round 1, Grok). The surviving log shows the *relaunch*, where step 1 climbs 60/120/240 and succeeds (`tape_quench/out.txt:112-125`) and soft-fails start at step 9 (`:252-258`) |
| Warm runs limp | 129 `-9` events across usermat runs 2-7, each rescued only by a Δt cut | `tape_quench_usermat/out.txt`, per-run counts 15/29/68/11/3/3 (python scan, 2026-09-01) |
| Raising the percentage does not close the gap monotonically | shortfalls 1 905 777 → 606 895 → 3 425 627 → 4 243 613 entries at 60/120/240/240 % | usermat run 2, `out.txt:1914-1924` (old binary, raw `INFOG(2)`) |
| The budget is unknown to the code but large | 8 ranks hold 4.6 GB RSS **including** the 240 % workspace; 44 GB available of 62 GB | `ps`/`free` measured 2026-09-01 13:28 |
| Diagnostics regression | soft path prints no shortfall; raw stream silenced | `cl_SolverMUMPS.cpp:1755-1760`, `mumpstools.f90:360-401` |
| Persist ratchet disarms the retry | once at 240, every later solve gets one attempt | `:566-569`; `todo/closed/mumps_workspace_cap_raise.md` F2 |

**Bottom line:** the ladder alone is the wrong instrument — MUMPS's guide says `ICNTL(14)` "does
not lead to a bound on the total memory allocated … we recommend the use of `ICNTL(23)`" (MUMPS
5.9.0 user guide §2.11, `tmp/userguide_5.9.0.txt:446-450`; installed library 5.9.1 per
`/opt/scls/mkl/include/dmumps_c.h`). But the same page says `-9` "may still occur … and as before,
ICNTL(14) should be relaxed" (`:437-443`). So: cap first, ladder behind it, give up on `-19`.

## 2. Architecture: Cap, Then Ladder, Then Give Up

MUMPS 5.9.0 guide, `ICNTL(23)` (`:2288-2324`): "maximum size of the working memory in MegaBytes
that MUMPS can allocate per working process … If `ICNTL(23)` is greater than 0 then MUMPS
automatically computes the size of the internal workarrays such that the storage for all MUMPS
internal data does not exceed `ICNTL(23)`." Accessed "at the beginning of the factorization
phase", so a `JOB = 5` retry picks it up without re-analysis. **A positive value on every rank is
"interpreted locally on each MPI process"**; the guide's lower bounds (`INFOG(16)` full-rank
in-core, `INFOG(36)` BLR factors, `INFOG(44)` low-rank CB) are stated for the host-only case. Since
`INFOG(16)` is the max over ranks of the per-process `INFO(15)`, a rank-uniform value ≥ `INFOG(16)`
satisfies every local bound — conservative, and now stated (jury round 1, Codex).

**`ICNTL(23)` is a cap, not a fill** (`:434-436`: "MUMPS will no longer try to allocate all the
memory authorized"). The earlier draft's claim that setting it "allocates the budget" was an
overstatement (Grok). Consequence for O5 below.

**The policy** (`next_workspace_action`, a pure function):

| state on a `-9` / `-8` / `-19` | action |
|---|---|
| budget known, `ICNTL(23)` slot still 0, budget ≥ bound | set the cap, retry `JOB 5` |
| budget known, slot 0, budget **<** bound | give up (the machine cannot hold the estimate; climbing the ladder raises the request further — Grok P1) |
| slot set, code `-9` or `-8`, ladder below ceiling | ladder rung (double `ICNTL(14)`), retry `JOB 5` — the guide's residual-`-9` instruction |
| slot set, code `-19` | give up — the cap cannot be met |
| ladder at ceiling | give up |
| budget unknown (0) | existing ladder, unchanged (now also on `-8`) |

Ladder ceiling `gMaxMemoryRelaxation`: 240 → **480** (one more rung; the defensible minimum from
`mumps_workspace_cap_raise.md` O1, decided 2026-09-01, Christian).

**Rejected alternatives** (unchanged in substance): raise the constant alone (bounds nothing,
non-monotone); derive a percentage from the probe (lossy re-encoding of `ICNTL(23)`).

## 3. Gap Table

| # | State / behaviour | Needed for | Handled today? | Class | Citation / rationale |
|---|---|---|---|---|---|
| 1 | Available physical memory, Linux | budget | no | (c) | `MemAvailable` from `/proc/meminfo` — **in kB, × 1024** (Grok P1). NOT `sysinfo().freeram` (free = 2 GB vs available = 44 GB on this box) |
| 2 | cgroup memory limit, Linux | budget under SLURM | no | (c) | v2: `/sys/fs/cgroup/memory.max` − `memory.current` (walk `/proc/self/cgroup` path; `max` = unlimited); v1: `memory.limit_in_bytes` − `memory.usage_in_bytes`. Budget = min(host available, cgroup headroom). Same class of bug the CPU probe closes (`cl_SolverWrapper.cpp:288-295`) (Codex + Grok P1) |
| 3 | Available physical memory, Darwin | budget | no | (c) | `host_statistics64(HOST_VM_INFO64)`: `free_count + inactive_count` × page size; release the host port. Owed on a Mac |
| 4 | Unknown platform / probe failure | safe fallback | n/a | (c) | return 0; never `BELFEM_ERROR`; not under `#ifdef OMP` |
| 5 | Ranks share a node's pool | per-process budget | `gComm.node_size()` (`cl_Communicator.cpp:196-204`) | (a) | `/ node_size()` |
| 6 | Live MUMPS instances on this rank | no double booking | `gNumSolvers` is occupancy on the Fortran side (`mumpstools.f90:217-225`), **not exported** | (c) | new `mumpstools_num_solvers()` bind(c) accessor; `/ max(1, count)`. Rank-uniform because create/free are collective (`cl_SolverMUMPS.cpp:268`, `:363`). The second instance on the reference deck is the conditioning shift-invert solver (`cl_FEM_DofMgr_EigenValues.cpp:1208`), not thermal (PETSc) |
| 7 | Rank-uniform decision | collective safety | ladder keys on `INFOG` (`:538-544`) | (c) | reduce **`int_t` MB** (no `comm_type<size_t>`, `commtypes.hpp:151`) with a new `allreduce_min` in `commtools.hpp` (L-21); 0 on any rank → 0 everywhere |
| 8 | Lower bound under BLR | correct bound | `INFOG(36)` is beyond the 40 copied (`mumpstools.f90:484-486`, `:204-205`; `cl_SolverMUMPS.cpp:51-52`); arrays are 80 (`dmumps_c.h:110`) | (c) | copy 80 at both shim exits and `free_solver` (`f90:165`); `set_size( 80 )`; hpp comments (`mumpstools.hpp:42,63-64`). **0-based:** `INFOG(16)` = `mInfoG(15)`, `INFOG(36)` = `mInfoG(35)` |
| 9 | `ICNTL(23)` slot in the shim | set the cap | none; `aIParameters` `dimension( 13 )` (`f90:263`, unpack `:322-336`), enum `mumpstools.hpp:91-106`, `set_size( 13 )` (`cpp:46`) | (c) | `Parameter::MemoryBudget` (14th); guarded write like ICNTL(14) (`f90:421-423`); `mumpstools_create_solver` does not take the array (`f90:51`) |
| 10 | `-8` and `-19` reach the policy | correct retry | `escalate_workspace` ignores both (`:541-544`); `error_message()` decodes `-19` (`:1245-1250`); `print_soft_fail()` does not (`:1744-1790`) | (c) | policy accepts `{-9,-8,-19}`; soft-fail decode for `-19` and `-8` |
| 11 | Shortfall decode overflow | safe print | `error_message()` `:1246` does `-aInfo[1] * 1e6` into `int_t` (int32 unless `BELFEM_INT64`, `typedefs.hpp:47-53`) | (c) | print "N million entries" when negative, never multiply in `int_t`; fix the existing arm too (Codex P1) |
| 12 | Cap and ladder persist per instance | no repaying per step | ladder persists (`:557-566`), reset in `free()` (`:444`) | (a) | same for the slot |
| 13 | Deck override | operator control | no key | (c) | `memory budget : 4000 ;` (integer, **MB**, dimensionless key via `get_int` — `unit_to_si` has no `MB` token, `stringtools.cpp:407-470`; `get_real` is fatal on a non-real, `cl_Input_Section.cpp:291-308`). Explicit value wins over the probe and is set **before the first factorization** (O5). `SolverParameters`: field, setter/getter/`have_`, **copy-ctor** (`cl_SolverParameters.cpp:119-140`), `synchronize()` width 12 → 13 (`:323-341`). Schema + reference same session |
| 14 | Soft path prints the numbers | diagnosis | omitted (`:1755-1760`) | (c) | second box row: estimate, cap, shortfall |
| 15 | User-facing doc row | honesty | `src/sparse/doc/solver_memory_and_compression.md:249` documents the 240 % ladder as current | (c) | update in the same session; language sweep |

### 3.1 Cross-cutting findings

- **Every decision keys on `INFOG`, never `INFO`** (`:534-544`). The probe is the only rank-local
  input and is reduced before any branch. `memory_budget_mb()` is therefore **collective** —
  called from every rank inside `escalate_workspace()`, which already is.
- **MegaBytes are 10⁶ bytes** in every MUMPS statistic (§5.12). Convert with 1e6.
- **Retry is `JOB = 5`, never `JOB = 2`** (`:584-593`).
- **The policy function is free and testable**: `belfem::solver::mumps::next_workspace_action()`
  in `mumpstools.hpp`/`cl_SolverMUMPS.cpp`, not a private static (Codex P2, Grok).

## 4. Ordered Steps

- [x] **R1 — Probe.** `src/core/fn_available_memory.{hpp,cpp}`: `size_t belfem::available_memory()`
  → bytes, 0 = unknown. Linux: `MemAvailable` × 1024, min with cgroup v2/v1 headroom when a limit
  is set. Darwin: Mach, port released. `#else` 0. Added to `src/core/CMakeLists.txt` `SOURCES`.
  Test: `tests/core` (`fast`): `> 0` on Linux.
- [x] **R2 — MIN reduction.** `belfem::allreduce_min` beside `allreduce` (`commtools.hpp:231-260`),
  same static_assert / serial identity. Test in `tests/comm/test_CommMPI.cpp` (np 2/4, label `mpi`;
  **not** in `check-fast` — stated, not pretended).
- [x] **R3 — Instance count.** `mumpstools_num_solvers()` bind(c) returning `gNumSolvers`; declared
  in `mumpstools.hpp`. (independent)
- [x] **R4 — Budget.** `MUMPS::memory_budget_mb()`: probe → `/ node_size()` → `/ max(1, instances)`
  → `× gMemoryBudgetSafety (0.5)` → `/ 1e6` → `int_t` → `allreduce_min`. 0 if any rank probed 0.
  (after: R1, R2, R3)
- [x] **R5 — INFO/INFOG width 80.** Both shim exits, `free_solver`, `set_size( 80 )`, hpp
  comments. (independent; smallest diff, do first)
- [x] **R6 — Shim slot.** `Parameter::MemoryBudget`, `dimension( 14 )`, `set_size( 14 )`, guarded
  `ICNTL(23)` write. (after: R5)
- [x] **R7 — Policy.** `next_workspace_action( aInfoG1, aRelax, aBudgetSlot, aBudgetMB, aBoundMB )`
  → `{ Ladder, Cap, GiveUp }` per the §2 table; `escalate_workspace()` calls it; bound =
  `mInfoG(15)` for `ICNTL(35) ∈ {0,3}`, `mInfoG(35)` for `{1,2}`; ceiling 480 with the three
  comment sites (`:28-31`, `:437`, `:558`); slot reset in `free()`. Box rows for the cap and for
  each rung. (after: R4, R6)
- [x] **R8 — Diagnostics.** `print_soft_fail()`: arms for `-19` and `-8`; second row with
  estimate / cap / shortfall in the millions convention. `error_message()` `-9`/`-19`/`-8` arms:
  no `int_t` multiply. (after: R5)
- [x] **R9 — Deck override.** `memory budget` (int, MB) in `linear magnetic` / `linear thermal`;
  `SolverParameters` field + setter (validates ≥ 0) + `have_` + copy-ctor + `synchronize()` 13;
  explicit value → slot set at construction (O5); `doc/input_schema.yaml` + `doc/input_file_reference.md`;
  Codex prose pass on the touched section. (after: R6)
- [x] **R10 — Docs.** `solver_memory_and_compression.md` §7 rows (ICNTL(14) ceiling, ICNTL(23),
  INFOG copy 80); language sweep. (after: R7)
- [x] **R11 — Tests.** `tests/sparse/test_Solver.cpp`: `next_workspace_action` table (all six rows
  of §2; `-8`; `-19`). Plus R1/R2 tests. (after: R7)
- [ ] **R12 — Gate.** Rebuild (Christian); `make check`; `tape_quench` with `library : mumps`,
  `compression scheme : off`, `timestep { restart : false ; }` (cold start without deleting the
  dump — `cl_FEM_Controller.cpp:5172-5183`), `-v 4`. Score: step 1 no longer needs three failed
  factorizations; the later "240 exhausted" single-shot soft fails drop; per-rank RSS recorded.
  Conditioning off for the gate run (the shift-invert instance walks the same policy). (after: R11)

### 4.0 Implementation Progress (updated 2026-09-01)

**Implemented and audited (code jury 2026-09-01, Codex + Grok, both high):**
- Probe `available_memory()` — `MemAvailable` × 1024, cgroup v2/v1 headroom walked leaf→root,
  v1 controller lists split per token, no exceptions (`strtoull`), a declared-but-unreadable
  hierarchy reports *unknown* rather than the host figure; Darwin Mach branch uncompiled (owed).
- `allreduce_min`; `mumpstools_num_solvers()`; INFO/INFOG 80; `Parameter::MemoryBudget` → `ICNTL(23)`.
- **The budget is measured ONCE in `MUMPS::initialize()`**, after the create and before any
  factorization (`mMeasuredBudgetMB`) — a design change from the plan's "on the first `-9`": the
  guide defines `ICNTL(23)` as the *total* an instance may hold, and a probe after a failed
  factorization would have measured the remnant beside the workspace MUMPS still held (Grok P1,
  = O4). `escalate_workspace()` now contains no collective.
- Policy `mumps::next_workspace_action()` per §2: Cap → Ladder (ceiling 480) → GiveUp on `-19`,
  GiveUp when the measured budget is below the estimate, 1 MB floor for a measured-but-full machine.
- Diagnostics: soft-fail reason distinguishes "budget below estimate" / "exhausted" / plain; second
  row with estimate, cap, shortfall; every `int_t × 1e6` in `error_message()` replaced (nine arms).
- Deck key `memory budget` (whole MB, `> 0`, `≤ int_t max`, integrality enforced), copy-ctor,
  sync width 13, schema + reference + module doc §7 (retitled).
- Tests: probe (core, fast), policy table, parameter + boundary + copy, `allreduce_min` np 2/4.

**Defect tracker (code round):**
- [x] **D1 (P1, Grok).** Probe could throw. Fixed — `strtoull`.
- [x] **D2 (P1, Grok).** Input docs omitted the below-estimate give-up. Fixed.
- [x] **D3 (P1, Grok).** Soft-fail said "exhausted" on an untried ladder. Fixed.
- [x] **D4 (P1, Grok).** Probe after failure measured the remnant. Fixed — measure at `initialize()`.
- [x] **D5 (P2, Codex).** Nonstandard cgroup mount → host figure. Fixed — unresolved → unknown.
- [x] **D6 (P2, Codex).** Setter accepted values that narrow negative in `int_t`. Fixed + test.
- [x] **D7 (P3, Codex).** `get_int` rounded `1.6` to 2. Fixed — integrality required.
- [x] **D8 (P2, Grok).** v1 co-mounted controllers missed. Fixed.
- [x] **D9 (P2, Grok).** Over-limit sentinel lost in MB truncation. Fixed — 1 MB floor.
- [x] **D10 (P2, Grok).** Cap row 73 > 71 chars. Fixed.
- [x] **D11 (P2, Grok).** Six more `× 1e6` arms. Fixed.
- [x] **D12 (P2, Grok).** §7 heading. Fixed.
- [x] **D13 (P1, Christian's build).** Three `snprintf` box rows tripped `-Werror=format-truncation`
  at `-Og` — invisible to the `-fsyntax-only` gate used before the jury. Fixed: buffers sized for
  the worst case, prints clip to the field (`%-71.71s`, `%-53.53s`); every touched object and the
  shim now compiled for real (`-c -o /dev/null`, tree flags).
- [x] **D14 (P1, gate run 2026-09-01 15:2x).** Step 2 of the first capped run soft-failed on
  `-20` (MPI reception buffer too small; guide: raise `ICNTL(14)`), a code the policy did not
  recognise, so the controller cut Δt. Cause: with `ICNTL(23)` set the relaxation is applied to the
  buffers *first* and the cap squeezes the rest, so buffer shortages become likelier. Fixed: `-17`
  and `-20` are retryable — ladder only (a cap cannot widen a buffer), give-up at the ceiling;
  soft-fail decode names the buffer. Tests extended. Same-session docs.
- **FALSE POSITIVE (Grok, retracted by verification).** `ICNTL(35) = 3` bound: the guide gives
  `INFOG(16)` for `ICNTL(35) ∈ {0, 3}`; the code is right.

**Still owed:** R12 gate (build, `make check`, cold `tape_quench` under MUMPS with `-v 4`,
`restart : false`, conditioning off); Darwin compile; an executable `-9 → cap → JOB 5` test and a
`synchronize()` test for the new field (Codex P3 / Grok residual).

## 5. Open Design Questions

- **O1 — budget first or ladder first?** **RESOLVED 2026-09-01 (jury + Christian):** neither as a
  binary — cap first, ladder behind it for residual `-9`/`-8`, give up on `-19`.
- **O2 — Safety factor.** **RESOLVED 2026-09-01 (Christian, on recommendation):** 0.5 of
  cgroup-aware available, split by ranks on node and by live instances. A named constant.
- **O3 — Two live instances.** **RESOLVED 2026-09-01:** divide by the shim's occupancy via a new
  accessor (R3). Rank-uniform because create/free are collective.
- **O4 — Does MUMPS release S before returning `-9`?** **RESOLVED BY DESIGN 2026-09-01 (code
  jury, Grok):** moot for the measured path — the budget is now measured in `initialize()`, before
  MUMPS holds anything, which is also the quantity `ICNTL(23)` bounds (a total, not headroom).
- **O5 — Cap before the first factorization?** **RESOLVED 2026-09-01 (Christian, on
  recommendation):** yes for an explicit deck value, no for the probe (the "cap not fill" reading
  is from the upgrade notes, not from the allocator; nobody has read `dfac_mem_*`).
- **O6 — Ceiling.** **RESOLVED 2026-09-01:** 480. Supersedes `mumps_workspace_cap_raise.md` O1.

## 6. Interface Design

```cpp
// src/core/fn_available_memory.hpp
namespace belfem { size_t available_memory(); }          // bytes; 0 = unknown

// src/comm/commtools.hpp
template < typename T > void allreduce_min( const T * aSend, T * aRecv, const int aCount );

// src/sparse/mumpstools.hpp
extern "C" int_t mumpstools_num_solvers();                // live instances (occupancy)
enum class Parameter { ..., ComputeDeterminant, MemoryBudget /* ICNTL(23), MB, 0 = off */ };
enum class WorkspaceAction { Ladder, Cap, GiveUp };
WorkspaceAction next_workspace_action( int_t aInfoG1, int_t aRelax, int_t aBudgetSlot,
                                       int_t aBudgetMB, int_t aBoundMB );   // pure, tested

// src/sparse/cl_SolverMUMPS.hpp (private)
int_t memory_budget_mb();                                 // collective; 0 = unknown
```

Deck (R9):
```
linear magnetic
{
    library       : mumps ;
    memory budget : 4000 ;        // MB per process; optional; wins over the probe; set before the first factorization
}
```

## 7. Definition-of-Done Checklist

- [ ] Every gap-table row mapped to a step or an On.
- [ ] Citations re-read at implementation time.
- [ ] `allreduce_min` exercised at np 2 and 4 (`make check`, label `mpi`).
- [ ] Darwin branch of R1 compiled on a Mac before this plan closes (owed).
- [ ] Schema + reference updated together (R9); `scripts/check_doc_claims.py` run.
- [ ] R12 gate run by Christian, numbers in the devlog.
- [ ] `todo/closed/mumps_workspace_cap_raise.md` moved to `closed/` with a pointer here.

## 8. Audit Trail

- Code jury round 1 (2026-09-01, same exchange file): no P0; D1–D12 above; collective safety of the
  new path confirmed by both auditors independently, then simplified further by D4's fix.
- Plan jury round 1 (2026-09-01): `tmp/ai_exchange/review_mumps_memory_budget.md`. P0 (Codex +
  Grok, independently): give-up on residual `-9` contradicted the guide. Eleven P1 confirmed by
  Claude against the tree and the guide; none refuted. Full verification and reconciliation table in
  the exchange file.
- Prior plan: `todo/closed/mumps_workspace_cap_raise.md` (2026-08-30 jury). D2/D4 there are gap rows 9/6 here.
- Literature: MUMPS 5.9.0 user guide (`tmp/userguide_5.9.0.txt`): §2.11 (`:431-450`), §5.12,
  `ICNTL(23)` (`:2288-2324`), errors `-8`/`-9`/`-19` (`:5672-5679`, `:5707-5710`).
