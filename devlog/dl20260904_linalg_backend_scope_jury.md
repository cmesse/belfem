# Backend Verification Plan: Scope Jury, Probe Declined, `USE_PARDISO` Preset

**Date:** 2026-09-04
**Purpose:** Record the jury round that settled the scope of `todo/linalg_backend_verification.md`,
Christian's ruling, and the edits that closed the plan the same day.
**Module:** `config/system`, `config/linalg`, `CMakeLists.txt` (build system only; no `src/` code)
**AIs involved:** Claude (pre-registration, verification, edits), Codex `gpt-5.6-terra`/`high` and
Grok `grok-4.6`/`high` (blind jury)

---

## 1. The question

Christian, reading the todo list: the idea is simple — if `$SCLS` is set, read the flavor; if it
is `mkl`, default `USE_MKL` and `USE_PARDISO` to ON, otherwise not. Is the plan's configure-time
link probe (R2–R5, with a `belfem_probe_link()` helper) the right path at all?

## 2. The round

Pre-registration, both audits and the verification pass are in
`tmp/ai_exchange/review_linalg_backend_scope.md`. Every `file:line` was re-read; no citation
failed, two lacked a directory component (`config/compiler/finalize_compiler.cmake`) and one was
off by one (`USE_OPENMP` at `CMakeLists.txt:78`).

All three voices drop the helper. Grok and Claude drop the probe entirely; Codex would keep one
stock `CheckCXXSourceCompiles` placed after `finalize_compiler.cmake`. Confirmed plan defects:

- The R7 row `debug` + `USE_MKL=ON` expected a `FATAL_ERROR` that R3's symbol probe could not
  produce — `$MKLROOT` is valid on this host, so `dgemm_` links. The row needed a
  flavor-versus-knob check, not a probe (Grok).
- D3's `cblas_dgemm_batch` requirement, applied to the netlib branch, would have refused the
  `debug` flavor, whose SLATE has no `cblas_` references (Grok). BELFEM never calls the symbol;
  only the `mkl` flavor's SLATE build has it undefined, and MKL provides it.
- R6 could not name which BLAS won the link under netlib — nothing resolves a path — while O4
  leaned on it as the mitigation (Codex).
- R5 offered `-DBLAS_DIR=…`, which R4 forbids (both).
- `find_scls_flavor.cmake:35-37` promised a "backend check in config_mkl.cmake" that did not
  exist (Grok).
- `LINK_OPTIONS` on `try_compile` is 3.14, not 3.25 as D5 said; the conclusion held on the 3.11
  floor (Codex).

Write-safety: two `config/doc/*.cmake` files changed during the jury window. Grok's
`sandbox-events.jsonl` shows an enforced landlock read-only profile and the Codex transcript never
names them; a second Claude session and CLion were open on the checkout. Not an auditor.

## 3. Christian's ruling

Presets only. `USE_PARDISO` defaults in parity with `USE_MKL`. No link probe, no library report:
every vendor other than MKL is assumed to behave like reference BLAS/LAPACK, and the group does
not need to know which library won. D3 is not needed. Remove the stale comment, strike the stale
plan lines, close the plan.

A first cut also added a flavor-versus-knob `FATAL_ERROR` at the top of `config_mkl.cmake`
(`mkl` + `USE_MKL=OFF`, or a netlib flavor + `USE_MKL=ON`), which the jury's R7 finding had
suggested. Christian tested it and withdrew it the same hour: the knobs are **defaults**, nothing
more. A flavor whose name contains `mkl` presets `USE_MKL` and `USE_PARDISO` ON; every other
flavor, and no `$SCLS`, presets them OFF; the user overrides with `-D` or in ccmake and is never
refused. The unrecognized-flavor notice went with it — a flavor is just a name.

## 4. What landed

| file | change |
|---|---|
| `config/system/find_scls_flavor.cmake` | `BELFEM_DEFAULT_USE_PARDISO` in parity with `USE_MKL`; the match is `MATCHES "mkl"` on the flavor name; comment table rewritten; the "backend check catches it" sentence and the unrecognized-flavor branch removed |
| `CMakeLists.txt` | `option( USE_PARDISO … ${BELFEM_DEFAULT_USE_PARDISO} )`; header comment names all three knobs and that `-D` wins |
| `config/linalg/config_mkl.cmake` | unchanged in the end — the refusal block was added and removed within the session |
| `CLAUDE.md` | flavor-default paragraph states the name rule, the PARDISO parity and that nothing is refused; `check_doc_claims.py` 38/38 |
| `src/sparse/doc/sparse_usage_guide.md` | `USE_PARDISO` default line |
| `todo/linalg_backend_verification.md` → `todo/closed/` | banner, Status, R0/R0b ticked, R2–R7 and D3/D5/D9–D11 struck with dated reasons, §7 boxes struck |
| `todo/README.md` | active count 5 → 4, closure sentence |

`cmake -P` on the flavor file alone: `mkl` presets ON/ON/OFF, `debug` presets OFF/OFF/ON, `cea`
and `gcc` preset all OFF silently, and a `SCLS_FLAVOR=mkl` over `/opt/scls/debug` warns and
trusts the variable. Christian ran a real configure against the first cut, which is what
surfaced the refusal as unwanted.

Existing trees are untouched: `option()` never overwrites a cached `USE_PARDISO:BOOL=OFF`, so
both current trees keep PARDISO off until reconfigured from an empty directory. On a fresh `mkl`
tree the change compiles `pardisotools.f90` and `PARDISOSolveTridiagonal` appears in the suite;
`gDefaultSolver` stays STRUMPACK (`en_SolverEnums.hpp:196-208`).

## 5. Owed

- One fresh configure per flavor to see the presets land. Christian's to run.
- A fresh `mkl` configure and build to confirm the PARDISO path still compiles (INC-527 rot
  pattern).

## 6. Same session: `USE_BELFEM_OPENMP` documented in `src/sparse/doc`

Christian asked under what circumstances `USE_BELFEM_OPENMP` should be ON. The answer lived only
in the 2026-09-01 devlog, `doc/parallel_execution.md` and the kernel comment, and the sparse
module's own docs did not mention the switch at all. Added a subsection "BELFEM's Own OpenMP
Kernels: `USE_BELFEM_OPENMP`" under "Solver Parallelism" in `sparse_usage_guide.md`, and a
paragraph plus a contents line in `src/sparse/doc/README.md`. Content: the two-switch split
(`USE_OPENMP` for the solvers, `USE_BELFEM_OPENMP` for the three Fortran kernels), the stack
reduction mechanism and the 262,144-row limit, why the kernels are master-only noise, that an MKL
build never reaches them, the two legitimate ON cases (a per-phase measurement with
`OMP_STACKSIZE` raised, and a toggle proof via `make check`), and what ON does not change. Codex
`luna`/`medium` language sweep run over both pieces (`tmp/ai_exchange/sweep_sparse_belfem_openmp.md`).

## 7. Same session: mechanical US-spelling sweep of the user-facing docs

Christian's standing rule, restated today: BELFEM uses American spelling throughout (the
`Neighbour` → `Neighbor` identifier sweep was his). A Codex prose sweep had returned US spellings
that I was about to revert; corrected, and the leftovers in the documentation were then swept
mechanically. Scope: `doc/`, every `src/**/doc/`, `CLAUDE.md`, `README.md`, `AGENTS.md`,
`METHODOLOGY.md`. About 200 tokens in 45 files, prose only. Guards that mattered:

- `doc/literature_references.md` excluded, and any capitalized token on a line carrying a year or
  "et al" skipped: journal and paper titles (*Mathematical Modelling and Numerical Analysis*,
  "Modelling the E–J relation …") keep their spelling.
- `aluminium` never touched: it is an accepted material-label alias in `cl_MaterialFactory.cpp`
  and appears as such in `doc/input_schema.yaml` enum lists.
- Tokens inside backticks or fenced code that also occur in code were left alone; path segments
  (`labelled_subsections`, `default_behaviour`, `behaviour_change` — schema keys) untouched.
- `cancelled` skipped throughout: `todo/cancelled/` is a directory name.

`check_doc_claims.py` 38/38 after the sweep. Not swept, reported with counts: `devlog/` (~200
hits) and `todo/` (~340 hits), and two runtime message strings in `src/` that still say
"labelled" (`cl_MaterialFactory.cpp:465-466`, `cl_MaxwellFactory.cpp:1075`) — source edits, so
not part of a docs sweep.
