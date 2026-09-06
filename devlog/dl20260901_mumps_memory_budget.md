# MUMPS memory budget: ICNTL(23) from a machine probe, ladder behind it

**Date:** 2026-09-01
**Purpose:** Session record — from "the intrinsic-material tape_quench deck won't launch under
MUMPS" to a jury-audited change that hands MUMPS a measured per-process memory cap.
**Module:** `src/sparse`, `src/core`, `src/comm`, `doc/input_*`

## Trigger

`cmake-build-debug/tape_quench` (builtin materials, B-spline jc/n tables) first crashed in
`ThinShellFactory::create()` — the layer stack named `rebco` while the materials section had been
renamed `ybco`; the map lookup dereferenced `end()` (fault address 0x28 = `second` in the
`unordered_map` node). Deck fixed; a `BELFEM_ERROR` guard now names the layer and the stack
(`cl_ThinShellFactory.cpp`, validated over the whole stack, outside the ghost-facet branch).

The same deck under `library : mumps` then could not cold-start: `-9` at every rung of the
`ICNTL(14)` ladder (30 → 60 → 120 → 240), every Δt cut down to 0.031 ms failing identically. The
sibling `tape_quench_usermat` directory "fired up" only because it warm-starts from a `memdump.hdf5`
at t = 22 ms (`belfem.cpp:207` loads it unconditionally); its runs 2–7 under MUMPS carried 129 `-9`
events rescued by Δt cuts, and its only clean long run (run 8, 1103 steps) is the one after the
deck moved to STRUMPACK. Christian's STRUMPACK relaunch of the intrinsic deck reached t = 20 ms in
56 steps with Δt at the 1 ms deck maximum and zero solver trouble — not an A/B (solver and
materials both changed), recorded as such.

## What landed

Plan `todo/mumps_memory_budget_plan.md` (jury round 1: P0 — my give-up-on-residual-`-9` policy
contradicted the MUMPS 5.9 guide, which says `-9` "may still occur" with `ICNTL(23)` set and still
wants a larger `ICNTL(14)`; `-19` is the give-up code). Code, jury round 2, fixes applied:

- `belfem::available_memory()` (`src/core/fn_available_memory`): `MemAvailable` × 1024 capped by
  the tightest cgroup limit (v2 `memory.max`/`memory.current`, v1 with controller lists split per
  token, walked leaf→root at the standard mount; declared-but-unreadable → unknown), no exceptions;
  Darwin Mach branch written, not compiled.
- `belfem::allreduce_min`; `mumpstools_num_solvers()`; INFO/INFOG copy 40 → 80 (library width;
  `INFOG(36)` is the BLR estimate); `Parameter::MemoryBudget` → guarded `ICNTL(23)` write.
- Budget measured **once in `MUMPS::initialize()`**, before any factorization (the guide defines
  `ICNTL(23)` as the total an instance may hold; a probe after a failed factorization measures the
  remnant — Grok's P1, which also closed the plan's O4): available / ranks on node / live instances
  × 0.5, MIN over ranks, 1 MB floor for a measured-but-full machine.
- Policy `mumps::next_workspace_action()` (pure, table-tested): on `-9`/`-8` with the slot empty,
  cap if the budget covers MUMPS's estimate, give up if it does not; behind a cap, the ladder to a
  new ceiling of 480; give up on `-19`. Soft-fail box now says which of the three happened and
  prints estimate, cap and shortfall; all nine `int_t × 1e6` overflows in `error_message()` gone.
- Deck key `memory budget : <whole MB> ;` (MUMPS only, applies from the first factorization,
  wins over the probe); schema, reference, module doc §7 updated in the same session.
- Tests: probe (core, fast), 13-row policy table, parameter/boundary/copy, `allreduce_min` np 2/4.

Syntax-gated with the tree's `-Wall -Werror -pedantic-errors` flags (C++ and gfortran); not built,
not run — Christian builds. **Gate owed (R12):** `make check`; cold `tape_quench` under MUMPS with
`compression scheme : off`, `timestep { restart : false ; }`, `-v 4`, conditioning off.

## Also this session

- `todo/mumps_workspace_cap_raise.md` (2026-08-30) superseded; its O1 resolved at 480.
- Monitor tooling for the running deck (`scratchpad/watch_tape_quench.sh`): milestones,
  relaunches, rate-limited trouble, exit — plus one false-alarm fix (it re-announced unchanged
  counts).
- Memory correction: both build trees read `USE_DEBUG=ON` today; the 2026-08-30 inversion note was
  stale. Read `CMakeCache.txt`, never the name.
- Diagnostics regression noted: the Sep 1 binary silences MUMPS's raw `INFOG(2)` lines below
  `-v 4`; the soft-fail box now carries the number instead.

Exchange: `tmp/ai_exchange/review_mumps_memory_budget.md` (plan pre-registration, two auditor
entries per round, verification and reconciliation tables).
