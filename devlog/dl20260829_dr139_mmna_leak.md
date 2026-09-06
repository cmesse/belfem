# DR-139: `~ElectricalCircuit` mMNA Leak — Destructor Fix

**Date:** 2026-08-29
**Purpose:** Fix DR-139 (register-cleaning session): `~ElectricalCircuit` deleted `mJ` but never `mMNA`
**Module:** circuit

## What landed

One guarded `delete mMNA ;` in `~ElectricalCircuit`, mirroring the existing `mJ` arm
(`src/circuit/cl_ElectricalCircuit.cpp:71-74`). Nothing else in the source tree was touched.

The DR-124 precondition the register row demanded was re-checked before the edit: the sole
`mMNA` allocation is the `mMNA == nullptr` lazy-init in `compute_MNA_matrix()`, no other free
or assignment exists anywhere in the tree, and the destructor never nulls the pointer on a
live object — the one-shot guard cannot re-enter. Confidence: high, confirmed by both vendors.

## Audit round

Code audit per the standing three-vendor rule; pre-registration, both vendor verdicts and the
reconciliation are in `tmp/ai_exchange/dr139_mmna_leak.md`.

- **Codex:** C1 (leak) and C2 (guard preserved) confirmed high; C3 (no double-delete) refuted
  *as an unconditional invariant* — the class never deletes copy/move, so an implicit shallow
  copy would double-free. No copy exists in the tree; verdict approve with the qualification.
- **Grok:** do-not-refute on all three claims, with tightenings adopted verbatim: the leak was
  one `SpMatrix` per lifetime that *entered* `compute_MNA_matrix()` (a never-stamped circuit
  leaked nothing), and the guard is safe *because* the free is destructor-only, not because a
  dangling pointer would read as null. Also traced factory/test ownership as single-owner on
  every live path, and noted the pre-existing `mSolver`-after-`mJ` declaration order as
  SuperLU-relevant but not a DR-139 defect (wrapper aliases only `mJ`).
- **Convergent by-catch → DR-141 filed:** `ElectricalCircuit` and `Solver` are implicitly
  copyable while owning raw pointers (`cl_ElectricalCircuit.hpp:107-109`,
  `cl_Solver.hpp:49-55` with `delete mWrapper` in `cl_Solver.cpp:121-124`). Latent UB,
  unreachable today; the DR-139 line added `mMNA` to that blast radius but did not create it.

Status: **reviewed, not verified.** The row's valgrind/ASan gate is owed and needs a rebuild.

## Warning for whoever lands the circuit-cluster WIP

The circuit-cluster session ("Clean up circuit cluster debt issues") exited during this
session, leaving ~18 uncommitted files across `src/circuit/` and `tests/circuit/` — its
DR-138 fix and the memdump-v2 state-I/O work, which per `dl20260829_circuit_restart_v2.md`
is complete with its unit gates run green (full suite 15/15), not abandoned mid-flight.
This DR-139 hunk sits **inside that uncommitted set** in `cl_ElectricalCircuit.cpp` — two
unrelated changes now share the file. Stage explicit paths at commit time; do not attribute
the destructor hunk to the state-I/O work or vice versa.

## Register bookkeeping

- DR-139 retagged `[CODE]` → `[RUN]`, fix column filled, gate unchanged.
- DR-141 filed (implicit-copy double-free hazard, P3, `[CODE][P]`). The eigen/ARPACK session
  then found the identical shape in `fem::dofmgr::EigenValues` (four raw owned pointers, no
  deleted copy/move — verified here at `cl_FEM_DofMgr_EigenValues.cpp:46-66`); by agreement
  that sibling is split to DR-142, filed and fixed (`= delete` copy/move, syntax-gated) by
  that session the same day, and DR-141 carries the pointer. Ownership split: this session's
  row keeps `Solver`, that session's row keeps `EigenValues`.
- Header `[P]` count recounted mechanically: 26 (the DR-140 filing had left it one stale).
