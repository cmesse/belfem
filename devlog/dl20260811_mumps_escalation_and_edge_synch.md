# MUMPS out-of-range escalation (DR-45) and quadratic edge/face synch (DR-31)

**Date:** 2026-08-11
**Purpose:** Record two solver/postproc defects closed in source, the jury round behind the
first, and one cross-cutting finding about how both were missed for so long.
**Modules:** `src/sparse`, `src/fem/kernel`

---

## 1. DR-45 — MUMPS `INFO(1) = +1` is now a hard error

`+1` means MUMPS found indices in `IRN`/`JCN` outside the matrix, **dropped those entries** and
factorised anyway. The solve then succeeds on a different matrix than the one BELFEM assembled,
and nothing downstream can notice: the residual is measured against the assembled operator, so
no tolerance ever sees the missing entries.

Christian's reading was that this should be fatal. I agreed before any auditor was involved,
which is the specific reason a `--jury` round was run rather than a straight implementation —
two parties agreeing is not evidence. The pre-registration named the two findings that would
have reversed the decision, and both were declared as *wanted* results.

**Verdict (Codex and Grok, blind, independently): escalate.** Neither could construct a benign
producer of `+1`. The decisive evidence is structural rather than the absence of a grep hit:

- Fixed dofs are **not** discarded through a sentinel index. They live in a separate matrix —
  `mDirichletMatrix` is free × fixed, `mSystemMatrix` free × free
  (`cl_FEM_DofMgr_SolverData.cpp:350-405`). There is no discard idiom for the escalation to break.
- With `ICNTL( 5 ) = 0` / `ICNTL( 18 ) = 0` and `irn`/`jcn`/`A` bound only on the master
  (`mumpstools.f90:321,340-348`), `+1` is **host-local by construction** — a worker has no
  indices to be out of range.

Implemented in `78ea534d`: both `solve()` overloads route their positive-`INFO(1)` branch
through a new private `MUMPS::check_warnings()`, which raises `BELFEM_ERROR` on the bit and
reports the remaining warnings per rank. `free()` deliberately keeps its own report-only loop —
it is teardown, reached from `~MUMPS()` after `JOB = -2`, where aborting would mask whatever
caused the teardown.

### The scope guard, which is the real hazard

Both auditors raised it independently: the check tests the **bit** (`mInfo( 0 ) & 1`), never
`mInfo( 0 ) > 0`. MUMPS returns warnings as a sum of flags, and `+8` (iterative refinement did
not converge) is *expected* on ill-conditioned HTS systems. Aborting on any positive `INFO(1)`
would have reproduced the failure from earlier the same day, where an always-active guard on a
normal distributed state aborted every parallel run (`61bbfa11`).

`BELFEM_ERROR` rather than `BELFEM_ASSERT`: the latter compiles out under `NDEBUG`, which would
leave the release build quietly solving the wrong matrix. `+1` is never routed through
`soft_fail()` either — a retry re-factorises the same wrong matrix.

### What the round corrected in my own reasoning

**A premise, refuted.** I argued `+1` was effectively silent because the message is
`InfoLevel::Minimal`. That is backwards: `Minimal = 1`, `Default = 2`, and the gate is
`aInfoLevel <= mInfoLevel` (`cl_Logger.hpp:33-38,102`), so `Minimal` is the tier *most* likely
to print. Only `Silent` would have hidden it. What survives is rank gating, not log level — and
under `ICNTL(18) = 0` the host is the only rank that can raise it, so an absorbed `+1` would have
been printing on rank 0 all along and read by nobody. **Visible but unread, not silent.**

**A disagreement between auditors, settled by source.** Codex wanted an `MPI_Allreduce`
rank-uniform verdict before raising, citing `cl_SolverSTRUMPACK.cpp:339-347`. That precedent is
the *soft-fail* path, which returns into further collectives and must stay uniform; a hard error
terminates via `MPI_Abort` and needs no reduction. Adding a collective per Newton iteration to
decide whether to abort buys nothing.

**A defect I had shipped that morning.** `5a2ddf81` claimed to fix "both `INFO` reads". There are
**three** sites — `MUMPS::free` and *both* `solve` overloads — and the Matrix-RHS one was still
`rank() == 0` gated. Raised by Grok in the round, confirmed by direct read, closed in the same
commit as the escalation.

---

## 2. DR-31 — quadratic edge and face fields were synched with the wrong stride

The row had sat at `unverified` with the note "needs a layout check before calling it a bug".
The check was done end to end and **the suspicion was correct**.

`DofManager::create_fields` sizes an edge field as
`edge_multiplicity() * number_of_edges()` (`cl_FEM_DofManager.cpp:215`), and `IWG_Maxwell` sets
that multiplicity to **2** whenever the element order is 2 (`cl_IWG_Maxwell.cpp:63-73`), storing
the two dofs of edge *i* at `2i` and `2i+1` — exactly how
`Calculator::nedelec_data_quadratic_2d/3d` reads them back.

`Postprocessor::synch_source_field` copied **one** value per entity index, and
`mAllEdgeIndices`/`mMyEdgeIndices` hold *entity* indices, not dof positions. For a quadratic
field it therefore read slot *i* out of a `2N`-long array: the second dof of every edge was
never transferred at all, and the value that *was* transferred belonged to a different edge.
`face_h` (multiplicity 2 over `number_of_faces()`) had the identical defect.

**It does not crash.** The indices stay in bounds, so this is silent wrong data on aura
entities — wrong recovered H/B/J near partition boundaries, in a run that otherwise looks
converged.

**Reachability is not opt-in.** `mElementOrder = mMesh->max_element_order()`
(`cl_MaxwellFactory.cpp:87`) and the flag is `mElementOrder == 2` (`:560`), so loading a
second-order mesh and running on more than one rank is enough. `MaxwellPostprocessor` puts
`edge_h`/`face_h` in the source list whenever `face_h` exists
(`cl_MaxwellPostprocessor.cpp:99-106,131-138`).

Fixed in `0671b29a`: the stride is recovered from the field's own length rather than assumed,
and `multiplicity` values are copied per index on both the send and the receive side. It is
derived from the field rather than queried from the IWG so the two cannot drift apart; asserts
cover the length dividing evenly and the received count matching.

Severity raised P2 → P1 on the wrong-answer class. The blocking-1.0 flag is left open — it turns
on whether order-2 meshes are in scope for the release, which is not a call this work can make.

---

## 3. The cross-cutting finding: references confirmed by recognition

Three distinct instances surfaced in one day, and they are the same defect:

| pointer | said | was |
|---|---|---|
| `cl_SolverMUMPS.cpp:238` / `:408` / `:563` | current line anchors | correct against `5a2ddf81`, stale the moment `78ea534d` rewrote the region |
| "DR-24" as the label for the parallel-run abort | that incident | DR-24 is the periodic free-cut `mPrescribedCurrents` row, closed 2026-08-09; the incident has **no row** — cite `61bbfa11` |
| "both `INFO` reads (`MUMPS::free`, `MUMPS::solve`)" | an exhaustive enumeration | three sites; `solve` is overloaded |

In every case the pointer was **correct when written** and was never re-followed. None is a
reasoning error; none would have been caught by thinking harder. The DR-24 mislabel survived a
jury round, a row edit and a second reader, because a plausible ID reads right at a glance and
confirming it costs a lookup nobody spends.

The convention already in force — anchor by a searchable token, never a line number — was
written about `file:line` in the input contract. A register ID, a commit hash and a prose
enumeration are the same object: a reference confirmed by recognition rather than lookup. The
line number is simply the instance that rots fastest. Whether to widen the written rule is
Christian's call and is recorded as such.

---

## 4. Open gates

Neither fix has executed. Both are **reviewed, not verified** — each compiles clean under
`-DNDEBUG` and `-DDEBUG` with `-Werror`, and that is all.

- **DR-45:** helix serial and on ≥ 2 ranks, completing *without* the new hard error firing. The
  regression risk is a false abort on a healthy deck, not a missed `+1`.
- **DR-31:** an order-2 mesh on ≥ 2 ranks, recovered H/B compared against the serial run. Note
  the serial path never calls `synch_source_field`, so a serial run proves nothing here.

The register row for DR-45 was struck on Christian's instruction; struck is not verified, and
the reproducer column carries the gate.
