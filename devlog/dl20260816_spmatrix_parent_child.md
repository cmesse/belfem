# Devlog 2026-08-16 — SpMatrix Parent/Child Structure Sharing

**Date:** 2026-08-16
**Topic:** Memory-sharing parent/child mechanism for SpMatrix pairs with identical sparsity; three cross-review rounds; int_t consistency campaign through the solver layer
**AIs involved:** Claude (primary), Codex + Grok (jury audits, 3 rounds), Codex (prose pass)
**Claude Confidence:** high on the landed mechanics; the two open questions are DR-80/DR-81
**Codex Audit Confidence:** high (round 3: "no blocking runtime defect")
**Grok Audit Confidence:** high on ownership/indexing/COO control flow; explicitly NOT signing off the linked-load path as verified until it has a test
**Literature References:** N/A (ownership/lifetime, not a formulation)
**Verification:** **focused regression — `make check` on Christian's machine,
`8/14 Test #8: sparse ... Passed 0.58 sec`**, covering the full 34-test parent/child +
HDF5 battery in `tests/sparse/test_SpMatrix.cpp`. That is protocol §11 rung 2, so the
mechanism is **verified**, not merely reviewed. (Prior to that gate the evidence was a
`mpicxx -fsyntax-only` pass with the real build flags on every touched TU.)
Not covered by the gate: the `transpose()` defect below (no test can currently fail on
it) and linked `multiply`, which has no test yet.
**Correction (round 4):** an earlier draft of this entry named `make check-fast` as the
gate. That is wrong and the auditors caught it — `check-fast` runs `ctest -L fast`
(`CMakeLists.txt:313`), `tests/sparse/CMakeLists.txt` sets no `TESTLABELS`, and `sparse`
is absent from the `check-fast` dependency list (`CMakeLists.txt:320-325`). The sparse
suite has never been part of the fast gate; making it so would need
`set( TESTLABELS fast )` plus `sparse` on that list — a deliberate change, not made here.

## Summary

Christian's idea: `mSystemMatrix`/`mJacobianMatrix` and `mFullMassMatrix`/`mFullStiffnessMatrix`
in `SolverData` are built from identical graphs, so each pair can share its pointer/index
(and MUMPS COO) arrays. Saving ≈ `(nnz + n + 1) * sizeof(int_t)` per pair. Verdict of all
three reviewers across three rounds: the idea is sound; the initial implementation had 7 P0s
(inverted HDF5 check, stale child index function on parent-initiated base flips, COO
fall-through/dangling aliases, half-initialized child ctor, missing lifecycle unlinking,
unsafe move/copy); all were fixed via the frozen plan+audit → code+audit method. The final
jury found no blocking defect.

## Design as landed (src/sparse/cl_SpMatrix.{hpp,cpp})

- `explicit SpMatrix( SpMatrix * aParent )`: guards (null, chain, second child, empty
  parent), copies structure via `update_from_parent()`, owns a zero-filled value array,
  inherits the indexing base, auto-links. Copy and move CONSTRUCTION are `= delete`
  (Rule-of-5 closure); deep copies go through assignment onto an empty matrix.
- Either destruction order is safe: child-first unlinks the parent's `mChild`;
  parent-first transfers structure ownership to the child (the `SolverData` teardown
  order is parent-first for both pairs).
- `set_indexing_base` propagates up (guarded) and down (UNCONDITIONAL — the downward
  guard would read the already-converted shared `mPointers[0]` and always skip).
- COO: child always re-aliases from the parent in both creation orders; frees null the
  child's alias in both orders (`free_coo_indices` propagates down unconditionally);
  `deallocate()` now clears `mHaveCooIndices`.
- `load()` on a linked matrix (either end) is verify-only: format + dims + full
  base-normalized pointer/index comparison against the existing arrays, values-only
  replace, no sort, no reallocation — keeps `SolverData::load_system` working on a
  linked parent. Unlinked load keeps the old rebuild+sort path. Scratch buffers are
  freed before any `BELFEM_ERROR` (test-throw mode safe); mallocs and HDF5 statuses
  checked. `sort_entries`/`set_type`/`transpose`/copy-onto/move-from-or-onto refuse
  linked matrices (`BELFEM_ERROR`, setup tier).
- `memory()`: values always owned; pointer + structural + optional COO arrays charged to
  the parent only (fixed a double-count that also dropped the pointer term).
- The dead friend `distribute( SpMatrix*, proc_t )` was deleted — a whole SpMatrix is
  never MPI-sent; MUMPS is host-centralized, STRUMPACK/PETSc scatter raw slices via
  `sparse::DistMatrix` (rank 0 only). DistMatrix untouched, per Christian.

## int_t consistency campaign (Christian-initiated, same session)

- `SpMatrix::indexing_base()` returns `int_t` (was `int`, silently narrowing
  `mPointers[0]` on INT64 builds); `tOldBase` locals widened.
- Christian's PARDISO change (`Vector<int_t>` mParameters/mInfo) confirmed correct — the
  `pardisotools` Fortran shim already takes/returns `integer(int_t)`; the old
  `Vector<int>` was a latent INT64 ABI bug. Follow-through fixes: the `initialize`
  vtable unified to `int_t` across Wrapper + all six solver overriders (MUMPS already
  had it; the split would not compile under INT64), `pardisotools.hpp` return types
  corrected from `int` to `belfem::int_t` (Fortran returns `integer(int_t)`), and the
  four PARDISO `%i` format sites converted to `%ld` + `(long)` casts (vararg UB under
  INT64 — Grok round 3).

## Tooling fixed on the side

- `scripts/scls_env.sh`: Linux-only pin now early-outs on any other OS via `uname -s`
  (was aborting every non-interactive script on macOS).
- `.claude/scripts/ask_codex.sh` / `ask_grok.sh`: BSD `mktemp` does not substitute
  `XXXXXX` before a suffix — the old templates created LITERAL files and concurrent
  auditor runs collided. Templates now end in `XXXXXX` and honor `TMPDIR`.
  (Lesson recorded: never edit a wrapper while a jury is running — bash reads scripts
  incrementally and the running leg died on a phantom syntax error.)

## Docs

`src/sparse/doc/sparse_usage_guide.md`: new "Parent/Child Structure Sharing" section;
the three fictional `distribute(&A, 0)` examples and the "Replicated/Distributed"
pattern list replaced with the host-centralized reality; copy-construction examples
rewritten to assignment form. README pointer added. Codex prose pass applied (8 edits).

## Round 4 — DR-80 settled, DR-81 tests written

**DR-80 REFUTED by source read (closed).** The three SuperLU routines the wrapper calls are
all read-only on the input arrays: `get_perm_c()` writes only `perm_c` and builds A'A /
A'+A in its own `SUPERLU_MALLOC`'d arrays; `sp_preorder()` never writes
`Astore->{nzval,rowind,colptr}` and forms AC by ALIASING `nzval`/`rowind` while allocating
only its own `colbeg`/`colend`; `dgstrf()` only scatters out of them
(`dense[asub[k]] = a[k]`), with no equilibration (that lives in the `dgssvx` driver, which
BELFEM does not call). Two facts kept: AC aliases the VALUE array, so a `SamePattern`
refactorization sees re-assembled values automatically; and `SUPERLU::solve` re-runs
`symbolic()` when `&aMatrix != mMatrix` (`cl_SolverSUPERLU.cpp:48-51`), so alternating
system/Jacobian solves re-wrap rather than reuse a stale alias.
*Note on provenance:* `tmp/superlu` was a saved GitHub landing PAGE (HTML), not the source
tree — the routines were read from the upstream raw sources instead.

**DR-81 CLOSED — the suite ran green** (`sparse ... Passed 0.58 sec`). Incidental datum for
the open `check-fast` question: the entire sparse binary, `test_Solver.cpp` and its
MUMPS/PARDISO initialization included, takes 0.58 s — comfortably inside the "<~10 s per
test" criterion, so labeling the suite `fast` is now an evidence-backed option rather than
a guess.

**DR-81: 10 HDF5 tests written** (`§6` in `tests/sparse/test_SpMatrix.cpp`): unlinked
round-trip (regression anchor for the branch that was restructured), child-loads-matching-
file, linked-parent values-only load (the `load_system` mechanism), cross-base load,
four rejections (index / pointer / dimension / format mismatch, each isolated by
same-nnz patterns), save-from-child, orphaned-owner rebuild. The register row stays OPEN
until the suite actually runs.

**One source bug found by the plan audit (Codex) and fixed:** the linked-load index scratch
guard `BELFEM_ERROR( tSwap != nullptr )` would fire spuriously for a zero-nnz matrix, since
`malloc(0)` may return `nullptr` — and zero nnz is treated as valid elsewhere in the class
(`create_coo_indices`, `allocate_values`). Both scratch guards now carry the
`|| size == 0` escape.

**Plan-audit corrections folded into the tests** (all verified against the tree before
applying): T10 must heap-allocate or the parent cannot die first; T4 must save from a
SEPARATE unlinked Fortran-base matrix, because flipping the pair converts the shared arrays
and leaves `tShift == 0`; T1/T9 need elementwise structure compares, not pointer identity;
file names must not contain `/`; the mismatch builders must place actual nonzeros or they
collapse into the dimension test. The T8 justification in the plan was wrong (the guard is
`tType == mType`, and `tComp` follows `mType`, so nothing null-derefs) — the test is kept
for the right reason: the tridiagonal matrix is symmetric, so its CSR and CSC dumps carry
identical arrays and the format string is the only discriminator.

## Open Questions

- **DR-81** stays open: the tests exist and are syntax-gated, but have not been run.
- The sparse suite is outside `make check-fast` (see the Verification correction above) —
  worth deciding whether it should be labeled `fast`; the tests are all tiny, but
  `test_Solver.cpp` pulls in MUMPS/PARDISO/PETSc initialization.
- `indexing_base()` callers outside sparse (PARDISO `mParameters(1)` cast) reviewed;
  EigenValues `mOriginalBase` was already `int_t`.

## Files Updated

- src/sparse/cl_SpMatrix.hpp / cl_SpMatrix.cpp (mechanism + all fixes)
- src/fem/kernel/cl_FEM_DofMgr_SolverData.cpp (two pairing lines by Christian; comment fixes)
- src/sparse/cl_SolverWrapper.{hpp,cpp}, cl_SolverSTRUMPACK.{hpp,cpp},
  cl_SolverUMFPACK.{hpp,cpp}, cl_SolverSUPERLU.{hpp,cpp}, cl_SolverPETSC.{hpp,cpp}
  (initialize → int_t), cl_SolverPARDISO.{hpp,cpp} (Christian's int_t members + format fixes),
  pardisotools.hpp (return types)
- tests/sparse/test_SpMatrix.cpp (14 parent/child tests)
- src/sparse/doc/sparse_usage_guide.md, src/sparse/doc/README.md
- scripts/scls_env.sh, .claude/scripts/ask_codex.sh, .claude/scripts/ask_grok.sh
- todo/debt_register.md (DR-80, DR-81)

Exchange thread: `tmp/ai_exchange/review_spmatrix_parent_child.md` (3 frozen rounds +
reconciliations). Plan: `tmp/ai_exchange/spmatrix_fix_plan.md` (+ amendment A1).
