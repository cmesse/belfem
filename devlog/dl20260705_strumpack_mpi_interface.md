# STRUMPACK MPI Interface Review + Low-Risk Fix Batch

**Date:** 2026-07-05
**Purpose:** Fresh-eyes review of the STRUMPACK interface in `src/sparse` (MPI clunkiness), Codex-audited, followed by the approved low-risk fix batch.
**Module:** sparse

> **REVERTED same day, RE-LANDED same evening with S1 default flipped.**
> With MC64 matching off (S1), the serial coupled h-φ run failed with
> ZERO_PIVOT at Timestep 1 / Newton 17; Christian reverted the whole batch.
> Root cause (verified in STRUMPACK 7.2.0 sources): `FrontDense::factor_phase2`
> runs getrf FIRST and records ZERO_PIVOT on *exact* singularity —
> `replace_tiny_pivots` patches the diagonal only afterwards and never
> clears the error (it handles tiny pivots, not exact zeros). Equilibration
> always runs for unsymmetric matrices, so MC64's **permutation** is the
> load-bearing part, not its scaling; no cheaper substitute exists (all
> matching jobs except CombBLAS gather to rank 0; SCLS builds lack
> CombBLAS). Cost re-framed: the gather is once per initialize, not
> per Newton step. The batch was re-landed from `./tmp/backup/` with
> `mUseMatrixMatching = true` (STRUMPACK default preserved; `matching : off`
> is a verified opt-out), comments/error hint updated. S2/S3/S5/S6
> unchanged; S4 still open. Codex's "needs-measurement" grade on S1 was
> the correct call — recorded in memory as load-bearing.

## Context

Christian asked for a fresh look at the STRUMPACK wrapper (works, but "clunky with MPI"),
with the STRUMPACK 7.2.0 sources available in `tmp/STRUMPACK` for cross-checking.
Claude reviewed read-only, Codex audited all findings (full thread:
`tmp/ai_exchange/strumpack_mpi_interface.md`, ephemeral).

**Architecture verdict:** the hand-rolled `DistMatrix` + `set_distributed_csr_matrix`
/ `update_matrix_values` design is correct and justified — STRUMPACK's
`broadcast_csr_matrix` (root-only input, internal Scatterv) resets `reordered_`,
so it would lose symbolic-factorization reuse across Newton steps.

## Findings (all Codex-confirmed)

- **S1 — MC64 matching silently on.** STRUMPACK defaults to
  `MatchingJob::MAX_DIAGONAL_PRODUCT_SCALING`; in the MPIDist solver this
  *gathers the entire matrix to rank 0*, runs sequential MC64, and broadcasts
  (`CSRMatrixMPI.cpp:733-790`; STRUMPACK's own docs recommend disabling).
  BELFEM never called `set_matching`. This was the largest hidden serial
  bottleneck + rank-0 memory spike in the distributed setup.
- **S2 — CompressionMethod parameter was decorative.** The N-based BLR
  heuristic in the wrapper unconditionally overwrote the user's OFF/BLR choice
  (explicit OFF re-enabled at N≥50k; explicit BLR disabled below). Only the
  command-line override survived. Identical 20-line block duplicated
  serial/parallel.
- **S3 — initial guess ignored on the MPI path.** Serial passed the caller's
  `aLHS`; parallel passed `mMyLhs`, never filled from `aLHS`
  (`DistMatrix::distribute_lhs` was called by PETSc but not STRUMPACK).
  STRUMPACK semantics verified: `use_initial_guess` matters for iterative
  Krylov only; DIRECT ignores it (harmless).
- **S4 — per-solve value scatter stages a persistent full-nnz copy on rank 0**
  (`mAllValues`), although per-rank CSR value slices are contiguous and
  `mOffsets` was already computed but never used. The zero-staging
  `distribute(tData, mOffsets)` pairs correctly with `receive(mMyValues)`
  (Codex-verified chunk/tag pairing). **Deferred** — graded medium-high
  (communication plumbing, wants a 2/4/8-rank smoke test); prepared as a
  follow-up.
- **S5 — closed by decision (Christian):** BELFEM-side metis_ndp
  pre-permutation was an experiment; the solver library should handle
  reordering itself. Codex confirmed no checked-in call path activates it
  (SolverParameters default AUTOMATIC), but the input key `reordering scheme`
  could — hence the guard below. Not to be confused with mesh partitioning or
  `KernelParameters`' DOF-side METIS default.
- **S6 — minors:** `reinterpret_cast` downcasts, `hatch_small_turtle`
  threshold/typo inconsistencies, `StrumpackSparseSolver` now a plain alias of
  `SparseSolver` in 7.x, redundant barriers (left alone), and the 7.1.2
  upstream "fix hang for small problems" note (check what SCLS links).

## Changes applied (user-approved batch)

- `en_SolverEnums.hpp/.cpp` — new `CompressionMethod::AUTOMATIC` (+ `to_string`
  "automatic"; parse comes free via the loop-based `compression_method()`).
- `cl_SolverParameters.hpp/.cpp` — compression default BLR → **AUTOMATIC**;
  new `mUseMatrixMatching` (default **true** after the same-day revert — see
  header note) with `set_matrix_matching`/`use_matrix_matching`, input key
  `matching` (bool, solver section), copy-ctor + `synchronize()` extended
  (6 → 7 uints).
- `strumpacktools.hpp/.cpp` — `set_strumpack_options` gains `aNumRows`;
  sets `MatchingJob::NONE` only when `matching : off` is requested (S1,
  default keeps STRUMPACK's MC64); the N-based BLR heuristic now lives in
  the `AUTOMATIC` case only, explicit OFF/BLR honored (S2); ZERO_PIVOT
  message points at the matching opt-out as the likely cause.
- `cl_SolverSTRUMPACK.cpp` — both duplicated heuristic blocks removed (options
  configured via the shared helper); parallel solve distributes the caller's
  LHS when `use_initial_guess()` (S3); 4× `reinterpret_cast` → `static_cast`
  (S6a); `hatch_small_turtle` severity tier now compares **per-proc** DOFs
  against Ncrit, call-site comment aligned with the 80k gate, "give" → "given"
  (S6b); serial solver spelled `<real,int>`.
- `cl_SolverDistMatrix.cpp` — `set_permutation_switch` returns false for
  `SolverType::STRUMPACK` (S5 follow-up): STRUMPACK reorders internally, the
  pre-permutation would only reorder twice. PETSc behavior unchanged.
- `tests/sparse/test_Solver.cpp` — `DefaultParametersValid` expects
  `CompressionMethod::AUTOMATIC`.

## Behavior changes to note

- **MC64 matching stays ON by default** (STRUMPACK behavior preserved;
  the first landing defaulted it off and was reverted — see header note).
  `matching : off` in the solver section is an opt-out for verified cases.
- Explicit `compression scheme : off/blr` is now honored verbatim; the size
  heuristic applies only for the (new) default `automatic`.
- An input file with `reordering scheme : metis` + STRUMPACK no longer
  triggers the BELFEM-side pre-permutation (STRUMPACK still gets
  PARMETIS/METIS as before).

## Verification

- All edited sources + test pass `g++ -fsyntax-only` with the real
  `flags.make` flags (BELFEM_STRUMPACK defined).
- **Confirmed working by Christian 2026-07-06:** the re-landed batch
  (matching ON default) builds and runs well; the ZERO_PIVOT abort did
  not recur.

## Open

- **S4 zero-copy value scatter** — ~~implement~~ **implemented 2026-07-06**
  (`DistMatrixCSR`: `mAllValues` staging deleted; `distribute_values`
  scatters straight from the matrix data via the previously-unused
  `mOffsets`; rank-0 `values()` returns a pointer into the live data —
  safe because STRUMPACK's `CSRMatrixMPI` ctor copies on set/update).
  Saves the per-solve O(nnz) pack loop and nnz×8 bytes on rank 0.
  Pending: serial + 2/4/8-rank smoke test by Christian (git backup in
  place).
- S6e — confirm which STRUMPACK version SCLS actually links (≥7.1.2 has the
  small-problem hang fix).
- Optional: revisit PARMETIS vs sequential METIS mapping for the MPI solver
  (STRUMPACK's own default keeps sequential METIS for ordering quality).
