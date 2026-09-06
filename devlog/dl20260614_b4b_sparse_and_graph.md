# Devlog: B4b Sparse Linear Algebra and Graph Theory Whitepaper

**Date:** 2026-06-14
**Purpose:** Read-only documentation pass for Tier B4b (sparse matrix class + solver/partitioner
backend wrappers + graph theory module).
**Module:** src/sparse, src/math/graph

## Context

Source-grounded whitepaper for the linear-algebra-and-graph substrate beneath FEM assembly/solving.
Scope: the matrix CLASS and the backend WRAPPERS, not their FEM use (assembly B3a, solver-driving B3c,
mesh partitioning B3.0 are call-site seams only). Read-only; `./archive` and `./nonfree` positively not
accessed.

**Path correction:** the sparse-matrix code is under `src/sparse/`, not `src/linalg/` as the brief
assumed (`src/linalg/` is dense `Vector`/`Matrix`). Graph is `src/math/graph/`.

## Method

Independent multi-agent cross-check (read-only subagents over the same files): one auditing solver
wrappers, one the graph module, one test coverage, one a from-scratch verification of the two
highest-risk facts (wired-vs-provisioned backends; UMFPACK->SuperLU + PARDISO-MKL), one gathering
in-repo license evidence. The independent verifier received no prior findings; it agreed on both
high-risk facts and surfaced the PARDISO contradiction independently.

## Key findings

- **SpMatrix:** CSC/CSR (`SpMatrixType`), switchable C++/Fortran indexing (`SpMatrixIndexingBase`),
  graph/dense/raw/HDF5 construction, runtime index-function-pointer dispatch (4 variants), COO indices
  for MUMPS. Tested: CSC/CSR, 0/1-based, dense-ctor, multiply, transpose, COO. Not tested:
  construction-from-graph (the production path) and save/load.
- **Solver wrappers:** abstract `solver::Wrapper` base + per-backend subclass + factory-by-switch
  (not pimpl); `BELFEM_<X>` guards with two-layer stubbing (factory throws when OFF). Six backends all
  have real implementations. Default build enables **only SuperLU** (USE_SUPERLU ON; UMFPACK/MUMPS/
  PARDISO/PETSc/STRUMPACK OFF). Every backend solve test is `#ifdef`-gated on a default-OFF macro, and
  **SuperLU - the only default solver - has no test**. `tests/old/` is dead; `USE_TEST` default OFF.
- **UMFPACK->SuperLU migration:** documented and license-driven (BSD SuperLU replacing GPL UMFPACK,
  `config_suitesparse.cmake:3-9`, `dl20260608_...:9`). Soft migration: UMFPACK still wins the default
  `#elif` chain when both are linked; SuperLU is default only because USE_SUITESPARSE defaults OFF.
- **PARDISO:** build-gated to Intel MKL (`FATAL_ERROR` otherwise) but the wrapper carries
  standalone-PARDISO license-error strings (-10/-11/-12) - contradiction flagged, author to confirm.
- **Eigensolver:** ARPACK-NG confirmed (`arpacktools.cpp:79` "ARPACK-ng" + `dnaupd`/`dneupd`),
  default-ON, untested. **SLEPc absent** from src/ and config/ -> roadmap only.
- **Graph:** BFS, DFS, RCM (symrcm), pseudo-peripheral x2, connected partitions - unconditional,
  tested. Partitioners METIS/ParMETIS/SCOTCH/PT-SCOTCH wrapped + `BELFEM_<X>`-guarded; METIS+SCOTCH
  default ON; parallel variants toggled in lockstep with serial via USE_METIS/USE_SCOTCH (no separate
  options). Only METIS is tested (smoke-level).
- **Licensing:** BELFEM = BSD-3-Clause by text. In-repo license statements exist only for SuperLU
  (BSD-3) and UMFPACK/SuiteSparse (GPL); the explicit GPL-avoidance narrative is in
  `config_suitesparse.cmake`. No in-repo license text for the other backends.
- **Seam:** STRUMPACK -> B1/comm: `BELFEM_STRUMPACK` makes comm request `MPI_THREAD_MULTIPLE`.

## Output

- Added `tmp/whitepaper/B4b_sparse_and_graph.md` (in-scope files first; confidence tags; solver +
  partitioner STATUS tables with license column; per-clause claim verdicts; independent cross-check;
  open questions; assumptions). Pure ASCII.

## Verification

- Confirmed `./archive` and `./nonfree` not accessed.
- Confirmed whitepaper is ASCII-only.
- No source edited, nothing compiled, no tests run (read-only documentation task).
