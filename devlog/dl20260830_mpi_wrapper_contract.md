# The Wrapper Contract: allreduce, comm_abort, and L-21

**Date:** 2026-08-30
**Purpose:** Record the hard policy ruling on third-party access, the jury round on the two new MPI
wrappers, and the four-site replacement that closed the violations
**Module:** `comm`, `sparse`, `core`

## The ruling

Christian, 2026-08-30, hard policy, now L-21 in `doc/lessons_learned.md`: **third-party libraries
are NEVER accessed directly — only through the dedicated wrapper layer**, in the fashion the
codebase applies elsewhere. A missing operation extends the wrapper first.

The trigger was a sweep finding `MPI_Allreduce` open-coded in `cl_SolverPETSC.cpp` and
`cl_SolverSTRUMPACK.cpp`, and `MPI_Abort`/`MPI_Initialized`/`MPI_Finalized` in `assert.cpp` — plus
a fourth site the sweep's own author supplied: `test_commmpi_main.cpp`, written earlier the same
day. Three of the four violations were fresh code by sessions that knew the tree well. That is what
makes this a tripwire, not a style note: the violations were nobody's ignorance, they were the
default behaviour of writing "the normal way".

## The jury round on the draft wrappers

Christian drafted `allreduce` and `abort` in `commtools.hpp` and asked specifically whether the
integer datatypes were right, given `int_t` splits 32/64-bit by suite. Claude pre-registered seven
findings sealed from the auditors; Codex (terra/high) and Grok (4.6/high) audited blind. Both:
**revise**. Full three-way convergence on the core defects, and the union was strictly better than
any one voice:

- the reduce call had the wrong arity for the MPI function it named — latent, because a template
  body is only type-checked at instantiation and nothing instantiated it yet
- `aRecv` was `const T*`; MPI's recvbuf is writable
- the serial branch returned without copying send→recv
- the count was `int_t`, which narrows silently above 2^31-1 on the 64-bit suite — the answer to
  Christian's question: `comm_type<T>` makes the *element* type safe on both widths; the *count*
  was the hazard, and it is now MPI's own `int`
- **the abort wrapper used `gComm.world()`** — re-introducing, hours later, the exact
  `error_abort` contract fixed that same day (`mComms` empty between `MPI_Init` and the push;
  `MPI_Abort` erroneous outside the init/finalize window)

Unique catches: **Codex** — the serial abort branch returning 0 lets `assert::error()` *resume
after a failed check* in a no-MPI build; complex `comm_type<>` types make `MPI_MAX` meaningless,
now a `static_assert`. **Grok** — `comm_check` on the abort path is recursion
(`comm_check` → `BELFEM_ERROR` → `error_abort` → abort); the name `comm_abort`, matching the
`comm_barrier`/`comm_check` siblings; keep `MPI_Op` OUT of the signature, since a parameter
re-leaks the vendor token the wrapper exists to hide.

## What landed

- `commtools.hpp`: `allreduce` — `void`, `comm_check`-wrapped `MPI_Allreduce`, `T*` recv, `int`
  count, arithmetic-only by `static_assert`, serial identity copy. MAX is deliberate and
  documented, not a parameter.
- `cl_Communicator.hpp` + `commtools.cpp`: `comm_abort` — plain function, `[[noreturn]]`,
  `MPI_COMM_WORLD` with lifecycle guards, unconditional `std::abort()` fall-through, never routed
  through `comm_check`. Declared in `cl_Communicator.hpp` rather than `commtools.hpp` for a
  layering reason found during implementation: core's include path is core/comm/containers only,
  and `commtools.hpp` drags `cl_Vector.hpp` from linalg — `assert.cpp` could not have included it.
- Four sites replaced: PETSC verdict fold, STRUMPACK failure fold, `error_abort` (now a one-line
  delegate; its contract comment moved to the wrapper so every future caller inherits it), and the
  Tier 2 test main (whose `#ifdef BELFEM_MPI` and `<mpi.h>` include disappeared entirely — the
  wrapper's serial branch carries the no-MPI case).
- Residual `MPI_` hits outside `src/comm`: comments, and `MPI_Comm` *type* plumbing in the
  petsctools glue, which is itself wrapper layer.

## Evidence

`mpicxx -fsyntax-only` with each TU's own `flags.make` flags (CXX_FLAGS + DEFINES + INCLUDES):
`cl_SolverPETSC.cpp`, `cl_SolverSTRUMPACK.cpp`, `assert.cpp`, `commtools.cpp`,
`test_commmpi_main.cpp` all clean — and the solver TUs instantiate `allreduce<int>`
unconditionally, so the template body is compiled, not merely parsed, closing the exact latency
that hid the draft's arity error. The no-MPI shape compiles too, instantiating the serial branch.
**Gate RAN 2026-08-30 (Christian): rebuild + `make check` green.** That verifies: the wrappers
compile and link in every registered TU; `allreduce<int>` executes for real at np=2 and np=4
through the Tier 2 exit fold (the all-pass case — MAX over zeros — on live communicators); the
solver TUs carrying the replaced sites build. Deliberately NOT claimed: `comm_abort` executing
(a green suite cannot reach the abort path, by definition — it is compiled and its contract is the
verbatim move of code that ran for months), and the failure-fold through the new wrapper (the
rank-3 probe measured that mechanism this morning against the direct call the wrapper replaced;
re-running the probe against the wrapper would lift that last inference to a measurement, at the
cost of one build cycle — noted, not owed).

Exchange: `tmp/ai_exchange/commtools_allreduce_abort.md` (pre-registration sealed before dispatch,
reconciliation appended).
