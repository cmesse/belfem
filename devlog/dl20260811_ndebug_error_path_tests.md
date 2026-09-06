# Error-Path Tests in a Production Build: `BELFEM_ERROR` Now Has a Runtime Reaction

**Date:** 2026-08-11
**Purpose:** Record why 5 of 13 test binaries aborted under `USE_DEBUG=OFF`, the three-AI
round that chose the fix, and what landed
**Module:** `src/core` (assert), `tests/` (all 13 mains), `tests/doc`

## Symptom

`./test/test_ode` announced 43 tests and died on test 15; `./test/test_mesh` announced 46 and
died on test 2. Both printed the error dragon and then `MPI_ABORT`. The natural reading — a
flaky or false-positive test — was wrong in a way that mattered: `MPI_Abort` kills the process,
so gtest never reached the remaining tests. 28 and 44 tests respectively did not run. A
truncated suite reports as one failure, which is why this had gone unnoticed.

## Root cause

`BELFEM_ERROR`'s *check* is compiled into every build; its *reaction* was not. Both macros
funnel through the `assert::error` function template, which chose at compile time:

```cpp
print_errorbox( aLocation, aTask, aCheck, tMessage );
#if !defined( NDEBUG ) || defined( DEBUG )
    throw aException;          // debug
#else
    error_abort();             // release -> MPI_Abort( gComm.world(), 1 )
#endif
```

So `EXPECT_THROW` on a `BELFEM_ERROR` was a debug-only construct. Seventeen tests across five
binaries did exactly that. The tree in use compounded it: `cmake-build-debug` carried
`USE_DEBUG=OFF` as a leftover cache override (`USE_DEBUG-MODIFIED:INTERNAL=ON`), so a directory
named "debug" was producing `-O2 -DNDEBUG` for libraries and tests alike.

Two test-plan documents had written the wrong inference down as a rule —
`tests/doc/tests_12_ode.md` §12 and `tests/doc/tests_11_io.md` §4 both said `BELFEM_ERROR` is
"always active" and therefore must **not** be guarded. True of the check, false of the throw.
The tests were correct by their documentation; the documentation was wrong. `doc/coding_philosophy.md`
had the mechanism right all along, which is why this survived: the authoritative document
agreed with reality and the derived one did not.

A second, quieter defect: `assert::error` is a function template instantiated in the *calling*
translation unit. Library errors instantiate in library `.cpp` files, so the throw/abort choice
is fixed with the library's flags. Adding `-DDEBUG` to a test target would not have restored
throwing for `cl_BDF.cpp:27`, and a `#ifndef NDEBUG` guard in a test file reads the wrong TU's
flags. It only worked because the whole tree is built with one flag set.

## Decision

Five options were put to a blind three-AI jury round (Claude pre-registered, Codex and Grok
audited independently). All three converged on a runtime switch, and both auditors rejected the
alternatives on the same grounds: an always-throw design moves the fatal path onto exception
machinery that `doc/coding_philosophy.md` explicitly tells callers not to rely on, and a throw
escaping `main` on one rank terminates that rank while its peers block in a collective. gtest
death tests were rejected because they `fork()` after `MPI_Init`. Guarding the tests out of
release was rejected by the requirement — release is the configuration that ships, and
`BELFEM_ERROR` covers exactly the failures that survive into it.

**Unanimity was not the evidence.** Agreement is the weakest rung of the ladder, and
process-death semantics under MPI is a safety boundary; the choice was ratified by Christian.

The one place the auditors disagreed was load-bearing and changed the design. The proposal as
pre-registered defaulted the flag to `false`. Codex refuted it: `USE_DEBUG` defaults `ON`, so
debug builds throw today, and 52 `#ifndef NDEBUG` blocks across `tests/` depend on that. A
blanket `false` would have broken all 52. The flag is therefore initialized to the build's
compile-time behaviour, and only test binaries override it. Confidence: high — the 52-site
count was measured, not estimated.

## What landed

`src/core/assert.hpp`
- `BELFEM_ASSERTIONS_ACTIVE` exported next to the macro definitions, so tests guard on one
  owned symbol instead of restating `#ifndef NDEBUG`. The two are not equivalent when both
  `NDEBUG` and `DEBUG` are defined.
- `throw_on_error()` / `set_throw_on_error()` declared; the `#if` inside `error()` replaced by
  a runtime branch on a path that is already fatal.
- Kept as two separate predicates deliberately. After this change a release binary in throw
  mode can catch a `BELFEM_ERROR` while `BELFEM_ASSERT` still expands to nothing — "assertions
  exist" and "errors are catchable" are different questions and must not share a macro.

`src/core/assert.cpp`
- The flag is a single definition here, initialized to `BELFEM_ASSERTIONS_ACTIVE`. It must not
  be a `static` in the header: `error()` is a template, so a header static would give every
  translation unit its own copy and a test binary could not change the reaction of an
  already-compiled library.

`tests/*/test_*_main.cpp` (13 files)
- `set_throw_on_error( true )` after `gComm.init`. Every rank sets it — MPI is multi-process
  and each process has its own copy.

`tests/doc/tests_11_io.md`, `tests/doc/tests_12_ode.md`
- The two rules corrected, and the io audit checklist item with them.

## Verification

| Configuration | Result |
|---|---|
| `USE_DEBUG=ON`, before the change | 13/13 pass — first run in this tree where the ASSERT tier compiled in at all |
| `USE_DEBUG=OFF`, before the change | 5/13 fail: `math`, `sparse`, `mesh`, `io`, `ode` |
| `USE_DEBUG=OFF`, after the change | see Status below |

The five release failures were predicted by name before the run and matched exactly, including
the triggering test and abort site in each: `cl_Graph_Vertex.cpp:84`, `cl_Solver.cpp:56`,
`cl_Element_Factory.cpp:197`, `hdf5_tools.hpp:168`, `cl_BDF.cpp:27`. No unexpected failures —
the `-O2` versus `-Og` tolerance risk flagged as a possible confounder did not materialise.

That prediction also settled the round's one source-trace-only finding. Codex had caught an
omission in the original sweep: `tests/io/test_HDF5.cpp:195` is a real `BELFEM_ERROR` path
(`FileMode::NEW` truncates via `H5F_ACC_TRUNC`, then `load_data` hits
`BELFEM_ERROR( hdf5::dataset_exists(...) )`). The `io` binary aborted at `hdf5_tools.hpp:168`,
inside the cited range and in exactly the named test.

## Known gap, stated rather than papered over

Under this design a test verifies that the check fires with the right message, **not** that the
process aborts. Release and test runs take different reactions by construction. That is the
price of not forking after `MPI_Init`, and it is the right trade — but no one should read a
green release suite as "the shipped abort behaviour is tested". `error_abort()` is three lines
and separately testable without MPI if that ever matters.

Related: release throw-mode test runs should stay single-rank unless a test is explicitly
rank-synchronised. Rank-asymmetric error paths exist (`cl_HDF5.cpp` opens files on rank 0 only),
and in throw mode those no longer take the whole job down.

## Not done here

- The 52 bare `#ifndef NDEBUG` guards in `tests/` still work correctly, since CMake never sets
  both macros. Migrating them to `#if BELFEM_ASSERTIONS_ACTIVE` is hygiene, not a fix, and a
  52-site mechanical edit was out of scope for an unsupervised session.
- `cmake-build-debug` still carries the `USE_DEBUG=OFF` override. Confirmed a leftover; both a
  debug and a release test tree are wanted, since the point of this change is that the suite
  must pass in both.
- `src/sparse/pardisotools.f90` emits two rank-mismatch warnings. gfortran has no
  `diagnostic push/pop`, so the fix is to declare the placeholder dummies as size-1 arrays
  rather than suppress. Left pending a decision. Compiling a patched copy showed it also
  surfaces a masked warning: the solve call passes 17 arguments ending in `gDPARM` while the
  other three `pardiso` calls pass 16, which may be a real ABI mismatch and deserves its own look.
