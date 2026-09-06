# Devlog 2026-08-25 — Gas Models: const Accessors, mutable Scratch

**Date:** 2026-08-25
**Topic:** Const correctness across the gasmodels hierarchy — property accessors made
`const`, memoization caches and preallocated scratch made `mutable`, with the compiler
used as the oracle for both directions
**Module:** `src/physics/gasmodels`
**AIs involved:** Claude
**Claude Confidence:** high on the mechanical result (compiler-proven), high on the
`mutable` inventory (each entry proven necessary by removal), medium on the claim that no
*further* method could be const (only the ones reached were tested)

## The question

Christian asked whether `const` + `mutable` is cleaner than the previous all-non-const
signatures, or whether the old version was "more honest", and which is better practice in
scientific code.

The answer taken: `const` on the evaluators is both cleaner **and** more honest. The C++
meaning of `const` on a member function is *logical* constness — "this call does not change
the observable state" — which is exactly true of `cp( T, p )`. A memoization cache is the
textbook justification for `mutable` (Meyers, *Effective C++*, Item 3; Stroustrup §16.2.9.3).
The all-non-const version was not more honest, it was less informative: it said "this may
change the gas" of `cp` and of `remix` alike, and the reader could not tell them apart.

One caveat had to be written down rather than assumed away: `const` normally implies
"safe to read concurrently" to anyone coming from modern C++, and here it does not, because
the const evaluators write the shared cache. That is consistent with the framework policy
(MPI, not threads) but it needed to be stated at the class.

## What was done

The split now runs through the whole hierarchy — `Gas`, `EoS`, `EoS_Idgas`, `EoS_Cubic`,
`AlphaFunction` and its four subclasses, `Helmholtz`, the H2/O2/CH4/N2 species equations,
`HelmholtzTransport` and the methane transport model:

| Group | `const` |
|---|---|
| evaluates a property of a fixed mixture | yes |
| changes what the gas is (`remix`, `set_mass_fractions`, equilibrium) | no |
| builds or rebuilds internal tables (splines, reference point) | no |

63 members carry `mutable`, in two classes: memoization caches (`Statevals`,
`mHelmholtzVals`/`mHelmholtzBits`, `mCubicStatevals`, the methane `mVals`/`mBits`) and
preallocated scratch (mixture work matrices, Cardano work vectors, the `delta^d` and
`tau^t` tables of the residual terms, the alpha-function help values).

Three changes were made instead of reaching for `mutable`, because a const overload or a
const reference was the honest fix:

- `Gas::heat_spline()`, `viscosity_spline()`, `conductivity_spline()` gained
  `const`-returning overloads rather than making the splines mutable.
- `prandtlmeyer::Wave` holds `const Gas &` instead of `Gas &`; it only evaluates.
- `EoS_Cubic::hdep`/`cpdep` bind `const Matrix< real > &` to the departure coefficients.

## Method: the compiler as the oracle

The propagation was not done by reading the call graph by hand. `const` was added to the
declarations and definitions of the evaluator set, then each translation unit compiled with
`-fsyntax-only`; every "assignment of read-only location" or "discards qualifiers" named a
member or a callee that had to be decided. That catches both directions: an over-applied
`const` fails on the mutation, an under-applied one fails as "no declaration matches".

**`mutable` cannot be checked that way** — an unnecessary `mutable` is silent. So each of
the 76 candidates was removed one at a time and all 14 translation units recompiled. Ten
turned out to be unnecessary and were dropped:

- `mComponents`, `mElements` — a `const Cell< RefGas * >` still yields a non-const pointee,
  so a const method can call through them without help. Marking the *container* mutable
  would have claimed that const methods may swap the components out.
- `mWorkMu`, `mWorkLambda` — written only by the viscosity/conductivity spline build.
- `mWorkVectorRAND0..2`, `mWorkMatrixRAND`, `mPivotRAND` — the RAND equilibrium solver
  changes the composition, so it was never a const path.
- `mFormationTable` — written once at setup, only *read* by the const `Gibbs`/`Hf`.

That last pair is the sharpest result of the sweep: `mFormationTable` is not mutable while
`mFormationWork` next to it is, and the difference is real — the table is built once, the
work vector is per-call scratch.

## Defects found on the way

- **`cl_GM_EoS_Nitrogen.hpp` constructor initializer list was malformed.** The paren closing
  `mPhi( ... )` was misplaced, so the new `mVapN` coefficient vector was being passed as a
  second argument to `mPhi`'s constructor. Fixed.
- **`cl_GM_EoS_Nitrogen.cpp` does not compile, and did not at `HEAD` either.** It is not in
  `CMakeLists.txt`, which is why this has gone unnoticed. Three pre-existing defects:
  the constructor never calls the `Helmholtz` base constructor; `compute_phir_t`/`_tt`
  contain `mJ( k ) *  - mBeta( kk )` and `mTauPowJ( k ) *  * mF( k )` with a missing operand;
  and `kk` is used without declaration in two scopes. Left alone — unfinished work, not a
  regression. `mVapN` is likewise declared but not yet read.
- The declarations of `evaluate_viscosity_interaction`, `evaluate_conductivity_interaction`
  and `dhdp_differential_quotient` had drifted out of sync with their definitions in the
  working tree; reconciled.

## Verification

Not "verified" in the executable-gate sense — no test was run.

- **Reviewed, compiler-proven:** all 14 gasmodels translation units and all 19 test
  translation units in `tests/physics/gasmodels` compile clean with `-fsyntax-only` under
  the flags taken from `cmake-build-debug/compile_commands.json`.
- **Not run:** `make check`. Neither build tree has `USE_GASMODELS=ON` (both caches say
  `OFF`), so the module and its tests are not currently built. Running them needs a
  reconfigure with `-DUSE_GASMODELS=ON`.
- Adding `const` cannot break a caller — a non-const object still calls a const method — and
  nothing outside `src/physics/gasmodels/` derives from these classes or takes a
  `Gas::*` member pointer. The only hits for those patterns are in
  `nonfree/physics/gasmodels/`, an untracked stale copy from 2026-07-21 that the live
  build does not use.

## Open

- `gastables::RefGas` was done later the same day — see
  `dl20260825_gastables_const_correctness.md`. The guess recorded here, that it "caches the
  same way", turned out to be **wrong**: `RefGas` does not memoize at all and needed no
  `mutable`. `mComponents` was nevertheless kept as `Cell< RefGas * >` by choice.
- `EoS::parent()` returns a non-const `Gas *`; left non-const since nothing in the module
  calls it.
- `nonfree/physics/gasmodels/` and `src/physics/gasmodels/` both define a CMake target named
  `gasmodels`, and `nonfree/physics/CMakeLists.txt` still adds its copy. Only reachable with
  `USE_NONFREE` and `USE_GASMODELS` both on; not investigated.
