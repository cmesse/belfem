# YBCO golden test fails on Ubuntu: unqualified `abs()` on doubles is the integer `abs`

**Date:** 2026-09-13
**Purpose:** Record why `YBCOThermalConductivity.RepresentativeTemperaturesAreFinitePositive`
fails on the `ubuntu_patch` branch (release build, GCC 13.3, libstdc++), what the executable gate
showed, and the fix applied after Christian's approval.
**Module:** src/physics/materials (`cl_Material_Metal.cpp`, `cl_Material_YBCO.cpp`,
`cl_Material_Alloy.cpp`, `cl_Material_HastelloyC276.cpp`); the same pattern in
`src/fem/kernel/cl_FEM_Tmatrix.cpp`, `src/mesh/fn_Mesh_ratio.cpp`, `src/executables/belfem.cpp`,
`src/executables/hphiTrun.cpp` (retired, not built)

## Symptom

`build/test/test_physics` (Release, Armadillo, gfortran/g++ 13.3.0-6ubuntu2~24.04.1, x86_64):
λ(20 K) = 37.7252 against golden 38.0820 (−0.94 %), λ(77 K) −0.11 %, λ(100 K) −0.07 %,
λ(300 K) −0.018 %. Tolerance is 1e-9 relative. Reproduces with the test run alone
(`--gtest_filter`), so it is not a test-order effect.

## Root cause

`Metal::set_RRR()` (`cl_Material_Metal.cpp:247`) solves `( rho_i( T_ref ) + rho_0 ) / rho_0 = RRR`
for `rho_0` with a secant loop whose exit test is `while ( abs( f ) > 1e-12 )`. The call is
unqualified and `f` is a `real`. In every affected translation unit on this toolchain the name
resolves to the C `int abs(int)` from glibc's `stdlib.h`, because libstdc++ puts the
floating-point overloads only into `namespace std` (`<bits/std_abs.h>`); the global
`::abs(double)` appears only when libstdc++'s own `<math.h>`/`<stdlib.h>` wrapper is included
somewhere in the chain, and here it is not. `abs( -0.37 )` is therefore `0`, the residual
(−0.485 at entry) is truncated to zero, the loop body never runs, and `rho_0` keeps the bisection
midpoint `1.0 · rho_i_ref / ( RRR − 1 )` = 1.2367e-8 Ω·m instead of the converged 1.2245e-8.
The electronic term `L0·T / ( rho_i + rho_0 )` dominates λ, and `rho_0` dominates the
denominator at low T, which is exactly the observed decay of the error with temperature.

The goldens were baked 2026-08-25 on a toolchain where the same expression converged (libc++ on
macOS exposes `::abs(double)`; on libstdc++ it depends on the include chain — which of the two the
baking machine was is not recorded, and the nightly CI has not been checked for this test).

Verified, not just reviewed: `mpicxx -fsyntax-only` of each listed `.cpp` with its own
`flags.make` plus an appended `static_assert( is_same< decltype( abs( -0.5 ) ), double > )`
fails in all of `cl_Material_Metal.cpp`, `cl_Material_YBCO.cpp`, `cl_Material_Alloy.cpp`,
`cl_Material_HastelloyC276.cpp`, `cl_FEM_Tmatrix.cpp`, `fn_Mesh_ratio.cpp`, `belfem.cpp`,
`cl_SimplicialComplex.cpp`.

## Executable gate

A scratch program linked with `mattest`'s link line (`CMakeFiles/mattest.dir/link.txt`)
constructs `YBCO`, prints the residual of the RRR equation, then reruns the `set_RRR` secant
iteration with `std::abs` from a subclass and rebuilds the λ spline with the converged `rho_0`:

| | rho_0 [Ω·m] | RRR residual | λ(20) rel. to gold | λ(77) | λ(100) | λ(300) |
|---|---|---|---|---|---|---|
| as built | 1.236735e-08 | −4.85e-01 | −9.37e-03 | −1.13e-03 | −7.16e-04 | −1.78e-04 |
| loop with `std::abs` (7 iterations) | 1.224490e-08 | 0 | −1.9e-12 | 4.8e-12 | −1.2e-12 | 2.7e-10 |

The as-built row reproduces the failing test's printed λ values to every printed digit; the
refitted row is inside the 1e-9 tolerance at all four temperatures. The golden values are
therefore correct and the Ubuntu build is wrong.

## Other sites with the same defect (all reviewed, only the YBCO path verified)

| Site | What the loop does | Effect under integer `abs` |
|---|---|---|
| `cl_Material_Metal.cpp:247` | RRR → `rho_0` secant | skipped; every `Metal` on libstdc++ carries the initial guess (YBCO −0.49 residual; other metals unmeasured) |
| `cl_Material_Metal.cpp:499` | `invert_debye` secant | skipped unless residual ≥ 1; only metals that derive θ_D from c_p |
| `cl_Material_YBCO.cpp:264` | Newton for the α plateau at ~350 K | skipped; plateau stays at 350 K exactly |
| `cl_Material_YBCO.cpp:369` | bisection for the c_p switch `T2` | skipped; `T2 = exp( 5.8 )` |
| `cl_Material_Alloy.cpp:670`, `cl_Material_HastelloyC276.cpp:136` | same kind of root finds | not read in detail |
| `cl_FEM_Tmatrix.cpp:50,88` | sparsity test `abs( a_ij ) > BELFEM_EPSILON` on `Matrix< real >` | entries with |a_ij| < 1 are dropped from the T-matrix — potential correctness defect, not verified |
| `fn_Mesh_ratio.cpp:36` | secant, `tF` starts at `BELFEM_REAL_MAX` | double→int conversion out of range is UB; loop likely skipped |
| `belfem.cpp:279` (`hphiTrun.cpp:167`, retired) | thermal/magnetic time sync `abs( t − t_thermal ) > 1e-12` | lags below 1 s read as zero; the thermal sub-stepping loop can exit early |
| `fn_Smith.hpp`, `cl_SimplicialComplex.cpp` | integer chain coefficients | integer `abs` is the intended overload; cohomology core, closed to AI |

## Fix (applied, approved by Christian: "replace all remaining abs() calls that refer to pure C")

Every floating-point call at the sites in the table is now `std::abs(`: eleven sites in eight
files (`cl_Material_Metal.cpp` ×2, `cl_Material_YBCO.cpp` ×2, `cl_Material_Alloy.cpp`,
`cl_Material_HastelloyC276.cpp`, `cl_FEM_Tmatrix.cpp` ×2, `fn_Mesh_ratio.cpp`, `belfem.cpp`,
`hphiTrun.cpp`). The homology sites stay: integer coefficients, integer `abs` intended, closed
core. No behaviour change on a platform where the double overload was already found. Each
edited unit that has an object in `build/` was re-checked with `mpicxx -fsyntax-only` under its
own `flags.make` (with `-Werror`) plus a `static_assert` that `std::abs( -0.5 )` is `double`:
all seven pass, so `<cmath>` reaches every site transitively. Not rebuilt or re-run here — the
build is Christian's; `test_physics` after the rebuild is the gate. Whether the Rocky and nightly-CI builds were on the integer path too is
worth one `static_assert` probe there, because if they were, every cryogenic `rho_0` in
production runs has been the initial guess.

Candidate tripwire for `doc/lessons_learned.md`: an unqualified `abs`, `sqrt`, `pow`, `exp` on a
floating-point argument inside `namespace belfem` is platform-dependent; `std::` or `<cmath>`
`using` is required, and `-Wall -Werror` does not catch it (GCC's `-Wabsolute-value` is C only).

## Status

Diagnosis verified (executable gate above, this host only). Fix applied at eleven sites,
syntax-checked per unit, not yet built or tested in the tree: reviewed, not verified, until
`make` and `test_physics` run. Scratch probes live in the session scratchpad, not the tree.
