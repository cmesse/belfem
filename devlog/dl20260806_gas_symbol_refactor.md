# Gas Modules: Physical-Symbol Naming Refactor

**Date:** 2026-08-06
**Purpose:** Drop the `a`-prefix from physical-quantity arguments across
`src/physics/gastables` and `src/physics/gasmodels`, so the API and kernels read
like the thermodynamics they implement (`shock( T1, p1, u1, alpha, … )`).
**Module:** physics/gastables, physics/gasmodels
**Decision:** Christian, 2026-08-06. Convention recorded in `CLAUDE.md`
(Canonical thermophysical symbols).

## What changed

**2590 identifiers across 39 of 81 files**, applied mechanically with
word-boundary matching (script kept in the session scratchpad):

| | | | | | |
|---|---|---|---|---|---|
| `aT`→`T` (1273) | `aP`→`p` (581) | `aV`→`v` (162) | `aH`→`h` (14) | `aS`→`s` (4) | `aU`→`u` (8) |
| `aT1`→`T1` (75) | `aP1`→`p1` (70) | `aU1`→`u1` (52) | `aA1`→`A1` (13) | `aT0`→`T0` (14) | `aP0`→`p0` (14) |
| `aT2`→`T2` (56) | `aP2`→`p2` (63) | `aU2`→`u2` (54) | `aA2`→`A2` (29) | `aAlpha`→`alpha` (17) | `aBeta`→`beta` (20) |
| `aMu`→`mu` (17) | `aLambda`→`lambda` (7) | `aCp`→`cp` (12) | `aC1..3`→`c1..3` (35) | | |

Plus two readability follow-ups the rename exposed: `fn_GT_idgas_mu.cpp` had a
local `mu` holding the **dipole moment** (renamed `dipole` — in this scheme `mu`
is the viscosity), and a `@param T : T` doxygen stub in `cl_GM_Helmholtz.hpp`
gained real text.

## Why it is safe — the collision survey that preceded it

The hazard is that **every target name is already a method** on these classes
(`Gas::T`, `::p`, `::v`, `::h`, `::s`, `::u`, `::cp`, `::mu`, `::alpha`,
`::beta`, …), so a parameter of the same name shadows the accessor. Four checks
were run before any edit:

1. **Unqualified self-calls: zero.** All 82 self-calls in `cl_Gas.cpp` use
   `this->`; a module-wide scan found only two bare hits, both inside string
   literals. Shadowing therefore cannot turn a call into a compile error.
   *This is now a load-bearing invariant — see `CLAUDE.md`.*
2. **Duplicate declarations: none.** Only nine locals in the two modules carry a
   target name, and none sits in a function taking the matching argument. Note
   `cl_Gas.cpp:2758` already had `const real s = this->s( T1, p1 );` and
   `prandtl_meyer` already hand-rolled `real & T = aT2;` aliases — the codebase
   was doing this by hand, which is what prompted the refactor.
3. **Injective map.** No two source names share a target, so no signature can
   collapse to duplicate parameters (verified after the fact as well).
4. **String literals: zero hits.** No error message or format string changed.

`aA`, `aB` and `aC` were **excluded** — they are comparator pointees
(`cl_GT_ComparisonObjects.hpp`), species-label strings
(`cl_GT_RefGasFactory`), and a coefficient output vector
(`fn_GT_create_glue_poly.hpp`), none of them physical quantities. The area
arguments `aA1`/`aA2` in `expand()` were renamed; the bare `aA` was not.

## Open / carried

- **Not compiled.** Per the standing convention the build is Christian's. The
  edit is mechanical and the four checks above are static, but nothing here is
  compiler-verified. This lands on top of the same day's uncommitted
  correctness fixes, so a build failure could come from either — the symbol
  refactor is a separate commit for exactly that reason.
- **`t`-prefixed temporaries were left alone.** Dropping those too
  (`tT`→`T`, `tV`→`v`, …) is the natural second pass, but it is the one
  direction that *can* collide with the now-renamed arguments inside a single
  function body, so it wants a compiler in the loop rather than a static scan.
- `u` remains context-bound: velocity in `total`/`expand`/`shock`/
  `prandtl_meyer`, internal energy in `Gas::u( T, p )`. Documented, not
  "resolved" — both are standard notation.
