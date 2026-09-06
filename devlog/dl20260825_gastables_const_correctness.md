# Devlog 2026-08-25 — Gas Tables: const Evaluators, and the M()/R() Reference Audit

**Date:** 2026-08-25
**Topic:** The const-correctness exercise carried one level down from `gasmodels` into
`gastables`, plus a read-only audit of the archived channel/combustion/engine code for
callers that depend on `Gas::M()` / `Gas::R()` returning a live reference
**Module:** `src/physics/gastables`
**AIs involved:** Claude
**Claude Confidence:** high on the gastables change (compiler-proven, no `mutable` needed),
high on the archive audit being a true negative (search terms recorded below)
**Related:** `dl20260825_gasmodels_const_correctness.md`

## Summary

`RefGas` is now const-correct, 100 signatures across `cl_GT_RefGas.{hpp,cpp}`. The
module needs **no `mutable` whatsoever** — see below, this is the interesting result.
Separately, the archived combustion-era code was searched for the one usage pattern that a
future `M()`/`R()` change would silently break. None was found, and the search terms are
recorded here so the negative can be trusted later without redoing the work.

## gastables needed no `mutable`

The prediction going in was that `RefGas` would look like `Gas`: caching evaluators needing
`mutable` scratch. That was wrong, and worth writing down because it explains why this
module was so much cheaper to fix.

`HeatPoly`, `TransportPoly` and `GasData` were **already** fully const-correct — every
evaluator on them was const before this session. And `RefGas` does not memoize at all:

```cpp
RefGas::Cp( const real T ) const
{
    return this->find_heat_poly( T )->Cp( T );   // cl_GT_RefGas.cpp:1029
}
```

`find_heat_poly` (`cl_GT_RefGas.cpp:770`) is a plain interval lookup with no memo. The
class holds `mHeatPolys`, `mViscosityPolys`, `mConductivityPolys`, three `Spline`s, a
`GasData` and a set of flags — all of it built once by the factory and read-only thereafter.
There is nothing a property call would want to write. So `Cp`, `H`, `S`, `cp`, `h`, `s`,
`mu`, `lambda`, all their derivatives, the twelve `poly_*` and twelve `spline_*` dispatch
targets, `zero`, and the three `find_*_poly` helpers all became const with zero `mutable`
added. Confirmed by compiling all 21 module TUs and both test suites.

The contrast with `gasmodels` (65 `mutable` members) is the whole point: caching lives in
the mixture layer, not the species layer.

Two further changes:

- The twelve `RefGas::*` member function pointers (`mFunctionCp` … `mFunctiond2LambdadT2`,
  `cl_GT_RefGas.hpp:116-149`) took `const` in their types, since `set_mode`
  (`cl_GT_RefGas.cpp:834`) assigns the now-const `poly_*` / `spline_*` / `zero` members
  to them.
- `find_heat_poly` / `find_viscosity_poly` / `find_conductivity_poly` now return
  `const HeatPoly *` / `const TransportPoly *`. They were const methods handing out a
  writable pointer — a const hole. All 18 call sites only use const poly methods, so
  closing it cost nothing.

What stays non-const is exactly the build-time surface: `add_*`, `create_*`, `delete_*`,
`finalize*`, `fix_*`, `set_*`/`unset_*`, `set_mode`. Plus `data()` and the three
`*_spline()` getters, which deliberately hand out writable handles; `data()` already had a
`const` overload returning `const GasData *`.

## The M() / R() reference audit

`Gas::M( T, p )` and `Gas::R( T, p )` return `const real &` into the state cache. We are
considering returning `real` by value instead. Christian recalled deliberately wanting a
reference somewhere in the combustion-era channel model, so that M and R would track a
remix automatically — the binding

```cpp
const real & tR = mGas->R( T, p );   // bound once, follows later remixes
```

which a by-value return would silently freeze into a snapshot. That is the one pattern the
change could break, and it would break quietly.

**Result: the pattern does not occur.** `archive/channel`, `archive/combustion`,
`archive/boundarylayer` and `archive/engine` contain no reference or pointer binding to
`Gas::M()`/`Gas::R()` that outlives its statement.

| Directory | `Gas::R(T,p)` sites | `Gas::M(T,p)` sites | live-ref bindings |
|---|---|---|---|
| `archive/channel` | 9 | 4 (+1 commented out) | 0 |
| `archive/combustion` | 0 | 0 | 0 |
| `archive/boundarylayer` | 0 | 0 | 0 — uses `h`, `Pr`, `hd`, `shock`, `total`, `isen_T`, never `M`/`R` |
| `archive/engine` | 3 | 1 | 0 |

Every site consumes the value immediately in arithmetic, copies it into a `real`, streams
it, or passes it to a setter. The closest thing to caching is a **value** copy, at
`archive/engine/cl_EN_State.cpp:48-49`:

```cpp
this->value( BELFEM_ENGINE_STATE_R ) = mCombgas.R( aT, aP );
this->value( BELFEM_ENGINE_STATE_M ) = mCombgas.M( aT, aP );
```

which already freezes a snapshot and is therefore unaffected.

### What was searched, so the negative can be trusted

Over every `.cpp`/`.hpp` in the four directories:

- `(->|\.)(M|R)\s*\(` — full call-site inventory, each hit classified by hand as
  `Gas` / `RefGas` / geometry
- `const real\s*&\s*<ident>\s*=` and the non-const `real &` variant
- `const real\s*&\s*m[A-Z]` — reference **data members**: there are none in any header
  under these four directories
- member-initialiser lists `^\s*[,:]\s*m<ident>\(` filtered for `.M(` / `->M(` / `.R(` / `->R(`
- multi-line bindings (`real & x =` with the call on the next line)
- `=\s*&...(\.|->)(M|R)\s*\(` and `real\s*\*\s*<ident>\s*=\s*&` — address-of / pointer capture
- `\bauto\b` — deduced reference bindings; the only `auto` tokens in these directories are
  inside comments in `archive/boundarylayer/fn_BL_VanDirest.cpp:237,344`

**False positives that had to be excluded**, recorded so a future search does not trip on
them again:

- `belfem::channel::Geometry::R( const real & aX )` — a nozzle *radius*, unrelated to the
  gas constant. Hits in `cl_CH_Geometry.cpp`, `cl_CH_GeometryNozzle.hpp`,
  `cl_CH_GeometryCombustor.cpp`, `cl_CH_GeometryCylinderCombustor.cpp`, `itlr_combustor.cpp:79`.
- `gastables::RefGas::M()` — the zero-argument molar mass (`cl_GT_RefGas.hpp:189`). It does
  return `const real &`, but it is constant per species, so it carries no remix hazard.
  Used in `archive/channel/combustor.cpp`, `cl_CH_ChannelODE.cpp:555`,
  `archive/combustion/cl_CN_Scheme.cpp:89,100`, `cl_CN_Injector.cpp:133,135`,
  `archive/engine/cl_EN_Parameters.cpp:384,394` — all value uses.

### The one reference binding that does exist

`archive/channel/cl_CH_Boundarylayer.cpp:1645`:

```cpp
const real & pc = mGas.component( 0 )->data()->p_crit() ;
```

`GasData::p_crit()` returns `const real &` (`cl_GT_GasData.hpp:213`), but it is an immutable
per-species table constant, not a remix-dependent cache, and it is used inside the same
function body. Different, non-hazard class.

## Correction to the gasmodels session

The warning comment added to `cl_Gas.hpp:339` in the companion session claims that the
reference returned by `M()`/`R()` is overwritten by the next property call at a different
state. **That is wrong.** `Statevals::reset()` (`cl_GM_Statevals.hpp:224`) clears only the
bitset, and `update_Tp()` (`:168`) resets bits and writes only the T and p slots. The M and
R slots are written in exactly three places, all remix/init paths (`cl_Gas.cpp:678,684`,
`753,759`, `789,790`). The reference is stable for the lifetime of the `Gas` and changes
only on remix. Corrected in this session: the comment now states the actual behaviour and
notes that both methods currently ignore T and p, the parameters being there for a derived
model whose mixture dissociates.

## Verification

Not "verified" in the executable-gate sense — no test was run.

- **Reviewed, compiler-proven:** 21 gastables TUs, 14 gasmodels TUs, 7 gastables test TUs
  and 19 gasmodels test TUs all compile clean with `-fsyntax-only`.
- **Not run:** `make check`. `USE_GASMODELS=OFF` in both build trees, so neither module nor
  its tests are currently built.

## The API changes that followed (same session, decided by Christian)

With the audit in hand, the deferred questions were settled:

**`Gas::M( T, p )` and `Gas::R( T, p )` now return `real` by value.** The live-view
behaviour was real and worked, but it was undiscoverable: a caller who wanted it had to
know, and a caller who did not had no way to tell. A `real` is a register return, so the
copy costs nothing. All 17 in-tree call sites already copied immediately; the archive
audit above is the evidence that nothing outside relied on the reference. The header
comment written earlier in the day was rewritten again to describe the new contract.

**Paired const / writable accessors.** `Gas::data( index )` gained a const overload
returning `const gastables::GasData *`, matching the pair `RefGas::data()` already had.
`RefGas::heat_spline()`, `viscosity_spline()` and `conductivity_spline()` gained const
overloads returning `const Spline *`. Overload resolution picks by the constness of the
handle.

**Pointers stay stored writable.** `Gas::mComponents` remains `Cell< RefGas * >` rather
than becoming `Cell< const RefGas * >` — Christian's convention, and it keeps the build
paths simple. Where a const method only reads a component, it binds its own read-only
handle at the point of use:

```cpp
const gastables::RefGas * tComponent = mComponents( aIndex );
```

Applied in `Gas::h`, `cp`, `dcpdT` (the component-wise trio) and in both interaction
loops, `evaluate_viscosity_interaction` and `evaluate_conductivity_interaction`, where it
also hoists the `Cell` indexing out of the inner loop. Not applied in `cea_mu`,
`cea_lambda` and `print`, which make a single const call through the pointer inside a
condition or a printf argument — a named local there would be noise, not enforcement.

**Note the scope of that guarantee, honestly:** the local handle enforces read-only use
*within the function*. It does not stop the next person writing `mComponents( k )->finalize()`
in some other const method. Only `Cell< const RefGas * >` would enforce class-wide. That
trade-off was made deliberately in favour of simple build paths.

## Verification of the API changes

Beyond recompiling all 61 TUs, two compile-only probes were written:

- a **positive** probe driving the whole read-only surface through `const Gas &` and
  `const gastables::RefGas &` — caloric, transport, coefficients, `M`/`R` by value, the
  const `data()` and const spline overloads, the component-wise accessors and `total()`.
  Compiles clean.
- a **negative** probe asserting the guarantee bites. All four attempts are rejected:
  `remix()` on a const `Gas`, `finalize()` on a const `RefGas`, assigning `data( 0 )` to a
  writable `GasData *`, and assigning `heat_spline()` to a writable `Spline *`.

The second probe is the one that matters: it shows the const is load-bearing rather than
decoration.

## Open

- `Gas::component( index )` still returns `gastables::RefGas * &` and `components()` returns
  the writable `Cell`, with no const twin. Same family as `data()`; not touched because
  nothing needed it yet.
- `EoS::parent()` returns a non-const `Gas *`; left as is, nothing in the module calls it.
- `nonfree/physics/gasmodels/` and `src/physics/gasmodels/` both define a CMake target named
  `gasmodels`. Only reachable with `USE_NONFREE` and `USE_GASMODELS` both on; not investigated.
- `cl_GM_EoS_Nitrogen.cpp` is still unbuildable and still absent from `CMakeLists.txt`
  ( see the companion devlog ).
