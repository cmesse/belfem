# Builtin Alloy Materials Never Received Their RRR

**Date:** 2026-08-12
**Purpose:** Record the crash of `hphirun` on an input file declaring a builtin alloy with an `RRR` key, its four causes, and the fix
**Module:** physics/materials

## Symptom

Running the `tapestack3d` example aborted during `MaxwellFactory::create_materials()`:

```
terminate called after throwing an instance of 'std::runtime_error'
  what():  Property is not constant
```

with the backtrace ending in `Material::constant_property( RRR )` →
`Alloy::set_components()` → `Alloy::Alloy( "pb40sn60", …, aRRR = NaN )` →
`MaterialFactory::create_material()`.

The triggering section was well-formed and matches the documented contract
(`doc/input_schema.yaml`, `materials` section: `builtin` also accepts an alloy
formula, `RRR` applies to pure-metal-shaped materials, and `Alloy` reports
`MaterialType::PureMetal`):

```
solder
{
    builtin : pb40sn60 ;
    RRR : 5 ;
}
```

The input file was not at fault. Four defects in the material path were.

## Causes

1. **RRR arrived too late.** `MaterialFactory::MaterialFactory` constructed the
   material first and only then read the `RRR` key, applying it through
   `Material::set_RRR()`. For a pure metal that works — `Metal::set_RRR()`
   re-solves `rho_0` by secant and rebuilds the rho database. For an alloy it
   cannot: `Alloy::set_components()` derives `rho_0` from RRR *inside the
   constructor*, so the constructor is entered with `aRRR = NaN` and fails
   before the factory ever gets to `set_RRR()`.

2. **`Alloy` has no `set_RRR()`.** It derives from `SplineLookupTable` →
   `Material`, not from `Metal`, so even a successful construction would have
   hit the base `Material::set_RRR()`, which raises
   *"not implemented for this material"*. RRR has to go in at construction.

3. **A NaN comparison hid the intended error.** The constructor guarded the
   fallback with `aRRR != BELFEM_QUIET_NAN`. NaN compares unequal to
   everything, so the guard was always true and the pair `("RRR", NaN)` was
   pushed into the component list. `set_constant()` stored the NaN, and
   `is_constant()` — which is just `! std::isnan` — then reported the property
   as absent, so the debug assertion in `constant_property()` fired instead of
   the intended `BELFEM_ERROR( "requires a valid RRR > 1" )`. Under `NDEBUG`
   the assertion is compiled out and the clean message appears; the two build
   types therefore failed differently on the same input.

4. **Case mismatch on the in-formula RRR.** The constructor detects an RRR
   token in the formula case-insensitively (`string_to_lower(…) == "rrr"`),
   but `Alloy::set_components()` recognised it only as the exact string
   `"RRR"`. A lowercase formula such as `pb40sn60rrr5` therefore passed the
   first check and then reached `create_component( "rrr" )`, aborting with
   *"Unknown metal: rrr"*.

## Fix

`cl_MaterialFactory.cpp` — read the `RRR` key before construction and forward
it, then drop the `MaterialType::PureMetal` branch that applied it afterwards:

```cpp
const real tRRR = tMatSection->key_exists( "RRR" )
        ? tMatSection->get_real( "RRR" ) : BELFEM_QUIET_NAN ;

tMat = create_material( tMatType, tRRR ) ;
```

Every constructor reachable through that branch (`Copper`, `Silver`,
`Indium`, `Lead`, `WhiteTin`, `Iron`, `Alloy`) already applies a non-NaN RRR
itself, so nothing is lost — and the pure metals now build their rho database
once instead of twice, since `Metal::set_RRR()` populates it at the end.
`HastelloyC276` is `MaterialType::LookupAlloy` and never entered that branch,
so its behaviour is unchanged.

`cl_Material_Alloy.cpp`:

- constructor: `aRRR != BELFEM_QUIET_NAN` → `! std::isnan( aRRR )`
- `set_components( Cell< pair > )`: the RRR token compare is now
  case-insensitive, matching the constructor's detection
- `set_components( Cell, Vector, real )`: the fallback read of the stored RRR
  is guarded by `is_constant()`, so a genuinely missing RRR produces the
  intended `BELFEM_ERROR` naming the material in debug and release alike

## Status

Both translation units pass `g++ -fsyntax-only` with the tree's own flags
(`-Wall -Werror -pedantic-errors -std=gnu++17`). **Reviewed, not verified** —
no build and no solver run were performed in this session; the `tapestack3d`
example still has to be run to confirm the material path now completes.
