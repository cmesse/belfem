# Tape Example: User Material Recognized as Buffer — ABI-Stale Plugin

**Date:** 2026-07-28
**Purpose:** Root-cause analysis of why the `ybco` layer in the tape
example is tagged `DomainType::Buffer` instead of being treated as an
HTS conductor.
**Module:** physics/materials, fem/kernel (ThinShellFactory)
**Participants:** Claude (primary), Codex (audit), Grok (third voice)

## Symptom

`cmake-build-debug/tape/input.conf` defines a user material `ybco`
loaded from `MatData/build/libmat.so` (init symbol `hts_init`). At
runtime the ybco thin-shell layer is classified as a buffer layer, and
the `resistivity type : power-law` key is silently ignored.

## Root cause (confirmed by disassembly, git history, and both audits)

`libmat.so` was built 2026-04-07; the `Material` ABI changed afterwards
in two independent ways, both baked into the compiled plugin:

1. **Vtable drift (the decisive one).** The plugin registers jc, n, and
   rho through virtual calls to vtable slot 22 (offset 0xb0). In the
   April layout, slot 22 was
   `set_user_defined_function(Property, Dependency, MatFunc1*)`. Eleven
   virtual slots were inserted before it since then (`set_bh_curve`
   became virtual in `ae3a6754`; six T-derivative virtuals in
   `b524fe08`; four B/beta-derivative virtuals in `82a1efc1`), so
   slot 22 now resolves to `mu(const real H, const real T)`. The three
   registration calls execute `mu()` with garbage arguments (the enum
   ints travel in integer registers, `mu` reads uninitialized xmm
   registers) and discard the result — **nothing is registered, with no
   error**.
2. **Enum drift (secondary).** The plugin baked in `ec=30, jc=31, n=32`
   (April values). `n_bloch_gruen = 27` was inserted 2026-07-14,
   shifting `layer_thickness=30, ec=31, jc=32, n=33`. The plugin's one
   surviving direct call, `set_constant(30, 1e-4)`, therefore sets
   `layer_thickness` instead of `ec`. Note: enum drift *alone* would
   not produce the silent-Buffer symptom — `rho=6` never moved, and a
   healthy vtable would have hard-errored in `set_custom`'s default
   case on the mis-mapped property (Grok's refinement).

Downstream: `ThinShellFactory::create_buffers()` tags a layer block
`DomainType::Buffer` iff its material lacks
`have(MaterialProperty::rho)` — the sole criterion
(`cl_ThinShellFactory.cpp:1758-1762`, blocks default to
`DomainType::ThinShell` at `:261`). With nothing registered, ybco fails
the check and becomes a buffer. Likewise `have(jc)` is false, so
`MaterialFactory` (`cl_MaterialFactory.cpp:186`) skips the
user-superconductor branch and the power-law key is ignored.

The plugin is stale because the current `matlib.cpp` cannot compile:
line 86 reads `mat->set_type()` — no such member exists on `Material`
(the API is the getter `type()`), and the statement lacks a semicolon.
Source mtime 2026-04-30 postdates the .so build by three weeks.

## Fix for the example

Rebuild `MatData` against the current tree with a corrected `hts_init`.
Since the example's jc/n/ec are constants, register them as constants —
the 1-arg T-only `set_user_defined_function` path must NOT be used for
jc/n (see pitfalls):

```cpp
extern "C" void hts_init(Material* mat)
{
    mat->set_constant( MaterialProperty::jc, 4e10 );
    mat->set_constant( MaterialProperty::n,  25.0 );
    mat->set_constant( MaterialProperty::ec, 1e-4 );
    mat->set_user_defined_function( MaterialProperty::rho,
                                    MaterialDependency::T, &hts_rhon );
}
```

For field-dependent jc/n, use the 2-arg `(normB, angleNxB)` or 3-arg
`(normB, angleNxB, T)` overloads — these are the only paths that create
the `JcFunctionUserDefined` object consumed by assembly
(`cl_Material_UserDefined.cpp:93-118, 170-211`).

## Latent framework issues found (not fixed; open)

- `Material::set_custom` (`cl_Material.cpp:416-424`): case `jc` assigns
  `mFunctionRhoI = &Material::jc_custom`, case `n` assigns
  `mFunctionDebye = &Material::n_custom` — jc/n hijack the rho_i and
  debye function-pointer slots; both target functions are unconditional
  `BELFEM_ERROR` stubs.
- A plugin registering jc/n via the 1-arg T-only overload sets
  `have(jc)/have(n)` but leaves `mJcFunction`/`mNFunction` null;
  assembly then routes through `jc_eval`/`n_eval`
  (`powerlaws.hpp:98-115`) into the `constant_property` fallback, which
  is `BELFEM_QUIET_NAN` after `set_custom` — debug builds assert in
  `constant_property` (`cl_Material.hpp:1442`), release builds silently
  propagate NaN into the solve.
- **No ABI guard on user material plugins.** Any change to the
  `Material` virtual interface or the `MaterialProperty` enum silently
  breaks every previously compiled `.so`. Suggested cheap guard: an
  `extern "C" int belfem_material_abi_version()` symbol emitted by the
  user-material template and checked by
  `MaterialFactory::create_material(file, label)` at `dlopen` time.

## Diagnostic technique (reusable)

`nm -D libmat.so` (undefined symbols reveal which calls are direct PLT
vs virtual) + `objdump --disassemble=hts_init` (vtable offsets and enum
immediates are visible as literals) + counting `virtual` declarations
in declaration order (Itanium ABI: destructor occupies slots 0-1)
against the header at the .so's build date from git.
