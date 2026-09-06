# Devlog 2026-08-05 — Gas Models: Departure Convention, Own Fluid Property Tables, Deuterium

**Date:** 2026-08-05
**Topic:** Reference pressure convention of the real gas functions, documentation sweep, a
generator for the cryogenic NASA-9 and CEA transport correlations, and the tables BELFEM
now builds for itself
**AIs involved:** Claude, Codex, Grok
**Claude Confidence:** high on the generator, high on the convention analysis
**Codex Audit Confidence:** high
**Literature References:** NASA RP-1311; Chemical engineering thermodynamics §8.8, Eqs. 8.9,
8.14, 8.18, 8.19 (pages supplied by Christian)

## Summary

The session set out to fix defects in the gas models and ended up establishing what the
`*dep0` terms actually mean. They are **not** a defect. They are a deliberate convention,
and the code now says so. Separately, a generator for the low temperature correlations was
written and verified, and then grew into a build for the whole table: BELFEM no longer
redistributes the vendor CEA files but extracts what it needs from them. Deuterium was added
on top of that, for fusion fuel handling.

## The departure convention

`realgas_cp`, `realgas_h`, `realgas_s` and their derivatives subtract the departure
evaluated at `gPref` = 1 bar. Since `idgas_cp` is the bare heat spline, that subtraction
makes **SRK and PR reduce exactly to the ideal gas model at the reference pressure**. The
departure splines are rebuilt when the gas model changes so this holds for each. Helmholtz
does not use the path at all, taking caloric properties from the equation of state.

The key structural insight, from Christian: the subtraction would be the **complete and
correct** assembly if the heat spline held the *real* gas at 1 bar, because

    h( T, p ) = h_real( T, gPref ) + hdep( T, p ) - hdep( T, gPref )

is an identity. `RefGas::create_splines` sets `RefGasMode::POLY` and samples the CEA
polynomial, so the spline holds the *ideal* gas at 1 bar. Two self consistent designs exist
and the code combines the spline of one with the correction of the other:

| | spline holds | dep0 | |
|---|---|---|---|
| A | ideal gas at 1 bar | none | correct |
| B | real gas at 1 bar | subtract | correct |
| current | ideal gas at 1 bar | subtract | differs by the departure at 1 bar |

Model continuity and real gas accuracy cannot both be had: A and B give the true property
but produce a step when switching from IDGAS to a cubic at 1 bar, and that step is the real
departure from ideality, not an artefact.

Measured for nitrogen against VDI Heat Atlas D2.3 and NIST, which agree to 0.02 %: cp is
6.6 % low at 80 K, 3.0 % at 100 K, under 0.2 % above 250 K; cryogenic enthalpy differences
are off by 1 to 2.5 %. So the convention costs nothing where it is relied upon (combustion,
high temperature) and matters in the cryogenic range.

**Nothing was changed.** Terms were removed, audited, and restored. The decision is
Christian's and is recorded in `todo/gasmodels_open_source_migration.md`.

## Corrections to earlier reasoning, for the record

- The tabulated CEA state is the **ideal gas at 1 bar**. Describing it as "0 bar" is
  misleading even though cp and h are numerically the same at both, and it obscured the
  distinction that actually matters, which is ideal versus real *at* 1 bar.
- Equations 8.18 and 8.19 are **identities**, obtained by writing the departure definition
  at two states and subtracting. The second bracket appears only because a difference is
  being taken; it is not a reference offset. An earlier reading of them as licensing a
  "fixed reference state" subtraction was wrong and was withdrawn.
- The claim that the component path omitted the subtraction was false. It is carried inside
  `EoS_Cubic::hdep( aIndex, ... )` and `cpdep( aIndex, ... )`. Found by Grok, confirmed by
  Codex.
- The offsets `b1`, `b2`, `mSref` and `Sref/M` carry the integral below the start of the
  fit, because the polynomials do not reach 0 K. They are load bearing for the equilibrium
  constants and must not be touched.

## Defects fixed

- **D1** isobutane CAS `75-28-50` → `75-28-5`; the hardcoded CCR entry had never matched.
- **D2** eight spreadsheet mangled CAS numbers, in the data files **and** in the factory
  written against them.
- **D3** `uint` underflow computing the last field width, superseded by giving the dipole
  field a fixed width so a ninth provenance column can follow it.
- **D4** `std::string( std::getenv( "SCLS" ) )` was undefined behaviour whenever `SCLS` was
  unset, reached from every `Gas` construction.
- **D5** an unreachable alpha function fallback, verified dead by both auditors.
- **D6** a column overflow giving 2,5-lutidine a molar mass of 257.65, and a literal tab in
  the formaldehyde row that made its dipole parse as `.270` instead of `2.331`.
- Column bounds asserts on the five `Spline` overloads that take a column index, and a grid
  invariant assert on `Gas::spline_col`, whose cross spline reuse was correct but undefended.

## The low temperature generator

`scripts/fluidprop/nasa9_lowT.py`, with `scripts/fluidprop/README.md`.

Samples ideal gas cp from CoolProp, fits the NASA-9 form by equality constrained linear
least squares through the KKT system, and solves the two integration constants exactly for
continuity of h and s at the junction. Transport uses the same method against the dilute
gas limit, with several intervals when four parameters cannot span the range.

Two findings worth keeping:

- Fitting transport intervals sequentially, each constrained to the one above, **propagates
  edge error downward** and can be worse than a single interval; hydrogen conductivity went
  to 18.9 %. Solving all intervals simultaneously with continuity as constraints gives 2.8 %.
- The `--dilute-density` sampling matters: the CEA transport data describes the dilute gas,
  so samples are taken at low density rather than at a finite pressure.

Verified: all junction jumps at round off, and the generated records parse back through the
same reader and recover the right physical limits (5/2 R monatomic, 7/2 R for N2 and O2,
5/2 R for H2 at 14 K where rotation is frozen out).

Not in CoolProp: neon and xenon transport, and ozone entirely.

## Data work

`share/fluidprop/` now holds a 62 species `thermo.inp` and `trans.inp` built by
`scripts/fluidprop/build_tables.py`, a 59 species `gasdata.inp` generated from CoolProp with
per record provenance in three source columns, and a 26 species `cubicalpha.inp` restricted to
what the reader can reach. `gasdata_ppoc.inp` stays private.

The CAS joins were checked end to end: 15 species reach the PM alpha data, including the Hill
ordered labels `NH3` to `H3N` and `CH3OH` to `CH4O` and both hydrogen isotopes. The other 40
are ions, atoms and radicals with no critical point, which is correct rather than a gap.

## Own tables, and deuterium

The plan to ship the vendor `thermo.inp` was dropped in favour of **building our own**, which
resolves the RPA redistribution question by construction rather than by argument: the third
party additions in the circulating CEA file are all solid and propellant materials, so a table
restricted to the species BELFEM needs contains none of them. `scripts/fluidprop/build_tables.py`
extracts those species and merges the generated low temperature intervals into each record, so
there is one table per property, no overlay and no load order contract. 1.3 MB becomes 51 kB.

Deuterium was then added for fusion fuel handling: `D2`, `D`, `HD`, `OD`, `D2O` and the ions,
62 species in all. Four findings came out of it.

- **CoolProp refuses `Cp0molar` below the melting line**, although `Cp0` is a function of
  temperature alone and does not move across four decades of pressure. It is a guard on the
  state, not on the property, and it made deuterium unreachable below 19.74 K. The sampler now
  retries at a lower pressure, which is what makes the fit down to 19.2 K possible.
- **The heat capacity needed more than one interval for the first time.** Normal deuterium
  overshoots the classical 7/2 R, peaking at 3.62 R near 98 K before settling to 3.51 R, because
  the ortho and para forms of a boson pair weight the rotational levels differently than in
  hydrogen. One interval gives 4.66 %, two give 0.08 %. The fit now escalates the way the
  transport fit does, solving all intervals simultaneously with continuity as constraints. It
  fires only for deuterium; hydrogen is monotonic over the same range and one interval holds it.
- **CEA's deuterium transport is derived, not measured.** Its D2 viscosity is exactly sqrt(2)
  times its H2 viscosity to within 0.05 % at every temperature. NIST disagrees with CEA by 4.6 %
  in conductivity at the junction, so the 5.85 % residual there is a conflict between the
  sources, not a fitting failure. Continuity was kept: a step at 250 K would be worse.
- `D`, `HD` and `OD` have no transport in any database, so theirs is derived from `H`, `H2` and
  `OH` by the same mass scaling NASA used. Derived records say so in their header instead of
  inheriting the parent's citation, which would credit a measurement of hydrogen for data that
  is not one.

`D2O` gets no generated interval: heavy water freezes at 276.97 K, above the junction. Tritium
is in neither CEA nor CoolProp.

Verified: all interval edges continuous to 1e-9 or better, the deuterium fit tracks CoolProp to
0.06 % across 19 to 250 K, the nine previously fitted species are byte for byte identical to the
version before this change, and every file is pure ASCII within the fixed column layout.

Two presentation defects were then caught by Christian and fixed. The per species notes had
been written inline on the record lines, pushing the longest to 125 characters in a file whose
whole point is fixed columns; and one of those notes described **transport** on a record in
`thermo.inp`, which carries caloric properties only. Notes are now footnote references, `@1`
for a generated interval and `@2` for a mass scaled transport record, with the text and the per
species detail in the file header. Neither file comments on the other's property any more, and
the build warns if any line passes column 80 rather than assuming it did not. Verified: zero
numeric change to either table, 62 of 62 records reparse.

One fix worth recording: interval edges were built as `exp( linspace( log( a ), log( b ) ) )`,
whose ends are not exactly `a` and `b`. The junction constraint was therefore imposed a fraction
of a kelvin off the junction. Harmless for well conditioned fits, but CO2 is fitted over 217 to
250 K only and its coefficients moved in the third digit. The ends are now restored exactly.

## The move, and untangling the build flags

Christian moved `gastables` and `gasmodels` into `src/physics`; the CMake was adjusted to
match and `USE_GASMODELS` no longer implies `USE_NONFREE`.

- `src/physics/CMakeLists.txt` adds both modules under `USE_GASMODELS`.
  `nonfree/physics/CMakeLists.txt` is trimmed to `atmosphere` and a `combustion` gated on a
  new `USE_COMBUSTION`, which requires `USE_NONFREE` **and** `USE_GASMODELS` and says which
  one is missing.
- `BELFEM_NONFREE_SOURCE_DIR` became `BELFEM_SOURCE_DIR` in the two moved CMakeLists and in
  `atmosphere` and `combustion`, which reach across trees for the gas headers. Only
  combustion's self-include still points into the nonfree tree, which is right.
- Nothing outside `physics` includes a gas header, so no other module needed touching.

Six option combinations were configured. The four legal ones produce exactly the right
library set, and the two illegal ones fail with the intended message. The case that matters
is `USE_GASMODELS=ON` with `USE_NONFREE=OFF`, which was a hard error before and now builds
gastables and gasmodels alone — that is the open-source configuration the whole campaign was
for.

**This was also the first compile of any of this work**, and it found a real break.
`cl_GM_EoS_Cubic.cpp:636` called `Spline::update_data( matrix, values, gTref, Sref/M )`, but
the signature had since gained two `SplineBC` parameters ahead of the reals, so the two
reals were landing on the boundary condition arguments. `real` does not convert to
`SplineBC`, so this had been a hard compile error ever since — the module simply had not
been built. The sibling call at `cl_GT_RefGas.cpp:1141` already passed the full list and the
fix matches it.

With that fixed, `gastables`, `gasmodels`, `atmosphere`, `combustion` and the four example
executables all build and link with no errors, in a scratch build tree so the working one
was left alone.

**Behaviour change worth announcing:** combustion used to build unconditionally whenever
gasmodels and nonfree were both on. It is opt-in now, and existing build trees carry
`USE_COMBUSTION=OFF`.

## Reviving the gas tests

`tests/old/physics` was dormant and is now `tests/physics/{gastables,gasmodels}`, wired into
`make check` and labelled `fast`; the whole set runs in under a second. **15 of the 19 tests
pass.** Getting there meant fixing five real defects, none of which were test problems.

- **`Spline`'s two vector constructors were mutually ambiguous.** Both had every trailing
  parameter defaulted, so *any* three argument construction failed to compile. Nothing in
  `src` used that form, so it went unnoticed; the tests do. The legacy overload's `aXref`
  lost its default, which leaves four and five argument calls alone and sends the three
  argument form to the boundary condition overload. The two agree there: the legacy body is
  the other one with `NoCurvature` and zero slopes.
- **`create_glue_polys_heat` could not terminate on a smooth junction.** It searches for a
  half width whose glue polynomial has no sign change in curvature, and gives up at 50 K.
  Across the new handover junction there is nothing to repair — cp and its slope match to
  1e-10 because they are constrained to — so the search was hunting a sign change in round
  off. It had never been exercised before: the lowest junction used to be 1000 K, where a
  `tTmid >= 500` escape bypasses the test entirely. Junctions that are already smooth are
  now skipped.

  **Which species: the five noble gases and deuterium** — `Ar`, `Ne`, `He`, `Kr`, `Xe`, `D2`
  — established by sweeping all 68 through `gastable --gas` with the skip disabled. Nothing
  else in the table fails, though `N2`, `O2`, `CO`, `CO2`, `H2` and `CH4` all have the same
  250 K junction and glue across it successfully.

  For the noble gases the reason is exact: cp is constant, so the curvature the test examines
  is pure floating point noise. Helium's straddles zero at the 1e-8 level and argon's and
  xenon's underflow to identically zero, and the condition
  `( min < 0 && max < 0 ) || ( min > 0 && max > 0 )` excludes zero by construction, so it can
  never be satisfied however wide the window grows. Zero curvature is the ideal case for a
  glue polynomial, not a failure. Deuterium is the odd one out and its 250 K junction is what
  fails, not the 69 K one; it sits on the descending tail of the rotational hump, and why the
  quintic glue wiggles there has not been run down.

  **Correction to a claim made earlier in this session.** The threshold was justified as
  sitting between two populations "six orders of magnitude apart". That was measured on too
  small a sample and is wrong: deuterium's junction between its *two generated* intervals
  agrees only to 4e-6, above the smallest tabulated kink found (argon at 1000 K, 4.5e-7), so
  the populations overlap. The constraint is satisfied in the fit; evaluating either seven
  term polynomial at 69 K costs several digits to cancellation. What the 1e-7 threshold
  actually isolates is the single case of a fitted interval handing over to tabulated data,
  which agrees to 1e-10 or better. Everything else still goes through the search, including
  that deuterium junction, and succeeds there. The comment in the source now says this.
- **`Spline::update_data` hardcoded `SolverType::UMFPACK`.** Every build without SuiteSparse
  — including Christian's — failed there. It now uses SuperLU unconditionally, per his call:
  SuiteSparse is legacy, is not BSD-3 clean, and is therefore off by default and cannot be
  relied on. `initialize` still selects by `#ifdef`, but that chain prefers SuperLU first, so
  the two agree in practice.
- **`Gas::remix_heat` fills its heat spline through `matrix_data()`**, bypassing the paths
  that declare what the extra coefficient row holds, so `entropy()` tripped its assert even
  though the row had been mixed correctly from components that do carry entropy tables. A
  `set_extra_mode()` setter now says so.
- **Synthesized transport was unreachable through `mu()`.** For a species with critical data
  but no tabulated transport, `create_splines` synthesizes viscosity and then calls
  `idgas_lambda`, which needs viscosity back. The object is still in `POLY` mode at that
  point, so `mu()` routed to the empty placeholder polynomial and threw. It now switches to
  spline evaluation first, where both branches have left the viscosity.

Six species were added to the table for this: `Kr` and `Xe`, because `Gas::Gas()` names both
in its default air composition and a default constructed gas could not otherwise be built at
all; `N2H4`; and `BrF3`, `CF2ClBr` and `CH2Cl2`, which earn their place by shape — a
composition line holds five element/count pairs and nothing else in the table fills more than
three. 68 species now.

One expectation was stale rather than wrong: `BrF3`'s dipole is 1.19 D, and the private table
carried 1.1, having truncated instead of rounded exactly as it did for water and hydrogen
sulphide. Every other field of that record matches.

Four tests are parked with `DISABLED_` and a written reason rather than deleted or quietly
loosened. `Gibbs` returns NaN across a sweep that starts at 100 K, below the lowest tabulated
interval of every species in it, and its reference data was taken against a locally modified
`thermo.inp` that is no longer shipped — its own comment still says so. `AlphaFunction`
compares an analytic second derivative against the second derivative of a cubic spline, which
is only piecewise linear, and reaches r2 = 0.66 to 0.73 against a required 0.99.
`Cubic_Departure` and `Cubic_State` are the same shape of comparison and miss by one decade,
0.99999 against 0.999999; those two most likely need the tolerance revisited rather than the
code, but that is a call to make deliberately.

## Jury round on the ExtraMode change

Slug `review_gas_spline_extra_mode`, target `src/physics/gasmodels/cl_Gas.cpp`, mode `--jury`
(parallel, blind). Findings pre-registered before dispatch and left unedited.

**Both auditors independently confirm the diagnosis and the fix.** Codex: "I do not see a
defect in the reviewed heat-spline extra-mode path." Grok: root cause confirmed, and the
question that actually worried me — whether mixing row 4 linearly is legitimate when entropy
of mixing is not linear — resolves in favour of the change. The mixing term is not in row 4:
`update_mixture_entropy()` holds `-sum x ln x` and `idgas_s` applies it separately. Row 4
carries standard entropy, which is mole-fraction linear. Verdict rests on source trace rather
than on agreement, since all three traces are independent.

Two corrections to my own material, both from Grok, both accepted:

- My pre-registration cited `idgas_s` at `:1116`. That line is inside `realgas_s`; `idgas_s`
  is at `:963-973`. The claim held, the citation did not.
- "The data were always right, only the flag was missing" **over-claims**. It is right for the
  pure standard-state row. It says nothing about the absolute real-gas level: `realgas_s`
  adds `mSref`, which itself contains `entropy( gTref )`, and whether that double counts is
  the open absolute-anchor question. The passing tests cannot settle it either — they offset
  the literature table by their own `s( Tref, Pref )` and accept on r2, which is blind to a
  constant offset. **This change must not be read as closing the `*dep0` work.**

Codex found five pre-existing defects in `cl_Gas.cpp`, none introduced by the change, all
confirmed against source. Two are wrong arithmetic: `T_from_h`'s bisection fallback evaluates
`h( tT1 )` where it means the midpoint `aT`, so the bracket updates on the wrong sign; and
`shock()` seeds `aP2 = pow( tT2/aT1, ... )` with no `aP1 *` factor, giving a dimensionless
number where a pressure is wanted, with the same omission repeated in the low temperature
branch of `prandtl_meyer()`. The other three are `Gas` being copyable while its destructor
deletes raw pointees, public composition input validated only by `BELFEM_ASSERT` which
compiles out in release, and an element count held as `real` with a reuse guard that tests a
different quantity than it allocates. Full record and reconciliation table in the exchange
file until it is swept; no fixes applied.

## Open

- Which departure convention BELFEM should carry. The jury explicitly declined to treat the
  ExtraMode fix as bearing on it.
- The five `cl_Gas.cpp` defects above, `T_from_h` and the pressure seeds first.
- Whether `set_extra_mode` should be replaced by a single fill-and-declare call so the
  invariant cannot be separated from the data. Both auditors call the current shape an
  acceptable minimal fix and a residual footgun.
- The four disabled tests above, `Gibbs` first: NaN is a defect somewhere, not a tolerance.
- `atmosphere` stayed in `nonfree/physics`. It is ISA-1976 and has no data or export control
  problem of its own, so whether it follows the other two is a decision, not a blocker.
- The runtime data path has still never been exercised: no binary has been run.
- **`gasdata.inp` has no generator script.** It was produced from CoolProp by hand, so it cannot
  be reproduced or extended the way the tables now can. This should be closed before the move.
- **Acetone is ideal gas only, and cannot be fixed by data alone.** Its CEA label
  `C3H6O,acetone` is 13 characters against an 11 character label field in `gasdata.inp`, so no
  record can be keyed to it. `fn_GT_fix_label.cpp` already strips the descriptor for
  `C2H3,vinyl` and needs the same case for acetone. The factory carries a hardcoded CCR entry
  for its CAS `67-64-1`, so the code expects it to work.
- Nothing in this session has been compiled.

## Files Updated

- src/physics/gasmodels/cl_Gas.cpp, cl_Gas.hpp
- src/physics/gastables/fn_GT_data_path.{cpp,hpp}, cl_GT_InputData.cpp,
  cl_GT_InputAlpha.cpp, cl_GT_RefGasFactory.cpp, CMakeLists.txt,
  (belfem_gastables_config.hpp.in was introduced and then removed: the data path
  now comes from the gBelfemDataPath global, set from $BELFEM_DATA)
- src/physics/gasmodels/cl_GM_EoS_AlphaFunctionFactory.cpp
- src/physics/gasmodels/doc/README.md, nonfree/physics/combustion/doc/README.md
- nonfree/physics/tables/ vendor data files (inputs to the build, not shipped)
- src/numerics/spline/cl_Spline.hpp
- scripts/fluidprop/nasa9_lowT.py, build_tables.py, species.txt, README.md
- share/fluidprop/
- todo/gasmodels_open_source_migration.md

Later in the session, for the move, the build flags and the tests:

- CMakeLists.txt (USE_COMBUSTION, the check/check-fast fan-out)
- src/physics/CMakeLists.txt, nonfree/physics/CMakeLists.txt, nonfree/CMakeLists.txt
- src/physics/{gastables,gasmodels}/CMakeLists.txt,
  nonfree/physics/{atmosphere,combustion}/CMakeLists.txt
- src/physics/gastables/cl_GT_RefGas.cpp (glue skip, synthesized transport mode)
- src/physics/gasmodels/cl_GM_EoS_Cubic.cpp (stale update_data call)
- src/numerics/spline/cl_Spline.{hpp,cpp} (ctor ambiguity, set_extra_mode, SuperLU)
- tests/physics/CMakeLists.txt and the revived tests/physics/{gastables,gasmodels}/;
  tests/old/physics removed
