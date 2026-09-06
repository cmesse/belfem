# Copper thermal conductivity ignores the field: dependency flag wiped by `set_RRR`

**Date:** 2026-09-03
**Purpose:** Root-cause a user report that `material copper --rrr 50 -b 0` and `-b 20` print identical `lambda(T)` columns while `rho` differs
**Module:** `src/physics/materials`, consumer in `src/fem/kernel`

## Symptom

Both tables show the same thermal conductivity at every temperature; resistivity moves from 0.031 to 0.054 ×1e-8 Ω·m at 4 K. Physically the Wiedemann–Franz coupling demands that λ drops with the magnetoresistance.

## Finding

The physics is implemented and correct: `Metal::lambda(T,B,β) = λ(T)·ρ(T)/ρ(T,B,β)` (`cl_Material_Metal.hpp`). It is never selected. Both consumers of thermal conductivity — the `material` tool (`main.cpp`, `lambda_at`) and the FEM calculator (`cl_FEM_Calculator.cpp`, `tLambdaFieldDependent`) — choose the three-argument overload only if `depends( lambda, normB )` is set, and for a constructed pure metal with a finite RRR it is not:

1. `Copper::create_kohler()` → `Metal::set_kohler_dependencies()` sets `{T, normB, angleBxJ}` on both `rho` and `lambda`.
2. `Copper` then calls `Metal::set_RRR()`, which rebuilds the lambda spline through `SplineLookupTable::create_spline` → `set_spline`, whose first act is `mPropertyDependencies( tIndex )->reset()` followed by `set_dependency( lambda, T )` only.
3. `rho` keeps its flags because `set_RRR` does not re-spline `rho`.

The bitset reset was introduced with the `SplineLookupTable` refactor (bc4a525c, 2026-07-02); the spline rebuild in `set_RRR` is older (present at 72c03c9f, 2025-11-18). Every pure metal (Aluminum, Silver, Nickel, Iron, Chromium, Indium, Lead, WhiteTin, Copper) reaches the same wipe; alloys restore the flags themselves after their table build (`cl_Material_Alloy.cpp`) and are unaffected; YBCO sets RRR before it declares lambda and is unaffected.

**Where it sits:** the material module. The FEM side obeys the flag correctly, so it is a victim: for every field-aware pure metal the coupled thermal solve uses `compute_lambda_bulk` (λ(T)) and `compute_dlambdadT_bulk`, i.e. the stabilizer conductivity ignores the magnetic field while its Joule heating does not.

## Evidence

Scratchpad probe linked against the prebuilt `cmake-build-debug/lib/libbelfem.a`, no `make`:

```
depends(rho,normB)=1 depends(rho,angleBxJ)=1
depends(lambda,T)=1 depends(lambda,normB)=0 depends(lambda,angleBxJ)=0
     T    lambda(T)  lam(T,20,0) lam(T,20,90)     rho(T)e8 rho(T,20,0)e8
   4.0      310.613      182.212       88.965       0.0314       0.0536
  20.0     1350.243      796.485      393.829       0.0323       0.0547
  77.0      510.860      470.991      390.209       0.2256       0.2447
 293.1      396.980      395.807      385.442       1.7081       1.7132
```

λ(4 K)/λ(4 K, 20 T, 0°) = 1.70 = ρ ratio, as the scaling demands. Codex (gpt-5.6-terra, medium) and Grok (grok-4.6, high) both confirmed the chain independently; Codex corrected the commit history, Grok mapped the other wipe paths and the YBCO/Alloy side cases. Every citation of theirs that was checked held. Exchange: `tmp/ai_exchange/copper_lambda_field_dependence.md`.

## Fix (applied, on Christian's approval — "we need all Metals to work properly")

`Metal::set_RRR` now ends with a gated `set_kohler_dependencies()`: when `depends( rho, normB )` is
true — the marker that the metal carries a magnetoresistance model — the lambda flags wiped by the
spline rebuild are re-registered. One edit in the base class covers all nine pure metals. `set_spline`
stays T-only by contract (alloys rely on the wipe and restore their own flags), the calculator is
untouched, and YBCO is unaffected because it calls `set_RRR` before declaring lambda, so the gate does
not fire there.

Executable evidence, again a scratchpad probe with the edited translation unit compiled under the
library's own flags and linked ahead of the prebuilt archive, all nine metals through
`MaterialFactory::create_material( label, 50.0, false )` at 4 K, 20 T:

```
material    rhoB  rhoA  lamB  lamA |     lam(4)  lam(4,20,0) lam(4,20,90) |   rho(4)e8 rho(4,20,0)e8
copper         1     1     1     1 |    310.613      182.212       88.965 |     0.0314       0.0536
aluminum       1     1     1     1 |    197.332      132.775       68.995 |     0.0495       0.0736
silver         1     1     1     1 |    326.354      180.678       95.250 |     0.0299       0.0540
nickel         1     1     1     1 |     93.946       70.506        9.513 |     0.1041       0.1387
iron           1     1     1     1 |     62.229       34.430        3.995 |     0.1570       0.2838
chromium       1     1     1     1 |     39.390       30.462        8.857 |     0.2477       0.3203
indium         1     1     1     1 |     59.269       42.612        7.520 |     0.1631       0.2269
lead           1     1     1     1 |     24.582       16.276       10.602 |     0.3942       0.5954
tin            1     1     1     1 |     47.316       12.036       24.381 |     0.2062       0.8106
```

Not examined: tin's transverse conductivity exceeds its longitudinal one at 20 T, which is a property
of the WhiteTin Kohler curves, not of this change.

**Regression test:** `tests/physics/test_MetalFieldDependence.cpp`, registered in
`tests/physics/CMakeLists.txt`. Two cases over the nine metals: the flags survive `set_RRR`, and
`lambda( T, B, beta ) · rho( T, B, beta ) = lambda( T ) · rho( T )` with `lambda( T, B, beta ) < lambda( T )`.
Built as a gtest probe against both libraries: passes with the fix, and the flag case **fails** against the
unfixed prebuilt archive — the test discriminates. Not yet run through `make check`; the shared tree was
not built.

**Docs:** one blockquote paragraph added to §6 of
`src/physics/materials/doc/materials_contracts_and_invariants.md` recording the wipe-and-restore
contract. Static gates passed on the edited TU: Blaze debug compile with `-Werror`, Armadillo
`-fsyntax-only` with the tree's flags.

**Owed:** `make check` (physics suite) after the next build; rerun the two `material` commands; a coupled
deck with a copper stabilizer to see λ follow the field.

## Second finding: the lambda spline's start tangent is the reciprocal of the slope

Asked whether the resistivity parameters are linked correctly, three links were probed on copper RRR 50.
`rho( 273.15 ) / rho( 0 )` is 50.000000; `lambda( T, 0, β )` agrees with `lambda( T )` to 2e-5 (the
Kohler zero-field branch returns the analytic ρᵢ, the plain accessor reads the ρᵢ spline). The third link
is wrong: `Metal::set_RRR` passes `x / constant::L0` = ρ₀/L₀ as the start tangent of the λ spline, and
the spline's `Tangent` boundary row pins the *first derivative*. The Hust curve's slope at 0 K is
L₀/ρ₀ (77.7 W/(m·K²) here), so the spline is pinned to a near-zero slope at the origin and rings through
the first knots. Spline vs analytic: −79 % at 0.5 K, −32 % at 2 K, exact at the 4 K knot, +3.3 % at 5 K,
+2.8 % at 6 K, −0.45 % at 10 K, negligible above ~14 K. With the reciprocal the error is below 0.01 %
everywhere. The 5 K and 6 K entries of the reported table are the wrong-BC numbers. The same expression
sits at `cl_Material_YBCO.cpp:123`; what the correct slope is there depends on the superconducting
electronic term at T → 0 and was not examined. Not fixed — reported for a decision. Cosmetic: the tool's
resistivity column is labelled `1e-8 A/m²`; the unit is Ω·m.

## Round 2: spline tangent and unit label fixed (Christian: "Of course, let's fix the unit label too")

`Metal::set_RRR` now passes `constant::L0 / x` as the start tangent of the λ spline; the YBCO
constructor passes `L0 / rho_0` for the same reason — its `lambda_custom` is `L0·T/(ρᵢ+ρ₀) + κ_ph(Callaway)`,
whose electronic term has slope L₀/ρ₀ at the origin while the phonon term starts flat. The `material` tool's
seven resistivity headers now read `1e-8 Ω·m` (was `A/m²`, and three of them `10-8`); display width kept.

A third gtest case, `LambdaSplineMatchesHustBelowFirstKnots`, compares the spline against the
`Metal::lambda_custom` expression at 0.5–6 K for all nine metals inside a 1 % band. With the fix seven
metals sit below 0.1 %; indium and lead reach 0.6 % and 1 % *between* the 4 K knots, which is knot
resolution, not the tangent (lead still misses by 1 % at 10 K, outside the test window). Against the
unfixed archive the case fails (copper −79 % at 0.5 K). Static gates on the edited TUs: Blaze debug compile
with `-Werror`, Armadillo `-fsyntax-only`, both clean. YBCO not exercised at run time (its construction
needs the data path); the change there is the same expression and rests on the derivative argument above.

## jc / n wiring check (Christian: "double check that the jc and n functions are wired correctly")

Read-only trace plus one probe. **The solve path cannot suffer the metals' mistake:** `rho_powerlaw( normJ, T, normB, angleNxB )` calls `jc_eval` / `n_eval` (`powerlaws.hpp`), which route on the *JcFunction's own* dependency bitset (`mJcFunction->depends_on( T )` picks the 3- or 2-argument `eval`) with a `constant_property` fallback when no function is attached. The material-level `depends( jc, … )` flags are consulted by nobody on the assembly path; their only consumer is the J/Jc post-processing branch for `UserDefined` materials. Each JcFunction sets its own flags (`JcFunctionDatabase`: normB, angleNxB, T; `JcFunctionModifiedKim`: normB, angleNxB; `JcFunctionUserDefined`: per overload arity). Argument order is consistent end to end: calculator `( normJ, T, normB, beta )` → `jc_eval( T, normB, angle )` → `eval( normB, angle, T )` → `Database::evaluate( T, log10 B, wrapped angle )`; bulk HTS passes the π/2 dummy beta, thin shells the unfolded `bn_angle`. Nothing in construction re-splines or `set_custom`s jc or n, and the one `reset_dependencies` (UserDefined 3-argument path) re-sets all three flags in the next lines. No spline-tangent analogue: the tables are tensor databases, not clamped splines.

**Latent defect found (not fixed):** `Material::set_custom( jc )` writes `mFunctionRhoI = &jc_custom` and `set_custom( n )` writes `mFunctionDebye = &n_custom` — copy-paste of the `rho_i` / `debye` cases (`cl_Material.cpp`, `set_custom` switch). No built-in material calls them, but `UserDefinedMaterial::set_user_defined_function( jc | n, MaterialDependency::T, MatFunc1 )` does, unguarded. Probe (Copper subclass calling `set_custom( jc )`, `set_custom( n )`): `rho_i( 77 )` then throws "jc_custom not implemented", `debye( 77 )` throws "n_custom not implemented", and `jc_eval` / `n_eval` throw "Property is not constant" (NaN in a release build) because no JcFunction is attached. The usage guide already says the `jc_custom( T )` path is not for user superconductors; the API does not enforce it. Proposed: reject jc and n in the one-argument overload with a `BELFEM_ERROR` pointing at the JcFunction route, and turn the two `set_custom` cases into errors. The post-processor's `jc_custom( Temp )` fallback becomes unreachable and can go.

## Round 3: the latent jc / n callback defect fixed (Christian: "Let's follow your proposed fix")

Three edits. `Material::set_custom` now rejects `jc` and `n` with a `BELFEM_ERROR` **before** it raises the
have-flag, so a rejected call leaves the material untouched; the two copy-paste cases that wrote into the
`rho_i` and `debye` dispatch slots are gone. `UserDefinedMaterial::set_user_defined_function( …, T, MatFunc1 )`
rejects `jc` and `n` with a message naming the field-dependent overloads. The J/Jc post-processing branch for
user materials branches on `depends( jc, T )` alone, since every user jc is now a JcFunction of
`( normB, angleNxB[, T] )`; the `jc_custom( Temp )` fallback is gone from both sites (bulk and thin shell).
Docstring of the one-argument overload and §7 of the usage guide updated. New gtest
`tests/physics/test_JcCustomGuard.cpp`: `set_custom( jc | n )` throws, `have( jc )` stays false, `rho_i` and
`debye` still evaluate — passes with the fixed TUs linked ahead of the archive, fails against the unfixed archive.
Gates clean on both backends for `cl_Material.cpp`, `cl_Material_UserDefined.cpp`, `cl_MaxwellPostprocessor.cpp`.

**Kept on purpose:** `Material::jc_custom` / `n_custom` and their `UserDefinedMaterial` overrides. I removed them
as dead code first and the compile refuted that: the reduced overloads `rho_powerlaw( normJ, T )`,
`rho_piecewise( normJ, T )` and their defect and derivative twins in `powerlaws.hpp` still call them. Whether
those overloads have any caller is a question for the audit round; removing them is out of scope here.

**Incident, my error.** To undo the dead-code removal I ran `git checkout -- src/physics/materials/cl_Material.hpp`
without reading the diff first. The tree was being edited concurrently by another session's doxygen-contradiction
sweep (appliers across ~90 files), and that checkout discarded its uncommitted comment edits to the header
along with mine: 16 insertions, 7–8 deletions. Recovered in full from the sweep's own records — the applier
log `tmp/ai_exchange/doxygen_sweep_reports/materials_applied.md` names the findings (F11 ×3, F37, F107 ×2;
F01–F06 were already in HEAD), and the subagent transcripts hold the verbatim old/new texts — and re-applied;
the header diff is again +16 lines. Rule broken: look at the target before overwriting, and `git status`
before edits in a shared checkout. Rule to keep: never `git checkout -- <file>` on this tree; revert by
inverse edit of the hunk you own.

Round-3 audits (Codex terra/high, Grok 4.6/high): ship. Both refuted my "every user jc is a JcFunction" —
`set_constant( jc )` is public on a user material too — and both confirmed the simplified branch survives
that case through the constant fallback in `Material::jc()`. Comments and the test adjusted accordingly
(`set_custom` is public; the test now also locks `have( n )`). Both confirm the reduced `( normJ, T )`
power-law family has no caller anywhere in the tree; with the T-only registration gone, those overloads and
`jc_custom` / `n_custom` are an unpopulatable trap. Deleting them is a decision for Christian; the usage guide
now says not to call them. Left as found (Grok): the 2-argument user jc/n overload does not validate that its
dependencies are normB and angleNxB.
