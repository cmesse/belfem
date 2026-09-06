# Three-AI Review of the Material-Model Commit + Iron Kohler Refit

**Date:** 2026-08-24
**Topic:** Jury cross-review of `acd4bbf0` ("more work on material model"), then the Iron Kohler
rebuild it triggered: transverse channel refit as a flat-ended Bezier, fixes applied on approval.
**Participants:** Christian, Claude (Fable), Codex, Grok
**Exchange record:** `tmp/ai_exchange/review_material_cryo_kohler.md` (pre-registration,
both audits, verification pass, reconciliation table)

## 1. Review round (jury mode, read-only)

Target: the full diff `200f0cba..acd4bbf0` — cryo-expansion refactor
(`create_alpha_table()` per-material duplicates folded into
`SplineLookupTable::create_cryo_expansion()`), 5-segment cp extension + `d2cpdT2_custom`,
new Debye curves (Al, Cr), new Kohler magnetoresistance (Al, Fe, Ni), new λ/ρ data,
`mTAlphaSwitch` → `mTAlphaPlateau` rename.

All three reviewers converged independently on the same core findings; every citation was
re-verified against source. Highest evidence level anywhere in the round: source trace
(protocol §11 rung 5) plus one build-tree probe. Reviewed, not verified.

**Held up cleanly:** the cryo-expansion fold, the 5-segment cp bookkeeping, the
`d2cpdT2_custom` chain rule (re-derived by hand, twice independently), Silver/Indium's
θ₀-from-β derivation, Nickel's Kohler implementation, the rename, all constructor reorders.

**Confirmed P0 (all latent — no live path until wired):** `Iron::create_kohler()` was a
copy of the Aluminum version with three paste errors (`mKohlerTransPolys(2)` written before
`set_size`; `q1` aliased `mKohlerLongPolys(1)` instead of `mKohlerTransPolys(1)`, leaving
the transversal mid-field poly empty and clobbering the longitudinal one;
`mKohlerLongPolys(2)` read in the high-field branch but never written) — and the function
was never called from the constructor, while the `kohler()` override had already removed
`Metal::kohler`'s `BELFEM_ERROR` backstop. Iron is factory-registered; the FEM path was
safe only because the unset `depends(rho, normB)` bit gates `cl_FEM_Calculator` off.

**Confirmed P1:** Aluminum (and the Iron copy) computed the longitudinal Kohler saturation
as `exp( polyval( p1, x1 ) )` with `x1 = exp(w1)` — the polynomial lives in `w = ln(B·S)`,
so the value underflowed to exactly 0 instead of ≈0.57, discontinuous at the switch.
Confirmed numerically by all three reviewers independently.

Christian triaged: the paste errors and the mattest findings are trivia (mattest is a
scratch stub, not shipped); the real question is the Iron curve model itself.

## 2. Iron Kohler refit

Christian supplied the raw Kohler data (ln(B·S) vs ln(Δρ/ρ_0T), 37 longitudinal points
from three series, 19 transverse). Findings from the fit work:

- **The committed quadratics ARE the least-squares optima** — `polyfit` reproduces both
  `p1` and `q1` to all twelve digits. Nothing better exists in the quadratic family; the
  defect was structural.
- The transverse quadratic, fitted on data ending at w = 6.16, has its vertex at
  w = 11.17 — inside the gap before the shared saturation knot w₁ = 14.21 — so the
  tangent-line tail anchored at w₁ (the established convention: longitudinal vertex sets
  the knot for both channels) inherited slope −0.63: transverse MR falling with field.
- An intermediate proposal (vertex-constrained quadratic, pinned to w₁) was superseded by
  **Christian's design: take the peak from the transverse dataset itself (the vertex of
  its LSQ quadratic: w = 11.1728, Δρ/ρ = 95.36) and fit a flat-ended Bezier there, in the
  Copper/Nickel idiom.** `Bezier::xi_by_x` clamps past the last control point
  (`cl_Bezier.cpp:97-98`), so the mid-field branch evaluates the Bezier up to w₁ and it
  holds the peak automatically — no transverse tail polynomial at all.
- Fitted control net (x₀ = 3.0 fixed at the low-field switch, end pinned at the peak,
  y₂ = y₃ for the flat end, y₁ ≤ g_peak enforced after the free fit overshot to 5.39):

  ```
  x = { 3.0,       7.132373, 10.138093, 11.172802 }
  y = { -2.324512, 4.549139,  4.557684,  4.557684 }
  ```

  rms 0.0568 vs 0.0565 for the unconstrained quadratic — no measurable fit cost.
- Full piecewise model validated in Python (same Hermite low-field construction as the
  C++): C⁰ at both knots to ≤1e-9, monotone in both channels, transverse saturation 95.4
  above longitudinal 34.6 everywhere — a more physical picture than either the fold-over
  or the 706 the vertex-pinned quadratic implied.

**Honest caveat, logged for the future:** the longitudinal saturation (34.6 at w₁ = 14.2)
is an extrapolation from data ending at w = 9.8, and the three measurement series fitted
separately put it at 1.6, 5.2, and 93.5 — the pooled value is dominated by the 24-point
series. Inside the realistic operating envelope (w ≲ 11.5) the curve is data-anchored.

## 3. Fixes applied (on explicit approval)

- `cl_Material_Iron.cpp` — `create_kohler()` rebuilt: longitudinal quadratic + saturation
  at its vertex (`polyval(p1, w1)`, the coordinate fix), transverse flat-ended Bezier +
  constant `mKohlerTransPolys(1)` above w₁; container layout now
  `mKohlerLongPolys{cubic, quadratic, saturation}`, `mKohlerTransPolys{cubic, saturation}`;
  all three paste errors gone by construction. `kohler()` mid/high-field branches updated.
  `create_kohler()` wired into the constructor (before `set_RRR`, so the ρ-database build
  sees the dependencies). Destructor deletes the new `mKohlerTransBezier`.
- `cl_Material_Iron.hpp` — `Bezier * mKohlerTransBezier = nullptr` member.
- `cl_Material_Aluminum.cpp` — `polyval( p1, x1 )` → `polyval( p1, w1 )` in the
  longitudinal saturation (one token + comment).

Hygiene sweep, second approval, same session:

- `mattest.cpp` — stray `#include <dense/DenseMatrix.hpp>` removed. It resolved only
  through `$SCLS/include` (confirmed: `find_scls.cmake` is gated on `DEFINED ENV{SCLS}`),
  and mattest is an unconditional CMake target — a build-breaker for any non-SCLS machine.
- Stale class headers rewritten: `cl_Material_Aluminum.hpp` and `cl_Material_Chromium.hpp`
  both claimed "the class does not construct yet" while their constructors run; now they
  state what is fitted, that only `create_mech()` is open (commented out), and for Cr that
  no Kohler curve exists yet.
- Unused `fn_create_beam_poly.hpp` includes removed from Chromium AND Iron (the Iron one
  was flagged by clangd during today's edits — same class as the review's Chromium find).
- DOI typo `0.1098/rsta.1959.0004` → `10.1098/…` fixed at both Chromium sites.
- `doc/thermal_expansion_from_heat_capacity.md` — `cp_from_debye` citation `:512` → `:537`.
- EOF newlines added to `cl_Material_Nickel.cpp` and `cl_Material_SplineLookupTable.cpp`;
  double semicolon in `cl_Material_Aluminum.cpp` removed.
- `Metal::create_cp` hardened: the re-entry guard now also checks `mCpBezierMedium`, and
  the `n == 5` block guards `Rx` against fold-back. Deliberately NOT the strict
  control-polygon check used for `Px`/`Qx` — Iron's R polygon legitimately dips
  (7.0214 → 7.016) while x(t) stays monotone — but the sharp Bernstein condition on the
  control-point differences: x'(t) ≥ 0 on [0,1] iff d0 > 0, d2 > 0, d1 > −√(d0·d2).
  Iron passes with d1 = −0.0054 against a bound of −0.154; Chromium's R is strictly
  increasing. A future refit that folds now aborts loudly at construction instead of
  silently mis-inverting in `xi_by_x`.

Deliberately NOT touched (design call / scratch, per triage): the Al/Ni
`set_kohler_dependencies()` wiring question, the exact-float `Qx(3)==Rx(0)` connection
contract, mattest's Chromium-value rho_ref and mid-main `exit(0)`, and clangd's
transitive-include warnings elsewhere (several are false positives, e.g.
`fn_create_beam_poly.hpp` in `cl_Material_Metal.cpp` is genuinely used).

## 4. Chromium Kohler (added by Christian later the same day, reviewed + validated)

Chromium got its Kohler curves: cubics in w for both channels (Kozlova & Kondorskii 1963
longitudinal, Arajs & Dunmyre 1965 transverse), with a new switch design — low-field
switch at the longitudinal cubic's inflection point (w = 4.274, BS ≈ 72), saturation at
its local maximum (w₁ = 6.544, BS ≈ 695, g′(w₁) = 0 to machine precision, so continuity
holds by construction). Review against the Iron-bug checklist: all four traps avoided
(sizing order, aliasing, unwritten vectors, w-vs-BS coordinate). Numeric validation of
the full piecewise model: C⁰ at both knots, monotone and non-negative in both channels,
longitudinal saturating at 0.293, transverse at 3.32 with a gently rising tangent tail
(slope +0.104 — the transverse maximum at w = 6.697 sits just above the shared knot);
transverse above longitudinal everywhere.

Cleanup applied on approval: the class doc's "No Kohler magnetoresistance curve yet"
sentence replaced (stale within a day of being written); the three Kohler members moved
from the public section (accidental placement) to the private data block; citation
comments corrected — "Arays and Dunmore, 1964" is actually **Arajs and Dunmyre, 1965**
(J. Appl. Phys. 36, 3555, verified against the DOI landing page), plus
Soviet/Strength/Kondorskii typos in the JETP line.

## 5. Status

Reviewed and numerically validated in Python; **not verified** — no build was run
(shared build tree, user-runs-builds convention). Executable gate owed at next rebuild:
construct `material::Iron` with RRR set (the constructor now builds the ρ database when
`aBuildTables` and the Kohler dependencies are set) and spot-check `kohler()` continuity
at both switch knots, e.g. via mattest. Since Iron now takes the `depends(rho, normB)`
path in `cl_FEM_Calculator`, the first FEM run with iron in field is part of the gate.
Elastic data (create_mech for Al/Cr) remains open work, per plan.
