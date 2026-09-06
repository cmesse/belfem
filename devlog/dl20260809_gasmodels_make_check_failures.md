# Gas Modules: `make check` Failures Root-Caused

**Date:** 2026-08-09
**Purpose:** Record the diagnosis and fix of the gasmodels failures in `make check`
**Module:** `src/physics/gasmodels`, `tests/physics/gasmodels`

## Summary

Two rounds, five failures, and **two of them are genuine physics defects**: a sign typo in
methane residual term 34, which breaks the Setzmann & Wagner equation of state in the critical
region only and survived the 2026-08-06 literature campaign; and the cubic entropy reference,
where commit `1c04933f` removed a term it read as a double count but which was in fact
cancelling a standard-state entropy baked into the departure spline. The remainder are
test-side: an overload ambiguity that stops the Gibbs test compiling, a Prandtl-Meyer premise
that assumed argon is calorically perfect over a wider temperature range than the NASA table
delivers, and a finite-difference step too coarse for a pseudo-critical cp spike.

Method note: CoolProp 8.0.0 is installed in the system python and carries the same published
equations (Setzmann & Wagner for methane, Lemmon for air). Used twice here as an independent
oracle — once to find the methane coefficient by dumping its EoS JSON and diffing term by term,
once to prove a failing test was measuring its own truncation error rather than a BELFEM defect.

## 1 — Methane: residual term 34 sign (CRITICAL, live defect)

**Symptom.** `GASMODELS.Methane_Cp_FD_NearCritical` aborts inside `Helmholtz::v()`:

```
Too many iterations for T= 190.595 K, p= 50.000 bar, rho= 107.217 kg/m^3, relax= 0.264
```

The density Newton is not fragile — it has no unique root to find.

**Cause.** `cl_GM_EoS_Methane.cpp:90` carried `mN(33)` = `+1.387292044e-2`. Setzmann & Wagner
1991 (J. Phys. Chem. Ref. Data **20**, 1061), Table 35, term 34 (d = 4, t = 22, c = 4) is
**−0.1387292044·10⁻¹**. It is the only one of the 40 residual terms that disagrees with the
published set; every n, d, t, c, α, β, γ and the Gaussian centres are otherwise exact.

**Consequence.** Eq. (5.3) violates the critical constraints it was fitted under (§5: first and
second partial derivative of pressure with respect to density zero at the critical point):

| quantity at (Tc, ρc) | as coded | required |
|---|---|---|
| p | 4.599174 MPa | 4.5992 MPa |
| ∂p/∂ρ | −16129 Pa·m³/kg | 0 |
| ∂²p/∂ρ² | −496 | 0 |

p was right because term 34 is even in nothing useful — the constraint that fixes p at one
point survived while the derivatives did not. Every isotherm up to ≈ 197 K therefore carried a
spurious van der Waals loop: at Tc the pressure runs over 3.94…4.82 MPa across
135 < ρ < 206 kg/m³, so at 5 MPa there are three roots and a stretch of ∂p/∂v > 0 between them.
Against the paper's own ancillary equations (3.2) vapour pressure and (3.4)/(3.5) saturation
densities, the coded EoS was exact below ~170 K and off by −4.5 % in p at 180 K, −16 % at 190 K.
Hence the invisibility: the whole D1 campaign validated the far field.

**Verification.** A numerical mirror of the coded coefficient set was checked against CoolProp
8.0.0 (whose methane is the same Setzmann & Wagner equation). Before the fix, agreement was
~1e-5 everywhere except the critical region, where it reached 14 %. After the fix, agreement is
**5.7e-6 uniformly across liquid, vapour, near-critical and supercritical states** — and that
residual is not an error but BELFEM's `Rm/M` = 518.2653 J/(kg·K) against the paper's rounded
518.2705. The critical constraints now hold (∂p/∂ρ and ∂²p/∂ρ² vanish to round-off), isotherms
above Tc are monotone, and ρ(190.595 K, 5 MPa) = 234.633 kg/m³, matching CoolProp to 5e-6. The
Newton in `Helmholtz::v()` was never the problem.

Tracked as **D16** in `../todo/closed/gas_correctness_fixes.md`.

## 1b — Cubic entropy: `mSref` was compensation, not a double count (CRITICAL)

Surfaced by the first run after the fixes above: `GASMODELS.Cubic_Entropy` fails with
r2 = −56.4 for PR **and** SRK — an RMS residual of 6648 J/(kg·K) against VDI D2.2 Table 7.

**Cause.** The departure spline is anchored on the component's *standard* entropy
(`cl_GM_EoS_Cubic.cpp:628`, `update_data( ..., gTref, data(g)->Sref()/data(g)->M() )`), so
`sdep0(T)` — an entropy *departure* — carried +Σ yᵢsᵢ°/Mᵢ = 6702 J/(kg·K) for air.
`Gas::realgas_s` subtracts `sdep0`, and the `+ mSref` term retired by commit `1c04933f` was
exactly what cancelled it: `mSref` = `idgas_s(gTref,gPref)` = 6702 (spline) + 162 (mixture).
The commit added the mixture term explicitly — correct — but dropped the 6702 cancellation with
it, so every cubic gas returned an absolute entropy one standard state too low. The residual
matches: 6702 minus the ~50 J/(kg·K) by which the cubic models genuinely miss the VDI table.

**The test data is sound**, checked independently: CoolProp's air EoS reproduces VDI Table 7
(relative to 298.15 K, 1 bar) to +50…65 J/(kg·K), so the test's `tS += tAirID.s(298.15,1e5)`
construction is the right anchoring. Note this also means the test has only ~0.005 of r2 budget
against ~60 J/(kg·K) of honest model-vs-table difference — it sits close to its tolerance by
design.

**Fix.** Anchor the departure spline's entropy at 0.0, so `sdep0` means what its name says and
the cleaned-up `realgas_s` becomes the pre-`1c04933f` expression term for term. The
cancellation is exact, not approximate: `spline_entropy(gTref)` is the same mass-weighted
Σ yᵢsᵢ°/Mᵢ that `sdep0` carried.

The commit message's "s(Tref,gPref) came out at 13.6 instead of 6.86 kJ/(kg·K)" was reasoned
from the code on the assumption that `sdep0` is a pure departure — it was not measured.
Tracked as **D17**. Related: `DISABLED_Gas_Entropy_RefState` is the natural regression gate,
but at `gPref` the cubic legitimately differs from the ideal gas by sdep(gTref,gPref) ≈
−0.5 J/(kg·K) ≈ 7e-5 relative, against its 1e-4 tolerance — re-tolerate before enabling.

## 1c — Methane FD consistency: the step, not the equation

With the sign fix in, `Methane_Cp_FD_NearCritical` no longer aborts but misses r2 = 1 by
1.68e-4. That is **not** a BELFEM defect: running the identical formula on CoolProp gives
1 − r2 = 1.682e-4, against BELFEM's 1.679e-4. The test's ±0.05 % temperature step truncates by
1.2 % at the 193.25 K sample, where cp spikes to 72 kJ/(kg·K) on the pseudo-critical line, and
that one point dominates the sum. Step tightened to ±0.01 %: truncation drops to 2.8e-7 in r2
while the enthalpy difference stays four orders above the noise the density solve leaves in h
(|Δp/p| < 1e-8 → 0.023 J/kg, against 2785 J/kg of signal). The agreement between BELFEM and
CoolProp on the *same* flawed measurement is itself the strongest confirmation that the EoS is
now right.

## 2 — Prandtl-Meyer perfect-gas limit: the premise was too wide

**Symptom.** `GASMODELS.PrandtlMeyerPerfectGasLimit` fails its own `ASSERT_NEAR` premise:
cp(200 K)/cp(4000 K) = 1 + 6.78e-6 against a 1e-12 tolerance.

**Cause.** Not a defect — argon's NASA table only carries cp = 2.5 R *exactly* on the
250–1000 K interval, where all polynomial coefficients but a₃ are zero. Below 250 K the table
hands over to a fitted low-temperature interval (84.306–250 K) whose seven-term polynomial
overshoots by 6.8e-6 at 200 K after heavy cancellation; above 1000 K the next interval carries
a₃ = 2.500069401 with small compensating powers of T (net +1.1e-8 at 4000 K, +4.9e-7 at 6000 K).

The test picked T₁ = 4000 K and Mach targets up to 10, which puts the downstream state at
~156 K — deep inside the fitted interval. No choice of T₁ makes the old sweep exact: expansion
to Ma = 10 costs a factor 25 in temperature, so the isentrope cannot stay inside a window that
spans a factor 4.

**Fix.** T₁ = 950 K with targets to Ma = 3.0, which puts the coldest state at 318 K — the whole
isentrope inside the exact window. The premise assertion now checks the coldest state against
the upstream state. (The junction at 1000 K does get a glue polynomial over [995, 1005]; 950 K
is clear of it.)

The measured ratio at 318 K/950 K is 1 − 8.2e-10, not 1 to machine precision: the caloric data
is a spline through the tabulated *enthalpy*, so cp is a spline derivative and carries ~1e-9 of
interpolation noise even where the underlying polynomial is exactly constant. Premise assertion
therefore set to 1e-8 and the two comparisons to 1e-6 — near the sonic point the inverse map
ν → T amplifies that noise, since dν/dMa vanishes at Ma = 1.

The high-Mach coverage is deliberately dropped rather than bought with a loosened tolerance:
the point of this test is that the exact method reproduces the closed form *exactly* when the
gas is exactly perfect.

## 3 — Gibbs test: `r2` overload ambiguity (compile error)

`tests/physics/gasmodels/cl_GM_Gas_Gibbs.cpp:188` called `r2( tValues.col(j), tExpect.col(j) )`.
`Matrix::col()` returns a backend view, and **both** `Vector<real>` and `Matrix<real>` have
implicit converting constructors from a column-vector expression (`cl_BZ_Vector.hpp:95`,
`cl_BZ_Matrix.hpp:97`) — so with two column views neither `r2` overload in `fn_r2.hpp` is
better. The calls above it compile because their second argument is already a `Vector`.
Fixed by binding both columns to `Vector<real>` first; `fn_r2.hpp` untouched. The test itself
stays `DISABLED_` for the reasons in its header comment (NaN sweep, stale reference data).

## Gates

Build and `make check` are Christian's. Expected: the methane near-critical FD test and the
hydrogen/methane cp-consistency neighbors pass, and the Prandtl-Meyer suite runs its exact
limit at 1e-8. Not compiled here.
