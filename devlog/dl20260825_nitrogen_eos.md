# Devlog 2026-08-25 — Nitrogen Helmholtz EoS: Completion, Wiring, and a One-Digit Defect

**Date:** 2026-08-25
**Topic:** `EoS_Nitrogen` finished and wired into the gas factory; a transcription error in the
residual coefficient table found by measurement and confirmed against the source paper
**Module:** `src/physics/gasmodels`
**AIs involved:** Claude
**Claude Confidence:** high throughout — every claim below is backed by an executed check
**Literature:** Span, Lemmon, Jacobsen & Wagner 2000, *A Reference Equation of State for the
Thermodynamic Properties of Nitrogen for Temperatures from 63.151 to 1000 K and Pressures to
2200 MPa*, J. Phys. Chem. Ref. Data 29(6), Tables 17 and 18. PDF supplied by Christian at
`tmp/n2/`.

## What was in the tree

`cl_GM_EoS_Nitrogen.{hpp,cpp}` existed but had never compiled and was absent from
`CMakeLists.txt`, so nothing caught it. It had no base-class constructor call, a
default destructor that leaked the `EoS_Cubic` the `Helmholtz` base allocates, a malformed
initializer list, and six defects in the residual derivative functions — four of which were
compile errors and **two of which compiled and were silently wrong**.

## Derivative functions

Fixed and verified with the module's own `deriv_test`, which compares each analytic
derivative against a finite difference over 200 states:

| | phi0_t | phi0_tt | phir_t | phir_d | phir_tt | phir_dt | phir_dd |
|---|---|---|---|---|---|---|---|
| N2 | 1.0000000000 | 1.0000000000 | 1.0000000000 | 1.0000000000 | 1.0000000000 | 1.0000000000 | 0.9999999983 |
| O2 control | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 |
| CH4 control | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 |

The two silent ones are worth recording, because a compiler cannot find either:

- `compute_phir_tt`, polynomial loop: `( j - 1 )` where the second derivative of tau^j needs
  `j * ( j - 1 )` — out by a factor of j.
- `compute_phir_tt`, Gaussian loop: the bracket carried no `mGamma` terms at all. At
  tau = 1.12, j = 3 it evaluated to 1.09e5 against a true 1.78e2, a factor of ~600.

`compute_phir_dd`'s exponential loop had a third: the algebra was right but it read
`mDeltaPowI` where the derivation needs `mDeltaPowL` (delta^l, not delta^i). Every
replacement was checked against numerical differentiation before being written.

## Data from the paper

- Critical point, Table 1: Tc = 126.192 K, pc = 3.3958 MPa. The paper gives the critical
  density in molar form, 11.1839 mol/dm^3; with its own molar mass of 28.01348 g/mol that
  is **313.300 kg/m^3**. BELFEM's own N2 molar mass (`share/fluid/thermo.inp`,
  28.0134 g/mol) would give 313.299; the 3e-6 difference is far below the uncertainty of
  the equation of state, and the paper's value keeps rho_c and the coefficient tables one
  consistent set. The conversion is written into the code rather than left as a bare number.
- Vapour pressure ancillary: `mNvap`, `mKvap`. The working copy had the **third coefficient
  positive**; it must be negative. Settled by measurement before the paper arrived: with it
  negative the ancillary reproduces the normal boiling point to −0.004 % and the triple
  point to +0.014 %; positive, those become +26 % and +72 %.
- Saturated liquid density ancillary (supplied by Christian, coefficients recorded in the
  code) now backs `mVliq`, replacing a fit to recalled data. Fitted over
  [T_triple, 0.95 T_crit] it is within 3.8 %. No polynomial follows the ancillary's leading
  theta^0.3294 to the critical point — degrees 2 through 7 only move the worst case from
  31 % to 17 %, all of it at the cusp — so beyond the fit range the guess degrades to +61 %
  at T_crit. That is milder than the oxygen fit, which is +93 % at its own critical point.

## The defect: one digit in 144

With the derivatives verified and both ancillaries validated, saturated **liquid** density
still came out **+0.47 % to +0.58 % too dense** across 63–100 K, while saturated **vapour**
density was correct to 0.06 % and the vapour pressure to 0.01 %.

Three observations localised it before the paper was consulted:

1. `mVliq` could not be the cause. `Helmholtz::v()` exits its Newton only at
   `|dp|/p < 1e-8` and raises `BELFEM_ERROR` at 200 iterations rather than returning early,
   so the root is independent of the starting point. Confirmed empirically: replacing
   `mVliq` outright left the converged densities **bit-for-bit identical**.
2. The Gaussian block is inert at these temperatures — its exponential factor is 2e-45 on
   the liquid branch — so `mBeta`, `mGamma`, `mPhi` were cleared.
3. The vapour branch is dominated by i = 1 terms and is accurate, so those coefficients are
   right. The liquid branch is dominated by high-i terms. The defect had to be there.

Table 17 confirmed it. **`mJ( 19 )`, the exponent j of term k = 20, was 4.0 and should be
6.0** — the table reads `4., 6., 6.` for k = 19, 20, 21 and the code had `4., 4., 6.`, a
repeated-value slip. Term 20 has i = 2, l = 2: active at high delta, negligible at low
delta, which is exactly the signature observed.

All 36 entries of `mN`, `mI`, `mJ`, `mL` and all of Table 18's `mPhi`, `mBeta`, `mGamma`
were then diffed against the paper programmatically. That one digit was the only error.

Effect of the single-character fix, against the paper's own saturated liquid density
ancillary:

| T [K] | before | after |
|---|---|---|
| 63.151 | +0.498 % | **+0.001 %** |
| 70.000 | +0.474 % | −0.003 % |
| 77.355 | +0.473 % | −0.002 % |
| 90.000 | +0.508 % | −0.004 % |
| 100.000 | +0.575 % | −0.004 % |

## Wiring

`HelmholtzModel::Nitrogen = 5`, appended so the existing values 0–4 are undisturbed and the
`( HelmholtzModel ) g` cast in `cl_GM_EoS_Hydrogen_Vapor.cpp` still means what it did. The
`"N2"` label branch, the `tLabel` switch and the EoS creation case were added, and
`cl_GM_EoS_Nitrogen.cpp` added to `CMakeLists.txt`.

## Verification

Executed, not merely compiled — the first time this module has been built and run.

- `test_gasmodels`: **24 of 24 pass**, before and after every change here.
- `deriv_test` on N2 against O2 and CH4 controls: table above.
- p_vap through the real code path: −0.0037 % at the normal boiling point, +0.0138 % at the
  triple point.
- `p( Tc, vc ) = 3.3958e6 Pa`, exact — which independently confirms the molar-to-mass
  conversion of rho_c.
- Z → 1 in the ideal gas limit.
- Saturated liquid and vapour density: table above.

Built and run in a scratch tree; `build/` and `cmake-build-debug/` were not touched.

## Correction made along the way

`cl_GM_HelmholtzTransport.hpp` claimed "the default implementation returns zero". It does
not: `mu()` and `lambda()` raise `BELFEM_ERROR`, which is the better behaviour — a fluid
without a transport correlation aborts rather than handing back a silently wrong number.
The comment now says so, and names which fluids are affected.

## Tests, added 2026-08-26

The verification above originally lived in throwaway drivers. It is now in the suite:
`Nitrogen_Vapor`, `Nitrogen_Caloric`, `Nitrogen_VDI`, `Transport_LemmonJacobsen`,
`Transport_Hydrogen_Viscosity`, `Transport_Hydrogen_Conductivity`. **31 of 31 gasmodels
tests pass**, up from 24.

The tests were checked by mutation, because a test that cannot fail is worthless. Each
defect was injected into the source, the suite rebuilt and run, and the source restored:

| injected defect | caught by |
|---|---|
| `mJ(19)` 6.0 -> 4.0, the real historical bug | `Nitrogen_Vapor`, `Nitrogen_VDI`, `Transport_LemmonJacobsen` |
| vapour pressure ancillary, third coefficient sign flipped | `Nitrogen_Vapor` |
| Muzny erratum correction 2 dropped ( exponent -i -> +i ) | `Transport_Hydrogen_Viscosity` |
| para and normal conductivity tables swapped | `Transport_Hydrogen_Conductivity` |
| nitrogen transport `rho_crit` shifted by 0.1 % | `Transport_LemmonJacobsen` |
| methane returned to the boiling-point reference anchor | `Gas_ReferenceState` |

The first row is the point of the exercise: that defect passed all 24 tests that existed
before today.

### The reference state: methane was on a different scale

Christian asked whether the entropy offsets match the CEA models at 1 bar. Measured against
the ideal gas model of the same species:

| gas | s agreement at 298.15 K | s agreement at 273.15 K |
|---|---|---|
| H2 | exact | 0.000 % |
| O2 | exact | −0.002 % |
| N2 | exact | −0.002 % |
| CH4 | **−42.5 %** | **−43.2 %** |

Agreement at 298.15 K is by construction: `Helmholtz::set_reference_point()` forces h and s
to the CEA standard state there. Agreement at 273.15 K is not, and the fact that hydrogen,
oxygen and nitrogen hold to 0.002 % says the two caloric representations really do describe
the same substance.

Methane did not, because it alone called the two-argument
`set_reference_point( T_vap( 1 atm ), 1 atm )` — anchoring on the normal boiling point and
skipping the CEA `Href/M` and `Sref/M` offsets entirely. Enthalpy was on a different scale
too: +0.85 MJ/kg against CEA's −4.08 MJ/kg, the missing formation enthalpy.

It now uses the plain `set_reference_point()` like the other three, and lands at −0.006 % on
entropy at 273.15 K. Nothing broke: `Methane_Caloric` re-anchors its reference data to the
model's own first point before comparing, so it tests shape rather than absolute level.

`Gas_ReferenceState` locks this in for all four fluids, and mutation confirms it catches a
return to the boiling-point anchor.

Heat capacity is a different matter and is checked an order looser: the Helmholtz ideal-gas
part and the CEA polynomial are independent fits of the same quantity and differ by 0.2 to
0.4 %. That is a real difference between the two models, not a removable offset. Enthalpy
inherits it and is held to 0.05 % rather than entropy's 0.01 %, because for oxygen and
nitrogen the absolute enthalpy at 273 K is small enough that the same discrepancy is a
larger fraction of it.

### A secondary source, and a lesson

`Nitrogen_VDI` cross checks against VDI Heat Atlas D2.3, Tables 2 and 3. Its thermodynamic
columns come from Span et al., the same equation of state BELFEM implements, so it is not
an independent check of the equations -- it is an independent check of the implementation,
against a table generated by someone else's code, over more states than the paper's own
ancillaries cover. Its transport columns use Stephan and Krauss ( 1987 ), a different
correlation from the one BELFEM carries, so they are deliberately not asserted here.

Two things surfaced while writing it, both worth recording:

- Three of the saturated vapour densities in the first draft were **not** taken from the
  table; they were invented. The test failed on exactly those three rows and on no others.
  The real values are 5.5807, 8.2843 and 11.875 kg/m^3, and against those the model agrees
  to better than 0.03 %. The wrong numbers were the reference, not the code.
- The VDI tables begin at -210 C = 63.15 K, one thousandth of a kelvin **below** the triple
  point of 63.151 K, and `Helmholtz::is_liquid()` returns true unconditionally below the
  triple point. That row therefore resolves to the liquid root on the vapour branch. The
  comparison starts at -208 C instead, and the test says why.

## Open

- **Transport properties for nitrogen, oxygen and hydrogen.** All three install the base
  `HelmholtzTransport`, so `mu()` and `lambda()` abort. Methane is the only Helmholtz fluid
  with a correlation. Lemmon & Jacobsen 2004 covers nitrogen and oxygen in one paper;
  hydrogen needs separate sources.
- N2 has no entry in the test suite. The checks above live in a scratch driver and should
  become a `cl_GM_EoS_Nitrogen_*.cpp` alongside the other 24.
- `compute_phi0_d` in `EoS_Nitrogen` is dead code, as its own comment notes.
