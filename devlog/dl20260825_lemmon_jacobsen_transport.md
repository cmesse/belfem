# Devlog 2026-08-25 — Transport Properties for Nitrogen and Oxygen (Lemmon & Jacobsen)

**Date:** 2026-08-25
**Topic:** Viscosity and thermal conductivity for the nitrogen and oxygen Helmholtz fluids,
one class for both; verified against the source paper's own verification table
**Module:** `src/physics/gasmodels`
**AIs involved:** Claude
**Claude Confidence:** high — reproduces Table V of the paper to the printed digits
**Literature:** Lemmon & Jacobsen, *Viscosity and Thermal Conductivity Equations for
Nitrogen, Oxygen, Argon, and Air*, Int. J. Thermophys. 25:21–69 (2004). Critical
enhancement after Olchowy & Sengers, repeated in that paper as Eqs. (7)–(11). PDF supplied
by Christian at `tmp/lemmon.pdf`.

## Why this was next

`mu()` and `lambda()` on the base `HelmholtzTransport` raise `BELFEM_ERROR`. Before this
session methane was the only Helmholtz fluid with a real correlation; hydrogen, oxygen and
the newly finished nitrogen all installed the base class, so asking any of them for a
transport property aborted the run.

## One class, two fluids

The starting suggestion was to copy `HelmholtzTransport_Methane`. The **architecture**
transfers — dilute gas plus residual plus a critical enhancement, served from a state cache
keyed on the last (T, p) — but the **equations do not**. Methane follows Friend et al.
(1989); this is a different correlation with a different collision integral, a different
residual form and a different critical enhancement. Nothing was carried across
symbolically, and the header says so, so that nobody later assumes the two are variants of
one model.

What the paper does give is better than a per-species class: **one functional form serves
nitrogen, argon, oxygen and air**, with only the coefficient tables differing. So this is a
single class carrying the fluid as a constructor parameter, the way `EoS_Hydrogen` carries
its three spin isomers, rather than two near-identical classes:

```cpp
new HelmholtzTransport_LemmonJacobsen( *this, HelmholtzModel::Nitrogen )
new HelmholtzTransport_LemmonJacobsen( *this, HelmholtzModel::Oxygen )
```

Nitrogen and oxygen are implemented, those being the two fluids BELFEM has a matching
Helmholtz EoS for. Argon and air are the same code plus a coefficient table if an EoS for
them ever lands.

Files: `cl_GM_HelmholtzTransport_LemmonJacobsen.{hpp,cpp}`, added to `CMakeLists.txt` and
wired into both branches of the EoS creation switch in `cl_Gas.cpp`.

## Units, which is where this kind of port usually goes wrong

The paper works in micro Pa s and milli W m^-1 K^-1; the accessors return SI. Three details
were easy to get wrong and are commented in place:

- `eta_0` keeps M in g/mol and sigma in nm, the units Eq. (2) is written for, and yields
  micro Pa s.
- Eq. (5) takes the dilute gas viscosity **divided by one micro Pa s**, so the unscaled
  `eta_0` feeds the conductivity, not the SI value.
- The critical enhancement, Eq. (7), is dimensionally SI on its own — density in kg/m^3,
  cp in J/(kg K), viscosity in Pa s, correlation length in m give W/(m K) directly. Only
  `lambda_0` and `lambda_r` carry the 1e-3; scaling `lambda_c` as well would be wrong.

`constant::kB` is the exact SI-2019 value where the paper quotes 1.380658e-23. The 6.5e-6
relative difference sits far below the digits of the verification table.

## Verification against Table V

The paper publishes calculated values expressly for program verification. Both fluids
reproduce them to the printed digits.

**Nitrogen**

| T [K] | rho [mol/dm^3] | eta paper | eta BELFEM | err | lambda paper | lambda BELFEM | err |
|---|---|---|---|---|---|---|---|
| 100 | 0 | 6.90349 | 6.9035 | 0.000 % | 9.27749 | 9.2775 | 0.000 % |
| 300 | 0 | 17.8771 | 17.8771 | 0.000 % | 25.9361 | 25.9361 | 0.000 % |
| 100 | 25 | 79.7418 | 79.7418 | 0.000 % | 103.834 | 103.8342 | 0.000 % |
| 200 | 10 | 21.0810 | 21.0810 | 0.000 % | 36.0099 | 36.0099 | 0.000 % |
| 300 | 5 | 20.7430 | 20.7430 | 0.000 % | 32.7694 | 32.7694 | 0.000 % |
| 126.195 | 11.18 | 18.2978 | — | — | 675.800 | — | — |

**Oxygen**

| T [K] | rho [mol/dm^3] | eta paper | eta BELFEM | err | lambda paper | lambda BELFEM | err |
|---|---|---|---|---|---|---|---|
| 100 | 0 | 7.70243 | 7.7024 | 0.000 % | 8.94334 | 8.9433 | 0.000 % |
| 300 | 0 | 20.6307 | 20.6307 | 0.000 % | 26.4403 | 26.4403 | 0.000 % |
| 100 | 35 | 172.136 | 172.1358 | −0.000 % | 146.044 | 146.0437 | −0.000 % |
| 200 | 10 | 22.4445 | 22.4445 | 0.000 % | 34.6124 | 34.6124 | −0.000 % |
| 300 | 5 | 23.7577 | 23.7577 | 0.000 % | 32.5491 | 32.5491 | 0.000 % |
| 154.6 | 13.6 | 24.7898 | 24.7890 | −0.003 % | 377.476 | 378.0997 | **+0.165 %** |

The oxygen critical-point row is the only one that is not exact. It is also the single most
sensitive state in the table: `lambda_c` contributes some 90 % of that number, and it is
driven by `(d rho / d p)_T` from the equation of state right where that derivative
diverges. A further contributor is that the transport class carries the paper's Table I
critical constants while the derivatives come from `EoS_Oxygen`, whose critical point is
Schmidt & Wagner's; the two are close but not identical. 0.165 % at that state is well
inside the correlation's own stated uncertainty.

`test_gasmodels`: **24 of 24 pass**, unchanged.

## Found on the way: v(T,p) cannot invert at the critical point

The nitrogen row of Table V could not be evaluated at all. `Helmholtz::v()` inverts
p(T,v) = p by Newton, and at the critical point the isotherm is flat, so the step diverges
and the 200-iteration guard fires:

```
Too many iterations for T= 126.195 K, p= 33.963 bar, rho= 309.748 kg/m^3, relax= 0.264
```

Measured reach along the critical isochore of nitrogen:

| T / T_crit | 1.0001 | 1.0010 | 1.0050 |
|---|---|---|---|
| v( T, p ) | fails | fails | converges, lambda = 89.15 mW/(m K) |

So the interface is usable from roughly half a percent above T_crit outward and unusable
within about a tenth of a percent. This is **pre-existing** and affects every Helmholtz
fluid, not just the new transport model — but it bites hardest here, because the critical
enhancement is precisely the term that matters in the region the inversion cannot reach.
A bracketed fallback for the near-critical isotherm would fix it; that changes
`Helmholtz::v()` for all fluids and is left as a decision rather than slipped in.

## Corrected along the way

`cl_GM_HelmholtzTransport.hpp` claimed the base implementation "returns zero". It raises
`BELFEM_ERROR` instead, which is the better behaviour — a fluid without a correlation
aborts rather than handing back a silently wrong number. The comment now says so and names
which fluids are still on the placeholder.

## Hydrogen, added the same session

`HelmholtzTransport_Hydrogen` implements the **viscosity** of normal hydrogen after Muzny,
Huber & Kazakov, J. Chem. Eng. Data 58:969-979 (2013), **with the 2022 erratum applied**
(J. Chem. Eng. Data 67:2855). The erratum is not optional — it changes three things, and
without it the correlation is simply wrong:

1. Eq. (6) is missing Avogadro's number.
2. In Eq. (7) the exponent on T* is **-i**, not +i.
3. The density scale of Eq. (9) is rho_sc = 90.909090909 kg/m^3; the body text says 90.5.

Two further traps, neither of them stated in the erratum:

- **Table 3 does not extract cleanly.** Its two columns interleave, so indices 3 and 5 land
  on the rows of 2 and 4. The reading b = { -0.1870, 2.4871, 3.71513, -11.0972, 9.09655,
  -3.8292, 0.5166 } was confirmed numerically, not by eye — see below.
- The erratum says N_A belongs "in the denominator", but the combination that actually
  reproduces the test values is `B* sigma^3 N_A / M`, giving the second viscosity virial in
  m^3/kg so that it multiplies a density in kg/m^3 directly. This was determined by
  measurement rather than from the printed wording.

The arbitration was the erratum's own three test values. Solving for the second viscosity
virial from the rho = 50 and rho = 100 points **independently** gave 8.928603e-3 and
8.928698e-3 m^3/kg — agreeing to 1e-5 relative, which simultaneously confirms the Table 3
reading, the symbolic-regression term and rho_sc. `B* sigma^3 N_A / M` evaluates to
8.929172e-3, matching to 6e-5.

Verified through the real code path:

| T [K] | rho [kg/m^3] | eta erratum | eta BELFEM | err |
|---|---|---|---|---|
| 40 | 0 | 1.9772 | 1.9772 | +0.0006 % |
| 40 | 50 | 5.9905 | 5.9906 | +0.0009 % |
| 40 | 100 | 49.034 | 49.0341 | +0.0002 % |

The correlation is for **normal** hydrogen. BELFEM's EoS distinguishes the para, normal and
ortho isomers but no transport correlation for the spin isomers exists, so this one is
installed for all three, as REFPROP does. The density entering the correlation still comes
from whichever EoS the gas carries, so the isomers do differ in the result.

`constant::NA` is the exact SI-2019 value where the erratum quotes 6.022137e23; the 1.2e-7
relative difference is far below the digits of the test values.

## Hydrogen thermal conductivity, also the same session

`HelmholtzTransport_Hydrogen::lambda()` implements Assael, Assael, Huber, Perkins & Takata,
J. Phys. Chem. Ref. Data 40:033101 (2011). Unlike the viscosity paper, this one **does**
distinguish the isomers and publishes separate coefficient tables for normal and
parahydrogen, so `lambda()` is isomer aware where `mu()` is not. Orthohydrogen has no table
of its own and is served by the normal one.

Verified against the paper's Table 7:

| T [K] | rho [kg/m^3] | normal paper | normal BELFEM | err | para paper | para BELFEM | err |
|---|---|---|---|---|---|---|---|
| 298.15 | 0 | 185.67 | 185.671 | +0.001 % | 192.38 | 192.379 | −0.000 % |
| 298.15 | 0.80844 | 186.97 | 186.972 | +0.001 % | 192.81 | 192.803 | −0.003 % |
| 298.15 | 14.4813 | 201.35 | 201.353 | +0.002 % | 207.85 | 207.848 | −0.001 % |
| 18 | 0 | 13.875 | 13.875 | +0.002 % | 13.643 | 13.643 | −0.003 % |
| 18 | 75 | 104.48 | 104.476 | −0.004 % | 100.52 | 100.510 | −0.010 % |
| 35 | 30 | 75.594 | 77.014 | **+1.9 %** | 70.335 | 74.460 | **+5.9 %** |

Five of six rows to 0.01 %. The sixth is the only state where the critical enhancement
matters, and the discrepancy has a specific, documented cause rather than being a
transcription error: **the enhancement is inversely proportional to the background
viscosity**, and Assael et al. fitted it in 2011 against the viscosity correlation REFPROP
carried then — which predates Muzny et al. 2013. BELFEM uses the newer viscosity.

Decomposing that row confirms it. The dilute-gas plus excess base is exact, as the other
five rows show; the whole discrepancy sits in the enhancement:

| | base | enhancement, paper | enhancement, BELFEM | error in the enhancement |
|---|---|---|---|---|
| normal | 60.72 | 14.88 (20 % of total) | 16.30 | +9.5 % |
| para | 59.22 | 11.12 (16 % of total) | 15.24 | +37.1 % |

Parahydrogen is far worse for a compounding reason: its background viscosity is the
*normal* hydrogen correlation, so the enhancement inherits that error too. This is the
concrete cost of the isomer approximation, and it is now written into the module README
rather than left as a code comment.

## Open

- The transport correlations have no entry in the test suite. All three papers publish
  verification tables — Lemmon & Jacobsen Table V, the Muzny erratum Table 1, Assael Table 7
  — which together make an exact, cheap, self-contained oracle covering all four fluids and
  both properties. They belong in gtest rather than in scratch drivers.
- Neither the transport model nor the nitrogen EoS has an entry in the test suite. Table V
  is an exact, cheap, self-contained oracle and should become a gtest case rather than
  living in a scratch driver.
- Argon and air are the same code plus one coefficient table each, if an equation of state
  for them is ever added.
