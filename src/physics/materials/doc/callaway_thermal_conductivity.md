# Callaway Phonon Thermal Conductivity — Usage Guide {#physics_materials_callaway_thermal_conductivity}

**Date:** 2026-08-25
**Purpose:** How to use `callaway_conductivity`, the Fortran kernel that evaluates the phonon
thermal conductivity λ_ph(T) of the Callaway model: the 20-parameter interface, the scattering
channels it wires, and how YBCO drives it
**Module:** `src/physics/materials`

The parameter indices, conversions, and mechanism wiring below come from `debye.f90` (the
kernel) and the two callers in `cl_Material_YBCO.cpp`. Where the header prose in
`debye.hpp` and the kernel disagreed in the past, the kernel is the ground truth.

> **Fixed 2026-08-25, refit degenerate — read §8 before touching YBCO's parameters.** The optical
> phonon–electron channel had a unit defect that switched it off; the kernel now evaluates it
> correctly. Refitting YBCO's parameters against the Sommerfeld et al. 2003 data showed that the
> class's *electronic* term reproduces the whole curve on its own (10 % rms), so the phonon
> parameters are not identifiable from that data; the values in the class are an interim set that
> restores the previous agreement, not measurements.

---

## 1. Purpose and location

`callaway_conductivity` computes the **phonon (lattice) contribution** λ_ph(T) to thermal
conductivity with Callaway's single-relaxation-time model. It is a Fortran routine with C linkage:

```c
void callaway_conductivity( const double * params,   // [20] inputs, §3
                                  double * result,   // λ_ph [W/(m·K)]
                                  int    * status ); // 0 = ok, 1 = error
```

Declared in `debye.hpp`; implemented in `debye.f90:207`.

**Callers.** The production entry point is the `YBCO::lambda_custom( const real T )`
override (`cl_Material_YBCO.cpp:475`). It packs the parameter array, calls the kernel, and returns
λ = λ_e + λ_ph — the Wiedemann–Franz electronic term from the normal-state resistivity plus the
Callaway phonon term. `YBCO::test_callaway( T,
aParams )` (`cl_Material_YBCO.cpp:406`) is the fitting harness; it packs the
same array, with the coupling taken from its argument. No other material uses the kernel:
`Metal::lambda_custom` is the Hust form (`fn_hust.hpp`), and there is no `lambda_phonon` method
anywhere in the tree — a `debye.hpp` comment once named one, and that comment was wrong.

The composition-driven inputs, molar mass `M` and impurity parameter `Γ`, are computed by
`material::Abundance` during construction (§6).

---

## 2. The physics in brief: four rates and one modifier

The kernel integrates

```
λ_ph = ( kB / (2π² v_g) ) · ( kB·T / ℏ )³ · ∫₀^(θ/T) τ(z) · z^n · e^z / (e^z − 1)² dz
```

(final assembly `debye.f90:428`; integrand `scatter_function`,
`debye.f90:171`), with z = ℏω/(kB·T), ω = K(1)·z, K(1) = kB·T/ℏ
(`debye.f90:303`), and one relaxation time τ from four additive rates
(Matthiessen), plus a superconducting reduction that multiplies the phonon–electron rate:

| # | Channel | Rate | Set from params (C index) | Where |
|---|---|---|---|---|
| 1 | Umklapp (phonon–phonon) | τ_U⁻¹ = K(2)·ω², with exp(−θ/(b·(1 + d·T/θ)·T)) when b ≠ 0 | θ(1), G(4), γ(5), b(11), d(12); V₀ from M(6), ρ(3) | `debye.f90:309` |
| 2 | Mass-difference impurity | τ_M⁻¹ = K(3)·ω⁴ | Γ(7), M(6), ρ(3), v_g(2) | `debye.f90:314` |
| 3 | Boundary | τ_B⁻¹ = K(4) = v_g / L₀ | v_g(2), L₀(9) | `debye.f90:317` |
| 4 | Phonon–electron | τ_phe⁻¹ = K(5)·ω | n_e(15), ε(13), m_eff(14), λ_pe(17), λ_opt(18), ω_opt(19) | `debye.f90:319` ff. |
| 5 | Superconducting reduction | multiplies K(5) | T_c(8), Δ₀ ratio(16) | `debye.f90:343` ff. |

The rates are summed in `scatter_function` (`debye.f90:181`–`debye.f90:190`).

---

## 3. The 20-parameter interface

Fortran `params(k)` is C `params[k-1]`. The "YBCO" column gives the value `YBCO::lambda_custom` passes.

| C | Fortran | Symbol | Meaning | Unit as passed | Converted to | Feeds | YBCO |
|--:|--:|---|---|---|---|---|---|
| 0 | 1 | T | temperature | K | — | everything | runtime T |
| 1 | 2 | θ | Debye temperature | K | ω_D = θ·kB/ℏ (`debye.f90:307`) | Umklapp; integral limit θ/T | `debye(T)` |
| 2 | 3 | v_g | phonon group velocity | m/s | — | τ_M, τ_B, τ_phe, prefactor | `group_velocity(T)` |
| 3 | 4 | ρ | density | kg/m³ | — | V₀ = M/(ρ·N_A); τ_phe | `density(T)` |
| 4 | 5 | G | shear modulus | Pa | — | Umklapp | `G(T)` |
| 5 | 6 | γ | Grüneisen parameter | – | — | Umklapp | `grueneisen(T)` |
| 6 | 7 | M | molar mass | kg/mol | — | V₀ | `Abundance` (§6) |
| 7 | 8 | Γ | mass-difference impurity parameter | – | — | τ_M | 1.0 effective (§6, §8) |
| 8 | 9 | T_c | critical temperature | K | — | SC branch gate; Δ₀ | `T_crit` (92.5) |
| 9 | 10 | L₀ | boundary length (REBCO layer thickness) | m | — | τ_B | `layer_thickness` |
| 10 | 11 | n | integrand exponent | – | — | `scatter_function` | 4 |
| 11 | 12 | b | Umklapp correction strength | – | — | Umklapp exponential (skipped if 0) | 26.1 |
| 12 | 13 | d | Umklapp temperature correction | – | — | Umklapp exponential | 5 |
| 13 | 14 | ε | deformation potential | **eV** | × e → J (`debye.f90:274`) | τ_phe (acoustic / fallback) | 10 |
| 14 | 15 | m_eff | effective mass | **× m_e** | × m_e (`debye.f90:275`) | τ_phe fallback | 5 |
| 15 | 16 | n_e | conduction-electron density | 1/m³ | — | τ_phe | 2e27 |
| 16 | 17 | Δ₀/(kB·T_c) | d-wave gap **ratio** | – | × kB·T_c → Δ₀ [J] (`debye.f90:278`) | SC branch | 2.1 |
| 17 | 18 | λ_pe | acoustic e–ph coupling | – | — | selects acoustic model | 0 |
| 18 | 19 | λ_opt | optical e–ph coupling | – | — | selects optical model | 1.06 |
| 19 | 20 | ν̃_opt | optical (Raman) band | **1/cm** | 100·ν̃·h·c → energy [J] (`debye.f90:288`) | τ_phe optical | 501 |

In practice, `v_g` is not a free parameter. `Metal::group_velocity(T)`
(`cl_Material_Metal.cpp:669`) computes it from the
Debye temperature and the atomic volume, v_g = (4πM / (3 q N_A ρ))^{1/3} · kB θ / h — about
2.6 km/s for copper at θ = 343 K.

---

## 4. Phonon–electron sub-models

K(5) is set by exactly one of three paths:

**Optical model** — when λ_opt > 0 (`debye.f90:321`):
```
U = ℏ·ω_opt/(kB·T);   V = exp(U)
A = n_e·kB·T/(ρ·v_g²) · U²·V/(V − 1)
K(5) = A·λ_opt ;  compute_density_function = .true.
```
Intended for YBCO's oxygen Raman band. **This is the branch with the unit defect of §8.**

**Acoustic model** — when λ_pe > 0 (`debye.f90:333`):
```
A = n_e·ε²/(ρ·v_g²·kB·T)
K(5) = A·λ_pe ;  compute_density_function = .true.
```

**Mutual exclusion:** if both couplings are positive, the kernel sets `status = 1`,
`result = 0` and returns (`debye.f90:296`). Always check
`status`.

**Fallback** — neither coupling is set (`compute_density_function` false): the kernel uses a default
deformation-potential rate from n_e, ε, m_eff, guarded by all three being positive
(`debye.f90:393`), else K(5) = 0:
```
A = n_e·ε²/(ρ·v_g²·kB·T);   U = m_eff·v_g²/(2·kB·T)
K(5) = A·√(π U)/exp(U)
```
with a simple BCS-like (not d-wave) reduction 2/exp(Δ₀/(kB·T)) if Δ₀ ≠ 0
(`debye.f90:400`).

**Choosing:** YBCO uses the optical model. A material whose dominant coupling is acoustic sets
λ_pe > 0, λ_opt = 0. The fallback needs n_e, ε and m_eff only.

---

## 5. Superconducting branch (d-wave)

This branch runs only when T < T_c, Δ₀ > 0, and an e–ph model is active
(`debye.f90:343`); otherwise K(5) keeps its normal-state value (the fallback path
has its own reduction, §4).

It multiplies K(5) by Y/(kB·T) (`debye.f90:389`), where Y integrates a
density-of-states-weighted Fermi window N(E)·f·(1 − f) over the reduced energy E_red = E/Δ₀:

- subgap, E_red < 1: N_ratio = (2/π)·E_red — the linear nodal density of states of a d-wave gap
  (`debye.f90:371`);
- above the gap: N_ratio = E_red / √(E_red² − 1), the BCS edge; the integration step is chosen so
  the square root is never evaluated at the singularity
  (`debye.f90:379`).

Below T_c the electrons condense and phonon–electron scattering weakens; the nodal
quasiparticles make the suppression partial and temperature-dependent rather than an s-wave
exponential freeze-out — the phonon-conductivity rise below T_c seen in cuprates.

---

## 6. Composition: where M and Γ come from

`M` and `Γ` are computed from the stoichiometry by `material::Abundance` in the YBCO constructor — although YBCO then overrides Γ with the effective disorder value of §8 —
(`cl_Material_YBCO.cpp:65`):

```cpp
Abundance tAbundance ;
real M, Gamma ;
tAbundance.compute_molar_mass_and_impurity_from_volumes(
        { "Y", "Ba", "Cu", "O" }, { 1.0, 2.0, 3.0, 7.0 }, M, Gamma );
this->set_constant( MaterialProperty::M,     M );
this->set_constant( MaterialProperty::Gamma, Gamma );
```

Γ = Σ_i Σ_j X_i·A_j^i·[(M_j^i − M̄_i)/M̄]² is the isotope and alloy mass-variance parameter
(`cl_Material_Abundance.hpp`, formulation after Zou & Balandin 2001; isotope data from the CRC
Handbook). `compute_molar_mass_and_impurity_from_masses(...)` takes weight fractions, and
`compute_molar_mass( element )` handles a single element. At call time both constants go into the array
as C indices 6 and 7.

---

## 7. Worked examples

### 7.1 YBCO — production path

From `YBCO::lambda_custom` (`cl_Material_YBCO.cpp:449` ff.):

```cpp
params[ 0] = T ;
params[ 1] = this->debye( T );                                    // θ  [K]
params[ 2] = this->group_velocity( T );                           // v_g [m/s]
params[ 3] = this->density( T );                                  // ρ  [kg/m³]
params[ 4] = this->G( T );                                        // G  [Pa]
params[ 5] = this->grueneisen( T );                               // γ
params[ 6] = this->constant_property( MaterialProperty::M );      // Abundance
params[ 7] = this->constant_property( MaterialProperty::Gamma );  // effective 1.0, see §8
params[ 8] = this->constant_property( MaterialProperty::T_crit ); // 92.5 K
params[ 9] = this->constant_property( MaterialProperty::layer_thickness );
params[10] = 4.0 ;     // n
params[11] = 26.14 ;   // b   [fitted, interim]
params[12] = 5.0 ;     // d   [fitted, interim]
params[13] = 10.0 ;    // ε   [eV]
params[14] = 5.0 ;     // m_eff / m_e
params[15] = 2e27 ;    // n_e [1/m³]
params[16] = 2.1 ;     // Δ₀/(kB·T_c) [fitted]
params[17] = 0.0 ;     // λ_pe  — acoustic off
params[18] = 1.06 ;    // λ_opt — optical on [fitted]
params[19] = 501.0 ;   // oxygen Raman band [1/cm]

callaway_conductivity( params, &k_ph, &status );
real k_e = constant::L0 * T / ( rho_i + rho_0 );   // Wiedemann–Franz, normal state
return k_e + k_ph ;
```

The optical model is active, the acoustic model is off, and the d-wave branch is active below T_c. Sources: the fit follows
Sommerfeld et al. 2003 (10.1103/PhysRevB.67.174520); m_eff and n_e from
10.1103/PhysRevLett.62.2317; the Raman band from 10.1103/PhysRevB.80.064505.

### 7.2 A normal metal — illustrative

No in-tree material drives the kernel for a normal metal; the metals use the Hust form. To
exercise the acoustic, non-superconducting path with copper-like numbers:

```cpp
double params[20] = { 0.0 };
params[ 0] = T ;
params[ 1] = 343.0 ;    // θ  [K], copper
params[ 2] = 2600.0 ;   // v_g [m/s], what group_velocity(T) gives for copper
params[ 3] = 8960.0 ;   // ρ  [kg/m³]
params[ 4] = 48e9 ;     // G  [Pa]
params[ 5] = 2.0 ;      // γ
params[ 6] = 0.06355 ;  // M  [kg/mol], or Abundance::compute_molar_mass( "Cu" )
params[ 7] = 0.0 ;      // Γ from Abundance for the real composition
params[ 8] = 0.0 ;      // T_c = 0: SC branch never fires
params[ 9] = 1e-3 ;     // L₀ [m]
params[10] = 4.0 ;
params[11] = 2.0 ;      // b
params[12] = 0.0 ;      // d
params[13] = 5.0 ;      // ε [eV]
params[14] = 1.0 ;      // m_eff / m_e
params[15] = 8.5e28 ;   // n_e [1/m³], copper
params[16] = 0.0 ;      // no gap
params[17] = 0.5 ;      // λ_pe: acoustic on
params[18] = 0.0 ;      // λ_opt: optical off
params[19] = 0.0 ;
```

The values of ε, b and λ_pe are illustrative, not fitted to a copper dataset.

---

## 8. Gotchas

1. **Acoustic and optical are mutually exclusive** (§4). Check `status`.

2. **`d` and ν̃_opt are adjacent but distinct slots** — C 12 / Fortran 13, and C 19 / Fortran 20.
   A 2026-07 fix separated them; the header, kernel and callers agree today.

3. **Silent unit conversions on input:** ε in eV, m_eff as a multiple of m_e, Δ₀ as a ratio to
   kB·T_c, the optical band in 1/cm (§3).

4. **The optical channel's unit chain, and what its repair revealed.** Until 2026-08-25
   `omega_opt = 100·params(20)·h·c` (`debye.f90:288`) — already an *energy*, ħω ≈ 1.0×10⁻²⁰ J for
   501 cm⁻¹ — was fed into `U = hbar * omega_opt / (kB*T)`, a second ħ, so U ≈ 8×10⁻³²/T and the
   channel evaluated to nothing. The kernel now forms `U = omega_opt / (kB*T)` (≈ 9.4 at 77 K).

   Refitting YBCO afterwards (Sommerfeld et al. 2003, 4–131 K; fit harness in the 2026-08-25
   devlog) found the problem is not the phonon parameters: the electronic term
   λ_e = L₀T/(ρ_i(T) + ρ₀) that `YBCO::lambda_custom` adds — a normal-state Wiedemann–Franz
   estimate carried below T_c — reproduces the measured curve by itself to 10 % rms, peak
   included. Every fit therefore converges to that same floor by suppressing λ_ph to ≲ 1 W/(m·K),
   with parameters at their bounds. Two pre-existing choices cause this: YBCO's resistivity anchor
   (60 µΩ·cm at 273 K, marked "[fitted]") is roughly four times below literature, which inflates
   λ_e; and the old Umklapp strength b = 8.17 had been calibrated against a Grüneisen input of
   ≈ 800 at 2 K (the linear-α artifact, since fixed), so with a sane γ ≈ 2.1 the phonon channel is
   wide open. The interim set in the class — b = 26.1, d = 5, Δ₀/(k_BT_c) = 2.1, λ_opt = 1.06, and
   an *effective* Γ = 1.0 replacing the isotopic 1.1×10⁻⁵ — restores the previous 10 % agreement
   with physically sane gap and coupling values. Making the phonon term identifiable requires a
   physical electronic term first (a literature ρ_ab(T) and a superconducting-state λ_e); that
   is a modeling decision, not a fit.

5. **`result` is λ_ph only.** Callers add the electronic term themselves.

---

## 9. References

- J. Callaway, *Model for Lattice Thermal Conductivity at Low Temperatures*, Phys. Rev. 113, 1046 (1959)
- J. Zou & A. A. Balandin, *Phonon heat conduction in a semiconductor nanowire*, J. Appl. Phys. 89, 2932 (2001) — τ formulation and the mass-variance parameter
- Sommerfeld et al. 2003, 10.1103/PhysRevB.67.174520 — YBCO fit reference
- 10.1103/PhysRevLett.62.2317 — effective mass and carrier density; 10.1103/PhysRevB.80.064505 — oxygen Raman band
