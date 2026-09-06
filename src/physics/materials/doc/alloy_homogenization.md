# Alloy Homogenization: Formula Alloys from the Pure-Metal Roster {#physics_materials_alloy_homogenization}

**Date:** 2026-08-25
**Purpose:** What `material::Alloy` mixes and how, the validation on record, and — the part to read first — for which alloys the result can be trusted
**Module:** `src/physics/materials`
**Implementation:** `cl_Material_Alloy.{hpp,cpp}`, invoked by `MaterialFactory::create_material` for any label that parses as a composition

This document replaces the February 2026 design notes `alloy_mechanical_mixing.md` and
`alloy_transport_mixing.md`; their proposals are now implemented, and this document describes
the current behavior.

---

## 1. What it is

An `Alloy` is built from a composition formula:

```
material Sn60Pb40 --rrr 10
material Fe71Cr19Ni10 --rrr 1.13
```

```cpp
MaterialFactory tFactory ;
Material * tSolder = tFactory.create_material( "Sn60Pb40", 10.0 );   // label, RRR
```

The numbers are **integer mass percent**: `to_pair` in `src/core/stringtools.hpp` parses letters followed by digits. Every element must be one of the pure metals in the roster
(Al, Cr, Fe, Ni, Cu, Ag, In, Sn, Pb). The factory instantiates one `Metal` object per component
and hands them to the alloy, which owns them. The **last component is the balance**: its
percentage is recomputed as 100 minus the others (`Alloy::set_components`,
cl_Material_Alloy.cpp:149), so `Fe70Cr19Ni10` silently becomes 70/19/11.

RRR is **mandatory and must exceed 1**; otherwise `BELFEM_ERROR` fires. It is the only alloying
effect the model knows about, and it is fixed at construction: the inherited `Material::set_RRR`
is not overridden and raises an error. The
`material` executable defaults it to 10 for a formula and to 100 for a pure metal.

Three object details matter to callers. It is typed `MaterialType::PureMetal`, so
FEM code treats it like a metal. Its components are constructed without lookup tables of their
own. Unlike `Metal`, whose `rho( T, B, beta )` can evaluate Kohler's rule on the fly,
`Alloy::rho( T, B, beta )` and `lambda( T, B, beta )` read the tabulated database and assert
its presence; construct the alloy with tables enabled if you will ask for field dependence.

The alloy is a `SplineLookupTable`: every property is tabulated once at construction on a 4 K
grid from 0 K up to the first 4 K multiple at or above the lowest `T_max` of its components, and
served from cubic splines afterwards.
The reference density is the mass-weighted harmonic mean of the components' room-temperature
densities; the molar mass is the molar-fraction average.

## 2. What is mixed and how

Mass fractions Y_k are fixed. At each temperature, the components' densities are used to recompute volume fractions v_k(T) = (Y_k / ρ_k(T)) / Σ_j (Y_j / ρ_j(T))
(`Alloy::update_volume_fractions`).

### 2.1 Elastic moduli: self-consistent scheme

Each component contributes isotropic stiffness from its own E_k(T), ν_k(T):
K_k = E_k / (3(1 − 2ν_k)), G_k = E_k / (2(1 + ν_k)). The effective K and G are found by the
self-consistent (Eshelby–Hill) scheme for spherical inclusions:

1. Voigt and Reuss bounds, K_V = Σ v_k K_k, 1/K_R = Σ v_k / K_k (likewise for G), give the
   starting point K = (K_V + K_R)/2 at the first temperature; every later temperature starts
   from the previous converged value.
2. With the current (K, G) as the reference medium, the Eshelby tensor for a sphere is
   S = γ δδ + δ' (sym), with γ = K / (3K + 4G) and δ' = 3(K + 2G) / (5(3K + 4G)) (Qu & Cherkaoui
   2006, Eq. 4.3.40). The Hill polarization tensor is P = S : M, with M the compliance of the
   reference medium.
3. For every component, compute the strain-concentration tensor A_k = [I + P : (C_k − C)]⁻¹
   and update C = Σ v_k C_k : A_k. K and G are read off the isotropic result
   (C₁₁₁₁ = K + 4G/3, C₂₁₂₁ = G), under-relaxed with ω = 0.99.
4. Iterate until the relative change in (K, G) falls below 1e-8 (`Alloy::self_consistent_KG`).

E and ν follow from the converged pair: E = 9KG / (3K + G), ν = (3K − 2G) / (6K + 2G). The
arithmetic uses GPa for scaling; the result is stored in Pa.

Spherical inclusions describe an equiaxed microstructure — cast or sintered alloys. Directional
or heavily worked microstructures would need an ellipsoidal Eshelby tensor and extra parameters,
which the class does not have.

### 2.2 Thermal expansion: Levin (Schapery) with the self-consistent K

With α_V = Σ v_k α_k and α_R = K_R Σ v_k α_k / K_k, the effective coefficient is

α = α_V + (1/K − 1/K_V) / (1/K_R − 1/K_V) · (α_R − α_V),

using the self-consistent K (cl_Material_Alloy.cpp:372). The
formula is exact for two phases and a tight estimate for more; when all K_k coincide the
interpolation is undefined and the code takes the mean of the two bounds. At 0 K, α, c_p and λ are
set to zero.

### 2.3 Specific heat and Sommerfeld coefficient: mass-fraction averages

c_p = Σ Y_k c_p,k(T) — exact for an extensive property. The electronic coefficient
γ = Σ Y_k γ_k is stored as a constant and sets the slope of the c_p spline at 0 K.

### 2.4 Electrical resistivity: Matthiessen with one residual term

ρ(T) = ρ₀ + Σ v_k ρ_i,k(T)

where ρ_i,k is each metal's phonon (Bloch–Grüneisen) resistivity from `Metal::rho_i_custom`, and
ρ₀ is a single residual resistivity for the whole alloy, fixed from the caller's RRR at room
temperature:

ρ₀ = Σ v_k ρ_i,k(293.15 K) / (RRR − 1)

(`Alloy::set_components`). One ρ₀ for all phases is a deliberate simplification: it reproduces
the macroscopic RRR by construction and needs no phase-diagram data; it loses the distinction
between the residual scattering inside a Pb-rich and a Sn-rich lamella.

Linear (volume-fraction) mixing of the intrinsic parts was chosen over a Bruggeman effective
medium after a comparison with the Hariharan et al. 1979 data for 60Sn–40Pb: linear mixing
tracks the measured ρ(T) from 4 K to 300 K, Bruggeman's spherical-inclusion geometry suppresses
the slope by 40–50 % on a lamellar eutectic. Nordheim's rule is not used — it applies to random
solid solutions below the solubility limit, not to the two-phase solders the class was written for.

### 2.5 Thermal conductivity: linear mixing of thermal resistivities

λ(T) = 1 / (w₀ + Σ v_k w_i,k(T)), with w₀ = ρ₀ / (L₀ T) the residual (Wiedemann–Franz) term
and w_i,k = hust(P_k, T) each metal's intrinsic phonon-limited thermal resistivity from its Hust
coefficients (`fn_hust.hpp`). The pure-metal Hust form has a third, interaction term w_i0 whose
empirical constant is unknown for a mixture; it is dropped, which slightly underestimates λ where
w₀ and w_i are comparable and is immaterial when the alloy's large ρ₀ dominates.

### 2.6 Field dependence: per-phase Kohler, mixed on the alloy's own baseline

`Alloy::populate_rho_database` (cl_Material_Alloy.cpp:529)
tabulates ρ(T, B, β) exactly as `Metal` does, but evaluates each component's Kohler curve at the
similarity parameter it would see inside the alloy:

S_k = ρ_ref,k / (ρ₀ + ρ_i,k(T)),   δ_k = kohler_k(B, S_k, β),   ρ(T, B, β) = ρ(T) · (1 + Σ v_k δ_k).

Because ρ₀ of an alloy is far larger than any pure metal's, S_k is small and the
magnetoresistance is weak — the expected behavior of a dirty conductor, though not one this
path has been checked against field data. Every component
must provide a Kohler curve; all nine roster metals do. λ(T, B, β) follows by Wiedemann–Franz
scaling, λ(T) ρ(T) / ρ(T, B, β), the same closure `Metal` uses. The table is cached as
`<label>_RRR<n>.hdf5` in the run directory like a metal's, probed for the current format on
load and rebuilt if it predates it.

### 2.7 Stored curves

| Property | 0 K start condition | End condition |
|---|---|---|
| E, ν | clamped, zero slope | parabolic |
| c_p | clamped, slope γ | parabolic |
| ρ | clamped, zero slope | parabolic |
| α, λ | parabolic | parabolic |

## 3. What the homogenizer can and cannot do

The model is a **two-phase composite rule plus one residual-resistivity knob**. It is
trustworthy exactly when the real alloy is one of those two things, and it fails when alloying
creates something no constituent has: a new crystal structure, a magnetic state, an
intermetallic. The dividing line in one sentence: *the method holds when the alloy keeps its
constituents' crystal structure and magnetic state; the residual resistivity is the only alloying
effect it models.*

**Trustworthy**

- **Eutectic solders — Sn–Pb, Sn–Ag, In–Sn, Sn–Pb–In.** Pb and Sn are nearly insoluble in each
  other, so a solder *is* a fine composite of two almost-pure phases; the mixture rule is the
  right physics, not an approximation. For Sn63Pb37: α ≈ 24.6×10⁻⁶ from the mix against ≈ 25
  measured, E ≈ 30–35 GPa against ≈ 30, ρ ≈ 1.5×10⁻⁷ Ω·m against 1.45×10⁻⁷. This is the case
  the class was written and validated for.
- **Cu–Ag.** Negligible mutual solubility: a two-phase composite, same reasoning.
- **Cu–Sn bronze, Al–Cu.** Solid solutions in the majority element's own fcc lattice; the matrix
  keeps its phase, and the residual resistivity is exactly what RRR ≈ 1.1 (bronze) or ≈ 2–3
  (2xxx aluminium) describes.
- **Cu–Ni, Ni–Cr, and Ni-base superalloy matrices (e.g. Ni72Cr16Fe8).** Fully miscible fcc
  solutions; Cr, or a copper majority, suppresses nickel's ferromagnetism, so the class's lack of
  a magnetic channel is correct here rather than accidental. Expect chromium's small α to pull
  the mixture ≈ 15 % low.

**Right mechanically, wrong magnetically**

- **Ferritic Fe–Cr (e.g. Fe83Cr17).** A bcc solution in bcc iron, so α, E, c_p mix well — but
  the alloy is ferromagnetic and the class returns μ₀. `Alloy` has no magnetic channel at all,
  whatever its components are.

**Not to be trusted**

- **Austenitic stainless steels (Fe71Cr19Ni10 ≈ 304).** The alloy is fcc; Fe and Cr are bcc.
  Whatever depends on the crystal structure rather than on the atoms is systematically off, see
  the worked example below.
- **Invar (Fe64Ni36) and other magnetovolume alloys.** Measured α ≈ 1×10⁻⁶ against a mixture
  value of ≈ 12×10⁻⁶ — a factor ten.
- **9 %-Ni cryogenic steels, intermetallics (Nb₃Sn).** Martensitic or ordered phases; not
  mixtures of anything in the roster.

### Worked example: Fe71Cr19Ni10 with RRR 1.13 against NIST 304

Reference values from the NIST cryogenic material fits
(https://trc.nist.gov/cryogenics/materials/304Stainless/304Stainless_rev.htm):

| T [K] | k [W/(m·K)] | c_p [J/(kg·K)] | E [GPa] | ΔL/L₂₉₃ [%] | α [10⁻⁶/K] |
|---|---|---|---|---|---|
| 4 | 0.27 | 2.1 | 210 | −0.297 | — |
| 77 | 7.9 | 205 | 214 | −0.280 | 7.0 |
| 293 | 15.1 | 471 | 200 | 0 | 15.4 |

What the mixture gives, and why:

- **E** lands near 215–220 GPa: ≈ 10 % high, the Voigt–Reuss–Hill average of bcc Fe, bcc Cr and
  fcc Ni. Ballpark.
- **c_p above ~50 K** is fine (Dulong–Petit does not care about the phase). **Below ~20 K it is
  low by a factor of 5–6**: 304's electronic Sommerfeld coefficient is ≈ 0.5 J/(kg·K²), the
  mass-weighted γ of the components ≈ 0.08. The austenitic density of states is not a mixture
  property. For a quench study at 4 K this is the number that matters most.
- **α** comes out ≈ 10.7×10⁻⁶ at room temperature against 15.4, and the 293→4 K contraction
  ≈ 0.18 % against 0.30 %: 30–40 % low, the fcc-versus-bcc difference again.
- **ρ and k** are within a factor two: the RRR fixes ρ ≈ 7–8×10⁻⁷ Ω·m (304: 7.2×10⁻⁷), and
  Wiedemann–Franz then gives the electronic part of k; the missing remainder is the lattice
  conductivity that dominates in a dirty alloy.
- **μ = μ₀** — correct for 304, by accident.

## 4. Validation on record

- 60Sn–40Pb: ρ(T) against Hariharan et al. 1979, 4–300 K (the linear-versus-Bruggeman decision
  above). No other alloy has been checked against measurement through this path; the 304 table
  above is an estimate from the pure-metal values, not a run.
- No test in `check-fast` exercises the alloy path. A test constructing `Sn60Pb40` and
  `Fe71Cr19Ni10` and asserting loose bounds against the values above would turn "ballpark" into an
  executable statement and would catch a regression in any component metal at once.

## 5. Known limitations

- One ρ₀ for all phases (§2.4); per-phase residual terms would need solubility limits and
  Nordheim coefficients.
- The Wiedemann–Franz closure assumes L = L₀; in alloys at low temperature L/L₀ ranges roughly
  0.6–1.9. Inherited from the pure-metal model.
- The w_i0 interaction term of the Hust form is dropped (§2.5).
- No magnetic channel: every alloy is non-magnetic.
- No alloy-specific Debye temperature; the phonon resistivities are those of the components.
- Spherical inclusions in the elastic scheme; the transport mixing was chosen for lamellar
  eutectics on the strength of the one Sn–Pb comparison in §4, the elasticity remains a
  spherical-inclusion approximation for them.

## References

- Eshelby, J. D. (1957), Proc. R. Soc. London A 241, 376–396
- Hill, R. (1965), J. Mech. Phys. Solids 13, 213–222
- Qu, J. & Cherkaoui, M. (2006), *Fundamentals of Micromechanics of Solids*, Wiley — Eshelby tensor for spheres, Eq. (4.3.40)
- Levin, V. M. (1967), Mech. Solids 2, 58–61; Schapery, R. A. (1968), J. Compos. Mater. 2, 380–404
- Voigt, W. (1928), *Lehrbuch der Kristallphysik*; Reuss, A. (1929), Z. Angew. Math. Mech. 9, 49–58
- Hust, J. G. & Lankford, A. B. (1984), NBS Special Publication 260-90 — thermal conductivity form
- Hariharan, Y. et al. (1979), Pramana 13, 117–125 — 60Sn–40Pb resistivity
- Bruggeman, D. A. G. (1935), Ann. Phys. 416, 636–664; Nordheim, L. (1931), Ann. Phys. 401, 607–640
- NIST Cryogenic Material Properties, 304 stainless steel — reference fits in §3
