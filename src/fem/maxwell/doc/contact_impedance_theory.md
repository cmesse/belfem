# Contact-Impedance Interface: Theory and Derivation {#fem_maxwell_contact_impedance_theory}

**Date:** 2026-04-09
**Purpose:** Derive the surface bilinear form for a thin resistive contact layer (rint) collapsed to a zero-thickness interface in the H-formulation, and clarify what this buys us relative to the existing PENTA6TS thin-shell treatment.
**Module:** fem/maxwell
**Status:** Theory note. No implementation yet.
**Related:** Schnaubelt et al. 2023 §III

---

## 1. Physical setup

The rint layer sits between the silver coating and the YBCO film in a REBCO tape stack. Its role is to model the electrical contact resistance between the two conductors.

| Property | Value |
|---|---|
| Resistivity `rho_rint` | ~1e-3 Ohm*m |
| Thickness `h_rint` | ~0.1 um = 1e-7 m |
| Contact resistance per unit area `R_ct = rho * h` | ~1e-10 Ohm*m^2 |
| Through-thickness stiffness `rho / h` | ~1e4 Ohm |

The dominant physics is **through-thickness current transfer**: current flows from Ag through rint into YBCO (normal to the layer plane). In-plane currents within the rint are negligible because the layer is too thin and too resistive for significant tangential current flow.

## 2. Bulk weak form in the rint layer

Starting from the H-formulation weak form in the rint volume `Omega_rint` (thickness h, resistivity rho):

```
integral_Omega_rint  mu_0 (dH/dt) . v  dV   +   integral_Omega_rint  rho (curl H) . (curl v)  dV  =  0
```

The first term (mass) is O(mu_0 * h) ~ 1e-13 and negligible for thin layers. The stiffness term dominates.

## 3. Thin-shell approximation (N = 1)

Following Schnaubelt 2023 §III (Equations 6-8), the thin-shell approximation (TSA) collapses the rint volume to a surface `Gamma` (the rint midplane). The H field is discretized through the virtual thickness using first-order Lagrange basis functions:

```
H(u, v, w) = H^-_t(u, v) * Xi_0(w)  +  H^+_t(u, v) * Xi_1(w)
```

where `w in [0, h]`, `Xi_0(w) = (h - w)/h`, `Xi_1(w) = w/h`, and `H^-`, `H^+` are the tangential H traces on the bottom (Ag) and top (YBCO) sides of the interface.

### 3.1 Through-thickness contribution (dominant)

The curl through the thickness gives the tangential current:

```
(curl H)_u  ≈  -dH_v/dw  =  -(H^+_v - H^-_v) / h
(curl H)_v  ≈   dH_u/dw  =   (H^+_u - H^-_u) / h
```

These are O(1/h) and dominate the curl. Integrating the stiffness through the thickness with the 1D Lagrange FE matrices:

```
integral_0^h  (dXi_i/dw)(dXi_j/dw) dw  =  (1/h) * [1  -1]
                                                     [-1   1]
```

The resulting surface integral is:

```
a_through(H, v)  =  integral_Gamma  (rho / h) * [[H_t]] . [[v_t]]  dS
```

where `[[H_t]] = H^- - H^+` is the tangential H jump across the interface.

**This is the contact-impedance term.** It scales as `rho / h`.

### 3.2 In-plane contribution (negligible for rint)

The tangential curl (in-plane derivatives) gives the normal current component. Integrating through the thickness:

```
integral_0^h  Xi_i * Xi_j  dw  =  (h/6) * [2  1]
                                            [1  2]
```

The resulting surface integral scales as `rho * h` ~ 1e-10, which is 14 orders of magnitude smaller than the through-thickness term and is negligible for the contact resistance case.

### 3.3 Mass contribution (negligible for rint)

```
a_mass(H, v)  =  integral_Gamma  mu_0 * (h/6) * [2  1] * (dH/dt) . v  dS
                                                   [1  2]
```

The mass contribution scales as `mu_0 * h` ~ 1e-13 and is negligible.

### 3.4 The simplified contact-impedance interface

For the rint application, retaining only the dominant through-thickness term:

```
a_rint(H, v)  =  integral_Gamma  (rho_rint / h_rint) * [[H_t]] . [[v_t]]  dS
```

This is the **simplified surface bilinear form** for the collapsed rint interface. The no-consistency-term assumption rests on treating this as a natural (Robin-type) boundary condition derived from the energy functional. Juntunen & Stenberg 2009 ("Nitsche's method for general boundary conditions") supports this for scalar elliptic problems, and the physical derivation (§3.1) suggests it extends to the H(curl) case. A formal proof for BELFEM's specific H(curl) conductor-conductor interface is still missing. Either a proof or a numerical convergence study confirming optimal order is an implementation gate. Until one of the two exists, treat the absence of consistency terms as a **plausible inference, not a settled result**.

## 4. Scaling: the matrix entry goes as rho/h

Two quantities with reciprocal thickness scaling must not be confused:

- `rho * h` is the **contact resistance per unit area** `R_ct`, the physically meaningful, measurable parameter (§10).
- `rho / h` is the **through-thickness impedance**, the reciprocal of the layer's conductance per unit area and the coefficient that enters the matrix.

Both the bulk discretization and the collapsed interface produce the same matrix scaling:

- **Bulk rint layer:** K_entry ~ (rho / h) * area. Scaling: `rho / h`.
- **TSA surface interface:** a_rint ~ (rho / h) * [[H_t]] * [[v_t]] (§3.1). Scaling: `rho / h`.

The scaling is the **same**, so the TSA does not change the conditioning ratio. For rint with `rho = 1e-3` and `h = 1e-7`, the matrix entry scales as `1e4` in both forms, and the conditioning ratio against adjacent Ag (`rho/h ~ 2e-3`) remains `~10^7`.

## 5. What the TSA surface form provides

The surface form does not improve conditioning. Its advantages are:

| Benefit | Bulk rint layer | TSA surface interface |
|---|---|---|
| DOFs in rint | 6 edge_h per element (3 bottom + 3 top) | 0 (collapsed to interface) |
| Coupling interfaces | 2 (rint-Ag + rint-YBCO) | 1 (Ag-YBCO directly) |
| Nitsche penalty parameter | Required (alpha must bridge rint-Ag contrast) | Not needed (the impedance is the coupling) |
| Matrix extreme values | Entire rint block has K ~ 1e4 | Only interface DOFs have ~ 1e4 contribution |
| Conditioning ratio | ~10^7 (rint block vs Ag block) | ~10^7 (interface DOFs vs Ag interior DOFs) |
| Localization | Extreme values distributed across all rint DOFs | Extreme values localized to interface DOFs only |

The key practical advantage is the third table row: eliminating the separate penalty parameter. The current `h_ghost()` kernel in `mt_maxwell_h.cpp` already bounds alpha via a regularized harmonic mean with `k_reg`, so the interface does not produce runaway off-diagonal entries. The regularized harmonic mean is a *stabilization choice* that trades accuracy for conditioning. The TSA contact-impedance form eliminates the issue at the source: the physical impedance `rho/h` is the coupling coefficient, exactly right by construction, with no tunable parameter and no accuracy-conditioning tradeoff.

## 6. Comparison with existing code paths

### 6.1 vs. the existing thin-shell conductor model

The existing PENTA6TS thin-shell model (the `h_picard` / `h_newton_*` kernels driving the thin-shell branch of `calculator::MaxwellData`) is the N=1 TSA. It assembles the full surface bilinear form: mass (§3.3), in-plane stiffness (§3.2), and through-thickness stiffness (§3.1). For the rint application, only the through-thickness term (§3.1) is numerically significant.

The difference between the existing thin-shell model and the collapsed interface:
- Thin-shell model: rint has its own PENTA6TS block with 6 edge_h DOFs per element. The through-thickness and in-plane terms are assembled as a bulk element. Two coupling interfaces (rint-Ag, rint-YBCO) each need Nitsche or condensation.
- Collapsed interface: rint has no block, no DOFs. The through-thickness term `(rho/h) [[H_t]] [[v_t]]` is assembled directly on the Ag-YBCO interface facets. One coupling interface, no separate penalty.

### 6.2 vs. h_ghost() (continuity Nitsche)

`h_ghost()` in `mt_maxwell_h.cpp` weakly enforces tangential H continuity:

```
[[H_t]] = 0     (continuity, not impedance)
```

The penalty form is: `alpha * [[H_t]] . [[v_t]]` plus SIPG consistency/adjoint terms.

The contact-impedance interface enforces a **different** condition:

```
(rho_rint / h_rint) * [[H_t]] = E_t     (impedance / Robin)
```

The coefficient is the physical contact impedance `rho / h`, not a tunable stabilization parameter. Consistency and adjoint terms remain an open question (see §3.4 above); current evidence suggests they are not needed, but this has not been proven for the H(curl) case.

The matrix forms look structurally similar because both are `coefficient * [[H_t]] . [[v_t]]`, with these differences:
- `h_ghost()` uses a tunable `alpha` (now a regularized harmonic mean with `k_reg` in `h_ghost()`, `mt_maxwell_h.cpp`) that must be "large enough" for stability
- The contact-impedance form uses `rho / h`, the physical impedance, as the coupling coefficient
- `h_ghost()` adds SIPG terms (`-{rho curl H} . [[v]]` and adjoint) for optimal convergence of the continuity condition
- The contact-impedance **likely** does not need SIPG terms (pending the derivation in §3.4), because it enforces a Robin condition rather than a Dirichlet-like continuity

## 7. Matrix structure for BELFEM implementation

The contact-impedance kernel operates on an **interface sideset** between two conductor blocks (Ag and YBCO). "No DOFs in rint" means no rint-owned block DOFs because the rint block does not exist. The H+ and H- traces at the interface still carry edge_h DOFs from the Ag and YBCO conductor blocks. These DOFs must remain distinct, rather than merged, so the jump `[[H_t]] = H_m - H_s` is representable. This is the same distinct-trace requirement as in Schnaubelt 2023 §III and Alves 2024.

At each integration point on the interface, the kernel assembles:

```cpp
// H_m, H_s : tangential edge-h DOFs from master (Ag) and slave (YBCO)
// E_m, E_s : edge basis functions evaluated at the integration point
// rho_h    : rho_rint / h_rint (contact impedance coefficient)
// wdS      : integration weight * surface Jacobian

// The jump is [[H_t]] = E_m * H_m - E_s * H_s
// The bilinear form is: rho_h * [[H_t]] . [[v_t]]

Kmm +=  rho_h * trans(Em) * Em * wdS ;   // master-master
Kms += -rho_h * trans(Em) * Es * wdS ;   // master-slave (coupling)
Ksm += -rho_h * trans(Es) * Em * wdS ;   // slave-master (coupling)
Kss +=  rho_h * trans(Es) * Es * wdS ;   // slave-slave
```

This is assembled into the stiffness matrix `K()` (not `M()`), because it represents resistive dissipation in the contact layer. The mass-matrix contribution is negligible for thin rint and is omitted.

This is structurally identical to the penalty part of `h_ghost()`, except that it omits the consistency and adjoint terms (the `Dm`/`Ds` curl-operator terms in `h_ghost()`). A new kernel function `h_contact_impedance()` in `mt_maxwell_h.cpp` would be cleaner than reusing `h_ghost()` with modified alpha.

## 8. Relation to the buffer phi-formulation

The buffer layer (phi formulation) and the rint contact impedance solve **different problems**:

| Aspect | Buffer (phi formulation) | Rint (contact impedance) |
|---|---|---|
| Physics | True insulator, J = 0 everywhere | Resistive contact, J_through ≠ 0 |
| DOFs | Nodal phi (scalar potential) | None (collapsed to interface) |
| Bulk operator | mu_0 * grad(phi) . grad(v) (Laplace) | None |
| Interface coupling | Static condensation (edge-to-node) | (rho/h) * [[H_t]] . [[v_t]] (impedance) |
| Current transfer | Zero | Finite, determined by rho/h |
| Cut topology | Rerouted through buffer node (see `buffer_cut_topology.md`) | Not affected **for the local same-loop case** (rint sits between two conductors in the same current loop). In NI-coil geometries where radial contact currents change the current partition between turns, Schnaubelt 2023 shows an additional free cut coefficient IC2 is needed; that case is outside our current scope. |

The two can coexist in the same tape stack:
```
Cu - Ag - Has - [buffer, phi] - YBCO - [rint, contact impedance] - Ag - Cu
```

The buffer splits the conductor into two halves (requiring cut rerouting). The rint provides the resistive bridge for current transfer between Ag and YBCO (if present). Without rint, the YBCO half is only connected to the terminals via the cut constraint.

## 9. Implementation path

The implementation is more than adding one kernel: the kernel itself is small (§9.2), but the Maxwell machinery does not yet provide the support code for a collapsed interface (§9.1), and that support code dominates the effort.

### 9.1 Plumbing prerequisites (the larger effort)

A new domain type `InterfaceTsCond` already exists in `en_DomainType.hpp:62`, but it is not wired through any dispatch tables. The following sites all need a new case or branch:

| Site | File:line | What's missing |
|---|---|---|
| String parser | `en_DomainType.cpp:97` | No "contact impedance" or "interfacetscond" string entry |
| Topology production | `cl_Topology.cpp:270` | Not produced by `sideset_type(master, slave)` |
| Domain parsing | `cl_FEM_Domain.cpp:28` | Not handled in the input-section dispatch |
| IWG dispatch | `cl_IWG_Maxwell.cpp:558` | Not dispatched to any kernel function |
| DOF allocation | `cl_Maxwell_FieldList.cpp:411` | Not given DOFs in the sideset DOF table |

Additionally:

- **Material/thickness on an interface sideset.** Maxwell's material assignment is block-based (`cl_MaxwellFactory.cpp:2226`), and thin-shell material/thickness metadata is indexed to existing layer blocks (`cl_ThinShell.hpp:45`, `cl_ThinShellFactory.cpp:109`). A collapsed rint (no block) needs a **new interface-property path**: either a dedicated sideset-property container, or a forward reference from the interface sideset to the (now-absent) rint layer's material and thickness. This is a design decision that must be resolved before the kernel can be written.

- **ThinShellFactory must skip PENTA6TS extrusion for rint.** Currently, every layer in the stack gets a block of PENTA6TS elements. A contact-impedance layer must not be extruded. Instead, the ThinShellFactory creates an interface sideset between the adjacent conductor layers (or reuses the existing ghost-facet sideset if one exists).

### 9.2 The kernel itself (the smaller effort)

1. **New kernel function:** `h_contact_impedance()` in `mt_maxwell_h.cpp`. It takes the interface facet's master/slave edge bases, computes `(rho/h) [[H_t]] [[v_t]]`, and assembles into Kmm/Kms/Ksm/Kss. Whether consistency/adjoint terms are needed is pending derivation (see §3.4). The `rho` and `h` come from the rint layer's material and thickness via the new interface-property path from §9.1.

2. **No changes to the cut or terminal machinery** for the local same-loop case (rint sits between two conductors in the same current loop). For NI-coil geometries where radial contact currents change the current partition between turns, an additional free cut coefficient would be needed (Schnaubelt 2023); that is outside current scope.

### 9.3 Effort estimate

| Component | Effort | Risk |
|---|---|---|
| Plumbing (§9.1): domain type wiring, topology, DOFs | 1-2 days | Medium (many files, all small) |
| Plumbing (§9.1): interface-property path design + implementation | 1 day | High (new design decision) |
| Plumbing (§9.1): ThinShellFactory skip-extrusion for rint | Half a day | Medium |
| Kernel (§9.2): `h_contact_impedance()` | Half a day | Low (simple bilinear form) |
| Consistency-term derivation (§3.4 gate) | 1 day | Medium (may not be needed, but must be checked) |
| Testing: analytical test problem + tape-stack validation | 1 day | Low |
| **Total** | **4-6 days** | |

## 10. Measured contact resistance values

The working value ρ = 1e-3 Ω·m used in the conditioning estimates above sits at the **extreme high end** of the measured range. The actual values depend critically on which interface is being modeled.

### 10.1 Internal Ag/YBCO interface (within a single tape)

This is the metallurgical bond between the silver cap and the YBCO film, the "rint" in our tape model. Measured contact resistance per unit area R_ct:

| Source | R_ct [Ω·m²] | Notes |
|---|---|---|
| Hayasaka & Ito 2019 | 7e-12 to 3.7e-11 | Contact-probing current transfer length method, temperature-dependent |
| Soldered overlap measurements | 2.5e-12 to 5e-12 | Indium-soldered overlap |
| As-supplied coated conductors | ~2.5e-12 | ~25 nΩ·cm², Wimbush & Strickland data |
| Cu-stabilized REBCO lap joints | ~3.6e-12 | ~36 nΩ·cm², dominated by tape interface itself |

**Typical range: R_ct ≈ 1e-12 to 1e-10 Ω·m²**

Converting to bulk resistivity for a modeled layer thickness h = 0.1 μm:
- R_ct = 1e-12 → ρ = R_ct / h = 1e-12 / 1e-7 = **1e-5 Ω·m**
- R_ct = 1e-11 → ρ = 1e-11 / 1e-7 = **1e-4 Ω·m** (typical)
- R_ct = 1e-10 → ρ = 1e-10 / 1e-7 = **1e-3 Ω·m** (worst case: degraded / aged contact)

### 10.2 Turn-to-turn contact in NI coils (Schnaubelt's T2TCL)

This is mechanical contact between tape surfaces under winding pressure, a completely different physical interface. Typical values R_ct ≈ 1e-9 to 1e-7 Ω·m² (10-100 μΩ·cm²), i.e. **100-10000× higher** than the internal Ag/YBCO interface. This is the application where the contact-impedance TSA (Schnaubelt 2023) is essential.

### 10.3 Implications for conditioning

The conditioning ratio ρ/h against adjacent Ag (ρ_Ag/h_Ag ≈ 2e-3):

| Scenario | ρ_rint [Ω·m] | ρ/h | Contrast vs Ag | Assessment |
|---|---|---|---|---|
| **Typical internal contact** | 1e-4 | 1e3 | ~5e5 | MUMPS handles this routinely |
| Moderate contact | 1e-3 | 1e4 | ~5e6 | Challenging, may need regularization |
| NI coil T2TCL | 1e-1 to 1 | 1e6 to 1e7 | ~5e8 to 5e9 | Requires TSA or resistivity cap |

## 11. Assessment: is this worth implementing now?

### The resistivity cap matches the typical contact resistance

The resistivity cap at ρ_eff = 1e-4 Ω·m is both Alves 2022b's literature limit for what the pure H-formulation can handle and the physical value for typical REBCO tapes. It coincides with the measured Ag/YBCO interface resistivity for a typical coated conductor (R_ct ≈ 1e-11 Ω·m², h ≈ 0.1 μm → ρ ≈ 1e-4 Ω·m).

This means:
- The "cap" is the right material property for a standard tape, not an approximation
- The conditioning ratio at ρ = 1e-4 is ~5e5, which MUMPS handles without difficulty
- The regularized harmonic mean in `h_ghost()` (`mt_maxwell_h.cpp`) is a reasonable stabilization at this contrast level
- Junie's empirical test (stable convergence at -156 dB residual) confirms the solver copes

### When ρ = 1e-3 applies

A resistivity of 1e-3 Ω·m corresponds to R_ct ≈ 1e-10 Ω·m², the **upper extreme** of the measured range. This occurs in:
- Degraded or aged contacts
- Mechanically damaged interfaces (delamination)
- Deliberately high-resistance interfaces (e.g., certain buffer-stabilized architectures)

For these cases, the conditioning ratio rises to ~5e6, which is more challenging but may still be manageable with the existing regularized h_ghost() approach.

### When the contact-impedance formulation becomes necessary

The collapsed-interface formulation (§9) is worth the 4-6 day implementation effort only for:

1. **NI coil T2TCL (Schnaubelt's application).** R_ct ≈ 1e-9 to 1e-7 Ω·m², giving ρ/h ≈ 1e6 to 1e7. At this level, neither the resistivity cap nor the regularized h_ghost() is expected to be adequate, so the TSA is essential. This is the application where Schnaubelt developed the method.

2. **Multi-tape stacks with many contact layers.** Even if each individual rint has manageable conditioning, 50-100 rint layers in a coil stack multiply the DOF and coupling overhead. The collapsed interface eliminates this.

3. **Quench propagation studies** where the time-dependent current redistribution through the contact layer is the primary quantity of interest and the cap value is too coarse.

### Recommendation

**For the current work (single-tape transport current at 77 K):** the resistivity cap at ρ_eff = 1e-4 Ω·m is sufficient and physically correct for typical REBCO tapes. The existing h_ghost() with regularized harmonic mean handles the conditioning. For the internal Ag/YBCO contact, the regularized h_ghost() path may already solve the ghost problem.

**For future work (NI coils, quench, multi-tape stacks):** the contact-impedance formulation becomes necessary. The theory in §2-7 of this note and the implementation scope in §9 provide the starting point.

**Immediate action:** verify that the regularized h_ghost() at ρ = 1e-4 Ω·m produces the correct current-transfer behavior in the single-tape test case. If convergence is clean and the physics is right, the collapsed-interface implementation (§9) can be deprioritized indefinitely.

### References for measured contact resistance values

- Hayasaka, K. & Ito, S. (2019). "Evaluation of Interface Resistance in a REBCO Tape at Different Temperatures by Contact-Probing Current Transfer Length Method." IEEE Trans. Appl. Supercond. **29**(5).
- "Interface properties and failures of REBCO coated conductor tapes: Research progress and challenges." Superconductivity **7** (2023), 100055.
- Fleiter, J. et al. (2017). "Contact resistance between two REBCO tapes under load and load-cycles." arXiv:1701.00447.
- Lee, J. et al. (2021). "Investigation on the electrical contact resistance of soldered metal insulation REBCO coil." IEEE Trans. Appl. Supercond. **31**(5).

## 12. References

- Schnaubelt, E. et al. (2023). "Electromagnetic simulation of no-insulation coils using H-phi thin shell approximation." IEEE Trans. Appl. Supercond. **33**(5), 4900906. Section III (TSA derivation, Equations 6-8).
- Juntunen, M. & Stenberg, R. (2009). "Nitsche's method for general boundary conditions." Math. Comp. **78**(267), 1353-1374. (Robin boundary conditions are natural for scalar elliptic problems.)
- Messe, C. et al. (2023). "BELFEM: ..." Supercond. Sci. Technol. **36**, 114001. (BELFEM thin-shell architecture.)
- Alves, B. de Sousa et al. (2022b). "A thin-shell H-phi-formulation..." Supercond. Sci. Technol. **35**(2), 024001. (Multilayer thin-shell, max rho = 1e-4.)
