# Selective Nitsche Coupling at Thin-Shell Layer Interfaces

**Date:** 2026-03-14 (updated 2026-03-16 after three-AI review)
**Purpose:** Long-term architecture improvement for high-contrast thin-shell layers
**Status:** Design concept — needs literature verification and prototyping
**Module:** `src/fem/maxwell/`, `src/mesh/cl_ThinShellFactory.cpp`
**Reviewed by:** Claude (Opus 4.6), Codex, Grok — consensus: feasible, penalty-only is insufficient

---

## 1. Problem Statement

In the h-φ thin-shell formulation, adjacent material layers with extreme resistivity contrast (e.g., YBCO at near-zero effective ρ next to hastelloy at ~1.2 μΩ·m) share DOFs at their layer boundary. This forces C⁰ continuity of H_t across a jump of 5+ orders of magnitude in transport properties. The resulting system is severely ill-conditioned, causing GMRES stagnation.

The broken edge key in merge_hh accidentally decoupled layer interfaces by creating private edge DOFs per element. This acted as a crude DG method and masked the conditioning issue. With the correct edge key (proper deduplication), the ill-conditioning is exposed.

The current workaround — manual sub-layering of thick materials — adds DOFs for resolution but doesn't address the fundamental conditioning problem at high-contrast interfaces.

---

## 2. Discrete Weak Form

### Volume terms per layer

Within a single thin-shell layer k with permeability μ_k and resistivity ρ_k, the two discrete matrix contributions are:

**Mass (μ term, multiplies ∂h̃/∂t):**

$$\mathbf{M}_k = \int_{\Omega_k} \mathbf{E}^\mathsf{T} \, \mu_k \, \mathbf{E} \; \mathrm{d}V$$

**Stiffness (ρ term, multiplies h̃):**

$$\mathbf{K}_k = \int_{\Omega_k} \mathbf{C}^\mathsf{T} \, \rho_k \, \mathbf{C} \; \mathrm{d}V$$

The semidiscrete form for a single layer reads:

$$\delta\tilde{\mathbf{h}}^\mathsf{T} \left( \mathbf{M}_k \, \frac{\partial \tilde{\mathbf{h}}}{\partial t} + \mathbf{K}_k \, \tilde{\mathbf{h}} \right) = 0$$

With shared edges between layers, the assembly accumulates M_k and K_k contributions from adjacent layers into the same DOF positions. When ρ_k differs by orders of magnitude between neighbors, the K rows become extremely ill-conditioned.

### Nitsche interface terms

At a layer interface Γ_I between layers k (side +, ρ⁺) and k+1 (side −, ρ⁻), we break the shared DOFs and introduce coupling terms. The edge DOFs on each side are independent: h̃⁺ (layer k top) and h̃⁻ (layer k+1 bottom). Let z denote the through-thickness (normal) coordinate.

**Harmonic weighting:**

$$\delta^+ = \frac{\rho^-}{\rho^+ + \rho^-}, \qquad \delta^- = \frac{\rho^+}{\rho^+ + \rho^-}, \qquad \{\rho\}_{\mathrm{harm}} = \frac{2\rho^+\rho^-}{\rho^+ + \rho^-}$$

**Jump and weighted average operators:**

$$[\tilde{\mathbf{h}}_t] = \mathbf{E}^+ \tilde{\mathbf{h}}^+ - \mathbf{E}^- \tilde{\mathbf{h}}^- \Big|_{\Gamma_I}$$

$$\left\{\rho \, \frac{\partial \mathbf{h}_t}{\partial z}\right\}_\delta = \delta^+ \rho^+ \frac{\partial (\mathbf{E}^+ \tilde{\mathbf{h}}^+)}{\partial z}\bigg|_{\Gamma_I} + \delta^- \rho^- \frac{\partial (\mathbf{E}^- \tilde{\mathbf{h}}^-)}{\partial z}\bigg|_{\Gamma_I}$$

In the thin-shell geometry, (curl h)·n_s reduces to ∂h_t/∂z — the through-thickness derivative of the tangential field. This is the quantity that the existing 1D auxiliary FE problem discretizes.

**The Nitsche interface terms contribute to K (the stiffness operator acting on h̃, not on ∂h̃/∂t).** They have the same ρ-based, curl-curl structure as the volume stiffness.

#### Full symmetric Nitsche (three terms)

1. **Consistency** — arises from integration by parts; ensures the exact solution satisfies the discrete equations. Involves flux of the TRIAL function and jump of the TEST function:

$$a_1 = -\int_{\Gamma_I} \left\{\rho \, \frac{\partial \mathbf{h}_t}{\partial z}\right\}_\delta \cdot [\delta\tilde{\mathbf{h}}_t] \; \mathrm{d}S$$

2. **Adjoint consistency** — added to make the bilinear form symmetric. Involves flux of the TEST function and jump of the TRIAL function:

$$a_2 = -\int_{\Gamma_I} \left\{\rho \, \frac{\partial (\delta\mathbf{h}_t)}{\partial z}\right\}_\delta \cdot [\tilde{\mathbf{h}}_t] \; \mathrm{d}S$$

3. **Penalty** — stabilizes the discontinuity:

$$a_3 = \alpha \int_{\Gamma_I} [\tilde{\mathbf{h}}_t] \cdot [\delta\tilde{\mathbf{h}}_t] \; \mathrm{d}S$$

where the penalty parameter α accounts for both layers' thicknesses independently:

$$\alpha = \eta \left( \frac{\delta^+ \rho^+}{\Delta z^+} + \frac{\delta^- \rho^-}{\Delta z^-} \right)$$

and η > C_inv² (η ≈ 4–10 as a starting range for lowest-order edge elements on PENTA6TS; may need larger values for high aspect-ratio layers).

#### Nonsymmetric Nitsche (recommended for Phase 1)

Drop the adjoint consistency term a₂, keeping only consistency + penalty:

$$a_{\mathrm{Nitsche}} = a_1 + a_3$$

This is cheaper (one fewer interface integral), usually more robust on high-aspect-ratio elements, and is the standard choice in modern CutFEM/immersed boundary implementations. The trade-off: suboptimal convergence order (h^{p+1/2} instead of h^{p+1}), but for a first prototype this is sufficient. The symmetric version can be added later if needed.

**Penalty-only (a₃ alone) is NOT sufficient** — it loses consistency and requires excessively large η, reintroducing conditioning problems. The consistency term a₁ is essential even in the nonsymmetric variant.

#### Key simplification: this is 1D scalar DG, not H(curl) DG

The through-thickness derivative ∂h_t/∂z is a finite difference of **E** evaluations, and the jump [h_t] is also a difference of **E** evaluations. The curl operator **C** never appears in any Nitsche term. The in-plane curl (which gives the out-of-plane current j) is handled entirely within each element's volume K and does not interact with the inter-layer interface.

This means the well-established 1D interior penalty DG theory for scalar diffusion (Brenner & Scott 2008, Ch. 10) applies directly. The more complex H(curl) DG theory (Houston, Perugia & Schotzau 2004-2005) is not required, though it provides additional reassurance. **This significantly de-risks the approach.**

#### Single YBCO layer is sufficient

A single YBCO layer with DG on both sides (hastelloy below, silver above) does not need to be split into three sub-layers. The harmonic weighting naturally protects the low-ρ side: when ρ_YBCO ≈ 0 (superconducting), the penalty α ≈ η·ρ_YBCO/Δz_YBCO is very small on both interfaces. The YBCO element is essentially free to carry its own field, weakly coupled to neighbors. This is physically correct — the superconductor screens the external field and should not be tightly constrained by the substrate.

### Mass matrix at DG interfaces

The mass term ∫ E^T μ E dV does **not** need Nitsche coupling. With duplicated edges, each layer simply gets its own mass contribution. No consistency or penalty terms are needed for the L² inner product, since μ ≈ μ₀ everywhere in HTS tapes.

(If ferromagnetic substrates are modeled in the future, B_n = μH_n continuity would need additional treatment.)

### Coercivity

For the combined bilinear form (volume stiffness + Nitsche interface) to be coercive, η must satisfy η > C_inv². The harmonic weighting ensures the penalty doesn't reintroduce the contrast: for extreme contrast (ρ⁺/ρ⁻ ~ 10⁵), {ρ}_harm ≈ 2ρ⁻, dominated by the smaller resistivity.

### Why there is no right-hand-side contribution

In a classical penalty method for imposing a constraint g = g₀, one adds α(g − g₀) to the equations, which produces both a matrix term (α·g) and a RHS term (α·g₀). Here, the constraint we enforce is **weak continuity**: [h_t] ≈ 0 (the tangential field jump should vanish). Since g₀ = 0, all Nitsche terms are purely bilinear — they couple h̃ with δh̃ and contribute only to the stiffness matrix K. There is no load vector contribution. The physics is unchanged: we are not adding or removing any source. We are simply replacing the **strong** enforcement of H_t continuity (shared DOFs, which forces both sides to have identical values algebraically) with a **weak** enforcement (separate DOFs coupled through Nitsche terms, which penalizes and minimizes any discontinuity). In the limit η → ∞, the penalty drives [h_t] → 0 and recovers the original C⁰ formulation. For finite η, a small controlled discontinuity is permitted, which is what relieves the conditioning problem.

---

## 3. Concrete Implementation

### What `h_ts_hts()` currently does

`h_ts_hts()` (`mt_maxwell_h.cpp:454`) assembles a 6×6 element matrix for a single PENTA6TS element. The 6 DOFs are edge DOFs: edges 0–2 (bottom face) and edges 3–5 (top face). At each integration point k:

```cpp
aMatrices->M() += trans( E ) * E * ( w(k) * mu0 * dV(k) );     // mass
aMatrices->K() += trans( C ) * C * ( w(k) * rho * dV(k) );     // stiffness
```

The `E` matrix (3×6) from `EF_PENTA6TS::E()` has through-thickness interpolation built in via factors `mF(0, k)` (bottom weight) and `mF(1, k)` (top weight):

- Columns 0–2: bottom edge contributions, scaled by `mS[0–2] × mF(0, k)`
- Columns 3–5: top edge contributions, scaled by `mS[3–5] × mF(1, k)`

At the **bottom face**: `mF(0,·) = 1, mF(1,·) = 0` → only columns 0–2 are active.
At the **top face**: `mF(0,·) = 0, mF(1,·) = 1` → only columns 3–5 are active.

### Why E, not C

The penalty term penalizes the jump of the tangential field VALUE [h_t], which is interpolated by **E** (edge interpolation). The **C** operator (curl) does NOT appear in any Nitsche term. The through-thickness derivative ∂h_t/∂z is computed as a DIFFERENCE of E evaluations divided by layer thickness — it's an algebraic operation on edge DOFs, not a curl.

- **Penalty:** E^T E (field values at the interface)
- **Consistency:** (E_top − E_bot)/Δz (through-thickness gradient of field values)
- **No bare C anywhere** in the Nitsche terms

### The Nitsche penalty as a 12×12 interface matrix

At a DG interface between layers k and k+1, edges 3–5 of element k are SEPARATE objects from edges 0–2 of element k+1. Let z denote the through-thickness coordinate, with Δz⁺ the thickness of layer k and Δz⁻ the thickness of layer k+1 (these are generally DIFFERENT).

The penalty parameter must account for both sides independently:

$$\alpha = \eta \left( \frac{\delta^+ \rho^+}{\Delta z^+} + \frac{\delta^- \rho^-}{\Delta z^-} \right)$$

The penalty couples the 12 DOFs from both elements:

```
                      element k (6 DOFs)    element k+1 (6 DOFs)
                    ┌───────────────────┬───────────────────────┐
                    │  +E_top^T E_top   │  -E_top^T E_bot       │
K_penalty =  α · dS├───────────────────┼───────────────────────┤
                    │  -E_bot^T E_top   │  +E_bot^T E_bot       │
                    └───────────────────┴───────────────────────┘
```

Where `E_top` = E evaluated at the top face of element k (columns 0–2 zero, columns 3–5 active), and `E_bot` = E evaluated at the bottom face of element k+1 (columns 0–2 active, columns 3–5 zero).

Since E_top only has entries in columns 3–5, the upper-left block `E_top^T E_top` has only the 3×3 sub-block (rows/cols 3–5) nonzero. Similarly `E_bot^T E_bot` has only rows/cols 0–2 nonzero.

### Assembly via faces — no new sidesets needed

The existing face infrastructure provides master/slave element relationships at layer boundaries. By creating faces on the DG interface blocks, each face gives:
- **Master element** (layer k) and its local face index → access to edges 3–5 (top face)
- **Slave element** (layer k+1) and its local face index → access to edges 0–2 (bottom face)

The full 12×12 Nitsche matrix (diagonal + off-diagonal blocks) can be assembled in a single loop over these faces. Each face provides the element pair and face indices needed to evaluate E_top and E_bot and assemble all four blocks into the global matrix.

### Sparsity pattern

The off-diagonal blocks (master-slave, slave-master) introduce NEW couplings between DOFs that are not connected in the current CG assembly. The sparsity pattern computation in `SolverData` must account for the DG face element pairs when building the matrix graph — each DG face adds connectivity between the 6 DOFs of the master element and the 6 DOFs of the slave element.

### Through-thickness derivative for the consistency terms

The consistency terms involve ∂h_t/∂z at the interface. For a PENTA6TS element with linear through-thickness interpolation:

```
∂h_t/∂z |_k = (E_top h̃_k − E_bot h̃_k) / Δz_k
```

This is computable from each element's own 6 DOFs (difference of top and bottom face values divided by layer thickness). The weighted average uses one-sided derivatives from each side, with DIFFERENT Δz values:

```
{ρ ∂h_t/∂z}_δ = δ⁺ ρ⁺ (E_top h̃_k − E_bot h̃_k)/Δz⁺  +  δ⁻ ρ⁻ (E_top h̃_{k+1} − E_bot h̃_{k+1})/Δz⁻
```

Both sides are available through the face's master/slave elements.

### Surface integral dS

The Nitsche terms are surface integrals over the triangular face at the layer boundary. dS is the triangle area — NOT the volume dV (which includes layer thickness). For prismatic elements, dS ≈ dV / Δz. Alternatively, dS can be computed directly from the face Jacobian or from `mSurfaceIncrement` in the Calculator.

### Complete interface assembly pseudo-code

```cpp
// Loop over faces at DG layer boundaries
for ( Face * tFace : tDGFaces )
{
    Element * tMaster = tFace->master();   // layer k PENTA6TS
    Element * tSlave  = tFace->slave();    // layer k+1 PENTA6TS

    // evaluate E at the interface from each side
    // E_top: E of master element at its top face (mF(0)=0, mF(1)=1)
    // E_bot: E of slave element at its bottom face (mF(0)=1, mF(1)=0)

    // layer thicknesses (generally different!)
    real dz_plus  = master_layer_thickness;
    real dz_minus = slave_layer_thickness;

    // material properties (lagged ρ for Picard)
    real rho_plus  = master_material->rho( ... );
    real rho_minus = slave_material->rho( ... );
    real delta_plus  = rho_minus / ( rho_plus + rho_minus );
    real delta_minus = rho_plus  / ( rho_plus + rho_minus );

    // penalty parameter (per-element Δz)
    real alpha = eta * ( delta_plus * rho_plus / dz_plus
                       + delta_minus * rho_minus / dz_minus );

    // --- Penalty term (12×12, assembled into 4 blocks) ---
    // master-master:  +E_top^T E_top * alpha * dS
    // master-slave:   -E_top^T E_bot * alpha * dS
    // slave-master:   -E_bot^T E_top * alpha * dS
    // slave-slave:    +E_bot^T E_bot * alpha * dS

    // --- Consistency term ---
    // one-sided derivatives (per-element Δz):
    // ∂h_t/∂z|⁺ = (E_top h̃_k − E_bot h̃_k) / Δz⁺
    // ∂h_t/∂z|⁻ = (E_top h̃_{k+1} − E_bot h̃_{k+1}) / Δz⁻
    // weighted average: δ⁺ ρ⁺ (∂h_t/∂z|⁺) + δ⁻ ρ⁻ (∂h_t/∂z|⁻)
    // multiplied by jump [δh̃_t] and integrated over dS

    // --- Adjoint consistency term ---
    // transpose of consistency (ensures symmetry)

    // assemble all into global K matrix using DOF indices
    // from tMaster and tSlave elements
}
```

---

## 4. CG-DG Hybrid Strategy

### Within a layer (CG, unchanged)
- Adjacent PENTA6TS elements on the same layer share edges
- Well-conditioned (same material on both sides)

### At high-contrast layer boundaries (selective DG)
- Duplicate the boundary edges so each layer owns its own copy
- Create faces at these boundaries (provides master/slave element pairs)
- Couple through the full Nitsche bilinear form (all three terms)
- Use harmonic-weighted averages in both the flux terms AND the penalty

### Contrast threshold

Classify interfaces based on worst-case material contrast (from material definitions, not instantaneous state):
- Any YBCO-metal interface: always DG (contrast can reach 10⁸ below Jc)
- Metal-metal interfaces with |log(ρ₁/ρ₂)| > 3 (factor of 1000): DG
- Low-contrast interfaces (copper-silver): CG (shared edges)

**Keep the classification fixed within a time step.**

---

## 5. Mapping to the Thin-Shell Architecture

### Current structure

```
Layer 0 (copper)     → bottom edges: hanging on volume
Layer 0-1 boundary   → SHARED edges (top of block 0 = bottom of block 1)
Layer 1 (silver)     → interior edges: free DOFs
Layer 1-2 boundary   → SHARED edges
...
Layer N (copper)     → top edges: hanging on volume
```

### Proposed structure

```
Layer 0 (copper)     → bottom edges: hanging on volume
Layer 0-1 boundary   → shared (copper-silver: low contrast)
Layer 1 (silver)     → interior edges: free DOFs
Layer 1-2 boundary   → DUPLICATED + face + Nitsche (silver-hastelloy: high contrast)
...
Layer k-k+1 boundary → DUPLICATED + face + Nitsche (hastelloy-YBCO: extreme contrast)
...
Layer N (copper)     → top edges: hanging on volume
```

### Implementation points

1. **ThinShellFactory** (`create_edges_on_layers`): At flagged boundaries, create TWO copies of the interface layer's edges. Currently `link_elements_with_edges()` assigns the same Layer::Edges to both blocks; the change is to assign separate copies. Then create faces connecting the PENTA6TS elements on each side, providing master/slave pairs for assembly.

2. **Element assembly** (`mt_maxwell_h.cpp`): The volume terms (M, K) in `h_ts_hts()` stay unchanged. A new Nitsche interface function loops over the DG faces and assembles the 12×12 coupling matrix using the face's master/slave elements.

3. **DOF manager**: Duplicated boundary edges become independent free DOFs. They are NOT hanging — the T-matrix framework only processes edges at the volume-shell interface (bottom of layer 0, top of layer N), unaffected by internal DG boundaries.

4. **Important limitation**: Only the tangential (in-plane) H_t components are resolved as thin-shell DOFs. The normal component H_n is obtained from the volume elements above and below. The Nitsche coupling therefore only stabilizes the in-plane field at layer interfaces.

---

## 6. Nonlinear Materials

YBCO's power-law resistivity ρ(J) = (e_c/j_c)|J/j_c|^(n-1) spans:
- Below Jc: ρ ~ 10⁻¹⁴ Ω·m → contrast with hastelloy ~ 10⁸
- At Jc: ρ ~ 10⁻⁶ Ω·m → contrast ~ 1
- Above Jc: ρ >> 10⁻⁶ → contrast inverted

**Recommended approach:** Always-DG with contrast-adapted penalty. Keep the DG structure fixed; update α at each Newton step using current ρ values. The penalty naturally adapts: large when contrast is low, small when contrast is high.

**First prototype:** Use lagged ρ (Picard-style) for all Nitsche terms. BELFEM's hybrid Picard/Newton strategy (Messe et al. 2023, §2.7) supports this. Full Newton linearization of the Nitsche terms can come later.

---

## 7. Stability Considerations

- **Inf-sup:** Selective Nitsche may affect the discrete inf-sup condition in the mixed h-φ formulation. Test per Dular et al. 2021. If the modified form fails, hierarchical enrichment (bubble functions) may be needed at DG interfaces.
- **Checkerboarding:** The Nitsche coupling introduces additional interface modes. The tight nonlinear tolerance (ε < 10⁻¹¹) already used in BELFEM should mitigate this.

---

## 8. Complementary Strategies

- **Lobatto sub-layering + Nitsche:** Complementary. Lobatto handles spatial resolution; Nitsche handles conditioning at material jumps. Together they replace both the manual sub-layering workaround and the accidental fix from the broken edge key.
- **HDG:** Overkill for this problem. Appropriate only if the entire thin-shell stack were treated as DG.
- **Domain decomposition:** Conceptually equivalent to FETI/BDDC with layer-aligned subdomains. Not directly exploitable with direct solvers, but relevant for future Krylov methods.

---

## 9. Implementation Roadmap

### Phase 1: Linear prototype with nonsymmetric Nitsche
1. Single interface: hastelloy-YBCO with constant ρ per layer
2. Nonsymmetric Nitsche (consistency + penalty, no adjoint consistency)
3. Create faces at DG layer boundaries for master/slave assembly
4. Validate against analytical through-thickness profile
5. Verify energy/symmetry of the assembled DG interface block numerically
6. Compare condition number: CG (shared) vs DG (Nitsche)
7. Calibrate η empirically (start with η = 10, sweep)

### Phase 2: Nonlinear extension
1. Lagged ρ (Picard-style) for Nitsche terms
2. Test with YBCO power-law on tapestack3d
3. Verify convergence of Newton iteration
4. Optionally upgrade to symmetric Nitsche (add adjoint consistency) if convergence order matters

### Phase 3: Production integration
1. Automatic contrast detection from material definitions
2. Integration with Lobatto sub-layering
3. Full Newton linearization (if needed)

### Phase 4: Nitsche at shell/volume interface (exploratory)
1. Replace hanging edges at thin-shell/volume boundary with Nitsche coupling
2. Eliminates T-matrix cascade logic — uniform DG+Nitsche for all shell interfaces
3. Requires H(curl) DG theory (Houston, Perugia & Schotzau 2004–2005) — harder than inter-layer 1D scalar DG
4. Motivation is code simplification, not conditioning (shell/volume is not a high-contrast interface)
5. Prerequisite: validated Phase 1 inter-layer Nitsche

---

## 10. Literature

### Directly applicable
- Ern, Stephansen & Zunino (2009) — Weighted interior penalty for high-contrast diffusion
- Hansbo & Hansbo (2002) — Nitsche for interface problems
- Houston, Perugia & Schotzau (2004-2005) — Interior penalty DG for Maxwell

### Stability analysis
- Dular et al. 2021 (paper0) — Inf-sup stability for h-φ in HTS
- Boffi et al. 2013 — Saddle-point analysis for mixed EM
- Arnold 2018 — De Rham complex, structure preservation

### Implementation patterns
- Brenner & Scott 2008 Ch. 10 — Interior penalty DG (directly applicable: the through-thickness Nitsche reduces to 1D scalar DG, not full H(curl) DG)
- Dawson & Proft (2002) — CG-DG hybrid methods
- Houston, Perugia & Schotzau (2004-2005) — H(curl) DG theory (provides additional reassurance but not strictly required)

### BELFEM-specific
- Messe et al. 2023 (paper1) — Architecture, static condensation, nonlinear strategy
- Arsenault et al. 2023 (paper3) — h-φ coupling, interface conditions
- Alves et al. 2022a/b (paper5/6) — Thin-shell interface conditions

**Note:** Specific section/equation numbers from AI reviews should be verified against the actual papers. No existing HTS paper uses Nitsche at thin-shell layer interfaces — this would be a novel application.

---

## 11. Related Documents

- `todo/thinshell_automatic_sublayering.md` — Lobatto/golden-ratio sub-layering
- `todo/convergence_investigation_resolution_2026-03-14.md` — investigation resolution
- `literature/doc/notation_guide.md` — FEM notation conventions
- `tmp/claude.md`, `tmp/codex.md`, `tmp/grok.md` — reviewer assessments
