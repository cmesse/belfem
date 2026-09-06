# Thin-Shell H–φ Formulation: Design Note

**Date:** 2026-04-07
**Author:** Claude (this trace), with corrections from two Codex reviews (same day, see `devlog/dl20260407_thinshell_hphi_formulation_review.md` and `devlog/dl20260407_thinshell_three_way_coupling_review.md`) and empirical confirmation from Junie (same day)
**Purpose:** Capture the motivation, theory, and implementation sketch for adding a per-layer H–φ formulation switch to BELFEM's thin-shell stack, so that insulator and high-contact-resistance interlayers can be modelled without the conditioning catastrophes that the current pure-H thin shell exhibits.
**Status:** **Proposal — not yet implemented.** First draft reviewed by Codex on 2026-04-07; **§4 implementation sketch revised** in response. Junie's audit confirmed the diagnosis empirically (see §1.1 below). Second Codex review (same day) flagged the contact-impedance case as a *distinct formulation milestone*, not a reinterpretation of `h_ghost()` — this is now §3.4 and §5.2 / §9.A.3. **Sections marked ⚠️ have known gaps and need re-derivation before implementation starts.**

## 1.1 Independent confirmation (Junie audit, 2026-04-07)

A second independent audit by Junie confirmed the diagnosis with concrete measurements:

- **Failure signature:** baseline pure-H simulation with the buffer layer enabled produced a residual of **+1544 dB** (i.e., `~ 10^77`, effectively unbounded — a numerical catastrophe consistent with NaN propagation through MUMPS factorization).
- **Stable bound:** when the resistivity used in the bulk K assembly was numerically bounded (the "Path 2" stopgap from §10), the residual converged stably to **−156 dB**, restoring the working envelope.
- **Mesh wiring is not the bug:** the audit verified that the air–shell coupling for reduced-edge PENTA6TS elements is correctly established (1469 hanging DOFs), ruling out any vertex-refactor regression in the hanging-edge wiring path. This is reassuring in two ways: (a) the trace work in `devlog/dl20260407_layered_patch_test_trace.md` was not chasing a wiring bug, the wiring is fine; (b) the H–φ refactor in this document doesn't need to revisit the hanging-edge plumbing for the *external* air–conductor case.
- **Codex's mathematical review of the H–φ direction was confirmed valid** by Junie's separate read. So three independent audits (Claude / Codex / Junie) converge on the same recommendation: H–φ is the right long-term answer for true insulators, and a numerical bulk-K bound is the right stopgap.

The 1544 dB → −156 dB delta is a 1700 dB swing, which is the most direct empirical evidence we have that the bulk-K conditioning (not the Nitsche penalty, not the hanging-edge wiring) is what kills the pure-H formulation when an insulator is in the stack.
**Cross-refs:**
- `src/fem/maxwell/doc/thinshell_postprocessor_node_sharing.md` (related but different layer-interface issue)
- `devlog/dl20260407_layered_patch_test_trace.md` (the trace that surfaced this)
- `literature/papers/fem/schnaubelt2023.md` (the formulation we want to follow)
- `literature/papers/fem/alves2022b.md` (the only multilayer thin-shell paper, and the limit of what pure H can do)

---

## 1. The problem

BELFEM's current thin-shell formulation uses Nédélec edge `H` DOFs in *every* layer of the stack. For a typical REBCO tape this is fine when all layers are good conductors (Cu, Ag, Hastelloy, YBCO) — the contrast in `ρ` is moderate (`~ 4e-9` to `~ 1.2e-6` Ω·m, four orders of magnitude) and MUMPS handles it.

It **breaks** when the layer stack contains:

1. **A true insulator interlayer** like the buffer layer (epitaxial oxide, e.g. `CeO₂` / `Y₂O₃` / `MgO` / `YSZ`) at `ρ ≈ 1 Ω·m` and `h = 0.15 μm`.
2. **A high-resistance contact layer** like the interface resistance between Ag and YBCO (`rint`) at `ρ ≈ 1e-3 Ω·m` and `h = 0.1 μm`.

For these layers the per-layer "stiffness scale" `k = ρ / h` becomes:

| Layer | ρ [Ω·m] | h [m] | `k = ρ/h` [Ω] |
|---|---|---|---|
| Cu @ 77 K | 4e-9 | 2.5e-6 | 1.6e-3 |
| Ag @ 77 K | 3.4e-9 | 1.5e-6 | 2.3e-3 |
| Hastelloy @ 77 K | 1.23e-6 | 5e-5 | 2.5e-2 |
| YBCO (deep SC) | ~1e-12 | 1e-6 | ~1e-6 |
| **Rint** | 1e-3 | 1e-7 | **1e4** |
| **Buffer** | 1.0 | 1.5e-7 | **6.7e6** |

The bulk PENTA stiffness term is

```
K_bulk_entry  ~  ρ · (1/h²) · (h · area)  =  k_layer · area
```

so the buffer's bulk K entries are about **10⁹ times** the copper's. The resulting matrix has entries spanning 9 orders of magnitude *before* the Nitsche penalty contributes anything, MUMPS hits its precision limit, and the factorization produces NaN. The rint case is "only" 10⁷ apart and survives the factorization, but the Picard iteration stalls in a wrong attractor (`residual ~ -40 dB`).

Capping the resistivity used in the bulk K assembly is a numerical band-aid that lets the simulation **run** but doesn't represent the right physics: a true insulator should not have a diffusion term at all, it should *enforce* `∇×H = 0` as a constraint. Under-capping it (treating an insulator as a "moderate metal") changes the physics in a way that's quietly wrong.

## 2. What the literature does

We surveyed every thin-shell paper in `literature/papers/fem/`. The pattern is consistent and unambiguous:

| Paper | Treatment of insulator interlayer in thin shell |
|---|---|
| **Alves 2022a** | Single SC layer only. Substrate and coating layers are *omitted*: "this assumption is usually valid (and frequently used) since the substrate and the coating layers of HTS tapes have resistivities several orders of magnitude higher than the superconducting material." |
| **Alves 2022b §5.2** | Multilayer thin shell (substrate + HTS + silver). Substrate = `σ = 10 kS/m → ρ = 1e-4 Ω·m`. Silver = `1.6e-8 Ω·m`. **No insulator interlayer; max ρ in any layer is 1e-4.** |
| **Alves 2024** | "For simplicity, only the superconducting layer of HTS tapes is considered in this article." Single SC layer. |
| **Messe 2023 §5.1** | Single SC layer (1 μm), 4 mm wide, no other layers. |
| **Riva 2023** | Stack of REBCO tapes; each tape modelled as a single SC layer; no internal layer detail. |
| **Schnaubelt 2023** | Uses **H–φ globally**. Conducting region: edge-`H` DOFs. Insulating region (the *air around* the coil, in their case): scalar-`φ` DOFs. Static condensation at the interface. Quote: "This choice ensures that no spurious current appears in insulating domains **without the need for artificially large resistivities as in a pure H formulation** solved for the magnetic field strength everywhere." |

The pattern is:

1. **No paper in this set models an insulator interlayer inside a thin shell with a pure H formulation.**
2. Schnaubelt 2023 explicitly **names the antipattern** ("artificially large resistivities in a pure H formulation") and avoids it by switching to H–φ in the insulator region.
3. The H–φ approach has been validated for the *external* air region (Schnaubelt 2023, Alves 2022b for the surrounding air, BELFEM's own implementation), but not for an *internal* insulator interlayer of a thin-shell stack — that part is genuinely outside the published envelope.

For BELFEM specifically, the H–φ machinery already exists and works for the external air–conductor boundary of a tape. The proposal in this note is to **extend the same machinery to the internal layer-to-layer boundaries inside the thin-shell stack**, so that an insulator interlayer can be flagged as a "φ-block" and statically condensed against its neighboring conductor blocks.

## 3. Mathematical formulation

### 3.1 Per-layer DOF kinds

Each thin-shell layer block becomes one of two kinds:

| Kind | DOFs | Bilinear form (per layer) |
|---|---|---|
| **`H` layer** | edge_h on top + bottom + (mid for higher order) edges | `∫ μ₀ ∂_t H · v dV  +  ∫ ρ (∇×H)·(∇×v) dV` |
| **`φ` layer** | nodal phi on top + bottom + (mid) nodes | `∫ ∇φ · ∇v dV` (Laplace, no time dependence) |

A layer's kind is determined by a flag set when the `Layer` struct is constructed in `ThinShellFactory`. The simplest selection rule:

- If the material has a non-superconducting `ρ` above a threshold `ρ_φ_threshold` (e.g. `1e-2 Ω·m`), tag as `φ`.
- Otherwise tag as `H`.

A more flexible rule allows the user to specify per-material in the input file, e.g.

```
materials
{
    buffer
    {
        usermat { ... }
        formulation : phi ;     // ← new per-material override
    }
}
```

### 3.2 ⚠️ Layer-to-layer interface coupling — REVISED after user decision (2026-04-07)

**Strategic decision (2026-04-07, after Codex and Junie reviews):** the design originally framed layer-interface coupling as an *extension of the Nitsche-DG-Ghost path*. The user pushed back: **static condensation is the more BELFEM-friendly default**, both because the infrastructure already exists (`hang_thinshell_edges_on_edges_*` for H–H, `hang_thinshell_edges_on_nodes_*` for H–φ, `DofData::create_dofwise_t_matrices_master` for the T-matrix construction) and because the physics doesn't actually need the weak-coupling features Nitsche was paying for.

The revised position is: **static condensation everywhere as the default**, with **Nitsche as an opt-in fallback** controlled by an input-file switch (see §9 item 5 for the input syntax).

The three layer-interface combinations under the new design:

1. **H–H** (current behaviour, simplified): conductor adjacent to conductor.
   - **New default:** static condensation. The "slave" side's edge DOFs are expressed as `±1` linear combinations of the "master" side's edge DOFs via T-matrix, exactly the same machinery as `hang_thinshell_edges_on_edges_*` already uses for the external air–conductor boundary. The slave side's bulk PENTA stiffness still gets assembled, but the T-matrix projection routes its contributions onto the master side's DOFs.
   - **Opt-in fallback:** Nitsche-DG-Ghost via the existing `h_ghost` kernel + ghost facet sideset machinery. This is what BELFEM does today and is preserved as-is. Selectable per-tape or globally via the input-file switch in §9 item 5.
   - Result: the regularized harmonic-mean penalty, the SIPG fix, the `Dm`/`Ds` operators, all the work we did this week becomes a *non-default but supported* code path, not the main path.

2. **H–φ** (new): conductor adjacent to insulator.
   - Static condensation only. **No Nitsche.** The φ side has nodal DOFs, the H side has edge DOFs, and the existing `hang_thinshell_edges_on_nodes_*` machinery already knows how to express a shell edge as `phi(node 0) - phi(node 1)` via the LINE2 T-matrix coefficients `[1, -1]` in `DofData::create_dofwise_t_matrices_master`. The plumbing change is *only* the dispatch — extending `hang_thinshell_edges_on_nodes_*` (or a sibling function with the same body) to fire on internal layer interfaces, not just external mesh facets.
   - **Codex review correction (preserved):** the existing `hang_thinshell_edges_on_*` functions hardcode the volume-tet master/slave roles via `aFacet->master()` / `aFacet->slave()`. For internal layer interfaces, both sides are PENTA layer elements, and the master/slave roles are the *layer indices* in the stack, not volume tet pointers. **The mathematical form is the same as the external case** (verified by inspection of `cl_FEM_DofMgr_DofData.cpp:3719-3727` — the `[1, -1]` LINE2 coefficient gives `edge_DOF = phi_0 - phi_1` regardless of whether the source is a volume tet node or a layer φ node), but the dispatch and facet wiring are new.

3. **φ–φ** (new): insulator adjacent to insulator (rare but possible).
   - No coupling needed at all. φ DOFs are nodal and the layer-interface nodes are shared between adjacent layer blocks (today's behaviour, since `hasDuplicates` is only true for material changes and would be false at φ–φ). Continuity of φ is automatic.

**Why this is a big simplification compared to the first draft:**

| First draft (Nitsche-extended) | Revised (condensation-default) |
|---|---|
| New `phi_ts_insulator` Nitsche kernel | Reuse existing `phi`/`phi_tet4` form, just on a PENTA6TS prism |
| New "internal H–φ coupling sideset" + new IWG branch | Extend existing `hang_thinshell_edges_on_*` dispatch |
| New variational derivation for the H–φ Nitsche penalty | The math is already in `DofData::create_dofwise_t_matrices_master` |
| H–H still uses Nitsche, with all the harmonic-mean / SIPG bookkeeping | H–H uses condensation by default; Nitsche becomes a debug fallback |
| Risk: bugs in the new H–φ Nitsche kernel | Risk: dispatch / facet wiring edge cases (more localized) |

**Why keep Nitsche at all:**

1. **Debugging.** A flip-the-switch alternative implementation is the strongest validation tool we have. The 20 dB patch test failure earlier this week was caught by exactly this kind of cross-check.
2. **Comparison testing.** Schnaubelt 2023, Alves 2022b et al. use condensation. Having Nitsche around lets BELFEM be the platform that can run *both* and report the difference.
3. **Future research.** Non-conforming meshes, intentionally weak coupling at degraded interfaces, etc. — Nitsche is the right tool for those if we ever want them.

The Nitsche infrastructure stays in tree as a maintained but non-default code path.

### 3.3 Cohomology cuts

A multiply-connected `φ` region (which is always the case when the layer wraps around, e.g. in a coil) needs cuts to make `φ` single-valued on closed loops. BELFEM already does this for the external air, and `CutFactory` already creates the cuts. The question is whether it also creates them inside an internal `φ` layer.

For a **planar buffer interlayer** in a flat tape, the buffer is *not* multiply connected (it's a flat slab with no closed loops in its interior), so no cut is needed. For a **wound coil** with a buffer interlayer, the buffer wraps around with the rest of the tape and inherits the same topology — cuts may be needed. This will need to be checked on a case-by-case basis when wiring it into `CutFactory`. For the initial implementation, **assume no cuts inside the buffer**, and add a runtime warning if a multiply-connected φ-block is detected without cuts.

### 3.4 Three-way coupling categorization (added 2026-04-07 after second Codex review)

The first revision of this note framed the layer-interface problem as a **two-way** split: continuity (H–H or H–φ) handled by condensation, true insulators handled by the φ formulation. After thinking about the rint contact-resistance case, the obvious next move seemed to be "use Nitsche for that, with a different α."

**Codex's second review caught the conflation.** BELFEM's existing `h_ghost()` is a **continuity-enforcing** Nitsche kernel — verified at `mt_maxwell_h.cpp:1838` and `:1901`. It weakly enforces

```
⟦H_t⟧ = 0      (tangential H continuous across the interface)
```

via a penalty term `α·⟦H_t⟧·⟦H_t⟧` plus the SIPG consistency / adjoint-consistency terms. Reinterpreting α does **not** turn this into a contact-impedance law. A true contact-impedance law is

```
n × ⟦E⟧ = ρ_sheet · J_n      (Robin / impedance jump, ρ_sheet = ρ · h)
```

which is a **different bilinear form**, with its own consistency / adjoint terms, derived from the asymptotic limit of a thin resistive layer (the T2TCL surface contribution in Schnaubelt 2023 §III). The two operators happen to live on the same kind of facet, but they enforce different physics. You cannot reach one from the other by tuning α.

The corrected three-way categorization is therefore:

| Interface physics | DOFs in the layer | Coupling mechanism | Status in this note |
|---|---|---|---|
| **Continuity** between two conductors (H–H), or between a conductor and an insulator (H–φ) | edge_h on conductor side; nodal φ on insulator side | **Static condensation** via existing `hang_thinshell_edges_on_*` machinery (default), with the existing `h_ghost()` Nitsche kernel preserved as an opt-in legacy fallback for H–H. | §4.2, fully scoped. Preconditions in §9.A.1 (scalar-φ operator) and §9.A.2 (internal H–φ topology). |
| **True insulator interlayer** (e.g. epitaxial buffer, ρ ≥ 1e-2 Ω·m) | nodal φ everywhere in the layer | The φ formulation makes the layer's bulk operator the elliptic Laplace `μ₀ ∫ ∇φ · ∇v dV`. Continuity at its boundaries is enforced by the H–φ static condensation above. | §3.1–§3.3, §3.5, §4.5, scoped. |
| **Contact-impedance interlayer** (e.g. rint, ρ ≈ 1e-3 Ω·m) | *None* — the layer is collapsed to a zero-thickness sideset between its neighbors | A **new dedicated kernel**, derived from the T2TCL surface contribution in Schnaubelt 2023 §III, that enforces `n × ⟦E⟧ = ρ_sheet · J_n` with `ρ_sheet = ρ_rint · h_rint`. **This is not a knob on `h_ghost()`. It is a new bilinear form.** | **NOT scoped in this note.** Promoted to a separate formulation milestone, gated on §9.A.3 below. |

**Why this matters now:** the rint case is the one that motivated the discussion in the first place ("can we still use Nitsche for the contact resistance?"), so the temptation was to claim it was already covered by the Nitsche fallback in §4.2. It is not. The Nitsche fallback in §4.2 is the *existing continuity kernel*, kept around for debugging the H–H continuity path. Treating rint correctly needs a separate derivation, a separate kernel, and ideally a separate test case before any code is written.

**Recommended path for rint while the contact-impedance formulation is being designed:** the **honest cap** (Path 2 in §10, Option A in §5.2). The cap is a numerical surrogate, but it is mathematically self-consistent (the operator that's solved and the operator that's reported agree) and the empirical evidence (Junie's 1700 dB swing from +1544 dB → −156 dB) shows it gives a working envelope for everything below the buffer. The cap is the right stopgap precisely because it does *not* pretend to be a faithful contact-impedance model.

### 3.5 Postprocessor

Currently the postprocessor recovers H from `edge_h` in conductor blocks and B from `phi` in air blocks. It would need:

- A new `MaxwellPostprocessorType::ThinShellInsulator` (or extend `ThinShellConductor`) for the φ-layer case.
- The recovery for the φ-layer evaluates `H = -∇φ`, `B = μ₀ H`, `J = 0` (no current in an insulator).
- Same caveat as the existing `J/Jc` leakage at the layer interfaces: the φ-layer's nodes are shared with the neighboring conductor's nodes, so the recovered field at those nodes will be the average. This is an *unrelated* problem documented in `src/fem/maxwell/doc/thinshell_postprocessor_node_sharing.md`; the H–φ change here doesn't fix it but doesn't make it worse either.

## 4. Implementation sketch

### 4.1 `ThinShellFactory` (`src/mesh/cl_ThinShellFactory.cpp`)

Today the `Layer` struct in `cl_ThinShellFactory.hpp:35-44` carries:

```cpp
struct Layer
{
    bool hasDuplicates = false ;
    Cell< Node *  > Nodes ;
    Cell< Edge *  > Edges ;
    Cell< Edge *  > EdgeDuplicates ;
    Cell< Face *  > Faces ;
    Cell< Face *  > FaceDuplicates ;
    Cell< Facet * > GhostFacets ;
};
```

Extend to:

```cpp
enum class LayerFormulation { H, Phi };

struct Layer
{
    bool hasDuplicates = false ;
    LayerFormulation formulation = LayerFormulation::H ;   // NEW
    Cell< Node *  > Nodes ;
    Cell< Edge *  > Edges ;
    Cell< Edge *  > EdgeDuplicates ;
    // ... unchanged ...
};
```

The `hasDuplicates` flag and `EdgeDuplicates`/`FaceDuplicates` containers are only meaningful for **H–H interfaces**. For **H–φ interfaces** the natural decoupling is "the φ side has nodes only, the H side has edges only", so duplicate edges aren't needed — but the *interface* between the two has to be set up via the static condensation path, not the ghost-facet path.

Selection of `formulation` happens in `create()` after the materials list is known. Roughly:

```cpp
for ( uint l=0; l<tNumLayers; ++l )
{
    Layer * tLayer = tLayers( l );
    real rho_l = material_resistivity_at_77K( aMaterials( l / tOrder ) );
    if ( rho_l > rho_phi_threshold )
        tLayer->formulation = LayerFormulation::Phi ;
}
```

Then `hasDuplicates` is recomputed: it's true only at H–H interfaces between layers of *different* materials. H–φ and φ–H interfaces don't get edge duplicates; they get a static-condensation hookup instead. φ–φ interfaces don't get anything special.

The `create_ghost_facets()` machinery only fires on `hasDuplicates == true`, so it automatically *doesn't* create ghost facets at H–φ or φ–φ interfaces. Good — the Nitsche penalty doesn't apply there.

### 4.2 ⚠️ `MaxwellFactory::hang_thinshell_edges_on_*` — REVISED after user decision (2026-04-07)

**The user's decision to default to static condensation simplifies this section substantially.** The first revision (post-Codex) called for new ghost-facet sidesets, a new IWG branch, and a new variational derivation. With condensation as the default, all of that goes away and the work reduces to **extending the existing `hang_thinshell_edges_on_*` dispatch to fire on internal layer interfaces**, not just external mesh facets.

#### What stays the same

- **The existing four functions** (`hang_thinshell_edges_on_nodes_bottom`, `_top`, `_on_edges_bottom`, `_top`) at `cl_MaxwellFactory.cpp:1364-1680` are *not modified*. They handle the external air–tape boundary correctly today and the corresponding test cases pass.
- **The T-matrix construction** in `DofData::create_dofwise_t_matrices_master` (`cl_FEM_DofMgr_DofData.cpp:3563-3760`) is *not modified*. The LINE2 case at line 3719-3727 already produces the right `[1, -1]` coefficients regardless of whether the source is a volume tet node or a layer φ node.
- **The `hang_thinshell_edges_on_edges_*` orientation matching by node index** (lines 1545-1559 and 1646-1660) works for any pair of edges that share two nodes. It doesn't care whether those nodes belong to a volume tet or to a layer PENTA.

#### What changes

The new work is a **second pass through the layer stack** in `MaxwellFactory::create_thinshells()`, after the existing external-boundary dispatch (lines 1220-1303), that walks each *internal* layer-to-layer interface and fires the appropriate hang function:

```cpp
// === existing pass: external air-tape boundaries ===
// (lines 1220-1303 today: bottom of stack hangs on volume master, top hangs on volume slave)
// → unchanged

// === NEW pass: internal layer-layer interfaces ===
for ( each pair of adjacent layer blocks (b, b+1) inside the stack )
{
    LayerFormulation fmt_b   = layers(b)->formulation;
    LayerFormulation fmt_bp1 = layers(b+1)->formulation;

    if ( fmt_b == H && fmt_bp1 == H )
    {
        // H-H internal interface.
        // Default: condensation. Use a new helper analogous to
        // hang_thinshell_edges_on_edges_top, but with the master being
        // layer b's TOP face edges and the slave being layer b+1's
        // BOTTOM face edges (both PENTA layers, both inside the stack).
        if ( use_nitsche_for_internal_hh )
        {
            // legacy / opt-in: do nothing here, the existing ghost-facet
            // path in ThinShellFactory::create_ghost_facets will handle
            // this interface via h_ghost. (See §9 item 5 for the input
            // file switch.)
        }
        else
        {
            hang_thinshell_edges_on_edges_internal( layers(b), layers(b+1) );
        }
    }
    else if ( fmt_b == H && fmt_bp1 == Phi )
    {
        // H-phi internal interface.
        // Layer b's TOP face edges hang on layer b+1's BOTTOM face nodes
        // via the LINE2 [1, -1] T-matrix coefficient — exactly the same
        // construction as the existing _on_nodes_top function but with
        // the slave being a PENTA phi layer instead of a volume tet.
        hang_thinshell_edges_on_nodes_internal( layers(b), layers(b+1), /*upper-is-phi=*/true );
    }
    else if ( fmt_b == Phi && fmt_bp1 == H )
    {
        // phi-H internal interface (mirror image).
        hang_thinshell_edges_on_nodes_internal( layers(b), layers(b+1), /*upper-is-phi=*/false );
    }
    else  // phi-phi
    {
        // No coupling needed. phi DOFs are nodal and the shared layer-
        // interface nodes are not duplicated for phi-phi pairs.
    }
}
```

The two new helper functions (`hang_thinshell_edges_on_edges_internal` and `hang_thinshell_edges_on_nodes_internal`) are the only genuinely new code. They are *near-copies* of the existing `_top` variants with three changes:

1. **Master and slave are both PENTA layer elements** (not a PENTA + a volume tet). The "master" is layer `b`'s top face (PENTA face index 1); the "slave" is layer `b+1`'s bottom face (PENTA face index 0).
2. **No `aFacet->master()` / `aFacet->slave()` calls.** The layer pair is passed in directly, so the master/slave PENTA elements are known by reference rather than fetched via the mesh facet.
3. **The shared "facet" between the two layer PENTAs is the layer-interface plane**, not a separate `Facet` object. The matching of node indices uses the temp-index-by-facet-position trick that the existing `_top` function already uses (verified in §3.1 of `devlog/dl20260407_layered_patch_test_trace.md`), so the geometry handling is already in place.

#### Estimated effort

- `hang_thinshell_edges_on_edges_internal`: **half a day**, mostly reading the existing `_top` function and copying it with the master/slave fetch removed.
- `hang_thinshell_edges_on_nodes_internal`: **half a day**, same shape.
- Dispatch in `MaxwellFactory::create_thinshells()`: **1-2 hours**, the new pass shown above.
- Input file switch parsing for the Nitsche fallback: **1 hour**.
- Tests: a regression test (Cu/Ag/Hastelloy/YBCO with no buffer, condensation default vs Nitsche fallback should give the same result to round-off) plus a new test with the buffer included.

**Total revised estimate: 1.5-2 days for §4.2** (compared to "medium-large, multi-day" in the first revision). This is the part of the design that benefited most from the user's decision.

### 4.3 `IWG_Maxwell::create_custom_vectors_and_matrices` (`src/fem/maxwell/cl_IWG_Maxwell.cpp:579-598`)

For φ-layer blocks, allocate scalar `phi` DOFs instead of `edge_h`/`face_h`. The `Group::dof_type` already supports per-block DOF assignments (this is what air vs conductor dispatch already uses); the φ-layer block would just be assigned the air-like DOF set.

### 4.4 `IWG_Maxwell::link_to_group` dispatch (`src/fem/maxwell/cl_IWG_Maxwell.cpp:430-520`)

Today the `DomainType::ThinShell` branch dispatches by material type:

```cpp
case( DomainType::ThinShell ):
{
    Material * tMat = aGroup->material();
    switch ( tMat->type() ) {
        case PureMetal: mFunMKF = &maxwell::h_ts_metal ;
        case LookupAlloy: mFunMKF = &maxwell::h_alloy ;
        case HTS: mFunMKF = &maxwell::h_ts_hts ;
        default: mFunMKF = &maxwell::h_ts ;
    }
}
```

Add a new branch *before* the material-type switch:

```cpp
case( DomainType::ThinShell ):
{
    if ( aGroup->layer_formulation() == LayerFormulation::Phi )
    {
        // Insulator layer: scalar phi formulation, no diffusion term.
        mFunMKF = &maxwell::phi_ts_insulator ;   // NEW
    }
    else
    {
        // Existing material-type switch for H-formulation layers.
        Material * tMat = aGroup->material();
        // ... existing switch ...
    }
}
```

`Group::layer_formulation()` is a new accessor that reads back the `formulation` flag set by `ThinShellFactory`.

### 4.5 ⚠️ New matrix kernel `phi_ts_insulator` — REVISED after Codex review

**Codex review correction:** the first draft of this section was wrong about which `aMatrices` slot the operator goes into. BELFEM's existing scalar-φ kernels assemble `μ₀ · trans(B) · B` into **`aMatrices->M()`**, not `aMatrices->K()`. Verified against `mt_maxwell_phi.cpp`:

- `phi()` at `cl_FEM_Calculator.cpp:24-44`: assembles into `M()`, then multiplies by `μ₀`.
- `phi_tri3()` at `:101-108`: assembles into `M()` with `μ₀ * dV / 2`.
- `phi_tet4()` at `:113-120`: assembles into `M()` with `μ₀ * dV / 6`.

This is BELFEM's convention for the magnetic scalar potential: the static elliptic operator `μ₀ ∫ ∇φ · ∇v dV` is stored in the *mass* slot rather than the *stiffness* slot. Why? Because for a static (time-independent) air block, BELFEM's time-stepping assembly adds `M + Δt·K` (or similar) into the global Jacobian, and putting the elliptic operator in `M` makes it survive the assembly without being scaled by `Δt`. The naming is "M" but the role is "elliptic stiffness for a static-in-time DOF".

The corrected kernel is:

```cpp
void
phi_ts_insulator( Calculator * aCalc, TimestepMatrices * aMatrices )
{
    const Vector< real > & w = aCalc->integration()->weights();

    for ( uint k=0; k<aCalc->num_intpoints(); ++k )
    {
        // gradient operator for the prism's nodal scalar basis
        const Matrix< real > & B = aCalc->B( k );

        // Laplace operator: -mu0 * div(grad phi) = 0.
        // BELFEM convention: assemble this into M(), not K() — that
        // matches phi(), phi_tri3(), phi_tet4() in mt_maxwell_phi.cpp.
        aMatrices->M() += w(k) * trans(B) * B * aCalc->dV(k);
    }

    aMatrices->M() *= constant::mu0;
    // K stays at zero. There is no diffusive (resistive) term in an insulator.
}
```

**Open question (carried to §9):** does the prism gradient operator `aCalc->B(k)` already exist for a PENTA6TS used as a *scalar* element? Today the PENTA6TS is only ever used as an `edge_h` (Nédélec) element. If we want it to also carry nodal scalar DOFs, the `Calculator` machinery may need a new dispatch path for the scalar gradient on a thin-shell prism. Estimated effort: small if the underlying `IF_PENTA6` already exposes `dN/dXi`, large if not.

### 4.6 ⚠️ `cl_Maxwell_FieldList` block-DOF assignment — REVISED after Codex review

**Codex review correction:** the first draft called this a "small field-list tweak". It is **not small**. Today the `ThinShell` mapping is hard-wired in **two** different switch tables that both need to learn about the new `ThinShellInsulator` kind:

- `cl_Maxwell_FieldList.cpp:255` — block-DOF dispatch: `case DomainType::ThinShell:` is a `case` *fall-through* with `Conductor`, both routed to `create_doftable( Conductor, ... )`. A new `ThinShellInsulator` case would route to `create_doftable( Air, ... )` instead, but care is needed because the existing case is a fall-through and the new case mustn't break the conductor path.
- `cl_Maxwell_FieldList.cpp:387` — sideset-DOF dispatch: `case DomainType::ThinShell:` routes to `create_doftable( ThinShell, ... )`. This mapping is for thin-shell *sidesets*, not blocks; it's used by the cohomology cut and ghost-facet paths. A new `ThinShellInsulator` sideset would need its own dedicated DOF list (likely `Air`-like, since all DOFs on the φ-layer are nodal).

In addition, `cl_Maxwell_FieldList.hpp` defines the `Conductor`, `Air`, `Ferro`, `Coil`, `ThinShell`, `Cut`, `BoundaryAir`, `Background`, `Symmetry` DOF lists as private members. A `ThinShellInsulator` list would need to be added there too.

**Other places that depend on `DomainType::ThinShell` and need to be updated for the new enum value:**

- `mesh::Block::set_domain_type` enum machinery (likely in `cl_Block.hpp` or similar)
- `MaxwellFactory::set_block_types_in_magnetic_equation()` switch list (around `cl_MaxwellFactory.cpp:1717`)
- The kernel dispatch in `IWG_Maxwell::link_to_group()` (the existing `case DomainType::ThinShell:` branch around `cl_IWG_Maxwell.cpp:430`)
- The postprocessor dispatch in `MaxwellPostprocessor::select_blocks_and_materials()` (`cl_MaxwellPostprocessor.cpp:152-264`)
- The `to_string` overload for `DomainType` (used in error messages and the new block label format we added on 2026-04-06)

**Estimated effort:** medium. ~10 file touches, all small but coordinated. Worth a dedicated commit *separate* from the kernel-level changes so the diff is reviewable.

**Backward compatibility:** existing input files don't use `ThinShellInsulator` and never have. The default selection rule (only flag a layer as `ThinShellInsulator` if its `ρ > ρ_φ_threshold`) leaves all current cases unchanged.

### 4.7 The Picard / Newton iteration

No change. The φ-layer's K is linear and constant; the Picard/Newton non-linear iteration only acts on conductor blocks. The φ-layer just sits in the global matrix unchanged from one iteration to the next.

## 5. Risks and unknowns

### 5.1 Backward compatibility

Existing input files that don't have an insulator interlayer should produce **bit-identical results** before and after this change. The selection rule (default `LayerFormulation::H` for any material with `ρ < ρ_φ_threshold`) needs to leave Cu/Ag/Hastelloy/YBCO untouched. Easy to verify with a regression test on the existing tape-stack test cases.

### 5.2 ⚠️ The "rint" case — REVISED after Codex review

The contact resistance layer (Rint, `ρ ≈ 1e-3`) is **not** an insulator — it's a bad conductor that should still carry some current. Putting Rint in the φ branch would lose the current-transfer physics entirely, which is wrong: the whole point of Rint is to model a current barrier.

So the threshold `ρ_φ_threshold` has to be set **above** Rint's resistivity, i.e. `> 1e-3`. Pick `1e-2` to be safe. The buffer at `1.0` is then comfortably above the threshold; rint at `1e-3` is comfortably below.

**Codex review correction:** the first draft proposed capping the bulk `ρ` for Rint while letting `mat->rho()` and the postprocessor see the *uncapped* `ρ`. **This is internally inconsistent:** the operator that's actually solved uses one ρ, the loss/J/Jc/postprocessor uses another, so the postprocessor's reported "physics" no longer corresponds to the solved state. That isn't a "small loss of physics fidelity for a numerical stability win" — it's a *quietly wrong* output.

The corrected position is one of these:

**Option A — honest cap:** if we cap `ρ` for the bulk K, we should *also* cap it for `mat->rho()`, the postprocessor, and any downstream physics. The simulation then uses `ρ_eff = min(ρ, ρ_cap)` consistently everywhere. The user is informed by an explicit `gLog` message of the form: "Layer X resistivity ρ=1e-3 capped to ρ_eff=1e-4 for numerical stability. Reported losses, J/Jc, etc. correspond to the capped value, not the input value." This is honest about the surrogate.

**Option B — no cap, accept the conditioning hit:** keep the real ρ, accept that Rint at 1e-3 with `h=0.1 μm` gives a 10⁷ conditioning ratio, and rely on MUMPS to grind through it. The "stalls at -40 dB" we observed today suggests MUMPS *can* factor it but the Picard iteration doesn't converge to the right attractor — so this option needs a better outer solver (Newton with line search? or the adaptive timestep work in `nonlinear_iteration_strategy_near_quench.md`).

**Option C — also put rint in the φ branch:** acknowledge that Rint physically is "almost an insulator", model it as one for the H-field, and add a *separate* current-transfer model (e.g. an explicit interface current source) to capture the current sharing. This is cleaner physics than A, but the current-transfer term is its own derivation.

**Option D — contact-impedance interface (Schnaubelt 2023 §III T2TCL):** collapse Rint to a *zero-thickness* sideset between Ag and YBCO (no extruded prism, no edge DOFs in the layer) and introduce a **dedicated surface contribution** to the weak form that enforces

```
n × ⟦E⟧  =  ρ_sheet · J_n      with     ρ_sheet = ρ_rint · h_rint
```

This is a Robin / impedance jump law on the interface, not a continuity penalty. It is *not* a reuse of `h_ghost()` with a different α — `h_ghost()` enforces `⟦H_t⟧ = 0` (continuity), and tuning α in `h_ghost()` cannot produce a Robin law because the underlying bilinear form is wrong. **A new kernel is required**, derived from the T2TCL surface contribution in Schnaubelt 2023 §III. The mathematical advantages are:

- The matrix entry sees `ρ_sheet = ρ · h ≈ 1e-10 Ω·m²` rather than `ρ / h ≈ 1e4 Ω`. The catastrophic `1/h` conditioning of the bulk K disappears.
- The current-transfer physics is captured exactly (the interface absorbs J_n through the sheet impedance, exactly as in a real contact film).
- The asymptotic h → 0 limit is exact: the formulation *is* the limit of the bulk-layer model with `ρ · h` held fixed. This is the standard TSA derivation.

**The mathematical disadvantages and unknowns:**

- Through-thickness variation of H inside the layer is lost (0th order in h). For a 0.1 μm contact film with no in-plane current variation, this is not a meaningful loss; for a thicker partial-contact layer it would be.
- The exact form of the consistency / adjoint-consistency terms (the analogues of the SIPG terms in `h_ghost()`) needs to be derived from the variational formulation, not guessed.
- The interaction with the existing `hang_thinshell_edges_on_*` static condensation at the *neighboring* H–H continuity interfaces (Ag-top and YBCO-bottom) needs to be specified — does the contact-impedance kernel replace those interfaces, sit between them, or live on a separate sideset?

**Codex's second review (2026-04-07) is explicit on this point:** Option D is plausible and literature-backed, but it is *not* equivalent to BELFEM's current `h_ghost()` and must not be described as such. It is a **separate formulation milestone**, gated on the precondition in §9.A.3 below. See `devlog/dl20260407_thinshell_three_way_coupling_review.md` for the full Codex finding and the file:line references that establish that `h_ghost()` is currently a continuity-only kernel.

**Recommended near-term path:** **Option A (honest cap)** as a stopgap, with the surrogate-warning gLog message and the threshold value (e.g. `ρ_cap = 1e-4`) exposed through the input file. The cap is mathematically self-consistent (solve and postprocessor agree on `ρ_eff`) and Junie's empirical 1700 dB swing confirms it gives a working envelope. **Option D is the right long-term answer** but it is a separate milestone with its own preconditions and derivations, and it should not be conflated with the H–H continuity Nitsche fallback in §4.2.

The recommended layer-treatment table becomes:

| Layer ρ range | Near-term treatment | Long-term treatment | Caveat |
|---|---|---|---|
| `ρ < 1e-4 Ω·m` (real conductors) | H-formulation, no cap. | Same. | None. |
| `1e-4 ≤ ρ < 1e-2 Ω·m` (high-resistance contacts, e.g. Rint) | H-formulation with **honest cap**: `ρ_eff = ρ_cap` everywhere — solve, rho-call, postprocessor, losses. **`gLog` warning emitted.** Rint goes here for now. | **Option D** (zero-thickness contact-impedance interface, T2TCL kernel from Schnaubelt 2023 §III). Gated on precondition §9.A.3. | Near-term: loss/J/Jc/heating values correspond to `ρ_eff`, not the input ρ. The user has to know this. Long-term: Option D removes the surrogate. |
| `ρ ≥ 1e-2 Ω·m` (true insulators, e.g. buffer) | **φ-formulation**. No diffusion, no resistive losses; H field is constrained to ∇×H = 0 by the elliptic operator. | Same. | No current flow at all in the layer. Loss is exactly zero. |

The cap value `1e-4` matches Alves 2022b's substrate — the literature limit for what pure H can handle. Above that, the H formulation breaks down regardless of solver, and the right answer is H–φ for true insulators or the contact-impedance interface for thin resistive films.

### 5.3 Cohomology cuts inside a φ-layer

For a planar tape this is fine (no cuts needed inside the buffer). For a wound coil, the buffer wraps around with the rest of the tape and inherits the multi-connectedness. The proposal: **assume no cuts inside the buffer** for the initial implementation, and add a runtime check that detects multi-connected φ-blocks and warns. The full cut machinery can be wired in as a follow-up if a real coil case demands it.

### 5.4 The `hang_thinshell_edges_on_*` dispatch

The four existing functions assume the original mesh facet has a master/slave that is either fully air or fully conductor on each side. With internal H–φ layers, the *layer-to-layer* interfaces inside the stack now also become "H to φ" boundaries that need static condensation. The existing functions can probably be reused as-is, but the dispatch logic in `MaxwellFactory::create_thinshells()` needs new branches for internal interfaces.

This is the most involved part of the implementation. Worth tracing carefully and writing dedicated tests for each combination of layer kinds.

### 5.5 Postprocessor node sharing

The existing issue documented in `src/fem/maxwell/doc/thinshell_postprocessor_node_sharing.md` (J/Jc leakage at layer interfaces because nodes are shared between adjacent layer blocks) is **not made better or worse** by this change. The H–φ split adds new interfaces but uses the same node-sharing convention. The eventual fix for both issues is to add `Layer::NodeDuplicates` and route adjacent blocks to disjoint node sets at H–φ interfaces (and at H–H interfaces with `hasDuplicates == true`).

### 5.6 Thermal coupling

If the thermal kernel is active, the φ-layer needs a corresponding thermal block. The simplest approach is to leave the thermal block as a `Conductor`-type with the insulator's thermal properties (heat capacity, conductivity), and just skip any Joule loss term (since `J = 0` in the insulator). No deep integration is needed beyond what BELFEM already does for air blocks in coupled thermal runs.

## 6. Effort estimate

Rough breakdown for a single-developer implementation:

1. `ThinShellFactory` Layer struct + per-layer formulation flag: small
2. `IWG_Maxwell` block-type dispatch + new `phi_ts_insulator` kernel: small
3. `cl_Maxwell_FieldList` DOF allocation per layer kind: small
4. `MaxwellFactory::create_thinshells` internal-interface dispatch (the trickiest part): medium-large
5. New `phi_ts_insulator` kernel implementation + tests: small
6. Regression test on existing input files (Cu/Ag/Hastelloy/YBCO stacks): small
7. New test case with buffer included, comparing against H-only-with-buffer-omitted as reference: small
8. Documentation: small

The non-trivial work is concentrated in step 4 — the geometry of the internal layer-to-layer H–φ interface. Worth bench-marking against a `mesh-aware` debug print of the assembled K rows to confirm the static condensation is doing what we think it is.

## 7. Test plan

1. **Regression**: existing `tape1` Cu/Ag/Hastelloy/YBCO stack (no buffer). Same K, same residual, same J/Jc as before. Bit-for-bit if possible.
2. **Buffer-omitted**: copy of `tape1` with buffer line commented out. Should match the regression case.
3. **Buffer-included with H-only (current behaviour)**: should reproduce the NaN we see today.
4. **Buffer-included with H–φ (new behaviour)**: should run cleanly. Compare K-field above the tape against case (2). If the buffer is truly inert physically, the field should be within discretization error of case (2).
5. **Buffer-included with rint at 1e-3**: should converge with the bulk-K cap on rint and the φ-formulation on buffer. Compare against (2) — the rint should produce a small but visible difference in the current ramp-up time.
6. **Wound-coil case (future)**: confirm cohomology cuts are correctly created/ignored inside the φ-layer. May need new test infrastructure.

## 8. References

- **Schnaubelt, E. et al. (2023)** "Electromagnetic simulation of no-insulation coils using H–φ thin shell approximation." *IEEE Transactions on Applied Superconductivity*, **33** (5), 4900906. Section II describes the H–φ formulation choice and explicitly names the "artificially large resistivities in a pure H formulation" antipattern.
- **Alves, B. de Sousa et al. (2022a)** "3-D finite-element thin-shell model for high-temperature superconducting tapes." *IEEE Transactions on Applied Superconductivity*, **32** (6), 7500411. Section II.A justifies omitting non-superconducting layers.
- **Alves, B. de Sousa et al. (2022b)** "A thin-shell H–φ-formulation for superconducting devices." *Superconductor Science and Technology*, **35** (2), 024001. Section 5.2 demonstrates the only published multilayer thin shell — substrate + HTS + silver, max ρ = 1e-4 Ω·m.
- **Alves, B. de Sousa et al. (2024)** "2-D H–φ thin-shell model for HTS tapes implemented in COMSOL." *IEEE Transactions on Applied Superconductivity*, **34** (5). Single SC layer.
- **Messe, C. et al. (2023)** "BELFEM: ..." *Superconductor Science and Technology*, **36**, 114001 (paper1). Section 5.1 single-tape benchmark uses single SC layer.
- **Riva, T. et al. (2023)** "H–φ Formulation in Sparselizard combined with DDM for superconducting tapes, stacks, and twisted wires." Single SC layer per tape in a stack.

## 9. Open questions for review — REVISED after Codex review

The first draft of this section had 6 open questions. Codex's review made it clear that **the first two questions are actually preconditions** that have to be answered before the implementation sketch is even meaningful. Restructured below.

### 9.A Preconditions — must be resolved before implementation begins

These are no longer "open questions" — they are *required proof steps* before §4 can be turned into code.

1. **What is the correct scalar-φ thin-shell operator?**
   - BELFEM's existing `phi`, `phi_tri3`, `phi_tet4` kernels assemble `μ₀ ∫ ∇φ · ∇v dV` into `M()`, not `K()`. (Verified at `mt_maxwell_phi.cpp:24, 107, 119`.) The first draft of §4.5 was wrong; the corrected version is in place but hasn't been **dimensionally checked**.
   - Required: write down the bilinear form for the φ-layer, including the `dV = thickness · area` weight. Confirm units. Confirm the time-stepping assembly does the right thing when both `H`-blocks (which use both `M` and `K`) and `φ`-blocks (which use only `M`) appear in the same global system.
   - Estimated effort: 1-2 hours, pen and paper, plus a quick read of `IWG_Timestep::compute_jacobian_and_rhs` to confirm the M-and-K combination logic.

2. **What is the exact internal H–φ interface topology?**
   - The first draft claimed reuse of the external `hang_thinshell_edges_on_*` path. Codex showed this doesn't fit. The corrected §4.2 says new code is needed. *What* new code is the question.
   - Required: pick one ghost facet on a hypothetical H–φ layer interface and walk by hand: which `Element*` is the master, which is the slave, what are the local face indices on each, what's the temp-index correspondence, what does the constraint look like as a row in the T-matrix or as a contribution to a `K_internal_hphi` block.
   - Compare against the external H–φ path (the `hang_thinshell_edges_on_nodes_*` functions) to identify what's actually reusable and what isn't.
   - Estimated effort: half a day. This is the core mathematical work.

3. **What is the variational form of the contact-impedance interface kernel? (added 2026-04-07 after second Codex review, gates Option D)**
   - The first sketch in §5.2 / §3.4 calls for a Robin / sheet-impedance jump law `n × ⟦E⟧ = ρ_sheet · J_n` derived as the asymptotic h → 0 limit of a thin resistive layer with `ρ · h` held fixed. Codex's second review (`devlog/dl20260407_thinshell_three_way_coupling_review.md`) confirms this is plausible and literature-backed but explicitly **not** equivalent to BELFEM's current `h_ghost()` continuity kernel.
   - Required: derive the bilinear form from Schnaubelt 2023 §III (the T2TCL surface contribution), including:
     - The interface penalty / impedance term: what does the test-function contraction look like? Is it `∫_Γ ρ_sheet · J_n · v_n dS`, or does it involve `⟦H_t⟧` directly?
     - The consistency / adjoint-consistency terms: do they exist for this operator at all? (For an SIPG-style continuity penalty they're required for optimal convergence; for a Robin BC they may not be needed because Robin is already a natural BC.)
     - The discretization: which test/trial spaces does the kernel act on — Ag-top edge_h, YBCO-bottom edge_h, both? How does it interact with the static condensation at the H–H continuity interfaces *neighboring* the rint sideset?
   - Required: dimensional check on `ρ_sheet = ρ · h` (units `Ω·m²`) and confirmation that the assembled K/M entries have the right units to combine with the surrounding bulk operators.
   - Required: a test problem that has an analytical solution (e.g. uniform current crossing a planar interface with known sheet resistance) so the kernel can be validated end-to-end before being wired into a real tape stack.
   - Estimated effort: **1 day** for the derivation and the analytical-solution test design, plus **1-2 days** for the kernel implementation and validation. Significantly more if the consistency/adjoint terms turn out to be non-trivial.
   - Until this precondition is resolved, **Rint stays on the honest cap (Option A in §5.2 / Path 2 in §10).** No code shortcut via "reinterpret α in `h_ghost()`" — the math is not the same.

All three preconditions should be settled (and ideally documented as their own short notes in `src/fem/maxwell/doc/`) before any code is written. Preconditions 1 and 2 are required for the φ formulation (§4); precondition 3 is required for Option D (the contact-impedance interface for Rint), which is a separate milestone.

### 9.B Open design questions — to be answered as part of the implementation

3. **Is the threshold `ρ_φ_threshold = 1e-2 Ω·m` for H vs φ the right number?**
   - The literature limit is `1e-4` (Alves 2022b substrate). Anything between `1e-4` and `1e-2` is the "high-resistance conductor" zone — gray area.
   - We want the threshold above Rint (1e-3) and below buffer (1.0), so `1e-2` is the obvious choice. Is there a more principled value? Probably not — this is a matter of where the user wants to switch from "resistive bad conductor" to "true insulator".

4. **Cohomology cuts inside a φ-layer in a wound coil:** how does the existing `CutFactory` machinery interact with a φ-block that lives inside a `ThinShell` sideset rather than in a regular volume air block? Worth a separate read of `cl_CutFactory.cpp` before committing.

5. **Per-material `formulation` override + per-tape `coupling` override** in the input file: yes/no?

   **Resolved (after Junie's review and the user decision on 2026-04-07 to default to condensation):** **both, with explicit override taking precedence over the auto rule.** Two separate switches:

   - **Per-material `formulation`**: which DOF kind the layer uses (`h`, `h_capped`, `phi`).
   - **Per-tape (or global) `coupling`**: which interface coupling method is used at H–H layer interfaces (`condensation` is the default; `nitsche` is the opt-in fallback that re-enables the existing `h_ghost` + ghost-facet path).

   This matches the design pattern of "sensible default + explicit escape hatch" used elsewhere in BELFEM. Concretely:

   ```
   materials
   {
       buffer
       {
           usermat { ... }
           // optional explicit override; if absent, auto rule applies
           formulation : phi ;
       }
       rint
       {
           usermat { ... }
           formulation : h_capped ;     // honest cap, see §5.2
       }
       copper
       {
           builtin : copper ;
           // no formulation key — auto rule kicks in,
           // sees rho < 1e-4, picks H (uncapped)
       }
   }

   topology
   {
       thinshell : tape1
       {
           sidesets : 5,6,7,8,9,10 ;
           // optional per-tape override of the layer-interface coupling
           // method. Default is `condensation`. Setting `nitsche` re-enables
           // the legacy h_ghost + ghost facet sideset path for H-H layer
           // interfaces (the work done in 2026-04 around regularized harmonic
           // mean penalty, SIPG fix, etc.). Selectable per tape; if absent,
           // the global default below applies.
           coupling : condensation ;
       }
       thinshell : tape2
       {
           sidesets : 11,12,13,14,15,16 ;
           coupling : nitsche ;       // tape2 uses the legacy Nitsche path
       }
   }

   solver
   {
       ...
       // optional global default for layer-interface coupling, overridable
       // per-tape via the `coupling` key above. Defaults to `condensation`.
       thinshell coupling : condensation ;
   }
   ```

   The auto rule (default behaviour, no input file change required for existing cases):

   - `ρ < 1e-4 Ω·m`  → `formulation : h` (pure H, no cap)
   - `1e-4 ≤ ρ < 1e-2 Ω·m` → `formulation : h_capped` (pure H, **honest cap** to `1e-4` per §5.2)
   - `ρ ≥ 1e-2 Ω·m`  → `formulation : phi`  (φ-formulation per §3 / §4)

   The explicit override lets a user *force* a specific treatment when the auto rule picks the wrong one for their physics — e.g. if their Rint is *physically* an insulator and they want it modeled as one (override to `phi`), or if they have an unusual ferromagnet they want treated as an H-block even though `ρ` is right at the threshold.

   The auto rule applies the moment the input file is parsed (in `MaxwellFactory::create_and_assign_materials` or similar), and the result is logged via `gLog` so the user sees which layer got which treatment without having to chase it in code. Logged at `InfoLevel::Detailed`, format like:

   ```
   Layer rint1 (rho = 1.00e-03 Ohm.m) auto-classified as h_capped (rho_eff = 1.00e-04 Ohm.m)
   Layer buffer (rho = 1.00e+00 Ohm.m) auto-classified as phi
   Layer copper1 (rho = 4.00e-09 Ohm.m) auto-classified as h
   ```

   The thresholds (`1e-4`, `1e-2`) live as named constants in `cl_IWG_Maxwell.hpp` or similar, and can themselves be overridden through the `iwg.psi()` plumbing if a user needs to retune them globally for a different problem class.

6. **Is there a *direct* weak coupling between H and φ that bypasses the static-condensation route entirely?**

   **Resolved (2026-04-07, user decision):** moot. The user decided to default to **static condensation everywhere** because the BELFEM infrastructure for it already exists (the existing `hang_thinshell_edges_on_*` functions and the `[1, -1]` LINE2 T-matrix in `DofData::create_dofwise_t_matrices_master`). A Nitsche-style direct H–φ coupling would still be valid math, but it would require new derivation, new IWG branches, and new test infrastructure, while condensation reuses what's already proven to work for the external air–conductor boundary. Static condensation wins on minimal-disruption grounds.

   The user agreed that an opt-in Nitsche fallback (per-interface or per-tape, controlled by the input-file `coupling : nitsche` switch in item 5 above) is worth keeping in tree for debugging, comparison testing, and future research where non-conforming meshes or intentionally weak coupling might be wanted.

### 9.C Operational / pragmatic questions

7. **Honest cap vs uncapped postprocessing for Rint?** Codex flagged that the first draft was internally inconsistent. Resolved: use the **honest cap** (Option A in §5.2), with a `gLog` warning explaining the surrogate. The user has to know that reported losses correspond to `ρ_eff`, not the input ρ.

8. **Should the bulk-K cap at `1e-4` be applied as a stopgap while §4 is implemented?**
   - Yes — but with the honest-cap caveat from §5.2 (it's a numerical surrogate, not faithful physics).
   - The cap alone runs the buffer case without NaN. It does NOT solve the underlying physics problem; it just keeps the matrix conditioning under control.
   - The stopgap should be implemented in a separate commit from the H–φ refactor so the diff is small and easy to revert if needed.

## 10. Decision: defer to next session — REVISED twice on 2026-04-07

This is a non-trivial structural change. The first draft of this section recommended jumping straight from "design note" to "implement path 3". **Codex's first review made it clear that's premature** — there are two preconditions in §9.A.1–§9.A.2 that need to be settled first. **Codex's second review** (after the rint discussion) added a third precondition in §9.A.3 for the contact-impedance interface.

The revised plan separates the work into **three independent milestones**, each with its own preconditions:

### Milestone A — H–φ for true insulators (the buffer case)

**Preconditions:** §9.A.1 (scalar-φ operator) and §9.A.2 (internal H–φ topology).

**Steps for the next session:**

1. Resolve §9.A.1: pen-and-paper derivation of `μ₀ ∫ ∇φ · ∇v dV` on a PENTA6TS used as a scalar element, dimensional check, confirmation that `IWG_Timestep`'s M-and-K assembly does the right thing for mixed H/φ blocks. **1-2 hours.** Documented as a one-page note in `src/fem/maxwell/doc/`.
2. Resolve §9.A.2: hand-trace one ghost facet on a hypothetical H–φ layer interface, identify what's reusable from the existing `hang_thinshell_edges_on_*` functions vs what needs new dispatch. **Half a day.** Documented as a second one-page note with concrete `Element*`-level wiring.
3. Implement §4 (the H–φ refactor) — `ThinShellFactory` Layer struct, `MaxwellFactory` internal-interface dispatch, `phi_ts_insulator` kernel, field-list updates. **2-4 days.**

### Milestone B — Static condensation as the H–H default (the cleanup pass)

**Preconditions:** none beyond §9.A.2 (which doubles as the H–H internal interface topology).

**Steps:**

1. Implement the new `hang_thinshell_edges_on_edges_internal` helper described in §4.2.
2. Wire the `coupling : condensation|nitsche` input switch from §9 item 5.
3. Regression test: existing tape stacks should give bit-identical results in condensation mode and Nitsche mode (to round-off).
4. **Total: ~1.5-2 days, can run in parallel with Milestone A** since both touch `MaxwellFactory::create_thinshells` but in disjoint dispatch paths.

### Milestone C — Contact-impedance interface for Rint (Option D)

**Preconditions:** §9.A.3 (T2TCL Robin derivation). This is the new precondition added after Codex's second review and the most uncertain of the three.

**Steps:**

1. Resolve §9.A.3: derive the contact-impedance bilinear form from Schnaubelt 2023 §III, including consistency / adjoint terms (or proof they're not needed), dimensional check, and an analytical-solution test problem. **1 day.**
2. Implement the new `h_contact_impedance` kernel (or whatever name) as a separate matrix module — *not* by extending `h_ghost()`, which is mathematically a different operator. **1-2 days.**
3. Validate against the analytical test problem from step 1.
4. Replace the honest cap on Rint with the new interface kernel in the production input files.

**Milestone C is gated on the precondition being met.** If §9.A.3 turns out to be hard or the consistency-term derivation is non-trivial, this milestone slips while milestones A and B proceed.

### For the running case in the meantime

- **Path 1 (immediate, no code change):** user updates the input file to omit the buffer line, following Alves 2022a/2024 / Messe 2023 / Riva 2023 convention. The case runs as-is.
- **Path 2 (small code change, small commit):** add the bulk-K cap as an **honest cap** (per the revised §5.2): cap `ρ` to `ρ_cap = 1e-4` *everywhere* — the K assembly, the `mat->rho()` callsites, the postprocessor — with a `gLog` warning explaining the surrogate. The Rint case converges; loss/J/Jc values correspond to the capped ρ; the user is informed.

**Path 2 is the stopgap for Milestone C.** Once Milestone C lands, Path 2 is removed and Rint goes through the contact-impedance interface kernel instead.

**Junie's empirical result confirms Path 2 is sufficient as a stopgap:** unbounded ρ gives a +1544 dB residual (NaN); bounded ρ gives a stable −156 dB residual (working). 1700 dB swing. Path 2 is *not* a guess — it's a measured working configuration.

Path 1 + Path 2 together get the immediate test cases running. Path 3 (this document) is the long-term physically-faithful answer for the buffer case and should land as its own milestone after the preconditions in §9.A are resolved.

---

**End of design note.**

- **First draft** 2026-04-07 (Claude).
- **Codex review #1** 2026-04-07 (devlog/dl20260407_thinshell_hphi_formulation_review.md). Caught four issues; revisions in §3.2, §4.2, §4.5, §4.6, §5.2, §9, §10.
- **Junie audit** 2026-04-07. Empirical confirmation of the diagnosis (1544 dB → −156 dB swing, see §1.1) and resolution of the auto-vs-explicit open question (§9 item 5). The diagnosis and direction are now **independently confirmed by three audits**.
- **User decision** 2026-04-07. Default to static condensation everywhere; keep `h_ghost()` Nitsche as opt-in fallback. Revisions in §3.2, §4.2, §9 item 5, §9 item 6.
- **Codex review #2** 2026-04-07 (devlog/dl20260407_thinshell_three_way_coupling_review.md). Flagged the conflation between the H–H continuity Nitsche fallback and a true contact-impedance Robin law. The contact-impedance case is **not** a knob on `h_ghost()` — it requires its own derivation from Schnaubelt 2023 §III T2TCL. Revisions: new §3.4 (three-way categorization), new Option D in §5.2, new precondition §9.A.3, new Milestone C in §10. Path 2 (honest cap) explicitly retained as the stopgap until Milestone C lands.
