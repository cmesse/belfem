# Ghost-Penalty Stabilization for Thin-Shell Layer Interfaces {#fem_maxwell_ghost_penalty_stabilization}

**Date:** 2026-08-24
**Purpose:** Explain what the ghost-penalty kernel `maxwell::h_ghost()` solves, how it works (with literature references), why its two parameters have their current defaults, and which deck block can override those defaults.
**Module:** src/fem/maxwell

---

## 1. The problem it solves

A stacked thin-shell conductor (REBCO tape, CORC cable) is modeled as several stacked
`PENTA6TS` layers — superconductor, buffer, metal stabilizer — each carrying its own
tangential-`H` (edge) field. Across the interface between two adjacent layers the
tangential `H` should be **continuous** (Ampère: no surface current sheet sits *between*
the modeled layers), but nothing in the per-layer edge spaces enforces this: the two
layers have independent DOFs that only meet weakly.

Without a coupling term the interface is under-constrained. Two failure modes appear:

- **Spurious tangential jumps / checkerboarding** at the layer interface, especially
  where one layer is superconducting (`ρ → 0`, near-singular stiffness) and the
  neighbor is resistive — the classic discontinuous-coefficient pathology.
- **Loss of robustness** during current sharing / quench, when the resistivity contrast
  across the interface sweeps over many orders of magnitude as the front moves.

`h_ghost()` adds a **weak interface coupling** (an interior-penalty / Nitsche term) that
enforces tangential-`H` continuity across the layer interface and stabilizes the
discrete operator against the resistivity contrast. It is registered for
`DomainType::Ghost` groups (`cl_IWG_Maxwell.cpp`, the `DomainType::Ghost` case) and runs on the "ghost" facets
that sit between stacked thin-shell layers.

## 2. How it works

The kernel (`h_ghost()` in `src/fem/maxwell/matrices/mt_maxwell_h.cpp`) assembles a
**symmetric weighted interior-penalty (Nitsche)** bilinear form on each ghost facet,
coupling the master (`+`) and slave (`−`) layer DOFs through four blocks
`K++, K+-, K-+, K--`. At each integration point, each block carries:

- a **penalty term** `α · ⟨Eᵀ E⟩` that penalizes the tangential-`H` jump, and
- the **consistency and adjoint-consistency terms** `∓ ρ_harm · ⟨Eᵀ D + Dᵀ E⟩` (the
  flux terms that make the penalty *consistent* — they vanish for the exact solution),

where `E` is the edge interpolation operator, `D` the per-layer normal-derivative
operator scaled by layer thickness `h`, and `ρ_harm` the harmonic mean of the two layer
resistivities. The adjoint-consistency terms mirror the consistency terms block by
block, so the assembled operator is **exactly symmetric**: `K-+ = (K+-)ᵀ`
term-for-term. The `GhostElementContract` test
(`tests/fem/test_InterfaceOrientation.cpp`) measures `max|K − Kᵀ| = 0.0` on real edge
functions. The kernel began as the nonsymmetric variant (consistency terms only) and
was deliberately symmetrized in April 2026 by adding the adjoint terms; the resulting
symmetric-with-harmonic-weighting scheme is the SWIP method of Ern, Stephansen &
Zunino (2009). Ruled intentional 2026-08-28: the symmetric form is the formulation of
record.

The **penalty coefficient** is:

```
k = ρ / h                       per-layer interface stiffness
k_reg_side = k + k_reg          regularized per-layer stiffness
k_pen = 2 · km_reg · ks_reg / ( km_reg + ks_reg )    regularized harmonic mean
α = eta · k_pen
```

The plain harmonic mean `2·km·ks/(km+ks)` is the standard heterogeneous-DG weighting
(**Burman & Zunino 2006**), robust to the coefficient contrast. It is bounded above by
`2·min(km, ks)` (well-behaved for an insulating neighbor, `ρ → ∞`), but it **collapses
to zero** when one layer is superconducting (`ρ → 0`). Adding the regularization offset
`k_reg` to each stiffness before taking the mean keeps the coefficient bounded in both
singular limits:
`k_reg → 0` recovers Burman–Zunino; both `k` tiny → `α ~ 2·eta·k_reg` (no collapse); one
`k` huge → `α ~ 2·eta·(other_k + k_reg)` (no blow-up).

**Literature.**
- Symmetric weighted interior penalty with harmonic flux weighting (the scheme
  implemented here): Ern, Stephansen & Zunino (2009), *IMA J. Numer. Anal.* 29:235–256.
- Interior-penalty family and coercivity: Brenner & Scott §10.5
  (`literature/books/brenner.txt`). Unlike the nonsymmetric variant, which is coercive
  for any `eta > 0`, the symmetric form is coercive only for `eta` above a mesh- and
  element-dependent threshold set by the trace inequality.
- Harmonic / weighted averaging of discontinuous diffusion coefficients:
  Burman & Zunino (2006), *SIAM J. Numer. Anal.* 44:1612–1638.
- Coefficient regularization (floor/ceiling on resistivity, same *principle* as `k_reg`):
  Dular et al. (2021), `literature/papers/fem/dular2021.txt:1109`.

**Restriction.** The kernel supports only first-order `PENTA6TS` facets: LINEAR
interpolation is asserted in `h_ghost` in `mt_maxwell_h.cpp`, and block sizes
come from `number_of_nedelec_dofs`. That gives 6 DOFs per TS element, 3 per
facet. Higher order requires adapting the `Dm/Ds` assembly.

## 3. The two parameters

Both constants live in `mPenalty` (`cl_IWG_Maxwell.cpp`, constructor) and are
user-settable since 2026-08-24 through the optional deck block
`nonlinear magnetic { nitsche ghost penalty { eta : … ; k_reg : … Ohm ; } }`
(see `doc/input_file_reference.md` §4.2.1). Since 2026-09-01 the ghost is
**opt-in**: when the block is absent, or states `eta : 0`, the thin-shell factory
creates neither the duplicate interface dofs nor the ghost facets, the layers
share their interface edges, and nothing in this document is assembled. With
`eta > 0` the values below apply (`eta` must then be stated; `k_reg` keeps its
default when omitted). The switch is read once for the factory, the controller
and the mesh cache tag (`src/fem/kernel/fn_FEM_ghost_switch.hpp`):

| Symbol | Slot | Default | Role |
|--------|------|---------|------|
| `eta`   | `penalty(0)` | **4.0** | dimensionless Nitsche stabilization multiplier |
| `k_reg` | `penalty(1)` | **1e-3 Ω** | stiffness-regularization offset (SC-collapse floor) |

**How the defaults were chosen.** Both are deliberate engineering defaults, *not* values taken from a specific
paper. An exhaustive search of curated `./literature`, web results, and
independent Grok/Codex scans found no published value or formula for either default. The *method* and
the *harmonic-mean weighting* are literature-grounded; the specific multiplier and offset
are calibrated to the HTS thin-shell problem class:

- `eta = 4.0` — for the symmetric scheme, coercivity requires `eta` above a mesh- and
  element-dependent threshold (the trace-inequality constant), and no published value
  exists for `PENTA6TS` ghost facets. `eta = 4` is therefore an engineering default,
  low-to-moderate compared with the usual SIPG default (~10), validated in practice on
  the supported problem class (point 3 below).
- `k_reg = 1e-3 Ω` — sits at the low end of the typical HTS-tape interface stiffness
  `k_conductor ~ 1e-3..1e-2 Ω` (Cu/Ag/Hastelloy at the few-micron scale): invisible at
  metal–metal interfaces, large enough to keep a superconducting layer from collapsing
  the harmonic mean, small enough to be dominated by real metal physics wherever that
  physics is present.

Both values are empirical and have no published derivation behind them. Two calibration
routes remain open: a mesh-aware `eta` taken from the element trace-inequality constant
à la Shahbazi 2005, and a physics-tied `k_reg = ρ_floor / h`. Neither is covered by a
contrast-sweep regression test yet.

**Why the defaults should normally stand.** The deck keys are for calibration
work and controlled experiments, not routine tuning:

1. **They are not physical inputs.** Unlike resistivity, `Jc`, or layer thickness, `eta`
   and `k_reg` are discretization-stabilization constants. Their *correct* values depend
   on the element type and mesh (the trace-inequality constant), not on the device being
   modeled — a physical argument alone gives no basis for changing them.
2. **Bad values degrade silently.** `eta` too small → weak jump control and
   checkerboarding; `k_reg` too large → over-stabilization that smears the SC↔normal
   front; both wrong → poor conditioning that surfaces only as a slow or failing solve
   (cf. the t69 breakdown). None of these raise an error. The parser rejects only the
   hard misconfigurations: `eta ≤ 0`, `k_reg ≤ 0` (an SC–SC facet would evaluate `0/0`
   in the harmonic mean with the unfloored power law), a missing unit on `k_reg`, and
   dimensioned values for the dimensionless keys.
3. **The defaults are validated for the supported problem class.** They are tuned for
   first-order `PENTA6TS` HTS thin-shell stacks — the only configuration the kernel
   currently supports.

If they ever need systematic recalibration (e.g. a new element order, or a mesh-aware
`eta`), the clean path is to derive the values from the element's trace constant inside
the kernel — the two calibration routes above — rather than relying on the raw deck
knob.

## See also

- `contact_impedance_theory.md` — the TSA contact-impedance alternative that eliminates the Nitsche penalty parameter entirely
- `mt_maxwell_h.cpp` — the `h_ghost()` kernel (search `h_ghost`)
- `tests/fem/test_InterfaceOrientation.cpp` — `GhostElementContract`, the executable contract for this kernel
