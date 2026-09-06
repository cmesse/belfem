# The Bearing Gauge and the Near-Null φ Eigenmode {#fem_kernel_bearing_gauge_eigenmode}

**Date:** 2026-08-06
**Purpose:** Why a single-node gauge pin ("bearing") leaves the Jacobian of an h-φ
system near-singular, why that paralyzes Newton while Picard is immune, and what the
remedies are.
**Module:** `src/fem/kernel` (`cl_FEM_Bearing`, `cl_FEM_DofMgr_BearingData`), with
consequences for the nonlinear Controller (see `nonlinear_controller_theory.md`)

---

## 1. The gauge freedom and the bearing

In the h-φ formulation the magnetic scalar potential φ in the air region enters the
physics only through ∇φ: the solution is determined up to an additive constant per
connected air component. BELFEM pins that constant strongly at a single mesh vertex —
the **bearing** (`bearing { nodes : <vertex id> ; }` in the input deck), implemented as
a Dirichlet fix of the node's φ dof (`Bearing::impose_dirichlet`). An id that does not
resolve is **no longer silent**: an empty bearing is legitimate on a worker rank, since
linking is master-only, but on the master it raises `BELFEM_ERROR`
(`cl_FEM_Bearing.cpp:104-120`). When the id resolved to nothing at all it names the deck's
`nodes` list without echoing the id; when the point exists but carries no dof of the requested
field it names the point and node. The old silent return is what let a dead bearing run for
months (DR-126).

The pin is *correct*: in exact arithmetic one point condition removes the constant
nullspace. The problem is its *strength*. A point Dirichlet is a measure-zero
constraint for a 3-D Laplacian: the discrete mode

> φ = c everywhere, forced to 0 at the bearing node

costs only O(h)·c² of energy, where h is the mesh size at the bearing (often the
coarsest region — the far boundary). The stiffness/Jacobian therefore keeps one
**near-zero eigenvalue** whose eigenvector is essentially the constant-φ mode, and the
condition number grows without bound under mesh refinement or coefficient contrast.

## 2. Why Newton dies and Picard does not

**Newton** solves J·Δx = r for a *correction*. In exact arithmetic r is orthogonal to
the nullspace; in floating point it carries a roundoff-level component along the mode,
which the solve amplifies by 1/λ_min. The direct solvers do not fail — their
tiny-pivot replacement (STRUMPACK) or pivot patching (MUMPS) returns SUCCESS — they
return a Δx dominated by a huge, physically meaningless constant-φ content. Observed
on the Garber CORC deck (2026-08-06, `BELFEM_PROBE_NEWTON_DX`):

- ‖Δx‖ ≈ 102 for ‖r‖ ≈ 5.7·10⁻⁶ (amplification ~2·10⁷);
- Δx components ≈ 0.43, *near-identical* across tens of thousands of φ nodes
  (the constant mode), largest entries in the finest-meshed region;
- the residual frozen to within ±0.01 dB while ω swept 0.5 → 0.004 — the mode is
  invisible to the lagged operator, so the line search can neither use nor reject
  the step, and the ill-conditioning destroys the *useful* part of Δx as well
  (the whole factorization loses its digits, not just the null component).

The result is the "Newton does nothing" signature: bit-flat accepted iterates at any
relaxation, until the stagnation guard falls back to Picard. Note the accepted
iterates *do* drift along the mode — harmless if the mode is pure gauge.

**Picard** solves A·x_G = b with the *physical* right-hand side. The solution's
component along the near-null mode is set once by the factorization (an arbitrary but
stable constant offset) and cancels from every physical quantity; the fixed-point
iteration converges as if the mode did not exist.

This behavior is independent of thin-shell edge fusing and identical on all branches;
it is a property of the gauge pin, not of any particular formulation detail.

## 3. What does NOT work: surface gauging with net transport current

The tempting fix — Dirichlet φ = 0 on the outer boundary sidesets — is only admissible
for problems with **zero net current** (background-field/screening problems). With net
transport current I, the far field is H_θ = I/(2πr) and any loop around the conductor
carries ∮H·dl = I: φ is multivalued out there (jump I across the cut, φ ≈ −I·θ/2π
along the boundary). A constant-φ surface condition clamps that circulation to zero
and fights the current constraint. The point bearing exists precisely because it pins
only the constant and makes no field statement.

## 4. Remedies

| Remedy | Status | Notes |
|---|---|---|
| **Run pure Picard** (`algorithm : Picard`) | RECOMMENDED for net-current h-φ decks (2026-08-06) | Immune to the mode; the hybrid controller's Newton stage buys nothing when every Newton solve is blind |
| Mean-value gauge row (Σφ = 0 as one dense constraint row, replacing the point pin) | future work | Pins the same constant with O(1) stiffness; same infrastructure class as the free-cut λ rows; physically identical to the bearing |
| Deflation (project the known constant mode out of Δx after the Newton solve) | future work | Cheap, but the factorization has already lost precision at solve time — inferior to fixing the operator |
| Absolute magnetic tolerance (`absolute tolerance`; **opt-in — consumed when set, disabled at 0.0 by default**, `cl_FEM_Controller.hpp:203-211`) | partly landed | Would stop iterating when ‖Ax−b‖ reaches machine floor; addresses a *second*, benign Newton-flat regime (solver returns Δx = 0 on machine-epsilon residuals), not the near-null mode itself |

**Caveat (open):** the observed mode has not been fully discriminated between (a) the
global constant against the point pin, and (b) a *bore* mode — the air inside a CORC
former couples to the pinned outer air only through the narrow tape gaps, so its own
constant is nearly floating regardless of the bearing. Both produce the same Newton
signature and both are gauge-like for Picard; they differ in remedy ((a) mean-value
row; (b) a second bearing inside the bore, at an error of the order of the neglected
gap coupling). A one-shot dump of the full Δx mode shape (unimodal vs bimodal in φ)
decides; the probe used on 2026-08-06 has been removed, so it would have to be
re-instrumented in the Newton branch of `cl_FEM_DofMgr_SolverData.cpp`.

## 5. Recognizing the mode

The signature, should it need re-confirming (the diagnostic probe used in the
2026-08-06 investigation has been removed again): print ‖Δx‖ of the Newton solve next
to ‖r‖. A ‖Δx‖ orders of magnitude above ‖r‖ (observed: 2·10⁷×), carried by
near-identical values on the φ dofs of most of the mesh, with the reported residual
frozen under any relaxation sweep, is this mode. Newton iterates in that state are
accepted-but-inert; Picard from the same state converges normally.

## 6. Literature

- Messe et al. 2023 — h-φ formulation and cut constraints; the gauge constant
  and its pinning are implicit in the formulation.
- Alves et al. 2022b — cohomology cuts; the multivalued φ around net-current
  conductors that rules out surface gauging (Section 3 above).
- Bathe, §8.2.4 — conditioning of near-singular systems and the loss of accuracy of
  direct solves.
