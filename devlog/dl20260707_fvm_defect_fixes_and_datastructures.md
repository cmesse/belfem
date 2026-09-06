# Devlog 2026-07-07 — FVM Defect Fixes + Data-Structure Design

**Date:** 2026-07-07
**Topic:** `src/fvm` MPFA-O — defect fixes D1/D3/D4/D6/D7, flux-measure verification, and the
global-solve data-structure design (§2.1/§2.2 of the plan)
**AIs involved:** Claude (implementation + design), Codex (two audit rounds, both high confidence)
**Claude Confidence:** high (fixes), high (measure `l = L/2`), medium→high (design, after Codex audit)
**Codex Audit Confidence:** high (both rounds)
**Literature References:** Aavatsmark 2002 (F0) Eq. 67, Fig. 9(a), `aavatsmark2002.txt:619-624,834-852`;
Agélas & Masson 2008 (F1) §5, `agelas2008.txt:257-264,305-308`; Messe et al. 2023 (paper1)
(static-condensation argument); MPFA.pdf (tmp/, meeting slides 2025-10-23)

## Summary

Execution start of `todo/fvm_implementation_plan.md`. All audited defects except D5 fixed in code
and Codex-verified; two previously unknown MPI bugs found and fixed in the same block; the
`S = t·l` flux measure confirmed against two independent sources; the global-solve data
structures decided with Christian and recorded in the plan (§2.1, §2.2). Design continues
2026-07-08.

## Key Findings

- **Two new latent MPI bugs** in the non-root branch of `Factory::compute_facet_normals`
  (neither in the 2026-06-26/07-05 audits): (1) the own-normals copy loop iterated *all* facets
  while the local arrays only hold owned columns — column misalignment past the first not-owned
  facet; (2) `tMyIDs` (sized for owned facets) was reused for the not-owned request list without
  resize — OOB write when not-owned > owned. Both fixed; Codex confirmed (high).
- **`S = t·l` measure question (Christian):** `l` is **half the full mesh-edge length** (`L/2`),
  *not* `|r|`. The `r` vectors live only in the gradient operator `G`. Canonical MPFA-O:
  sub-interface = vertex-to-midpoint half-edge (Aavatsmark 2002 (F0), Eq. 67, Fig. 9(a)); Agélas
  variant: face-measure split `m^s_σ = m_σ/card(V_σ) = L/2` regardless of the 1/3–2/3 continuity
  points. The code already stores the right value (`0.5·tLength` in `create_subedges`).
- **TPFA through-thickness "exact" was overclaimed** (Codex): exact only for flat/parallel-normal
  extrusions; on curved normal-extruded shells it is a K-orthogonal thin-shell approximation
  (Aavatsmark 2002 (F0), `aavatsmark2002.txt:834-852`). Plan wording corrected.
- **Strong-local Dirichlet confirmed** (decides plan O2): boundary-face elimination exactly as
  Agélas & Masson 2008 (F1) (`agelas2008.txt:257-264`). Guards recorded: BC tags per subedge
  (not per node), never eliminate a `Tc` and keep its flux row, always scatter the affine local
  RHS.

## Changes Made

- `src/fvm/cl_FVM_Factory.hpp` — D6: injective facet key `tA·N + tB`.
- `src/fvm/cl_FVM_Factory.cpp` — D1: slave normal negated in both branches (serial + distributed,
  applied by Christian); D4: distributed tangent `mC = mQ − mP`, non-root `set_size`, owned-index
  copy loop, `tMyIDs` resize; D7: both distributed lookups → `mMesh->facet(id)`; malformed
  `compute_surface_adjacencies` deleted; unused `tM` removed; `solve_local_systems` reads the
  stored `SubEdge` normals.
- `src/fvm/cl_FVM_SubEdge.hpp` — `master_normal()`/`slave_normal()` accessors (single source of
  truth for the D1 convention).
- `src/CMakeLists.txt` — `add_subdirectory( fvm )` (D3).
- `src/math/tensor/cl_Tensor.hpp` — 3rd-order ctor + accessors (Christian; resolves the D3
  `G(m,n,3)` mismatch).
- `src/fvm/cl_FVM_Dof.hpp/.cpp`, `cl_FVM_Resistor.hpp`, `fvmtest.cpp` — new scaffolds
  (Christian, evening): cell Dof carrying element/index/material; two-Dof resistor
  `R = (ρ₀l₀ + ρ₁l₁)/A` for layer conductances + ESATAN-style connectors.
- `src/fvm/doc/fvm_mpfao_theory.md` §2.3 — rewritten: Agélas 1/3–2/3 continuity points with
  citation (closes the D2 residual) + D1 normal-sign convention (slave normal computed in the
  slave's own tangent plane, then negated).
- `todo/fvm_implementation_plan.md` — §2.1 data-structure decisions (factory-moves-to-MPFAO,
  per-layer subcells mapping PENTA6, TPFA-in-thin-shell-limit, direct scatter instead of stored
  `C`, subedge-local BCs, no DofManager reuse) + §2.2 evening amendment (`fvm::Dof` +
  `fvm::Resistor`, material on the Dof, sparsity = node cliques + resistor pairs, two-target
  constraint Maxwell/nonfree-satellite); R3/R4 ticked, R5/R7/R8/R12/O2 reconciled, Codex prose
  pass applied.
- `todo/fvm_module_next_steps.md` — defect ledger: D1/D3/D4/D6/D7 ticked with fix notes.

## Open Questions

- Clean `-Wall -Werror` build of the module (Christian builds; R2 gate).
- `fvmtest` driver content (R1; file scaffold exists in `src/fvm/`).
- Data-structure design continues 2026-07-08: Dof/Resistor finalization, SubCell→Dof link,
  environment-Dof flavor, radiative-coupling sibling, MPFAO assembly interfaces.
- D5/R5 (local solve in `MPFAO`) — next code milestone after the design settles.

## Files Updated

- src/fvm/cl_FVM_Factory.hpp / .cpp
- src/fvm/cl_FVM_SubEdge.hpp
- src/fvm/doc/fvm_mpfao_theory.md
- src/CMakeLists.txt
- src/math/tensor/cl_Tensor.hpp (Christian)
- src/fvm/cl_FVM_Dof.hpp/.cpp, cl_FVM_Resistor.hpp, fvmtest.cpp (Christian, new)
- todo/fvm_implementation_plan.md
- todo/fvm_module_next_steps.md
