# Quadratic Thin-Shell Nédélec Edge Functions (EF_QUAD9TS, EF_PENTA18TS)

**Status:** Deferred (linear shells are the active production path).
**Severity:** Latent — fires only when quadratic thin-shell meshes are run with the Maxwell kernel. The factory currently throws a clear `BELFEM_ERROR` naming the missing class and pointing to this document.
**Effort estimate:** 1-2 focused sessions per element class (derivation + implementation + unit-circulation tests).

## Context

After the Option-B canonicalization (April 2026), linear thin shells (`QUAD4TS`, `PENTA6TS`) share facet topology and orientation layout with their volume counterparts. Quadratic thin shells (`QUAD9TS`, `PENTA18TS`) are defined at the mesh and element-template level — their facet layouts, orientation tables, and integration-point paths all canonicalize through `GeometryType::QUAD` and `GeometryType::PENTA` — but no `EdgeFunction` implementation exists yet.

Consequence: trying to solve a Maxwell problem on a mesh with `QUAD9TS` or `PENTA18TS` blocks will throw from `EdgeFunctionFactory::create_edge_function` during Calculator setup. The error message points here.

## What's missing

Two classes, both in `src/fem/interpolation/nedelec/`:

### `EF_QUAD9TS`

Quadratic 2D thin-shell Nédélec basis on a 9-node QUAD with 4 CCW facets (LINE3 bottom/top curves + 2 side-line facets; see `cl_Element_QUAD9TS.hpp`). `fn_num_nedelec_dofs.hpp` declares **3 DOFs**, vs. QUAD4TS's 2. The ansatz needs to be derived from the hierarchical p=2 Nédélec space restricted to the thin-shell reduction (no side-edge DOFs).

### `EF_PENTA18TS`

Quadratic 3D thin-shell Nédélec basis on an 18-node PENTA with 5 facets mirroring volume PENTA18. `fn_num_nedelec_dofs.hpp` declares **16 DOFs**. Linear sibling `EF_PENTA6TS` uses 6 DOFs (3 bottom-triangle edges + 3 top-triangle edges). The 16 for PENTA18TS is not an obvious factoring — it needs to come from the formulation document, not be guessed at.

## Derivation sources

- **Messe et al. 2023 (paper1)** — BELFEM's core architectural paper. §2 covers the thin-shell H-φ formulation and the hierarchical basis reduction that produces the `Q*TS` DOF counts. This is the primary reference for both classes.
- **Arsenault et al. 2023 (paper3)** — magnetodynamic H-φ coupling. Cross-reference for the edge-DOF structure consistent with how volume elements couple.
- Existing `src/fem/interpolation/nedelec/cl_EF_PENTA6TS.cpp` and `cl_EF_QUAD4TS.cpp` are the templates to mirror. Their code structure (precompute, link, E, C) is the contract any new class must honor.

## Implementation checklist

1. Derive the QUAD9TS basis from Messe §2.X with the 4-facet topology now in place (bottom LINE3 curve + 2 side lines + top LINE3 curve). Verify 3 DOFs match `num_nedelec_dofs`.
2. Write `EF_QUAD9TS.hpp/.cpp` following `EF_QUAD4TS`'s contract. Verify unit circulation per edge DOF.
3. Register in `EdgeFunctionFactory::create_edge_function` (replace the current "not implemented" error branch with `return new EF_QUAD9TS();`).
4. Repeat for PENTA18TS.
5. Unit tests: add cases to `tests/fem/test_FacetIntegrationPoints.cpp` (if applicable at the edge-function level) and any edge-circulation test infrastructure.

## Related

- `src/fem/interpolation/cl_EdgeFunctionFactory.cpp:62-80` — current error branch pointing here.
- `devlog/dl20260422_thinshell_facet_renumbering_break.md` — original audit that kicked off the Option-B canonicalization.
- `src/fem/interpolation/nedelec/fn_num_nedelec_dofs.hpp` — DOF counts that the new classes must honor.
