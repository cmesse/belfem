# Free periodic cut dof: dense graph row, uint16 overflow, and the physics verdict

**Date:** 2026-08-04 (investigation started 2026-08-03)
**Purpose:** Root-cause and resolution of the `populate_graph` crash in the
greg CORC problem; formulation decision on the surplus periodic cohomology
generator; graph-layer hardening.
**Modules:** math/graph, fem/kernel, fem/maxwell, homology
**AI collaboration:** Claude (primary), Codex + Grok (three audit rounds via
`tmp/ai_exchange/greg_corc_graph_overflow_*.md`)

## Symptom

`hphirun` on `cmake-build-debug/greg` (3-tape CORC, twisted periodic BC,
thin-shell tapes, thin cuts) aborted at the first timestep:

```
Expected 265727 vertices but got 3583 (MatrixType=1)
cl_FEM_DofMgr_SolverData.cpp (line 861), SolverData::populate_graph()
```

## Root cause

Three layers, innermost first:

1. **Messenger:** `graph::Vertex::mVertexCounter` was `uint16_t`;
   `insert_vertex()` had no overflow guard. 265727 mod 2^16 = 3583 — the
   counter wrapped four times while inserting a genuinely 265,727-entry
   adjacency row.
2. **The row:** dof 358252 = abstract node 179126 × 2 dof types — the
   **4th cohomology generator** (`corc.bfm /nodes/abstract` =
   179123…179126). Its cut lies on the periodicity plane; through the cut
   duplicate sources (`cl_CutSet.cpp`, duplicates list their cuts' abstract
   nodes with weight 1) and hanging-dof condensation, the lambda appears in
   the dof list of every element adjacent to the plane → a mesh-scale
   adjacency row. Abstract nodes are mesh-bare by design (no elements, no
   sources) — a probe showing all zeros there is *not* an orphaned node.
3. **Why 4 generators but 3 currents:** `IWG_Maxwell::set_currents()` fixes
   the first `aI.length()` abstract dofs (the terminal-linked generators come
   first in the cohomology transformation, `cl_Cohomology.cpp:1240-1272`);
   generators beyond the prescribed currents remain free unknowns.

## Physics verdict: the 4th generator is legitimate and must stay free

A periodic geometry has one more generator than its non-periodic
counterpart: across the periodic identification it is ∇φ = H that must be
continuous, not φ itself. The φ-jump per cell is ∮H·dl along the axial
cycle — the axial magnetomotive force — and the 4th generator's cut carries
exactly that dof. Homotopic longitudes share the *linked current*, which for
3 same-handed tapes at `num_pitch = 0.01` scales as ~3·I·num_pitch ≈ 2.7 A
peak (I = 90 A) — small but not zero, so prescribing 0 A would make the bore
axial field first-order wrong.

An intermediate attempt that fixed surplus cut dofs to 0 A was tested and
**reverted**: besides the physics objection, this problem also builds the
full (n+m)² matrices (`mUseFullForce` path), whose graph includes fixed
rows — the dense row reappeared in `FullMass` regardless. With an electrical
circuit attached, a cut current can legitimately be a free unknown, so no
warning is emitted for unprescribed cuts either.

The condensed equation a free abstract dof receives (verified independently
by both auditors) is the Faraday/flux balance of the cut cycle, assembled as
Tᵀ M T over all cut-adjacent elements — consistent with the corrected
air-domain weak form (Arsenault et al. 2026 erratum to Arsenault et al. 2023
(paper3); BELFEM's air operator is the magnetodynamic mass term). Known
caveats of the free-cut path: the air block contributes only M, so the
lambda equation degenerates in the magnetostatic/large-Δt limit, and a
mesh-spanning cut produces an inherently dense matrix row (root-separator
material for the direct solver; symrcm bandwidth is meaningless for that
row).

## Changes

- `src/math/graph/cl_Graph_Vertex.hpp/.cpp` — `mVertexCounter` widened
  `uint16_t` → `uint32_t`; overflow asserts in `insert_vertex()`,
  `increment_vertex_counter()`, `init_vertex_container()` check the 32-bit
  limit; `<limits>` included.
- `src/fem/kernel/cl_FEM_DofMgr_SolverData.cpp` — the two `tOffsets`
  lookups in `compute_dof_dof_connectivity()` switched from
  `Map::operator[]` (which default-inserts 0 on a missing key — a silent
  corruption path: offset 0 reads the element count as a dof count) to the
  asserted `Map::operator()`.
- `src/math/graph/doc/graph_usage_guide.md` — counter type updated.
- `src/fem/maxwell/cl_IWG_Maxwell.*` — net unchanged (tail-fix and
  diagnostic added, then reverted by decision).

## Validation

Run proceeds past graph construction, matrix allocation, and into the
solve. Expected check once periodic steady state is reached: the value on
abstract node 179126's dof (the axial MMF per cell) should oscillate at
~2.7 A peak in phase with the 90 A drive.

## Open items

- `IWG_Maxwell::mPrescribedCurrents` is declared and read in
  `collect_abstract_node_dofs()` but never filled anywhere in-tree; the
  accessor `prescribed_currents()` has no callers. Either wire it or remove
  it.
- `set_currents()` maps prescribed currents to abstract dofs by prefix
  order; this silently assumes the abstract-node list order matches the
  cohomology transformation row order (terminals first). Worth a one-time
  consistency assertion.
- `mPairVerdict` in `CutProcessor` is built but never consumed; the
  "Step 6c" comments in `cl_CutSet.cpp` describe unimplemented registration
  of fragmentation extra duplicates. Unrelated to this crash; still debt.
- Free-cut solver economics on larger meshes: one dense row per free cut is
  tolerable for a direct solver, but monitor STRUMPACK memory; a
  border/Schur treatment (free cut lambdas outside the adjacency graph)
  remains the clean long-term architecture if this becomes a bottleneck.
- Magnetostatic/large-Δt limit leaves a free cut lambda without a
  constraining equation (air contributes only M); guard if such runs appear.
