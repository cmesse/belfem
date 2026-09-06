# Fix Plan: Parallel Source-Closure Gap in Mesh Distribution

**Date:** 2026-06-20
**Purpose:** Fix the distributed-mesh inconsistency that causes the z-periodic CORC model to abort in `MUMPS::solve()` on a worker rank, while the same model solves correctly in serial.
**Module:** src/mesh (Distributor), src/fem/kernel (Kernel::distribute_mesh)
**Status:** PLAN — tri-AI audited (Claude + Codex + Grok). Pending user approval. **5a (confirm-first) is mandatory before any selection-logic change.**

---

## 1. Symptom

4-proc CORC (z-periodic HTS cable, cohomology thin cuts, thin-shell layers) reaches `MUMPS::solve()` and aborts with `INFO(1) = -1` ("error on processor `INFO(2)`") on a worker. The identical model solves in serial. The MUMPS wrapper inspects only rank 0's INFO (`cl_SolverMUMPS.cpp:298-308`), so the worker's true error code is not surfaced.

## 2. Root cause (confirmed by three independent scans)

Claude (trace), Codex, and Grok — each working without seeing the others' conclusions — independently converged at high confidence (~98%) on a flag-ordering bug in `Distributor::select_entities()` (`src/mesh/cl_Mesh_Distributor.cpp`):

| line | action | element flag state |
|------|--------|--------------------|
| 290 | `unflag_all_elements()` | all cleared |
| 311 | owned elements → `mElementBitset` only (no `flag()`) | still cleared |
| 404 | edge `select_sources(edge)` loop, gated `if ( tElement->is_flagged() )` | **runs on zero flagged elements → dead** |
| 424/428 | face `select_sources(face)` loop, same gate | **dead** |
| 497 | `tElement->flag()` (first and only), from the owned+aura `mElementBitset` | too late |
| 518 / 543 | edge/face bit loops set `mEdgeBitset` / `mFaceBitset` + nodes | no `select_sources` |
| 560-577 | fixpoint closure loop calls `select_sources` on NODES (and duplicates) only | edges/faces/facets never source-closed |

The precise defect is **incomplete transitive / off-element source closure**:

- A source entity that lives **on a ghosted element** still reaches the worker, because the post-aura edge/face bit loops (`:518-558`) add every edge and face of each flagged element to the bitsets. This fallback path keeps distribution from aborting outright.
- A source entity that is **not** an edge/face/facet of any ghosted element — i.e. an off-element or transitively-introduced source — is never selected, because the only code paths that would select it (`select_sources` on the edge/face/facet that owns the source, or a transitive re-pass) never run.

(An earlier draft said edge/face sources are "never ghosted" and the dead loop is "the only mechanism"; Grok refuted both. The on-element fallback above is why the run survives distribution and reaches the solver at all.)

## 3. Why this can produce the MUMPS failure (plausible, not yet confirmed)

The causal links are verified, but the conclusion is not; confirmation is gated on 5a:

- Periodicity attaches EDGE/FACE sources to slaves (`cl_Mesh_Periodicity.cpp:100-141`); thin-shell attaches EDGE→EDGE sources (`cl_MaxwellFactory.cpp:1562-1581, 1656-1676`). CORC exercises both.
- `expand_hanging_basis_sources()` deliberately skips edges/faces (`cl_Mesh.cpp:2135`), so there is no source-flattening fallback for edges as there is for nodes.
- The DOF layer dereferences the source edge's DOF directly, with no NODE-style "all sources found" guard (`cl_FEM_DofMgr_DofData.cpp:3561, 3673-3693`). An incomplete or wrong edge constraint on a worker can yield an inconsistent / structurally-singular distributed matrix that surfaces only at `MUMPS::solve()`.

**Why this is not yet "confirmed" (Grok):** distribution *completes* and reaches MUMPS. If a worker were missing a source entity outright, `ProtoMesh::create_t_matrices()` would abort earlier via `basis(EDGE,id)` → `Mesh::edge(id)` → `BELFEM_ERROR` (`cl_ProtoMesh.cpp:872`, `cl_Mesh.hpp:1236-1240`). Reaching MUMPS therefore points to **partial** closure (entities present, constraints wrong/incomplete) rather than a wholly missing entity — or, less likely, a non-structural MUMPS error (`-9` workspace, `-13` allocation). The 5a diagnostic plus one-time surfacing of proc 3's true `INFO(1)` will decide which.

## 4. Secondary gaps (same class)

- **Aura facets are not source-closed.** `select_sources(tFacet)` runs only for *owned* facets (`:321-331`); facets pulled in by facet aura (`:464-467`) or introduced as sources inside `select_sources` get no source-closure pass (Codex + Grok).
- **Element / control-point target sources are not selected.** `populate_t_matrices` supports element and control-point targets, but selection never calls `select_sources` on the entities bearing those sources.
- **ELEMENT vs CELL type mismatch (Codex).** `select_sources()` has an `EntityType::ELEMENT` case, but `Element::entity_type()` returns `EntityType::CELL` (`cl_Element.hpp:788-792`); there is no `CELL` case in `select_sources()` and none in `Mesh::basis()` (`cl_Mesh.hpp:1768-1801`). A genuine element/cell source would therefore hit the `default:` branch → `BELFEM_ERROR`, not be ghosted. This is latent for CORC (no element-typed sources today), but it must be handled or the ELEMENT/CELL convention normalized before the plan can claim complete source closure.

The unified closure in 5b covers the first two at no extra cost; the CELL mismatch needs an explicit case or a guard.

## 5. Proposed fix

### 5a. Confirm first — MANDATORY (temporary diagnostic)

Before changing any selection logic, prove the gap empirically:

- [x] **5a-1** Source-closure check added at the end of `select_entities(aTarget)` (`cl_Mesh_Distributor.cpp`). **RESULT (confirmed):** procs 1 and 3 each report **29 missing sources**, proc 2 reports 0. Every miss is a **NODE source of a hanging EDGE** (e.g. `edge 565254` → missing nodes `19006`/`19026`) — the thin-shell/interface edges that hang on surrounding nodes. This directly confirms the dead-loop gap: `select_sources(edge)` never runs, so a hanging edge's sources are never ghosted. Refinement vs the original wording: the missing sources are NODE-typed (edges hanging on nodes), not edge-typed; the 5b fix is unchanged (calling `select_sources(edge)` ghosts these node sources).
- [ ] **5a-2** (deferred) Per-rank MUMPS `INFO` print added (`cl_SolverMUMPS.cpp`), but **not yet captured**: with the kernel volume-exchange fix reverted, `compute_element_volumes` aborts at element-id-0 *before* the solve, so MUMPS is never reached. Deferred until the run reaches the solver (needs the element-0 path handled too). Not blocking — 5a-1 is direct structural evidence.
- [x] **5a-3** Ran on 4 procs. The 5b gate ("proceed only if a worker shows missing sources") is **satisfied** — 29 missing on procs 1 and 3. Diagnostics still in place; remove after 5b lands.

### 5b. Fix: joint source-closure + topology fixpoint

- [x] **5b-1** DONE — replaced the node-only closure with a **source-closure-only** fixpoint that calls `select_sources` on every set node (+ duplicates), edge, face, and facet until the summed bitset population stops growing.
  - **Correction to the plan:** ~~run the geometric aura (node→element→edge/face, facet→facet) *inside* the fixpoint~~ — **wrong**: this grows the halo one element-layer per iteration → it walks the whole mesh → hang (observed). The geometric aura must stay a **single pass before** the fixpoint. A source entity only needs to *exist* as a basis (plus its sub-nodes, which `select_sources` adds), not its element neighborhood. **Codex + Grok confirmed** this for the CORC path; it is a documented general limitation (an off-element edge/face source discovered post-aura could lack DOF context — does not occur in CORC).
  - **Convergence** — sums node+edge+face+facet+**control-point** bitset counts (the CP term prevents early termination if a CP source is added); monotonic `set()` over a finite mesh guarantees termination.
- [x] **5b-2** DONE — deleted the dead pre-aura edge/face source loops.
- [ ] **5b-3** (deferred, latent — no CORC trigger) Closing the sources *of* elements/control points: the fixpoint does not call `select_sources` on elements or control points, and `select_sources` lacks a `CELL` case (`Element::entity_type()` returns `CELL`, not `ELEMENT`). Convergence already counts the CP bitset; what remains deferred is iterating element/CP carriers + the `CELL`/`Mesh::basis()` case. Safe to defer until a mesh actually carries element- or control-point-typed sources.

### 5c. Why a joint fixpoint rather than a single post-aura sweep

A single post-aura sweep is insufficient (both auditors): it iterates a fixed snapshot and misses sources introduced mid-pass — a source edge may itself be hanging, or may pull in a node whose element carries further hanging edges (the EDGE→EDGE cascade at `cl_FEM_DofMgr_DofData.cpp:3673-3706`). Source-closure and aura are mutually dependent, so they must interleave until the bitsets stabilize. Keeping the two phases distinct inside one loop, rather than blending them, is clearer and less error-prone.

## 6. Acceptance criteria

- [x] The 5a check reports **zero** missing sources for every proc (was 29 on procs 1/3). Now a debug-only `BELFEM_ERROR` invariant (`select_entities`), compiled out in release.
- [x] 4-proc CORC runs the full pipeline and many timesteps (58+), and **parallel tracks serial**: both converge cleanly without the buffer, and both stagnate identically *with* the buffer + ghost-stabilization config. So the distribution is validated — the residual stagnation is a **separate buffer/ghost-stab formulation issue** (serial-reproducible), not a distribution defect. (The original `-1`/`-9`/element-0 stops were a buffer double-count + MUMPS workspace, all resolved.)
- [x] Serial CORC is unchanged (the closure only affects the parallel `select_entities` path).
- [ ] No regression in an existing parallel Maxwell/thin-shell test — not formally run; covered by the tri-AI audit + the live multi-timestep CORC run.

## 7. Risks and mitigation

- **Over-ghosting.** The closure adds only entities genuinely referenced as sources, so growth is bounded by the real source graph and cannot pull unrelated mesh regions. Compare per-rank ghost counts before/after.
- **Termination / convergence metric.** Bits are only ever set, never cleared, so the loop terminates — *provided* convergence is measured across all bitsets, not node count alone (see 5b).
- **Aura ordering.** The topology/aura phase must run inside the fixpoint, after new sources are added, including facet-on-facet expansion; otherwise newly ghosted source nodes/facets miss their elements and neighbors.
- **Rollback.** The change is localized to `select_entities`; reverting restores prior behavior. No serialized-format change, so workers and master stay compatible.

## 8. Audit outcomes and remaining open items

Resolved by the tri-AI audit:

- [x] **Root cause (flag-order dead loop): verified high** by all three.
- [x] **Fix shape (joint fixpoint over all bitsets, aura inside the loop): endorsed** by Codex and Grok; a single sweep is insufficient.
- [x] **Convergence must sum all bitsets**, not node count (Codex + Grok).
- [x] **Facet-on-facet aura (`:464-467`) must be inside the fixpoint** (Codex + Grok).
- [x] **Fix the element/control-point gap now**, including the CELL/ELEMENT case (Codex), rather than deferring while claiming complete closure.
- [x] **Causality framing softened:** the defect is incomplete transitive/off-element closure, not "never ghosted"; the MUMPS link is plausible, not confirmed, until 5a runs.

Remaining open:

- [ ] Whether CORC actually exercises the control-point / element-source paths today, or whether a guard there suffices (5a will show).
- [ ] Whether proc 3's true `INFO(1)` is structural (confirming the matrix-inconsistency theory) or a resource code (which would re-scope the fix).
