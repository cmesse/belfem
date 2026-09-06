# Thin-Shell Break Points After PENTA/HEX Facet Renumbering

**Date:** 2026-04-22
**Purpose:** Trace the thin-shell model (PENTA6TS / HEX8TS) and `mt_maxwell_h.cpp` to find where the new PENTA/HEX edge/facet numbering and outward-normal convention have broken the H-φ coupling.
**Module:** fem/maxwell, mesh
**Audience:** BELFEM developers working on the thin-shell or Nédélec basis.
**Mode:** Read-only trace. No source edits made.

---

## 1. What the user changed

Confirmed by diff against HEAD (`git diff HEAD`).

### 1.1 Solid hexes and prisms — edge slot re-ordering

HEX8, HEX20, HEX27, HEX64 and PENTA6, PENTA15, PENTA18 had their edge slot numbering rearranged. For HEX8 (`src/mesh/cl_Element_HEX8.hpp:125-200`):

| Slot | Old `get_nodes_of_edge` | New `get_nodes_of_edge` |
|------|---|---|
| 0-3 | bottom face, **reversed** direction (1→0, 2→1, 3→2, 0→3) | bottom face, forward (0→1, 1→2, 2→3, 3→0) |
| 4-7 | verticals (0→4, 1→5, 2→6, 3→7) | top face (4→5, 5→6, 6→7, 7→4) |
| 8-11 | top face (4→5, 5→6, 6→7, 7→4) | verticals (0→4, 1→5, 2→6, 3→7) |

For PENTA6 (`src/mesh/cl_Element_PENTA6.hpp:116-175`):

| Slot | Old | New |
|------|---|---|
| 0-2 | bottom, reversed directions (1→0, 2→1, 0→2) | bottom, forward (0→1, 1→2, 2→0) |
| 3-5 | verticals (0→3, 1→4, 2→5) | top (3→4, 4→5, 5→3) |
| 6-8 | top (3→4, 4→5, 5→3) | verticals (0→3, 1→4, 2→5) |

`get_edges_of_facet` was updated in lockstep for each solid element, so the mapping slot↔facet is internally consistent again.

### 1.2 Thin-shell elements — facet count and normal convention

`PENTA6TS` template parameters changed from `<6,6,6,2,2>` to `<6,6,6,5,1>` (`src/mesh/cl_Element_PENTA6TS.hpp:31-70`). Facet count went 2 → 5, matching volume PENTA6:

| New facet | Nodes | Shape | Geometric role |
|---|---|---|---|
| 0 | {0,1,4,3} | quad | side |
| 1 | {1,2,5,4} | quad | side |
| 2 | {0,3,5,2} | quad | side |
| **3** | **{0,2,1}** | **tri** | **bottom, normal points DOWN (outward)** |
| **4** | **{3,4,5}** | **tri** | **top, normal points UP (outward)** |

Old facet numbering (pre-commit):

| Old facet | Nodes | Normal |
|---|---|---|
| 0 | {0,1,2} | UP (CCW seen from above) |
| 1 | {3,4,5} | UP (CCW seen from above) |

That is the "parallel normals" convention the user explicitly called out. Bottom facet is no longer at index 0 — it is at index 3, and its winding has been reversed so the outward normal flips direction.

`HEX8TS` went `<8,8,8,2,1>` → `<8,8,8,6,1>` (`src/mesh/cl_Element_HEX8TS.hpp:24-165`). Now:

| New facet | Nodes | Geometric role |
|---|---|---|
| 0-3 | side quads {0,1,5,4}, {1,2,6,5}, {2,3,7,6}, {0,4,7,3} | lateral (side) |
| 4 | {0,3,2,1} | bottom, normal DOWN |
| 5 | {4,5,6,7} | top, normal UP |

Old HEX8TS had facet 0={0,1,2,3} (bottom, normal UP) and facet 1={4,5,6,7} (top, normal UP) — same-direction normals.

### 1.3 Connector element type swap

`ThinShellFactory::create_side_elements` now builds connectors as `ElementType::HEX8TS` (8 edges, 6 facets) instead of solid `ElementType::HEX8` (12 edges, 6 facets). `create_vertical_edges` is commented out. Edge slots 0-3 receive the bottom-layer horizontal + longitudinal edges, slots 4-7 receive the top-layer equivalents (`src/mesh/cl_ThinShellFactory.cpp:2545-2621`).

---

## 2. Break points in the consumer code

Confidence: **high** on everything in this section — the lines are grepped and the logic is re-derived from current source.

### 2.1 MaxwellFactory — four hardcoded `0`/`1` facet indices on PENTA6TS

`src/fem/maxwell/cl_MaxwellFactory.cpp`:

| Function | Line | Current call | Intent |
|---|---|---|---|
| `hang_thinshell_edges_on_nodes_bottom` | 1367, 1380 | `aElement->get_{nodes,edges}_of_facet(0, …)` | "PENTA6TS bottom triangle" |
| `hang_thinshell_edges_on_nodes_top` | 1438, 1447 | `…get_{nodes,edges}_of_facet(1, …)` | "PENTA6TS top triangle" |
| `hang_thinshell_edges_on_edges_bottom` | 1494, 1510 | `…get_{nodes,edges}_of_facet(0, …)` | "PENTA6TS bottom triangle" |
| `hang_thinshell_edges_on_edges_top` | 1601, 1614 | `…get_{nodes,edges}_of_facet(1, …)` | "PENTA6TS top triangle" |

Three of those sites carry a comment stating explicitly:
> "always face 0 because aElement is a PENTA6TS whose bottom face is at index 0 by construction"
> "always face 1 because aElement is a PENTA6TS whose top face is at index 1 by construction"

Both comments are now stale.

With the new PENTA6TS layout:
- `get_nodes_of_facet(0, …)` returns the **4-node side quad** {0,1,4,3} instead of the 3-node bottom triangle {0,1,2}. The caller then iterates `aFacet->number_of_nodes()` (3 for a TRI3 mesh facet) into that container, so three of the four stored nodes are *some* three of the side-quad nodes — not the intended bottom-triangle nodes. The hanging-source wiring is therefore silently wrong, not crashing.
- `get_edges_of_facet(0, …)` on the new PENTA6TS runs through the `default` branch of the switch (`cl_Element_PENTA6TS.hpp:192-211` defines only cases 3 and 4), which calls `throw_facet_error(0)` and fails loudly via `BELFEM_ERROR`. This is almost certainly the first visible failure mode when running anything that reaches `MaxwellFactory::create_hanging_edges_and_facets()`.

Fix needed (not applied): replace `0` → `3` and `1` → `4` at all four sites; update the three stale comments.

### 2.2 ThinShellFactory::create_ghost_facets — same issue on PENTA6TS

`src/mesh/cl_ThinShellFactory.cpp:1839-1840`:

```cpp
tGhost->set_master( tMasterBlock->element( f ), 1 );
tGhost->set_slave( tSlaveBlock->element( f ), 0, 1 );
```

Old reading: the lower block's TOP triangle is master-facet 1, the upper block's BOTTOM triangle is slave-facet 0, both with the old same-direction-normal convention.

New reading: master-facet 1 is a side quad `{1,2,5,4}`, slave-facet 0 is a side quad `{0,1,4,3}`. When `Mesh::update_facet_nodes()` (`src/mesh/cl_Mesh.cpp:683-702`) later calls `tFacet->master()->get_nodes_of_facet(1, tNodes)` on a PENTA6TS master it copies 4 nodes into a ghost TRI3 element that has room for 3 — silently overwriting. The ghost facet ends up wired to four lateral-quad nodes instead of three top-triangle nodes.

Fix needed (not applied): `set_master(…, 4)` and `set_slave(…, 3, …)` in the new convention.

### 2.3 Calculator::normal_penta_ts — wrong cases + missing sign flip

`src/fem/kernel/cl_FEM_Calculator.cpp:1857-1912`. The switch handles only `mMasterIndex ∈ {0, 1}` and returns the **same** formula `(dr/dξ) × (dr/dη)` for both. That was consistent with the old thin-shell convention (both facets' normals parallel and pointing up in ref space).

With the new PENTA6TS:

- The master index reaching this function is now 3 (bottom) or 4 (top), not 0 or 1 — so the switch falls through to `default` and fires `BELFEM_ERROR("Invalid master index for facet")` as soon as `compute_bn` is called. This will be the crash site for any h-φ thin-shell kernel (`h_ts_metal`, `h_ts_hts`, all the `h_ts_*` variants in `mt_maxwell_h.cpp`).
- Even after cases 3 and 4 are added, the formulas have to differ in sign. The outward bottom-triangle normal is `-(dr/dξ × dr/dη)`; the outward top-triangle normal is `+(dr/dξ × dr/dη)`. Compare `normal_penta` (volume) cases 3 and 4 at `cl_FEM_Calculator.cpp:1816-1845`, which already implement exactly that sign pattern.
- The stale comments on the existing cases are misleading:
  - `case 0 : // corresponds to negative normal of facet 4 of volume penta` — formula is NOT negated; appears to be a note-to-self from an earlier rewrite that was never acted on.
  - `case 1: // corresponds to normal of facet 5 of volume penta` — volume PENTA has no facet 5.

Fix needed (not applied): add cases 3 (negate the cross product) and 4 (keep formula) to mirror `normal_penta`.

### 2.4 `compute_bn` semantics — check sign convention downstream

`src/fem/maxwell/matrices/mt_maxwell_h.hpp:130-173`. The helper computes

```
bn = 0.5 · (bm + bs);       bn ← (bn · n) · n
```

where `bm` and `bs` are the φ-side B-fields on the two volume neighbors and `n = tCalc->normal(k)`. The old code implicitly assumed `n` points from slave to master (or vice versa) consistently across master/slave facets because both PENTA6TS facet normals were parallel.

In the new convention:
- For the thin-shell sideset facet whose master is a volume PENTA on the UP side, `mMasterIndex = 3` (volume PENTA bottom tri {0,2,1}, outward is DOWN, i.e. pointing toward the shell from above).
- For the thin-shell sideset facet whose master is the volume PENTA on the DOWN side, `mMasterIndex = 4` (top tri, outward is UP).

If `compute_bn`'s callers assumed a fixed sign for `n` (e.g. "always from master side toward slave side"), the new outward-on-each-side convention will introduce a sign flip whenever master/slave are on opposite sides from the old assumption. The code path (see `h_ts_metal` at `mt_maxwell_h.cpp:91-161`, `h_ts_hts` at `:467-552` and the `_t` / piecewise variants) multiplies by `n` twice so `bn` itself is sign-invariant in magnitude, but `tape_normal = n` is used directly to compute β = acos(H · n / |H|); a flipped `n` flips β across the tape, which changes `Ic(B, θ)` for all piecewise HTS paths.

Fix needed (not applied): once `normal_penta_ts` is corrected, audit the thin-shell callers to confirm whether `n` still satisfies the orientation contract they rely on. If the contract is "always from master to slave", then after adding cases 3/4 the `n` for new-master-index-3 will point from master INTO the shell, which is opposite to what the slave-index-4 case does. That's the concrete downstream check the user is worried about.

### 2.5 HEX8TS side facets have no `get_edges_of_facet` entry

`src/mesh/cl_Element_HEX8TS.hpp:181-205`: only cases 4 (bottom) and 5 (top) are implemented; cases 0-3 fall through to `throw_facet_error`. This is fine as long as only the bottom / top faces are queried, but any future h-φ coupling through the **lateral** HEX8TS face (the shell-connector fold described in `src/fem/maxwell/doc/shell_connector_coupling.md`) would hit this.

Not a current break, but worth a note — the comment on side-connector outer-face BCs in `shell_connector_coupling.md §4.2` may need to cross-reference this once someone adds the fold coupling.

### 2.6 Typo in commented-out connector code

`src/mesh/cl_ThinShellFactory.cpp:2609`:

```cpp
//Element->insert_edge( tVerticalEdgesOuter( k+1 ), 7 );
```

Missing `t` — should be `//tElement->insert_edge(...)`. Harmless because the line is commented out, but it will break the moment anyone un-comments the block.

---

## 3. Why is the error sometimes silent vs loud?

Given the four classes of defect above, three failure signatures are possible depending on where execution reaches first:

| Signature | Source |
|---|---|
| `BELFEM_ERROR: "Invalid master index for facet"` | §2.3, first call to `compute_bn` after the new master indices 3 or 4 arrive (any `h_ts_*` kernel). |
| `throw_facet_error(0)` from `cl_Element_PENTA6TS.hpp` | §2.1, as soon as `MaxwellFactory::create_hanging_edges_and_facets` runs `get_edges_of_facet(0, …)` on a PENTA6TS. |
| Wrong hanging-edge sources, garbled-looking fields | §2.1 / §2.2, if node-only paths run first (`get_nodes_of_facet` writes 4 nodes where only 3 are read). |

A silent failure of the third form would be consistent with "equations in `mt_maxwell_h.cpp` look wrong" as opposed to an outright crash.

---

## 4. Is the root cause the normals, face numbering, or edge directions?

Confidence: **high**.

**Facet numbering is the primary driver.** The new `5`-facet PENTA6TS pushes its top/bottom triangles from indices 0/1 to 3/4. Four sites in `MaxwellFactory` and one in `ThinShellFactory` still hardcode the old indices. `normal_penta_ts` additionally needs new cases 3/4 with a sign flip on case 3.

**The normal-convention change is a smaller, dependent effect.** Once the facet index is correct, the bottom facet's outward normal naturally flips (because its winding was reversed). `normal_penta_ts` currently has the "both normals up" formula baked in; one-half of the cases needs a minus sign. `compute_bn`'s callers must also be rechecked for the "master→slave" vs "always outward" assumption.

**Edge direction changes are not the culprit for the thin-shell break.** The only edge-direction change was on the bottom faces of HEX8/PENTA6 (reverse → forward) and the reshuffle of top-vs-vertical slots for HEX8/PENTA6. `ThinShellFactory::link_elements_with_edges` at `cl_ThinShellFactory.cpp:1573-1636` indexes by `tFacet->element()->edge(e)->index()`, and `Element::compute_edge_directions()` (via per-edge node-ID comparison) resolves local-vs-stored orientation. So the Nédélec edge basis's scalar signs flip in the local frame automatically, and the new convention on solid hexes/prisms is internally consistent. `cl_EF_HEX8.cpp` was updated in the same commit and the prior CODEX audits at `todo/ai_exchange.md:46-109` confirm its twelve reference-edge circulations now sum to `+1`.

---

## 5. Minimal repair list (for the next session)

Listed in dependency order so that each stage compiles. Not applied here — read-only trace.

1. `src/mesh/cl_ThinShellFactory.cpp:1839-1840` — change `set_master(…, 1)` → `set_master(…, 4)` and `set_slave(…, 0, 1)` → `set_slave(…, 3, 1)`.
2. `src/fem/maxwell/cl_MaxwellFactory.cpp` lines 1367, 1380, 1438, 1447, 1494, 1510, 1601, 1614 — replace `0` → `3` and `1` → `4` in the four `hang_thinshell_edges_on_{nodes,edges}_{top,bottom}` helpers, and refresh the three stale "face 0/1 by construction" comments.
3. `src/fem/kernel/cl_FEM_Calculator.cpp:1857-1912` — add `case 3` (bottom, negated cross product, triangular weight) and `case 4` (top, cross product, triangular weight) to `normal_penta_ts`, matching the `normal_penta` pattern at `:1816-1845`. Delete the stale `// facet 4 / facet 5 of volume penta` comments. Keep or remove cases 0/1 depending on whether any old call-site still feeds them.
4. Re-run a thin-shell patch test; confirm the `h_ts_*` kernels see a consistent `n` (§2.4 contract check).
5. Low priority: `src/mesh/cl_ThinShellFactory.cpp:2609` — fix `//Element->` → `//tElement->` typo in the commented-out `aSign < 0` HEX8TS block.

## 6. Open questions for the user

1. Was the new outward-normal convention meant to apply to both `normal_penta_ts` cases, or only the bottom triangle? The `normal_penta` (volume) pattern has opposite signs on bottom/top; `normal_penta_ts` currently has identical formulas. I've assumed the former was the intent.
2. Is the plan to keep `mMasterIndex ∈ {0, 1}` as a legacy alias for any code path, or can those cases be removed outright?
3. Does `compute_bn` need to encode "which side is master" explicitly (e.g. store an explicit sign), or will the two sign flips (facet normal flip + master choice flip) cancel? This is §2.4 and I think it's worth a paper-and-pencil check before running.

---

## References

- Prior Codex audits on the solid-hex edge direction / EF_HEX8 fixes: `todo/ai_exchange.md:46-183` (2026-04-20 entries).
- `src/fem/maxwell/doc/shell_connector_coupling.md` — background on HEX8TS side-connector design.
- Messe et al. 2023 (paper1) §2 — thin-shell formulation overview and static-condensation interface design.
