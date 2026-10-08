# Thin-shell sideset alignment after the same-type walk

**Date:** 2026-10-03
**Purpose:** record the fix for mixed-winding thin-shell sidesets in 3D and the updated orientation checker
**Module:** fem/maxwell

## Trigger

NI pancake test deck (`NI_coi/demo_clean_dummy_contact_test`): a 4.5-turn REBCO spiral as one
thin shell (14 sidesets), copper between the turns and a copper liner inside and outside, so
every tape sideset is conductor|conductor. Stock 0.9.2 (release build) aborted in setup:

```
SEGV / "corrupted double-linked list"
  mesh::Vertex::add_facet
  mesh::ThinShellFactory::update_node_facet_tables
  mesh::ThinShellFactory::create
  fem::MaxwellFactory::create_thinshells
```

## Diagnosis

- `toolbox/check_shell_orientation.py` (stock model) reported sideset 130 with 14 of its 300
  facets mastered by the other block, and sideset 177 with 4 of 228. All of them sit at the
  180° seam.
- Tracing which facet the walk came from showed two sources:
  - the liner's split face on the 180° plane (sidesets 209/210, copper|copper, perpendicular
    to the tape), seeded from the liner|air faces 201/191;
  - near 0°, the coil's top and bottom faces (copper|air), passing through sideset 124.
- Where a perpendicular face meets the tape along an edge shared by three or more facets,
  `check_facet_orientation` compares facets that are not neighbours on one surface. So one
  sideset is reached from two directions with opposite windings.
- No deck sign can repair this, because a sign flips the whole sideset.
- Removing the copper|copper split faces (simulated in the checker) left 6 bad edges in 130,
  so no mesh change fixes it reliably.
- Crash mechanism: with a mixed master side, boundary nodes are duplicated on both sides.
  `update_node_facet_tables` sizes the facet arrays over `tNode`'s duplicates
  (`cl_ThinShellFactory.cpp:2770`) but fills them over `tOrg`'s (`:2754`, `:2788`). The guards
  in `Vertex::add_facet` are `BELFEM_ASSERT` only, so a release build corrupts the heap.
  Confidence medium: this follows from reading the code, not from a debugger trace.

## Change

`MaxwellFactory::align_thin_shell_facets()`, called at the end of `fix_facet_masters` in 3D
while the facet-to-facet connectivity still exists:

- For every sideset of every `Protoshell`, a breadth-first pass over the edge neighbours that
  belong to the same sideset records each facet's winding relative to the seed of its patch.
  Facet flag slots 1 to 3 hold this state; they are cleared afterwards.
- **Vote:** the cross-type facets decide when the patch has any, otherwise all facets do. A tie
  keeps the seed's winding. The minority is flipped with `Facet::flip()`.
- **Hard error:** a facet to be flipped has no slave, or the patch is not orientable
  (a parity conflict, for example a non-manifold or self-crossing patch).
- **Messages:** facets aligned per sideset; a sideset that splits into several patches; a
  patch whose cross-type facets disagree.
- **Unchanged:** uniform sidesets (so existing deck signs keep their meaning), the walk itself,
  and 2D, which never runs the walk.
- Why this option is the most robust: it is independent of element numbering and of the order
  the walk visits facets. Leaving shell facets out of the walk would have fallen back on the
  lower-element-ID rule, which is uniform only when gmsh numbers each volume as a block.

## Files changed

| File | Change |
|---|---|
| `src/fem/maxwell/cl_MaxwellFactory.hpp` | declaration of `align_thin_shell_facets()` |
| `src/fem/maxwell/cl_MaxwellFactory.cpp` | call at the end of `fix_facet_masters` (3D); new method after it |
| `toolbox/check_shell_orientation.py` | models the new pass (step 2b) by default; `--legacy` models stock 0.9.2; reports per sideset what the pass aligns |
| `CHANGELOG.md` | `[Unreleased]` → Fixed entry |
| `devlog/dl20261003_thin_shell_alignment.md`, `devlog/README.md` | this record |

The same script was also copied to `belfem-0.9.2/toolbox/` (outside this tree). With
`--legacy` it gives byte-identical output to the previous version on the toolbox deck.

## Evidence (evidence ladder)

- **Static:** `cl_MaxwellFactory.cpp` passes `-fsyntax-only` under the module's own
  `compile_commands.json` flags (`-Werror`, MPI compiler from `/opt/scls/mkl`), rc 0.
- **Checker:**
  - test deck, patched model: 130 and 177 aligned (14 and 4 facets), every seam relation
    clean; suggested line `-120, -130, 124, 133, 139, 147, 153, 161, 167, 177, -119, 182, 175, -103`,
    the same as with `--radial 0,0,0,0,0,1` (all normals outward);
  - original toolbox deck: "already uniform", output otherwise identical to the stock model.
- **Not run:** no rebuild, no belfem run of the patched source, no `make check`. The fix is
  **reviewed, not verified**. The gate is a rebuild followed by the test deck: setup must pass
  `create_thinshells` and print the two "aligned" lines.

## Open

- `update_node_facet_tables` loops over the duplicates of `tNode` in one place and of `tOrg`
  in another (`:2770`, `:2802`). That inconsistency is still in the code; this change only
  removes the trigger the deck hit. Worth unifying on `tOrg`.
- `Vertex::add_facet` overflow is guarded only by `BELFEM_ASSERT`; a `BELFEM_ERROR` capacity
  check on the setup path would turn the next case into a message.
- Whether alignment cures the 180° current jump seen in `demo_clean_dummy_contact` (the
  air-facing turns) is a separate question; the deck with liners is the test for it.
