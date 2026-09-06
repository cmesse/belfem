# Terminal Curve Orientation Ordering Fix (Tapestack Abort)

**Date:** 2026-07-02
**Purpose:** Root-cause and fix for `BELFEM_ASSERT( tFacet != nullptr, "No surface found on tape/boundary" )` in `CutFactory::orient_terminal_curves_sub()` on the 2-tape stack problem
**Module:** homology

## Symptom

`hphirun` on `tapestack3d.msh` (2 tapes, fresh Gmsh path) aborted at
`cl_CutFactory.cpp:2285` while the single-tape corc problem ran fine. A user
diagnostic confirmed the failing sideset (10) held 48 facets — not empty.

## Root Cause (confidence: high, confirmed by Codex and Grok audits)

Sideset 10 is the small inter-tape connector boundary at z=50 (not a domain
cap); its master volume is the thin slab sandwiched between the two tapes.
Inside `create_thin_shell_cuts()` the old call order was:

1. `duplicate_nodes_on_face_sidesets()` — duplicates all tape nodes
2. `relink_slave_elements_with_duplicate_nodes()` — slave-side volume elements
   get duplicate nodes
3. `relink_non_thinshell_facets()` — boundary facets whose master element was
   altered are rewritten to the master's (now duplicated) nodes
4. `orient_terminal_curves()` — flags the ORIGINAL nodes of `segments()(0)`
   and searches sideset facets for ≥2 flagged corners

When the inter-tape slab is the slave side of a tape, step 3 relinks the
connector-boundary facets to duplicate node pointers, so step 4's
original-node search finds no matching facet → abort. corc never hits this
because its terminal boundary belongs to the master-side air block.

## Rejected fix

First proposal (flag original + duplicate via `mThinShellDuplicates(index())`)
was refuted in audit: `append_move( mMesh->nodes(), mThinShellDuplicates )` at
line 1636 clears the source Cell (`cl_Cell.hpp:519`) *before* orientation ran,
so the lookup would index an empty container (debug assert / release OOB).

## Applied fix

Reordered `create_thin_shell_cuts()` (`cl_CutFactory.cpp:1613-1642`): the 3D
`orient_terminal_curves()` now runs immediately after
`create_curves_for_thinshells()`, before any node duplication or relinking.
At that point every facet and segment still holds original nodes, so the
existing search works unmodified. The 2D branch
(`orient_terminal_curves_2D()`) stays at its original position.

Safety checks done before the move:

- `orient_terminal_curves()` is self-contained (unflags all nodes on entry,
  cleans up its own flags); no dependency on the duplication steps.
- `Curve::reverse()` before `close_terminal_loops()` is coherent: the
  arclength fixup is monotonicity-agnostic, and the only `mS` consumers are
  CurveFactory, the fixup itself, and BfmFile serialization — the constant
  offset relative to the old ordering is inert.
- Do NOT switch to `original()`-identity matching: thin-shell duplicate nodes
  are created without `set_original()` (`cl_CutFactory.cpp:~1802`).

## Verification

Build and run handed to user: tapestack3d must pass the previous abort point;
corc must be re-run as regression. The temporary `#DIAG` cout in
`orient_terminal_curves_sub` is still in place and should be removed after a
successful run.

## AI collaboration

Exchange thread: `tmp/ai_exchange/terminal_curve_orient_dupnode.md`.
Codex (high confidence) confirmed the mechanism, refuted the first fix with
the `append_move` clear, and proposed the reorder. Grok (third voice)
independently confirmed the mechanism, the necessity of a `gNoIndex` guard in
the rejected variant, and rated residual wrong-facet risk low.
