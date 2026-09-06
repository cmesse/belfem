# Ferro-Air Interface Removal and InterfaceCondFerro Deactivation

**Date:** 2026-08-12
**Purpose:** Retire the dead ferro-air interface kernel and give `InterfaceCondFerro`
the same treatment its siblings already have, closing R9 of
`todo/ferro_undulator_mu_and_interface_bugs.md`.
**Module:** `src/fem/maxwell`

## Why

Tracing D3 (`InterfaceCondFerro` is admitted to the magnetic equation but has no case in
`IWG_Maxwell::link_to_group`, so it aborts `"Not implemented"`) surfaced a second, larger
finding: **there is no ferro-air interface weak form, in 2d or in 3d.** Christian confirmed
the design intent — the interface nodes are duplicated for visualization only.

The tree agrees, twice over:

1. `phi_phi_2d` printed its own name, dumped the dofs and called **`exit(0)`**. Everything
   below that was the bubble-enrichment body behind `#ifdef BELFEM_ENRICH_PHIA`, a macro
   defined nowhere in the tree. `phi_phi_3d` was `BELFEM_ERROR( false, "not implemented" )`.
2. `InterfaceFerroAir` is set `Inactive` under h-phi in `create_hanging_edges_and_facets`,
   so the group never links and the stub was never dispatched.

That kernel is residue of a dropped experiment: enriching the phi elements of the iron.
It was abandoned, and the enrichment space was the wrong one — bubble functions rather
than hierarchical (Dular et al. 2021 for why identical FE orders at an interface need
hierarchical enrichment). `cl_MaxwellFactory.cpp` already recorded this in a comment:
*"the bubble functions in BELFEM are not the right space for it"*. `mUseEnrichment` holds
its `false` initializer and has no writer.

### Unreachability, established rather than assumed

`phi_phi_2d` / `phi_phi_3d` could not be reached in any configuration:

- **Under h-phi** — the sideset is `Inactive` before the dof manager builds its groups, and
  both live sideset loops skip inactive groups (`cl_FEM_DofManager.cpp:688`, `:750`).
- **Under any other formulation** — the `else` backstop in the `InterfaceFerroAir` case
  raised before dispatch.
- **The one unguarded `link_to_group( tSideSet )`** (`cl_FEM_DofManager.cpp:510`) sits behind
  `if ( mIWG->compute_jacobian_on_sideset() )`, and `mComputeJacobianOnSideset` is `false`
  at `cl_IWG.hpp:103` **with no writer anywhere in the tree**.

So the `exit(0)` could not fire today. It was one bypassed deactivation away from firing,
with no MPI finalize and an exit code of 0 — a driver script would have read it as success.

## What was done

1. **`cl_IWG_Maxwell.cpp`** — the `InterfaceFerroAir` case collapses to a single
   `BELFEM_ERROR( false, "Ferro-Air Interfaces must be disabled!" )`, matching
   `InterfaceCondAir`. The `if ( HPhi ) / else` structure went with it: the `else` was a
   backstop added earlier the same day, and the branch it guarded turned out to be the dead
   half. The `mt_maxwell_phi_phi.hpp` include is removed.
2. **`InterfaceCondFerro` gains a case** beside `InterfaceCondAir`, with the same
   "must be disabled!" error, so a bypassed deactivation fails with the rule rather than
   with `"Not implemented"`.
3. **`matrices/mt_maxwell_phi_phi.{cpp,hpp}` deleted**, with the `CMakeLists.txt` entry.
   No `phi_phi` reference remains in `src/`.
4. **`cl_MaxwellFactory.cpp`, `create_hanging_edges_and_facets`** — `InterfaceCondFerro` is
   set `Inactive` once its hanging-edge condensation is done, mirroring `InterfaceFerroAir`.
   H-phi is implied by the enclosing branch condition. **This is the D3 fix (R9).**
5. **`src/fem/kernel/doc/dof_manager_usage_guide.md`** — two stale references to
   `phi_phi_3d` as live φ-φ coupling corrected, and the "Disabled Interfaces" list now names
   `InterfaceCondFerro` and `InterfaceFerroAir` with the visualization rationale.

### Why deactivate rather than de-admit

The first proposal was to drop `InterfaceCondFerro` from the whitelist in
`set_block_types_in_magnetic_equation`, making it match `InterfaceCondAir`. Reading the tree
changed the recommendation: the whitelist runs *before* the condensation pass, and the fem
groups take the **live** mesh type (inherited at `cl_FEM_SideSet.cpp:47-50`, re-stamped at
`cl_MaxwellFactory.cpp:736-740`), which is exactly how `InterfaceFerroAir` and `ThinShell`
are already made dormant. Deactivating follows the established idiom and does not depend on
reasoning about admission ordering.

| type | whitelisted | deactivated after condensation | net |
|---|---|---|---|
| `InterfaceCondAir` | no | no | never a group |
| `InterfaceFerroAir` | yes | yes (h-phi) | dormant |
| `ThinShell` | yes | yes | handled separately |
| `InterfaceCondFerro` | yes | **now yes** (h-phi) | dormant — was an abort |

`Maxwell_FieldList` needs no change: it has cases for both the pre- and post-deactivation
type, so either snapshot resolves.

## Verification

`g++ -fsyntax-only` on both edited translation units with the exact
`CXX_FLAGS`/`CXX_DEFINES`/`CXX_INCLUDES` from `cmake-build-debug`. **Reviewed, not verified:**
no build, no run. The gate is a deck with a conductor block touching a ferro block — which
could not run at all before this change, since it aborted at link time.

## Open

- Nothing from this work. **Note for anyone reading the two 2026-08-12 devlogs together:** the
  `tPropagate = ( dim == 3 )` guard that appeared in `cl_MaxwellFactory.cpp` while this session
  was editing is **R13** of the parallel ferro-undulator session, not a competing fix — it is
  Christian's option B, taken after the R8 run cleared D4's abort and hit D5. R10's chain rule
  survives for 3D; 2D no longer propagates at all. See
  `devlog/dl20260812_ferro_undulator_d1_fix.md`.
- **`mComputeJacobianOnSideset` is `false` with no writer**, so the entire sideset branch of
  `DofManager::compute_jacobian` (`cl_FEM_DofManager.cpp:496` and its element loop) is dead
  for every IWG in the framework. Not touched here — bigger blast radius, deserves its own look.
- **The bubble machinery remains**: `src/fem/interpolation/bubble/` (14 files), the
  `mEnrichmentData` paths in `cl_FEM_Block.cpp:145` and `cl_FEM_SideSet.cpp:503`,
  `cl_IWG_Maxwell.cpp:518` and `cl_CutFactory.cpp:722`, all gated on a permanently-false flag.
  Retire, or keep as scaffolding for a future hierarchical enrichment.
