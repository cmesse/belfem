# DR-194 and DR-197: Two One-Hunk Fixes in the Thin-Shell Setup Path

**Date:** 2026-10-08
**Purpose:** Close the two latent `[W]` rows left by the NI-coil campaign. DR-194:
`ThinShellFactory::update_node_facet_tables` walks the duplicates of `tNode` in two of its five
loops and of `tOrg` in the other three, which aliases today and overflows the facet container the
day it does not. DR-197: `CutFactory::~CutFactory` deletes the facets of each temporary thin-shell
sideset and then the sideset, whose own destructor deletes them again — a double free on any throw
between facet duplication and restoration. Mechanism: walk `tOrg` everywhere and guard the
container in setup code; delete only the sideset, which owns its facets.
**Module:** `src/fem/kernel` (`cl_ThinShellFactory`), `src/mesh` (`cl_Vertex.hpp`, one accessor), `src/homology` (`cl_CutFactory`, outside the closed core), `tests/fem`, `tests/homology`
**AIs involved:** Claude Fable 5.1 (plan, pre-registration, reconciliation); Codex gpt-6-astra/high + Grok grok-4.7/high (juries); Claude Opus 5.5 (implementation)
**Status:** PLAN v2 — 2026-10-08; Codex astra/high folded in (test seams, guard per destination, no strike before the valgrind gate); Grok pending; no source modified

> **Scope guards.** In scope: the five loops of `update_node_facet_tables`, one `facet_capacity()`
> accessor on `Vertex`, one `BELFEM_ERROR` in the factory, the `~CutFactory` loop, two tests.
> Out of scope: the thin-shell node tables' design, the ownership of `mesh->facets()` (a non-owning
> view: "facets are deleted by sideset", `cl_Mesh.cpp:~244`), anything in the closed cohomology core.

## 1. Current Behaviour and How It Fails

| Failure | Mechanism | Evidence |
|---|---|---|
| DR-194 heap overflow when `tNode != tNode->original()` | loops 1 (reset), 2 (count), 4 (fill) iterate `tOrg->duplicate( d )`; loop 3 (allocate) and loop 5 (unflag) iterate `tNode->duplicate( d )` (`cl_ThinShellFactory.cpp:2728-2805`). A duplicate counted in loop 2 but not allocated in loop 3 is filled in loop 4 past a zero capacity; `Vertex::add_facet` guards with `BELFEM_ASSERT` only (`cl_Vertex.cpp:250-254`), so release writes past the end | benign today: `aNodes` are `mMasterNodes`, asserted non-duplicates by `reset_node_indices` (`:701`), so the lists alias |
| DR-197 double free on unwind | `~CutFactory` deletes every facet of each temporary sideset, then `delete tSideSet`; `SideSet::~SideSet` deletes `mFacets` again (`cl_CutFactory.cpp:62-72`, `cl_SideSet.cpp:34-40`). On the normal path `restore_thin_shell_sidesets` moves the facets back and resets the container before deleting (`:463-491`), so only a throw between `duplicate_and_relink_facets` (`:1486`) and the restore reaches it | SEGV seen in gdb on an early draft of the `ThinShellAlignment` fixture, debug, throw mode |

**Bottom line:** two ownership slips in setup code, each one hunk, each with a throw-path or
alias-path test that is red on HEAD.

## 2. Ordered Steps

- [ ] **R1** `cl_Vertex.hpp`: inline `uint16_t facet_capacity() const` returning `mFacetCapacity` (the counter accessor `number_of_facets()` exists at `:182`).
- [ ] **R2** `update_node_facet_tables`: loops 3 and 5 iterate `tOrg->duplicate( d )`; in loop 4, immediately before **each** insertion — `tOrg->add_facet( tFacet )` and `tOrg->duplicate( d )->add_facet( tFacet )` — a `BELFEM_ERROR( <destination>->number_of_facets() < <destination>->facet_capacity(), … )` naming the destination node id (two guard sites, one per destination; setup code, always active; the container primitive keeps its `BELFEM_ASSERT`). Codex P3: no `tNode` is in scope in that loop.
- [ ] **R3** `~CutFactory`: delete the temporary sidesets only; the facet loop goes. Header comment on `mTemporaryThinShellSidesets` states that the sideset owns the facets it holds.
- [ ] **R4** tests. (a) `tests/fem`: the public path cannot produce a duplicate-valued `aNodes` — the constructor requires original master nodes (`cl_ThinShellFactory.hpp:363`), `create()` asserts it (`:701`) and `collect_nodes()` draws from `mMasterNodes` (`:852`) (Codex P2; O1 resolved). So the test executes the production implementation directly through a test seam: the loop body moves into a static free helper `update_node_facet_tables( Cell< Node * > &, Cell< Facet * > & )` in a header the factory and the test both include (or the member becomes public static; the implementer picks the smaller diff), and the test feeds it (i) original-only nodes and (ii) nodes that are themselves duplicates of originals with further duplicates, asserting exact facet membership and counts on every original and duplicate, and flag slot 2 cleared on all of them afterwards — the flag check is what catches an unfixed loop 5. (ii) is red on HEAD. (b) `tests/homology`: a `CutFactory` on a minimal thin-shell mesh, call the public `duplicate_and_relink_facets()` (`cl_CutFactory.hpp:158`), then let the factory go out of scope **without** `restore_thin_shell_sidesets`, then destroy the mesh (Codex P2: the closed-rim refusal runs before duplication, `:2096`/`:1146`, so it cannot reach this destructor state). Red on HEAD as a double free under valgrind; clean with the fix.
- [ ] **R5** gates: `make check` (17/17); the NI F deck of `dl20261008_seam_short_interface_duplicates.md` to its first save as a thin-shell smoke; **valgrind on test (b) and on the mesh destruction after it, in a non-MKL tree (`/opt/scls/debug` flavor, `test_homology` target only)** — this gate is part of the definition of done, not owed.
- [ ] **R6** strike DR-194 and DR-197 only when every R5 gate has run (Codex P2: the register forbids striking with residual work, `debt_register.md:51`); if the valgrind tree cannot be built tonight the rows stay open with the plan's state in their status cell; devlog.

## 3. Open Questions

- **O1** RESOLVED 2026-10-08 (Codex): it cannot; R4 (a) is a helper-level test of the production body.
- **O2** Whether any other owner holds pointers to the temporary sidesets' facets on the throw path (`mMesh->facets()` is a view; `mTemporaryThinShellSidesets( aID )->facets()` is handed out by `cl_CutFactory.hpp:236`). The jury should look for a consumer that deletes them.

## 4. Audit Trail

- Exchange thread: `tmp/ai_exchange/review_dr194_dr197.md`.
