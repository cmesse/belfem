# DR-194 and DR-197: Two One-Hunk Fixes in the Thin-Shell Setup Path

**Date:** 2026-10-08
**Purpose:** Close the two latent `[W]` rows left by the NI-coil campaign. DR-194:
`ThinShellFactory::update_node_facet_tables` walks the duplicates of `tNode` in two of its five
loops and of `tOrg` in the other three, which aliases today; the day it does not, a duplicate that
was reset and counted but never allocated receives `add_facet` on a null container (release:
a store through `nullptr`; debug: the assert throws). DR-197: `CutFactory::~CutFactory` deletes the facets of each temporary thin-shell
sideset and then the sideset, whose own destructor deletes them again — a double free on any throw
between facet duplication and restoration. Mechanism: walk `tOrg` everywhere and guard the
container in setup code; delete only the sideset, which owns its facets.
**Module:** `src/fem/kernel` (`cl_ThinShellFactory`), `src/mesh` (`cl_Vertex.hpp`, one accessor), `src/homology` (`cl_CutFactory`, outside the closed core), `tests/fem`, `tests/homology`
**AIs involved:** Claude Fable 5.1 (plan, pre-registration, reconciliation); Codex gpt-6-astra/high + Grok grok-4.7/high (juries); Claude Opus 5.5 (implementation)
**Status:** ✅ COMPLETE 2026-10-08 — landed on devel the same evening (Christian's autopilot ruling). R1 `Vertex::facet_capacity()`; R2 the `tOrg` walk in loops 3 and 5 plus a `BELFEM_ERROR` on each `add_facet` destination; R3 `~CutFactory` deletes only the owning sidesets; R3b `restore_thin_shell_sidesets` nulls each slot after delete; R4 (a) `test_ThinShellNodeFacetTables` (duplicate case throws on HEAD) and (b) `test_CutFactory_unwind` (SEGV / 4 valgrind contexts on HEAD, 0 with the fix); R5 gates: new gtests green, valgrind 0 errors in the `/opt/scls/debug` tree (suite + unwind), `make check` 17/17 zero warnings, NI F deck smoke unchanged; R6 DR-194 and DR-197 struck. Juries: plan astra/high + grok-4.7/high (every finding folded, no P0), code astra/high (approve) + grok-4.7/high (no defect; three text-level points applied). Residue: no test can provoke a throw inside the private restore loop (R3b read-verified only); O3 limit stands

> **Scope guards.** In scope: the five loops of `update_node_facet_tables`, one `facet_capacity()`
> accessor on `Vertex`, one `BELFEM_ERROR` in the factory, the `~CutFactory` loop, two tests.
> Out of scope: the thin-shell node tables' design, the ownership of `mesh->facets()` (a non-owning
> view: "facets are deleted by sideset", `cl_Mesh.cpp:~244`), anything in the closed cohomology core.

## 1. Current Behaviour and How It Fails

| Failure | Mechanism | Evidence |
|---|---|---|
| DR-194 null store when `tNode != tNode->original()` | loops 1 (reset), 2 (count), 4 (fill) iterate `tOrg->duplicate( d )`; loop 3 (allocate) and loop 5 (unflag) iterate `tNode->duplicate( d )` (`cl_ThinShellFactory.cpp:2728-2808`). A duplicate reset in loop 1 (`mFacets == nullptr`, capacity 0, `cl_Vertex.cpp:214-222`) and counted in loop 2 but not allocated in loop 3 receives `add_facet` in loop 4: `BELFEM_ASSERT` only (`cl_Vertex.cpp:250-254`), so release stores through null; the wrong unflag in loop 5 also leaves flag slot 2 set on the originals' duplicates | benign today: `aNodes` are `mMasterNodes`, asserted non-duplicates by `reset_node_indices` (`:701`), so the lists alias |
| DR-197 double free on unwind | `~CutFactory` deletes every facet of each temporary sideset, then `delete tSideSet`; `SideSet::~SideSet` deletes `mFacets` again (`cl_CutFactory.cpp:62-72`, `cl_SideSet.cpp:34-40`). On the normal path `restore_thin_shell_sidesets` moves the facets back and resets the container before deleting (`:463-491`), so only a throw between `duplicate_and_relink_facets` (`:1486`) and the restore reaches it | SEGV seen in gdb on an early draft of the `ThinShellAlignment` fixture, debug, throw mode |

**Bottom line:** two ownership slips in setup code — DR-194 is two loop edits plus one accessor,
DR-197 one destructor site plus a two-line hardening of the restore — each to be pinned by a test
that is red on HEAD (none exists yet).

## 2. Ordered Steps

- [x] **R1** `cl_Vertex.hpp`: inline `uint16_t facet_capacity() const` returning `mFacetCapacity` (the counter accessor `number_of_facets()` exists at `:182`).
- [x] **R2** `update_node_facet_tables`: loops 3 and 5 iterate `tOrg->duplicate( d )`; in loop 4, immediately before **each** insertion — `tOrg->add_facet( tFacet )` and `tOrg->duplicate( d )->add_facet( tFacet )` — a `BELFEM_ERROR( <destination>->number_of_facets() < <destination>->facet_capacity(), … )` naming the destination node id (two guard sites, one per destination; setup code, always active; the container primitive keeps its `BELFEM_ASSERT`). Codex P3: no `tNode` is in scope in that loop.
- [x] **R3** `~CutFactory`: delete the temporary sidesets only; the facet loop goes. Header comment on `mTemporaryThinShellSidesets` states that the sideset owns the facets it holds (a contract comment, no ids or dates).
- [x] **R3b** (Grok) `restore_thin_shell_sidesets`: set the cell slot to `nullptr` right after each `delete tTemp`, and let `~CutFactory` skip null slots, so a throw later in that loop (the `new` at `:469`, the `sideset()` lookup) cannot leave a freed pointer for the destructor.
- [x] **R4** tests. (a) `tests/fem`: the public path cannot produce a duplicate-valued `aNodes` — the constructor requires original master nodes (`cl_ThinShellFactory.hpp:363`), `create()` asserts it (`:701`) and `collect_nodes()` draws from `mMasterNodes` (`:852`) (Codex P2; O1 resolved). So the test executes the production implementation directly through a test seam: the loop body moves into a static free helper `update_node_facet_tables( Cell< Node * > &, Cell< Facet * > & )` in a header the factory and the test both include (or the member becomes public static; the implementer picks the smaller diff), and the test feeds it (i) original-only nodes and (ii) nodes that are themselves duplicates of originals with further duplicates, asserting exact facet membership and counts on every original and duplicate, and flag slot 2 cleared on all of them afterwards — the flag check is what catches an unfixed loop 5. (ii) is red on HEAD by the thrown `std::runtime_error` from the `add_facet` assert (assertion build throws, `assert.cpp:310`); with the fix (ii) completes and the expectations hold. (b) `tests/homology`: a `CutFactory` on a minimal thin-shell mesh (a sideset of one or two facets between two blocks, a terminal curve as the existing fixtures build them), call the public `create_thin_shell_cuts()` (or `duplicate_and_relink_facets()`, `cl_CutFactory.hpp:147-158`) and then let the factory go out of scope **without** `run()`, then destroy the mesh (Grok: the window opens at the first push `:1495`; `restore` is private and only `run()` reaches it; the closed-rim refusals run before duplication and never enter the window). Red on HEAD as an invalid free under valgrind (no throw needed); clean with the fix. The test's contract is a clean destroy of factory then mesh, not a continued solve.
- [x] **R5** gates, in this order of weight: (1) the new gtests run and green, with their HEAD-red controls recorded verbatim (the ctest rollup `make check` 17/17 is a non-regression rollup, not a gtest count — Grok); (2) **valgrind on `test_homology` in the non-MKL `/opt/scls/debug` tree** (built: `wt2-build-nomkl`), the (b) test and the mesh destruction after it clean, against the HEAD baseline of the same binary; (3) `make check` in the MKL debug tree; (4) the NI F deck to its first save as a thin-shell non-regression smoke — HEAD passes it too, so it is not closure evidence.
- [x] **R6** strike DR-194 and DR-197 only when every R5 gate has run (Codex P2: the register forbids striking with residual work, `debt_register.md:51`); if the valgrind tree cannot be built tonight the rows stay open with the plan's state in their status cell; devlog.

## 3. Open Questions

- **O1** RESOLVED 2026-10-08 (Codex): it cannot; R4 (a) is a helper-level test of the production body.
- **O2** RESOLVED 2026-10-08 (Grok): no second owner — `thin_shell_facets` (`cl_CutFactory.hpp:233-237`) has no callers; `mFacets` and `mFacetMap` are non-owning (`cl_Mesh.hpp:87`, `cl_Mesh.cpp:243`). After R3 the one `delete tSideSet` is the single free.
- **O3** (Grok, limit, not fixed) after a throw inside the window the mesh is not usable: the live sideset holds tA/tB and `mFacets` may dangle until `~Mesh`, which does not delete through it. Out of scope: the contract is a clean destroy.

## 4. Audit Trail

- Exchange thread: `tmp/ai_exchange/review_dr194_dr197.md`.
