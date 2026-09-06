# 2D h-φ Thin-Shell Pipeline — Todo Checklist

**Date:** 2026-07-01 (v1) / 2026-07-10 (v2 — second iteration)
**Purpose:** Dependency-ordered work items to bring the 2D thin-shell path up to the 3D
reference. Derived from `todo/2d_thinshell_gap_analysis.md` (file:line evidence and
confidence labels live there; v2 line references re-baselined to HEAD of 2026-07-10).
**Status:** EXECUTION STARTED 2026-08-04/05, **PAUSED — the sprint calendar below has
expired and should be read as an ordering, not a schedule** (corrected 2026-08-11).
Sprint goal was first-order (linear) 2D thin shells within two weeks from 2026-08-04; that
window closed, and the day-by-day plan ("Days 1-2", "Day 13") no longer describes anything.
Order-2 items remain out of scope.

**2026-08-11 — Christian reports the 2D thin shells are WORKING**, on a deck that is not in
the repository. That most likely settles B8 (its fix was applied 2026-07-10 and marked
"pending build+run verification"), E4, C3 (explicitly "M1-contingent: fix on first failure"),
and M1 — but none of that is recorded here, because nobody else can reproduce the run.
**T0 is therefore now the single highest-value item on this list:** the code works, and what
stands between that and it *staying* working is a deck checked into `examples/` with a line in
the test docs. 2D periodic (C1, D5, M4) is **deferred** on Christian's instruction — it has not
been run and is not being pursued for now.

**Where it otherwise stands (2026-08-11, spot-checked against the tree):** 10 of 53 items
ticked. The blocker is unchanged and is the first item on the list — **T0, the reference
deck, is still unwritten**, so M1 has never been reachable and nothing downstream of it has
been exercised. The code fixes that did land (E5/`4a42d982`, D2/`43474c9f`, D1, A3) are real
and verified in tree.

**Two health warnings for anyone picking this up.** Every `file:line` citation is baselined
to 2026-07-10 and has drifted — believe the code, not the number. And the tick state has
lagged reality at least once: B6 was fixed on 2026-07-15 and sat unticked for four weeks, so
**verify an item against the tree before working it**, the way A3/A7/B6 were verified on
2026-08-11.

> **2026-08-09 currentness sweep.** Two of the sprint's hardest items landed on 2026-08-04,
> out of sprint order, driven by Gregory's multi-layer alternation report rather than by the
> T0-first plan:
> - **E5 (EF_QUAD4TS) — landed in `4a42d982`.** The `-mS[1]` top-dof flip is gone, the curl
>   parameters are the true thickness derivatives (`B(0,0) = -0.5`, `B(0,1) = +0.5`,
>   `cl_EF_QUAD4TS.cpp:143-144`), ∇η comes from real element geometry instead of an
>   axis-aligned assumption, and `mDetJ = area/4`. This also dissolves E5's premise: the
>   basis no longer *depends* on the `mS[1] = -mS[0]` invariant, so the missing
>   edge-direction writer stops being a correctness question (it stays a T1 hygiene
>   question). Root-cause record: `devlog/dl20260804_2d_thinshell_layer_alternation.md`.
> - **D2's code half — landed in `43474c9f`** (`get_top_nodes`, `src/mesh/meshtools.cpp`),
>   fixing the flipped top-surface edge-to-φ tie.
> - **E5(c)'s analytic patch test — landed as the permanent battery**, not as a one-off:
>   `tests/fem/test_InterfaceOrientation.cpp` (10 cases under `make check-fast`, both
>   backends) covers unit circulation, the sign sweep, inter-layer continuity and the N=4
>   Ampère telescope, and was mutation-checked red/green against the historical `-mS[1]`
>   flip. Physics expected values still pend Christian's sign-off
>   (`falsification_tooling.md`).
>
> **What this does NOT mean:** M1 has not been reached. T0 (the reference deck) is still
> unwritten, so nothing below has been exercised end-to-end; the greg2 4-layer smoke rerun
> is the outstanding evidence (`debt_register.md` DR-15).
>
> **Line-number caveat:** all `file:line` citations below are baselined to HEAD of
> 2026-07-10 and have drifted — spot-checked examples: D4's stale error message is now at
> `cl_MaxwellFactory.cpp:1356` (was :1295-1297). Re-locate by symbol, not by line.
**Second iteration (2026-07-10):** all five verify-first items were resolved by dedicated
verification passes (Claude subagents, Grok as a third voice, Codex audit; exchange thread
`tmp/ai_exchange/2d_thinshell_todo.md`, distilled into
`devlog/dl20260710_2d_thinshell_plan_v2.md`). Key changes: B1 already fixed on HEAD;
E1/E2 gate multi-layer only (M2, not M1); D2 and A4 upgraded from suspected to
confirmed; B3 downgraded to defensive validation; E5 escalated (missing
edge-direction initialization); new day-1 task T0 (reference deck).
**Labels:** A = Stage 1 (spatial cuts), B = Stage 2 (cohomology), C = Stage 3 (thin cuts),
D = Stage 4 (shell insertion), E = Stage 5 (solver), T = sprint tasks, M = milestones.

## Sprint scope (v2)

**In (two-week, first-order):** T0, A1, A4, D1 (conditional), D2, D3, E1, E2, E4, E5 —
plus C3 on first failure during bring-up, and the M1–M3 milestones (M4 stretch).
**Out (explicitly deferred):** the order-2 block (A3, C5, D6, D7, E8), all CLEANUP
items (A6, A7, B5, B6, B7, D8, D9, D10), B2 (unless the sprint model uses symmetry),
C2 (unless cut asserts fire), E7.

### Suggested two-week sequence

- [ ] **Days 1–2 (S):** T0 reference deck + A1 validation; B3 deck check; instrumented
      debug run until the first failure; D3 one-liner while working in `cl_MaxwellFactory.cpp`.
- [ ] **Days 3–5 (M):** D2 pairing fix + A4 orientation sign (one numeric sanity run
      each); re-run hanging-edge verification.
- [ ] **Days 5–8 (M/L):** T1 + E5 edge-direction mechanism + EF_QUAD4TS patch test;
      E4 sizing fix with the E6 `h_ts_*` sweep in the same pass.
- [ ] **Days 8–9:** M1 green (serial, linear, single-layer, air both sides); D1 if the
      deck has a conductor above the tape; C3 if cut asserts fired during bring-up.
- [ ] **Days 10–12:** E1 + E2 together → M2 green (multi-layer).
- [ ] **Day 13:** M3 parallel smoke test (B1 already fixed).
- [ ] **Day 14:** buffer; if green early, start M4 (C1, D5).

## Stage 0 — Sprint infrastructure (new in v2)

- [ ] **T0 — HIGH (day 1) — author a known-good 2D reference deck.** No 2D thin-shell
      input deck exists in the tree (`tmp/examples/Validation/Parallel_Planes` is 3D;
      `CCT_2D` has no `thinshell` section). Author a minimal serial deck — one straight
      tape in air, single layer, terminal curves, sigmoid current (reuse the BC block from
      `tmp/examples/Conductor_Symmetry_2D/input.conf`) — plus a two-layer variant for M2.
      This is the top schedule risk (Grok concurs): without it, every fix below is
      unverifiable and integration discovery consumes the schedule.
- [ ] **T1 — MEDIUM — cross-check how 3D PENTA6TS survives the missing edge-direction
      initialization found under E5.** Either the audit missed a writer, or
      3D runs with default directions and its EF tolerates them — determine which before
      copying any 3D convention into EF_QUAD4TS. (medium confidence that no writer
      exists; 3D works in production, so one of the two explanations must hold.)

## Stage 1 — Spatial cut creation

- [ ] **A1 — HIGH — validate 2D terminal-curve association.** Add a `BELFEM_ERROR`
      (rank 0) in `MaxwellFactory::read_thin_shell_data` when a 2D Protoshell ends up with
      no terminal curves (cl_MaxwellFactory.cpp:2500-2530), so a missing or misordered
      `topology/curves` entry fails loudly instead of silently tearing open the tape endpoint
      (= gap G1 in the analysis). Also document in the input-format docs that, in 2D, the curve's
      first sideset ID must be the tape sideset.
- [ ] **A2 — MEDIUM — relax the 2D terminal-curve match rule.** Decide whether the 2D
      match should also accept the tape id anywhere in the curve's sideset list
      (mirroring the 3D a/b swap, cl_MaxwellFactory.cpp:2514-2530), rather than `sideset_a`
      only (:2507); implement the chosen rule.
- [x] **A3 — HIGH (order-2, one-liner) — LINE3 missing `break` in `close_terminal_loops`.**
      ~~Fix the fall-through after the node-container `ElementType::LINE3` case
      (cl_CutFactory.cpp:2015-2030) — it falls into `default:` and aborts with "Invalid
      element type at terminal loop" (= gap G2). v2 note: the arclength LINE3 switch
      (:1962-1987) already has its `break`; only the node-container switch is affected.
      It is off the sprint's critical path (LINE3 = order-2); fix opportunistically — this also
      unblocks 3D order-2.~~
      **FIXED by Christian 2026-08-11.** Verified in tree: the node-container `LINE3` case
      now carries its own `break` before `default:` (`cl_CutFactory.cpp:2030-2040`), so an
      order-2 terminal loop no longer falls through into
      "Invalid element type at terminal loop".


- [x] **A4 — HIGH — 2D terminal-curve orientation sign.** ~~verify-first~~ **CONFIRMED
      2026-07-10 (high):** `orient_terminal_curves_2D` computes `cross(shell, tape)`
      with hard-coded shell normal +z, while the 3D reference computes
      `cross(boundary, shell)` after the a/b swap — anticommuted operand order, so the
      2D curve direction is mirrored relative to the validated 3D convention
      (cl_CutFactory.cpp, 2D ~:2140-2220 vs 3D ~:2039-2096 at HEAD;
      doc/thin_shell_facet_orientation.md documents no explicit winding convention).
      Remaining work (in sprint, days 3–5): one numeric sanity run on the T0 deck
      checking imposed current direction, then flip the operand order. *(Box ticked for the
      verification half; the code fix is tracked in the sprint sequence.)*
- [ ] **A5 — LOW — assert exposure in `orient_terminal_curves_2D`.** Revisit
      `BELFEM_ASSERT(tFacet != nullptr, "No surface found on tape/boundary")`
      once A1/A2 settle the terminal-curve semantics (the assert
      fires when the first segment's nodes don't lie on ≥2 corner nodes of a tape facet).
- [ ] **A6 — CLEANUP (deferred) — dead spatial-cut machinery.** Remove or clearly quarantine
      the dead code in CutFactory — `create_sidesets_2d/3d`, `create_cut_sideset_2d/3d`,
      `flag_nodes_and_facets_of_tape_sidesets` (contains a divergent second 2D endpoint
      heuristic), `create_facet`, `connect_facet_to_slave`, `check_element_types`,
      `save_edges` — plus `SideSetFactory` and `CutProcessorManual`
      (no callers repo-wide) (= gap G5).
- [x] **A7 — CLEANUP (deferred) — botched rename.** ~~Rename `tProtoQUAD4TS` back to a
      sensible name in the 3D branch of `duplicate_nodes_on_face_sidesets`.~~
      **DONE 2026-08-11:** renamed to `tOtherShell`, which is what the loop is — the inner
      pass over `mProtoshells` that skips the current shell by id, exactly as the comment
      above it says ("look if other thin shells have side nodes flagged"). Three occurrences,
      all in `cl_CutFactory.cpp`; zero left in the tree; TU syntax-checks clean.

## Stage 2 — Cohomology computation

- [x] **B1 — HIGH — 2D branch for `compute_element_adjacencies`.** **FIXED on HEAD
      (2026-07-02 CW-mesh session, verified 2026-07-10):** the 2D branch now resets
      per-element containers and calls `Mesh::reset_connectivity(ElementToElement)`
      (cl_CutFactory.cpp:2761-2773; new API cl_Mesh.hpp:841) instead of stamping an
      empty adjacency — the next kernel recomputes the graph from scratch, so 2D
      parallel partitioning is no longer poisoned. M3 becomes a verification run.
- [ ] **B2 — MEDIUM (deferred unless sprint model uses symmetry) — 2D symmetry
      consistency sweep.** Add the 2D analog of the face-consistency sweep in
      `CutFactory::unflag_symmetry_sidesets`: after unflagging symmetry edges/nodes,
      unflag top ELEMENTS whose edges were all unflagged (the existing face loop is a
      silent no-op in 2D), so 2D symmetry models don't feed an inconsistent complex to
      Smith/coreduction.
- [x] **B3 — ~~MEDIUM~~ LOW — 2D terminal convention in `suggest_Homology`.**
      ~~verify-first~~ **LARGELY REFUTED 2026-07-10 (high):** the parser doubles 2D
      terminal lists (`tDomainGroupsOut = tDomainGroupsIn`,
      cl_MaxwellBoundaryConditionFactory.cpp:125, concatenated :144-151), so
      `size()/2 ≥ 1` for every parsed deck and the halving at cl_Homology.cpp:259-267
      picks exactly the input half. The empty-generator death only occurs if the
      doubling is bypassed. Residual (out of sprint): add a defensive validation error
      for undoubled lists. Low-cost in-sprint action: one deck check on T0 (day 1).
- [ ] **B4 — LOW — element-type guard in `CutProcessor`.** The ctor infers TRI3/TRI6 from
      dim+order (cl_CutProcessor.cpp:33-35); add an explicit check that 2D φ-blocks are
      triangles so a QUAD block errors out instead of silently running TRI edge kernels on
      4-edge elements.
- [ ] **B5 — CLEANUP/DOC (deferred, but document during M1 bring-up) — `mCuts` no-op.**
      Either populate `CutFactory::mCuts` or delete `cuts()` and the no-op
      `DomainType::Cut` loop (cl_MaxwellFactory.cpp:837). v2 note (Codex): live cut
      sidesets are added at cl_CutData.cpp:180, and
      `set_block_types_in_magnetic_equation` excludes `Cut`/`Default` sidesets from
      assembly (cl_MaxwellFactory.cpp:1737) — the no-op may be deliberate (cuts as
      topology-only sidesets). Document the intent, since it also determines whether
      C3's master-only facets can ever reach an assembly consumer (relevant to C4).
- [x] **B6 — CLEANUP (deferred) — stray pasted declaration.** ~~Delete
      `void collect_thin_cut_edges(...)` inside `determine_cut_case_3d`
      (cl_CutData.cpp:460, still present at HEAD).~~
      **ALREADY DONE, ticked 2026-08-11:** zero occurrences of `collect_thin_cut_edges`
      remain anywhere in `src/`. `devlog/dl20260715_spfa_clean_implementation.md` records the
      removal in `6ee74493` (2026-07-15); the box was simply never ticked.


- [ ] **B8 — HIGH — facet-flag leak from `fix_facet_masters` tears all 2D tape tips
      (FIX APPLIED 2026-07-10, pending build+run verification).** Full trail in the
      exchange thread (2026-07-10 entries) and
      `devlog/dl20260710_2d_thinshell_t0_flag_leak.md`. Symptom chain, all measured
      on the T0 11-tape stack: `fix_facet_masters` leaves every equal-domain-type
      facet FLAGGED when its propagation queue is empty (all-air 2D model has no
      domain-type-contrast seed facets) → the tip test in
      `duplicate_nodes_on_face_sidesets` counts connector-line facets as flagged
      (tCount 2–3, never 1) → zero tips protected → all tape tips duplicated (tapes
      silently torn) → `close_terminal_loops` closes through disconnected node sets
      (4 nonzero-boundary nodes per chain) → open suggested chains pair with the
      cohomology basis only by accidental edge overlap (Grok+Codex confirmed the
      cycle/representative math; Pellikka pairs terminals as ∂Σ cycles) →
      rank-deficient tTemp (difference/anchor structure) → Matrix OOB in
      `updatekGeneratorsFromHomology` (cl_Cohomology.cpp:560). `suggest_Homology`
      and `close_terminal_loops` are CORRECT — the closure failed upstream.
      Fixes applied: (1) `unflag_all_facets()` after the propagation loop in
      `fix_facet_masters`; (2) `unflag_all_facets()` before tape-flagging in the 2D
      branch of `duplicate_nodes_on_face_sidesets`; (3) `BELFEM_ERROR` when an open
      terminal curve finds ≠ 2 tips. Verify: bracketed [1]..[11] deck passes
      cohomology (rank 11), then re-test the unbracketed union form (closed chains
      make the union a valid cycle — may work now, like 3D corc). Residual
      (deferred): cosmetic rank `BELFEM_ERROR` in `updatekGeneratorsFromHomology`;
      note `Chain::getBoundary()` is nullptr for isBound chains and
      `addSimplexToChain` never maintains boundaries — future cycle invariants need
      a manual boundary walk.
- [ ] **B9 — MEDIUM (found in passing during B8) — `fix_facet_masters` orientation
      propagation never runs on all-air 2D thin-shell meshes.** The queue is seeded
      only by domain-type-contrast facets; with every block Air, no facet is ever
      queued, so the normal-alignment sweep (cl_MaxwellFactory.cpp:~1105-1128) is a
      no-op and tape facet orientation stays raw gmsh order. Interacts with A4
      (terminal-curve sign) and master-side selection; decide whether 2D thin-shell
      tapes need an explicit orientation pass before M1 physics is trusted.
- [ ] **B7 — CLEANUP/DOC (deferred) — empirical sign hack.** Resolve (or document) the
      sign hack in `Homology::reorient_generators` — "should be +1.0 … another sign
      error elsewhere in 3-D"; 2D currently avoids the hack. Revisit together with the
      A4 sign fix, which may be the "elsewhere".

## Stage 3 — Thin-cut computation

- [ ] **DEFERRED 2026-08-11 (Christian: 2D periodic is not being run for now) — C1 — HIGH (blocks 2D periodic; M4) — periodic edges in 2D `collect_facets`.**
      Extend the 2D branch (cl_CutProcessor.cpp:307-414) to use
      `mPhiBoundariesAndPeriodic` (as the 3D branch does) instead of `mPhiBoundaries`
      only (:376-383).
- [ ] **C2 — MEDIUM (deferred unless cut asserts fire) — 2D dangling-cut pruning.** Port
      the 3D iterative pruning loop to the 2D branch of `collect_facets` — 1D analog:
      iteratively drop thin-cut edges whose endpoint is not shared with another kept cut
      edge or a protected boundary, so dead-end cut chains don't reach
      `CutData::add_thin_cut_sidesets_to_mesh` and trigger the cl_CutData.cpp:89-90 asserts.
- [ ] **C3 — HIGH (M1-contingent: fix on first failure during bring-up) — symmetrize 2D
      `add_thin_cut_sidesets_to_mesh`.** Align the 2D branch (cl_CutData.cpp:55-141)
      with 3D: set BOTH master and slave on 2D cut facets (currently master-only at
      :139), replace the coordinate-epsilon side identification (:98-134, comment
      "something is off with the original pointers") with topological side selection,
      and add a boundary fallback equivalent to the 3D Step-5a slave-only path so cuts
      touching non-φ boundaries do not assert in debug or corrupt state in release.
      v2 note: interior cuts on a clean T0 deck may pass with master-only facets if no
      consumer dereferences `slave()` — therefore this is contingent, not day-1; it stays HIGH
      because the failure mode in release is silent.
- [ ] **C4 — MEDIUM (after C3) — audit Cut-sideset consumers.** Audit 2D consumers of
      Cut-sideset facets for `slave()`-dereference assumptions (IWG cut path, DofManager)
      and remove any workarounds that assumed master-only facets.
- [x] ~~**C5 — (order-2, deferred) — second-order T-matrix weights.** Fix the weights
      flagged at cl_CutProcessor.cpp:97 (`// todo: fix weights for second order!`) —
      needed for 2D TRI6 (and 3D TET10).~~ **RULED STALE 2026-08-27 (Christian, at
      DR-103's strike): there is nothing to fix.** The cut constraint is a nodal
      T-matrix tie (`CutSet::create_duplicates`, all weights 1.0); the jump
      [phi] = I is constant over the cut and phi is Lagrange at every order, so
      the unit weights are exact for TRI6/TET10 too — no facet quadrature is ever
      formed. The source todo was removed and a rationale comment added at the
      weights in `cl_CutSet.cpp`. The genuine order-2 residue in this region
      (diagonal-case midside coverage, missing hanging-source cascade in
      `create_duplicates`, no TRI6/TET10 homology fixture) is tracked as
      `debt_register.md` DR-122.

## Stage 4 — Thin-shell insertion

- [x] **D1 — HIGH (conditional on geometry) — LINE-facet case in `to_master_orientation`. FIXED 2026-08-09** (`fn_to_master_orientation.cpp:596-604`, tracked as `debt_register.md` DR-14). The `case 1` is the identity: a LINE2 facet carries exactly one edge, so there is nothing to permute. The "relative sign from the master/slave node order" half of the original prescription turned out **not** to belong here — the caller already derives the ±1 weight from the node-index comparison (`cl_MaxwellFactory.cpp:1758-1770`) after the Node overload has reversed the LINE2 nodes (`:32-36`). Still unexercised: no deck yet puts a conductor above a 2-D shell. Original text:
      Add a facet-edge-count == 1 case to
      `to_master_orientation(Facet*, Cell<Edge*>&, ...)`
      (fn_to_master_orientation.cpp:594-681, default error :676-681): identity index plus
      relative sign from the master/slave node order. Without it,
      `hang_thinshell_edges_on_edges_top` hard aborts whenever a 2D shell has a
      conductor slave above it. v2 note (Grok): the air-above path uses the node
      overload, which already handles LINE2 (:32-36) — D1 gates conductor-above
      geometries only, not every M1. Schedule when the deck needs it.
- [x] **D2 — HIGH — 2D top-node pairing. CODE FIX LANDED 2026-08-04 (`43474c9f`, `get_top_nodes` in `src/mesh/meshtools.cpp`)** — the box previously covered only the verification half. ~~verify-first~~ **CONFIRMED CROSS-WIRED
      2026-07-10 (high, two independent traces):** QUAD4TS layer insertion order is
      n2 = top(F1), n3 = top(F0) (cl_ThinShellFactory.cpp:1218-1221);
      `get_top_nodes(QUAD)` returns facet 2 = {n2, n3} (meshtools.cpp:1313-1316,
      cl_Element_QUAD4TS.hpp:74-78); combined with the `to_master_orientation(LINE2)`
      flip and the top-edge storage order {n3, n2}, the lookup at
      cl_MaxwellFactory.cpp:~1535 pairs top(F0) with the volume node at F1 and vice
      versa — position-misaligned, unlike the PENTA path
      (cl_Element_PENTA6TS.hpp:100-106). Remaining work (in sprint, days 3–5): pick the
      fix location — reverse `get_top_nodes(QUAD)` to return {n3, n2}, or skip the
      orientation flip on the 2D top path — and confirm numerical impact with the
      two-element run. *(Box ticked for the verification half; the code fix is tracked
      in the sprint sequence.)*
- [ ] **D3 — MEDIUM (both dims, one-liner) — `break` vs `continue`.** Change `break` to
      `continue` at cl_MaxwellFactory.cpp:1292 so one empty ThinShell does not silently
      skip edge-hanging for all remaining shells.
- [ ] **D4 — LOW — stale error message.** Fix cl_MaxwellFactory.cpp:1295-1297 ("Only
      linear triangular elements are supported" — LINE2 is in the accepted list).
- [ ] **DEFERRED 2026-08-11 (Christian: 2D periodic is not being run for now) — D5 — HIGH (blocks 2D periodic; M4) — QUAD4TS branch in `create_periodic_sideset`.**
      The PENTA-only assert and the `f<3` side-facet loop
      (cl_ThinShellFactory.cpp:2305 ff., called :310-311) make 2D + periodic assert in
      debug and mis-enumerate facets in release.
- [ ] **D6 — (order-2, deferred) — transposed indexing in `process_nodes_line3`.** Fix
      cl_ThinShellFactory.cpp:654 ff.: normals are accumulated row-per-node into a
      `(3, numNodes)` matrix that the read side consumes column-per-node into a
      length-2 vector — out of bounds for node index ≥ 3. Adopt the column-per-node
      layout of `process_nodes_line2`/`_tri3`.
- [ ] **D7 — (order-2, deferred) — `QUAD9TS::get_nodes_of_edge`.** Give QUAD9TS a proper
      specialization (cl_Element_QUAD9TS.hpp currently inherits `get_nodes_of_facet`
      semantics via cl_ElementTemplate.hpp) so edge-based code cannot silently
      misread its three curve-edges.
- [ ] **D8 — CLEANUP (deferred) — dead legacy layers-Element path.** Delete the
      `Element(SideSet*, DofManager*, Facet*, Cell<Facet*>& aLayers, ...)` ctor
      (cl_FEM_Element.cpp:90-111), `link_dofs_thin_shell` (:902 ff.), and
      `compute_edge_directions_thinshell` (:1415-1431, the "tape roller" physical_tag
      mechanism that no longer exists). Verified no callers; its `dims==2` assert can
      mislead audits into treating it as a live 2D path. v2 note: do NOT delete before
      E5/T1 resolve where layer-element edge directions should come from — this dead
      code is the only in-tree precedent for setting them.
- [ ] **D9 — CLEANUP (3D, deferred) — `process_nodes_tri6` tM bug.** Fix the
      triple write of `tM(0)` instead of `tM(0)/tM(1)/tM(2)`
      (cl_ThinShellFactory.cpp:914-916, still present at HEAD).
- [ ] **D10 — CLEANUP (deferred) — complete side-connector removal.** Remove the
      remaining side-connector residue: `EF_HEX8TS`
      (cl_EdgeFunctionFactory.cpp, cl_EF_HEX8TS.cpp), HEX8TS cases in
      meshtools.cpp and cl_FEM_Calculator.cpp, fn_num_nedelec_dofs.hpp,
      cl_Element_Factory.cpp, `ThinShellFactory::compute_binomial_vectors`
      and `compute_side_indices`, `ThinShellFactory::check_inside` (dead; also
      makes the hardcoded `mElementMapper->set_dimension(3)` moot).

## Stage 5 — Main solver + preprocessors

- [ ] **E1 — CRITICAL for M2 (release-mode memory corruption; NOT on the M1 path) —
      2D variant of `maxwell::h_ghost`.** Replace the hard-coded `n=6; m=3; d=3`
      (mt_maxwell_h.cpp:1878-1880) and the TRI3/PENTA6TS asserts (:1771, :1873-1875)
      with sizes taken from the group's edge function (`ndofs()`; QUAD4TS has 2 edge
      dofs per side, fn_num_nedelec_dofs.hpp:28). v2 re-scope (E3 resolved): a
      single-layer 2D shell creates NO Ghost sideset (early return `nb==1`,
      cl_ThinShellFactory.cpp:1837; sideset gate `g>0` :289), so this is the
      out-of-bounds path for every MULTI-layer 2D shell only (mt_maxwell_h.cpp:1767 ff.).
- [ ] **E2 — CRITICAL for M2 (pairs with E1) — type-dependent Ghost scratch sizing.** In
      `IWG_Maxwell::create_custom_vectors_and_matrices` (cl_IWG_Maxwell.cpp:607-614):
      `e = 6` and `D± (3, e)` are PENTA6TS-specific; use `e = 2` and 2D-row `D±` for
      QUAD4TS groups.
- [x] **E3 — INFO — single-layer Ghost check.** **ANSWERED 2026-07-10 (high):** a
      single-layer shell (1 thickness → 2 layers → 1 block for order 1) hits the
      `nb == 1` early return in `create_ghost_facets` (cl_ThinShellFactory.cpp:1837)
      and the `g > 0` sideset gate (:289) — no Ghost sideset exists, `h_ghost` is never
      bound. E1/E2 block ONLY multi-layer 2D shells; M1 is unblocked without them.
- [ ] **E4 — HIGH (M1 path) — dimension-consistent thin-shell work vectors.**
      ~~suspected~~ **CONFIRMED 2026-07-10 (high): silent wrong-physics in BOTH debug
      and release.** `"bt"`, `"bn"`, `"n"` are hard-coded size 3 while `"b"`, `"bm"`,
      `"bs"`, `"j"`, `"h"` are sized `mesh::dimension(type)` (= 2 for QUAD4TS)
      (cl_IWG_Maxwell.cpp:585-598). Armadillo expression assignment silently RESIZES
      the target (cl_AR_Vector.hpp:375-381 delegates to `arma::Mat::operator=`;
      BELFEM's asserting `op_VectorPlus` is bypassed by expression templates), so
      `compute_bn` round-trips `bn` 3→2→3 (mt_maxwell_h.hpp:146-164) and `b = bt + bn`
      (mt_maxwell_h.cpp:137/:2271) silently grows `b` to size 3 — no assert fires
      anywhere. Pick one convention (size-by-dimension, or pad all to 3) and align
      `compute_bn`/`h_ts_*` with it.
- [x] **E5 — HIGH (M1 path) — EF_QUAD4TS edge directions + validation. FIXED 2026-08-04 (`4a42d982`) + battery.** See the status block; (b) and (c) are done, and (a) is dissolved rather than implemented — the corrected basis does not rely on the edge-direction invariant. Original text:
      ~~verify-first~~ **ESCALATED 2026-07-10 (medium-high):** the verification pass
      found NO live writer of layer-element edge directions — the legacy tape-roller
      path (cl_FEM_Element.cpp:90-111 → :1415-1431) is dead, and no replacement was
      located, so `edge_directions()` returns the default bitset → mS = {-1, -1},
      violating the `mS[1] = -mS[0]` invariant assumed by `EF_QUAD4TS::E()`
      (cl_EF_QUAD4TS.cpp:185-187; `link()` :71 copies whatever the element reports).
      Caveat (hence medium-high, see T1): 3D PENTA6TS works in production, so either a
      writer was missed or 3D tolerates the defaults — resolve T1 first. Then:
      (a) implement/confirm the edge-direction mechanism for QUAD4TS layer elements;
      (b) derive/verify the curl coefficients `B(0,0)=B(0,1)=+0.5` (:157-167) against
      the PENTA6TS reduction (Messe et al. 2023, paper1); (c) analytic patch test.
- [ ] **E6 — MEDIUM (pairs with E4, same sprint block) — audit 2D `h_ts_*` bindings.**
      Confirm every `h_ts_*` variant reachable from a QUAD4TS ThinShell block in
      `IWG_Maxwell::link_to_group` (cl_IWG_Maxwell.cpp:437-527) behaves with 2D-sized
      Calculator data (bodies were not audited). v2 note (Codex): do this while fixing
      E4 — all reachable `h_ts_*` bodies share the same vector-size risk, so one pass
      over them settles both items.
- [ ] **E7 — LOW (deferred) — BC dimensionality cross-check.**
      `MaxwellBoundaryConditionFactory` stores `mNumDimensions` but never branches on
      it; 2D dimensionality is inferred from missing "output terminal" keys, so a
      3D-style deck on a 2D mesh fails late or scales silently (voltage `length`). Add
      an explicit check.
- [ ] **E8 — (order-2, deferred, last) — 2D second-order enablement block.** Implement
      `EF_QUAD9TS`, a 2D analog of the TRI6/TET10-only `TMatrix` + PART 3 facet
      weights, and extend the facet-source branch in `DofData` to 2D; then lift the
      order gates at cl_MaxwellFactory.cpp:1176 and cl_FEM_Calculator.cpp:637. Do this
      only after A3, C5, D6, and D7.

## Verification milestones (v2 exercise lists corrected)

- [ ] **M1 — Milestone A (serial, linear, non-periodic, single-layer):** the T0 tape
      deck runs end-to-end in debug — exercises A1, A4, D2, E4, E5 (+ D1 only if the
      deck has a conductor above the tape; + C3 if cut asserts fire).
      ~~exercises A1/A3, D1/D2, E1/E2/E4~~ (v1 list was wrong: E1/E2 are multi-layer
      only, A3 is order-2).
- [ ] **M2 — Milestone B (multi-layer):** two-layer T0 variant (Ghost sidesets active) —
      exercises E1/E2 under load, then E6.
- [ ] **M3 — Milestone C (parallel):** 2D parallel run — B1 is already fixed, so this is
      a verification/smoke run, not a fix.
- [ ] **DEFERRED 2026-08-11 (Christian) — M4 — Milestone D (periodic, stretch):** 2D periodic model — exercises C1 and D5.
- [ ] **M5 — Milestone E (order-2, out of sprint):** 2D second order — the (order-2)
      block (A3, C5, D6, D7, E8).
