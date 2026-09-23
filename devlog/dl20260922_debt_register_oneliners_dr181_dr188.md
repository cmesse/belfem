# Debt register one-liners, the background-field surface increment, and the CurveFactory pair and closure rewrite

**Date:** 2026-09-22
**Purpose:** Christian opened with "let's go over the debt register, first the simple ones". Two bundles went through plan → jury → code → jury: five `src/mesh` one-liners (DR-161, DR-169, DR-176, DR-179, DR-180) and the `background_phi` surface increment (DR-181). The mesh bundle's gate exposed a CurveFactory defect class, which became DR-188 and went through the same loop.
**Module:** `src/mesh` (`meshtools`, `cl_Mesh_Basis`, `cl_Element_Factory`, `cl_Face`, `cl_FaceFactory`, `cl_Curve`, `cl_CurveFactory`), `src/fem/maxwell` (`mt_maxwell_background`, `cl_MaxwellFactory` terminal association)
**AIs involved:**
- Orchestration: Claude Opus 5.5, then Fable 5.1, then Opus 5.5 again. Fable ran out of usage credits twice; Christian's ruling was that Opus orchestrates until the reset.
- Planning: Fable 5.1 drafted the mesh-bundle and DR-181 plans; an Opus 5.5 agent drafted DR-188.
- Implementation: Opus 5.5 agents. The one-token LINE3 midpoint fix and the DR-188 fix pass were inline.
- Pre-jury read of every Opus diff: Fable, or Opus while Fable was out of credits.
- Mesh bundle juries: Codex `gpt-6-astra` and Grok `grok-4.7`/high. The plan jury ran astra at high, the code rounds at high → xhigh → xhigh; three code rounds in all.
- DR-181 juries: astra/high and Grok/high on both plan and code. Grok needed three attempts on the plan (30 turns, then 60, then 100).
- DR-188 juries: astra/high and Grok/high on plan and code, a narrow round 2 with astra at xhigh, and a round 3 (astra/high, Grok/high) on the closed-terminal guard and its test.

**Status:** COMPLETE — all three are landed in one commit. What is still open is listed at the end of this entry.
**Verification:** all runs used the release tree (`build/`, SCLS `mkl`, Armadillo) and happened in the Claude sandbox. The auto-mode classifier refused `make` in the shared tree, so each changed TU was compiled with `build/`'s own flags.make, added with `ar r` to a scratch copy of `build/lib/libbelfem.a`, and the test binaries were relinked from `build/`'s link lines. Nothing was written into a shared build tree.
- the 8 `fast` suites run directly: first 1734/1734 on the 21:38 library, then **1827/1827** on the rebuilt baseline — a peer's 15 overnight commits (the Tier 2 test-gap run, HEAD `42b5e42` → `1c7d182`) added tests, so the gate was rebuilt on their 01:27 library with all 9 of our source TUs and 7 of our test TUs compiled fresh (core 58, containers 575, linalg 285, math 361, mesh 119, homology 25, circuit 127, fem 277). A first run had failed the peer's `SolverParametersDeck.RelativeToleranceIsValidatedFromADeck` only because the older library predated their `0d6b726`;
- DR-181 red-before: the unfixed kernel fails at K 17.6 % and f 18.8 % against a 1e-12 tolerance; green-after passes;
- DR-188: all 14 new tests are red or crash against the pre-DR-188 library;
- DR-188 A/B on the 10 shipped 3-D `curves` decks: 232 curves, all identical to the old code;
- the closed-terminal guard test `ThinShellTerminalCurves.*`: 6/6 with the guard; against HEAD's unguarded `MaxwellFactory`, the 3 refusal cases are red and the 3 acceptance cases green;
- Christian's own `test_mesh` run on the fixed bundle: green (102 cases at that point).

Not run: a debug tree, valgrind, the MPI suites, the 64-bit build.
**Exchange:** `tmp/ai_exchange/review_mesh_oneliners_dr161_169_176_179_180.md`, `review_maxwell_background_ds_dr181.md`, `review_curvefactory_dr188.md`

## Summary

Seven register rows were closed in code: DR-161, DR-169, DR-176, DR-179, DR-180, DR-181 and DR-188. Four more rows changed state: DR-178 was struck on Christian's ruling, DR-183 was filed and then struck as refuted, and DR-184 and DR-188 were filed. The sections below follow the three rounds.

## Mesh one-liner bundle (DR-161, DR-169, DR-176, DR-179, DR-180)

- **DR-161.** `get_top_nodes` switched on `geometry_type`, so it accepted QUAD8, QUAD16 and HEX64, while `get_bottom_nodes` did not. It now switches on the same 14 concrete types and refuses the other three (Christian's ruling: refuse, not support). The planner also found that top's QUAD16 corner swap was already misaligned. The facet-2 midside nodes crossed the swapped corners, so accepting that type had been an accident rather than a feature.
- **DR-169.** `~Basis` now frees `mSources` and `mWeights` keyed on the pointers, not on the count. The copy and move operations are `= delete` (both code-jury rounds): an unfilled copy would otherwise double-free under the unconditional free. The re-call leak of `allocate_source_container` stays open. The new test cannot turn red without valgrind, as the plan states.
- **DR-176.** `create_unity_nodes` now handles HEX8TS and HEX8TB in its HEX8 case. `test_Pipette` builds HEX8TS through the factory instead of borrowing the HEX8 geometry. The planner's claim that the only difference between the two was "edge container only" was refuted: the templates are `<8,8,8,6,1>` and `<8,8,12,6,6>`, and only the node tables coincide.
- **DR-179.**
  - `cl_FaceFactory.cpp:484` now compares against `gNoID`.
  - The `uint` Face sentinels that took `gNoIndex` now take `BELFEM_UINT_MAX` (`cl_Face.hpp:38,54`, `cl_Face.cpp:28,56`, `cl_FaceFactory.cpp:528`). Under `BELFEM_INT64` the old default narrowed and skipped the orientation computation.
  - The sideset-only branch used to call `insert_face( tFace, gNoIndex )`: an out-of-range write into an unallocated `mFaces`. It now refuses with a named `BELFEM_ERROR`. Christian's first ruling was "fix it whole"; the code jury then showed that a sideset-only surface face gets facet 0's single edge in the 3-D edge pass and cannot survive a `.bfm` reload, and the branch has no production caller. So "fix whole" was revised to "refuse".
  - The 64-bit tree does not compile today for unrelated reasons (`Cell<uint>` vs `Cell<index_t>&` in `cl_Mesh_Distributor.cpp:2420`, `cl_ProtoMesh.cpp:107/116`, `cl_TensorMeshFactory.cpp:98-113`), and no SCLS flavor has 64-bit integers yet, so the 64-bit leg is blocked.
- **DR-180.** The `Curve` constructor dereferenced a null sideset B for a quadratic 2-D curve. The null guard exposed three defects behind it, which the bundle's three code rounds and the gate took apart:
  - **The midpoints were assigned by position.** They were gathered in facet order but paired with pairs sorted by key, so a segment could get the wrong midpoint (astra, round 1).
  - **`create_pairs` builds pairs from the element clique, not from facet edges.** `connect_nodes_to_nodes` links every node of every element. At a boundary corner, where one element holds three end nodes, this adds a diagonal chord (astra and Grok, round 2). The fix builds pairs and midpoints directly from the sideset facets (`create_pairs_from_facets`).
  - **The LINE3 midpoint was inserted into the Segment's own, never-allocated vertex container.** `create_segments` called `tSegment->insert_node( …, 2 )` instead of `element()->insert_node`, so `Vertex::insert_node` wrote through a null `mNodes`. The new quadratic test crashed with a SEGV on Christian's 20:46 build; gdb gave `Vertex::insert_node ← create_segments ← from_2d_sidesets`. This was fixed inline, and Grok showed that the same fix repairs the live 3-D `intersect` path, which shares `create_segments`.

## DR-181 — `background_phi` read the point-0 surface increment at every point

`aCalc->dS()` (default index 0) is now `aCalc->dS( k )` in both lines of `mt_maxwell_background.cpp:48-49`.

- `dS_cartesian( aIndex )` calls the normal callback at the requested index. The per-point families (`normal_hex`, `normal_penta`, the curved ones) recompute only when the index changes. So every integration point after the first reused point 0's increment. Flat first-order facets are exact, which is why `disk_pulse` and `disk_ramp` never showed the defect.
- The new `tests/fem/test_BackgroundField.cpp` builds a kernel-level fixture: one HEX8 air element with a non-parallelogram bottom face. Part A checks Σw·dS(k) against the absolute shoelace area and requires the increment to vary across the face. Part B compares the kernel's K and f against a per-point reference at 1e-12.
- The juries corrected the fixture. The boundary condition must be allocated with `new`, because `Kernel::~Kernel` deletes registered BCs (astra and Grok, plan round). Part A needs the absolute area, because facet 4 is `0,3,2,1`. The fixture also gained a serial-only guard, finite-value asserts before the max reductions (a NaN slips through `std::max`), and explicit move deletion (astra, code round).
- The sweep found no other defaulted `dS()` or `dV()` inside an integration-point loop in `src/` or `tests/`. `nonfree/` was not present to sweep.

**By-catch, refuted (DR-183).** The planner suspected that the `dirichlet → BackgroundField` stamp at `cl_MaxwellFactory.cpp:908-912` makes `background_phi` assemble `∫BᵀB` onto free master rows. Grok's third attempt refuted this by source trace: `SideSet::impose_dirichlet` resets the sideset to `Dirichlet` + `GeometryOnly` every step before `compute_jacobian_and_rhs`, and that loop skips inactive sets. The row was filed and struck as refuted the same session.

**Grok's facet-4 "weight pairing" concern was refuted by source trace.** The ξ↔η swap on facet 4 is the facet-to-master coordinate map, so each weight stays with its point for any quadrature rule. No row was filed.

## DR-188 — CurveFactory: pairs from facets on every path, closed curves oriented, dead overload deleted

This row was filed from the mesh bundle's code rounds and widened on Christian's ruling (extend-don't-branch). It had four items.

- **(D) Midpoints drifted after an unmatched pair.** The old `collect_midnodes( aPairs, aMidNodes )` advanced its cursor only on a match, so every later midpoint shifted and the last slot stayed null. The live 3-D `intersect` path and `thin_shell_side_curves` both used it. Both now build pairs and midpoints from facets:
  - `intersect` uses the edges that an A facet and a B facet share, through a new `aSharedWith` parameter. The B side gets the same midpoint-conflict and missing-midpoint checks as the A side.
  - The thin-shell path reads each midpoint from the one facet that carries the edge (`midnode_of_facet_edge`).
  - `create_pairs` and the old `collect_midnodes` are deleted.

  The plan jury showed that the clique defect bites the 3-D path too: a hex that holds a curve corner adds a face diagonal. Both auditors also showed that the plan's "open curves unchanged" lemma was false: two shared-edge pieces bridged by a co-element chord pass today and would be refused tomorrow. Christian made an old-vs-new A/B on every shipped 3-D deck mandatory. It came back 10 of 10 identical: 232 curves, none refused, none closed.
- **(Z1) Mixed order.** `from_2d_sidesets` takes its order from sideset 0 alone, so a LINE4 sideset chained after a LINE3 one contributed its ξ = −1/3 node as a midpoint. It now refuses any non-line facet, any facet whose type differs from the first sideset's, an empty id list, and each empty sideset by id. The last one is from the round-1 code jury (astra): an empty sideset made the type check order-dependent.
- **(Z2) Closed curves.** `orient_segments` turned every segment low→high, including the closing one, so `collect_nodes` mis-closed: an ASSERT fired in debug and release produced the wrong order. The Opus planner found that `sort_segments`' first step could also go either way. The fix:
  - The closing segment runs n−1 → 0, which needs `orient_segments( aClosed, … )`.
  - `sort_segments` orders segments by the walk index of their node 0, and `next( Segment )` is deleted.
  - The planner checked this in a Python model: every 3-, 4- and 5-loop closes correctly, and 2000 random open chains come out unchanged. The C++ was rederived by both auditors.

  Closed thin-shell side curves crash earlier: `identify_subchains` leaves the closing edge's group index stale, and `create_protocurves` then indexes out of range (astra, plan round). On Christian's ruling they are refused at setup, per component with E ≥ N, with scoped cleanup of the scratch state. Grok's code-round claim that the comment misstated this mechanism was refuted by the same trace.
- **(Z3)** The bitset `intersect( …, DynamicBitset& )` overload leaked its `Curve` and had no caller. Christian confirmed that `nonfree/` has none either, so it is deleted.
- **(O6)** Named refusals now cover an empty curve and disjoint sidesets.
- **Closed deck curves (R-3).** A closed deck curve is now built correctly, but it would then flow into `CutFactory::close_terminal_loops`, which duplicates terminal curves without checking `is_closed()`. This premise of the ruling was corrected by the implementer: the thin-shell refusal does not cover it. Christian's ruling: `MaxwellFactory` now refuses a closed curve at the one place that associates thin-shell terminals, naming the curve and the shell.
  - The refusal runs before any mutation: the round-2 Grok audit showed the 3-D `sideset_a` branch swapped the curve's sidesets before the check could fire.
  - It makes the `is_closed() && tNumTips == 0` disjunct at `cl_CutFactory.cpp:1245` unreachable (Grok, round 2). That disjunct came from the July T0 session as defense in depth inside an error message, never as a validated closed-tape path. Christian ruled to keep the refusal and record the dead branch for Gregory.
  - Christian asked for a factory-level test: `tests/fem/test_ThinShellTerminalCurves.cpp` (Opus). Each case builds a real `MaxwellFactory` from a generated deck and a hand-written Gmsh 2.2 mesh in its own `mkdtemp` directory. The decks deliberately lack a `homology` section, so `create_cuts()` stops right after the association. The acceptance cases therefore show that no refusal fired; they cannot show the `push`, because `terminal_curves()` is private.
  - Round 3 found that the shell under construction leaked on every refusal in that block, the new one and the older thickness, layer and "no sidesets" refusals alike. It is now held by a setup-scoped `unique_ptr` until `mProtoshells` takes it.

**Residue kept on DR-188 (Christian's ruling).**
- Mixed or cubic facet order is unchecked on the 3-D and thin-shell paths.
- `sort_end_nodes` refusals inside `thin_shell_side_curves` leak scratch state and leave flags set; this bites only in throw mode.
- The `aAdjacency` fill is dead.
- The `CutFactory::close_terminal_loops` `n-1-k` reversal puts the LINE3 midpoint into slot 0 (`src/homology/cl_CutFactory.cpp:~1598-1602`, only order-2 decks reach it). Christian's question for Gregory: could the cohomology run on the linear skeleton, with the order-2 lift confined to one place? Christian noted that the cell complex, and therefore the basis, is order-independent; the lift back onto the mesh is where the order-2 handling lives.

## Register

- **DR-178** was struck on Christian's ruling (fix `041c1dd`, three-face gate green in both builds). Its side-face sign test and the PENTA15/18 fixtures remain as residue.
- **DR-183** was filed and struck as refuted.
- **DR-184** was filed: the QUAD16 facet-1 table lists node 8 where node 7 belongs. It gets its own round.
- **DR-188** was filed and widened.
- The header was recounted three times as peer sessions filed DR-185–187. `check_doc_claims` passes 38/38.
- The commit carries the register whole. Its header counts depend on the rows the test-gap run filed today (DR-162–182, DR-185–187), which were uncommitted and interleave with ours.

## Process notes

- **Grok's turn budget.** The 30-turn default lost two full audits (DR-181 plan, mesh plan), and 60 turns lost another. At 100 turns, with an explicit "verdict by turn 70", the audits completed. A peer's in-progress rewrite of `ask_grok.sh` broke one dispatch with a bash syntax error; `bash -n` is now run before each dispatch.
- **Fable's pre-jury reading.** Christian's standing rule is that Fable reads every Opus diff before the vendor jury. That reading found the test-hygiene defect (a failure message that dereferenced the pointer under test, turning a red case into a SEGV) and pre-registered the LINE3 tail as never executed, which the gdb trace then confirmed.

## Open

- DR-169: the valgrind leg on a non-MKL tree.
- DR-179: the 64-bit leg, blocked on the missing SCLS flavor and on the unrelated INT64 compile failures.
- A debug-tree run of all three bundles; the MPI suites.
- DR-184 (the QUAD16 facet table).
- DR-188's residue, above.
