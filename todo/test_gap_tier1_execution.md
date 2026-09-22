# Test-gap Tier 1 execution (overnight run)

**Date:** 2026-09-21
**Purpose:** Work off Tier 1 of `tmp/devlog_sweep/triage.md` v2 (22 rows) with the
`code-with-jury` loop, unattended: for each row re-read the mechanism against local `devel`,
plan, plan jury, implement the test (and the local fix a red test exposes), code jury, gate in
the debug tree, commit own files. Findings that need a ruling or a wider fix go to
`todo/debt_register.md` instead of an edit.
**Module:** `tests/` (all suites); `src/` only for jury-confirmed local fixes
**AIs involved:** Claude Fable 5.1 (loop), Codex (`gpt-6-astra`/high on critical rows,
`gpt-5.6-sol`/medium plan jury and `gpt-5.6-terra`/high code jury otherwise), Grok `grok-4.7`/high
**Status:** ✅ RUN CLOSED 2026-09-22 (afternoon wall clock after the jump) — 21 of 22 rows committed on `devel` (`b5a18cf` … `f2d34c2`, plus `8f317e8` comment), T-08 carried with plan v3 rulings; DR-161..182 filed, DR-178 fixed and closed; every src fix consumer-swept (8/8 LOW); the quiz answered O1–O21 and DR-163/168/170/171/172; follow-ups in `todo/test_gap_followups.md`; Gregory's items in `devlog/dl20260922_homology_report_for_gregory.md`; final `make check-fast` in both trees at close

> **Standing approvals (Christian, 2026-09-21, recorded in memory `test-gap-overnight-run`):**
> blanket approval for tests and for jury-confirmed local `src/` fixes; I run `make`/`ctest`
> in `cmake-build-debug` per row and confirm in `build`, full `make check` when it makes sense;
> one commit per row on `devel`, own files only; peer WIP (`git status` at start) never staged;
> on a red gate not fixed in one attempt, a vendor quota failure or an open ruling: file a
> `DR-` row and move on. Cohomology core stays closed.
>
> **Scope guards:** Tier 1 rows only, in `triage.md` order. A row's "live in tree" claim is
> re-read against local `devel` first (the vendors read published `main` at `ef6d47f2`).
> Register and README edits are appended but left unstaged (those files carry peer edits).

---

## 1. Environment facts (read 2026-09-21)

| fact | value |
|---|---|
| debug tree | `cmake-build-debug` (`USE_DEBUG=ON`, MKL ON, gasmodels OFF, tests ON) |
| release tree | `build` (`USE_DEBUG=OFF`, MKL ON, gasmodels OFF, tests ON) |
| cores | 1 (serial builds; Tier 2 MPI runs oversubscribed) |
| next free register id | DR-183 (DR-161..182 filed this run, appended unstaged to the peer-edited register; its `[P]`/`[W]` recount line is owed) |
| peer WIP not to stage | see `git status` list of 2026-09-21 23:58 (18 modified, 11 untracked); notably `tests/comm/test_CommMPI.cpp`, `todo/debt_register.md`, `todo/README.md`, `devlog/README.md`; **plus** `src/physics/materials/cl_Material_Silver.cpp` and `.../doc/material_property_sources.md`, modified 00:06 by a peer (Gerritsen 1952→1956 citation year), not by this run |

## 2. Row ledger

One entry per Tier 1 row, in `triage.md` order. Each carries the re-read verdict, the plan
slug, the jury outcomes, the gate, the commit, and the residue.

| # | id | re-read on `devel` | plan/jury slug | gate | commit | residue |
|---|---|---|---|---|---|---|
| 1 | T-02 | live: count/fill predicate mismatch + unguarded serial layer lookup | `review_tgap_t02` (plan v3; code jury astra + Grok re-run: approve, comments/moves fixed) | debug fem 187/187, release 185/186; mutation red 2 of 3 | `b5a18cf` | connector loop → DR-162 |
| 2 | T-01 | live: HEX20 case sizes 9 and reads `node( 21 )`; QUAD9 bottom midside reads `node( 5 )` | `review_tgap_t01` (plan + code juries approve, 5 P3 nits applied) | debug: red 3/15 → green; mesh 74/74 both trees after the nits | `c41c1c0` | O1 → DR-161; release cannot bound `node( 21 )` (R32), G4 on the debug tree |
| 3 | T-26 | three live defects (DR-163, DR-164, R5 fixed at six sites) | `review_tgap_t26` (code jury terra + Grok: approve with edits) | physics 34/34 both trees; R5 mutation SEGV; `DISABLED_` pins red as pre-registered | `99ac9d0` | DR-170 (`n = 4` tier), O5 (iron Kohler tangent) |
| 4 | T-12 | no live defect; test pins the lazy-inverse shape | `review_tgap_t12` (code jury approve; comment reworded) | debug 283/283, release 263/263; G4 probe: lazy 2×3 vs wrapper 3×2 | `90c7bfe` | storage order of `inv( product )` unpinned (residue) |
| 5 | T-11 | contract: same-class width converts, class refused | `review_tgap_t11` (code jury: Cell/Matrix paths + all-ones prefill applied) | io 57/57 both trees | `7d3dfff` | handle leak → DR-165; philosophy exception → O6; G4 mutation not run (header library-wide), reviewed by trace |
| 6 | T-13 | counters: sources uint16, dofs uint8, node/edge/face/element capacity uint8, facet uint16, all `BELFEM_ERROR`-guarded; graph counter uint32 assert-only (DR-166) | `review_tgap_t13` (code jury astra + Grok: fatal counts, RAII owner, fills to capacity applied) | debug mesh 95/95 + math 420/420 (round 7); mutants A, B, C red (B redone serially) | `e2f66ab` | unfilled source container leaks in `~Basis` → DR-169 |
| 7 | T-10 | no live defect; R1 reads `rho_i_custom` | `review_tgap_t10` (code jury astra + Grok: approve with edits) | fem 187/187 + physics 34/34 debug, 185/186 + 34/34 release; integer-truncation mutants red (the bare `abs` spelling is green under these flags) | `81ffe4e` | ferromagnet `rho_mag` outside the RRR identity → O7 |
| 8 | T-24 | code IS symmetric; unpinned coefficients now pinned | `review_tgap_t24` (plan v2; code jury astra + Grok approve) | fem 187/187 debug, 185/186 release | `73f8aa4` | DR-168 `[F]` (jump vs average, α floor three numbers); unequal-h fixture residue |
| 9 | T-32 | live R32 defect at `:743` | `review_tgap_t32` (plan v3; code jury astra + Grok: approve with edits) | mesh 81/81 both trees; mutation red by text | `4871e04` | DR-167 (leaks on refusal), O2 (off-axis rotation contract); the match count is hits, not distinct partners (INC-009 geometry invisible to a count) |
| 10 | T-25 | every public cp/λ/ρ is a 4 K spline of a PROTECTED raw law except Cu/Ag cp, YBCO ρ, Fe/Ni ρ (raw); Magnesia has no ρ; the concrete classes' switch members are private | `review_tgap_t25` (Grok refuted v1: probe does not compile, 0 K logs, exponent bound not a law property; plan v2 written; Codex re-run queued 03:37) | debug + release physics 100/101 (magnesia has no rho); the two windows red as disabled pins | `62b20c2` | DR-173, DR-174; O8 (envelope thresholds are a project threshold, ruling on the values); raw sweeps carry no step envelope |
| 11 | T-51 | deck path hardcodes `aBuildTables = true` → S1 seam landed; K-1/K-2/Q3/Q4 → DR-171 | `review_tgap_t51` (plan v2 approved by terra + Grok) | debug physics 98/101 (3 skips), release filtered green | `2076ec2` | DR-171 `[F]`; DR-175 (refusal-path leaks, both auditors); flag not forwarded for Hastelloy/YBCO/magnesia/curve (no current deck changes) |
| 12 | T-23 | derivation +1 (plan agent, Grok, Codex), row/core comment −1; the test's `gInnerRingSense = 1` is a calibration, disagreement → Gregory (O10) | `review_tgap_t23` (plan v2 approved) | debug homology 21/21, release green | `5834e8b` | O10 (sense vs core comment → Gregory); DR-177 (`getSimplicesMap` declared type, closed core → Gregory) |
| 13 | T-33 | contract exists; the row's "throws in debug" is wrong (`BELFEM_ERROR`); M5 not dead arithmetic (10 ε guard) | `review_tgap_t33` (plan v2 approved by astra + Grok) | debug fem 197/197 (round 8), release filtered green; M1-M5 red by construction (M5 via the near-unit case) | `8124a70` | DR-172 `[P]`; untested: null/hanging source refusal, duplicate bypass, weights between one ulp and 0.5 |
| 14 | T-34 | face key = sorted corner triple; third owner silently dropped (O9) | `review_tgap_t34` (plan v2 approved by astra + Grok) | debug mesh 95/95, release filtered green | `abd598f` | O9 (third owner); DR-179 (`BELFEM_INT64` sentinel mismatch, both auditors); quad key drops the fourth corner (coverage residue) |
| 15 | T-60 | no live defect on the normals; NEW derived findings: 2-D QUAD4TS ghost master normals point inward on the production CW winding (Q2), PENTA side-face `dS` short by ×2 (Codex, derivation) | `review_tgap_t60` (plan v2 after astra + Grok; implementation on Opus running with Christian's TS-numbering notes) | debug fem 198/198, release 196 + 1 skip; three side cases mutated separately, each red on its own face, restore green; consumer sweep LOW, O20 ruled none | `041c1dd` | DR-178 FIXED (closes on the gate); DR-176 (`create_unity_nodes` HEX8TS arm); Q2 (2-D ghost master normals inward on the CW winding) unpinned; side-face sign untested; the check compiles only under USE_MAXWELL |
| 16 | T-36 | domain-type skip in `update_facet_nodes`; release write unobservable from `tests/mesh` | `review_tgap_t36` (plan v2 approved) | debug mesh 95/95, release filtered green | `892f862` | O11 |
| 17 | T-38 | **live defect** at `cl_CurveFactory.cpp:1016/:1096`, R1 landed | `review_tgap_t38` (plan v2 approved) | debug mesh 95/95, release filtered green; G4: R1 removed → exactly the three refusal tests red with the pancake 'Row index 4294967295' crash, restore green | `f58524c` | `intersect( id, id, bitset )` never unflags, no open-tree caller (Grok); self-loop accepted (pre-existing) |
| 18 | T-44 | **row contract (1) not in the tree**; R1 landed in the open file | `review_tgap_t44` (plan v2 approved) | debug homology 21/21 (round 7), release green; G4: R1 else-branch removed → exactly the two ±2 refusal tests red, restore green | `842213c` | periodic bitset branch uncovered; the compacted production index space is not what the fixture uses; O18 to Gregory; `compute_edge_bitsets` dead-work by-catch |
| 19 | T-08 | three LIVE parser defects by reading (comment-eaten `;` dropped silently; unsigned depth decrement; no end-of-scan check), C4 fixed and unpinned | `review_tgap_t08` (round 1 terra + Grok on v1 hard errors; Christian ruled warnings; v2 by a Fable agent; round 2 terra/xhigh + Grok/xhigh: two load-bearing refutations — the `;`-statement header is a live legacy form, EOF close re-parents siblings — both ruled) | no gate (not implemented) | | CARRIED to the next session by Christian's ruling; plan v3 list appended to `review_tgap_t08_plan.md` (header-and-statement, close before the next header, InputFile forwarders, two ctest commands, Python must not skip silently) |
| 20 | N-17 | all three CONTRACT (integer fill typed since `8465284`; bitset copy assignment copies the index but nothing pins it; `collect_node_coords` self-sizing, protected, callerless) | `review_tgap_n17` (plan jury terra + Grok: approved with corrections; code jury: Codex approve, Grok's failure claim refuted — it read a scratch mutant, P5) | debug linalg 300, containers 606, fem 201; release 280, 570, 199+1; M5 red on the three IwgNodeCoords cases, the index-copy deletion red on IndexSurvivesCopyAssignment | `c2490a4` | optional IWG_Maxwell sentinel case dropped (no C++ Maxwell macro in the fem suite; would need its own Maxwell-listed file); `allocate_work_matrices` and `collect_node_coords` remain callerless; stray `.rej` removed |
| 21 | N-18 | all five CONTRACT (GEMM NT/TN untested; `eigen_sym` copies + 'U' on both backends; Bezier third derivative; Database axis steps; HEX27 entry 23) | `review_tgap_n18` (plan jury terra + Grok: approved with the mutation map corrected; code jury: Codex changes → applied, Grok approve with two comment fixes) | debug linalg 300, math 424, physics 111+1; release 280, 358, 111+1; M3, M7, M10b, M11 red on exactly their predicted tests | `a05c080` | M1/M2 debug-only (not run); M9/M9b/M12/M13 not run (budget); non-affine FD margins are an off-line measurement; Blaze complex-double padding depends on SIMD width |
| 22 | N-19 | ONE LIVE: two-argument user `jc`/`n` registration ignores its dependency flags (silent argument transposition, `cl_Material_UserDefined.cpp:101-126`); three CONTRACT pins | `review_tgap_n19` (plan jury terra + Grok: approved with the dlopen leak, citation and sentinel corrections; code jury: Codex blocker refuted (shared CMake), Grok three nits applied) | debug + release physics 111/112 + 1 skip incl. the plugin build; R3 refusals fire; R5 measured 4.84e-4 on both backends, band 0.5 % | `f2d34c2` | O21 ruled and applied (template prose corrected; the one in-tree plugin uses the three-argument order); handle-leak path has no gate; `debye.f90:234` implicitly SAVEd logical (register candidate); `debye.f90:25` default-real `c` literal |

Order note: T-01 is taken before T-02 because its fix and test are self-contained (mesh suite,
hand-built elements), which makes it the cheapest first pass through the whole loop; T-02 needs
a two-shell kernel fixture.

## 3. Defects found during the run (Dn)

- **D1 (T-01)** `get_bottom_nodes` HEX20 sized 9 and read `node( 21 )`; QUAD9/QUAD9TS read the right-edge midside. Fixed in `src/mesh/meshtools.cpp`, red-then-green.
- **D2 (T-02)** collector count/fill predicate mismatch; serial linker unguarded layer lookup. Fixed (seam + guard), gate running.
- **D3 (T-26)** piecewise `dJ` tangent inconsistent with the raw residual below `j1` → DR-163.
- **D4 (T-26)** Kim `deval_dB` missing → DR-164.
- **D5 (T-26)** 3-argument power-law overloads dereference a null `mNFunction` → R5 local fix.
- **D6 (T-11)** HDF5 loaders leak handles on the class-mismatch throw → DR-165.
- **D8 (T-26, test)** the plain h/(h/2) agreement at 1e-6 refused three smooth points (n = 27 curvature, spline `rho_n( T )`); the oracle now uses the Richardson-extrapolated value with a 3e-5 agreement gate. Test-side.
- **D10 (T-25)** lead's public cp negative below 3 K → DR-173; **D11 (T-25)** YBCO raw Callaway λ NaN below ~1 K → DR-174.
- **D7 (T-01, test)** gtest parameter names collided: `to_string( ElementType )` folds TS/TB variants; fixed in the test's name generator.

## 4. Open questions (On)

- **O1 (T-01, from Grok)** `get_top_nodes` accepts HEX64, QUAD8, QUAD16 through `geometry_type`
  while `get_bottom_nodes` rejects them (`meshtools.cpp:1389-1411`); the two helpers do not
  share one supported-type list. Outside T-01; a later row or a register note.
  **RULED 2026-09-22 03:47 (quiz):** one shared supported-type list for `get_top_nodes` / `get_bottom_nodes` (refusing TRI3/TET4, HEX20 per O14, and the higher-order hexes and quads alike); folded into D13b (src, jury).
- **O2 (T-02)** whether a strict subset of thin-shell layer blocks is a supported selection or
  should be refused at `create_blocks`; R1 makes it work consistently. → DR-162.
  **RULED 2026-09-22 03:35 (quiz):** supported: guard the serial side-connector loop (`cl_FEM_DofMgr_BlockData.cpp:481`, `:501`) with `key_exists` and add a side-connector fixture — queued as D12 (plan + jury), DR-162 → ruled.
- **O3 (T-26)** `mRhoMin` has no setter, so the "both sides of the floor" clause of the row is
  not testable without a `src/` change. Left out of the test; register note deferred.
  **RULED 2026-09-22 03:47 (quiz):** a setter is acceptable but `mRhoMin` is an EXPERT setting: zero is a fine default; HTS values are of order 1e-17 and 1e-16 is already too big; a user would not normally touch it. Queued low as D15 (setter + two-sided test; if it becomes a deck key, the input contract rule fires).
- **O4 (depth, all sol rows)** both auditors of T-12 and T-11 flag that `gpt-5.6-sol`/medium sits
  on no row of protocol §9.1 for a round-1 plan audit. This run follows Christian's 2026-09-21
  instruction ("astra for the highly critical tests, terra and sol otherwise") as a run-scoped
  ruling; whether §9.1 should gain a "test-only plan" row is Christian's call.
  **RULED 2026-09-22 03:42 (quiz):** terra/high for plan and code rounds on non-critical modules, astra/high on the critical list; sol retired from plan rounds. Standing table (memory updated; protocol §9.1 edit owed to a session that owns `doc/ai_collaboration_protocol.md`, which is peer WIP tonight).
- **O5 (T-26, from the gate)** iron's Kohler `drhodT( T, B, β )` at 7.3 K, 5 T differs from a
  Richardson central difference by 1.7e-5 relative (the accessor is itself a fixed-step FD at
  0.001 K). The test holds the field-dependent overloads to 1e-4 for that reason; whether a
  1e-5-class truncation on a Newton tangent matters is a physics call.
  **RULED 2026-09-22 03:47 (quiz):** "there is a global parameter for gDeltaT or something like that; check whether using it would cause a regression; hard-coded steps are not good practice" — the Kohler accessor's hard-coded 0.001 K step is to be replaced by the framework's global FD step after a regression check (D16: locate the global, compare the two derivatives on the T-26 fixture, plan + jury). Checked: the Kohler derivative `Metal::drhodT_kohler` (`cl_Material_Metal.cpp:758-768`) already steps by the framework global `gFinDiffDeltaT` — the hard-coded value lives in that global's definition, not in the accessor. D16 therefore reduces to: is the global's value the right step (the 1.7e-5 residual is its truncation), and should the T-26 test read the same global instead of its own 1e-4 tolerance rationale. No src change until that is answered; the test's comment is corrected to name the global.
  **RULED (follow-up, 03:52):** nothing to change; the T-26 test comment names `gFinDiffDeltaT` as the source of its 1e-4 tolerance.
- **O6 (T-11, from Codex)** the HDF5 loaders convert same-class widths on read
  (`hdf5_tools.hpp:185-202`), an I/O exception to the philosophy's "no automatic numerical
  conversions" (`coding_philosophy.md:103-111`). Document the bounded exception there, or not.
  **RULED 2026-09-22 03:47 (quiz):** document the bounded I/O exception in `doc/coding_philosophy.md` (same-class widths convert on read, never across classes); queued as a doc edit for a session that owns that file (peer WIP tonight).
- **O7 (T-10, from Grok)** for iron and nickel `Ferromagnetic::rho_custom` adds `rho_mag( T )`
  that `Metal::set_RRR` never includes, so the stored `rho_0` satisfies
  `( rho_i( T_ref ) + rho_0 ) / rho_0 = RRR` and not `ρ( T_ref ) / ρ( 0 K ) = RRR` as the
  header comment (`cl_Material_Metal.hpp:117-123`) states. Which definition is meant for the
  ferromagnets is a physics call; the test pins the code's residual and says so.
  **RULED 2026-09-22 03:52 (quiz):** the header is right: RRR = ρ(293 K)/ρ(0 K) with the magnetic term included; `set_RRR` for the ferromagnets must solve with `rho_mag` (src, jury, physics check against the RRR data source; this changes iron/nickel resistivity values → consumer inventory first). Queued as D19.
- **O9 (T-34, from both auditors)** three tetrahedra on one triangle are a legal non-manifold
  input; `FaceFactory` silently drops the third owner (`cl_FaceFactory.cpp:359-373`). If a guard
  is wanted it is a `BELFEM_ERROR` on the slave slot at the entry of the `aMasterIDs != gNoID`
  arm; the two-TET fixture cannot show it. Christian's call whether to add it.
  **RULED 2026-09-22 03:52 (quiz):** leave it; non-manifold meshes are not an accepted input.
- **O10 (T-23)** the first debug run calibrates the sign `s` of the CCW inner-ring pairing on
  today's tree; both auditors and the plan derive +1, the row and the core comment say −1. The
  calibration and the derivation-vs-comment disagreement go to Gregory.
- **O11 (T-36, from both auditors)** `Mesh::update_facet_nodes()` re-derives a plain facet's
  nodes with no check that the master's facet node count equals the wrapper's (`cl_Mesh.cpp:648`);
  a one-line `BELFEM_ERROR` there is allowed by the tier rule and would catch the short-list
  case the `insert_node` assert misses in release. Whether to add it is Christian's call; the
  row stays test-only.
  **RULED 2026-09-22 03:52 (quiz):** add a `BELFEM_ASSERT` (debug tier, not ERROR) on the facet node count in `update_facet_nodes`, and a regression test that runs in debug mode. Queued as D18.
- **O8 (T-25)** the no-jump envelope on the 0.5 K grid uses one named threshold pair
  (`p = 5` local exponent, `c = 0.01`); Grok showed neither is a property of the fitted laws
  (Fe/Ni Bloch-Grüneisen 4.5, uncapped Bézier/log-polynomial pieces). The test prints the
  measured maximum exponent and asserts the pair; whether the pair is the contract or only a
  smoke bound is Christian's call.
  **RULED 2026-09-22 03:42 (quiz):** smoke bound; keep the pair; a rejected fit is a ruling for the property owner.
- **O12 (T-08)** the deck parser refusals (orphan statement whose `;` a comment ate, stray `}`,
  unclosed section) are a behaviour change at the boundary every run crosses; decks that "parse
  badly but parse" today will stop with the file and line named. Recommended: hard error; the
  one recorded reason for deferring (release day) has lapsed. The plan's shipped-deck scan is a
  Python mimic and did not see `nonfree/`.
  **RULED 2026-09-22 03:35 (quiz):** the input format is deliberately failure tolerant: a new deck gets a WARNING naming file and line and proceeds when the intent is clear (orphan statement → taken as the statement; stray `}` → ignored; unclosed section → closed at end of file). Example decks must be flawless; a Python deck checker (the plan's mimic, promoted to `scripts/check_input_decks.py`) must ERROR on any of the three in shipped decks. T-08 plan → v2 after the running jury lands.
- **O13 (T-01 sweep)** `ThinShellFactory` builds QUAD9TS / PENTA18TS layers in full while
  `cl_MaxwellFactory.cpp:1711` refuses any shell whose facet is not TRI3/LINE2: is second-order
  thin-shell support meant to be live (then the T-01 fix becomes load-bearing and the rest of the
  second-order hang path needs its own audit) or are those branches dead?
  **RULED 2026-09-22 03:38 (quiz):** dead for now, mark it — a comment at the QUAD9TS/PENTA18TS extrusion branches of `ThinShellFactory` says the layers are built but the Maxwell factory does not assemble them (comment-only, no round). Queued as D13a.
- **O14 (T-01 sweep)** no thin-shell path can produce a HEX20 layer; should HEX20 join the
  refused set instead of keeping a `get_bottom_nodes` case?
  **RULED 2026-09-22 03:38 (quiz):** refuse HEX20 — move it to the refused set in `get_bottom_nodes`/`get_top_nodes` beside TRI3/TET4 (src, jury). Queued as D13b; the T-01 HEX20 test case then pins the refusal.
- **O15 (T-26 sweep)** the T-less power-law family (eight bodies, five still callerless) has now
  produced one defect; delete it or forward it to the T-carrying form with `gTbulk`? By-catch:
  `Material::n( normB, angle[, T] )` lacks the constants fallback `jc()` has.
  **RULED 2026-09-22 03:38 (quiz):** keep the T-less family as is; the fix and the tests cover it.
- **O16 (T-38 sweep)** `intersect( id, id, DynamicBitset & )` never unflags before sorting and has
  no caller: add the two-line unflag (jury) or delete the overload.
  **RULED 2026-09-22 03:38 (quiz):** delete the callerless `intersect( id, id, DynamicBitset & )` overload (src, jury). Queued as D14.
- **O17 (T-38 sweep)** keep the refusal of a multi-chain intersection over a silent split (a curve
  id is addressed by the deck); an explicit list-valued entry would be the contract change if the
  pancake geometry needs it.
  **RULED 2026-09-22 03:42 (quiz):** "please check the examples, this might be a valid case" — a per-deck check of every `curves` entry against its `.geo` is running (read-only agent, `tmp/ai_exchange/consumer_impact_t38_curves.md`); the definitive gate is `belfem` to the cut stage on each curve deck with a fresh gmsh mesh, offered to Christian. If a shipped intersection is legitimately two chains, the refusal must become a split (O17 reopens as a contract change). **Checked (04:00, `tmp/ai_exchange/consumer_impact_t38_curves.md`):** all ~550 `curves` entries of the 14 shipped decks resolve to one open chain each (zero loops, zero multi-chain, zero undecidable); pancake ships the fixed form its generator now emits; the decks already express a conductor crossing a plane twice as one curve id per chain grouped in the boundary condition (`corc_periodic_bc` `[1,2]`, `racetrack_usermat` separate ids) — the refusal stands. Runtime confirmation still owed on `pancake` (fixed deck never run) and `corc_twolayer` (python mesh); a `belfem` run past `create_curves` on a fresh gmsh mesh is the gate.
- **O18 (T-44 sweep, Gregory)** write `clean()`'s unit guarantee as a contract beside it and as a
  `@pre` on `CutData` (comment-only); the new error says "Cut N" (1-based) where `clean_spfa` says
  "generator N" (0-based); `compute_edge_bitsets()` fills two bitsets nothing reads.
- **O19 (T-32 sweep)** no opt-out for the periodic node-match refusal exists and none is
  recommended (a wider radius turns a loud refusal into the quiet double-hit mis-pairing).
  **RULED 2026-09-22 03:42 (quiz):** no opt-out.
- **O20 (DR-178 sweep)** production reach is Christian's own decks: a hand-built prism tape whose
  narrow edge face sits in a symmetry or background sideset changes by 2× on that term;
  PENTA15/18/18TS models share the callback and have no fixture.
  **RULED 2026-09-22 03:35 (quiz):** no such decks. DR-178 commits as a pure correction and closes on the gate + mutation evidence.
- **O21 (N-19)** R1 refuses a two-argument user `jc`/`n` plugin whose declared dependencies are
  not `( normB, angleNxB )`: a plugin that "works" by accident today will be refused — consumer
  inventory of the shipped plugin templates and example plugins before it lands.

### Register rulings from the quiz (04:15 PDT)

- **DR-168 [F]** Christian: "I think the jump is intended, these are layers of different materials
  with different resistivities. The magnesia layer is actually a phi layer without any electric
  conductivity. Let's double check if my reasoning seems sound to you." Claude's assessment is in
  the report (the physics argument for a flux-type interface term is sound; whether the kernel's
  jump of OUTWARD slopes equals the tangential-E continuity condition or the SWIP average is a
  sign-and-factor question only the unequal-ρ fixture can settle) — the fixture stays queued as
  D20; the row keeps `[F]` with the ruling recorded.
- **DR-170 [W]** Christian: avoid `BELFEM_ERROR` in hot paths; `n = 4` is legal and fractions such
  as 4.5 exist; "I don't understand why this is even a DR." → the row is ONLY about the abort tier
  in `drho_piecewise_dJ` at the degenerate knot, never about the legality of `n = 4`; reworded,
  and the fix (the same frozen-knot fallback as `drho_piecewise_dT`) queued as D21 (src, jury).
- **DR-171 [F]** Christian: contradictions are refused, redundancy that agrees is allowed
  (`builtin : pb40sn60rrr5 ; RRR : 5 ;` is fine, `RRR : 100` beside it is refused; `builtin`
  beside `curve` is refused). Queued as D22 (src, jury, `input_file_reference.md` +
  `input_schema.yaml` in the same session).
- **DR-172 [P]** Christian: "I don't know. Is there a use case? Refusing sounds reasonable, but it
  might cause a regression within the shipped examples." → consumer inventory first (do any
  shipped decks impose a Dirichlet condition on a multi-source hanging dof?), then D23 decides
  between refusal and acceptance.

- **DR-163 [P]** Christian: "I would argue that parallel is right, but this might need a jury
  because I am not sure." → a physics round (Codex + Grok, literature: Messe et al. 2023 §2.6,
  Rhyner 1993) on which side of the piecewise law carries the parallel combination BEFORE any fix;
  queued as D24. Consumer inventory applies: changing `rho_piecewise` changes every piecewise
  residual.
- **O10 / O18 → Gregory** by a report file `devlog/dl20260922_homology_report_for_gregory.md`
  (O10 sense calibration vs the core comment, O18 `clean()` contract + numbering + dead bitsets,
  DR-177 declared type) written at close-out; Christian forwards it.
- **T-08 S4** confirmed: a headerless `{` is the one refusal at the deck boundary.
- **D-queue scope**: finish the 22 rows tonight; D12–D24 become a registered todo plan for the
  next session (`todo/test_gap_followups.md`).

## 5. Process incidents

  **RULED 2026-09-22 03:35 (quiz):** land the refusal; the shipped example plugins and templates must be corrected in the same row (consumer inventory first).
- **P1 (commit sweep, 01:28 PDT).** The first T-02 commit swept six files a peer session had
  STAGED in the shared index (`.claude/commands/cross-review.md`, two `.claude/skills/*/SKILL.md`,
  `CLAUDE.md`, `doc/ai_collaboration_protocol.md`, `devlog/dl20260922_workflow_skills.md`).
  `git add <own files>` does not protect against a non-empty index. Repaired the same minute:
  the four row commits were rewritten from `7479cb2` with explicit own-file staging only
  (`b5a18cf`, `81ffe4e`, `73f8aa4`, `99ac9d0`); the peer's six files are back in the working
  tree as unstaged modifications with their content intact, but their staging is lost — the
  peer has to `git add` them again. Nothing was pushed. Rule from here: `git diff --cached
  --stat` must be empty before every `git add`, and the staged list is read back before commit.

- **P2 (Codex quota, 01:46 PDT).** The Codex CLI reports "You've hit your usage limit … try again
  at 3:30 AM" from about 01:45; the T-23 plan jury lost its Codex voice (Grok delivered). Per
  the standing rule the missing voice is reported here, not substituted: the T-23 plan is
  re-dispatched to Codex once Christian topped up the credits at 02:00 (the sleeping 03:35 jobs
  were cancelled); T-25, T-51, T-33 and T-34 got their Codex voice the same way. Rows whose code round would need
  Codex before then (T-13 is already twice juried; nothing else is at the code stage) are not
  blocked. Christian may top up credits.

- **P3 (two makes in one tree, ~02:28 PDT).** The second T-13 mutant chain was launched while
  the first chain's restore build was still running in `cmake-build-debug` (the completion
  notice had not arrived; the log tail was misread). The two `make`s raced for ~40 minutes and
  the restored `test_math` binary still carried the narrowed counter. Repaired: second chain
  stopped, its scratch edit restored by hand (`git diff` clean on all three files), the three
  files touched, a serial clean rebuild + rerun (round 7) queued as the gate. Rule from here:
  one `make` per tree, verified with `pgrep -fa "make test_"` before every launch.

- **P4 (Codex usage window exhausted again, 02:42 PDT).** Every code round of the second wave
  (T-33, T-51, T-34, T-23, T-36, T-38) came back `Audit FAILED` on the Codex side ("try again at
  3:30 AM"); the Grok halves landed. Reported; Christian bumped the credits at 02:58 and all six
  Codex halves were re-run from `tmp/ai_exchange/tgap_codex_rerun.sh` (same prompt shape as
  `cross_review.sh --jury`, one prompt per row) and appended to their threads. Side quest from
  the same message: `.claude/scripts/ask_copilot.sh` written as a third seat, smoke-tested
  (read works, a write attempt is refused by the CLI's permission rules); on this account only
  `auto` (→ gpt-5.6-luna) is policy-enabled, every stronger model is refused before billing until
  the github.com policy toggles are flipped, so it is a luna-class voice for now. The T-33 Grok
  half also hit its 60-turn budget re-deriving the fixture; re-run at 120 turns with a brief
  that named the three functions to read, landed.

- **P5 (a jury read a scratch mutant, ~14:35 PDT).** The N-17 code jury ran while the N-17 and
  N-18 G4 mutation chains were applying scratch edits to `src/`; Grok read
  `cl_DynamicBitset.cpp` with the index-copy line deleted (and saw the `cl_Database.hpp` marker of
  the other chain) and reported the new test as failing against the tree. The mutants were restored
  byte-identical minutes later; the finding is REFUTED as a tree defect and the audit is otherwise
  usable. Rule from here: no G4 chain while a jury is reading the same tree — run mutations after
  the round lands, or jury a frozen diff whose src files the chain does not touch.

- **Retrospective consumer-impact sweep (Christian, ~03:20 PDT).** Rule: before a fix that can
  change a value or reject an input lands, inventory its consumers and the shipped decks that
  reach them (memory `feedback_consumer_impact_before_behavior_change`). Applied retroactively to
  every src fix of the run: T-01, T-02, T-26, T-32, T-38, T-44, DR-178 (one read-only Opus agent
  each, reports in `tmp/ai_exchange/consumer_impact_<id>.md`) and T-51 (done by hand: one
  production caller, default unchanged, LOW). Verdicts and the reach rulings they leave are
  collected in §4 of this file once the agents return.

## Consumer-impact sweep (retrospective, Christian's rule)

| fix | verdict | reach | ruling left |
|---|---|---|---|
| T-01 `c41c1c0` (QUAD9 midside, HEX20 size) | LOW (high) | both branches unreachable from every shipped deck and production path: `cl_MaxwellFactory.cpp:1711` refuses any thin shell whose facet is not TRI3/LINE2, so no QUAD9TS shell reaches the four callers (`:1856`, `:1992`, `:1930`, `:2099`); no shipped `.geo`/`.msh`/deck raises the order; `OrderConverter` is never instantiated. Old HEX20 branch was an over-read never dereferenced | O13: is second-order thin-shell support meant to be live (the shell factory builds QUAD9TS/PENTA18TS, the Maxwell factory refuses them)? O14: should HEX20 be in the refused set rather than a case? |
| T-51 `2076ec2` (factory flag) | LOW (high) | one production caller, default unchanged | none |
| DR-178 (penta side dS) | LOW (high) | a prism volume block is refused at calculator init (only `PENTA6TS` passes, `cl_FEM_Calculator.cpp:1436-1448`); thin-shell side faces are mastered only by side connectors (`dV`) and the unassembled periodic sidesets; ghosts are the caps; `TMatrix` is tet-only; the thermal film BC is a `dS` consumer but `StaticHeatConduction` is never requested by any deck; no shipped mesh has prisms and no `.geo` extrudes with `Layers{}` | O20: Christian's own decks — a hand-built prism tape whose narrow edge face sits in a symmetry or background sideset changes by 2× on that term; PENTA15/18/18TS models share the callback and have no fixture. By-catch filed: DR-181 (`background_phi` dS at point 0), DR-182 (`Periodic` stamp removes a kernel-bearing sideset) |
| T-32 `4871e04` (periodicity refusal) | LOW (high) | premise corrected: the old release did not continue with unmatched nodes, it left null targets and segfaulted in `create_facet_map` / `match_edges` (`:826-828`, `:891`) — the change is a diagnostic, not a new refusal class; the check sits behind two pre-existing 1e-9 gates (`select_sidesets:1464`, `map_facets:532-538`, compile-time `BELFEM_MESH_EPSILON`, no deck key); all eight periodic decks are node-for-node images (`Periodic Surface` in the five gmsh decks, one triangulation written twice at 12 decimals of a mm in the three python decks); two observed to pass the cut stage on 2026-09-19 | O19: no opt-out today and none recommended (a wider radius makes the double-hit mis-pairing at `:729-737` a quiet wrong answer); cosmetic: the message's `m` unit assumes the scaled mesh |
| T-44 (cut coefficient refusal) | LOW (high) | `Cohomology::clean()` already carries two release-active unit gates (`cl_Cohomology.cpp:605-613`, `:883-890`) after the one non-unimodular producer (Step 7 of `updatekGeneratorsFromHomology`, followed by `clean()` at `cl_CutFactory.cpp:350-352`), so the refusal is a third line of defence; no record of a non-unit reaching `collect_coefficients` on a shipped deck (non-units occur upstream: corc_twolayer, helix, all rectified inside `clean()`). The old path was worse than the doc said: 3-D aborted with the wrong diagnosis ("not sign-coherent"), 2-D release ran a cut of the wrong sign silently | O18 (Gregory): write `clean()`'s unit guarantee as a contract next to it and as a `@pre` on `CutData` (comment-only); the new error says "Cut N" (1-based) where `clean_spfa` says "generator N" (0-based); by-catch: `compute_edge_bitsets()` fills `mCohomologyEdgesPlus/Minus` that nothing reads (`cl_CutProcessor.cpp:219`, `:223`) — dead work, medium-high |
| T-38 `f58524c` (curve refusal) | LOW (high) | the cut factory's `thin_shell_side_curves` runs `identify_subchains` first (`cl_CurveFactory.cpp:453-500`), one component per curve → structurally immune (a seam-broken tape side yields several curves, not a refusal); the reachable callers are the deck `curves` paths (`cl_MaxwellFactory.cpp:608`, `:613`); 14 of 20 example decks build curves, all read as single chains, and every one has been run in a debug tree where the old code aborted on a multi-chain curve. Smoke order if wanted: pancake, rutherford, gantry/undulator2d, corc_periodic_bc | O16: `intersect( id, id, DynamicBitset & )` never unflags and has no caller — add the unflag (two lines, jury) or delete the overload; O17: keep the refusal over a silent split (a curve id is addressed by the deck); by-catch DR candidate: `from_2d_sidesets` passes a null sideset B into `Curve( id, A, nullptr )` whose initializer reads `mSideSetB->element_type()` behind a short-circuit that only holds for first-order A → second-order 2-D mesh would dereference null (unreachable today) |
| T-02 `b5a18cf` (thin-shell facet selection) | LOW (high) | a strict narrowing on the single funnel `DofManager::set_equation` → `create_blocks` (`cl_FEM_DofManager.cpp:99`); both equations select by `DomainType` and a shell's layer blocks are only `ThinShell`/`Buffer`, which both lists contain, so every shipped multi-shell deck (gantry 928 shells … corc_garber_effect 3) selects every shell and gets identical ids; no deck key can produce a partial selection (`input_schema.yaml:1237-1247`) | DR-162 unchanged in evidence, option (a) now cheaper (one `key_exists` in the serial side-connector loop `cl_FEM_DofMgr_BlockData.cpp:481`, `:501`); corner: a field selecting only a shell's coating blocks skips the facet link. Agent's low-confidence by-catch `ThermalFactory.cpp:121/127` (share Vector, receive Cell) checked by Claude: `share( Cell )` is documented "same contract as share( Vector )" and both chunk a contiguous `T` array with the same tags → compatible, no defect |
| T-26 `99ac9d0` (constant-n overloads) | LOW (high) | no production caller selects the T-less overloads (all sixteen `MaxwellData` entry points pass `T` → the 4-arg / 8-arg forms, `cl_FEM_Calculator.hpp:3798`, dispatcher `:272-440`); the factory makes jc-function + constant-n mutually exclusive; a deck that had reached them would have SEGV'd (virtual `eval` on null), never miscomputed. Eight value sites + four assertion texts in the diff, five still callerless | O15: delete or forward the T-less family (eight bodies nothing calls)? by-catch: `Material::n( normB, angle[, T] )` has no constants fallback (`powerlaws.hpp:977-996`), same shape one level up, harmless today |

## 5. Audit trail

Exchange threads: `tmp/ai_exchange/review_tgap_<id>.md`, one per row.
