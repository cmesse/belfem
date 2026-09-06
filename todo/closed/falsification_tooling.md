# Falsification Tooling — Regression Battery, Evidence Hierarchy, Campaign State

**Date:** 2026-08-05
**Purpose:** Institutionalize the falsification side of the workflow: a permanent
interface/orientation regression battery (`make check-fast`, both backends), an evidence
hierarchy wired into the collaboration protocol and `/cross-review`, and a compression
layer over the devlog (campaign state pages + a 1.0-lens debt register). Complements the
committed cross-review tooling (`todo/closed/cross_review_tooling.md`).
**Module:** `tests/fem`, `doc/`, `devlog/campaigns/`, `todo/`
**AIs involved:** Claude (plan + implementation + builds/tests, newly sanctioned), Codex +
Grok (audit of the battery's expected-value derivations and the protocol diffs)
**Status:** EXECUTABLE CORE COMPLETE 2026-08-05 — battery green on both backends
(`make check-fast`: Armadillo 0.93 s / Blaze 0.50 s / `env -i` pass), mutation check
red(4)/green(10) demonstrated, two defects found by the battery itself (DR-47 HEX8TB
factory case, DR-48 test_math broken). AWAITING CHRISTIAN: (1) EXPECTED sign-off block
(devlog dl20260805_falsification_tooling), (2) D2 protocol/command diffs approval,
(3) correction pass on campaign pages + debt register, (4) commit slicing. Deferred
follow-ups: ghost Calculator fixture (DR-46), dedicated 18-case pairing sweep (R5◐),
Codex/Grok wording audit of the D2 diffs.

**2026-08-09 currentness sweep — one of the four AWAITING-CHRISTIAN items is already
satisfied, and both day-one defects are closed:**
- **(2) D2 protocol/command diffs — the edits are IN, not pending.**
  `doc/ai_collaboration_protocol.md` now carries **§11 Evidence Hierarchy** with the
  evidence ladder (`:354-397`), the "reviewed ≠ verified" Verification field in the devlog
  template (`:220`), and the registry rows for `devlog/campaigns/<name>.md` and
  `todo/debt_register.md` (`:347-348`). What is still owed is Christian's *review* of that
  wording, not its application.
- ~~**DR-47 is FIXED:** ... the disabled `Hex8TbUnitCirculation` test is unblocked and should
  be re-enabled — that is a concrete, small next action for this file.~~ **Done 2026-08-10:**
  the factory case landed 2026-08-08 (`cl_IF_InterpolationFunctionFactory.cpp:198`) and the
  `DISABLED_` prefix is now removed. The test's own `EXPECTED: pending Christian sign-off`
  line stays — it belongs to item (1) below, not to DR-47.
- **DR-48 was closed 2026-08-05** (already reflected in the register).
- Items (1) EXPECTED sign-off, (3) correction pass on the campaign pages + register, and
  (4) commit slicing are still open. Note that (3) has grown: the six campaign pages under
  `devlog/campaigns/` were seeded 2026-08-05 and at least
  `controller_anderson.md` is now factually stale (it still says the Anderson bundle is
  UNCOMMITTED; it was committed in `6b2a2b98`).
- The one open confirmation inherited from `closed/cross_review_tooling.md` — the
  `scls_env` g++ guard being existence-only rather than under `/opt/scls/` — is parked here
  so it does not get lost.

> **Scope guards (from the task brief):**
> - Physics-load-bearing expected values (telescope sums, sign conventions) are PROPOSED
>   with written derivations and flagged `// EXPECTED: pending Christian sign-off` —
>   reviewer agreement does not settle physics.
> - Battery stays fast: minutes, reference elements + code-built micro-meshes, no MPI
>   solve, no gmsh input. Full-model smoke runs remain Christian's.
> - No CI system, no new frameworks, no YAML. Campaign pages are markdown; the register
>   is one markdown table.
> - Everything runs through `scripts/scls_env.sh`; `make check-fast` must pass from
>   `env -i`.
> - D3 seeding is timeboxed to ONE focused pass over the 117 devlogs since 2026-06-01
>   (Open-items sections only); if that is not enough, stop and ask.
> - Mutation acceptance: the battery must go red on the reintroduced `-mS[1]` flip —
>   a battery that cannot catch its target bug class is decoration.

---

## §1 — Recon facts the design rests on

- Test infra: gtest per module via `config/scripts/Add_Test.cmake` (one executable
  `test_<module>` + one ctest entry per module); `make check` runs everything
  (`CMakeLists.txt:278-290`). `tests/fem/` exists (4 suites, kernel+interpolation+mesh
  in LIBLIST). Reference elements come from `mesh::ElementFactory::
  create_reference_element` (pattern: `tests/fem/test_FacetIntegrationPoints.cpp:45-48`).
- Edge functions live in `src/fem/interpolation/nedelec/`: `EF_QUAD4TS`, `EF_PENTA6TS`,
  `EF_HEX8TS`, `EF_HEX8TB` share the `EdgeFunction` API `precompute(aXi)` /
  `link(Element*)` / `E(i)` / `C(i)`.
- **Feasibility crux:** `EF_QUAD4TS::link()` needs a full `fem::Element` — facet nodes,
  `element()` nodes, `edge_directions()`, and the thickness via
  `aElement->parent()->parent()->mesh()->block(id)->thickness()`
  (`cl_EF_QUAD4TS.cpp:92`, NaN falls back to geometric thickness but the parent chain
  is still dereferenced). Standalone EF tests need either a minimal real
  Group/DofManager or a seam. → O1.
- Mutation target (greg2 alternation bug, fixed in `4a42d982`): top-edge dof sign at
  `cl_EF_QUAD4TS.cpp:166-167` — `mE(·,1) = mS[1]·…`; the historical defect was `-mS[1]`.
  The telescope identity it broke: "per-layer currents telescope to the outer-trace
  difference (Ampère)" (dl20260804_2d_thinshell_layer_alternation.md:107).
- Orientation family: "18 orientation cases (12 quad + 6 tri)"
  (dl20260320_vertex_capacity_refactoring.md:87); PENTA-TS facet-normal sign pattern
  documented in dl20260422_thinshell_facet_renumbering_break.md:121-196
  (`normal_penta_ts` cases vs `normal_penta`, `cl_FEM_Calculator.cpp:1816-1912`).
- Ghost facets: inter-layer facets exist only at material-interface layers; the
  `h_ghost` trace and the deferred dof-count check are in
  dl20260804_2d_thinshell_layer_alternation.md:117-135. The exact "12-dof" contract
  source gets pinned in R4 (the devlogs defer it by description, not by class name).
- Backends: `cmake-build-debug` is Blaze (`USE_MATRIX_BLAZE=ON`). No Armadillo tree
  exists → R6 creates one. Shared-tree rule applies: Christian may be running `make`
  in `cmake-build-debug`; I build the new Armadillo tree freely but touch the debug
  tree only for targeted test targets, and stop on any collision symptom.
- HEX8TB status: `cl_EF_HEX8TB.cpp` exists (268 lines, exact-cuboid design) but phase-2
  FEM wiring is blocked on the edge-function derivation
  (`todo/hex8tb_phase2_fem_wiring.md` R1). → O4.

---

## §2 — Design

### D1 — `tests/fem/test_InterfaceOrientation.cpp` + `make check-fast`

**Fixture (R1):** a small test-support header `tests/fem/support/fn_make_thinshell_stack.hpp`
that code-builds micro-meshes (no gmsh, no file I/O): an N-layer 2-D QUAD4TS stack with
LINE2 facets and controlled node ordering/edge directions, plus single reference
PENTA6TS / HEX8TS / HEX8TB elements. Whether the fixture can hand the EF a usable
`fem::Element` decides O1.

**Test groups (R2-R5):**

1. **Circulation identities (R2).** For each EF and each edge dof k: Gauss-integrate
   `∫_edge_j E_k · dl` on the reference element → `s_k δ_jk` (unit circulation on its
   own edge, zero on others); E-column face activity (QUAD4TS bottom columns vanish at
   η=+1, top at η=−1; HEX8TS/PENTA6TS analogues at ζ=±1); sign consistency against
   `edge_directions()` for both orientations of each edge.
2. **Ampère telescope (R3).** N=4 QUAD4TS stack (greg2 configuration): prescribe edge
   dofs q for (a) a uniform tangential field — per-layer curl must be ZERO per layer
   (the alternation bug made it oscillate), and (b) a linear through-thickness profile —
   per-layer currents `I_l = ∫ C·q dA` must telescope: `Σ_l I_l = H_t(top trace) −
   H_t(bottom trace)`. This is the minutes-scale test that would have caught greg2.
3. **Ghost contract (R4).** Pin the ghost facet element's dof-count (the deferred
   12-dof check) and its element-matrix sparsity pattern from the `h_ghost` path;
   material-interface stack variant of the fixture (ghosts exist only there).
4. **Orientation sweep (R5).** All 18 facet-pairing cases (12 quad + 6 tri) through
   `to_master_orientation` + the thin-shell facet pairing; PENTA-TS facet-index/normal
   sign table against the `normal_penta` reference pattern.

**Expected values:** every physics-load-bearing constant carries a derivation comment +
`// EXPECTED: pending Christian sign-off`; all of them are ALSO collected in one
sign-off block in the devlog (acceptance 3). Structural expectations (δ_jk shape,
sparsity, dof counts, orientation-table involution) are not physics and not flagged.

**`make check-fast` (R6):** `Add_Test.cmake` gains an optional `TESTLABELS` variable →
`set_tests_properties(... PROPERTIES LABELS ...)`; the fast set (proposed: containers,
linalg, math, mesh, fem) gets the `fast` label; new root target `check-fast` = build
those test targets + `ctest -L fast --output-on-failure`. Verified from `env -i` through
`scls_env.sh`, on Blaze (`cmake-build-debug`) and a new `cmake-build-armadillo` tree.

**Mutation check (R7):** re-introduce `-mS[1]` at `cl_EF_QUAD4TS.cpp:166-167` locally →
`make check-fast` must go red (telescope + circulation sign tests); revert → green.
Recorded with the exact commands in the devlog. Never committed.

### D2 — Evidence hierarchy (protocol/command/template diffs, shown as diffs)

1. `doc/ai_collaboration_protocol.md`: new section **"Evidence Hierarchy"** —
   (a) vocabulary rule: *reviewed* (static audit done) ≠ *verified* (executable gate
   passed); devlogs must not say "verified" for syntax-only work;
   (b) the ladder, strongest first: end-to-end reproducer > focused regression >
   compile/link > numeric probe > static source trace > literature consistency > AI
   agreement — lower levels support, never replace, higher ones;
   (c) **review stop condition**: a new audit round requires material code change, new
   experimental evidence, reviewer disagreement on a load-bearing point, or a safety
   boundary (ownership, parallelism, persistence, ABI, formulation) — otherwise the
   next step is the executable gate;
   (d) same-session rule for D3: a session resolving/creating an open item updates the
   campaign page + register in that session; the devlog links, it does not restate.
2. `.claude/commands/cross-review.md`: reconciliation table gains an **evidence** column
   (highest ladder level backing the verdict); step 6 gains the stop-condition note.
3. Devlog template (protocol §6): one-line `**Verification:**` header field — highest
   evidence level reached; for probe level or above also branch/commit + the command
   (one line, not a form).

### D3 — Campaign pages + debt register

1. `devlog/campaigns/<name>.md`, ACTIVE campaigns only. Proposed list (→ O3):
   `side_coating_wall_element`, `2d_thinshell_validation`, `controller_anderson`,
   `bfm_persistence`, `gasmodels_migration`, `release_1.0`. Each ≤1 page: current
   accepted design, branch, last passing reproducer, open P0/P1, superseded approaches,
   links to dated entries. Seeded from June+ devlogs; every seeded claim tagged
   `[seeded — confirm]` for Christian's correction pass before commit.
2. `todo/debt_register.md`: ONE table `ID | area | severity | status | reproducer |
   blocking-1.0?`, seeded from Open-items sections of the 117 June-onward devlogs in a
   single pass (grep-assisted: `## Open`-family headings), blocking-1.0 column filled
   as proposal for the September lens. Register IDs are stable (`DR-nn`).
3. The same-session update rule ships in the D2 protocol diff (2d above).

---

## §3 — Steps

- [x] R1 — 2026-08-05: `tests/fem/support/cl_TS_TestStack.hpp` — TS_TestStack2D (N-layer
      QUAD4TS, rotation + edge-flip options), TS_TestPrism (PENTA6TS/HEX8TS), TS_TestWall
      (HEX8TB + recovery facet). O1(b) WORKED: Mesh → Kernel → DofManagerBase → Block →
      aura `fem::Element` + `set_facet`, zero production seams. Gotcha: edges must be
      created AFTER the Kernel (MeshChecker rejects pre-existing edges) and are
      fixture-owned
- [x] R2 — 2026-08-05: QUAD4TS circulation/sign/activity + Stokes E–C consistency;
      PENTA6TS 6×6 and HEX8TS 8×8 circulation matrices = identity (canonical edges),
      flip isolation. ~~HEX8TB~~ DISABLED: no Lagrange factory case for HEX8TB
      (cl_IF_InterpolationFunctionFactory.cpp:215) — the hex8tb_phase2 R4 gap, now with
      executable evidence (DR-47)
- [x] R3 — 2026-08-05: N=4 Ampère telescope (uniform ⇒ I_l=0; linear ⇒ I_l = h_l−h_{l+1},
      Σ = h_0−h_N) + inter-layer continuity (the greg2 mechanism), rotated variants
- [◐] R4 — Ghost contract pinned (h_ghost 4-block Nitsche, 2×6 dofs,
      mt_maxwell_h.cpp:1805-1937) but needs a Calculator fixture — scaffolded DISABLED,
      DR-46
- [◐] R5 — Orientation: `get_top_nodes`/`get_bottom_nodes` positional-alignment test
      (the 43474c9f class) for QUAD4TS/PENTA6TS/HEX8TS; the full 18-case facet-pairing
      sweep partially covered by the existing test_FacetIntegrationPoints orientation
      sweep — a dedicated thin-shell pairing sweep remains follow-up work
- [x] R6 — 2026-08-05: TESTLABELS + `check-fast`. ~~math OUT → DR-48~~ DR-48 CLOSED same
      day (Christian's ask): Quaternion value-type reconcile, Tensor matching-shape
      rule, Blaze-safe reference loops, + the fn_ddot/fn_kelvin_christoffel
      Armadillo-only-signature source fix. Fast set = containers/linalg/math/mesh/fem,
      green: Armadillo 0.70 s, Blaze 0.65 s, `env -i` pass. New trees
      `cmake-build-arma-test` / `cmake-build-blaze-test`; `cmake-build-debug` untouched.
      New DR-49: spline tests SuiteSparse-gated, compiled out in all current configs
- [x] R7 — 2026-08-05: mutation `-mS[1]` → 4 tests RED (circulation, sign, Stokes,
      continuity; C-only telescope green as designed — Stokes carries the E↔C tie),
      revert → 10 GREEN. Commands in the devlog
- [◐] R8 — D2 edits in place (protocol §11, template Verification field, cross-review
      evidence column + stop condition); diffs presented for approval — Codex/Grok
      wording audit deferred until after Christian's approval
- [◐] R9 — Six campaign pages seeded, all claims `[seeded — confirm]`; awaiting
      Christian's correction pass (NOT committed)
- [◐] R10 — Register seeded (DR-01…DR-48) in ONE pass over the 56 June+ Open-items
      sections (timebox held); blocking-1.0 proposals await the correction pass
- [x] R11 — 2026-08-05: devlog `dl20260805_falsification_tooling.md` (mutation record,
      consolidated EXPECTED block, Verification field exercised); READMEs updated
- [x] R12 — 2026-08-14: **the volume-element half of the battery landed** (Christian's ask,
      the day the TET defects surfaced): `tests/fem/test_EdgeFunctions.cpp` + fixture
      `support/cl_EF_TestVolume.hpp` (TRI3/TRI6/TET4/TET10 over the TS_TestStack kernel-chain
      pattern, canonical edges + self-mastered faces). Circulation identity on reference AND
      distorted elements, C-vs-central-difference curl (exact for these degrees), EXODUS
      detJ > 0 guard (pins the node-map booby trap), edge-flip isolation, TET10
      curved/straight path equivalence. Probe-run 22/22 green (46/46 with the Lagrange suite)
      against the fixed `EF_TET10`; against the pre-fix library it reproduces the defect
      signature — the battery catches D1/D2/D3 of `todo/nedelec_edge_function_defects.md`,
      and its curved-path test is what FOUND D3. Interior point deliberately at column 0
      (link's Jacobian point), because a boundary point there masks the D3 class.

## §4 — Open questions

- O1 — RESOLVED 2026-08-05 (Christian): probe (b) minimal-Group `fem::Element` first;
  fall back to (a) code-built micro-mesh over a minimal real Kernel — build (a) even if
  heavy (the fixture amortizes over every future interface test). **(c) is OFF the
  table for link():** thickness recovery is a designed data path (the HEX8TB
  recovery-facet route); a test seam there could mask exactly the defect class the
  battery targets.
- O2 — RESOLVED 2026-08-05 (Christian): module list approved as SEED only; the criterion
  is per-test wall time — label individual tests FAST (<~10 s each, whole set <~3 min).
  A slow test in a fast module stays out.
- O3 — RESOLVED 2026-08-05 (Christian): add `gasmodels_migration`. Periodic and
  double-corc are dormant as campaigns — fold their remaining open items into the
  `release_1.0` page as validation/debt rows, no standalone pages.
- O4 — RESOLVED 2026-08-05 (Christian): approved — battery covers the current
  exact-cuboid EF_HEX8TB only; wall-term tests stay with hex8tb_phase2_fem_wiring.md
  behind the WIP stop.
- O5 — RESOLVED: `cmake-build-armadillo/` beside `cmake-build-debug/` (mine); debug tree
  only for targeted test targets; the scls_env generator pin applies to both.
- O6 — RESOLVED 2026-08-05 (Christian): approved, with the addition that a signed-off
  value's pending flag is REPLACED by a one-line provenance comment (derivation source:
  paper+section or hand-check reference) so the test file carries its own adjudication
  trail. Post-sign-off values count as verified at focused-regression level.

## §5 — Risks

- fem::Element constructibility may force option O1(c) (source seam) — flagged, not
  assumed approved.
- `check-fast` wall time: mesh module tests may exceed "minutes" on cold build; the
  label set is adjustable (O2) without touching the battery.
- The 117-devlog seeding pass is large; the timebox guard applies (stop and ask rather
  than a shallow-but-wrong register).
- Concurrent builds in `cmake-build-debug` (shared-tree rule): targeted targets only,
  stop on missing/partial `.o` symptoms.
