# Cut-Aware Side-Edge Fusing — Single-Authority Branch Resolution (D8)

**Date:** 2026-08-06
**Purpose:** Make the thin-shell rim fuse ("fuse-to-g") cut-consistent: today, different
interface levels of the SAME rim station are constrained to different branches of the
multivalued air trace g (with/without the free-cut λ), producing inter-level trace jumps
of order λ(t) that curl into the rim-strung, time-oscillating J blobs on Gregory's CORC.
Mechanism of the fix in one sentence: resolve each station's source nodes ONCE (a single
authority, cut-composed like the outer-interface anchors already do) and hand the same
sources to every interior level and to the node fusing.
**Module:** `src/fem/kernel` (`cl_ThinShellFactory`), touching the anchor machinery in
`src/fem/maxwell` (`cl_MaxwellFactory`) read-only
**AIs involved:** Claude (diagnosis + plan, this file), executor = fresh session; Codex
audit recommended before implementation (protocol default)
**Status (2026-08-26):** the **debt** on this plan is closed — DR-53 was struck by ruling, and
O1's default question resolved with it: **keep both fuse flags, classified EXPERIMENTAL, hardcoded
`false`, never deck-exposed.** What remains here is a *research* question, not debt: the physics
half of O1 (graded two-branch authority vs the HEX8TB wall element). **Neither flag may be
re-defaulted to `true`** — the Δt-collapse reproducer is `sidecoatings` at t = 2.0825 ms, and
fuse-off passes that cliff warm and cold.

**Earlier status:** R3 IMPLEMENTED AND BUILT 2026-08-06 (committed since — rode in with
the `backup` commits); first post-fix
physics run (`greg/garber2.png`, 23:07): the alternating time-oscillating rim
blobs are GONE — the D8 incoherence is fixed at the constraint level — but the
rim shows coherent edge-current DEPRESSION where Norris-type peaking is expected.
Cause measured (§4a e): 90/348 rim stacks have mismatching λ branches above vs
below the tape (Christian's insight, confirmed same evening), so **O1 is REOPENED**
— single authority is coherent but over-constrains; the through-thickness branch
transition is real physics. R2/R3 artifacts stay (they fixed the incoherence and
are the substrate for any refinement). **NEXT: (i) decision on O1
(free-rims-for-ship vs graded authority vs HEX8TB wall element — Christian +
Prof. Sirous, data in §4a e); (ii) formal R4 probe gate still pending a POST-fix
dump** (`BELFEM_PROBE_FUSED_ROWS=1`, setup-only; the 20:13 dump on disk is
pre-fix; expect residual "mixing" at exactly the 90 λ-mismatch stacks — that is
now the O1 signal, not a bug signature). **`mFuseEdges` switched OFF in source
(`cl_ThinShellFactory.hpp:280`) 2026-08-06 late, per the agreed O1 option (a)** —
the tree had carried `true` since the fusing was made unconditional (b42044da),
which is why Gregory's runs were fused; free rims are the ship state again. Flip
it locally for O1/R4 investigation runs. **Update 2026-08-08 (Christian's
disposition): D8 CLOSED in `hex8tb_phase2_fem_wiring.md`** — the incoherence is
fixed and verified end-to-end; R4/R5 below are demoted to `mFuseEdges` re-enable
preconditions. They gate nothing in the ship configuration and nothing in the
HEX8TB campaign (connector and fuse paths are mutually exclusive branches of
`ThinShellFactory::create()`).

> **2026-08-09 currentness sweep — a second, independent confirmation of the same
> conclusion arrived, from a different direction.** During the side-connector viz-decoupling
> round (`devlog/dl20260809_sideconnector_viz_regression_round.md`), the sidecoatings run
> collapsed its timestep at t = 2.0825 ms against a state-independent ~2e-6 magnetic
> residual floor. Root cause: the **wall-side** fuse
> (`mFuseEdgesWhenHavingSideConnectors`, `cl_ThinShellFactory.hpp:348`) hangs OUTER wall
> edges on tape temp-edge nodes while the new decoupled inner sheets hang 1:1 on tape edges
> — the two constraint families clash at exactly the rim stations where the λ branches
> differ above vs below, i.e. **the same 90/348 stacks this plan's §4a e census identified**.
> Christian's A/B (fuse off passes the cliff warm *and* cold) is reproducer-tier evidence.
>
> Two consequences for this file:
> 1. The "ship = free rims" conclusion now rests on two independent failure modes, not one.
>    **Neither fuse flag should be defaulted back to `true`** without O1 being resolved
>    first. Both are `false` in source today: `mFuseEdges` and
>    `mFuseEdgesWhenHavingSideConnectors` (`cl_ThinShellFactory.hpp:347-348` — note the
>    anchors below say `:276` / `:280`; they have drifted).
> 2. The scope guard claiming the HEX8TB campaign "does NOT depend on the fuse" needs a
>    caveat: it is true of the *connector* path, but the **wall-side** fuse flag is
>    reachable from the side-connector configuration and is what produced the collapse.
>    Read the guard as "the connector element does not need the fuse", not as "the fuse
>    cannot affect a connector run".
>
> O1 (graded two-branch authority vs the HEX8TB wall element as the proper physics) remains
> the open decision, with Christian and Prof. Sirous; ~~the default question (false vs.
> removing the flag) was routed to them on 2026-08-09.~~
> **The default question is RESOLVED 2026-08-26 (Christian's ruling, recorded when DR-53 was
> struck): KEEP the flag, classified EXPERIMENTAL — not for a user audience, and not to be
> exposed.** Both flags stay hardcoded `false`; neither is deck-settable (verified — no matching
> string literal in `src/`, no entry in either input-contract artifact), so reaching the fuse
> requires editing and rebuilding source. That is what makes it an experiment rather than a
> supported configuration, and it is why DR-53 closed as debt. **The physics half of O1 is NOT
> settled by that ruling** and stays open here as a research question, not as debt.
>
> **Physics position (2026-08-13, Christian and Prof. Sirous, recorded in the debt
> register DR-53 and at the flag site):** fusing is believed to be the mathematically
> more correct continuity statement, yet the code converges slower with it and produces
> no significantly better result — the suspicion is that the fuse **overconstrains** the
> problem, analogous to weakly enforcing a B·n = 0 that the formulation already fulfills
> at the boundary. Both flags stay `false`; remove-vs-keep of the flag is still open.

> **Scope guards:**
> - Fix the BRANCH MIXING only. No formulation redesign, no enrichment work, no
>   Alves-style non-duplication rework (rejected 2026-08-06, Christian: collides with
>   the HEX8TB wall element's per-level rim slots, the `.bfm` layout, and the D1/D7
>   investments).
> - The HEX8TB wall element R4/R5 wiring (`todo/hex8tb_phase2_fem_wiring.md`) is a
>   separate campaign and does NOT depend on the fuse — the connector path is a
>   mutually exclusive branch; this plan only makes the fuse itself correct.
> - `mFuseEdges` stays `false` (source-level switch, `cl_ThinShellFactory.hpp:280`)
>   until R5 gates pass. Gregory's production runs are NOT gated by this plan.

---

## 1. Context for a cold session

The long arc lives in `tmp/ai_exchange/side_edge_fusing_handoff.md` (ephemeral — read it
if it still exists) and devlog `dl20260806_matfix_controller_port_executed.md`
(Addendum 2 + the D7/D8 updates). Compressed:

- A stacked h-φ thin shell duplicates rim ("side curve") entities per interface level.
  Published methods (Alves 2022a, `alves2022a.txt:215-217`) have exactly ONE rim trace
  by non-duplication; BELFEM must therefore tie its duplicated rim dofs to the air
  trace g ("fuse-to-g", edge-on-node hanging with LINE2 ±1 coefficients, converted at
  `cl_FEM_DofMgr_DofData.cpp` — search `EntityType::EDGE` in the hanging-dof
  conversion). Free (untied) rims are the validated legacy approximation and the
  current ship state.
- With the fuse ON, Gregory's CORC (deck `cmake-build-debug/greg/input.conf`, periodic
  twisted cell, 3 tapes × 10 layers, cuts terminate on the side curves) shows
  element-scale alternating J blobs strung along the gap-facing rims. The blobs are
  algorithm- and solver-independent (persist under Picard+MUMPS) — they are a property
  of the constraint system, not of the iteration.
- Two prior suspects were ELIMINATED 2026-08-06 by direct evidence: D7 (missing
  `original()` normalization in `connect_side_nodes` — fixed, was real but not the blob
  cause) and the ±1 sign/orientation defect class (refuted: 141,901/141,901 clean).
  The "Newton does nothing" phenomenon was resolved separately and is NOT this bug:
  near-null constant-φ mode of the single-node gauge pin, see
  `src/fem/kernel/doc/bearing_gauge_eigenmode.md`; remedy pure Picard.

## 2. The defect (D8) — verified evidence

Constraint-row dump of the fused Garber setup (probe `BELFEM_PROBE_FUSED_ROWS`, v2,
210,066 rows; analysis 2026-08-06, Claude):

- Row-level integrity CLEAN: every two-source row has +1 on the edge's own first hang
  node; zero internally inconsistent hang-pair groups; zero-weight λ entries are exact
  +λ−λ cancellations of edges whose both endpoints sit on one cut side.
- **61 of 151 rim stations mix cut branches ACROSS levels**: coincident rows of one
  station read `h_e = Δφ` on some levels and `h_e = Δφ ± λ` on others (λ = free-cut
  axial-MMF dof, observed basis id 179126 on tape 1), partly with different node
  representatives outright (e.g. original 3533 vs duplicates 35356 / 157464).
- Consequence: adjacent levels' tangential traces differ by λ(t) ≈ amps → through-
  thickness curl → rim J blobs, oscillating with the drive. 61/151 ≈ 40% of stations
  matches "blobs strung along the whole rim".
- Root cause: `cl_ThinShellFactory.cpp` `connect_side_edges` sources the ORIGINAL
  mid-surface curve nodes (`aEdges( e )->node( s )` — pure Δφ after composition), while
  the outer-interface anchors (`hang_thinshell_edges_on_nodes_bottom/top`, wired from
  `cl_MaxwellFactory.cpp:1469/:1529` — line numbers from the 2026-08-06 jury round,
  re-verify) resolve through the sheet/cut duplicate nodes, whose hanging composition
  correctly picks up the cut jump (Δφ ± λ). One station, two branches of g.

Tracked as **D8** in `todo/hex8tb_phase2_fem_wiring.md` (defect tracker) — keep the two
files in sync when resolving.

## 3. Design: single authority (decided), branch convention (open)

**Single authority (decided 2026-08-06, Christian: "cut-aware method").** Do NOT teach
`connect_side_edges` its own cut logic — a second implementation that must agree with
the first is the failure mode we are fixing. Instead, resolve each station's source
node pair ONCE and reuse it everywhere:

1. For each side-curve station, determine the authoritative source nodes — the node
   representatives (original or cut duplicate) that the chosen authority sheet's rim
   edge actually references at that station.
2. Hand exactly these sources to: all interior-level rim edges
   (`connect_side_edges`, both `Edges` and `EdgeDuplicates`), the interior-level rim
   nodes (`connect_side_nodes` — same authority, same representatives), and verify the
   outer-interface anchors resolve to the same (they define the authority, so this is
   a consistency assert, not a change).
3. Keep the D1/D7 `original()->index()` container-indexing pattern — authority changes
   WHICH node is the source, not how containers are indexed.

**Branch convention at cut-terminating stations (O1, RESOLVED 2026-08-06 — see §9).**
Where a cut terminates on
the side curve, g is genuinely two-valued (branches differ by λ). Fusing all levels to
one branch deliberately flattens that jump at those stations — a localized rim-model
error of the same order as the thin-shell rim neglect. The convention (which sheet is
the authority — lean: bottom sheet, fixed) and the acceptability of the flattening need
Christian + Prof. Sirous sign-off. The physical counterpart of the jump is the
through-thickness bridging of the trace at the cut station; single-branch fusing
suppresses it.

## 4. Ordered steps

- [x] **R1 — Archaeology (read-only).** DONE 2026-08-06 (Claude, this session).
  Findings in §4a below; D9 mechanism = key-collision hypothesis (medium confidence,
  definitive assert folds into R3), D10 classified as multi-generator composition
  (final verdict via the R4 station-consistency gate).
- [x] **R2 — Design the authority API.** A small `ThinShellFactory` helper that, given
  a side-curve station, returns the authoritative source node pair (cut-composed).
  ORDERING CONSTRAINT from R1: the anchors run AFTER the fuse
  (`create_cuts()` at `cl_MaxwellFactory.cpp:444` → `create_thinshells()` :470 →
  anchor PART 2 :1310-1407), so "read what the anchors built" is impossible at fuse
  time. But the anchors contain NO cut logic of their own — their cut-awareness is
  inherited entirely from reading `aFacet->master/slave()->get_nodes_of_facet(...)`
  on volume elements the CutFactory already relinked. The authority helper must
  therefore read the SAME input: the bottom-side volume element's facet nodes at the
  rim station (available inside `ThinShellFactory::create`, since cuts precede it).
  That is still one source of truth — element→node linkage as rewired by CutFactory —
  not a second cut implementation. **Design DRAFTED → §4b (2026-08-06); Codex audit
  before R3 per protocol.** (after: R1; O1 resolved 2026-08-06)
- [x] **R3 — Implement.** DONE 2026-08-06 (Claude, this session; code compiles
  pending Christian's build). `compute_side_authority` added
  (`cl_ThinShellFactory.cpp`, between `compute_side_edge_indices` and
  `connect_side_edges`); `connect_side_edges( Layer* )` takes the flat
  authority-pair container; `connect_side_nodes` takes the per-station authority;
  call site builds both once per curve. `original()->index()` indexing preserved
  (D1/D7). Per design v2: `SideLayer` untouched, anchors untouched. Deviation from
  the draft: the D9 detector is a LOUD LOG (`message( InfoLevel::Default, ... )`),
  not an assert — the dump evidence says it will fire on the greg deck, and an
  abort would block the R4 probe run. Station-match and both-sides-conductor
  checks are `BELFEM_ERROR` (data conditions, always active). (after: R2)
- [ ] **R4 — Constraint regression gate (setup-only, no physics; since 2026-08-08 a
  `mFuseEdges` RE-ENABLE precondition, not a D8 gate).** Fused Garber setup
  run with `BELFEM_PROBE_FUSED_ROWS=1` (flip `mFuseEdges` locally); rerun the
  station-consistency analysis (§6). REVISED criterion (2026-08-06 addendum, made
  formal at the 2026-08-08 D8 closure): **0** stations with differing composed
  constraints OUTSIDE the 90 λ-mismatch stacks of §4a(e) — residual mixing at
  exactly those stacks is the O1 physics signal, not a failure (pre-fix baseline:
  61/151 mixed everywhere). Pairing test stays 141,901/141,901 clean, and the
  R1(c)/R1(d) items are either fixed or explained-and-documented. v2 addition
  (D13): explicitly examine the periodic seam-station rows (z ≈ 0 clusters) for
  stale node-fuse compositions — the census script's multi-λ station list
  localizes them. (after: R3)
- [ ] **R5 — Physics gates (`mFuseEdges` re-enable precondition).** (a) Flat 3-layer homogeneous tape, no cuts/periodicity,
  fused, vs Norris (the long-planned experiment — now as confirmation, not discovery);
  (b) Garber deck fused-vs-free A/B: blobs gone, fields match free-rim away from the
  rim; (c) net-current sum rule per tape (∮h·dl = I(t)) unperturbed by the cut-station
  flattening. Christian runs; Claude analyzes. (after: R4)
- [ ] **R6 — Close-out.** Update `side_coating_wall_element.md` §6 (the orientation
  gate this campaign kept deferring is subsumed by R4/R5 evidence — say so explicitly),
  ~~tick D8 in `hex8tb_phase2_fem_wiring.md`~~ (done 2026-08-08, ahead of R4 — see the
  Status update), devlog, decide whether `mFuseEdges`
  becomes deck-selectable or stays source-level. Probe hygiene: `BELFEM_PROBE_FUSED_ROWS`
  is the regression tool — keep, but note it in the DR-25 probe register.

## 4a. R1 findings (2026-08-06, Claude; line numbers verified this session)

**(a) Anchor resolution path.** `hang_thinshell_edges_on_nodes_bottom`
(`cl_MaxwellFactory.cpp:1469`) sources each sheet rim edge from
`aFacet->master()->get_nodes_of_facet( index_on_master )` — the air volume element's
nodes; `_top` (:1529) same via `slave()` + `to_master_orientation`. **No cut logic
exists in the anchors.** Their cut-awareness is inherited: `create_cuts()` (:444) runs
before everything and the CutFactory relinks air elements to duplicate nodes, so the
volume facet nodes ARE the correct branch representatives. The Δφ ± λ composition
happens later, in the DofData conversion (`cl_FEM_DofMgr_DofData.cpp:3562` EDGE
branch): a source node that `is_hanging()` expands into its own sources/weights
(:3576-3585, :3622-3628), which is where λ enters.

**(b) Fuse resolution path.** `connect_side_edges` (`cl_ThinShellFactory.cpp:2931`
`Layer*`, :2964 `SideLayer*`) sources from `tOrg->node( s )` — the mid-surface curve
edge's own nodes, which the CutFactory never relinked → always the no-λ branch.
`connect_side_nodes` (:3002) sources from the curve node itself (forwarding its
hanging composition if present) but indexes the target by
`tOrg->original()->index()` — so at a duplicate curve station the layer node at the
ORIGINAL's position is fused from whichever representative comes last in the curve
container: last-in-wins branch selection. Fuse call site: :360-388, interior levels
`l = 1 … tNumLayers-2`; outer levels are anchor territory. Ordering:
cuts (:444) → thinshells/fuse (:470) → anchors (:1310) → DofData composition.

**(c) D9 mechanism — 6 of 156 edges missing per interior block.** Census on the
existing dump (`cmake-build-debug/greg/fused_edge_rows_rank0.txt`, 210,066 rows):
9 interior blocks of 156 consecutive edge ids, each missing EXACTLY the same six
offsets **{1, 50, 53, 102, 105, 154}** — three reversal-symmetric pairs
(1+154 = 50+105 = 53+102 = 155), all near the periodic seam / cut-termination
stations (missing-edge neighbors sit at z ≈ 0…1.7·10⁻⁴, λ-carrying rows adjacent).
Leading hypothesis (medium confidence, REVISED during R2 — see §4b): the twin
collapse happens upstream, in `create_temporary_edges`' `unique()` over
original-normalized keys (:1727) — `create_edge_map` (:2821) itself cannot collide
because the temp container it indexes is already unique per original pair. Where the
side curve carries duplicate stations (periodic seam crossings / cut terminations),
two curve positions resolve to ONE temp/layer edge; the survivor is fused, other
layer-edge ids in the block are never referenced → free levels inside an otherwise
fused rim. Definitive test = the repeated-index detector in §4b (fold into R3).

**(d) D10 classified — 425 rows with |weight-sum| ≥ 2.** There are FOUR cohomology
generator dofs in the composed rows: bases 179123, 179124, 179125, 179126 (=the
free axial-MMF λ). The 425 rows are compositions through up to three generators
(e.g. `±(λ₁₂₃+λ₁₂₄+λ₁₂₅)` with clean ±1 φ pairs; no weight anywhere has magnitude 2;
sum = net λ content since φ contributions cancel pairwise). They cluster at ~55
stations, z ∈ {0 … 2·10⁻⁴} — the periodic seam / multi-cut termination zone. This is
structurally LEGITIMATE multi-cut composition; whether the branch combos are
mutually consistent per station is exactly what the R4 station-consistency gate
checks, so D10's verdict rides on R4.

Full weight-sum histogram: 0: 160,609 · ±1: 49,032 · ±2: 351 · ±3: 74. Rim-edge
stations normally carry 11 coincident rows (all levels of the through-thickness
stack). Census script: `tmp/ai_exchange/side_edge_fusing_r1_census.py` (rerunnable
against any new dump; ephemeral — R4 should re-home it if it becomes the gate).

**(e) POST-R3 MEASUREMENT (2026-08-06 evening, Christian's question): the upper and
lower cuts do NOT match along the rim.** Method: per coincident rim stack (≥10
rows, pre-fix dump), take the bottom-sheet anchor row (lowest edge-node ids) and
the top-sheet anchor row (highest) — both are anchor territory, untouched by R3 —
and compare composed constraints. Result: **348 full stacks; 129 identical; 90
differ in λ branch** (almost all ±λ₄ = 179126); 129 differ in φ node
representatives with equal λ. The λ-mismatching stacks form extended stretches
along the rim, not isolated stations. Physical reading: along a stretch where
exactly one cut sheet passes between the below-path and the above-path around the
rim, g_top − g_bot = ±λ(t) is REAL (the MMF of the transport current wrapping the
tape edge); its through-thickness transition profile is solution physics, not a
constraint. Single authority forces all interiors onto the bottom branch, so the
entire λ jump collapses into the topmost interface — consistent with the first
post-fix physics run (`greg/garber2.png`): the old alternating time-oscillating
blobs are GONE (the D8 incoherence is fixed), but the rim now shows a coherent
edge-parallel current DEPRESSION where Norris-type edge peaking is expected.
Refutes the O1 premise that the flattening is localized — see reopened O1, §9.

## 4b. R2 design — the authority API (draft 2026-08-06, Claude; Codex audit pending)

### Additional structural facts (verified this session)

- **The temp-edge layer is branch-blind by construction.** `create_temporary_edges`
  (`cl_ThinShellFactory.cpp:1695`) keys facet edges by original-normalized node pairs
  (:1718-1721), `unique()`s the keys (:1727) — collapsing any cut-side or periodic
  twins into ONE temp edge — and rebuilds each edge's nodes from the positional node
  container (:1737-1738), discarding whatever branch representative the facet element
  actually referenced (:1716 reads it, only the key survives). No fix at the temp-edge
  level can recover the branch; the authority must come from the facets/volume side.
- **Layer::EdgeDuplicates** are material-interface twins: created only when
  `mCreateGhostFacets` is on and adjacent layer materials differ (:175-188). Both
  copies tie to the same g pair; the current identical fusing of both is correct and
  stays.
- **`SideLayer` (HEX8TB side-connector wrap) pre-fuses its outer edges in its own
  constructor** (`cl_ThinShellFactory.hpp:206-213`, sources = reference-edge nodes —
  same no-λ branch defect) and `connect_side_edges(SideLayer*)` then resets and
  re-fuses them. R3 must route BOTH through the authority (or drop the ctor pre-fuse
  as redundant).
- **Layer bottom nodes carry `original()` links to the mid-surface originals**
  (used by the anchor DEBUG check, `cl_MaxwellFactory.cpp:1521-1524`) — this enables
  identity-based station matching with zero orientation logic.

### Design

**Principle (v2 after Codex audit — reconciliation in
`tmp/ai_exchange/side_edge_fusing_r2_audit.md`).** The authority evaluates the *same
expression* the anchors evaluate later —
`tFacet->master()->get_nodes_of_facet( index_on_master )` — and matches stations by
`original()` identity **directly on the volume facet nodes** (cut duplicates carry
`original()` links normalized to the true original, `cl_CutFactory.cpp:2657-2658`;
plain nodes return themselves). No slot conventions, no CW/CCW reasoning, and — v2
correction (audit item 2, D11) — **no layer-node hop**: layer nodes are fresh
`new Node` + `set_index` only (`cl_ThinShellFactory.cpp:1404-1405`), their
`original()` returns `this` (`cl_Node.hpp:319`), so the originally proposed
`get_bottom_nodes` matching can never fire. Consistency with the anchors is by
construction; R4 stays the empirical gate.

**Side selection (v2, audit item 3, D12).** Authority side = `master()` unless the
master block is `DomainType::Conductor`, in which case fall back to `slave()` —
identity matching is orientation-free, so the fallback needs no
`to_master_orientation`. `BELFEM_ERROR` if both sides are conductor (no φ authority
exists there). This mirrors the PART 2 routing (`cl_MaxwellFactory.cpp:1355-1376`).

**New helper** (private, `cl_ThinShellFactory`):

```cpp
void compute_side_authority(
    Curve                 * aCurve,          // stations, in traversal order
    Cell< Facet * >       & aFacets,         // mid-surface facets (create() scope)
    Cell< Edge * >        & aEdges,          // temp edges
    const Cell< index_t > & aEdgeIndices,    // curve position -> temp-edge index
    Cell< Cell< Node * > >& aEdgeSources,    // out: per curve position, source PAIR,
                                             //      ordered like the temp edge's own nodes
    Cell< Node * >        & aNodeSources );  // out: per curve station, ONE authority node
```

Per curve position `p` with temp edge `tOrg = aEdges( aEdgeIndices( p ) )`:
1. Find the adjacent mid-surface facet + slot via the existing
   `create_edge_to_face_map` (:2834) — already handles flag hygiene and asserts
   one side edge per facet.
2. Pick the authority side (master unless conductor, see above); read
   `get_nodes_of_facet( ... , tVolumeNodes )` on that volume element; for each
   temp-edge node `s`, find the `k` with
   `tVolumeNodes( k )->original()->index() == tOrg->node( s )->original()->index()`
   → `aEdgeSources( p )( s ) = tVolumeNodes( k )`. `BELFEM_ASSERT` both stations
   match exactly once (an unlinked "decoupled duplicate",
   `cl_InterfaceProcessor.cpp:673`, would fail loudly here — that is wanted).
3. `aNodeSources`: station authority = the volume node resolved from the facet of the
   FIRST curve position adjacent to that station (traversal order — deterministic).
   This is the ONLY place the O1 flattening enters: edges are side-resolved exactly
   like the anchors (adjacent facets across a cut termination legitimately use
   different representatives — the λ jump lives BETWEEN midpoints, never within one),
   while the single per-station node dof must pick one branch.

**Rewiring (R3, v2).**
- `connect_side_edges( Layer * ... )` (:2931): add the `aEdgeSources` parameter;
  `tDup->add_source( aEdgeSources( p )( s ) )` replaces
  `tDup->add_source( tOrg->node( s ) )`. Source order stays positionally aligned
  with the edge's own nodes — the DofData ±1 convention
  (`cl_FEM_DofMgr_DofData.cpp:3635-3641`) is untouched; the R4 pairing test
  (141,901/141,901) re-verifies this.
- `connect_side_nodes` (:3002): source = `aNodeSources( station )`; keep the
  existing hanging-forwarding (if the authority node `is_hanging()`, copy its
  sources/weights — same pattern as today; pointer semantics are NOT an option,
  DofData's NODE branch asserts non-hanging sources,
  `cl_FEM_DofMgr_DofData.cpp:3506`). Target indexing UNCHANGED
  (`original()->index()`, D1/D7).
- **SideLayer: OUT OF SCOPE (v2, audit item 7, D14).** The
  `connect_side_edges( SideLayer * )` overload (:2964) is never called; the
  SideLayer construction site lives in the side-connector path that hard-aborts
  (`BELFEM_ERROR` at :358). Leave the ctor pre-fuse (hpp :206-213) and the unused
  overload untouched; wiring the authority into the HEX8TB wrap belongs to
  `todo/hex8tb_phase2_fem_wiring.md` (note added there at R6).
- Call site (:360-388): build the edge-to-face map + authority once per curve,
  before the `l`-loop.

**D13 — known, pre-existing, unchanged (audit item 4).** Periodicity `update()` +
`set_entity_dependencies()` run AFTER `create_thinshells()`
(`cl_MaxwellFactory.cpp:483-488`) and reset periodic-slave node source containers
(`cl_Mesh_Periodicity.cpp:60`), so a node-fuse snapshot taken at fuse time from a
periodic-slave authority node can be superseded. Today's `connect_side_nodes` has
the IDENTICAL exposure (it snapshots the curve node); `set_entity_dependencies`
then reconciles periodic layer pairs by copying the master's expansion (:62-67).
The authority change neither widens nor fixes this. Disposition: R4 gets an
explicit seam-station node-row check; any fix is a separate campaign.

**D9 detector (R3, costs one loop).** `BELFEM_ASSERT` (or loud `gLog` + count) that
`aEdgeIndices` contains no repeated temp-edge index per curve. Revised D9 hypothesis
after the `create_temporary_edges` finding: the twin collapse happens in the
`unique()` at :1727 (periodic seam crossings and/or cut-side twins whose
original-pairs coincide), so two curve positions can map to ONE temp/layer edge —
the survivor is fused (with the authority fix: per-position sources now make the
LAST-writer position win explicitly), and some other layer-edge ids in the block are
never referenced → the 6 missing rows. If the detector fires, the residual question
is whether the collapsed positions are periodic-identified anyway (then one row is
legitimately shared, document and keep) or genuinely distinct (then the layer needs
per-branch edge entities — escalate to Christian before widening scope).

### What this deliberately does NOT do

- No new cut logic, no reading of cut/cohomology data structures — the branch
  information is consumed exclusively through element→node linkage that CutFactory
  already rewired (cuts precede thinshells, `cl_MaxwellFactory.cpp:444` vs :470).
- No change to which entities exist (`.bfm` layout, HEX8TB per-level rim slots,
  D1/D7 container indexing all untouched).
- No touching of the anchor machinery (it defines the authority; it is read-only
  reference behavior).

## 5. Defects

- [ ] **D8 (HIGH, driving defect):** cross-level cut-branch mixing at rim stations —
  full evidence in §2 and in `todo/hex8tb_phase2_fem_wiring.md` D8. Fix = R2/R3.
- [ ] **D9 (open, mechanism hypothesized — see §4a c):** 6 of 156 interior-block rim
  edges absent per level, same six offsets {1, 50, 53, 102, 105, 154} in all 9
  blocks, at seam/cut-termination stations. Hypothesis (revised in R2): twin collapse
  in `create_temporary_edges`' `unique()` (:1727) maps two curve positions onto one
  temp/layer edge → some layer copies are never fused (FREE level inside a fused
  rim). Confirm with the §4b repeated-index detector in R3; disposition depends on
  whether collapsed positions are periodic-identified (document) or distinct
  (escalate — needs per-branch edge entities, scope decision by Christian).
- [ ] **D10 (classified benign-pending-R4 — see §4a d):** the 425 multi-λ rows are
  compositions through up to three of the four generators (179123-179126) at ~55
  seam-zone stations; no anomalous weights. Verdict = R4 station consistency; no
  separate fix expected.
- [x] **D11 (design defect, CLOSED in design v2):** the drafted station matching went
  through layer-node `original()` links that do not exist (audit item 2, verified) —
  v2 matches on the volume facet nodes directly. Never reached source.
- [x] **D12 (design gap, CLOSED in design v2):** unguarded `master()` read at
  conductor-master facets — v2 adds the DomainType side selection with slave
  fallback (audit item 3).
- [ ] **D13 (pre-existing, documented — audit item 4, verified):** node-fuse
  snapshots can be superseded by the post-fuse periodicity reset
  (`cl_MaxwellFactory.cpp:483-488` → `cl_Mesh_Periodicity.cpp:60`); identical
  exposure exists in today's code; authority change is neutral. Checked empirically
  in R4 (seam-station rows); any fix is out of scope here.
- [x] **D14 (design overreach, CLOSED in design v2):** the proposed SideLayer ctor
  pre-fuse deletion would have changed behavior — the `SideLayer*` connect overload
  is never called and the whole path hard-aborts (audit item 7, verified). SideLayer
  untouched; deferred to `hex8tb_phase2_fem_wiring.md`.

## 6. The regression analysis (runnable, cold)

Probe: env `BELFEM_PROBE_FUSED_ROWS=1`, fused build, any run of the deck — the file
`fused_edge_rows_rank<r>.txt` is complete at setup (kill after the first residual
line). v2 line format:

```
edge <id> nodes <n0> <n1> orig <o0> <o1> pdc <p0><d0> <p1><d1> hang <s0> <s1> src <basis>:<type>:<w> ... x <x> <y> <z>
```

`hang` = nodes the edge hangs on, in the edge's OWN node order; `src` = the final
composed constraint row (solver order). Checks (python, as run 2026-08-06):

1. *Pairing:* rows where both hang nodes appear directly in `src` must have
   +1 on `hang[0]`, −1 on `hang[1]`. Expected: 100% clean.
2. *Station consistency:* strip zero-weight `src` entries; cluster rows by midpoint
   (grid hash, ±2·10⁻⁶ m, merge neighbor cells — layer offsets are ~10⁻⁷ m); within
   every cluster of ≥8 rows, the set of composed constraints must have size 1.
   Expected after fix: 0 violating stations (baseline: 61 of 151).
3. *Census:* weight-sum histogram (±1 = single-λ rows, fine; ±2/±3 = D10), rows per
   edge-id block (156 expected per interior level; 150 seen = D9).

## 7. Definition of done

R1–R4 landed and green (revised criterion: 0 mixed stations outside the 90 λ-mismatch
stacks), R5 physics gates pass on Christian's runs, O1 signed
off and documented (branch convention + flattening note in
`src/fem/kernel/doc/` — extend `hanging_dofs_static_condensation.md` or a new
`thin_shell_rim_fusing.md`), D9/D10 dispositioned (D8 closed 2026-08-08 in
`hex8tb_phase2_fem_wiring.md`; the wall-element campaign is NOT gated on the fuse —
connector and fuse paths are mutually exclusive branches of `ThinShellFactory::create()`).
Plan moved to `todo/closed/` with DONE summary.

## 8. Audit trail

- `tmp/ai_exchange/side_edge_fusing_r2_audit.md` — Codex audit of the §4b design
  (brief, 8 findings, Claude verification pass + reconciliation table; log:
  `_codex_r2_authority.log`). Ephemeral; dispositions distilled into §4b/§5.
- `tmp/ai_exchange/side_edge_fusing_handoff.md` — the full campaign state (two jury
  rounds + the 2026-08-06 D7/D8 updates); ephemeral, distilled here.
- `tmp/ai_exchange/review_side_edge_fusing_physics.md`,
  `review_greg_corc_regression.md` — the verified jury rounds (ephemeral).
- Devlogs: `dl20260806_matfix_controller_port_executed.md` (Addendum 2 + updates),
  `dl20260806_greg_corc_convergence_regression.md`.
- `src/fem/kernel/doc/bearing_gauge_eigenmode.md` — the separately-resolved Newton
  paralysis (NOT this bug; do not re-conflate them).

## 9. Open questions

- [ ] **O1 — REOPENED 2026-08-06 evening (was: closed same day).** The closure
  rested on the premise that the g-jump flattening is LOCALIZED to cut-terminating
  stations. Christian's post-run question ("the cuts on the upper and lower side
  don't necessarily match") is confirmed by measurement (§4a e): 90 of 348 rim
  stacks carry different λ branches above vs below, in extended stretches — the
  two-valuedness is a property of whole rim segments, and the through-thickness
  λ transition is real transport-current physics that single-branch fusing
  suppresses (garber2.png: edge current down instead of up). Options on the table:
  (a) ship with free rims (`mFuseEdges = false` — already the default; validated
  legacy), keep the R3 single-authority code as the now-coherent fuse
  infrastructure; (b) graded authority — interior level k hangs on BOTH branches
  with complementary weights (imposes a linear through-thickness transition;
  needs the DofData edge-on-node conversion to accept 4 sources — real surgery);
  (c) proper rim physics via the HEX8TB wall element (R5 r′ theory — the
  standing campaign). Needs Christian + Prof. Sirous with the §4a(e) data.
- [ ] **O2** — Should the R4 station-consistency check become a checked-in regression
  test (battery-style, small fused fixture) rather than a manual probe run? (Lean:
  yes, after R5 — fits `make check-fast`.)
