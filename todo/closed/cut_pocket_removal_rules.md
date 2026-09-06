# Cut Pocket Removal: Rules of Engagement (pre-implementation agreement)

**Date:** 2026-07-16
**Purpose:** Agree on the rules for a pocket-removal / representative-canonicalization post-pass on rectified cohomology generators, before any code is written. The pass removes closed, support-only cut fragments by firing exact quotient-node coboundaries through the already-audited `cohomology::fire_node_coboundary()`.
**Module:** `src/homology` (+ flood-fill primitive from `src/math/graph` / revived `./archive/graph`)
**AIs involved:** Claude (draft), Codex (audit + prose), Christian (rule decisions)
**Status:** ✅ IMPLEMENTED 2026-07-16 (`Cohomology::remove_cut_pockets`, census + Tier-A
firing live in `clean_spfa()`; conservative O1 predicate = any-sideset node contact;
O2 → Tier B via that predicate; O3 → pure-only as collapsed; O4 → local BFS, revived
Tarjan kept as future cross-check). Field outcome on double-layer corc: the census
became the decisive diagnostic — it exposed the SPFA **density** defect (support 226k
vs 5.4k) that motivated the certify-then-greedy pivot, and after the sparse rerun it
proved the `match_edges` abort is the Rule-7 seam-emission defect (mismatch scales
and flips sign with cap-touching cut mass), NOT a removable pocket — zero Tier-A
candidates existed on corc, exactly as the rules predicted for a cap-footprint
symptom. Tier-B enablement remains blocked on validation cases. Devlog:
`devlog/dl20260716_census_density_pivot.md`.
>
> **Closed 2026-08-09** (todo currentness sweep). All seven rules and O1–O4 are
> ticked: `Cohomology::remove_cut_pockets( bool aFireTierA )` is in the tree and
> committed (`cl_Cohomology.hpp:155`, called from `clean_spfa()`), together with
> `fire_node_coboundary()` and `rectify_greedy_sweeps()`. The one item this file
> deliberately leaves open — **Tier-B enablement** — is not work-in-flight but a
> gate waiting on validation cases per structure type; it is carried as a future
> option here rather than as an active task. Rule 7's excluded defect (essential
> cut terminating on a periodic cap) was diagnosed and fixed separately under
> `periodic_cap_cut_emission.md`, exactly as this file predicted.

> **Scope guards:**
> - This pass is a canonicalization of an already-valid unit generator. It is not the fix for an essential cut that terminates on a periodic cap (seam-face continuation); that is a separate defect with its own plan (see Rule 7).
> - The Phase-4 min-cost optimal representative (dual-area weights) remains out of scope. This pass is a cheap monotone local cleanup.
> - Runs after `clean_spfa()` rectification only; the greedy path gets nothing.

---

## §1 Definitions

All statements below live on the quotient graph used by `clean_spfa()` and `fire_node_coboundary()`: slave periodic nodes fold onto their masters, self-paired nodes remain their own representatives, and the edge domain is `original_edges()`.

- **Support:** number of nonzero coefficients of the generator on the complex edge domain.
- **Zero-graph:** quotient representative nodes connected by current domain edges whose rectified coefficient is `c' = 0`. Components are always computed from the current generator state.
- **Component flip:** for a current zero-graph component `K` and sign `s ∈ {+1, -1}`, fire every representative node of `K` with weight `s`, adding the exact coboundary `s·δ(1_K)`. Edges with both endpoints in `K` are unchanged (this includes any nonzero edges interior to `K` — the in-`K` support is not claimed to be exactly `δ(1_K)`); only edges between `K` and its complement can change.
- **Pure pocket:** a current zero-graph component for which one flip sign zeroes every boundary edge `∂K`. In the rectified unit setting, this is the only admissible support-reducing component flip.

## §2 The Rules

- [x] **Rule 1 — Class safety is unconditional for quotient component flips.**
  Every accepted flip adds an exact coboundary, so the cohomology class is preserved (`dd = 0`). This is the same quotient firing convention already used by `clean_spfa()` and `fire_node_coboundary()` (`cl_Cohomology.cpp:346-383`, flags set at `:405-423`). The risk budget is spent on admissibility and downstream structure effects, not on the cohomology class.

- [x] **Rule 2 — Unit preservation is a hard admissibility guard.**
  Before firing, check every current boundary edge `e ∈ ∂K`. The proposed increment on `e` must not have the same sign as the current coefficient; otherwise `±1` would become `±2`. The release-active unit gate reruns after the pass. (Consequence, Codex-verified: for a *current* zero-component every `∂K` edge is nonzero, so an admissible flip zeroes all of them — **admissible ≡ pure**; a mixed-orientation boundary has no admissible sign.)

- [x] **Rule 3 — Quotient awareness is mandatory.**
  Pockets can cross the periodic plane. The flood-fill therefore runs on the quotient graph, and a crossing pocket is one component whose flip fires both geometric sides through `fire_node_coboundary()`. Raw-mesh flood-fill is forbidden: it can split a crossing pocket into one-sided cap traces — the same failure class that later appears as periodic edge-count/key asymmetry in `match_edges()` (`cl_Mesh_PeriodicityFactory.cpp:863-868`; field case 2026-07-16: 12458 vs 12498 with matching node counts, because `collect_nodes` folds through `original()` while `collect_edges` counts raw edge objects).

- [x] **Rule 4 — Structure taxonomy has two confidence tiers.**
  - **Tier A — air-interior pure pockets (flip by default):** the changed cut footprint is disjoint from conductor interfaces, thin-shell sidesets, terminal/cut sidesets, outer boundaries, periodic cap sidesets, and any other structure-sensitive sideset. These flips are class-safe and high-confidence because they remove only an isolated air-only closed fragment. (Not "byte-inert" downstream — `CutData`/`CutProcessor` lists do change — but the changed footprint never touches structure-sensitive paths.)
  - **Tier B — structure-adjacent or ambiguous pockets (count/report only initially):** the changed footprint touches interfaces, shells, boundaries, periodic caps, terminal curves, or another generator's nonzero cut support. Rule 1 still preserves the class, but changed cut faces may exercise the historical side-curve, interface, one-sided seam, and cut-case enumeration paths. Enabling Tier B flips requires validation cases per structure type.
  - The footprint predicate must classify the **changed cut footprint** (cut faces / duplicated nodes derive from surrounding elements, `cl_CutData.cpp:264-320`, `:443-665`), not merely the sideset membership of `∂K` edge IDs; periodic cap sidesets are named explicitly (`CutProcessor` handles them separately from ordinary boundaries, `cl_CutProcessor.cpp:21-28`, `:53-68`).
  - The per-generator log reports Tier-B candidates: count, `∂K` size, periodic-cap contact, and overlap with other generators. This is also the cheap diagnostic for deciding whether a periodic mismatch trace is a removable pocket or part of the essential cut.

- [x] **Rule 5 — Monotone, deterministic, and terminating.**
  Fire only strictly support-reducing pure-pocket flips. Process current components in ascending minimum representative-node index. **After every accepted flip, rebuild the zero-graph components from the updated coefficients** — zeroed boundary edges merge components, and scanning a stale component list can miss newly merged pure pockets or manufacture inadmissible fragment flips (Codex M1, required not optional). Termination is guaranteed because support is a finite non-negative integer and strictly decreases per accepted flip. For a nonempty pure pocket the admissible sign is unique; empty-boundary components are no-ops and are skipped.

- [x] **Rule 6 — Pipeline position and gates.**
  The pass runs inside `clean_spfa()` after the rectification unit gate (`cl_Cohomology.cpp:599-608`) and before `mMesh->unflag_everything()`, while the same periodic slave flags and `original_edges()` domain are still available. The unit gate reruns after pocket removal. Per generator, log: zero-graph component count, pure Tier-A pockets removed, Tier-B pockets skipped, periodic-cap contacts, and **other-generator overlaps** (generators are independent classes, but `CutProcessor` builds cut sets jointly, `cl_CutProcessor.cpp:43-103`; pockets overlapping another generator's support are Tier B until validated — Codex M4).

- [x] **Rule 7 — What the pass must not be expected to fix.**
  If, after pocket removal, the essential cut still terminates on a periodic cap, the `match_edges` mismatch should remain visible. A shared cap trace with a one-sided cut body is a seam-continuation defect in cut construction (one-sided periodic pair handling exists at `cl_CutProcessor.cpp:704-764`), not a representative-choice problem. That gets its own plan; this pass must not hide it.

## §3 Open questions (decisions before implementation)

- [x] **O1 — Code predicate for Tier A/B.** Define the exact source-level predicate for the changed cut footprint: interfaces, thin-shell sidesets, generated cut/terminal sidesets, outer boundaries, periodic cap sidesets, and other-generator overlap. File:line citations required during implementation.
- [x] **O2 — Outer-boundary pockets.** Draft recommendation: Tier B initially, because boundary sidesets carry BCs and cut motion against Dirichlet/symmetry surfaces touches BC handling. Revisit only with a validation case.
- [x] **O3 — Partial improvers.** ~~Allow flips that shrink but do not erase?~~ COLLAPSED by Codex audit: under the current zero-component + unit-preserving rules, admissible support-reducing flips are exactly pure pockets — partial improvers do not exist in this formulation. Any generalized mode (arbitrary subsets, weighted gains) is a separate future design with its own tie-break and validation, not a switch on this pass.
- [x] **O4 — Primitive reuse from `./archive/graph`.** The pass needs deterministic connected components on an edge-filtered quotient graph and a component-boundary scan. **Update 2026-07-15/16:** the parallel Fable session has rewritten and verified the archived Tarjan pockets algorithm (see memory `project_archive_graph_algorithms`; MV matching has pre-existing `find_path` aborts and is not needed here). Evaluate the revived primitive first; otherwise a small local BFS. Either way, components are rebuilt after each accepted flip (Rule 5).

## §4 Audit trail

- Exchange thread: `tmp/ai_exchange/corc_clean_divergence.md` (SPFA implementation + audits; pocket-removal rules audit under `# CODEX RULES AUDIT`, 2026-07-15 21:20 PDT — M1–M6 all incorporated above).
- Related plans: `todo/thin_cut_nonunit_rectification_implementation.md` (T1–T9; Phase 4 min-cost representative), future seam-continuation plan (Rule 7).
- Historical note: the year-ago flood-fill attempts predate rectification — with `|c| ≥ 2` coefficients, inside/outside labeling was ill-posed. Post-SPFA (`|c'| ≤ 1` + exact closedness) the flood-fill is well-defined; the idea was right, the prerequisite was missing.
