# Periodic Consistency of the Recovered Nodal Fields

**Date:** 2026-09-23
**Purpose:** On a periodic model, a node and its periodic partner carry the same physical state, but
the nodal recovery fits each from its own one-sided patch, so the Exodus fields H, B, J, J/Jc differ
across the periodic planes (visualization and the printed J/Jc_max; the solve ties the dofs). Mechanism: after each
postprocessor's gather, rank 0 replaces the two values of every periodic node pair by one combined
value — vector fields rotated through the periodic transform, scalars as they are.
**Module:** `src/fem/kernel` (`cl_FEM_Postprocessor`), `src/mesh` (read-only use of `Periodicity`)
**AIs involved:** Claude (plan), Codex `gpt-6-astra` + Grok `grok-4.7` (plan jury)
**Status:** PLAN — audited 2026-09-23 (Codex astra/high + Grok 4.7/high), amendments C1–C12 (§4.1)
supersede the conflicting text above them; O1 decided (plain mean); no source modified — ready
for implementation approval.

> **Scope guards:**
> - Postprocessing only; no change to the dofs, the solve, or the mesh.
> - One periodicity per mesh (`Mesh::periodicity()`, `cl_Mesh.hpp:1979-1990`).
> - The shared layer-interface nodes keep their last-writer semantics
>   (`thinshell_postprocessor_node_sharing.md`); this plan only makes each postprocessor's own
>   output periodic.

---

## 1. Current Behavior

Facts from a code trace (2026-09-23), citations from that report:

| Fact | Citation |
|---|---|
| node/edge/facet partners are one symmetric pointer each, set both ways | `cl_Mesh_PeriodicityFactory.cpp:731-732`, `cl_Node.hpp:346-368` |
| master/slave only from `Periodicity::master_nodes()/slave_nodes()` (index k pairs k); input `source` = master, `target` = slave | `cl_MaxwellFactory.cpp:625-646` |
| no transform stored after setup; the factory's plane bases imply X_t = T_tᵀ(T_s X_s + ΔD e_z): a rotation R = T_tᵀ T_s plus a shift along the target normal | `cl_Mesh_PeriodicityFactory.cpp:715-729`, `:1360-1432` |
| the solver slaves each slave node/edge dof to its master with weight +1, no rotation (edge integrals are invariant; nodal vector components are not) | `cl_Mesh_Periodicity.cpp:49-114`, `cl_FEM_DofMgr_DofData.cpp:2964-2975` |
| `Periodicity` and the partner pointers are set on rank 0 (the distributor ghosts sources, not partner elements) | `cl_Mesh_Periodicity.cpp:51`, `cl_Mesh_Distributor.cpp:305-336` |
| thin-shell layer nodes are paired too (`set_periodic` both ways) | `cl_ThinShellFactory.cpp:1387-1426` |
| after `synch_target_fields`, rank 0 holds every target field on the full mesh | `cl_FEM_Postprocessor.cpp:149-193`, `:381-392` |
| target fields are scalar components: `Hx, Hy, Hz, Bx, …, Jx, …, JJCx, JJCy, JJCz` | `cl_MaxwellPostprocessor.cpp:62-150` |
| decks with periodicity: translated `corc_periodic_bc`, `rutherford`, `tapestack3d(_quench)`, …; **twisted** `corc_garber_effect`, `corc_twolayer` (`examples/README.md:103-104`) | example `input.conf`s |

**Bottom line:** the partner pairs, the fields and the geometry all sit on rank 0 after the gather;
what is missing is the transform as a matrix, a rule for which components form vectors, and the
combination step.

## 2. Architecture

**Chosen: combine on rank 0 after the gather.** For each pair (m, s) selected by this
postprocessor, with R the rotation master → slave:

- vector field v (component triples `?x, ?y, ?z`): v̄ = (w_m R v_m + w_s v_s)/(w_m + w_s);
  write v_s = v̄, v_m = Rᵀ v̄
- scalar field: the same weighted mean, no rotation

with w the node's patch weight (Σ element volume over its patch, available as G(0,0) of the node
matrix — O1). For a field inside the fit's basis both one-sided values are already exact and equal,
so the mean changes nothing; otherwise it splits the difference by patch size.

**Rejected: the union patch across the periodic planes** (the sibling mechanism, with the partner's
samples moved through the transform and its vector values rotated). It is the more accurate fit,
but the partner's elements are not on the same rank in parallel (the distributor ghosts sources
only), the partner pointers exist on rank 0 only, and it would touch the patch walk that R9 just
pinned. The visualization goal does not need it.

## 3. Gap Table

| # | Piece | Needed for | Handled today? | Class | Citation / rationale |
|---|---|---|---|---|---|
| 1 | Rotation R and shift t, master → slave | vectors, validation | no | (a) from the three plane-node pairs `master_plane()` / `slave_plane()` | frames [P; Q; N] per plane as the factory builds them; R = F_tᵀ F_s |
| 2 | Validation of R, t on every node pair | safety | no | (c) | `‖R x_m + t − x_s‖ ≤ tol·L`; on failure a warning and scalar-only combining (no rotation) |
| 3 | Pair list with direction | combining | yes (rank 0) | (a) | `master_nodes()/slave_nodes()`; layer nodes via `periodic()` pointers, direction decided by the transform (O2) |
| 4 | Vector grouping of target fields | rotation | no | (c) | component triples by name suffix `x, y, z` on a common stem; everything else scalar |
| 5 | Patch weights per node | weighted mean | no | (b) → O1 | G(0,0) of the node matrix before inversion |
| 6 | Where it runs | ordering | — | (a) | end of `Postprocessor::recover_fields`, rank 0, after `synch_target_fields`, only pairs whose both nodes this postprocessor owns-or-selected |
| 7 | Duplicates on periodic planes (cut siblings) | consistency | no | (b) → O3 | a sibling of a periodic node may or may not carry its own partner pointer |
| 8 | 2D periodicity | generality | — | (a) | R is 2 × 2 (in-plane), vector triples become pairs |

## 4. Ordered Steps

- [ ] **R1** `src/fem/kernel/fn_FEM_periodic_combine.hpp` — header-only, testable without a mesh:
  given R, the pairs (index m, index s, weights), the scalar fields and the vector triples, write
  the combined values in place.
- [ ] **R2** Rank-0 setup in `Postprocessor::initialize` (once): rebuild R, t from the plane nodes,
  validate on all node pairs, build this postprocessor's pair list (both nodes among its selected
  nodes) and its vector grouping (after: R1).
- [ ] **R3** Call at the end of `recover_fields` on rank 0 (after: R2).
- [ ] **R4** Tests: R1 directly (rotation round trip, weighted mean, scalar untouched by R, 2D);
  a periodic fixture through the real recovery if the kernel chain allows it (translation and a
  90° twist, as `tests/mesh/test_Periodicity.cpp:372-390`); mutations: R replaced by the identity
  → the twisted case red; combine skipped → both red.
- [ ] **R5** Deck gates (Christian runs, or the scratch worktree with approval): `corc_periodic_bc`
  (twisted 1/3 pitch) and `rutherford` (translation): every master/slave node pair equal after
  rotation to round-off; fields away from the planes unchanged bit for bit.
- [ ] **R6** Docs: `postprocessor_recovery_theory.md` (a periodic section), the node-sharing note.


### 4.1 Amendments from the plan jury (2026-09-23, `tmp/ai_exchange/review_periodic_fields.md`)

- [ ] **C1** (Grok) Gates: the twisted deck is `corc_garber_effect` (or `corc_twolayer`), not
  `corc_periodic_bc`, which is translated (`examples/README.md:103-104`, its README). R5's
  identity-R mutation only discriminates on a twisted deck.
- [ ] **C2** (Grok) Transform: X_s = T_tᵀ(T_s X_m + (D_t − D_s) e_z), R = T_tᵀ T_s, t = (D_t − D_s) N_t,
  rebuilt from the stored plane nodes with the factory's frame construction; both frames are
  right-handed, so det R = +1. Validate only on **surface pairs** (those on the Hesse planes, within
  `BELFEM_MESH_EPSILON`); layer nodes sit off the plane by up to half the stack and must not veto R.
- [ ] **C3** (Codex, Grok) Failure contract: a degenerate frame, det R ≠ 1, a non-finite coordinate or a
  failed surface-pair check is a setup `BELFEM_ERROR` (runs once), not a warning with scalar-only
  combining, which would draw a smooth, wrong field on a twist.
- [x] **C4** *decided 2026-09-23, Christian: plain mean* (no communication). (Codex, Grok) Weights: `G(0,0)` is n_q Σ V_e over the node's patch (volume once per
  integration point), computed on the **owning** rank and discarded after inversion; rank 0 does
  not have it for remote nodes. Either a plain mean (no communication), or a snapshot plus a gather
  in `initialize` on every rank, empty ranks included, before any `recover_fields` (a collective
  after its early return would break the drain check). **Christian's call (O1).**
- [ ] **C5** (Codex, Grok) Membership on rank 0: a pair is combined only if both nodes are among the
  nodes this postprocessor wrote — rank 0's `mMyOwnedNodeIndices` ∪ `mAllOwnedNodeIndices( p ≥ 1 )`
  — captured as mesh indices in `initialize` (the node indices are repacked inside `recover_fields`).
- [ ] **C6** (Codex P1, Grok) Duplicates: iterate the authoritative `master_nodes( k ) / slave_nodes( k )`
  once (surface, both-sided cut and layer pairs are all there). One-sided cut duplicates have no
  pair (`cl_CutSet.cpp:148-158`); after combining a pair, copy the result to the selected co-located
  siblings of each endpoint that have no pair of their own, so a sibling does not keep the old value.
  Fixture with an unpaired duplicate required.
- [ ] **C7** (Grok) No 2D path: every periodic deck is 3D and the factory builds a 3 × 3; row 8 of the
  gap table is dropped. Component grouping: complete `?x ?y ?z` triples only; anything else is a
  scalar (a lone `Jz` in 2D stays scalar).
- [ ] **C8** (Codex, Grok) Order: the combine runs at the end of `Postprocessor::recover_fields` after the
  index restore, before `MaxwellPostprocessor::run`'s own tail (air `B = μ0 H` from the combined H;
  the side-connector seam copy from combined tape rims). Acceptance inspects the finished sequence.
- [ ] **C9** (Codex, Grok) Claims: for a twist the stored components differ, `R v_m = v_s` is what agrees;
  exactness holds only for fields the fit reproduces (no truncation, identifiable); the weighted
  mean is a new consistency projection, not an SPR accuracy result. Vector-averaging B (ferro μ(|H|))
  and J/Jc (J / jc(|B|)) is a visualization choice and does not commute with those maps.
- [ ] **C10** (Codex) Scope: the combine also changes the printed `J/Jc_max` (`cl_FEM_Controller.cpp:3797-3821`).
- [ ] **C11** (Codex, Grok) Tests through the real pipeline, not only the header: the R9 fixtures
  (`test_PostprocessorRecovery.cpp`) with a hand-made `Periodicity` — translation and a 90° twist,
  an unpaired cut duplicate, a non-periodic no-op; MPI: pairs with endpoints on different ranks, an
  empty rank, rank 0 owning neither endpoint (a Tier-2 rank-count suite). Expectations computed
  independently, not only pair equality.
- [ ] **C12** (Grok) "Unchanged away from the planes" is checked on a same-process snapshot before and
  after the combine; two deck runs differ by the solver's run-to-run floor.

**Revised effort:** medium. The combine itself stays small (~100 lines); the cost is in C4 (if
weighted), C5/C6 (membership and siblings), C3 (setup validation) and C11 (a periodic fixture
through the kernel plus an MPI suite).

## 5. Open Questions

- **O1** **RESOLVED 2026-09-23 → plain mean** (Christian). Original: plain mean or patch-weighted mean? Weighted needs G(0,0) kept per node (one real
  per node). Proposal: weighted.
- **O2** **RESOLVED (jury, verified)** → layer pairs are in `master_nodes()/slave_nodes()` via the backup (`cl_Mesh_Periodicity.cpp:283-297`). Original question: are they in the lists or only linked by pointers?
  If pointers only, the direction is decided per pair by testing `R x_a + t ≈ x_b`.
- **O3** **RESOLVED → C6.** Cut duplicates on a periodic plane: combine per sibling pair, or per original? To be read
  from how `set_periodic` treats duplicates before R2.
- **O4** Should J/Jc (`JJC*`) be combined as a vector (it is J scaled by jc) — yes by the grouping
  rule; the shared-node J/Jc leak is unaffected.

## 6. Effort

Small to medium: one header function (~80 lines), ~60 lines of setup and call in `Postprocessor`,
tests. The risk is in the details, not the size: the transform's direction and the plane frames,
layer and duplicate nodes, and the rank-0 assumption.

## 7. Audit Trail

- Plan jury: `tmp/ai_exchange/review_periodic_fields.md`.
