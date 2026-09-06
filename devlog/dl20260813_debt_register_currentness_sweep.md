# Debt Register Currentness Sweep — 24 Open Rows Verified Against the Tree

**Date:** 2026-08-13
**Purpose:** Sweep every unstruck debt-register row against the current source to find
rows whose debt is already paid, and to correct rows the tree contradicts.
**Module:** `todo/` (verification touched `src/fem`, `src/mesh`, `src/homology`,
`src/circuit`, `src/fvm`, `tests/`)

## Method

Six read-only verifiers ran in parallel, grouped by module, each given the row's own
claims and one rule: **the evidence is the tree, never the plan and never the devlog.**
Every finding below carries a `file:line` a verifier actually read. Nothing was compiled
or run; no source was edited.

## Closed

**DR-28** (periodic, thin-cut QB) — closed as already-fixed. The row asked for a cheap
step-1 diagnostic separating "feasible but greedy failed" from "genuinely infeasible".
That diagnostic is the **default path** and has been since the SPFA certifier landed:
`Cohomology::clean_spfa()` (`cl_Cohomology.cpp:491`) calls
`graph::spfa_difference_constraints` at `:578` *before* any greedy sweep. Infeasible
produces a hard error naming the negative cycle, its edge count, node pairs, nearest
element and an `error.exo` dump (`:658-664`); feasible-but-stalled produces the distinct
warning at `:682` plus a dense fallback with a re-solve guard (`:686-691`). All three
constructors bind `mFunClean = &Cohomology::clean_spfa`. The old undifferentiated message
survives only in the unselected `clean_greedy()`.

## Rows whose meaning changed

**DR-30 (postproc aura stale-T) — the row as written is refuted; a different, narrower
defect survives.** T on ghost nodes *is* written by every thermal solve
(`DofManager::solve()` → `FieldData::distribute`, `cl_FEM_DofManager.cpp:1032-1038`),
and the comm table's node set covers the aura layer
(`cl_Mesh_Distributor.cpp:306-317,353-359`); Anderson writes fields before that
distribute, so it does not reopen the hole. What actually survives: T is never in
`mPostprocessorSourceFields` (`cl_FEM_DofManager.cpp:1292` distributes only
`{edge_h, face_h}`/`{phi}`), so correctness rests entirely on a preceding thermal solve;
and `cl_FEM_Calculator.cpp:3003-3012` silently keeps the thermal peer's *previous*
element when `element_exists()` is false, while `compute_superconductor` re-checks only
`block_exists` (`cl_MaxwellPostprocessor.cpp:674,:747`). Severity drops from "wrong on
every boundary element" to "wrong where maxwell/thermal group membership diverges".

**DR-19 (side-connector) — the path is unblocked, serial AND parallel.** The
`edge coating : on` path was read end to end: no WIP stop, early return, stub or
todo-abort from input parse (`cl_MaxwellFactory.cpp:2904-2921`) through
`create_side_connectors` (`cl_ThinShellFactory.cpp:420-678`) to a real IWG dispatch
(`cl_IWG_Maxwell.cpp:421-428`) and full kernel bodies (`mt_maxwell_h.cpp:264-320`).
**This sweep first reported multi-rank as broken; that was a false positive, retracted
the same day — see the section below.**

**DR-02 (kernel collapse) — not purely a run.** The plan's R2d converged-at-clamp
diagnostic never landed: `Calculator::T_clamped()`/`rho_clamped()` exist
(`cl_FEM_Calculator.hpp:296,:299`) with zero consumers tree-wide — roughly ten lines of
Controller-side print. The Gate A amendment is otherwise true, but its cited coordinates
were `bc578b5e`'s; at HEAD the env-gated dump blocks are `:973-987` and `:1671-1684`.

DR-02 is the **second pattern-break** of the register's life — a row found *worse* than
recorded (the first was DR-45(a) on 2026-08-11).

## The sweep's own false positive — retracted the same day

This sweep initially reported a third finding: that **P9(a) was a live code gap and
multi-rank `edge coating : on` was broken**, on the evidence that
`Distributor::send_thinshell_data` / `receive_thinshell_data`
(`cl_Mesh_Distributor.cpp:1686,:1732`) contain zero occurrences of
`side_connector`/`coating`, so non-root ranks rebuild ThinShells without
`side_connector_blocks()`.

Christian challenged it within the hour — *"I am running a problem with multiple ranks.
How can edge coating be broken?"* — and the code refutes it. The grep premise is true;
the conclusion does not follow. The connector record is **deliberately** absent from the
shell record, and everything derived from it travels by two other explicit channels:

1. `cl_MaxwellFactory.cpp:2160-2200` broadcasts `tCoatingBlocks` and `tCoatingTypes`,
   carrying the walls' true `Left/RightCoating` domain types to every rank. Its own
   comment states the reason outright: *"the types must travel explicitly since non-root
   ranks rebuild ThinShells without the connector record"* (`:2162-2164`).
2. `BlockData::link_thin_shell_facets_parallel()`
   (`cl_FEM_DofMgr_BlockData.cpp:522-680`, dispatched whenever `comm_size() >= 2`,
   `:447-466`) builds `tSideData` from `side_connector_blocks()` on rank 0, `share`s it,
   and every other rank re-establishes each wall's recovery facet and its master
   reference — with aura-aware skips and always-active `BELFEM_ERROR` checks for owned
   walls.
3. P9(b), the wall thickness off-root, was already fixed on 2026-08-09 via
   `proto::GroupData::mThickness`.

Every consumer of `side_connector_blocks()` is root-only or serial-only:
`cl_FEM_Kernel.cpp:471,496` (root-only `partition_mesh`), `cl_MaxwellFactory.cpp:1085,1137`
(inside `if mCommRank == 0`), `:2176,2189` (root-only builder feeding the broadcast),
`cl_FEM_DofMgr_BlockData.cpp:491` (serial path), `:572,582` (root-only builder feeding the
share). So the empty accessor off-root is a cosmetic incompleteness of the ThinShell
object, not a functional gap.

**The lesson, which is worth more than the finding would have been:** the verifier proved
an *absence in one function* and the sweep promoted it to a *functional gap* without
tracing the alternative routes. That is the standing audit false-positive shape — check
that the thing is really unavailable before flagging it as missing — and a grep-shaped
absence is its most convincing costume. It survived because the plan file's own P9(a)
text says the same thing, which made the inference look corroborated when it was merely
repeated.

## Recommended for closure — Christian's ruling, not a fact

- **DR-16** — the empty-seed-queue mechanism is confirmed, but 2D no longer relies on it:
  `cl_MaxwellFactory.cpp:1260` sets `tPropagate = (dimensions == 3)` and `:1296` gates
  the BFS on it, with `:1252-1259` explaining that in 2D the thin-shell cut pipeline needs
  a uniform per-tape master side, so selective pairwise flips could only break it.
  Deliberate and documented — a by-design closure like DR-59.
- **DR-17** — both PENTA6TS circulation tests are written, enabled and fast-labelled
  (`tests/fem/test_InterfaceOrientation.cpp:440,:478`); only the run and sign-off remain.
  A DR-42/49-style strike.
- **DR-31** — `0671b29a` is an ancestor of HEAD and the working tree matches it exactly;
  the stride fix is entity-agnostic, so `edge_h` **and** `face_h` are covered by one path.
  Only the ≥2-rank order-2 run remains. A DR-52/64/66-style strike.

## Ready-to-execute cleanups (change set pinned in the row)

- **DR-26** — `mPairVerdict` has **zero reads** tree-wide (8 hits, all declaration,
  comment, clear or write). Full deletion set recorded in the row: `cl_CutProcessor.hpp:66-72`
  and `:159-163`, `cl_CutProcessor.cpp:82` and `:724-787`, `cl_CutSet.hpp:26-41`, plus the
  stale comment text at `cl_CutSet.cpp:105-107,112-119`.
- **DR-13** — one-line clamp or assert at `cl_Vertex.cpp:224` closes it. A caller audit
  showed no reachable overflow today (every site passes 0 or a `uint16_t`-derived size);
  the assert cited by the row is at `cl_Vertex.hpp:270`, not `.cpp:270`.

## Confirmed still open, no change in kind

DR-06, DR-11 (the PID controller is an **iteration-cost** controller — no ΔT-per-step cap
exists anywhere, so the proposal is absent rather than superseded), DR-15, DR-18, DR-22,
DR-27 (no guard exists; in fact no free-cut λ machinery exists in the FEM layer at all),
DR-29 (no worklist; #10's stale-adjacency case pinned to the conditional at
`cl_SimplicialComplex.cpp:1176` versus the unconditional removal at `:1221-1222`),
DR-37 (both halves absent; `UserMaterialTemplate.cmake`'s backend detection is dead code),
DR-39 (all three absent — and there is no netlist parser at all, so `.subckt` is a question
about something that does not exist), DR-40, DR-41 (two further build blockers found),
DR-46, DR-53, DR-65, DR-69.

**DR-23** was missed by the sweep's first filter (its ID is bolded in the table, which the
row-matching pattern did not catch) and was checked separately: both halves read exactly as
written — `tests/math/test_GraphSpfa.cpp` carries 16 `TEST`s, while a tree-wide grep over
`tests/` for `Cohomology`, `clean_spfa`, `rectify_greedy_sweeps`, `fire_node_coboundary`
and `remove_cut_pockets` returns zero files. The `Cohomology`-layer half is untouched.

## Register hygiene fixed in passing

Two main-table rows carried **unescaped `|` characters inside code spans**, which silently
breaks the markdown table for every reader: `netlist|subckt|spice` (added by this sweep)
and `acos( dot(n,b)/|b| )` in DR-69 (added earlier the same day by the DR-69 ruling — so
that row had been rendering broken since it was written). Both are now `\|`-escaped, and
every row below the table header was re-checked to carry exactly seven unescaped pipes.
