# Session Handoff: Double-Layer CORC — D1 Fix and First Timestep

**Date:** 2026-07-16
**Purpose:** Hand off to a fresh session. The periodic-seam campaign is complete; the
double-layer corc now reaches its first Jacobian assembly and dies at the KNOWN D1
defect of the Maxwell kernel-collapse plan. This file contains everything a new
session needs; the full trail lives in the documents listed at the bottom.
**Status:** ✅ SPENT / CLOSED 2026-08-09 — this brief has been consumed. Every
premise it hands off has since changed:

- **§3 (the immediate task) is done.** D1 of `maxwell_kernel_collapse_plan.md` is
  fixed: `Calculator::allocate` now builds `MaxwellData` only for non-Air blocks
  that exist in the peer dof manager (`cl_FEM_Calculator.cpp:1111-1137`), so the
  "MaxwellData has not been initialized" abort on air blocks is gone (the guard
  survives as a debug `BELFEM_ASSERT`, `cl_FEM_Calculator.hpp:3346`).
- **§2 is obsolete.** The whole 2026-07-15/16 campaign is committed
  (`ce6a0e8b` and neighbors); `src/math/graph`, `src/homology` and `src/mesh`
  are clean in the working tree.
- **§5.1 is done.** All seam probes (`#pocket*`, `#xor-once`, `#tie-*`,
  `probe_duplicate_symmetry`, `#probe-sym`) have been stripped — zero hits in
  `src/homology` and `src/mesh`.

Residual §5 items live on elsewhere and are NOT lost:
T8 SPFA test port + single-layer corc regression → `debt_register.md` DR-23 and
`thin_cut_nonunit_rectification_implementation.md` T8/T9; half-cut conductor watch
→ `periodic_cap_cut_emission.md` (closed, open-question section); Phase-4 L1
min-cost representative → `thin_cut_nonunit_rectification_implementation.md`
Phase 4; helix re-run → DR-09/DR-23 run gates; stale `mPairVerdict`/Step-6c
comments → DR-26.

The §4 seam validation gate numbers remain useful as historical reference values
for any rerun that recomputes cuts on the same corc mesh.

---

## §1 Where things stand (one paragraph)

`cmake-build-debug/corc` (double-layer CORC, refined mesh after the SPFA certificate
located the coarse throat) runs the ENTIRE cut + periodicity machinery cleanly:
SPFA certify-then-greedy rectification, pocket census, XOR-once quotient emission,
facet-mediated combo-III periodic edge ties, parity-aware trace-twin handling,
compacted pair lists, original-aware `set_entity_dependencies`. It then enters
`Controller::iterate_coupled → DofManager::compute_jacobian_and_rhs →
IWG_Maxwell::compute_mkf → maxwell::h_calc` and aborts:
**"MaxwellData has not been initialized"** — defect **D1** of
`todo/maxwell_kernel_collapse_plan.md` (null-material MaxwellData on air blocks,
diagnosed 2026-07-13, Codex-confirmed same day: air-specific, null-deref at first
initialize). The seam is no longer the blocker.

## §2 CRITICAL: uncommitted working tree

Nothing from the 2026-07-15/16 campaign is committed. Christian decides commit
granularity (suggested slices: solver / census / seam). Changed files:

- `src/math/graph/fn_Graph_spfa.{hpp,cpp}` (new) + `CMakeLists.txt` — SPFA
  difference-constraint solver (warm start, amortized negative-cycle scan,
  verified certificates)
- `src/homology/cl_Cohomology.{hpp,cpp}` — `clean()` → `mFunClean` scaffold
  (Christian), `clean_spfa()` certify-then-greedy hybrid, `rectify_greedy_sweeps`,
  `fire_node_coboundary`, `remove_cut_pockets` census + Tier-A pass
- `src/homology/cl_CutProcessor.cpp` — XOR-once quotient self-cancel
- `src/homology/cl_CutData.cpp` — T7 message split + stray-paste removal
- `src/mesh/cl_Mesh_PeriodicityFactory.{hpp,cpp}` — `match_edges` facet-mediated
  rewrite (combo-III ties, half-cut policy, parity-aware aliases, E1 compaction),
  `collect_edges` fallback comments, `probe_duplicate_symmetry`
- `src/mesh/cl_Mesh_Periodicity.cpp` — original-aware orientation asserts +
  authoritative-pointer pairing in `set_entity_dependencies` (pre-existing edits
  extended during the campaign)
- `todo/`, `devlog/`, `tmp/ai_exchange/corc_clean_divergence.md` — documentation

## §3 Immediate task: D1

- **Plan:** `todo/maxwell_kernel_collapse_plan.md` — D1 = air blocks construct a
  `MaxwellData` with a null material; first `initialize` call dereferences /
  asserts. D2–D5 are siblings (D2 thermal-side material binding, D3 β cache-slot
  aliasing, D4 UserDefined overload parity, D5 peer-element linking) — check the
  plan's R-step sequence and the §4.1 run protocol before editing; the plan may
  already prescribe the fix shape.
- **Reproduce:** `cd cmake-build-debug/corc && ../bin/hphirun` (dies in the first
  Picard/Newton assembly).
- **Fast iteration trick:** `corc.bfm.presave` (157 MB, in the corc dir) is a saved
  mesh WITH cuts and seam ties baked in. Rename to `corc.bfm` and the run reloads
  it, skipping the whole cohomology/cut/periodic phase (~minutes) and going straight
  to kernel + assembly — ideal for D1 iteration. Delete/rename away to exercise the
  fresh path again.
- **Related memory:** MaxwellData helper history in
  `todo/maxwell_kernel_collapse_plan.md` + memory notes (R5 confirmed, R6 LookupAlloy
  flip queued; B7 qhist fix closed).

## §4 Seam validation gates (for any rerun that recomputes cuts)

The `#`-probes are still in the tree (deliberately, for Christian's field
validation). Expected fresh-path numbers on the CURRENT corc mesh:

```
#xor-once   : 12+17+29+14+17+23+37 = 149 canceled quotient face pairs (7 cuts)
#tie-summary ( pre-cut )    : tied 9318, all other counters 0
#tie-summary ( solve-phase ): tied 8565 ( mixed 4333 ) halfcut 1678 alias 1305
                              facetMismatch 0  uncovered 0 + 0
#pocket-seam : per-generator support ~24k-53k ( sparse greedy ), NOT ~226k
```

Release-active invariants that must stay quiet: Whitney-twin alias assert,
facet-mismatch/uncovered abort, unit gate after rectification, tie-consistency in
`set_entity_dependencies`. If any fires on a NEW mesh, believe it — every one of
them found a real defect during the campaign.

## §5 Open items (tracked, not blockers)

1. **Probe cleanup** (after Christian's validation): `#pocket*`/`#pocket-seam`
   (cl_Cohomology.cpp), `#xor-once` (cl_CutProcessor.cpp), `#tie-*` +
   `probe_duplicate_symmetry` + `#probe-sym` (cl_Mesh_PeriodicityFactory.cpp).
   Gate behind a flag or strip.
2. **Half-cut conductor watch:** half-cut edges are untied by policy — valid while
   they stay in the φ-region (DOF-inert). If `#tie-halfcut` ever reports edges on
   conductor caps, the affine (jump-offset) tie needs a DofData extension (mixed
   Edge+Node sources). Census guards this.
3. **T8:** port the SPFA scratch tests (`test_spfa.cpp`, session scratchpad —
   recreate from `devlog/dl20260715_spfa_clean_implementation.md` if gone) into the
   repo suite.
4. **Phase-4 L1 min-cost representative:** deliberate deferral (both auditors);
   complexity table in the exchange thread. Delivers guaranteed-sparse cuts +
   weighted cap avoidance.
5. **Single-layer corc regression:** the archived pair at the build root is stale
   (input names vertex IDs absent from the mesh) — needs a matching mesh/input to
   re-validate; scaffolding in `cmake-build-debug/corc_single_regression/`.
6. **helix re-run** with the final seam code (its earlier success predates
   the facet-mediated matching): expect unchanged-or-better; cut pockets should be
   reduced by sparse rectification + Tier-A flips.
7. Stale `mPairVerdict`/Step-6c comments in cl_CutSet.cpp still mislead audits
   (memory: superseded by Step 9; recommend "8j" comment cleanup).

## §6 Key documents

- **Plans:** `todo/periodic_cap_cut_emission.md` (✅ seam, full mechanism +
  formulation record), `todo/cut_pocket_removal_rules.md` (✅ rules + field
  outcome), `todo/thin_cut_nonunit_rectification_implementation.md` (T1–T7 done,
  pivot recorded), `todo/maxwell_kernel_collapse_plan.md` (**next**)
- **Devlogs (the story, in order):** `dl20260715_corc_clean_divergence.md`,
  `dl20260715_spfa_clean_implementation.md`, `dl20260716_spfa_warmstart_cycle_scan.md`,
  `dl20260716_census_density_pivot.md`, `dl20260716_periodic_seam_phase_two.md`
- **Exchange thread (ephemeral, sweep after distilling):**
  `tmp/ai_exchange/corc_clean_divergence.md` — all Codex/Grok verdicts
  (GROK EDGE-TIES combo-III algebra and CODEX PHASE-TWO FINAL are the load-bearing
  entries if anything in the tie logic is ever revisited)

## §7 Conventions for the next session

Per `doc/ai_collaboration_protocol.md` / `AGENTS.md`: read-only until Christian
approves edits; Codex audits + Grok third voice for anything correctness-critical
(launch `codex exec` from the REPO ROOT — a cmake-build-debug cwd sandboxes it away
from `tmp/ai_exchange/`); calibrated confidence on non-trivial claims; devlog +
README index for meaningful sessions; Christian runs builds himself unless he
explicitly grants run permission (he granted it during this campaign for the corc
iteration loop — re-confirm, don't assume).
