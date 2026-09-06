# Periodic Seam Phase Two: XOR-Once + Facet-Mediated Combo-III Ties — COMPLETE

**Date:** 2026-07-16
**Purpose:** Execute Christian-approved phase two of `todo/periodic_cap_cut_emission.md`:
make the periodic edge rebuild duplicate-aware and fix the one-sided in-plane cut
emission, so the double-layer corc passes the solve-phase periodic rebuild.
**Modules:** src/homology (CutProcessor), src/mesh (PeriodicityFactory)
**Exchange thread:** `tmp/ai_exchange/corc_clean_divergence.md` (GROK EDGE-TIES,
CODEX SIGN-ASYMMETRY, phase-two entries)

## What landed (iterative, each step run-validated on corc)

1. **XOR-once quotient self-cancel** (`CutProcessor::collect_facets`): for each
   periodically identified face pair with BOTH copies claiming a thin-cut face,
   unflag both — the quotient analog of the interior same-Face self-cancel. A single
   claim stays where it fired (one-sided emission + straight ties = jump I once).
   **Validated by cancellation arithmetic: the seven cuts canceled 149 pairs —
   exactly the prior run's 149-edge count mismatch.**
2. **Facet-mediated periodic edge matching** (`match_edges` rewritten): the
   index-aligned facet lists from `map_facets` are the sector-resolved carrier;
   within a matched facet pair, edges correspond through their originals (unique per
   facet). Raw-profile classification per Grok's EDGE-TIES combo-III algebra:
   - **pure pairs** (each side both-dup or both-orig) tie DIRECTLY — the abstract
     current cancels in the Whitney relation `h = φ₀ − φ₁`, the jump stays in the
     node hanging (my earlier "direct ties erase the jump" was node-level reasoning
     wrongly lifted to edges — Grok corrected it);
   - **half-cut edges** (exactly one dup endpoint) are the genuinely affine class:
     untied by policy, censused (`#tie-halfcut`) — DOF-inert in the φ-region;
   - **trace twins** (1:2 realizations along the cut boundary curve, BOTH
     directions): Whitney-identical; first proposal tied, twin aliased — the
     symmetric alias also protects the first tie's direction alignment;
   - direction alignment through originals with a debug-tier topology assert;
   - hard aborts only for true inconsistencies (facet mismatch / uncovered edges).
3. **E1 compaction:** after matching, `master_edges()/slave_edges()` hold aligned
   tied pairs only — `set_entity_dependencies`, BFM save, and crosslink keep their
   matched-pair contracts.

## Final census (double-layer corc, after the parity-aware alias fix)

```
pre-cut:     tied 9318 | everything else 0                      (bijection)
solve-phase: tied 8565 ( 4333 through cut duplicates ) |
             halfcut 1678 | alias 1305 | facetMismatch 0 | uncovered 0 + 0
```

Post-audit hardening (Codex phase-two-final): the alias skip now carries a
release-active Whitney-twin invariant — and it FIRED on first contact: some trace
twins are half-cut twins (circulation differs by I), and with the half-cut twin
arriving first the pure twin had been silently skipped. The parity-aware,
order-independent dedup recovered 27 previously-lost ties and reclassified 48
half-cut twins into the censused affine class. Bonus validation: the intermediate
corc.bfm (saved after a successful periodic rebuild) reloads cleanly — the
compacted pair lists survive BFM serialization (the Codex-E1 consumer path);
parked as corc.bfm.presave.

`match_edges` passes. `set_entity_dependencies` passes (original-aware orientation
asserts hold for all mixed ties). Kernel creation and DofManager assembly pass.
**The run reaches the first Jacobian assembly** and stops at
`maxwell::h_calc: "MaxwellData has not been initialized"` — the KNOWN D1 defect of
`todo/maxwell_kernel_collapse_plan.md` (null-material MaxwellData on air blocks,
diagnosed 2026-07-13). The periodic seam is no longer the blocker; the baton passes
to the kernel-collapse plan.

> **Correction 2026-07-19 (Claude):** the D1 attribution above is wrong — the D1 air
> gates were already in tree and functional. The real cause is **D13** (see the
> kernel-collapse plan): the silently-landed R6 flip dispatched `h_calc` on the
> ThinShell·LookupAlloy branch, but thin shells are FEM sidesets and never build a
> MaxwellData. Mitigated 2026-07-19 by reverting that dispatch line to legacy
> `h_alloy`/`h_alloy_t`; the double-corc run should now clear this point.

## Formulation record (why the ties are what they are)

Grok EDGE-TIES (thread): combo (III) — straight node ties + direct edge ties — is
the unique self-consistent structure under XOR-once for pure pairs; combos (I)
affine-edges and (II) crosswise-same-sign are algebraically inconsistent; mirrored
emission + straight ties cancels the jump to ZERO (not 2I). Half-cut edges are the
only place the current survives on an edge DOF; expressing that affinely would need
a DofData extension (mixed Edge+Node sources) — avoided by the untied policy, valid
while such edges stay in the φ-region (the census guards this).

## Open

- Final Codex audit of the accumulated phase-two diff (in flight at devlog time).
- Half-cut closure verification if a conductor-side case ever appears.
- The diagnostics (`probe_duplicate_symmetry`, `#tie-*`/`#xor-once`/`#pocket*`
  probes) are still in the tree for Christian's field validation; strip or gate
  before production if desired.
- Next blocker: ~~D1~~ D13 in `todo/maxwell_kernel_collapse_plan.md` (corrected
  2026-07-19; mitigated same day by the TS-alloy dispatch revert).

## Files touched

- `src/homology/cl_CutProcessor.cpp` (XOR-once pass)
- `src/mesh/cl_Mesh_PeriodicityFactory.cpp` (`match_edges` facet-mediated rewrite,
  aliases, compaction, probes), `cl_Mesh_PeriodicityFactory.hpp` (signature)
- Earlier same-arc: `cl_Cohomology.cpp/.hpp` (SPFA hybrid + pocket census),
  `fn_Graph_spfa.{hpp,cpp}` — see `dl20260715_*` and `dl20260716_*` devlogs
