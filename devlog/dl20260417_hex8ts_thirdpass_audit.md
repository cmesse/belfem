# Devlog 2026-04-17 — HEX8TS Third-Pass Audit

**Date:** 2026-04-17
**Topic:** Read-only third-pass audit of Claude's H1-H5 hypotheses about HEX8TS side-connector failure modes, independently re-verifying Codex's two 2026-04-17 findings and adding one detail-level observation about *which* layer-edge view the wrap's inner edges hang to.
**AIs involved:** Claude
**Claude Confidence:** high on concurrence, medium (~60%) on the novel detail
**Codex Audit Confidence:** prior Codex passes (21:46 and 22:30) already on file
**Literature References:** none invoked; source-only audit

## Summary

User left a standing audit request in `todo/ai_exchange.md` (Claude 2026-04-17 16:30:00 PDT) and went offline for the weekend. Codex already posted two audit responses on 2026-04-17 (21:46 and 22:30). This session re-executes the audit as a fresh independent pass, so that three independent reads now exist on the same question. I concur with Codex on every H1-H5 verdict and add one observation at the detail level that both prior passes appear to have under-examined.

## Key Findings

### Concurrence with Codex (independently re-verified)

- `ThermalFactory::create_thermal_kernel()` omits `LeftCoating`/`RightCoating` from its block switch at `src/fem/thermal/cl_ThermalFactory.cpp:63-104`, so wrap nodes have no thermal DOFs and Joule heating is not coupled back.
- `h_tb_t` samples a vertical end-face, not the 4 inner-face nodes, at `src/fem/maxwell/matrices/mt_maxwell_h.cpp:1863-1876`. Additionally the Left and Right samplings are at *different* relative `k` positions (`k+1` vs `k`), which biases the averaged T whenever T has an along-curve gradient.
- H3 wrap material assignment is structurally correct (`cl_MaxwellFactory.cpp:885-888` plus `materials.first()==materials.last()` at `cl_ThinShellFactory.cpp:1938`).
- No mutual-exclusion guard between side-connector creation and a buffer-cut-rerouting path exists, but there is also no implemented buffer-cut-rerouting in the current `CutFactory`/`MaxwellFactory`. So §6.4 of `side_connector_modeling.md` is a future concern, not a present runtime one.
- `mCreateSideConnectors` is hard-true at `cl_ThinShellFactory.hpp:66` with no setter.
- Silent connector suppression (when `materials.size() < 3` or first/last labels mismatch) at `cl_ThinShellFactory.cpp:1933-1938` is real and unreported.

### Novel observation — inner-edge hanging target

Codex 21:46 verifies that `collect_inner_nodes_and_edges` plumbs `set_sources` correctly: one source, weight +1, node order preserved via `mNodeMap`. That is true in the narrow sense. But the question "which of the two layer-edge views does the hanging source *point at* when the layer has duplicates?" was not explicitly adjudicated.

From `link_elements_with_edges` (`cl_ThinShellFactory.cpp:1568-1594`) the convention is:
- At a `hasDuplicates` layer L, `Edges` (non-duplicate) = lower block's top view, `EdgeDuplicates` = upper block's bottom view.

From `collect_inner_nodes_and_edges` (`cl_ThinShellFactory.cpp:2282-2362`) and `create_side_elements` (`cl_ThinShellFactory.cpp:2453-2539`):
- Wrap row j's Lo inner horizontal = `tLo->InnerEdgesHorizontal` ← always `tLayer->Edges` → **block j-1** view.
- Wrap row j's Hi inner horizontal = `tHi->InnerEdgeDuplicates` if non-empty → **block j+1** view; else `Edges`.

Net effect: when layers have duplicates, the wrap row j's inner-face horizontal edges are tied to blocks j−1 and j+1, never to block j (its physically adjacent block). In a Hastelloy-buffer-HTS-buffer-Hastelloy tape, 6 of 10 inner-horizontal ties point at a non-matching (non-Hastelloy) material's view.

This is unlikely to fully explain "wrap carries no current," because vertical edges (1/3/5/7) are all free and axial `J` is driven by the ζ-gradient of `H_η` through those. But at the level of the continuous-model consistency it is a discrete H-tangential jump inside a uniform-material wrap. Worth empirical verification; full argument and diagnostics are in the third-pass entry I appended to `todo/ai_exchange.md`.

## Changes Made / Proposed

No source edits (read-only audit per `AGENTS.md`).

Artifact changes:

- Appended `CLAUDE 2026-04-17 22:01:31 PDT` entry to `todo/ai_exchange.md` with concurrence on Codex 21:46/22:30 findings and the novel inner-edge-target observation, plus four Monday diagnostics ordered by signal-to-effort.

## Open Questions

- Is the wrap's inner horizontal edge wiring deliberate? (i.e., is there a physical rationale for tying to blocks j-1/j+1 rather than block j, or a later pass that re-points one side?)
- Does the Hast-Hast-Hast sanity test in diagnostic (3) of the ai_exchange entry actually isolate the issue?
- Does the wrap's behavior change under the two confirmed fixes (thermal activation + `h_tb_t` node sampling) alone, before attempting any deeper changes?

## Files Updated

- todo/ai_exchange.md
- devlog/dl20260417_hex8ts_thirdpass_audit.md (this file)
- devlog/README.md
