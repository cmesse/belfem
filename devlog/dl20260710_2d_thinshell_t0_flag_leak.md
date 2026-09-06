# Devlog 2026-07-10 — T0 First Run: Facet-Flag Leak Tears 2D Tape Tips

**Date:** 2026-07-10
**Topic:** First-ever 2D thin-shell run (T0 tapestack deck): cohomology abort
root-caused to a facet-flag leak in `fix_facet_masters`; three fixes applied.
**AIs involved:** Claude (diagnosis, fixes), Codex (audit — key dissent), Grok
(math/literature confirmation), Christian (instrumented runs, probes, approvals)
**Claude Confidence:** high (root cause; every measurement explained), pending
build+run (fix verification)
**Codex Audit Confidence:** medium-high (dissent on the intermediate hypothesis
proved decisive)
**Literature References:** Pellikka et al. 2013 (topology) — terminal pairing via
∂Σ cycles / relative 1-cycles in H₁(S,∂S)

## Summary

The T0 deck (11-tape 2D stack, first exercise of the 2D thin-shell pipeline)
aborted in `Cohomology::updatekGeneratorsFromHomology` with a Matrix bounds error.
A four-probe instrumented debugging session (pairing-matrix dumps, chain boundary
walk, tip-count probe) traced the failure through five layers to a single scratch
flag leak. Fixes applied with Christian's approval; build and verification run
handed to Christian.

## The causal chain (each link measured)

1. `MaxwellFactory::fix_facet_masters` flags every interior facet whose master and
   slave share a domain type, and unflags only those reached by a propagation
   queue seeded from domain-type-CONTRAST facets. An all-air 2D thin-shell model
   has no contrast facets → empty queue → **function exits with every interior
   facet flagged**.
2. The 2D tip test in `CutFactory::duplicate_nodes_on_face_sidesets` assumes only
   tape facets are flagged; tape tips counted their (flagged) connector-line
   facets too (measured: node0-facets = 2 for outer tapes, 3 for middle — tape +
   1–2 connectors) → `tCount == 1` never fired → **zero tips protected**.
3. All 22 tape tips duplicated → every tape silently torn open (gap G1 made real).
4. `close_terminal_loops` then "closed" each terminal curve through disconnected
   node sets — measured: 208 distinct chain edges but **4 nonzero-boundary nodes**
   per chain (2 tips × 2 sides), the fingerprint of two open polylines.
5. Open suggested chains pair with the computed cohomology basis only by
   accidental edge overlap (`Cochain::operator()` is a shared-edge sum,
   representative-independent only for cycles — Grok+Codex confirmed; Pellikka
   pairs terminals as ∂Σ cycles). Measured tTemp: adjacent-pair difference loops +
   sporadic anchors + one pure-air zero row; component {4,5,6} unanchored
   (`c4+c5+c6 = 0`) → rank 10 of 11 → transform overflow → abort. ≤5 conditions
   passed by adjacency luck.

**Codex's dissent was the turning point:** it refused the intermediate hypothesis
("`suggest_Homology` lacks a closure step") by finding that `close_terminal_loops`
already implements the master-forward/slave-return closure — which forced the
contradiction (closed chains + genuine cocycles ⇒ provably full rank) that pushed
the hunt upstream to the flag state. `suggest_Homology`, `close_terminal_loops`,
and the cohomology module are all correct.

## Changes Made

- `src/fem/maxwell/cl_MaxwellFactory.cpp` — `fix_facet_masters`:
  `mMesh->unflag_all_facets()` after the propagation loop (worklist flags are
  scratch; previously leaked whenever the queue never reached a facet).
- `src/homology/cl_CutFactory.cpp` — `duplicate_nodes_on_face_sidesets` (2D):
  `mMesh->unflag_all_facets()` before flagging tape facets (defense in depth), and
  a `BELFEM_ERROR` when an open terminal curve finds ≠ 2 tip nodes (0 allowed for
  closed loops) — this failure class now announces itself in one line.

## Input-format findings (T0 deck, first 2D thin-shell deck in the tree)

- 2D `topology/curves` syntax: plain sideset list (`1 : 5, 6 ;`), NOT the 3D
  `A @ B` intersection form; first listed sideset must belong to the tape
  (association via `sideset_a`, cl_MaxwellFactory.cpp:2507).
- Current BC: `input curves : [1],[2],…` — each bracket group is its own
  condition (`get_id_groups`); omit `output curves` in 2D (parser doubles the
  input list). Deck: `cmake-build-debug/tapestack/input.conf`.

## Open Questions

- Verification run pending (Christian): bracketed [1]..[11] should pass cohomology
  with rank 11; then re-test the unbracketed union form — with closed chains the
  union is a valid cycle and may work like 3D corc does.
- **New B9:** `fix_facet_masters` orientation propagation never runs on all-air 2D
  meshes (empty seed queue) — tape facet orientation stays raw gmsh order;
  interacts with A4 sign conventions.
- Residual polish: rank-deficiency `BELFEM_ERROR` in
  `updatekGeneratorsFromHomology`; `Chain::getBoundary()` is nullptr for isBound
  chains and `addSimplexToChain` never maintains boundaries — future cycle
  invariants need a manual boundary walk.
- Christian's diagnostic prints in cl_Homology.cpp / cl_Cohomology.cpp can be
  removed once the fix is verified.

## Files Updated

- src/fem/maxwell/cl_MaxwellFactory.cpp
- src/homology/cl_CutFactory.cpp
- todo/2d_thinshell_todo.md (B8 rewritten with root cause; B9 added)
- cmake-build-debug/tapestack/input.conf (curves + BC sections, earlier today)
- devlog/dl20260710_2d_thinshell_t0_flag_leak.md (this file)
- devlog/README.md
