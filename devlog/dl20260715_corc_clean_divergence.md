# Cohomology::clean() Divergence on Double-Layer corc — Diagnosis (Read-Only)

**Date:** 2026-07-15
**Purpose:** Root-cause the greedy `clean()` stall on the double-layer corc case and
decide the path forward; no source edits (todo/devlog updates only).
**Module:** src/homology
**Exchange thread:** `tmp/ai_exchange/corc_clean_divergence.md` (Claude + Codex + Grok)

## Symptom

Christian's `#Probe` in `Cohomology::clean()` (`cl_Cohomology.cpp:195-322`) on
`cmake-build-debug/corc` (double-layered CORC): generators 0–5 clean to unit
coefficients in a few sweeps, but **generator 6 diverges** — non-unit count dips
4042 → 1869, then grows monotonically to a 55'511 plateau while the generator's support
balloons from 29'527 to 194'288 edges. The exact plateau trips the stall detector
`BELFEM_ERROR` at `cl_Cohomology.cpp:316` ("mesh is probably too coarse …").

## Diagnosis (three-AI concurrence, high confidence)

Two structural defects of the greedy algorithm, no outright bug:

1. **In-sweep cascading (the amplifier).** The sweep range-fors over the live
   `getSimplicesMap()` — an `OrderedMap` backed by `std::map`
   (`cl_OrderedMap.hpp:28`) — while `addSimplexToCochain` (`cl_Cochain.hpp:251-275`)
   inserts into the same map. Legal (`std::map` insertion preserves iterators), but any
   freshly minted `|c|≥2` edge with a **higher edge index** is processed in the *same*
   sweep: the sweep is order-dependent and self-feeding, matching the observed support
   flood.
2. **Sign-only endpoint rule (the misdirection).** The fired node is chosen purely by
   coefficient sign (`node(0)` if c>0, else `node(1)`, `cl_Cohomology.cpp:244`) — no
   degree/potential/progress measure. Each fire reduces the target edge by exactly 1 but
   perturbs every incident complex edge by ±1; nothing prevents excess from circulating
   around a cycle instead of cancelling.

**Ruled out** (both auditors, code-grounded): iterator invalidation on the
currently-visited entry (its own fire cannot zero it), periodic-partner double-firing of
one edge id (`:277-309` is the intended master∪slave quotient adjacency), and sign
errors (both passes use the same consistent orientation test, `:257`, `:287`).

**Secondary defects noted:** the stall detector at `:316` only catches an *exact*
count plateau (a period-2 oscillation would loop forever), its message over-claims
("too coarse") for an undecided regime, and non-unit coefficients are only peeled by 1
per visit.

## Feasibility remains undecided — do not read it off the probe curve

The probe signature is a greedy-cascade signature, **not** a feasibility bit. Grok's
prior exact-rule simulation saw transient coefficient growth to |c| = 7 on instances
that were *feasible* under SPFA, so the climb is compatible with Regime 1
(feasible-but-greedy-failed) as well as Regime 2 (genuinely infeasible throat between
the two layers). Only the difference-constraint SPFA verdict decides
(plan's standing 2026-07-01 caveat, reaffirmed three-way).

## Decisions

- **No greedy patches** (offender-list snapshot, one-fire-per-sweep, endpoint
  heuristics): unanimous "not worth it" — none yields a certificate, all compete with
  T1 effort. At most, snapshotting is a debug-only probe cleanup.
- **Execute T1+T2 of `todo/thin_cut_nonunit_rectification_implementation.md`** (shared
  `feasibility_solve()` SPFA core + read-only Phase-0 diagnostic) with the double-layer
  corc as reproducer, ending in the T3 feasible/infeasible verdict for generator 6.
  Implementation start awaits Christian's go.
- **Plan updated:** the Phase-0 "corc is unit-coefficient / does not qualify" caveat is
  obsolete — it described the single-layer corc (which stays the periodic regression
  case). The double-layer corc is the missing non-unit reproducer; the reproducer
  precondition is now satisfied. T2 and the greedy-termination open item updated
  accordingly.

## Files touched

- `todo/thin_cut_nonunit_rectification_implementation.md` — reproducer precondition
  ticked, Phase-0 caveat superseded, T2 retargeted, 2026-07-15 field observation added.
- `tmp/ai_exchange/corc_clean_divergence.md` — full audit thread (to be swept).
- No source edits.
