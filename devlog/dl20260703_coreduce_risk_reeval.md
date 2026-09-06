# Devlog 2026-07-03 — Coreduce Performance Findings: Literature-Grounded Risk Re-Evaluation + Low-Risk Fix Batch

**Date:** 2026-07-03
**Topic:** Re-evaluate the risk table in `todo/coreduce_complexPellikkaGeneralized_performance_findings.md` against new primary literature; apply the approved low-risk fixes
**AIs involved:** Claude, Codex
**Claude Confidence:** high on applied fixes; medium (~65%) on finding #11 reachability (partially refuted by Codex, accepted)
**Codex Audit Confidence:** high (C1–C7 confirmed with independent file:line evidence; C8 partially refuted)
**Literature References:** Mrozek & Batko 2009 (coreduction; §4 reduction pairs, §5 empty-cell trick, §6 Algorithm 6.1 / Thm 6.2), Pellikka et al. 2013 §2.2 — both added to `literature/papers/topology/` this session

## Summary

The user asked for a fresh risk pass over the 2026-06-22 performance findings before fixing
anything, since the Cohomology module was not written by him. Two papers (mrozek2009,
pellikka2013) were added to the literature library mid-session at Claude's request; they
materially changed several risk ratings. After a Codex audit of all claims
(`tmp/ai_exchange/coreduce_perf_risk_reeval.md`), the approved low-risk batch (#4, #7/#5
residual, #9) plus the audit-unblocked #3 were applied to source.

## Key Findings

- **#1 de-risked:** the proposed worklist IS Mrozek & Batko 2009 Algorithm 6.1 (queue of
  neighbors of the last removed pair, membership-checked); Thm 6.2 guarantees homology
  preservation for any valid removal order — order affects only reduction depth. The
  full-rescan `pCoreduce()` is the deviation, not the baseline.
- **#2 option (b) struck (unsafe):** Mrozek §5 empty-cell trick licenses omitting ONE 0-cell
  per connected component; batch omission leaves empty-boundary edges `pCoreduce`'s
  `size()==1` test can never remove → junk survives into cocombine → spurious generator risk.
- **#3 re-rated medium-low:** the "invariant" cleanup only runs when ≥1 neighbor passes the
  guard (naive hoist changes the zero-neighbor case), and at `tID2 == a` line 1173 erases `b`
  from the very map the neighbor loop iterates. Codex confirmed the live valid-neighbor set
  equals a pre-scan snapshot (loop body never inserts into `a`'s boundary).
- **#6 struck (premise wrong):** `mCochainsMap` is `Cell<Map<...>>`; `Cell::operator()` is
  bounds-checked array indexing, free in release. All three AIs missed this on 2026-06-22.
- **NEW #11 (latent correctness gap):** Mrozek §4 requires the pair coefficient κ(b,a) to be
  invertible (±1 over ℤ); `pCoreduce()` checks only boundary size == 1. Codex refuted the
  claimed reachability via `pGeneralizedCocombine(p)` → `pCoreduce(p+1)` (cocombine writes
  (p+1)-boundaries; that pCoreduce reads (p+2)-boundaries). Most reachable inside the
  non-generalized `pCocombine()` (used only by `coreduce_complexPellikka()`, not the
  generalized production path). Defensive `abs(coeff)==1` hardening left as a user decision.

## Changes Made

All in `src/homology/`, user-approved:

- `cl_SimplicialComplex.cpp` — #4: scan condition and `val1` now use `it2->second` instead of
  two redundant `getCoefficient()` re-searches of the same map. #7/#5: `mCochainsMap(p+1)(a)`
  hoisted into `tCochainA`, serving both the neighbor setup and the per-neighbor (p+2)
  cleanup. #3: `pGeneralizedCocombine()` restructured into snapshot (valid neighbors →
  `tNeighborIDs`/`tNeighborCochains`, function-scope buffers) → one conditional invariant
  cleanup pass → neighbor-dependent accumulation pass; exact pre-split semantics preserved
  (cleanup skipped when no valid neighbor, identical coefficient reads, commuting erases/adds).
- `cl_SimplicialComplex.hpp` — #9: `remove_kcochainFromMap()` now single `find()` + delete +
  erase-by-iterator (behavior-identical even for absent keys; the old `operator[]` path
  default-inserted a nullptr and erased it again).

`todo/coreduce_complexPellikkaGeneralized_performance_findings.md` gained a
"Risk Re-Evaluation (2026-07-03)" section, updated priority table, and a fix-campaign
checklist; finding #2(b) and #6 struck in place.

## Open Questions

- #11 **deferred by Christian same session** after cross-reading
  `todo/thin_cut_nonunit_rectification_implementation.md`: non-unit *generator* coefficients
  are a confirmed downstream reality handled by the planned rectify-or-certify solver (a
  distinct phenomenon from invalid non-unit *pair removal*, which #11 would guard). Adding
  the unit check now would change which cells reduce and shift the generator representatives
  under that plan mid-flight. Revisit once rectification T5 guards + T8/T9 tests exist.
- #1/#2(a): port `pCoreduce()`/`coreduceOmit()` to Algorithm 6.1 worklists — highest
  expected speedup; same sequencing argument: land the rectification plan first, then its
  guards/tests double as the regression instrumentation for the worklist port (plus Betti
  number and cut-count comparisons).
- #10 (stale adjacency cleanup) and the `addCochainToCochain()` structural merge remain open.
- Build + cohomology regression run handed off to Christian (shared build tree).

## Files Updated

- src/homology/cl_SimplicialComplex.cpp
- src/homology/cl_SimplicialComplex.hpp
- todo/coreduce_complexPellikkaGeneralized_performance_findings.md
- devlog/dl20260703_coreduce_risk_reeval.md (this file)
