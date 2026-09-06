# Devlog 2026-04-06 — normal_penta Follow-Up Review

**Date:** 2026-04-06
**Topic:** Follow-up review of prism normal conventions after inspecting `tmp/penta/facet.m`
**AIs involved:** Codex
**Claude Confidence:** N/A
**Codex Audit Confidence:** high
**Literature References:** N/A

## Summary

Reviewed the user's stated prism-normal convention and the Octave derivation in `tmp/penta/facet.m`. This resolves one earlier objection: equal TS normals on the two shell faces are consistent if `normal_penta_ts()` is meant to encode a common shell director rather than outward normals. Two issues remain: the branch still does not compile, and the implemented z-component for the top triangular face still disagrees with the user's own stated formula.

## Key Findings

- The new `GeometryType::PENTA` dispatch is now present in `Calculator::allocate()`.
- The TS convention in `tmp/penta/facet.m` supports `n0_ts = n1_ts` when TS facet `0` is treated as a flipped copy of volume face `3`.
- The code still fails to build because of `mFunInverJ` / `mFunInvertJ` and `mDeJ` / `mDetJ` naming mismatches.
- The z-component for volume case `4` and TS case `1` is still opposite to the user-stated formula for `n4`.

## Changes Made / Proposed

- No source changes made.
- Appended the follow-up audit to `todo/ai_exchange.md`.

## Open Questions

- Whether the generic documentation for `Calculator::normal()` should distinguish outward normals from shell-direction normals on TS prism sidesets.

## Files Updated

- todo/ai_exchange.md
- devlog/dl20260406_normal_penta_followup_review.md
