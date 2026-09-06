# Devlog 2026-05-22 - Thick-Cut Multiplicity Audit

**Date:** 2026-05-22
**Topic:** Audit of proposed thick-to-thin cut handling for tetrahedral edge coefficients greater than one
**AIs involved:** Codex
**Codex Audit Confidence:** high
**Literature References:** Alves et al. 2022b Eq. 10; Alves et al. 2024 thin-cut jump condition; local homology docs

## Summary

Read-only audit of the hypothesis that tetrahedral thick-cut edge coefficients can be transformed into thin-cut face jumps by `k = (1/3) A e`, and that duplicate-node identity across adjacent elements can be decided by sharing two of the three cut-causing edges.

## Key Findings

- The local rows in the proposed matrix match BELFEM's existing TET4/TET10 edge and facet conventions for the four vertex-cap cut cases.
- The extraction formula `k = (1/3) A e` is not a valid inverse. The row matrix satisfies `A A^T = 4 I - J`, so the four cap multiplicities are only determined up to a common gauge. A correct recovery must integrate the edge cochain and choose a gauge, for example `k3 = 0`, `k0 = e5`, `k1 = e3`, `k2 = e4`, with consistency checks on `e0`, `e1`, and `e2`.
- Current BELFEM code does not support coefficients with absolute value greater than one in cut processing. `CutData::collect_coefficients()` only records coefficients exactly `+1` or `-1`; higher magnitudes remain in the edge support but have zero weight during cut-case validation.
- Duplicate identity should be based on the full branch label, i.e. the vector of cut jump coefficients after gauge selection, not only on whether two local cut triangles share two edges. The two-shared-edges rule is at most a simple-case adjacency hint.

## Changes Made / Proposed

- No source-code changes made.
- Proposed direction: replace direct averaging with an integer potential/gauge recovery per tetrahedron, and carry signed integer jump weights into duplicate-node source weights instead of using only bitsets.

## Open Questions

- Decide whether BELFEM should continue reducing all cohomology generators to `{-1,0,1}` before cut processing, or intentionally support integer multiplicities in `CutProcessor`.
- If multiplicities are supported, define a canonical gauge policy that minimizes duplicate creation and is stable across adjacent elements.

## Files Updated

- devlog/dl20260522_thick_cut_multiplicity_audit.md
- devlog/README.md
