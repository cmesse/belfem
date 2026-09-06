# Devlog 2026-06-09 - CutSet Periodic Duplicate Audit

**Date:** 2026-06-09
**Topic:** Audit of periodic handling added to `CutSet::create_duplicates()`
**AIs involved:** Codex
**Codex Audit Confidence:** high

## Summary

Reviewed the user's draft changes to `CutSet::create_duplicates()` for Step 2 of the periodic thin-cut continuity plan.

## Findings

- The change is aimed at the right timing window, because `CutSet` duplicate maps are still available before `CutProcessor::collect_duplicates()` clears them.
- The implementation is not correct as written: a periodic pair present in `mNodeOriginals` is processed twice because only the current node is marked as processed, not its periodic partner.
- The new periodic duplicate pointers target the original partner nodes instead of the duplicate partner nodes.
- Both periodic duplicates are assigned the same original source node, so the second duplicate would later be linked to the wrong original by `CutFactory::link_node_duplicates_and_originals()`.
- The code uses node flag `1` as scratch, which conflicts with the planned need to preserve/copy periodic ownership flags `1/2`.
- The code assumes the periodic partner is present in `mNodeOriginals` and has a valid `mNodeBitset` index. That should be asserted explicitly before indexing the bitset.

## Files Updated

- `devlog/dl20260609_cutset_duplicate_periodic_audit.md`
- `devlog/README.md`
