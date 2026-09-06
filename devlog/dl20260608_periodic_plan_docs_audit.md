# Periodic Plan Docs Audit

Date: 2026-06-08

## Summary

Read-only audit of the updated periodic boundary condition planning docs:

- `todo/periodic_bc_fix_plan.md`
- `todo/periodic_thin_cut_continuity_fix.md`
- `todo/periodic_thin_cut_strategy.md`

The focused thin-cut continuity plan is directionally correct, but the broader fix plan still contains stale API/function references from an older periodicity design. The continuity plan also needs two implementation corrections before it is used as a coding checklist.

## Findings

- `periodic_bc_fix_plan.md` references several mechanisms that do not exist in the current tree: `create_dofwise_periodicities_master()`, `flag_periodic_entities_12()`, `collect_nodes_from_flags_12()`, and `Periodicity::match_faces()`. The current live rebuild is `Periodicity::update()` plus `Periodicity::set_entity_dependencies()`, and periodic hanging relationships are consumed through the ordinary `create_dofwise_t_matrices_master()` path.
- `periodic_thin_cut_continuity_fix.md` correctly identifies the thin-cut periodic continuity defect, but its proposed Step 1 must split `CutProcessor::mBoundaries` semantics. Periodic sidesets should be excluded from face/edge trimming, while still being available to mark boundary edges used by the peel-loop guard.
- The proposed Step 2 location is too late if it expects to use `CutSet::mNodeDuplicates`. `CutProcessor::collect_duplicates()` clears each `CutSet` duplicate map before `CutFactory::link_node_duplicates_and_originals()` runs.
- The old `periodic_thin_cut_strategy.md` is mostly superseded as an active implementation strategy, but it still contains useful material worth scavenging: the physical-lift acceptance contract, the post-clean periodic-edge diagnostic, and the empty/stale sideset-node guard idea for `PeriodicityFactory::select_sidesets()`.

## Recommendation

Keep `periodic_thin_cut_continuity_fix.md` as the active thin-cut plan after correcting the two implementation traps above. Clean up or archive `periodic_bc_fix_plan.md` sections that describe the removed flag-based and DofData-periodicity paths. Delete or archive `periodic_thin_cut_strategy.md` only after scavenging the useful diagnostics and acceptance criteria into the active plan.

## Tests

No code was modified and no tests were run.
