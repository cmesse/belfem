# Periodic facet-map crash — root cause (unpaired one-sided-jump seam duplicates)

**Date:** 2026-06-17
**Purpose:** Record the auditor-confirmed root cause of the Maxwell-pass
`create_facet_map` "Node N is not flagged" crash, and the fix options it implies.
**Module:** src/mesh (periodicity) ↔ src/homology (cut duplication)
**Status:** Diagnosis only — no source behavior changed. Confirmed by Codex + Grok (high confidence).

---

## Where we are

Christian parked the facet-based `match_nodes_and_edges` refactor (commented,
`cl_Mesh_PeriodicityFactory.cpp:372-411`), reverted `update_periodicity` to the
original **backup-restore + old `match_edges`** path, and added a read-only
`check_health()` DEBUG detector. With the periodic block enabled, the run now clears
cohomology and the thick-cut relink, and fails in the **Maxwell-pass** periodicity
update (`MaxwellFactory.cpp:479` → `Periodicity::update`):

```
what():  Node 19145 is not flagged   (create_facet_map, cl_Mesh_PeriodicityFactory.cpp:890)
```

`check_health` **passes** — the seam facets are correctly relinked to their masters.

## Root cause (confirmed by both auditors)

The `#DIAG 5j` probe (added to `create_facet_map`, dumps every facet node missing from
the collected `aNodes`) shows, in the Maxwell pass:

- **565 unflagged facet-node incidences = 178 distinct nodes** (pre-cut passes: 0).
- Every one: `is_dup=1` (a cut duplicate), `periodic=0` (no partner),
  index-consistent + registered in `mMesh->nodes()`, on **sideset 7 / conductor
  block 2** (thin-shell seam).
- 178 distinct sits between the 6a/6b counts (167 one-sided-jump + 35 fragmentation =
  202) → these **are** the one-sided-jump / unbacked cut duplicates.

The chain (Codex + Grok CONFIRMED Q1/Q2/Q3, REFUTED Q4):

1. The 178 duplicates are placed on the periodic master facets **upstream**, by the
   cut pipeline (`CutFactory::relink_slave_elements_with_duplicate_nodes`), *before*
   `update_periodicity` runs. `collect_nodes` (`:357`) does gather them.
2. But the active path takes `restore_node_pairs` (`:415`), which **overwrites**
   `mMasterNodes`/`mSlaveNodes` from the **paired-node backup only**
   (`cl_Mesh_Periodicity.cpp:267-289`) — discarding the collected facet nodes.
3. The backup contains only **paired** nodes. 6b's exactly-one-bit branch
   (`cl_CutSet.cpp:144-155`) creates a duplicate with **no `set_periodic` and no
   `add_node_pair_to_backup`**. (6c, which would register the fragmentation/extra
   dups, is **not implemented** — `mPairVerdict` is built but never consumed.)
4. So the 178 unpaired seam duplicates are absent from `master_nodes`, and
   `create_facet_map`'s invariant (every facet node must be in `aNodes`,
   `:836-890`) correctly fails.

**Corrections to the in-progress framing** (both auditors): it is a **list overwrite**
(`restore_node_pairs`), *not* a facet relink between `:357` and `:467`; nothing in that
window relinks facet nodes.

## What this means

This crash is the **facet-level manifestation of the long-deferred one-sided-jump
policy** (6b/6c). The auditors agree the `create_facet_map` invariant is **correct** —
the fix is **not** to weaken the check. The 178 dups are genuinely unpaired (they sit
on opposite sides of a period-wrapping cut), so they cannot enter a pair-ordered node
list cleanly.

## Fix options (Christian's design call)

- **(A)** Exclude periodic seam facets that carry unpaired period-wrapping-cut
  duplicates from the periodic facet/face matching — they have no slave partner anyway.
- **(B)** A cut-aware pairing space that admits unpaired seam duplicates.
- **(C)** Implement 6c so the one-sided-jump / fragmentation duplicates are handled
  consistently end-to-end, and decide whether such a duplicate belongs in
  `master_nodes` at all.

Recommend discussing **A vs C** — both tie back to the jump policy we kept deferring;
this crash finally pins exactly where it has to be decided.

## Probes left in place (removal: plan step 8i)

- `create_facet_map`: `#DIAG 5j` (unflagged-facet-node registry/dup/periodic survey).
- `match_edges`: `#DIAG 5f/5g/5i` (now on the dead `match_edges`, since the active path
  uses `match_edges` at `:454` — note: `5f/5g/5i` are in the *old* `match_edges` body
  which IS still called on the active path; verify before removal).
- `collect_nodes`/`check_health`: Christian's relinking checks (`#ifdef DEBUG`).

No source behavior was changed tonight. The parked facet-based refactor (`:372-411`)
and its fixes (1–4 + collect_edges index/id) remain commented for when the jump policy
is settled.
