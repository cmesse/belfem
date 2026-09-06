# Helix Periodic Example: Cross-Block InterfaceSet Pairing Bug

**Date:** 2026-07-10
**Purpose:** Root-cause analysis of the `hphirun` abort on the helix translational-periodic example; three-AI audit (Claude primary, Codex + Grok auditors). Fix implemented after Christian's approval (see "Applied Fix" below).
**Module:** homology (InterfaceProcessor), mesh (Periodicity)

## Symptom

`cmake-build-debug/helix` (`helix_translate_periodic.msh`, four helical
conductors, 90° twist, pure-translation periodicity along z) aborts during
`MaxwellFactory::create_cuts`:

```
cl_InterfaceProcessor.cpp:206  InterfaceSet::duplicate_nodes
"periodic partner 223 of interface node 328 was not duplicated"
```

Mesh generation and periodic node-pair detection work; the failure is in
interface node duplication.

## Root Cause (confirmed by Codex and Grok independently)

`InterfaceProcessor` creates one `InterfaceSet` per block-pair bitset
(`cl_InterfaceProcessor.cpp:325-397`): {cond0, air}, {cond1, air}, … are
distinct sets with private `mOriginals`/`mDuplicates` maps. The periodic
re-pairing inside `InterfaceSet::duplicate_nodes` (`:183-226`) looks up
`A->periodic()` only in **its own** `mDuplicates` — assuming a periodic
partner of an interface node always belongs to the same interface set.

The helix geometry breaks that assumption by design: with 90° twist +
4-fold symmetry, `Periodic Surface` maps the front face of conductor *k*
to the back face of conductor *k−1*. Verified from the mesh:

- Node 223: curve 9, (0.278, 1.416, 0) — conductor-1/air interface
  (surface 51 bounds volumes {2, 5}).
- Node 328: curve 25, (0.278, 1.416, 5) — conductor-0/air interface
  (surface 41 bounds volumes {1, 5}).

So partner 223 is (or would be) duplicated in the **sibling** set
{block2, air5}, invisible to set {block1, air5}. Every earlier periodic
test (corc etc.) had block-preserving periodicity, so the per-set
assumption never fired. Timing (Grok): either set order fails — if the
{1,5} set runs first, the partner is not duplicated anywhere yet; the fix
must re-pair only after all sets created their duplicates. Release builds
abort too: with `BELFEM_ASSERT` compiled out, `Map::operator()` raises
`BELFEM_ERROR` on the missing key.

Thin cuts are NOT implicated in this abort: `CutSet::create_duplicates`
(`cl_CutSet.cpp:121-155`) already branches both-bits / one-sided per set,
and `CutProcessor::classify_periodic_pairs` (`cl_CutProcessor.cpp:704-765`)
is cross-set aware.

## Applied Fix (Christian-approved decisions: implement now / hard
## BELFEM_ERROR / keep Decouple pairing)

Minimal-diff two-phase design, per audit constraints:

1. **Same-set pairing stays in `InterfaceSet::duplicate_nodes`** — the
   hard assert on the set-local map was replaced by a skip
   (`if ( ! mDuplicates.key_exists( B->id() ) ) continue;`), deferring
   block-crossing pairs. Preserves today's semantics exactly for
   block-preserving periodicity, including multi-set (triple-junction)
   originals, which each set still pairs with its own copies.
2. **New `InterfaceProcessor::pair_cross_set_periodic_duplicates()`**,
   called from `InterfaceProcessor::duplicate_nodes` after all sets
   created their duplicates (guarded by `mMesh->has_periodicity()`).
   For each set's periodic original `A` whose duplicate `C` is still
   untied (`!C->is_periodic()` — the natural "already paired" gate, also
   preventing double-processing from the partner's side), it searches the
   sibling sets for the partner's duplicate `D`, ties `C↔D`, and adds the
   pair to the periodicity backup once.
   Guards (all `BELFEM_ERROR`, always-on per decision): partner
   duplicated in no set; partner duplicated in more than one sibling set
   (ambiguous); partner's duplicate already tied to another node.
3. Decouple (coil) sets keep participating in periodic pairing —
   behavior-preserving; revisit if a periodic coil case misbehaves
   (Codex's `set_entity_dependencies` re-sourcing concern,
   `cl_Mesh_Periodicity.cpp:53`, stays on record).

Files: `src/homology/cl_InterfaceProcessor.cpp` (skip + new method),
`src/homology/cl_InterfaceProcessor.hpp` (declaration).
Syntax-verified with the real homology-target build flags
(`g++ -fsyntax-only` + flags.make); build and helix run handed to
Christian.

## Separate Debt Found (not this abort)

Step 6c is unimplemented: `mPairVerdict` from
`CutProcessor::classify_periodic_pairs` is write-only
(`cl_CutProcessor.cpp:764`; zero readers in `src/`). The comment at
`cl_CutSet.hpp:26-27` claiming it is "consumed by create_duplicates()"
over-claims, and FRAGMENTATION extra duplicates remain unregistered in the
periodic backup (`cl_CutSet.cpp:106-118` caveat). Also, the "Periodicity of
duplicates | handled" row in `todo/interface_node_duplication_coil_ferro.md`
should be demoted to cover cross-InterfaceSet partners.

## Artifacts

- Exchange thread: `tmp/ai_exchange/helix_periodic_interface_sets.md`
  (Claude claim + Grok audit + Codex audit + resolution)
- Example: `cmake-build-debug/helix/` (`helix_translate_periodic.geo/.msh`,
  `input.conf`)

## Correction (same day): "Step 6c debt" framing was wrong

Christian challenged the 6c claim after the helix ran. Re-reading
`todo/closed/periodic_thin_cut_continuity_fix.md` settles it:

- Step 6 (a/b/c/f) was **SUPERSEDED** by Step 9 (plan §Step 6, items
  6a-6c/6f). The corc `create_facet_map` "Node not flagged" blocker was an
  induced-periodic **keying** problem, not a duplication-policy problem.
- Step 9 (original-identity keying) **is implemented**: the periodic
  rebuild keys every facet/edge corner by `node->original()->index()`
  (`cl_Mesh_PeriodicityFactory.cpp:563, 729-730, 801-817, 854-860`);
  functional blocker cleared 2026-06-18 (corc full pipeline).
- Pairing the induced one-bit duplicates was explicitly **rejected as
  unsafe** (plan §"Why Pairing Was Rejected"): may have no image dup to
  pair with; would clobber the original's single `mPeriodic`;
  `set_entity_dependencies()` wipes slave sources from the master →
  T-matrix corruption. Implementing 6c could break the working path.
- The real debt is **cleanup, not function** (plan item 8j): remove or
  document `classify_periodic_pairs()`/`mPairVerdict` (write-only WIP),
  and fix two stale comments that misled all three AIs in this session's
  audit — `cl_CutSet.hpp:26-27` ("consumed by create_duplicates()", false)
  and `cl_CutSet.cpp:106-119` ("Step 6c ... must register them or the
  periodic rebuild drops them", pre-Step-9 framing).
- Still-live Step 6 residuals per the plan: 6e (promote relink/symmetry
  asserts to named errors), 9e (genuine geometric jumps with non-periodic
  originals), 5d (always-active coverage checks). 6d (release-mode
  asymmetry) appears resolved by the unconditional three-way branch
  ("correct in BOTH builds", cl_CutSet.cpp:103-104).
