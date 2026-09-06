# DR-109 / `mPairVerdict`: bisect, provenance, and the deletion that never landed

**Date:** 2026-08-30
**Purpose:** Answer whether `CutProcessor::mPairVerdict` is safe to delete; establish when and why it was added.
**Module:** `src/homology`

---

## The question

Christian: `mPairVerdict` appears to have no consumer. Bisect when it was added,
read the devlogs for its intended purpose, and decide whether it can go.

## Provenance (git, not devlog)

`git log --all -S mPairVerdict -- src/homology/` returns exactly two commits.

- **`6347dd7b`** "saving work in progress", 2026-06-16 — the birth. One commit
  introduces `mPairVerdict`, `classify_periodic_pairs()`, `CutPairVerdict`,
  `CutPairInfo` and `CutSet::cut_pattern()`. The same commit rewrites
  `CutSet::create_duplicates` from a debug-only `#DIAG 4c` survey into the live
  three-way both-bit / one-bit / neither branch, and writes the comment block
  that promises a "Step 6c ( consuming mPairVerdict )".
- **`69fe2a58`** "updating AI collaboration rules", 2026-06-18 — removes the
  `#DIAG 4h` probes. The pickaxe matches only because a header comment changed
  wording. The data flow is untouched.

**The table was therefore born without a consumer, and never acquired one.**
Step 6c was described in a comment in the same commit that created the producer,
and was never written. Two days later `69fe2a58` landed on the same branch, and
Step 9 (original-identity keying in `cl_Mesh_PeriodicityFactory`) became the
mechanism that actually solved the problem. `todo/closed/periodic_thin_cut_continuity_fix.md:111-129`
records that pairing the induced duplicates was **rejected as unsafe**, not
deferred: it would overwrite a node's single `mPeriodic` slot, and
`set_entity_dependencies()` would then corrupt slave T-matrix semantics.

## What the devlogs show: the same trap, four times

The comment has misled every audit that touched it.

| Date | Record | What happened |
|---|---|---|
| 2026-06-17 | `dl20260617_periodic_facetmap_rootcause.md:49` | Noted 6c "not implemented"; listed implementing it as fix option C |
| 2026-07-10 | `dl20260710_helix_periodic_interface_sets.md:88` | A three-AI round concluded 6c was debt that must be closed — **wrong**; INC-030 records that two stale source comments misled all three AIs |
| 2026-08-13 | `dl20260813_debt_register_currentness_sweep.md:119` | Zero reads re-confirmed; exact deletion set pinned |
| 2026-08-26 | `dl20260826_dr108_uint8_source_overflow.md:138` | DR-109 opened, asserting `create_duplicates` **consumes** the map — the trap again |

## Two register rows in direct conflict

- `todo/debt_register.md` DR-109 (open, P2) says "The Step 6 duplication branch
  in `CutSet::create_duplicates` consumes that verdict map", and rates the
  64-cut silent cap at `cl_CutProcessor.cpp:734` on that consumption.
- `todo/debt_register_closed.md` DR-26 says the whole machinery is dead and was
  **deleted on 2026-08-24**.

Both are wrong in different ways, and the code settles it.

## Audit round

Codex `gpt-5.6-terra`/`xhigh` and Grok `grok-4.6`/`xhigh`, blind and parallel,
from a brief that presented the DR-26/DR-109 conflict without naming a preferred
side. Claims pre-registered before dispatch in
`tmp/ai_exchange/dr109_pairverdict_deletion.md`. Read-only respected by both.

- **Zero reads.** 1 declaration (`cl_CutProcessor.hpp:78`), 2 writes
  (`clear()` `:731`, `operator[]=` `:783`), 0 reads. Both auditors independently
  produced the same counts and both checked the indirect paths — accessors,
  friends, `#if !defined(NDEBUG)`, `save_debug_meshes`, serialization, tests,
  `nonfree/`, and CMake-conditional sources. `cl_CutProcessor.cpp` and
  `cl_CutSet.cpp` compile unconditionally.
- **`classify_periodic_pairs()` is mesh-state pure.** It reads node
  periodicity/id/index and `DynamicBitset::test` (const), sets no flag, mutates
  no entity, changes no index or ordering. Its only write is the map. Removing
  the call at `:82` cannot change solver behaviour.
- **DR-109's mechanism is refuted.** `CutSet::create_duplicates` takes no
  processor and no verdict argument (`cl_CutSet.hpp:85-86`); its branch reads
  `mNodeBitset` and nothing else (`cl_CutSet.cpp:125-159`). The claimed
  consumption exists only as prose at `:107-123`.
- **The width concern does not migrate.** Both auditors swept `src/homology/`
  for a surviving cut-count mask. The only `unsigned long long` / `1ull` sites
  are `CutPairInfo::mUnionBits` (`cl_CutSet.hpp:40`) and the four locals at
  `cl_CutProcessor.cpp:750-768` — all inside the deletion set. The live
  cut-membership representation is a runtime-sized `DynamicBitset`, uncapped.

## The finding that changes what to do next

**DR-26's deletion was never lost — it is in `stash@{0}`.**

`git log --all -S` finds no commit removing the symbol, `git status
src/homology/` is clean, and the code is live at HEAD. But
`dl20260824_night_shift_dr02_dr102_dr103.md:243` listed DR-26's deletion as
uncommitted working-tree state, and `stash@{0}` ("bisect: session homology/mesh
changes (DR-26/102/103/13)") carries it in full: the member, the function
declaration + body + call site, `CutPairVerdict`/`CutPairInfo`,
`CutSet::cut_pattern()`, and the corrected `cl_CutSet.cpp` comment that cites
the closed plan instead of promising 6c.

It is **not** a clean cherry-pick. The same stash carries DR-103's
`check_midside` signature change (adds `Element *`, adds an `is_flagged` guard)
and edits to `src/mesh/cl_FaceFactory.cpp` and `cl_Vertex.cpp`. Extracting DR-26
alone means taking four hunks across `cl_CutProcessor.hpp/.cpp` and
`cl_CutSet.hpp/.cpp` and leaving the rest.

## Answer

**Yes, safe to delete** — *reviewed*, not verified. Three independent static
traces, no build and no periodic run against the deletion.

## Owed with the deletion, not optional

1. **Rewrite `cl_CutSet.cpp:107-123`.** Deleting the code and leaving the comment
   saying "Step 6c must register them" recreates INC-030 a fourth time. The
   stashed text already does this correctly.
2. **Strike `todo/test_hardening_campaign.md` step R3d** (`:291-294`) — an open,
   unticked test contract aimed at `classify_periodic_pairs` above 64 cuts.
3. **Close DR-109** as dissolved, not fixed. Its named mechanism is wrong; there
   is nothing left at `:743-746` once the function goes.
4. **Correct the DR-26 row** in `todo/debt_register_closed.md`. It records a
   landed deletion that never reached git. Left as is, the next sweep re-derives
   the same false closure.

## Adjacent, out of scope

`CutSet` allocates `mNodeBitset` at `cl_CutSet.cpp:31`; `~CutSet()` (`:35-38`)
deletes only `mBitset`. A real leak, confirmed by reading, independent of this
deletion. Do not fold it into the cleanup.

## Standing state

No source edited. Exchange thread `tmp/ai_exchange/dr109_pairverdict_deletion.md`.

---

## Execution (2026-08-31, Fable, Christian's "let's do that")

The deletion landed by hand — the stash's `src/homology/` diff was used as the
reference, **not** applied: `stash@{0}` also carries DR-103's `check_midside`
change and `src/mesh/` edits, all of which remain stashed and unclaimed.

Applied, exactly the DR-26 set:

- `cl_CutProcessor.cpp` — the `classify_periodic_pairs()` call (`:82`) and the
  full function body removed
- `cl_CutProcessor.hpp` — `mPairVerdict` member + comment, declaration block
  removed
- `cl_CutSet.hpp` — `CutPairVerdict` enum, `CutPairInfo` struct, `cut_pattern()`
  declaration + inline definition removed
- `cl_CutSet.cpp` — the Step-6c comment replaced with the stash's corrected
  account (rejection recorded in `todo/closed/periodic_thin_cut_continuity_fix.md`;
  the periodic rebuild keys by original identity in
  `cl_Mesh_PeriodicityFactory.cpp`)

All four owed items done: comment rewritten (item 1), R3d struck with the
campaign Status line updated (item 2), DR-109 struck as **dissolved** with the
`[W]` header recount 6→5 (item 3), and the DR-26 closed row carries a
**CORRECTION 2026-08-31** paragraph recording the false closure and the hand
re-landing (item 4).

Gates: tree-wide sweep finds zero remaining references to any deleted symbol in
`src/`/`tests/`/`nonfree/`; both touched TUs pass `g++ -std=gnu++17
-fsyntax-only` with the homology module's production flags (`-Wall -Werror`);
`scripts/check_doc_claims.py` 37/37 after the recount. Still **reviewed, not
verified** — no build ran, `make check` owed with the next build Christian runs.

The `mNodeBitset` leak stays out, as ruled above.
