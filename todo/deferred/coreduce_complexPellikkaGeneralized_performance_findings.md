# Performance Findings: `SimplicialComplex::coreduce_complexPellikkaGeneralized()`

**Date:** 2026-03-11
**Refreshed:** 2026-06-22 — re-audited against the current tree (Claude + Codex + Grok, 3-way agreement, all citations checked)
**Risk re-evaluation:** 2026-07-03 — literature-grounded risk pass (Mrozek & Batko 2009 and Pellikka et al. 2013 added to `literature/papers/topology/`); see the section below
**Purpose:** Merged performance analysis, evaluated against actual code
**Module:** homology

> **Status: DEFERRED 2026-08-11** (currentness sweep) — moved to `deferred/`, where the file's
> own text has been pointing since 2026-07-03. The move records a state, it does not create one:
> the remaining item (#1/#2a, the worklist port) was explicitly sequenced behind
> `todo/thin_cut_nonunit_rectification_implementation.md`'s T5 guards and T8/T9 tests, which
> supply the Betti/cut-count instrumentation the port needs to be validated safely — and T8 is
> only part done. Nothing is waiting on this file and nobody has owned it since 2026-07-03.
> **Revive when** the rectification tests land, or when a profile puts `coreduce` back on the
> critical path. Every finding below stands; re-locate by symbol, not by line.
>
> **Status (2026-08-09 re-check):** the findings still stand — nothing in
> `cl_SimplicialComplex` has been restructured since — but **the 2026-06-22 line numbers
> have drifted too** and the refresh table below is no longer literal. Current anchors:
> `pCoreduce()` is `cl_SimplicialComplex.hpp:332` (decl `:129`), `coreduceOmit()` is
> `cl_SimplicialComplex.cpp:1236` (callers `:1269`, `:1297`). Re-locate the neighbor-loop
> findings (#3–#7) by symbol inside those two bodies rather than by the ranges quoted below.
> This file is a **cold** performance ledger: it has had no owner since 2026-07-03, is
> tracked as `debt_register.md` DR-29 (deferred, sequenced after the rectification
> guards/tests), and is not on the 1.0 path.
>
> **Status (2026-06-22):** The code changed in commit `cec4fb7` (2026-05-13, "Fixing the
> clean function again"). Every line number in the original 2026-03-11 analysis is stale;
> the structural bottlenecks mostly remain. The table below is the authoritative current
> state — independently re-audited by Claude, Codex, and Grok with full agreement. The
> per-finding sections that follow carry corrected line numbers and a verdict tag.

## Status Refresh (2026-06-22)

| # | Finding | Verdict | Current location |
|---|---------|---------|------------------|
| 1 | `pCoreduce()` full rescan | **STILL-VALID** | `cl_SimplicialComplex.hpp:339-381` |
| 2 | `coreduceOmit()` one-at-a-time | **STILL-VALID** | `cl_SimplicialComplex.cpp:1214-1226` |
| 3 | invariant cleanup in neighbor loop | **STILL-VALID** | `cl_SimplicialComplex.cpp:1173, 1182-1188` |
| 4 | redundant `getCoefficient()` | **STILL-VALID** | `cl_SimplicialComplex.cpp:1135-1136, 1153` |
| 5 | `key_exists` + `operator()` double lookup | **PARTLY RESOLVED** | neighbor loop now `find()` (`1158-1172`); residual at `1152/1185/1187` |
| 6 | hoist `mCochainsMap(p/p+1/p+2)` refs | **STILL-VALID** | `cl_SimplicialComplex.cpp:1119-1187` |
| 7 | hoist `a`'s coboundary out of neighbor loop | **STILL-VALID** | `cl_SimplicialComplex.cpp:1184-1185` |
| 8 | `erase_key()` double lookup | **RESOLVED** | `cl_OrderedMap.hpp:147`, `cl_Map.hpp:204` |
| 9 | `remove_kcochainFromMap()` triple lookup | **CHANGED → double** | `cl_SimplicialComplex.hpp:318-319` |
| 10 | stale adjacency entries inflate search | **STILL-VALID** | `cl_SimplicialComplex.cpp:1135-1143` |
| — | `addCochainToCochain()` recursive cost | **STILL-VALID** | `cl_Cochain.hpp:281-307` (call `cpp:1165`) |

**What changed since 2026-03-11:**

- **#8 RESOLVED.** Both `erase_key()` implementations now call `mMap.erase(aKey)` directly,
  with no `key_exists()` guard (`cl_OrderedMap.hpp:150`, `cl_Map.hpp:207`). The old code
  snippet quoted in the original finding is now factually wrong.
- **#9 downgraded to double access.** Because `erase_key()` was simplified (see #8),
  `remove_kcochainFromMap()` is now `operator[]` (find) + `erase_key()` (erase) — a *double*
  access, not the old triple lookup. Concern reduced, not eliminated.
- **#5 partly resolved.** The neighbor loop's `key_exists` + `operator()` pattern was
  replaced with a single `find()` (`cpp:1158-1161, 1170-1172`). What remains is *unguarded*
  `operator()` membership access — `mCochainsMap(p+1)(a)` at `cpp:1152/1185` and
  `mCochainsMap(p+2)(tID2)` at `cpp:1187`. That is a different concern (it assumes
  membership and would trip a debug `BELFEM_ASSERT` on a stale ID, see #10), not the
  original double-lookup.
- All other findings hold unchanged in substance; only line numbers moved.

---

## Risk Re-Evaluation (2026-07-03)

Fresh risk pass requested by the user before fixing anything, since the Cohomology module
was not written by the current maintainer and semantic side effects are hard to judge from code alone.
The pass is now grounded in the primary literature: **Mrozek & Batko 2009** (coreduction algorithm,
`literature/papers/topology/mrozek2009.txt`) and **Pellikka et al. 2013** (Gmsh solver, §2.2,
`literature/papers/topology/pellikka2013.txt`). Codex audit thread:
`tmp/ai_exchange/coreduce_perf_risk_reeval.md`.

| # | 2026-06-22 risk | Re-evaluated risk | Reason |
|---|-----------------|-------------------|--------|
| 1 | Medium | Medium, **de-risked in principle** | The worklist matches Mrozek & Batko 2009 Algorithm 6.1 (§6). Each valid elementary coreduction preserves homology (Thm 6.2), so with dequeue-time validity checks, processing order affects reduction depth rather than correctness. Residual risk is implementation care (stale queue entries) and cut-shape reproducibility, not mathematics. `pCocombine` (cpp:1023) already uses a `Q` worklist in-module. |
| 2 | Medium | Option (a): medium. Option (b): **unsafe, struck**. Option (c): medium | Mrozek & Batko 2009 §5: the empty-cell trick licenses omitting ONE 0-cell per connected component; the rest must be *paired away* by coreduction. Batch omission leaves edges with empty boundaries that `pCoreduce` (size==1 test) can never remove → junk survives into cocombine → possible spurious generators. |
| 3 | Low | **Medium-low** | The hoist is not purely mechanical: (i) the cleanup currently runs only when ≥1 neighbor passes the guard at cpp:1159 — a naive hoist runs it unconditionally; (ii) at `tID2 == a`, cpp:1173 erases `b` from `a`'s boundary map, which is the very map the neighbor loop iterates (safe today only by std::map erase semantics). A split must snapshot valid neighbors and keep the cleanup conditional. |
| 4 | Low | **Low, confirmed** | Pure no-op: `it2` iterates the same map that `getCoefficient(it2->first)` re-searches; `val1` equals `it2->second`. |
| 5 res. | Low | Not a perf item | The unguarded `operator()` at cpp:1187 is an invariant assert — a feature, keep. The `a` lookups at cpp:1152/1185 fold into #7. |
| 6 | Low | **Premise wrong, impact ≈ nil** | `mCochainsMap` is `Cell<Map<...>>` (hpp:36); `Cell::operator()` is bounds-checked array indexing, free in release. The 2026-06-22 audit mistook it for a hash lookup. Cosmetic only. |
| 7 | Low | **Low, confirmed** | `a` is guaranteed present (while-scan at cpp:1135) and not deleted until cpp:1193; hoisting the `Cochain*` is safe. |
| 9 | Low | **Low, confirmed** | find + delete + erase-by-iterator is equivalent even for absent keys (current `operator[]` default-inserts nullptr; `delete nullptr` is a no-op). |
| 10 | Medium | Medium | Unchanged; interacts with #3(i) — the zero-valid-neighbor skip may itself leave stale entries. |
| — | High (addCochainToCochain) | High | Unchanged, defer. |

**NEW finding #11 — `pCoreduce()` never checks the unit-coefficient condition.**
Mrozek & Batko 2009 §4: a (co)reduction pair requires κ(b,a) *invertible* in R (±1 over ℤ);
`pCoreduce` (hpp:353) tests only `size() == 1`. Safe inside `coreduceOmit` (coefficients stay
±1 there), but in the Generalized pipeline `pCoreduce(p+1)` runs *after*
`pGeneralizedCocombine(p)`, which accumulates coefficients (cpp:1165, 1177) and itself guards
its own pair choice with `abs(...) == 1` (cpp:1136). If a size-1 boundary with |coeff| ≥ 2
reaches `pCoreduce`, removing it violates Thm 4.1's premise → homology/torsion may change.
Confidence: medium (~65%) that the state is reachable; correctness gap if it is. The
non-generalized `pCocombine` (cpp:1052) applies no unit check on `val1`/`val2` either.
Possible defensive fix: add `abs(coeff) == 1` to the `pCoreduce` pair test (behavior-narrowing,
literature-aligned).

**Codex audit outcome (2026-07-03):** the local gap is real, but the specific
`pGeneralizedCocombine(p)` → `pCoreduce(p+1)` path is NOT reachable: cocombine writes its
accumulated coefficients into (p+1)-boundaries, while that `pCoreduce` call reads
(p+2)-boundaries, and `pCoreduce(p)` never runs again afterwards. The gap is most reachable
*inside* the non-generalized `pCocombine()` (size==2 test, unchecked `val1`/`val2`, worklist
feeds candidates back within the same pass) — but that variant is only used by
`coreduce_complexPellikka()`, not the generalized production path.

**Deferred (user decision, 2026-07-03) — sequenced behind the thin-cut rectification plan.**
Context from `todo/thin_cut_nonunit_rectification_implementation.md`: non-unit coefficients in
the *generators* are a confirmed downstream reality (coarse-CCT reproducer; transient growth to
|c| = 7 observed in greedy `clean()`), handled by the planned rectify-or-certify solver — a
distinct phenomenon from an invalid non-unit *pair removal*, which is what #11 would guard.
Adding the unit check now would change which cells reduce → change generator representatives →
shift the inputs to that plan mid-flight. Revisit after T5's release-active guards and T8/T9
tests land, which will provide the regression instrumentation to validate this change safely.
**Status: deferred — latent, not reachable in the generalized pipeline (Codex).**

### Fix campaign (2026-07-03)

- [x] Apply #4 — use `it2->second` instead of the two redundant `getCoefficient()` calls (cpp:1136, 1153) *(applied 2026-07-03)*
- [x] Apply #7 (+#5 residual) — hoist `mCochainsMap(p+1)(a)` into a local `Cochain*` (`tCochainA`) serving cpp:1152 and cpp:1185 *(applied 2026-07-03)*
- [x] Apply #9 — rewrite `remove_kcochainFromMap()` as find + delete + erase-by-iterator (hpp:314-321) *(applied 2026-07-03)*
- [x] Codex audit of the re-evaluation claims C1–C8 (`tmp/ai_exchange/coreduce_perf_risk_reeval.md`) — C1–C7 confirmed; C8 partially refuted (see #11) *(2026-07-03)*
- [x] Decide on #11 (unit-coefficient check) — **deferred by Christian (2026-07-03)**: Codex confirmed the local gap but ruled out the `pGeneralizedCocombine(p)` → `pCoreduce(p+1)` path (cocombine writes (p+1)-boundaries, while that `pCoreduce` call reads (p+2)-boundaries). Sequenced behind the thin-cut rectification plan (T5 guards + T8/T9 tests provide the safety net); see the #11 section
- [x] Apply #3 via snapshot-based split — audit confirmed the aliasing analysis; implemented with `tNeighborIDs`/`tNeighborCochains` snapshot Cells and cleanup conditional on ≥1 valid neighbor *(applied 2026-07-03)*
- [ ] ~~#2 option (b): batch-remove surviving 0-cochains~~ (unsafe per Mrozek & Batko 2009 §5)
- [ ] #1/#2(a): port `pCoreduce`/`coreduceOmit` to Algorithm 6.1 worklists — separate task, needs regression meshes (compare Betti numbers and cut counts before/after). Best sequenced AFTER the thin-cut rectification plan (`todo/thin_cut_nonunit_rectification_implementation.md`) lands: the port changes which cells survive reduction and thus the generator representatives, and T5's release-active guards + T8/T9 tests are exactly the instrumentation needed to validate that safely

---

## Call Path (current line numbers)

```
coreduce_complexPellikkaGeneralized()          cl_SimplicialComplex.cpp:1261
├── coreduceOmit()                             cl_SimplicialComplex.cpp:1206
│   ├── pCoreduce(0..2)                         cl_SimplicialComplex.hpp:327
│   └── while(0-cochains remain):              cl_SimplicialComplex.cpp:1214
│       ├── remove one 0-cochain                          (cpp:1222)
│       └── pCoreduce(0..2)                    ← full rescan after each removal
└── for p = 0..2:
    ├── pGeneralizedCocombine(p)               cl_SimplicialComplex.cpp:1119
    └── pCoreduce(p+1)                          cl_SimplicialComplex.hpp:327
```

## Data Structures in Play

| Container | Underlying type | Used for |
|-----------|----------------|----------|
| `mCochainsMap(k)` | `Map` = `std::unordered_map` | Top-level: simplex ID → Cochain* |
| `mSimplicesMap` (in Cochain) | `OrderedMap` = `std::map` | Sparse coefficients (+1/−1) per cochain |

---

## Algorithmic Bottlenecks (High Impact)

### 1. `pCoreduce()` repeated full rescans — replace with worklist/frontier — STILL-VALID

**Identified by:** Codex, Claude | **Location:** `cl_SimplicialComplex.hpp:339-381`

`pCoreduce()` uses a `while(tRemoved)` loop that rescans the *entire* (p+1)-cochain map
after each pass. Late in coreduction, most cochains won't qualify (boundary size ≠ 1),
so most scan work is wasted.

**Fix:** Maintain a queue of candidate (p+1)-cochains whose boundary size may have
dropped to 1. Seed once from the initial map. After removing pair (a,b), enqueue only
the directly affected neighbors instead of rescanning.

**Expected impact:** Largest single improvement — turns O(passes × |map|) into
O(|removals| × avg_neighbors).

### 2. `coreduceOmit()` quadratic behavior from repeated rescans — STILL-VALID

**Identified by:** Codex, Grok, Claude | **Location:** `cl_SimplicialComplex.cpp:1214-1226`

After the initial `pCoreduce(0..2)`, the while-loop removes ONE 0-cochain (line 1222), then
reruns all three `pCoreduce(p)` passes from scratch (lines 1223-1225). This multiplies the
cost of finding #1 by the number of surviving 0-cochains.

**Fix options (in order of preference):**
- a) Keep per-dimension worklists alive across the entire omit phase (combines with fix #1).
- b) ~~Remove all orphan 0-cochains in a batch, update boundaries, then run `pCoreduce` once.~~
  **Struck 2026-07-03 — unsafe.** Mrozek & Batko 2009 §5: only ONE 0-cell per connected
  component may be omitted (the empty-cell trick); the rest must be paired away by
  coreduction. Batch omission leaves empty-boundary edges `pCoreduce` can never remove.
- c) At minimum, only rerun the dimensions whose neighborhoods were actually affected
  (note: a removal cascade still propagates p=0 → 1 → 2, so the saving is smaller than it looks).

**Expected impact:** High — eliminates the multiplicative cost on top of fix #1.

### 3. Invariant cleanup repeated inside neighbor loop in `pGeneralizedCocombine()` — STILL-VALID

**Identified by:** Codex | **Location:** `cl_SimplicialComplex.cpp:1167-1189` (neighbor loop opens at 1156)

Inside the per-neighbor loop, two pieces of work do **not** depend on the current neighbor
`tID` yet are re-executed once per neighbor:

- the `setCoefficient(b,0)` clear at line 1173 (does not use `tID` or `val2`), and
- the entire `p < 2` block at lines 1182-1188, which clears `a` from each (p+2) boundary.

Only `add_simplex_to_boundary(tID, …)` at lines 1176-1177 (and the `val2`/`addCochainToCochain`
work at 1162-1165) is genuinely neighbor-dependent.

**Fix:** Snapshot the valid neighbors for (a,b) first. If the snapshot is non-empty, run the
invariant cleanup for line 1173 and the 1182-1188 block once, then run the neighbor-dependent
accumulation over the snapshot. Codex and Grok concurred on the invariant/neighbor-dependent
split; the 2026-07-03 caveats below constrain the exact implementation.

**Caveats (2026-07-03 re-evaluation — re-rated medium-low):**
- The cleanup currently runs only when at least one neighbor passes the guard at cpp:1159;
  a naive hoist runs it unconditionally, changing behavior in the zero-valid-neighbor case.
- At `tID2 == a`, line 1173 erases `b` from `a`'s boundary map — the very map the neighbor
  loop at cpp:1156 iterates. The split must snapshot the valid neighbors first and keep the
  cleanup conditional on that snapshot being non-empty to stay exactly behavior-preserving.

---

## Redundant Lookups (Medium Impact)

### 4. Redundant `getCoefficient()` when iterator already has the value — STILL-VALID

**Identified by:** Claude, Codex | **Location:** `cl_SimplicialComplex.cpp:1135-1136, 1153`

```cpp
// Line 1136: it2 iterates the coboundary's getSimplicesMap()
// but then re-searches the same map for the same key:
abs( tCochainReduce->getCoboundary()->getCoefficient( it2->first )) != 1
// Should be:
abs( it2->second ) != 1
```

Similarly at line 1153, `tCochainReduce->getCoboundary()->getCoefficient(a)` re-searches
the coboundary map for `a`, which is exactly the key the scan terminated on (`it2->first`).
The value is already known.

**Fix:** Use `it2->second` directly. Cache the coefficient of `a` from the search result.

**Impact:** Eliminates O(log n) `std::map::find` per inner iteration.

### 5. `key_exists()` + `operator()` double lookup pattern — PARTLY RESOLVED

**Identified by:** Codex, Grok, Claude | **Location (resolved):** `cl_SimplicialComplex.cpp:1158-1161, 1170-1172` | **Location (residual):** `1152, 1185, 1187`

The originally-cited neighbor-loop double lookup is **gone** — the code now does a single
`find()` and checks the iterator against `end()`:

```cpp
auto it3 = mCochainsMap( p ).find( tID );        // single lookup
if ( it3 != mCochainsMap( p ).end() and tID != b )
{
    Cochain * tCochainAdd = it3->second;
```

What remains is a *different* pattern: unguarded `operator()` membership access at
`mCochainsMap(p+1)(a)` (lines 1152, 1185) and `mCochainsMap(p+2)(tID2)` (line 1187). These
assume the key is present rather than doing `key_exists` first — correctness depends on the
no-stale-entry invariant (see #10), and in debug builds a stale ID would trip a
`BELFEM_ASSERT` inside `Map::operator()` (`cl_Map.hpp:227`) rather than degrade silently.

### 6. Hoist `mCochainsMap(k)` Cell references — premise wrong, impact ≈ nil

**Identified by:** Codex, Claude | **Location:** throughout `pGeneralizedCocombine()` (`cl_SimplicialComplex.cpp:1119-1187`)

`pCoreduce()` already caches a `mCochainsMap(p+1)` Cell entry (line `hpp:346`:
`Map<...>& tMap1 = mCochainsMap(p+1)`). Doing the same in `pGeneralizedCocombine()` would be
harmless, but the 2026-07-03 audit found the premise was wrong: `mCochainsMap` is
`Cell< Map<...> >` (hpp:36), and `Cell::operator()` is bounds-checked array indexing that is
free in release builds. Hoisting is cosmetic; the repeated *inner* `Map` lookup of `a` is the
real cost, covered by #7.

### 7. Hoist `a`'s coboundary map out of the neighbor loop — STILL-VALID

**Identified by:** Codex | **Location:** `cl_SimplicialComplex.cpp:1184-1185`

`mCochainsMap(p+1)(a)->getCoboundary()->getSimplicesMap()` is fetched inside the
per-neighbor loop but doesn't depend on the neighbor.

---

## Container/Helper Costs (Lower Impact)

### 8. `erase_key()` double lookup in both `Map` and `OrderedMap` — RESOLVED (2026-06-22)

**Identified by:** Codex, Claude | **Location:** `cl_OrderedMap.hpp:147-151`, `cl_Map.hpp:204-208`

The `key_exists()` guard has been removed; both implementations now erase directly:

```cpp
void erase_key( const Key & aKey )
{
    mMap.erase( aKey );   // std::map / std::unordered_map already no-op if absent
}
```

The original finding (which quoted a guarded `if (key_exists(...)) mMap.erase(...)` body) no
longer applies. Verified in both container headers by Claude, Codex, and Grok.

### 9. `remove_kcochainFromMap()` — CHANGED: now double access, not triple

**Identified by:** Codex | **Location:** `cl_SimplicialComplex.hpp:314-321`

```cpp
delete mCochainsMap(k)[aID];      // operator[] → find (lookup #1)
mCochainsMap(k).erase_key(aID);   // erase only (lookup #2, since #8 removed the guard)
```

With #8 resolved this is now a **double** access (find + erase), not the previously-reported
triple lookup. **Fix:** use `find()` once, `delete it->second`, erase by iterator.

> Note (Codex): the helper uses `operator[]`, which would default-insert a `nullptr` if `aID`
> were absent; it is safe only because callers guarantee presence.

### 10. Stale adjacency entries inflate search cost — STILL-VALID

**Identified by:** Codex | **Location:** `cl_SimplicialComplex.cpp:1135-1143`

The backward scan through a coboundary's `OrderedMap` checks `mCochainsMap(p+1).key_exists()`
for each entry. Stale entries (simplices already removed from the complex but still
referenced in local coboundary maps) increase the number of failed lookups. This invariant
also underwrites the unguarded `operator()` accesses noted in #5.

More aggressive local cleanup during removal, or lazy compaction, could reduce this.

---

## `addCochainToCochain()` Recursive Cost (Structural) — STILL-VALID

**Identified by:** Codex, Grok, Claude | **Location:** `cl_Cochain.hpp:281-307` (called at `cl_SimplicialComplex.cpp:1165`)

Each call recursively updates the cochain, its boundary, AND its coboundary. Inside
`pGeneralizedCocombine()`, this is called once per neighbor, and each call chains through
three levels of `std::map` operations.

This is inherent to the algorithm's algebraic structure. Possible mitigation:
- A specialized non-recursive merge for this specific algorithmic path.
- Separate "update coefficients" from "repair reciprocal adjacency" so each map is
  touched once.

Not a quick fix — requires understanding the invariants carefully.

---

## Evaluation of Grok-Specific Suggestions (still valid)

Several Grok suggestions were evaluated and found to be **invalid or impractical**. None of
these verdicts are affected by the 2026-05-13 code change:

| Suggestion | Verdict | Reason |
|-----------|---------|--------|
| Replace `OrderedMap` with `Vector` + `DynamicBitset` for simplex maps | **Impractical** | Simplex maps store sparse coefficients (+1/−1) keyed by ID. After coreduction, IDs are highly sparse. Dense arrays would waste enormous memory. `DynamicBitset` cannot store coefficient values. |
| Object pooling for Cochains (claimed "2× speedup") | **Premature** | Cochain destructor deletes sub-cochains recursively (`cl_Cochain.cpp:48-51`); pooling requires careful lifecycle management. Worth measuring first, but the claimed speedup is unfounded. |
| `reserve()` on `mCochainsMap` | **Marginal** | Maps are *shrinking* during coreduction, not growing. |
| `try/catch` for map access | **Invalid** | BELFEM compiles release with `-fno-exceptions`. |
| Specific timing numbers ("6-10 s", "3-5× speedup") | **Fabricated** | No measurements were performed. |
| "The rest of the pipeline is very clean" | **Unfounded** | No evidence this was actually analyzed. |

**Note on `mCochainsMap` (top-level) replacement:** Replacing the `Map<index_t, Cochain*>`
with a flat `Cell<Cochain*>` indexed by element index IS potentially viable for the
top-level maps (element indices are bounded), paired with a `DynamicBitset` for existence
checks. This would give O(1) lookups instead of hash table lookups. However, this only
helps the top-level container, not the per-cochain `mSimplicesMap` which is the more
frequently accessed structure.

---

## Recommended Priority Order

(Risk column updated by the 2026-07-03 re-evaluation; see that section for rationale.)

| # | Fix | Type | Risk | Effort | Status |
|---|-----|------|------|--------|--------|
| 1 | Worklist/frontier for `pCoreduce()` (#1) — port to Mrozek & Batko 2009 Alg. 6.1 | Algorithmic | Medium (math de-risked; needs regression meshes) | Medium | open |
| 2 | Persistent worklist for `coreduceOmit()` (#2, option a; option b struck) | Algorithmic | Medium | Medium | open |
| 3 | Snapshot-based split of invariant cleanup (#3) | Algorithmic | Medium-low | Low | **done** (2026-07-03, audit-confirmed) |
| 4 | Eliminate redundant `getCoefficient()` lookups (#4) | Code cleanup | Low | Low | **done** (2026-07-03) |
| 5 | Single `find()` instead of `key_exists`+`operator()` (#5) | Code cleanup | Low | Low | mostly done; 1152/1185 fold into #7; 1187 assert kept as invariant check |
| 6 | ~~Hoist `mCochainsMap(k)` Cell references (#6)~~ | Code cleanup | — | — | struck: premise wrong, impact ≈ nil |
| 7 | Hoist `a`'s `Cochain*`/coboundary out of the neighbor loop (#7) | Code cleanup | Low | Low | **done** (2026-07-03) |
| 8 | Fix `remove_kcochainFromMap()` lookup (#9) | Container fix | Low | Low | **done** (2026-07-03) |
| 9 | Clean stale adjacency entries during removal (#10) | Algorithmic | Medium | Medium | open |
| 10 | Specialized non-recursive merge for `addCochainToCochain` | Structural | High | High | open |
| 11 | Unit-coefficient check in `pCoreduce()` (NEW #11) | Correctness | Low (behavior-narrowing) if confirmed | Low | **deferred** (2026-07-03) — sequenced behind thin-cut rectification T5/T8/T9 |

Fixes 3, 4, 7, and 8 were applied on 2026-07-03 after the Codex audit. Fixes 1-2 remain the
highest-impact open work and need regression meshes (Betti numbers and cut counts before/after).
Fixes 9-10 are deeper structural changes for later; #11 is pending user decision.

## Measurement Before Changing Code (still applicable)

The existing `PERFORMANCE_CHECK` in `pCoreduce()` (`cl_SimplicialComplex.hpp:334-383`) is too
coarse — it only logs loop count and aggregate cochain-map size. Useful counters:

- Number of full-map passes in `pCoreduce()` (reveals rescan waste)
- Cochains examined vs. actually removed per pass (reveals scan efficiency)
- Average/max backward scan length in `pGeneralizedCocombine()` candidate search
- Count of `setCoefficient(..., 0)` calls (reveals erase traffic)
- Time in `addCochainToCochain()` (reveals recursive merge cost)
