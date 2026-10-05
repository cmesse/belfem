# 0.9.3 ultra-review follow-ups: dangling mesh dof containers after a refusal, pure-metal range assert

**Date:** 2026-10-04
**Purpose:** Close the two findings of the `/code-review ultra` pass over the 0.9.3 branch that
Christian ruled to fix (2026-10-04, question form): `DofData::reset` deletes the dofs of a refused
field but leaves the mesh's node/edge/face/facet/element dof containers pointing at them, and the
pure-metal classification of `BuiltinName` is a range test guarded only by a comment. The third
finding (the duplicated constructor try/catch cleanup in three factories) is declined as style.
**Module:** `src/fem/kernel` (`DofData`), `src/physics/materials` (`MaterialFactory`); tests in `tests/fem`
**AIs involved:** Claude Fable 5.1 (plan, pre-registration, reconciliation, commit); Claude Opus (implementation agent); Codex `gpt-6-astra`/high + Grok `grok-4.7`/high (plan jury); Codex `gpt-6-astra`/xhigh + Grok `grok-4.7`/xhigh (code jury, round ≥ 2 of the subject)
**Status:** ✅ COMPLETE 2026-10-05 — committed on `devel`. Landed: `DofData::reset` disconnects before deleting; `disconnect_dofs_from_mesh` releases a basis if it holds the dof, is empty, or has no array; `connect_dofs_to_mesh` resets its bases before counting; `Basis::holds_dof` / `has_dof_container`; the `is_pure_metal` `BELFEM_ASSERT`; five `DofDataReset` tests (nine in all). Verified: `make check` 17/17 debug + 19/19 release in private trees on the final diff, valgrind clean, every positive control and mutation red then restored. Juries: plan (astra/high, grok-4.7/high), code ×3 (astra/xhigh, grok-4.7/xhigh; Codex round 2 lost to the usage limit, re-run on the final diff). Residue: `todo/test_gap_followups.md` R13 (hanging-edge fixture, `tBefore` leak); `staticheat.cpp` order if revived.

> **Scope guards:**
> - In scope: `DofData::reset`, `DofData::disconnect_dofs_from_mesh`, the `DofDataReset.*` tests, `is_pure_metal`.
> - Out of scope: the `tBefore` leak on the same refusal path (`cl_FEM_DofMgr_DofData.cpp:1950`, `:2027`) and the
>   full-path hanging-edge refusal fixture — both stay in `todo/test_gap_followups.md` R13; the factories'
>   try/catch idiom (declined); `static_assert` anywhere (Christian: `BELFEM_ASSERT`, never `static_assert`);
>   the disabled `staticheat.cpp` example, which deletes its mesh before its kernel (D2).
> - No `input.conf` key changes; no behavior change on an accepted deck.

---

## 1. Current behavior and how it fails

| Failure | Mechanism | Evidence |
|---|---|---|
| Use-after-free hazard (throw mode) — mesh entities keep pointers to freed dofs after a refused field setup | `create_dofwise_t_matrices_master` connects every dof to its mesh basis (`connect_dofs_to_mesh`, `cl_FEM_DofMgr_DofData.cpp:2977` → `Basis::insert_dof`) and can refuse later in the same call (hanging-edge weights other than ±1, `:3236`, `:3243`, `:3355`); the refusal unwinds into `Kernel::create_field` (`cl_FEM_Kernel.cpp:609-620`), whose owner destroys the manager → `DofData::~DofData` (`:62-66`) → `reset` (`:70-118`), which deletes the dofs but never calls `disconnect_dofs_from_mesh` (`:669-685`). Release builds abort inside `BELFEM_ERROR` and never reach `reset`; the hazard is a throw-mode (debug, gtest) one | `cl_FEM_DofMgr_DofData.cpp:62-118`, `:647-685`, `:2977`; `cl_FEM_Kernel.cpp:596`, `:609-620`; `cl_Mesh_Basis.hpp:301-308` |
| Silent misclassification if the enum is edited — a non-metal inserted before `Lead` would pass as an alloy component | `is_pure_metal` is `aName <= BuiltinName::Lead` (`cl_MaterialFactory.cpp:111-114`); the contract "the pure metals come first" lives in the comment at `:48-49` only | `cl_MaterialFactory.cpp:44-66`, `:111-114`, `:716`, `:732-735`; consumer `cl_Material_Alloy.cpp:179` |

**Where the containers are live (corrected by the plan jury, D1).** `DofManager::set_equation`
(`cl_FEM_DofManager.cpp:85-128`) connects and never disconnects; `consolidate_dofs` (`:1250-1253`)
and `link_bearings_with_dofs` (`:1345-1361`) read the containers inside it, `MaxwellFactory`'s
autopinning reads `tNode->dof( 0 )` after it (`cl_MaxwellFactory.cpp:947`), and the disconnect is
the last step of `DofManager::initialize` (`:298-325`, the call at `:320`). `Kernel::create_field`
returns without initializing (`cl_FEM_Kernel.cpp:625`), so a field sits in the connected window
until its factory calls `initialize` (thermal: `cl_ThermalFactory.cpp:236`) or disconnects by hand
(`cl_MaxwellFactory.cpp:1008`). A second field on the same mesh entering that window replaces the
containers: `allocate_dof_container` frees the old array and allocates a new one
(`cl_Mesh_Basis.hpp:285-294`), so the first field's links are gone from that moment — a pre-existing
property, not touched here.

Who reads the dangling pointers today: no shipped path found by three readers (me, Codex, Grok);
the next `connect_dofs_to_mesh` increments before it allocates (`:649-651`), so the stale slots are
freed, not loaded. Grok did not inventory every `Basis::dof` caller that could run after a caught
`create_field`; neither did I. Confidence medium-high that the hazard is latent.

**Bottom line:** in throw mode, `reset` is the teardown every refused setup reaches, so it must
release the mesh containers the setup filled — and, because a sibling field may own those
containers by then, it must release only the ones that still hold this manager's dofs. The enum
range contract should fail loudly in a debug build instead of silently.

## 2. Mechanism

**A — `reset` disconnects, ownership-aware, before it deletes.**

`DofData::disconnect_dofs_from_mesh` (`:669-685`) changes in two ways:

1. a null guard per slot (`create_dofs` pre-fills `mDOFs` with `nullptr`, `:421`;
   `collect_hanging_dofs` does the same for `mHangingDOFs`, `:1977`, `:2016`; a `new Dof` that throws
   inside a fill loop reaches `reset` with a null suffix);
2. a basis releases its container only if that container still holds the dof being disconnected:
   scan `aBasis->dof( i )` for `i < aBasis->number_of_dofs()` and compare pointers (no dereference
   of the stored pointers, so a container of already-deleted dofs is safe to scan). A container
   that a second manager replaced does not hold this manager's dof and is left alone (D1). The scan
   is at most 255 entries (`mNumberOfDofs` is `uint8_t`, `cl_Mesh_Basis.hpp:46`) per dof, teardown
   only. On the normal path (`initialize:320`, `MaxwellFactory:1008`) every basis holds the
   manager's own dofs, so every container is released as before.

   **Amended by the code jury (D4, D5):** the scan lives in `Basis::holds_dof` (null-guarded, so a
   counted-but-unallocated basis is skipped) and the release predicate is
   `number_of_dofs() == 0 || holds_dof( tDof )`: a count of 0 has no readable entry, so releasing
   it never drops a sibling's live links, and it frees the empty arrays a sibling's connect leaves
   behind. `Basis::has_dof_container` exists for the tests' oracle. The invariant the predicate rests
   on: a container with entries belongs to exactly one manager, because connect replaces before it
   inserts.

`DofData::reset` calls `this->disconnect_dofs_from_mesh()` immediately before the flag pass, while
every pointer is still live. After `initialize` has disconnected, the second call finds empty
containers (`free( nullptr )`, count 0 → no scan) and is a no-op. Before `initialize` — a field
destroyed in the connected window, which is exactly the refusal case and also a clean teardown of a
never-initialized field — the call is the fix. A dof always carries a mesh basis
(`cl_FEM_Dof.cpp:23,41,60,78,95`), so only the dof pointer is guarded. The mesh must outlive the dof
manager for the destructor to read the bases: `Kernel::~Kernel` deletes its managers (`:128-131`)
before an owned mesh (`:146-153`), and `HangingStrip` declares the mesh first; the contract is
written into the `DofData` class comment (D2).

Alternative rejected: a disconnect in `Kernel::create_field`'s catch block — it would cover one
owner and not the destructor; and an unconditional disconnect (plan v1) — it frees a sibling's
containers (D1).

**B — a `BELFEM_ASSERT` pins the range.** In `is_pure_metal`, before the comparison:
`BELFEM_ASSERT( static_cast< uint >( BuiltinName::Lead ) + 1 == 9, "the pure metals of BuiltinName must be the nine names Aluminum .. Lead; update is_pure_metal() with the enum" )`.
An insertion before `Lead` changes the count and trips it in every debug run that reaches the
function — the alloy-formula scan (`cl_MaterialFactory.cpp:704-716`, reached by `pb40sn60`-shaped
labels; plain names return from `builtin_kind` at `:664` first, D3) and `is_pure_metal_label`
(`:732-735`, every alloy component check). A reorder that keeps `Lead` at 8 passes, which the
comment still has to carry — the assert is a tripwire, not a proof. `static_assert` is excluded by
ruling. Debug-only is the right tier: the fact is a compile-time property, so a release message
buys nothing.

## 3. Steps

- [x] **R1** A: null guard and ownership scan in `disconnect_dofs_from_mesh`; `reset` calls it first; the mesh-lifetime contract in the class comment; the comment above the flag pass extended by one sentence. `cl_FEM_DofMgr_DofData.cpp` (+ `.hpp` for the class comment if it lives there).
- [x] **R2** A tests, `tests/fem/test_DofSelfContainment.cpp`, beside `DofDataReset.*` (the fixture's field has already filled every node's container through `create_field`, so each test replaces the container with the full trio `increment_dof_counter` → `allocate_dof_container` → `insert_dof`, exactly as `connect_dofs_to_mesh` does; a bare `insert_dof` would add to the field's count):
  - `MeshContainersAreReleased`: three counting dofs on three distinct nodes — one in both containers (the shared shape), one in `mDOFs` only, one in `mHangingDOFs` only (D3: the hanging-only basis must be cleared by the hanging loop, not by the `mDOFs` loop). Each node's count is 1 before `reset`, 0 after; each dof destroyed once.
  - `AnotherManagersContainerIsKept`: two `DofData` objects on the fixture's field. A's dof on node 1 is inserted; then B's dof replaces node 1's container (the trio). `reset` of A: node 1 still has count 1 and `dof( 0 )` is B's dof; A's dof destroyed once. `reset` of B: count 0, B's dof destroyed once.
  - The null-slot shape is already covered by `NullSlotsAreSkipped` (`:320-338`), which now exercises the guard; the "never-connected basis" case of plan v1 is dropped (cannot fail on this fixture, Grok).
  - Positive controls (R4), test hunk only on HEAD `77b4567`: `MeshContainersAreReleased` fails (counts stay 1); `AnotherManagersContainerIsKept` passes on HEAD (no disconnect at all) — its control is plan v1's unconditional disconnect, which the implementation agent is told not to write; I check the ownership scan by reading the diff and by a one-line mutation (scan disabled → that test red) in the private tree.
- [x] **R3** B: the `BELFEM_ASSERT` in `is_pure_metal`, `cl_MaterialFactory.cpp:111-114`; comment at `:48-49` points at it. No new test: the `pb40sn60` / `Pb38Sn62` cases of `test_MaterialFactoryDeck` and every alloy-component check reach it.
- [x] **R4a** Baseline `make check` on HEAD `77b4567` in both private trees (reused from the pre-check session; debug: SCLS `debug`, reference BLAS, Armadillo, 17/17, zero warnings; release: SCLS `mkl`, Armadillo, gas models, shared libs, 19/19). Done 2026-10-04.
- [x] **R4b** With the diff, three passes in Claude's own private trees (2026-10-04/05): pass 1 (pre-D4/D5) debug 17/17 + valgrind 6/6 clean, release 19/19, positive control both new tests red on HEAD, scan mutation red; pass 2 (post-D4/D5) 17/17 + valgrind 8/8, 19/19, positive control on the pre-fix predicate SIGSEGV + red, count-0 mutation red; pass 3 (post-D7/D8) and pass 4 (final diff) 17/17 + valgrind 8/8, 19/19, zero warnings both, no-array mutation red. Every mutation restored to green afterwards.
- [x] **R5** Code jury on the own-files diff: round 1 → D4–D6; round 2 (Grok only, Codex on the usage limit) → D7–D8; round 3 on the final diff, both voices, no defect, one P3 test gap closed by `AFullBasisIsReleasedBeforeItIsCountedAgain`. (`tmp/ai_exchange/review_release093_ultra_followups.diff`), same slug, round-2 depth (astra/xhigh, grok-4.7/xhigh); verification pass; reconciliation.
- [x] **R6** Close out (2026-10-05): devlog, `devlog/README.md`, this file's Status, `todo/README.md`, `todo/test_gap_followups.md` R13 (narrow it to the fixture, the `tBefore` leak and the full-path regression — the dangling slots are closed here), commit on `devel`.

## 4. Defects found by the juries

- **D1 (HIGH, Codex P1 + Grok 1–3, confirmed by source read)** Plan v1 claimed `set_equation` disconnects at `:320`; that line is `initialize`. Consequence: an unconditional disconnect in `reset` frees containers a sibling field may own (sequence: initialize A, create B on the same mesh, destroy A — no shipped caller does it, `Kernel::~Kernel` destroys all managers together and `create_field`'s catch destroys only the refused one, but the ownership hole is real). Fixed in plan v2 by the ownership scan (§2 A.2) and the two-manager regression (R2).
- **D2 (LOW, Codex P2)** The new destructor reads the mesh; a mesh deleted before its manager would be a use-after-free in `~DofData`. No shipped path (`Kernel::~Kernel` order; `HangingStrip` order); the disabled `staticheat.cpp` example (`:59`, `:133`; off in `src/fem/kernel/CMakeLists.txt:80`) does delete its mesh first and would need reordering if it were ever re-enabled. Plan v2: the contract goes into the `DofData` class comment; the example stays out of scope.
- **D3 (LOW, Codex P2 #3 + Grok 7–8)** R2 of plan v1 would have passed without the hanging loop and its "unconnected basis" half could not fail; "every `test_MaterialFactoryDeck` case" reaches `is_pure_metal` was wrong (plain names return from `builtin_kind` first). Fixed in R2/R3 wording.
- **D4 (HIGH, Codex P1 code round, confirmed)** The ownership scan indexed `Basis::dof( k )` on a basis whose count is positive but whose array is null — the state `allocate_dof_container` leaves when its `malloc` fails (the `BELFEM_ERROR` precedes the count reset, `cl_Mesh_Basis.hpp:283-296`) and the state `increment_dof_counter` alone produces. Fixed 2026-10-04 (Opus, second pass): the scan moved into `Basis::holds_dof`, null-guarded; test `ACountedContainerWithoutAnArrayIsSkipped`.
- **D5 (MEDIUM, Grok 1 code round, confirmed)** A sibling's `connect_dofs_to_mesh` runs `allocate_dof_container` over every basis and leaves count-0 arrays where it inserts nothing; the scan never released a count-0 container, so they were retained until `~Basis` (plan v1's unconditional disconnect freed them). Fixed 2026-10-04: the predicate is `count == 0 || holds_dof` — a count of 0 has no readable entry, so it is never a sibling's live links; test `AnEmptyAllocationIsReleased` with the new `Basis::has_dof_container` accessor as oracle.
- **D6 (LOW, Codex P3 + Grok 3–5)** Three comments overstated: "deleted pointers are safe to compare" (the scan is safe because `reset` disconnects before it deletes), the class comment's unconditional "releases the mesh's dof containers", the enum comment's "pinned by the assert" (it checks the range length only). Reworded 2026-10-04. Declined: `static_assert` (ruling), `k` → `tIndex` (file idiom, index exception).
- **D7 (MEDIUM, Grok defect 1 code round 2, confirmed)** `holds_dof` trusts the counter as the array length, but `connect_dofs_to_mesh` increments every counter before any allocation; a throw in that loop (the `uint8_t` ceiling, debug build) leaves earlier bases with a counter above a sibling's stale array, and the failing manager's `reset` would scan past it. Needs a 255-dof basis, which no mesh has; the API produces the state. Fixed 2026-10-05 (Opus, third pass): `connect_dofs_to_mesh` resets each of its dofs' bases before counting, so no basis carries a counter above an existing array and a refused connect leaves no counter residue. Success path unchanged (`allocate_dof_container` freed the same arrays a few lines later). No distinguishing oracle short of a 255-dof basis; reviewed by trace, exercised by every connect in `make check`.
- **D8 (LOW, Grok defect 2 code round 2, confirmed)** `ACountedContainerWithoutAnArrayIsSkipped` could not tell skip from release. Fixed 2026-10-05: predicate `! has_dof_container() || count == 0 || holds_dof` — a basis with no array has nothing to protect and gets its counter cleared; test renamed `ACountedBasisWithoutAnArrayIsCleared`, asserts `number_of_dofs() == 0` after `reset`.
- **Codex leg of round 2** failed on the ChatGPT usage limit (reset 00:40 PDT); re-run scheduled on the final diff with the identical prompt. Round 2 was one voice; round 3 (Grok on the final diff) plus the Codex re-run restore two.
- **Grok 4–6 (prose)** O1 misdescribed the worker path (`connect_dofs_to_mesh` has no rank test; it is master-only because `:1985` sits in the `mCommRank == 0` branch from `:1890`; a rank-0 refusal leaves the workers in the barrier at `:2192`, so worker `reset` never runs on the refusal path); "every exit" must say throw mode; "overwritten" must say freed. Fixed in §1.

## 5. What could break

- A: a caller that destroys a manager while a sibling holds the containers now keeps the sibling's links (the scan); the unconditional variant would have dropped them. A caller that leaves two fields in the connected window on one mesh has already lost the first field's links to `allocate_dof_container` — pre-existing, unchanged.
- A: the mesh-lifetime contract (D2).
- A: parallel refusal is a deadlock at barrier 1 (`:2192`) independent of this fix; R13's fixture must stay serial (Grok).
- B: `BELFEM_ASSERT` with a constant-true condition under the debug tree's `-Wall -Werror` (`CMakeLists.txt:110`, `config/compiler/config_gcc.cmake:86`); `-Wtautological-compare` covers self- and range-compares, not a constant equality — the debug build (R4b) decides.

## 6. Open questions

- ~~**O1** worker path~~ — RESOLVED 2026-10-04 by the plan jury (Grok 4): no worker reaches `reset` on the refusal path; on a successful worker teardown the bases were never filled on that rank and the call is a no-op. No decision needed.

## 7. Exchange

`tmp/ai_exchange/review_release093_ultra_followups.md` (plan round, verification, reconciliation; code round follows in the same thread).
