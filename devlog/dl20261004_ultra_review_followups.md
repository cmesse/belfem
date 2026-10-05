# Ultra-review follow-ups for 0.9.3: `DofData::reset` releases the mesh dof containers; pure-metal range assert

**Date:** 2026-10-04 (overnight run, closed early 2026-10-05)
**Purpose:** Christian tried `/code-review ultra` on the 0.9.3 branch for the first time (the whole `devel`–`main` span exceeded its 8 000-line limit; the `devel`-scoped run covered 77 files). It returned three findings; two were ruled to fix overnight through the standing loop, the third (the duplicated constructor try/catch cleanup in three factories) declined as style. Question: close the dangling mesh dof containers that `DofData::reset` leaves after a refused field setup, and pin the pure-metal range of `BuiltinName`, without changing any accepted deck's behavior.
**Module:** `src/fem/kernel` (`DofData`), `src/mesh` (`Basis`, two inline accessors), `src/physics/materials` (`MaterialFactory`); tests in `tests/fem`
**AIs involved:** Claude Fable 5.1 (plan, pre-registration, verification, reconciliation, gates, commit); Claude Opus (two implementation passes from written briefs); Codex `gpt-6-astra`/high + Grok `grok-4.7`/high (plan jury); Codex `gpt-6-astra`/xhigh + Grok `grok-4.7`/xhigh (code jury round 1 on the first diff, round 2 on the post-fix diff)
**Status:** COMPLETE — committed on `devel`; §Owed lists what stays open
**Verification:** executed in Claude's own private trees (not the shared `build/` or `cmake-build-debug/`): debug (SCLS `debug`, reference BLAS, Armadillo) and release (SCLS `mkl`, Armadillo, `-DUSE_GASMODELS=ON`, shared libs). Baseline on the committed HEAD `77b4567`: 17/17 debug, zero warnings; 19/19 release. First pass (pre-D4/D5 diff): 17/17 and 19/19, valgrind on `DofDataReset.*` clean (every leak record under OpenMPI init, no BELFEM-owned allocation, no invalid access), positive control red (both new tests fail on HEAD with only the test hunk), mutation red (scan answering true for every basis turns `AnotherManagersContainerIsKept` red). Second pass (post-fix diff): debug build zero warnings, `make check` 17/17, valgrind on `DofDataReset.*` 8/8 pass, 28 error contexts all under OpenMPI init, zero BELFEM frames outside `Communicator::init`, no invalid access; release build zero warnings, 19/19. Positive control on the pre-fix predicate (round-1 `basis_holds` + the new header + the new tests): `ACountedContainerWithoutAnArrayIsSkipped` dies with SIGSEGV (the null array indexed, Codex's P1), `AnEmptyAllocationIsReleased` fails on `has_dof_container()` (Grok's finding). Mutation (count-0 clause dropped): `AnEmptyAllocationIsReleased` red alone; restored, 8/8. Third pass (after D7/D8) and fourth pass (the final diff, after the round-3 comment fixes and the extra counter assertion): each time debug zero warnings, 17/17, valgrind 8/8 clean (28 contexts, all OpenMPI init, zero BELFEM frames); release zero warnings, 19/19; mutation (no-array clause dropped): `ACountedBasisWithoutAnArrayIsCleared` red alone; restored, 8/8. Fifth pass (the committed diff, with the pre-pass regression): debug zero warnings, 17/17, valgrind 9/9 clean; release zero warnings, 19/19; positive control (pre-pass loop removed): `AFullBasisIsReleasedBeforeItIsCountedAgain` fails with "Dof counter of basis 1 is full ( max 255 )"; restored, 9/9. Both trees at `e5a5873` + the committed diff. Christian's own `make check` in the shared trees not run.
**Exchange:** `tmp/ai_exchange/review_release093_ultra_followups.md` (plan round, code rounds 1 and 2, verifications, reconciliations), `review_release093_ultra_followups.diff`

## Summary

`reset` is the teardown every refused field setup reaches in throw mode, and until now it deleted the
dofs while the mesh's node/edge/face/facet/element containers still pointed at them. It now calls
`disconnect_dofs_from_mesh` first. The disconnect is ownership-aware: a basis releases its
container if the container still holds the dof being disconnected, or if the container is empty.
That predicate is what the two juries made of my plan — my first version disconnected
unconditionally on a wrong premise, the second scanned without a null guard and never released
empty arrays. The pure-metal range test gained a `BELFEM_ASSERT` on the range length (never a
`static_assert`, Christian's ruling). No accepted deck changes behavior: on the shipped paths every
container is already released at the end of `DofManager::initialize` or by `MaxwellFactory`, and a
second disconnect finds count 0.

## What the ultra review found, and what it did not

The three findings, verified against the source before the ruling: the dangling containers
(pre-existing, already known as test-gap R13 from the pre-check's code round — the ultra pass added
no new fact but confirmed the row), the comment-only enum ordering contract, and the duplicated
try/catch idiom. No P0, no P1; it agreed with the three-vendor pre-check of the same morning.

## What the juries changed

**Plan round (astra/high, grok-4.7/high).** Both refuted my lifecycle trace: the disconnect at
`cl_FEM_DofManager.cpp:320` is the last step of `initialize`, not of `set_equation`, so a field sits
with live containers from `set_equation` until `initialize` — and `Kernel::create_field` returns
inside that window. An unconditional disconnect in `reset` would free a sibling field's containers
(Codex P1: initialize A, create B, destroy A). Also: the new destructor reads the mesh, so the mesh
must outlive the manager (Codex P2 — `Kernel::~Kernel` and the test fixture satisfy it; the disabled
`staticheat.cpp` example does not and stays out of scope); my R2 test would have passed without the
hanging-dof loop (Codex, Grok); my O1 worker story and three prose claims were wrong (Grok). I
proceeded on plan v2 despite my own pre-registered "stop on P1" note, because the P1 was a defect
in the plan's trace with a correction both voices converged on, not a design ruling; recorded in
the exchange.

**Code round 1 (astra/xhigh, grok-4.7/xhigh) on the first diff.** Codex P1: the ownership scan indexed
`Basis::dof( k )` on a basis with a positive count and a null array — the state `allocate_dof_container`
leaves when its `malloc` fails (the `BELFEM_ERROR` precedes the count reset) and the state a bare
`increment_dof_counter` produces. Grok 1: a sibling's `connect_dofs_to_mesh` runs
`allocate_dof_container` over every basis and leaves count-0 arrays where it inserts nothing; the
scan never released those, so they were retained until `~Basis`. Both P3-level comment overclaims
(Codex, Grok) confirmed. Declined: `static_assert` (ruling), `k` → `tIndex` (the file's 38 index
loops all use a bare `k`; CLAUDE.md's index exception).

**Fixes (Opus, second brief).** The scan moved into `Basis::holds_dof` (null-guarded, pointers
compared, never dereferenced) with `Basis::has_dof_container` as the tests' oracle; the predicate
became `number_of_dofs() == 0 || holds_dof( tDof )` — a count of 0 has no readable entry, so
releasing it never drops a sibling's live links, and it frees the empty arrays a sibling's connect
leaves. Two tests: `ACountedContainerWithoutAnArrayIsSkipped`, `AnEmptyAllocationIsReleased`. The
invariant the predicate rests on — a container with entries belongs to exactly one manager, because
connect replaces before it inserts — is stated in the disconnect comment (Grok's open risk).

**Code round 2 (astra/xhigh, grok-4.7/xhigh) on the post-fix diff.** Codex's leg died on the ChatGPT
usage limit (reset 00:40 PDT; the jury's prompt was copied from the live temp file and a detached
runner re-issues it after the reset, against the final diff). Grok, alone this round, found one more
real gap: `connect_dofs_to_mesh` increments every counter before any allocation, so a throw inside
the counting loop (the `uint8_t` ceiling, debug build) leaves earlier bases with a counter above a
sibling's stale array, and `holds_dof`, which trusts the counter as the array's length, would read
past it in the failing manager's `reset`. No mesh has a 255-dof basis, but the API produces the
state, as D4 did. Grok also showed that the skip test could not tell skip from release. Fixes
(Opus, third brief): `connect_dofs_to_mesh` resets each of its dofs' bases before counting — a
no-op on the success path, where `allocate_dof_container` freed the same arrays a few lines later —
so no basis ever carries a counter above an existing array and a refused connect leaves no counter
residue; the predicate gained a no-array clause (`! has_dof_container()`: nothing to protect, the
counter is cleared) and the test, renamed `ACountedBasisWithoutAnArrayIsCleared`, asserts the
counter. `static_assert` was raised for the third time and declined for the third time.

**Code round 3 (grok-4.7/xhigh; astra/xhigh deferred to the reset) on the final diff.** Grok: no new
code defect; the three predicate arms, their short-circuit order, and the pre-pass hold, and each of
the four new tests locks one arm. It refuted two sentences of my own pre-registration — the pre-pass
makes the malloc exact-size rather than "identical", and sibling-only bases are handled by the
count-0 arm, not the pre-pass — and it showed that after the pre-pass a throw in the counting loop
drops the previous manager's arrays on every shared basis, where the old code dropped nothing
before allocate (reachable only at the 255-dof ceiling; accepted as the price of a safe scan). Three
comment overclaims and one missing counter assertion were fixed after the round; the diff and the
saved Codex prompt were regenerated so the deferred leg audits the final state, and a fourth gate
pass ran on it. Codex, after its reset and on the same final prompt: approve at review level, no blocking
defect, one P3 — the pre-pass had no distinguishing regression — with the oracle supplied: a basis
filled to the 255-dof ceiling, then a second manager connecting one dof there succeeds only because
connect resets the basis first (`increment_dof_counter` refuses at 255). Added as
`AFullBasisIsReleasedBeforeItIsCountedAgain` (Opus, fourth brief), with its positive control in the
fifth gate pass. Two voices on the final diff, in agreement.

## Tests

`DofDataReset.*` in `tests/fem/test_DofSelfContainment.cpp` grew from four to nine:
`MeshContainersAreReleased` (three dofs on three nodes — shared, owned-only, hanging-only — each
node's count 1 before, 0 after, each dof destroyed once), `AnotherManagersContainerIsKept` (two
`DofData` on one field; B replaces node 1's container; A's reset keeps B's link; B's reset releases
it), `ACountedBasisWithoutAnArrayIsCleared`, `AnEmptyAllocationIsReleased`, and
`AFullBasisIsReleasedBeforeItIsCountedAgain` (255 dofs on node 1, a second manager connects one;
red without the pre-pass). The fixture fact that shaped them: `HangingStrip` calls `create_field`
without `initialize`, so every node already holds the strip field's dofs, and a test must perform
the full `increment_dof_counter → allocate_dof_container → insert_dof` trio to replace a container
(a bare `insert_dof` would add to the field's count and the test would lie — Grok, plan round).

## Process notes

- I borrowed the pre-check session's two configured private trees to save an hour of building. That
  session closed mid-run and swept its scratchpad: both `src` worktrees vanished while my debug
  `make` and release `make check` were running. The positive control and the release build had
  completed; the rest was redone in my own trees. Lesson saved to memory: never reuse a peer
  session's scratchpad trees.
- The peer session also committed `e5a5873` (its devlog's commit line) on top of `77b4567` while
  this run was in flight; the diff was rebased onto it with no conflict.
- Codex's round-1 leg finished in six minutes; Grok's xhigh legs took about twenty each. Codex's
  round-2 leg died on the ChatGPT usage limit; the jury's prompt was copied from the live temp file
  and re-issued from a detached runner after the 00:40 reset, against the regenerated final diff, so
  the record has two voices on the state that was committed. Three code rounds in all; every
  post-round edit was a comment or a test, except the D7/D8 predicate change, which went back
  through a round.
- The ultra review's one substantive finding was already in R13; its value here was as a cheap
  breadth pass that agreed with the vendor jury. Its two nits were worth a `BELFEM_ASSERT` and a
  declined refactor.

## Owed

- `todo/test_gap_followups.md` R13, narrowed this session: the hanging-edge h-curl fixture for the
  full refusal path and the `tBefore` leak; serial only (a rank-0 refusal leaves the workers in the
  barrier of `collect_hanging_dofs`).
- `staticheat.cpp` (disabled target) deletes its mesh before its stack kernel; if it is ever
  re-enabled, the order must change (Codex, both rounds).
- `CITATION.cff` version and date with the 0.9.3 tag (unchanged from the pre-check).
