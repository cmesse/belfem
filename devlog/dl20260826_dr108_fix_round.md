# DR-108 fixed: the uint8_t hanging-source counter

**Date:** 2026-08-26
**Topic:** fix round for the source-counter overflow that flux-starved the gantry yoke arm.

## What was wrong

`mesh::Basis::mNumberOfSources` was `uint8_t`. `CutSet::create_duplicates` gives each
duplicated cut node a source list of `N + 1` entries (one abstract node per crossing cut,
plus the original node, last). The gantry has 464 cuts, so the trunk duplicates needed 465
sources; `set_sources` narrowed that to 209, `malloc`ed and copied only 209, and the
**original node — last in the list — was dropped**. The intended constraint
`phi_dup = phi_org + sum_c I_c` degenerated to `phi_dup = 209 * I`, the cut stopped
transmitting flux, and the yoke arm behind the trunk went dead.

Root cause was established by two independent measurements: the φ-ladder census (32 N-groups,
32 exact `(N+1) mod 256` hits) and a live write-site probe (797 calls above 200 sources;
426→170, 430→174, … 465→209, all wraps observed as they happened). The `465 = 464 cuts + 1
original` decomposition is what identifies the dropped entry.

## What landed

`src/mesh/cl_Mesh_Basis.{hpp,cpp}`, +58/-1:

- `mNumberOfSources` widened to `uint16_t`. `sizeof(Basis)` is unchanged — the member sits in
  the 8-byte slot ahead of `Basis** mSources`. `uint16_t` over `uint` deliberately: it keeps
  the allocation guard able to fire, and matches `mFacetCounter`, which had already been
  widened *because of cohomology cuts* — the same pressure, on the neighboring counter.
- `BELFEM_ERROR` width guards before every narrowing write: the three `set_sources` overloads,
  an allocation-size guard in `allocate_source_container`, and a no-new-member wrap guard in
  `add_source`.
- `mNumberOfDofs` stays `uint8_t` but is guarded in `increment_dof_counter` and `insert_dof`.
  Audit upgraded this from judgment call to required: a wrapped dof counter makes
  `allocate_dof_container` skip its `malloc`, and `insert_dof` then writes through `mDofs`,
  which has no `nullptr` initializer.
- Every limit via `std::numeric_limits<decltype(member)>::max()`; no hardcoded cliffs.
- No `.bfm` load guard, deliberately: a wrapped file is *below* any width check and cannot be
  detected there. Stale-`.bfm` deletion is a gate condition instead.

## Gates

- **G1** serial gantry, fresh mesh, stale bfm deleted: census **0 corrupted nodes of 2334**.
  The 183 trunk nodes formerly pinned to exactly 209·I regardless of position now span
  0.4621..159.6016·I with real spatial structure.
- **G2** np = 10, fresh build: zero master-not-found aborts, residual sequence **identical to
  serial** (0.00/−71.33, −4.77/−79.07, −4.77/−77.34, −4.77/−75.62 dB).
- **G3** `make check-fast` 9/9.

Unplanned corroboration: convergence went from 6 iterations per step (Picard 1–2 then Newton
3–6) to **2 Picard iterations** at −71…−85 dB. The wrapped constraint had been making the
solver fight an inconsistent system.

## Round notes

Plan and code both audited by Codex and Grok, 2/2 each. Their corrections shaped the result:
the "four set_sources overloads" in my plan was a miscount (there are three), the dof-counter
guard was upgraded to required, and the `.bfm` question was settled as a process gate rather
than a code guard. On the plan's one genuine vendor split — whether a `userdefined` section
could produce per-member values — I read the construction site myself: the `SourceFunction`
is built inside the per-member loop, so instances really are per-member (Codex right on the
fact), but all are initialized from the same `(file, label)` of the same section, so
divergence requires an impure user library (Grok right on the consequence).

Three claims from my own opening brief on this defect were refuted along the way and are
corrected in `todo/ferro_cut_flux_island.md`: the "inconsistent per-fragment MMFs" were the
wrapped constant plus a floating island's offset, the "mixed jump signs" were an artifact of
unordered pair enumeration, and "iron-piercing cuts" were never the mechanism — 205 of the 302
corrupted nodes are in air, and the trunk begins on the outer air boundary. Christian's
ferro-cut input flag would not have fixed this; it stands on its own merits.

## Also in this session

`cmake-build-claude/` was set up as a private build tree so gates stop queueing behind
Christian. It ran G1–G3 here and cleared the two gates that commit c6ef67ef had declared
owed (check-fast 9/9; fresh-build gantry np=10 distributing with zero master-not-found).
