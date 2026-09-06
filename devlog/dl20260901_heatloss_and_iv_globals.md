# Devlog 2026-09-01 — total ohmic heat loss as a mesh global; I/U published alongside the csv

**Date:** 2026-09-01
**Topic:** Two features Christian added between commits — a `heatloss` mesh global integrating
`∫ρ|j|²dV` over the h-domain, and per-terminal `I_nn`/`U_nn` mesh globals written next to the
existing `iv_results.csv` — put through a three-AI jury round, then the P0/P1 repairs that
followed, including a **0-byte-csv regression introduced by the first repair round**.
**AIs involved:** Claude (Opus 5), Codex `gpt-5.6-terra`/`high`, Grok `grok-4.6`/`high`
**Claude Confidence:** high on every defect below — each is a source trace, and the two that
mattered most were closed with standalone reproducers rather than by reading
**Verification:** regression (ladder level 2) — **`make check` green, 15/15** (8 fast, 2 mpi;
54.4 s), run by Christian on branch `main` @ `7e78fac0` + this working tree. Supported by two
standalone probes compiled and executed in the scratchpad (level 4): the `ofstream` double-open
reproducer and the label-agreement probe, and by a three-AI source trace with every citation
re-checked. **Be clear what the suite establishes here:** it proves these changes break nothing it
covers, and it contains **no heat-loss and no IV case** — a 0-byte `iv_results.csv` and a
monotonically growing `heatloss` both pass it. The discriminating evidence for the two P0s and for
the csv regression is the probes, not the suite. Still owed for the features themselves: one
two-rank run confirming `iv_results.csv` is non-empty with an `I_01 U_01 …` header, and that
`heatloss` in the Exodus frames does not drift upward step over step.
**Exchange:** `tmp/ai_exchange/review_heatloss_and_iv_globals.md` (pre-registration, both auditor
entries, verification pass, reconciliation table, round-2 re-check)

## Summary

The two features are small — one `+=` in five element kernels, one accumulator, one label loop —
and the round found **two P0s and four P1s** in them, all agreed 3/3 or confirmed by trace. The
recurring shape is the same in both features: *a diagnostic that is coupled to "whatever happens
to run" rather than to the thing it claims to measure.* `heatloss` accumulates from any code path
that assembles an element; the I/U globals are refreshed by whatever ran most recently before the
file was written. Neither is wrong in the kernel; both were wrong in the plumbing.

Two of my own round-1 claims did not survive verification, and one more was withdrawn under
Christian's pushback. Recorded below, because a round that only records the auditors' errors is
not an honest record.

## 1. What the features do

**`heatloss`** — `MaxwellFactory` creates the global on **every** rank
(`cl_MaxwellFactory.cpp:1035-1039`), deliberately outside the `if ( comm_rank() == 0 )` block that
governs every other mesh global, because workers must own a slot to `+=` into. Each h-kernel
accumulates `dotQ += rho * dot( j, j ) * wdV` over its integration points and calls
`save_heatloss` (`mt_maxwell_h.hpp:32-36`), an unconditional `+=` — five call sites,
`mt_maxwell_h.cpp:155/226/335/377/423`. The integrand is the same one the thermal Joule source
already uses (`mt_thermal_h.cpp:49`), which is the strongest single piece of corroboration the
round produced.

**`I_nn`/`U_nn`** — `Controller::save_IV` publishes one global per abstract dof beside the csv row
it already wrote (`cl_FEM_Controller.cpp:4186-4211`).

## 2. The two P0s

### 2a. The accumulator was reset and reduced on the minority solve path only

The reset/collect pair existed only in `iterate_magnetic`. But `belfem.cpp:234`,
`hphirun.cpp:111` and `hphiTrun.cpp:119` all enter through `solve_coupled()` → `iterate_coupled`,
which assembles at `:1123` (HEAD) and `:1232` (line-search retrial). Neither had it. On that path
`heatloss` was a monotonically growing rank-0 partial across every iteration of every timestep,
and a warm restart resumed from the accumulated value — `load_globals` writes into the existing
variable (`cl_Mesh.cpp:3118-3121`) and the memdump loads after factory construction
(`belfem.cpp:207`).

Grok's contribution here was the one that mattered: it predicted, *before* the fix existed, that
copying the block to the HEAD of `iterate_coupled` would not be enough, because the retrial at
`:1232` re-assembles at a **different** state and would add on top. That is exactly what the first
repair did, and the second one closed it.

Fixed by extracting `Controller::reset_heatloss()` / `collect_heatloss()`
(`cl_FEM_Controller.cpp:2323-2345`) and bracketing all three magnetic assemblies:
`:1122-1124`, `:1231-1233`, `:2149-2151`.

### 2b. `h_newton_mu0` double-counted, and it is the HTS kernel

`dotQ` was added unconditionally, then added again inside the `if ( mx->norm_j( k ) > BELFEM_EPSILON )`
tangent block — a copy-paste into the guard. The four sibling kernels accumulate once.

This is not a corner case: `cl_IWG_Maxwell.cpp:337-341` dispatches `Conductor`/`ThinShell` to
`h_newton_mu0` precisely when μ is **constant**, and HTS carries μ₀. Every integration point with
current reported 2× the Picard value at the same state. Fixed by deleting the second add.

## 3. The P1s

- **`save()` ran before `save_IV()`** in all three executables, and `ExodusWriter::save` snapshots
  the globals inside that call (`cl_Mesh_ExodusWriter.cpp:73-74`, `:661-718`). So the first frame
  carried no `I_*`/`U_*` at all — the Exodus global set was not constant across the `.e-s` series —
  and every later frame carried the previous output step's values. The memdump was unaffected
  (it runs last). Reordered in all three: `belfem.cpp`, `hphirun.cpp:123-124`,
  `hphiTrun.cpp:132-133` and `:193-194`.
- **`I_nn` and `U_nn` were created with the same explicit ID**, `tDof->mesh_basis()->id()`, and
  abstract dofs extracted per node (`cl_FEM_DofMgr_DofData.cpp:4161-4167`) mean several `I_*` can
  share one basis id. `create_global_variable` keys `mGlobalVariableIDMap` on that id
  (`cl_Mesh.cpp:591`), whose uniqueness a 2026-08-15 comment in that very function exists to
  protect. Grok correctly narrowed the blast radius — Exodus and `save_globals` walk the container,
  not the id map, so both values still reached the file; the damage was confined to
  `Mesh::global_variable( id_t )` and to `load_globals` collapsing the map. Fixed by dropping to
  the auto-id form, which is what every other call site in the tree uses.
- **`save_IV` runs a second full element pass into the accumulator — still open, by decision.**
  `:4123 compute_full_matrices()` sits outside the rank-0 guard and reaches the h-kernels through
  `cl_FEM_DofManager.cpp:874` → `cl_IWG_Maxwell.cpp:248` → `mFunMKF`. Christian disputed this;
  the chain was re-verified statically and stands, and a grep of `save_IV`'s body for `heatloss`
  returns nothing. Because the state is converged there, `reset_heatloss()`/`collect_heatloss()`
  around that call would land on the correct value rather than double it — a one-line-pair fix,
  **not applied**. Left as the session's one open defect.

## 4. The regression the repair round introduced, and how it was caught

The "empty filename skips the csv" feature arrived as:

```cpp
tFile = new std::ofstream( aFilename );        // already opens AND truncates
if ( mFirstIVSave )  tFile->open( aFilename, std::ios::trunc );
else                 tFile->open( aFilename, std::ios::app );
```

`basic_filebuf::open` returns a null pointer when `is_open()` is already true, and
`basic_ofstream::open` then calls `setstate( failbit )`. The stream is left open-but-failed, so
every subsequent `<<` is silently dropped — after the constructor already truncated the file.
**`iv_results.csv` became a 0-byte file on every run, both branches.** Probe (libc++, `-std=c++17`):

```
first  call: is_open=1 fail=1
second call: is_open=1 fail=1
--- resulting probe.csv (bytes: 0) ---
```

Worth stating plainly: this is a defect that *reads* correctly. Two `open`s, a truncate branch and
an append branch, a guard on the filename — the logic is right and the stream state is wrong. No
amount of source-trace review would have produced the `bytes: 0`; a fifteen-line probe did, in
about a minute. It is the cheapest gate in this session and it caught the worst defect.

Fixed at `cl_FEM_Controller.cpp:4153-4165`: a stack `std::ofstream` opened exactly **once** with a
ternary on `mFirstIVSave`, `is_open()` as the write guard throughout, and a `BELFEM_ERROR` on a
failed open (file I/O — always-active tier). The `new`/`delete` also leaked on any early return
and is gone.

## 5. csv header aligned to the mesh globals

The header wrote `I_0 V_0 I_1 V_1 …` — 0-based, letter `V`, unpadded — while the globals were
`I_01`/`U_01`. Anyone correlating a csv column with an Exodus global was off by one and looking
for the wrong letter. `tFormatI`/`tFormatU` are now hoisted next to `tAbstractDofs`
(`:4142-4146`) and the header prints through the same `sprint` the globals use (`:4173-4174`), so
one source of truth spells both. `tCount` changed `index_t` → `uint`, which is what the `%0Nu`
that `format_with_leading_zeros` emits actually expects — that also closes the varargs mismatch
the round flagged. Probe over the real `format_with_leading_zeros`/`sprint` bodies:

```
3 dofs   header: time I_1 U_1 I_2 U_2 I_3 U_3          globals: I_1 U_1 I_2 U_2 I_3 U_3
12 dofs  header: time I_01 U_01 … I_12 U_12            globals: I_01 U_01 … I_12 U_12
```

The padding width follows the dof count, as the globals already did.

## 6. Claims that did not survive — mine and the auditors'

- **Mine, retracted:** I asserted `index_t` is `uint64_t`, making the `%u` format UB. It is
  `uint64_t` only under `BELFEM_INT64` (`typedefs.hpp:47-52`), which only `USE_MKL_64BIT_API` sets
  (`CMakeLists.txt:143-145`), and that option **defaults OFF** (`:80`). Codex's narrower framing
  ("UB in ILP64 builds") was the correct one. Real but build-conditional — P3, and now moot.
- **Mine, retracted:** I raised as an open question whether the rank-sum double-counts aura/ghost
  elements. Grok refuted it: `cl_FEM_Block.cpp:63-80` fills `mElements` from
  `aOwnedElementIndices` only, aura goes to a separate container, and the assembly loops walk
  `elements()`.
- **Mine, withdrawn under pushback:** I flagged that `collect_heatloss` omits the `broadcast` its
  in-tree template has (`fn_Mesh_compute_volume.cpp:166`), leaving workers with partials.
  Christian pointed out `reset_heatloss()` zeroes on every rank before the next assembly, so a
  worker partial can never be read as anything else. Not a defect.
- **Mine, wrong in detail:** every `file:line` in my pre-registration was 2-3 lines low, read off
  diff-hunk offsets instead of the file. Both auditors' numbers were right and mine were not.
  The habit to keep: cite from the file, never from the hunk header.
- **Codex, refuted:** "`tDotQ` should be a `Cell`, not a `Vector`" — `fn_Mesh_compute_volume.cpp:151`
  is the in-tree precedent for exactly this reduction and `sum()` is defined only for `Vector`
  (`src/linalg/fn_sum.hpp:46-51`). The *allocation* half of the same finding stands and is
  unfixed: `Vector<real> tDotQ` is constructed on every `collect_heatloss` call.
- **Neither auditor** found the `save_IV`/`compute_full_matrices` pollution (§3) — a single-raiser
  finding, flagged as needing human adjudication and duly disputed.

## Still open

- `compute_full_matrices()` at `:4123` unguarded (§3) — by decision, not oversight.
- Cosmetic, unfixed: `I`, `U`, `Ilabel`, `Ulabel` drop the `t` prefix in framework code (both
  auditors, and the math exemption is scoped to `src/math`/`src/physics`); per-call `Vector` in
  `collect_heatloss`; the `heatloss` duplicate guard in `cl_MaxwellFactory.cpp:1035-1036` carries
  the BC guard's message verbatim ("label the sections to disambiguate"), which cannot be acted on
  for a hard-coded name.
- **Units are recorded nowhere.** `∫ρ|j|²dV` is instantaneous *power* — W in 3D, W/m in a 2D deck —
  and the name "heatloss" reads as an energy. `heatloss`, `I_*`, `U_*` are new user-visible output
  names with no `doc/` entry.
- Whether the coating/side-connector `dV` already carries the shell thickness, and whether coatings
  belong in the integral at all, is a **physics question for Christian** — the round explicitly did
  not vote on it.

## Note for whoever commits

`cl_MaxwellFactory.cpp` was **reindented wholesale** — 737/731 lines changed for a 6-line
substantive edit, and `namespace fem\n{` became `namespace fem {`, matching neither the file nor
the tree. Recommend reverting the reformat and re-applying just the `heatloss` block; as it stands
the diff is unreviewable and will conflict with the parallel session. The four
`examples/*/src/CMakeLists.txt` changes and the `todo/` edits in the same working tree belong to
other work.
