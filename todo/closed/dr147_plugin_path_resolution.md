# DR-147: One Search Order for Every Plugin `.so`

**Date:** 2026-08-31
**Purpose:** Three plugin `.so` keys in the deck spell the same `file :` but resolve **two**
different ways; the source-function one bypasses the resolver entirely. Move the search-order
primitive down to `src/io`, where every module can already see it, and route all four load paths
through it — closing the layering block that has kept DR-147 open.
**Module:** `src/io` (primary), `src/physics/materials`, `src/numerics/sources`
**AIs involved:** Claude (exploration + plan), Codex (audit), Grok (third voice)
**Status:** ✅ COMPLETE 2026-08-31 — DR-147 STRUCK and archived (Christian's ruling).

Filed 2026-08-29, unblocked, fixed, audited in four rounds and gate-verified by execution in one
session. **Landed:** `belfem::search_data_file()` in `io/filetools`, with `material::data_file()`
reduced to a forwarder and `SourceFunction::read_user_defined` routed through it — no `src/`
CMakeLists changed. Three by-catch defects closed (D1/D3, both plugin loaders leaked their `dlopen`
handle; D2, an error naming a file the loader never opened) plus D4, a pre-existing doc falsehood
about `$BELFEM_DATA` on an installed tree. Both input-contract artifacts and the plugin-author
template updated.

**Verified by execution**, not by suite result: `SourcePluginPath` 5/5 in `test_circuit`,
`SearchDataFile` 12/12 in `test_io`, each confirmed present in the built binary by name — this file
is `#ifdef`-wrapped and would otherwise have passed vacuously. Both halves of the red-then-green ran
in one binary: `EnvironmentIsNotContaminated` exercises the pre-fix path (empty root → bare name →
`dlopen` fails) and `BareNameLoadsFromTheDataDirectory` the fixed one.

**Residual follow-ups, deliberately not closed with the row:** (i) R15, the Codex prose sweep over
the rewritten `input_file_reference.md` §5/§9 sections; (ii) under a multi-config generator
`LIBRARY_OUTPUT_DIRECTORY` gains a per-config subdirectory and `BELFEM_TEST_PLUGIN_ROOT` would stop
matching the fixture — harmless while the generator is pinned to Unix Makefiles
(`CMakeLists.txt:277-286`), and it bites only on Xcode/VS; (iii) `material::data_path()` still has
zero production callers, so a broken join there would be silent — a pre-existing property this
change did not introduce.

> **Scope guards (Christian's ruling, 2026-08-31 — scope B):**
> - The source plugin searches the **`material` subdirectory**, same as the other two `.so` keys.
> - `gastables::data_path()` and `embed_python_guide()` are **explicitly OUT of scope** — see O1.
> - By-catch defects in the same subsystem ARE in scope (D1–D3), per extend-don't-branch.
> - **No `src/**/CMakeLists.txt` changes and no `config/scripts` include-path changes.** If the
>   implementation needs one, the design is wrong. *Qualified twice:* (i) both auditors flagged the
>   original blanket wording as self-contradictory with the test step, so a SOURCES append in a
>   `tests/**/CMakeLists.txt` is an allowed exception; (ii) **O3 resolved to (a) on 2026-08-31, so
>   building a MODULE `.so` fixture under `tests/` is now explicitly IN scope.** The invariant that
>   still holds without exception: no `src/**/CMakeLists.txt`, no `config/scripts` include-path
>   change.
> - `material::data_file()` keeps its name, signature and header contract. Zero caller churn.

---

## 1. Current Behaviour and How It Fails

| Path | Resolver | Citation |
|---|---|---|
| usermat `.so` | `material::data_file()` at the caller | `cl_MaterialFactory.cpp:558` → `cl_Material_UserDefined.cpp:33` |
| defect `.so` | `material::data_file()` in-function | `cl_Material.cpp:762,764` |
| `bhfile` / HTS `file` (HDF5) | `material::data_file()` | `cl_MaterialFactory.cpp:567,580` |
| **source-function `.so`** | **none — raw deck string to `dlopen`** | `cl_SourceFunction.cpp:221` |

> **Corrected after audit (Codex, confirmed).** This plan's first draft and the DR-147 register row
> both say "three plugin paths resolve three different ways". That is wrong: usermat and defect
> **both** go through `material::data_file()`. There are three plugin **keys** but **two**
> resolution behaviours. The register row has carried this error since it was filed and must be
> corrected with the fix.

| Failure | Mechanism | Evidence |
|---|---|---|
| Same key, different contract | `file :` under `materials` resolves one way; under a `userdefined` source another, and nothing documents the difference | `doc/input_file_reference.md:492-508` vs the source row at `:942`, which is silent |
| **D1** leaked loader handle (source) | `read_user_defined` never stores or `dlclose`s its handle | `cl_SourceFunction.cpp:221` vs `cl_Material_UserDefined.cpp:33,51` |
| **D2** error names a file the loader never opened | `read_defect` reports the resolved path on `dlopen` failure, the unresolved one on `dlsym` failure | `cl_Material.cpp:766` vs `:774` |
| **D3** leaked loader handle (defect) | `Material::read_defect` also keeps `tHandle` local and never closes it | `cl_Material.cpp:762-777` — **found by both auditors; the first draft missed it entirely** |

**Bottom line:** the resolver needs only `file_exists` (`src/io/filetools.hpp:42`) and
`gBelfemDataPath` (`src/core/globals.hpp:53`), so the layering block dissolves the moment the
primitive moves one layer down.

## 2. Architecture: Why the Primitive Goes to `src/io`

`config/scripts/Add_Library.cmake:4,15` puts `src/core` and `src/io` on **every** module's include
path unconditionally; `Add_Executable.cmake:3,15` and `Add_Test.cmake:2,14` carry the same block.
Both auditors independently confirmed no consumer is excluded and no production CMake change is
needed for R1–R4, R6.

**Rejected — resolve at the three `read_user_defined` call sites.** `fem/maxwell` and `fem/thermal`
already carry `physics/materials` (`CMakeLists.txt:35`, `:22`), but **`src/circuit/CMakeLists.txt`
has no physics include** (it has `numerics/sources` at `:29-30`). That branch needs a CMake edit
giving `circuit` a dependency it does not have.

**Rejected — move the resolver up and rename it.** The register row prices this as "renames
`material::data_file` and touches every caller". There are four callers, all inside
`physics/materials`; a forwarder makes the count zero. **That sentence in the register's open
column is now dead and must be struck when the row is updated** (both auditors).

## 3. Gap Table

| # | State / behaviour | Handled today? | Class | Citation / rationale |
|---|---|---|---|---|
| 1 | three-step search order | yes, only inside `physics/materials` | (c) | `fn_material_data_path.cpp:36-74` |
| 2 | subdirectory constant | hardcoded, **and it carries its own separator** | (c) | `:23` is `"/material"`, not `"material"` — see §3.2 |
| 3 | return-unchanged-on-miss | yes | (c) | must survive verbatim; `fn_material_data_path.hpp:36-42` |
| 4 | source `.so` resolution | **no** | (c) | `cl_SourceFunction.cpp:221` |
| 5 | source loader handle lifetime (D1) | **no** | (c) | O2 — RESOLVED |
| 6 | `dlsym` error names resolved file (D2) | **no** | (c) | `cl_Material.cpp:774` |
| 7 | defect loader handle lifetime (D3) | **no** | (c) | `cl_Material.cpp:762-777` |
| 8 | `gastables` / python-guide resolvers | separately, deliberately | (b) | O1 — out of scope |
| 9 | contract in both input artifacts | material only | (c) | `input_file_reference.md:492-508`; schema `:1573` has **no** `resolution:` |
| 10 | plugin-author-facing contract | silent on search order | (c) | `UserLibraryTemplate.cmake:13-17` — Grok blocks on this |
| 11 | directory-hit and empty-string quirks | yes, incidentally | (c) | §3.3 — must be preserved, not "cleaned up" |

### 3.1 The precedence change — corrected after audit

The first draft claimed this makes `$BELFEM_DATA/material` win over `$LD_LIBRARY_PATH` for source
plugins generally. **Grok narrowed that and is right:** POSIX `dlopen` consults `$LD_LIBRARY_PATH`
**only for a name containing no slash**. Every in-tree deck's source `file :` contains a slash —
`examples/disk_pulse/input.conf:101`, `sidecoating:152`, `2D_Undulator:311,320`,
`tape_quench_usermat:167` — so `$LD_LIBRARY_PATH` was never in play for them and step 1 finds them
all from their own run directory.

So the accurate statement is: **the change adds a lookup that did not exist for bare basenames, and
for those it takes precedence over `$LD_LIBRARY_PATH`.** For slashed names it adds two fallbacks
that previously did not exist at all. No in-tree deck changes behaviour.

### 3.2 The join rule must be written, not implied — BLOCKING (Grok C1, Claude C4)

`gMaterialsSubdir` is `"/material"` **with a leading slash** (`fn_material_data_path.cpp:23`), and
`data_path()` concatenates it directly (`:30-31`). Two failure directions follow, and R2's original
wording ("verbatim except `aSubDirectory` replaces `gMaterialsSubdir`") hits one of them:

- pass `"material"` into a body that concatenates directly → `$BELFEM_DATAmaterial` (**Grok**);
- pass `"/material"` into a body that joins with `/` → `$BELFEM_DATA//material`, which POSIX
  collapses so nothing breaks, but which surfaces in every error message naming a resolved path,
  and which no test would catch (**Claude**).

**The rule, pinned:**

```
if ( gBelfemDataPath.empty() )      -> skip steps 2-3 entirely
tBase = aSubDirectory.empty() ? gBelfemDataPath
                              : gBelfemDataPath + "/" + aSubDirectory
tCandidate = tBase + "/" + aFile
```

`data_path()` keeps `gMaterialsSubdir` exactly as it is (R3) — do not "tidy" the constant while
also changing the join, or the two changes will mask each other.

> **DEVIATION FROM THE AUDITED PLAN, 2026-08-31, deliberate and flagged for the code round.**
> The constant WAS changed: `gMaterialsSubdir` is now bare `"material"`, and both `data_path()` and
> `data_file()` join explicitly. Grok's round-1 advice above was to leave it alone precisely because
> changing the constant and the join together lets two errors cancel. I took the other side: leaving
> it would have meant either `gMaterialsSubdir + 1` pointer arithmetic at the call — obscure, and
> silently wrong the day someone drops the slash — or a second constant, i.e. two sources of truth
> for one directory name. The masking risk is covered by
> `SearchDataFile.JoinProducesNoDoubleSeparator` and by the three cases that resolve against real
> shipped files, which fail on any miscombination of constant and join. **The code audit was asked
> to test this judgement specifically rather than to discover it.**

### 3.3 Quirks that must be preserved, not cleaned up (Grok C2)

`file_exists` is `std::filesystem::exists` (`filetools.cpp:23-26`), which is **true for
directories** — it does not test `is_regular_file`. Consequences of today's code that a
"clean-room rewrite" would silently change:

| Input | Today |
|---|---|
| `""` | step 2 is `<data>/material/` — the directory — which exists, so the directory path is returned |
| `"foo/"` missing locally | step 3's basename is `""` → the directory again |
| missing absolute path | step 2 concatenates (POSIX `//` mid-path); step 3 still fires and tries the basename |
| `"name.so"` (no slash) | step 3 is skipped by the `find_last_of('/') != npos` guard (`:59-69`), but step 2 already tried `<data>/material/name.so`, so the guard is a no-op here |

R3's "verify observable behaviour is identical" means **these**. Rejecting a directory hit is a
behaviour change with its own gap row — it is not DR-147.

## 4. Ordered Steps

- [x] **R1** — `src/io/filetools.hpp`: declare `search_data_file( aFile, aSubDirectory )`. Doc block
      states the search order, the join rule from §3.2, the return-unchanged clause and *why* it is
      load-bearing, the O1 note that `gastables` deliberately does not use this, and names
      **`gBelfemDataPath`** — not `$BELFEM_DATA` — as the source of truth (see D4).
- [x] **R2** — `src/io/filetools.cpp` (after: R1): define it, preserving §3.3's quirks and
      implementing §3.2's join rule. Add `#include "globals.hpp"`.
- [x] **R3** — `fn_material_data_path.cpp` (after: R2): `data_file()` → `search_data_file( aFile,
      "material" )`. `data_path()` and the header unchanged.
- [x] **R4** — `cl_SourceFunction.cpp:218-234` (after: R2): resolve before `dlopen`; report the
      **resolved** path in both `BELFEM_ERROR`s. Needs `#include "filetools.hpp"` — the file has
      none today (`:12-20`). No CMake change (§2).
- [x] **R5** — **D1** (after: R4, O2): store `mHandle`, **nullptr-initialised**; `dlclose` in the
      destructor **guarded on non-null**; **delete copy ctor, copy assignment, move ctor and move
      assignment**; decide and document what a second `read_user_defined()` call does (today it
      would silently drop the first handle).
- [x] **R6** — **D2**: `cl_Material.cpp:774` reports `tPath`, not `aLibraryPath`.
- [x] **R7** — **D3** (after: O2): same lifetime treatment for `Material::read_defect`, or an
      explicit gap-table line saying why the defect loader keeps leaking while the source one does
      not. Do not leave it asymmetric and undiscussed.
- [x] **R8** — `doc/input_file_reference.md` §9 (after: R4): the `file`/`label`/`units` row gains
      the resolution rule and §3.1's corrected precedence statement.
- [x] **R9** — `doc/input_schema.yaml:1573` (after: R4): the source `file` key gains
      `resolution: material_data_path`. Widen that block's prose (`:874-886`) to cover all four
      keys. Note the circuit component `file` is documented by reference to the BC keys
      (`:1722-1726`), so `:1573` covers it — say so rather than duplicating.
- [x] **R10** — `UserLibraryTemplate.cmake:13-17` (after: R4): one sentence on the search order.
      **Grok blocks on this**; Codex judged it non-stale. Grok's reasoning wins: it is where a
      plugin author learns how `file :` is opened, and after R4 it is incomplete rather than wrong.
- [x] **R11** — **D4**: `doc/input_file_reference.md:505` and
      `src/physics/materials/doc/README.md:96-97` both say that when `$BELFEM_DATA` is unset only
      step 1 applies. **False for an installed tree** — `cl_Communicator.cpp:79-92` fills
      `gBelfemDataPath` from `BELFEM_INSTALL_DATADIR` in that case. Pre-existing defect, found by
      Codex, verified. Fix both.
- [x] **R12** — Gate built, per O3(a): `tests/circuit/plugin/test_source_plugin.cpp` (MODULE, the
      suite's first), `tests/circuit/test_SourcePluginPath.cpp` (the discriminator),
      `tests/io/test_filetools.cpp` (search-order regression). **Not yet RUN.**
- [x] **R13** — `make check` green (Christian, 2026-08-31), and the new cases confirmed **present
      in the built binaries by name and executed**: `SourcePluginPath` 5/5 in `test_circuit`,
      `SearchDataFile` 12/12 in `test_io`. The `#ifdef BELFEM_TEST_PLUGIN_ROOT` vacuous-pass trap
      was checked explicitly — the define is in `tests/circuit/CMakeFiles/*/flags.make` and the
      fixture `.so` is on disk under `plugindata/material/`.
      All seven changed/added translation units compile clean under the project's real flag set
      (`-std=gnu++17 -Wall -Werror -Wno-long-long -pedantic-errors`), which is still not a
      substitute for running the suite. An earlier pass used bare `-fsyntax-only`, reported clean,
      and the build failed on `-Werror=sign-compare` in `test_filetools.cpp` — bare syntax-only
      ignores warnings. Fixed, and the lesson recorded.
- [x] **R14** — Update the DR-147 register row: strike the dead "renames + touches every caller"
      pricing (§2), correct "three ways" → two behaviours (§1), record what landed.
- [ ] **R15** — Codex language sweep over the touched `input_file_reference.md` section
      (`gpt-5.6-terra`, `medium`).

### 4.1 The gate — first draft REFUTED by both auditors

> **The original §4.1 claimed a `tests/io/test_filetools.cpp` case would run red before the fix and
> green after. That is false and both auditors said so independently.** `search_data_file` does not
> exist pre-fix, so such a test does not fail — it does not *compile*. It tests that the function
> was added, not that DR-147 was fixed. `material::data_file("bhdata.hdf5")` is likewise green
> before *and* after R3, by construction.

Two separate things, and only the second is a discriminator:

**(i) Primitive unit tests — regression value, NOT red-before.** `tests/io/test_filetools.cpp`,
one line appended to `tests/io/CMakeLists.txt` SOURCES (currently only `test_HDF5.cpp`, `:5-7`);
`test_io` is already wired into `check` by the TESTDIR loop (`CMakeLists.txt:433-436`). Cases, with
both auditors' corrections applied:

1. step 1 wins — `aFile` must be a **CWD-relative** name. *(An absolute path makes this case
   vacuous — Grok. My "use absolute paths / Ascii-PWD trap" note was a misapplied memory: that trap
   is `cl_Ascii.cpp:36`'s `getenv("PWD")`, and `file_exists` does not go through Ascii.)*
2. step 2 — relative path below the subdir.
3. step 3 — basename alone.
4. miss — returns the argument unchanged.
5. empty resolver root — **must save, clear and restore `gBelfemDataPath` directly.** It cannot be
   done via the environment: `Add_Test.cmake:35` force-sets `BELFEM_DATA`, and
   `test_io_main.cpp:22` runs `gComm.init()` once before any test. Describe it as "empty
   `gBelfemDataPath`", never as "`$BELFEM_DATA` unset" (see D4).
6. empty `aSubDirectory` — needs an explicit fixture; `share/` root holds only `fluid/`,
   `material/`, `python/`, so there is nothing there to find as written.
7. §3.3 quirk cases: `""`, trailing slash, missing absolute.

**(ii) The actual DR-147 discriminator.** Call `SourceFunction::read_user_defined()` with a
**basename** `.so` present only under `gBelfemDataPath + "/material"`, absent from CWD and from
`$LD_LIBRARY_PATH`, then `compute()`. Red before R4, green after. Belongs in
`tests/circuit/test_SourceFunction.cpp` (which exists, and already holds DR-117's cases), not in
`tests/io`. **This requires a test plugin `.so` fixture, and there is no precedent:** no test in
the tree builds a MODULE library and none calls `dlopen` — the only `add_library` under `tests/` is
the OBJECT library at `tests/physics/backendfree/CMakeLists.txt:26`. That cost is O3.

## 5. Open Design Questions

- **O1 — should `gastables` and `embed_python_guide` adopt the primitive?**
  **RESOLVED 2026-08-31 → no (Christian, scope B).** They are not near-duplicates.
  `gastables::data_path()` is a *directory* resolver with marker-file validation (`gasdata.inp`)
  and a four-deep `../share/fluid` ladder (`fn_GT_data_path.cpp:26,63-66`); `embed_python_guide` is
  a fixed subpath with no search (`fn_embed_python_guide.cpp:30,42-52`). Two further reasons:
  `data_path()` **deliberately skips validation** when `$BELFEM_DATA` is set and says why at
  `:51-55` — a unification that tidies that away is a silent regression that reads like a
  simplification in a diff. And one caller is `nonfree/physics/combustion/main.cpp:48`, in the
  proprietary repo, where the protocol suppresses the exchange — an audit could not see half the
  blast radius. R1's header carries this note so the next duplication sweep does not refile it.

- **O2 — source plugin handle lifetime. RESOLVED 2026-08-31 → (a) store-and-close, hardened.**
  All three of Claude, Codex and Grok converged independently. Evidence:
  1. `_init` installs a pointer into the plugin's text — `set_user_defined()`
     (`cl_SourceFunction.cpp:208-214`), called by every in-tree plugin
     (`examples/disk_pulse/src/bgpulse.cpp:79-81`, `tape_quench_usermat/src/current.cpp:46-49`,
     `src/numerics/sources/example_user_source.cpp:108-111`), invoked via
     `cl_SourceFunction.hpp:92,432-434`. So an *early* `dlclose` is use-after-unload.
  2. Destructor-`dlclose` is nonetheless safe: owners `delete` the object
     (`cl_FEM_PhysicalBoundaryCondition.cpp:29`, `cl_CurrentSource.cpp:36`,
     `cl_VoltageSource.cpp:32`), after which nothing can reach `mUserFunction`.
  3. **The strongest argument for closing rather than leaking is Grok's:** the Maxwell factory
     builds **one `SourceFunction` per group member** from the same `(file, label)` — stated in the
     code's own comment at `cl_MaxwellBoundaryConditionFactory.cpp:256-259`, constructed at `:364`.
     That is N `dlopen`s of one `.so`. POSIX `dlopen` is refcounted, so N store-and-close pairs are
     correct, while a deliberate leak is N unreclaimed mappings.
  4. **`UserDefinedMaterial` is the wrong template to mirror.** It always `dlopen`s in its only
     constructor and `dlclose`s unconditionally (`cl_Material_UserDefined.cpp:33,49-52`).
     `SourceFunction` is a plugin for 1 of 8 types, so `mHandle` is usually null — and
     `dlclose(nullptr)` is **not** a POSIX no-op. That class also has a destructor without deleted
     copy/move; do not replicate its Rule-of-Three hole.
  5. `SourceFunction` is copyable today (`~SourceFunction() = default`, `cl_SourceFunction.hpp:101`;
     `mValues` is a copyable `Cell`). Nothing copies it in production — all six sites are
     `new SourceFunction()` — so the hazard is **latent**, which is exactly DR-141's shape.
     Deleting copy/move costs nothing and is what makes (a) safe.

- **O3 — is the discriminating gate worth building the suite's first `.so` test fixture?**
  **RESOLVED 2026-08-31 → (a), build the fixture (Christian.)** A MODULE library fixture under
  `tests/`, giving a genuine red→green gate on R4 inside `make check`, permanently. Rejected:
  **(b)** a hand-run deck gate — cheaper, and genuinely discriminating, but DR-140 records that a
  hand-linked scratch test cannot verify build wiring, and the project has repeatedly paid for
  treating a hand gate as equivalent; **(c)** primitive tests only with R4 recorded
  reviewed-not-verified — the register already carries several rows in that state and does not need
  another. The fixture is the suite's first MODULE target and first `dlopen`-based test; it is
  reusable for the material and defect plugin paths, which are equally ungated today.

## 6. The Interface

```cpp
// src/io/filetools.hpp
/**
 * Resolve a data file named in an input file, below aSubDirectory of the
 * shared data directory.
 *
 * The run directory always wins, so a local copy overrides the shared
 * database. Failing that the file is looked up below the data directory,
 * first under the same relative path and then by name alone.
 *
 * If nothing is found, aFile is returned UNCHANGED. This is load-bearing:
 * it makes the opener report the name the user wrote, and it keeps the
 * dlopen search path intact for plugin libraries, which resolve through
 * $LD_LIBRARY_PATH and need not exist as a file at all.
 *
 * The root is gBelfemDataPath, NOT $BELFEM_DATA directly: Communicator::
 * set_globals() fills it from the environment, and on an installed tree
 * from BELFEM_INSTALL_DATADIR when the environment is silent. An empty
 * gBelfemDataPath skips the two fallbacks. An empty aSubDirectory searches
 * the data directory itself.
 *
 * Note: gastables::data_path() deliberately does NOT use this. It is a
 * directory resolver with marker-file validation and a relative-path
 * ladder, and it intentionally skips validation when the root is set.
 */
std::string
search_data_file( const std::string & aFile, const std::string & aSubDirectory );
```

## 7. Definition-of-Done Checklist

- [x] Every gap-table row mapped to a step or an open question.
- [x] **O3 decided before R12.**
- [x] §3.2's join rule implemented exactly; §3.3's quirks preserved.
- [x] D1–D4 each either fixed or carrying an explicit written reason not to be.
- [x] `make check` green; new cases confirmed present by symbol **and executed** — 17/17.
- [x] Both input-contract artifacts updated in this session, plus `UserLibraryTemplate.cmake`.
- [x] `check_doc_claims.py` clean. **37/37** — it caught the `[P]`/`[W]` count drift from the retag.
- [x] No `src/**/CMakeLists.txt` and no `config/scripts` changes (test SOURCES append excepted).
- [x] `git status` reviewed after every Grok round.
- [x] DR-147 row corrected (two behaviours, not three ways) and updated with what landed.

## 8. Audit Trail

- Exchange thread: `tmp/ai_exchange/dr147_plugin_path_resolution.md`
- **Plan audit round 1, 2026-08-31.** Codex `gpt-5.6-terra`/`high`; Grok `grok-4.6`/`high`.
  Both returned "do not implement as written". Convergent blocking findings: the R9 gate does not
  test the fix; case 5 is untestable via the environment; the scope guard contradicts R9; O2 is
  store-and-close with copy/move deleted; `Material::read_defect` leaks its handle too (D3).
  Grok alone: the join rule (C1), the directory-hit quirks (§3.3), the `dlopen`-slash narrowing of
  §3.1, the per-group-member `dlopen` argument, `UserLibraryTemplate.cmake`. Codex alone: the
  installed-tree `BELFEM_DATA` doc defect (D4), the "three ways" factual correction.
  Claude pre-audit: the join-rule trap from the other direction, the O2 evidence chain.
  **Split verdict:** `UserLibraryTemplate.cmake` — Grok blocking, Codex non-stale; resolved for
  Grok. `share/material/README.md` — Claude said stale, both auditors said accurate; resolved for
  the auditors, no edit.
- Every auditor citation was re-read against the tree before inclusion. Two were wrong and are
  corrected here: Grok's schema line for the source `file` key is **`:1573`**, not `:1536`; its
  `CMakeLists.txt:435-436` is the `add_dependencies` pair, with the TESTDIR list at `:433-434`.

- **Code audit round, 2026-08-31.** Codex `gpt-5.6-terra`/`xhigh`; Grok `grok-4.6`/`xhigh` — §9.1's
  safety-boundary row, chosen because the diff is centrally about ownership and lifetime.
  **Split verdict: Grok "I would ship this"; Codex "do not ship as fully validated/documented yet".**
  Both agreed there is no production correctness blocker; they differ on whether two overclaims are
  ship-stops. Everything below was fixed rather than argued.

  Convergent findings:
  - **The gate was not hermetic.** `Add_Test.cmake:35` sets only `BELFEM_DATA` and does not clear an
    inherited loader path, so a contaminated environment could make the "red before" case green.
    Fixed by adding `SourcePluginPath.EnvironmentIsNotContaminated` — Grok's named control: with an
    empty root, the bare name must FAIL to load. If it loads, the fixture is reachable some other
    way and the discriminating case is announced as meaningless instead of silently passing.
  - **`SearchDataFile.MaterialDataFileContractIsUnchanged` never called `material::data_file()`.**
    Renamed to `SearchOrderBehindTheMaterialForwarder`, with the limit stated in place: `tests/io`
    cannot include `physics/materials`, so the one-line forwarder is covered by inspection only.
  - **Both refuted D3's stated rationale.** `Material` copy/move were **already deleted at
    `cl_Material.hpp:416-419`**, so the recorded reason for leaving them alone — "pre-existing
    hazard, larger blast radius" — was false on its premise. The decision stands; the reasoning
    given for it did not, and is corrected in the row.

  Grok alone:
  - `dlopen` of a slash-free name searches `DT_RPATH`/`DT_RUNPATH`, `ld.so.cache` and the system
    directories too, **not** `$LD_LIBRARY_PATH` alone as the test comment claimed. Comment corrected.
  - `LIBRARY_OUTPUT_DIRECTORY` on the fixture is **load-bearing for discrimination**, not tidiness:
    without it the module lands in `${CMAKE_BINARY_DIR}/lib`, which a `USE_SHARED_LIBS=ON` test
    binary carries on its RPATH — and the bare name would then load without the fix. Now commented.
  - `fn_material_data_path.hpp:23-25` still told the pre-R11 `$BELFEM_DATA` story. Fixed.
  - Quirk coverage was half-landed. Added `TrailingSlashMissResolvesToTheDirectory` and
    `MissingAbsolutePathStillFallsBackToTheName`.
  - The header comment claiming "nothing in the tree copies a SourceFunction — every site is
    `new SourceFunction()`" is false: the tests stack-allocate. Corrected.
  - Open risk, NOT fixed: under a multi-config generator `LIBRARY_OUTPUT_DIRECTORY` gains a
    per-config subdirectory and `BELFEM_TEST_PLUGIN_ROOT` would not match the artifact. The
    generator is pinned to Unix Makefiles (`CMakeLists.txt:277-286`), so this bites only if someone
    ships on Xcode/VS.

  Codex alone:
  - **The loader path was documented platform-neutrally while the template declares macOS support.**
    macOS uses `$DYLD_LIBRARY_PATH`. Corrected in six places.
  - The `<build>/plugindata/material/` comments were stale — the root is
    `${CMAKE_CURRENT_BINARY_DIR}`, so it is `<build>/tests/circuit/plugindata/material/`. Corrected.
  - Confirmed the §3.2 deviation against `HEAD`: `data_path()` byte-identical in output, `data_file()`
    matching in strings and branch order for every input including a trailing-slash root.

  Residual, stated rather than closed: Grok notes `material::data_path()` has **zero production
  callers**, so a broken join there would be silent. That is a pre-existing property of that
  function, not something this change introduced.
