# DR-147: one search order for every plugin `.so`

**Date:** 2026-08-31
**Purpose:** Session record — DR-147 fixed, three by-catch defects closed, and what two audit
rounds changed about the work
**Module:** `src/io`, `src/physics/materials`, `src/numerics/sources`

## 1. The block that was not a block

DR-147 had been parked since 2026-08-29 as "blocked on a layering choice": the source-function
plugin loader in `src/numerics/sources` could not call `material::data_file()`, because numerics
sits below physics. The row framed the choice as resolving at the three call sites versus moving
the resolver up and renaming it, "which touches every caller".

Both branches were unnecessary. `data_file()`'s body depends on exactly two things —
`file_exists()` (`io/filetools.hpp:42`) and `gBelfemDataPath` (`core/globals.hpp:53`) — and
`config/scripts/Add_Library.cmake:4,15` puts `core` and `io` on **every** module's include path
unconditionally. So the search order moved down into `io/filetools` as
`belfem::search_data_file( aFile, aSubDirectory )`, `material::data_file()` became a one-line
forwarder, and **no `src/` CMakeLists changed at all**.

Two of the row's own claims were wrong and are corrected in it:

- it said three plugin `.so` paths "resolve differently". Usermat and defect **both** went through
  `data_file()`. Three keys, **two** behaviours; only the source function bypassed it.
- it priced the move as touching every caller. There were four, all inside `physics/materials`, and
  a forwarder made the count zero.

The call-site branch was also not uniform, which is what killed it: `fem/maxwell` and `fem/thermal`
already carry `physics/materials` on their include path, but `src/circuit/CMakeLists.txt` does not.

## 2. What the audit rounds changed

Four rounds ran: plan audit (Codex `gpt-5.6-terra`/`high`, Grok `grok-4.6`/`high`), then code audit
at `xhigh` for both — §9.1's safety-boundary row, because the diff is centrally ownership and
lifetime.

**The plan round refuted the gate.** The original design was a unit test on `search_data_file`,
pre-registered red-then-green. Both auditors independently pointed out that it cannot be red: the
function does not exist pre-fix, so such a test fails to *compile*. It would only prove the function
was added. Both also caught that the plan's scope guard ("no CMakeLists changes") contradicted its
own test step.

That refutation is what produced **the suite's first MODULE target and first `dlopen`-based test**.
No test in the tree had ever built a `.so` or called `dlopen`; the only `add_library` under `tests/`
was an OBJECT library. Christian ruled for building the fixture (O3(a)) over a cheaper hand-run deck
gate, on the DR-140 precedent that a hand-linked test cannot verify build wiring.

**The plan round also found D3**, which the plan had missed entirely: `Material::read_defect` leaks
its `dlopen` handle exactly as `read_user_defined` did.

**Grok supplied the argument that settled ownership.** The Maxwell factory builds one
`SourceFunction` per group member from the same `( file, label )` — its own comment says so at
`cl_MaxwellBoundaryConditionFactory.cpp:256-259` — so D1 was N unreclaimed mappings, not one.
`dlopen` refcounts, so N guarded closes are exactly right. It also noted that `UserDefinedMaterial`
is the wrong template to mirror: it always `dlopen`s, whereas `SourceFunction` is a plugin for 1 of
8 types, so `mHandle` is usually null and `dlclose( nullptr )` is not a POSIX no-op.

**The code round found the gate was not hermetic.** `Add_Test.cmake:35` sets only `BELFEM_DATA` and
does not clear an inherited loader path, so a contaminated environment could turn the red case
green. Fixed with Grok's named control, `EnvironmentIsNotContaminated`: with an empty root the bare
name must fail to load, so a reachable-by-other-means fixture announces itself instead of silently
passing. Grok also showed the test's own comment was wrong — `dlopen` of a slash-free name searches
`DT_RUNPATH`, `ld.so.cache` and the system directories too, not `$LD_LIBRARY_PATH` alone — and that
`LIBRARY_OUTPUT_DIRECTORY` on the fixture is load-bearing rather than tidiness: without it the
module lands in `${CMAKE_BINARY_DIR}/lib`, which a `USE_SHARED_LIBS=ON` test binary carries on its
RPATH, and the bare name would load without the fix.

**Codex found the docs were wrong on macOS** — six places said `$LD_LIBRARY_PATH` while the plugin
template declares macOS support, where it is `$DYLD_LIBRARY_PATH`.

**Both refuted a stated rationale of mine.** D3's record said `Material`'s copy/move were left alone
because the hazard was pre-existing with a larger blast radius. `cl_Material.hpp:416-419` already
deletes both. The decision was right; the reason given for it was false on its premise, and the row
now says so.

**Codex also found a pre-existing doc defect** unrelated to this change: `input_file_reference.md`
and `physics/materials/doc/README.md` both claimed that an unset `$BELFEM_DATA` leaves only step 1.
False on an installed tree — `cl_Communicator.cpp:79-92` falls back to `BELFEM_INSTALL_DATADIR`
when that directory exists. Fixed in all three places that told the story, including
`fn_material_data_path.hpp`, which Grok caught still telling it after the first two were fixed.

**Split verdicts, both resolved rather than argued.** At plan stage,
`UserLibraryTemplate.cmake`: Grok blocking, Codex non-stale — resolved for Grok, because it is where
a plugin author learns how `file :` is opened. At code stage, Grok said ship, Codex said not yet;
both agreed there was no production correctness blocker, and every overclaim they named was fixed.

Two auditor citations were themselves wrong and were corrected against the tree: Grok's schema line
for the source `file` key (`:1573`, not `:1536`) and its `CMakeLists.txt:435-436`.

## 3. One deliberate deviation from the audited plan

The plan (§3.2) recorded Grok's advice to leave `gMaterialsSubdir` alone, because changing the
constant and the join together lets two errors cancel. The constant was changed anyway, to bare
`"material"`, with both readers joining explicitly. Leaving it meant either `gMaterialsSubdir + 1`
pointer arithmetic at the call — obscure, and silently wrong the day someone drops the slash — or a
second constant, i.e. two sources of truth for one directory name.

The deviation was recorded in the plan **before** the code round and the auditors were pointed at it
as item 1, so it was tested rather than discovered. Codex compared against `HEAD` and confirmed
`data_path()` byte-identical in output and `data_file()` matching in strings and branch order for
every input, including a root with a trailing slash.

## 4. Evidence

All seven changed or added translation units compile clean under the project's **actual** flag set —
`-std=gnu++17 -Wall -Werror -Wno-long-long -pedantic-errors` (`config/compiler/config_gcc.cmake:75`).
`check_doc_claims.py` 37/37; it caught the `[P]`/`[W]` count drift the `[W]`→`[P]` retag caused,
which is what it exists for.

**A process correction worth keeping.** The first pass through this used bare `-fsyntax-only` with
no warning flags and reported clean — and the build then failed, because bare syntax-only ignores
warnings while this build is `-Werror`. The defect was a gtest
`ASSERT_GT( gBelfemDataPath.size(), 0 )`: `size_t` against `int` is `-Werror=sign-compare`. So
"`-fsyntax-only` clean" is a weaker rung on the evidence ladder than it sounds, and was reported
here as if it were stronger. Recorded so the next session checks with the real flags. Note also that
adding `-Wextra` is the wrong correction — the project does not use it, and it fires on pre-existing
headers, which sends you chasing phantom defects in code you did not touch.

**`make check` RAN GREEN** (Christian, 2026-08-31). The new cases were confirmed **present in the
built binaries by name and executed** rather than inferred from the suite result — the DR-117/DR-140
lesson — because this file is wrapped in `#ifdef BELFEM_TEST_PLUGIN_ROOT` and would have passed
*vacuously* had the define not reached it. It did: the define is in
`tests/circuit/CMakeFiles/*/flags.make` and the fixture `.so` is on disk under
`plugindata/material/`.

- `SourcePluginPath` 5/5 in `test_circuit`
- `SearchDataFile` 12/12 in `test_io`

**Both halves of the red-then-green ran in the same binary**, which is stronger than the usual
pre-registration. `EnvironmentIsNotContaminated` sets an empty root, which makes `search_data_file`
hand back the bare name unchanged — the exact pre-fix code path, a raw
`dlopen( "libtest_source_plugin.so" )` — and asserts it FAILS. It does.
`BareNameLoadsFromTheDataDirectory` then loads that same bare name through the data directory and
evaluates the plugin's waveform. So the discrimination is **executed**, not asserted; the control
Codex and Grok asked for turned out to double as the red half.

## 5. Housekeeping note

This session ran in a shared checkout alongside at least two others — a riva-law n-source guard in
`cl_Material.cpp` (a file this work also edits, in separate hunks) and HDF5 fixed-width file types.
Their work was excluded from the audit prompts by name so it would not be reviewed as this change,
and nothing of theirs was touched. `git status` was diffed across both Grok rounds; the writes
observed during them were those sessions, whose devlogs appeared the same night.

Plan, gap table, both audit rounds and the deviation record:
`todo/dr147_plugin_path_resolution.md`.
