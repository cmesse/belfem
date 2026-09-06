# The Installed Public API: Replace the Whole-Tree Header Install With a Flat User API Set

> **DEFERRED 2026-09-04** (Christian's ruling): the plan was split. The cheap, layout-independent
> items landed the same day — `LICENSE` installs to the prefix root (R7/D4), the umbrella and
> `cl_SourceFunction.hpp` guards are `BELFEM_`-prefixed (D1 half, D9), the umbrella carries the
> prohibition note and the `gTbulk` sentence (R2 half), the shipped compile line is correct for the
> module layout (D10), and D8 was already gone from the tree. The **spine** — the flat
> `install( FILES )` list, the umbrella move, the isolated staged gate, the template rewrite and the
> end-to-end gate (R1, R2 move, R3, R4, R6, R8) — is parked: the 2026-08-31 repairs made an installed
> prefix build plugins, which was the pain this plan answered, and R1 would now have to rewrite the
> include-path blocks in four decks and two templates that those repairs just fixed. O1 and O5 stay
> unanswered. Revive only if the "application plus plugin SDK" boundary is wanted for a release;
> the revision-2 re-audit is still owed before any spine code is written.

**Date:** 2026-08-29
**Purpose:** `make install` currently deploys every `.hpp` under `src/` — 629 headers in module
layout — of which the documented plugin API is 14. Replace that glob with an explicit, flat
`include/belfem` containing exactly the user API plus the `belfem_user_api.hpp` umbrella, so a
plugin author needs one `#include` and one `-I`. The mechanism is an `install( FILES )` list pinned
by a compile gate that builds the shipped examples against a staged flat copy **with the inherited
include path removed**.
**Module:** `CMakeLists.txt`, `src/fem/kernel` (umbrella move), `src/physics/materials`
(templates), `tests/physics/backendfree`, plus the doc set in R6
**AIs involved:** Claude (exploration + plan), Codex + Grok (plan audit round 1, 2026-08-29)
**Status:** PLAN — **revision 2**, 2026-08-29. Revision 1 failed its audit on three P0s, all in the
step specification, none in the architecture; see §8. Revision 2 rewrites R1/R3/R4 and adds R6.
**Re-audit owed before any code is written.** Decision to replace rather than supplement taken by
Christian, 2026-08-29.
**Update 2026-09-04:** split by Christian — see the banner. R7 and the layout-independent halves
of R2 landed; the spine is deferred.
**Update 2026-08-31:** R5, D2 and D3 landed independently of the replan — they repair the
*existing* module-layout install rather than anticipating the flat one, so they neither
depend on nor prejudge R1. An installed prefix now builds plugins. R1–R4 and R6–R8 are
unchanged and still owed; R1 will have to rewrite the include-path block in all four decks
and both templates when the layout goes flat.

> **Scope guards:**
> - **Out of scope:** the `belfemConfig.cmake` / `find_package` question (open as **O4** in
>   `todo/closed/shared_library_and_install_plan.md:334-337`); RPATH; static-vs-shared `libbelfem`; any
>   change to solver or physics behaviour; Windows.
> - **Compatibility promise explicitly dropped (decided 2026-08-29, Christian):** an installed
>   BELFEM can no longer compile a program that *links* `libbelfem` as a general C++ library. The
>   install becomes an **application plus a plugin SDK**: 7 executables by default, one library,
>   and 16 header files. Anyone embedding BELFEM in their own solver needs the source tree. This is
>   policy, not a provable claim — Grok's falsification attempt found no in-tree consumer that
>   contradicts it, but an out-of-tree embedder cannot be ruled out by grep.
> - **In scope:** the header install rule, the umbrella's location/guard/licence, both plugin
>   templates, the backend-free compile gate, the installed Undulator example, and the four
>   documents that go false with R1.

---

## 1. Current Behaviour and How It Fails

One rule deploys headers (`CMakeLists.txt:341-346`):

```cmake
install( DIRECTORY ${BELFEM_SOURCE_DIR}/
         DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${LIBPREFIX}
         FILES_MATCHING PATTERN "*.hpp"
         PATTERN "doc" EXCLUDE PATTERN "CMakeFiles" EXCLUDE )
```

`BELFEM_SOURCE_DIR` is `src` (`:143`), `LIBPREFIX` is `belfem` (`config/globals.cmake:5`).
Measured 2026-08-29: **629 headers** across 14 module trees.

| Failure | Mechanism | Evidence |
|---|---|---|
| No API boundary — every header is a de-facto promise | directory glob; no public/private marker exists in the tree | `CMakeLists.txt:341-346`; no `PUBLIC_HEADER` property, no `BELFEM_API` macro |
| `include/belfem` alone resolves nothing | headers include each other by bare name; in-tree this works only because `config/scripts/Add_Library.cmake:3-25` puts every module dir on the path | consumer needs 6–13 `-I` entries |
| Headers ship for modules absent from the library | the glob has no option gate | `src/CMakeLists.txt:20-22` (`USE_VTK` OFF, `CMakeLists.txt:104`); `src/physics/CMakeLists.txt:11-14` (`USE_GASMODELS` OFF, `:97`) |
| Those headers are unusable as a public API | `cl_VTK_Curve.hpp:15-18` includes `cl_Matrix.hpp`; `cl_VTK_MeshView.hpp:15-17` includes `cl_Mesh.hpp` — neither is in the API set | *(Revision 1 claimed `vtktypes.hpp` had unguarded VTK includes. **False** — they sit inside `#ifdef BELFEM_VTK`, `vtktypes.hpp:15-30`, and `BELFEM_VTK` is set only under `USE_VTK`, `CMakeLists.txt:209-212`. Retracted; the row survives on the mechanism at left.)* |
| The include-dir list is triplicated by hand | `UserMaterialTemplate.cmake:78`, `UserLibraryTemplate.cmake:99-108`, `examples/2D_Undulator/lib/CMakeLists.txt:47-56` | which is why the third copy rotted — D2 |
| Neither working template is installed | the glob matches `*.hpp` only; `src/physics/materials/CMakeLists.txt` has no `install()` | `grep -rn "UserMaterialTemplate\|UserLibraryTemplate" --include=CMakeLists.txt --include=*.cmake` → one comment hit |

**Bottom line:** the install ships 629 headers to serve an API of 14, ships no artifact telling a
user how to consume them, and the one plugin `CMakeLists.txt` it does ship cannot compile against
it. The count is not the problem; the absence of a boundary is — and 14 is now a defensible
boundary because the three entry headers were made backend-free on 2026-08-29.

## 2. Architecture: Why a Flat Explicit List Is Now the Right Spine

A curated set was rejected in the 2026-08-29 jury round
(`tmp/ai_exchange/review_install_headers.md`) on two grounds, **both removed by Christian's API
change the same day**:

1. *"It breaks the moment a plugin uses `Vector`/`Matrix`."* — The supported API now excludes it.
   All three entry headers carry the same prohibition: `cl_Material.hpp:24`,
   `cl_Material_UserDefined.hpp:17`, `cl_SourceFunction.hpp:43` — *"NEVER include cl_Vector.hpp or
   cl_Matrix.hpp here, or any class that uses it."* `cl_Material_UserDefined.hpp` was decoupled the
   same day: closure **57 headers across 10 module dirs → 11 across 3**, no backend define.
   The notes are a developer contract and are treated as one — R3 does not police them. What R3
   proves is narrower and mechanical: that the list in §6 compiles as shipped, standalone.
2. *"Flattening collides."* — Not for this set. No basename among the 14 appears twice in `src/`.
   The known collision — `bernstein/cl_IFB_LINE3.hpp:13` and `lagrange/cl_IFB_LINE3.hpp:13`, same
   basename **and** same guard `BELFEM_CL_IFB_LINE3_HPP` — is in `fem/interpolation`, outside the
   set.

**Flat, not module layout.** The point is that `#include "cl_Material.hpp"` works with one `-I`.

**Explicit `install( FILES )`, not a glob.** A glob cannot express "these 14"; the list *is* the
boundary. R3 makes it mechanically checked.

**A shipped boundary, not access control.** All 16 files sit in one public directory, so nothing
stops a user including `globals.hpp` directly. Only `belfem_user_api.hpp` is the supported entry
point; the rest are shipped dependencies. §6 says so; CMake cannot enforce it.

**Rejected — ship both.** Keeps every header a promise, which is the failure being fixed. Decided
against by Christian 2026-08-29.

**Rejected — module layout for the 14.** Costs the single `-I` and buys nothing.

## 3. Gap Table

Measured with `g++ -H` / `-fsyntax-only`, gnu++17, against a staged **flat** directory containing
only the listed files. All gates rc=0 under `-DDEBUG` and `-DNDEBUG`, with **no** backend define
and **no third-party include directory**. *(Revision 1 said "TPL headers where needed" — that
clause belonged to the pre-decoupling 57-header measurement and contradicted §6. Removed.)*

| # | State | Needed for | Class | Citation / rationale |
|---|---|---|---|---|
| 1 | `core/typedefs.hpp`, `constants.hpp`, `assert.hpp`, `fn_sprint.hpp`, `globals.hpp` | every entry point | (c) explicit | `globals.hpp` and `powerlaws.hpp` also arrive transitively via `cl_Material.hpp:14` and `:2279`, and were umbrella-named 2026-08-29 (O2). Omitting either from the install list fails: `fatal error: globals.hpp: No such file or directory`. `powerlaws.hpp` needs `gTbulk` (`globals.hpp:44`) |
| 2 | `containers/cl_Bitset.hpp`, `cl_Cell.hpp` | `Material` members | (c) explicit | `cl_Material.hpp:21-22` |
| 3 | `physics/materials/cl_Material.hpp`, `powerlaws.hpp`, `cl_JcFunction.hpp` | materials + defects | (c) explicit | closure of `cl_Material.hpp` = **10** |
| 4 | `cl_Material_UserDefined.hpp`, `fn_Material_UserDefinedPolynomials.hpp`, `cl_JcFunction_UserDefined.hpp` | plugins reaching further | (c) explicit | each closes at **11** = the 10 plus itself; no fan-out |
| 5 | `numerics/sources/cl_SourceFunction.hpp` | source / current functions | (c) explicit | closure 5, four already above; callback is `typedef real ( UserFunc )( const real )` at `:78` |
| 6 | `belfem_user_api.hpp` umbrella | one-`#include` entry | (c) explicit | exists, wrong place — D1 |
| 7 | `belfem_version.hpp` | version inspection | (a) keep | `CMakeLists.txt:349-350`; already lands at the include root |
| 8 | Plugin CMake templates | telling a user how to build | (c) explicit | D3 |
| 9 | `LICENSE` | every installed header points to it | (c) explicit | D4 |
| 10 | `visualizer`, `gastables`, `gasmodels` headers | nothing in the user API | (a) dissolved | the replacement ships none of them |
| 11 | Programs linking `libbelfem` | not a supported install use case | (b) → **decided** | scope guards |

### 3.1 Cross-cutting findings

**The gate must remove the inherited include path, or it proves nothing.** This is the single most
important correctness item, and revision 1 got it wrong. `tests/physics/backendfree/` is
`add_subdirectory`'d from `tests/physics/CMakeLists.txt`, which has already run
`include_directories` for `math/tensor`, `math/tools`, `numerics/bezier`, `numerics/integration`,
`numerics/spline`, `physics`, `physics/database`, `physics/materials`, `mesh` (`:3-12`);
`config/scripts/Add_Test.cmake:3-19` adds `core`, `comm`, `containers`, `linalg`, `linalg/lapack`,
`linalg/{armadillo|blaze}`, `linalg/operators`, `io`, `math/graph`, `sparse`, `numerics/spline`,
`numerics/bezier`, `generated`. **`include_directories()` in a child directory appends; it does not
replace.** Staging a flat copy while that inheritance stands is a no-op.

**The prohibition covers 3 of the 14 headers.** `cl_Material.hpp:24`,
`cl_Material_UserDefined.hpp:17` and `cl_SourceFunction.hpp:43` carry the *"NEVER include
cl_Vector.hpp"* note. The other eleven — including `cl_Cell.hpp`, `powerlaws.hpp`,
`cl_JcFunction.hpp` and `globals.hpp` — carry nothing. A **direct** violation of that note is a
developer error and is not the gate's business (Christian, 2026-08-29). A **transitive** one is
different: an include added to one of the eleven pulls linalg in two hops past a note its author
never saw. That is what `cl_Material_UserDefined.hpp` was doing until 2026-08-29 — `cl_Vector.hpp`,
`cl_Database.hpp`, `fn_polyval.hpp`, 57-header closure — by accretion, not by defiance.

**Guard 2 is narrower than it looks.** `tests/physics/backendfree/test_MaterialBackendFree.cpp:34-39`
`#error`s on `BELFEM_CL_VECTOR_HPP`, `BELFEM_CL_MATRIX_HPP`, `BELFEM_SPLINE_HPP`,
`BELFEM_CL_SPMATRIX_HPP` and the four wrapper guards — so it catches the backend and spline
families but not `fn_polyval.hpp` or `cl_Database.hpp`. Worth knowing when reading a green gate;
not worth extending speculatively.

**In-tree and installed layouts differ permanently.** In-tree keeps the module `-I` list
(`Add_Library.cmake:3-25`); installed is flat. The templates must therefore keep two branches
forever. R3 tests only the installed shape; the source-tree branch is exercised only by R8.

**The umbrella's include list IS the ship list** (since 2026-08-29, Christian — O2). All 14 content
headers are named explicitly in `belfem_user_api.hpp:7-20`, so the R1 CMake list and the umbrella
can be cross-checked mechanically against each other; a disagreement is the drift this plan is
guarding against, and catching it does not require R3 to run. 14 named + umbrella + version header
= **16 installed files**.
*(Before that change the umbrella named 12 and relied on `globals.hpp` and `powerlaws.hpp` arriving
transitively through `cl_Material.hpp:14` and `:2279` — anyone deriving the install list from it
would have shipped a prefix that does not compile.)*

**Every header in the set is self-contained.** Measured 2026-08-29: each of the 14 compiles alone
against the flat directory, no backend define, under `-DDEBUG` and `-DNDEBUG`. So the umbrella's
include order carries no meaning, and a user who includes one header directly instead of the
umbrella still compiles. Not true of the tree in general — `cl_AR_Matrix.hpp:34` uses `arma::Mat<T>`
without including `armadillo.hpp`, relying on `cl_Matrix.hpp:15` having pulled `cl_Vector.hpp`
first. Worth a regression thought if the set ever grows.

## 4. Ordered Steps

- [ ] **R1 — Flat `install( FILES )` *and* both template probes, in one change.**
      Define the 14 paths once as a CMake list (visible to `tests/`, which is added at
      `CMakeLists.txt:409`, after the install block). Delete `CMakeLists.txt:341-346`; install the
      list plus `belfem_user_api.hpp` flat into `${CMAKE_INSTALL_INCLUDEDIR}/${LIBPREFIX}`; keep
      `:349-350` unchanged. **In the same commit**, rewrite the *install branch* of both templates
      to one `-I ${BELFEM_DIR}/include/belfem` probing a flat file, keeping the module `-I` list on
      the source-tree branch. Revision 1 split this and would have shipped a prefix whose own
      templates abort in `UserMaterialTemplate.cmake:86-91`.
- [ ] **R2 — Umbrella: move, guard, licence, and the strays it should take with it**
      *(after: R1)*. **Everything but the move landed 2026-09-04** (guard, licence block, note,
      `gTbulk` sentence, D8, D9); the move to the include root stays with R1. Move to the include root; guard `BELFEM_BELFEM_USER_API_HPP` →
      `BELFEM_USER_API_HPP`; add the LBNL/UC block (it has none — `:1-3` is an IDE stub); add the
      *"NEVER include cl_Vector.hpp"* note its three wrapped headers carry; add one sentence that
      `gTbulk` resolves from the host at load. Also: delete the dead
      `template< typename T > class Vector;` at `cl_Material.hpp:85` (**zero** uses of `Vector<` in
      that file), and rename `cl_SourceFunction.hpp:12`'s guard `CL_FUNCTION_HPP` to a
      `BELFEM_`-prefixed one — unprefixed and collision-prone in a flat public directory.
- [ ] **R3 — Rebuild the backend-free gate as a genuinely isolated staged compile** *(after: R1)*.
      Three requirements, all load-bearing:
      1. **Isolate.** Overwrite the directory/target `INCLUDE_DIRECTORIES` property so the object
         library sees *only* the staged directory. Appending is not enough (§3.1).
      2. **Stage with dependency tracking.** `add_custom_command( OUTPUT … DEPENDS <source header> )`
         per file, driven by the same list R1 installs. Configure-time `file( COPY )` goes stale on
         a header edit until CMake re-runs. Delete stale stage contents so a *removed* file cannot
         mask an omission.
      3. **Cover all three examples.** `example_user_material.cpp` (today's only subject),
         `example_user_defect.cpp` and `example_user_source.cpp` — the latter two have **never been
         compiled by any gate** — plus a TU including only `belfem_user_api.hpp`. Keep the existing
         `remove_definitions` (`:13-14`) and Guard 2, and keep the `check-fast` hook
         (`CMakeLists.txt:446-451`).
- [ ] **R4 — Rewrite `UserLibraryTemplate.cmake`'s backend machinery out** *(after: R3 proves it
      unnecessary)*. A **rewrite, not a deletion list**: revision 1 listed `:99-108` and `:111-121`
      but not `:125-128`, which still reads `BELFEM_HEADER_ROOT` and iterates `BELFEM_MODULE_DIRS`
      — both would be unset. Remove `BELFEM_BACKEND` (`:55`), `BELFEM_TPL_INCLUDE_DIRS` (`:65-66`),
      the validation (`:90-93`), the linalg/sparse module dirs (`:99-108`), the `cl_Vector.hpp`
      probe (`:111-121`), the TPL append (`:130-132`), `target_compile_definitions` (`:175`) and
      the summary line (`:209`). Also drop the `math/graph` + `io` "headroom" from
      `UserMaterialTemplate.cmake:78,94-95` — those subdirectories will not exist on the prefix.
      See **O1** on merging the two.
- [x] **R5 — Install the templates and the three `example_user_*.cpp`** ~~*(after: R4)*~~ into a
      **subdirectory** — `${CMAKE_INSTALL_DATADIR}/${LIBPREFIX}/templates/` — not the `share/` root,
      which `CMakeLists.txt:359-361` already occupies. Fixes D3.
      **Done 2026-08-31, ahead of R4** — the dependency on R4 was ordering convenience, not a
      requirement: the templates already configure and build against an installed prefix as they
      stand, and shipping them is what makes the prefix usable today. `install( FILES )` after the
      `examples/` rule in `CMakeLists.txt`, anchored on `BELFEM_SOURCE_DIR`. R4's backend cleanup
      is unaffected and still owed.
- [ ] **R6 — Update the documents R1 makes false.** `CLAUDE.md:409` ("`make install` deploys …
      headers …"), which `scripts/check_doc_claims.py` guards; a `HISTORICAL` note on the
      module-layout paragraph of `todo/closed/shared_library_and_install_plan.md`; the doxygen line
      `g++ … -I/path/to/belfem/include` at `cl_Material_UserDefined.hpp:120` — wrong today, and it
      ships as one of the 14; and `src/physics/materials/doc/materials_usage_guide.md:334-373`,
      which mixes plugin-side `#include "cl_Material.hpp"` with host-side `MaterialFactory` (not in
      the set). Not optional polish.
- [x] **R7 — Install `LICENSE`** to `${CMAKE_INSTALL_PREFIX}/LICENSE`, matching the "top-level
      LICENSE file" wording the headers use. Independent of the others. Fixes D4.
      **Done 2026-09-04**: `install( FILES ${CMAKE_SOURCE_DIR}/LICENSE DESTINATION . )` after the
      templates rule in `CMakeLists.txt`.
- [ ] **R8 — End-to-end gate (the definition of done).** Real `make install DESTDIR=<tmp>`, then
      build all four shipped plugin examples against the prefix **using the installed template**,
      on Linux and Darwin. Specify precisely, or it is not a gate (Codex): `BELFEM_DIR` is
      `${DESTDIR}${CMAKE_INSTALL_PREFIX}`, not `DESTDIR`; the templates hardcode `myalloy.cpp` /
      `defect.cpp`+`current.cpp` (`UserMaterialTemplate.cmake:50`, `UserLibraryTemplate.cmake:68`)
      and must be overridden to the shipped `example_user_*.cpp` names — a gate whose operator
      hand-edits installed files is not a gate. Also exercise the **source-tree** branch, which R3
      does not cover. Depth: see **O5**.

## 5. Open Design Questions (not silently decided)

- **O1 — Merge the templates into one `UserPluginTemplate.cmake`?** After R4 they still differ in
  variable names (`MATERIAL_*` vs `LIBRARY_*`), comments and default sources — revision 1's "differ
  only in the source list" was wrong. Merging is reasonable and not free. **Recommend merge**;
  three copies of the include list is what rotted the Undulator example. Christian's call.
- **O2 — RESOLVED 2026-08-29, Christian → both named.** `belfem_user_api.hpp:9,20` now include
  `globals.hpp` and `powerlaws.hpp` explicitly, so the umbrella's 14 `#include` lines are exactly
  the R1 install list. Re-verified after the change: umbrella alone and all four shipped examples
  compile against the flat set, one `-I`, no backend define, `-DDEBUG` and `-DNDEBUG`. Buys more
  than readability — the two lists are now mechanically comparable (§3.1).
- **O3 — RESOLVED 2026-08-29, Christian → moved.** `git mv src/physics/materials/example_user_source.cpp
  src/numerics/sources/example_user_source.cpp`, beside the `cl_SourceFunction.hpp` it demonstrates.
  No build impact: the file appears in no `CMakeLists.txt`. R3 and R5 must use the new path.
- **O4 — RESOLVED 2026-08-29 (Grok falsification attempt).** In-tree compile-time consumers of an
  installed prefix are exactly: the two templates, the Undulator example, and the backend-free gate
  (which uses the build tree, not the prefix). Runtime loaders are exactly three `dlopen` sites
  (`cl_Material.cpp:764`, `cl_Material_UserDefined.cpp:33`, `cl_SourceFunction.cpp:224`). Host code
  compiles against the source tree. No hidden consumer contradicts the scope guard; an out-of-tree
  embedder remains unprovable either way.
- **O5 — How deep is R8?** Compile-only, or `dlopen` a built plugin with the installed `belfem` /
  `material` binary? Grok argues for load-and-run, citing the stale-plugin vtable failure as
  precedent: a compile gate cannot catch an ABI break, and this plan publishes an SDK. Against:
  materially heavier, and DR-131 records that the Darwin `dynamic_lookup` path has never been
  compiled on a Mac at all, so Darwin load-and-run may not be reachable this release.
  **Raised for Christian; unanswered.**

## 6. The Installed Interface

`${CMAKE_INSTALL_PREFIX}/include/belfem/` — flat, **16 files**:

```
belfem_user_api.hpp          <- umbrella; the only SUPPORTED entry point
belfem_version.hpp           <- existing rule, unchanged

typedefs.hpp  constants.hpp  assert.hpp  fn_sprint.hpp  globals.hpp
cl_Bitset.hpp  cl_Cell.hpp
cl_Material.hpp  powerlaws.hpp  cl_JcFunction.hpp
cl_Material_UserDefined.hpp  fn_Material_UserDefinedPolynomials.hpp
cl_JcFunction_UserDefined.hpp
cl_SourceFunction.hpp
```

Source paths for the R1 list: `src/core/{typedefs,constants,assert,fn_sprint,globals}.hpp`,
`src/containers/{cl_Bitset,cl_Cell}.hpp`,
`src/physics/materials/{cl_Material,powerlaws,cl_JcFunction,cl_Material_UserDefined,fn_Material_UserDefinedPolynomials,cl_JcFunction_UserDefined}.hpp`,
`src/numerics/sources/cl_SourceFunction.hpp`.

### 6.1 The contract a plugin must satisfy

Backend-free is **not** host-configuration-free (Codex). Every macro below changes what the shipped
headers compile to, and the install exports none of them today:

| Macro | Effect | Source |
|---|---|---|
| `NDEBUG` / `DEBUG` | sets `BELFEM_ASSERTIONS_ACTIVE`; changes the **bodies** of inline `Material` and `Cell` methods, not layout | `src/core/assert.hpp:24-32` |
| `BELFEM_INT64` | switches `int_t` / `index_t` between 32 and 64 bit | `src/core/typedefs.hpp:47-53`; set from `USE_MKL_64BIT_API` at `CMakeLists.txt:135-136` |
| *(none for the backend)* | `BELFEM_ARMADILLO` / `BELFEM_BLAZE` must **not** be set for a plugin | `cl_Material.hpp:24` and siblings |

Both templates currently default a plugin to `Release` (`UserMaterialTemplate.cmake:141-143`) while
the host defaults to Debug (`CMakeLists.txt:90`, derived at `:115-123`) — a copy-paste plugin
against a default host is mismatched under a comment that says "match the host's build type".
R4/R1 should fix the default or the comment.

Also: `gTbulk` is `extern` in `globals.hpp:44` and used from `powerlaws.hpp` inlines, so a plugin
instantiating them takes an undefined symbol resolved from the host at load (`-rdynamic` on ELF,
`dynamic_lookup` on Mach-O). Pre-existing, now part of the published contract.

There is **no** version or ABI handshake at the `dlopen` boundary (DR-132, decided for 1.0). Out of
scope here; it remains the reason an ABI break is silent.

## 7. Definition-of-Done Checklist

- [ ] Every gap-table row mapped to a step or an open question.
- [ ] Each claimed gap backed by a citation or an explicit measurement.
- [ ] Ordered steps with dependencies; R1 and the template probes in one commit.
- [ ] R3's isolation demonstrated **once, during implementation**: delete `globals.hpp` from the
      stage and confirm the gate fails. If it still passes, the inherited `-I` was not removed and
      the gate is compiling against the source tree. A realistic failure — a transitive file
      dropped from the install list is the mistake that ships a broken prefix — and it needs no
      permanent fixture. *(Revision 2 first proposed injecting `#include "fn_polyval.hpp"` into
      `cl_Material.hpp`. Dropped 2026-08-29, Christian: the prohibition is written in that header,
      and a gate whose purpose is to catch a developer defying an explicit note is ceremony. R3 is
      a packaging gate — does the shipped list compile as shipped — not a discipline gate.)*
- [ ] O1, O2, O3, O5 logged; O4 resolved.
- [ ] **R8 passes** on Linux and Darwin, source-tree branch and install branch.
- [ ] Codex + Grok re-audit of this revision before code; a third audit on the diff.

## 8. Audit Trail

- `tmp/ai_exchange/review_install_headers.md` — 2026-08-29 jury round on the install question;
  origin of D1–D4.
- `tmp/ai_exchange/user_api_plan_audit.md` — 2026-08-29 plan audit, revision 1. **Verdict: changes
  required.** Grok raised three P0s (R1 before the template probes; R4 as a deletion list; R3 not
  isolating the inherited `-I`), all confirmed against the tree and all fixed above. Codex raised
  the missing `BELFEM_INT64` contract, the unexecutable R8, the staging-staleness risk, and the
  "shipped boundary ≠ access control" correction, all folded in. Two Codex claims from the earlier
  round (MPI and MKL "required" by installed headers) were refuted — both are `#ifdef`-guarded.
  One claim of my own was retracted: `vtktypes.hpp` VTK includes **are** guarded.
- Re-audit of revision 2: pending.

---

## Defect tracker

- [ ] **D1 — The umbrella is in `src/fem/kernel/`. HIGH.** Among 34 headers, none of the 14 API
      headers among them; neither template puts `fem/kernel` on the include path
      (`UserMaterialTemplate.cmake:78`, `UserLibraryTemplate.cmake:99-101`), so no plugin can find
      it. Guard is `BELFEM_BELFEM_USER_API_HPP` against the house `BELFEM_VERSION_HPP` style, and
      it carries no licence block and no prohibition note. *(Revision 1 justified this partly with
      "none of which is backend-free" — overstated; `en_FEM_SolverAlgorithm.hpp:16-25` is a bare
      enum. The placement is still wrong, for the include-path reason.)* Fixed by R2.
      Found by Claude 2026-08-29.
      **Half fixed 2026-09-04**: guard is `BELFEM_USER_API_HPP`, the LBNL block is present (added
      by Christian before 09-04), and the prohibition note is in. The placement is deferred with R1.
- [x] **D2 — the example decks' plugin `CMakeLists.txt` cannot compile against an install. HIGH.**
      Originally raised against `examples/2D_Undulator/lib/CMakeLists.txt`: `:58-61` set the include
      path to `${BELFEM_DIR}/include` — wrong root, no module dirs. Also `:75` C++14 vs project
      C++17 (`CMakeLists.txt:62`); `:91` `SHARED` vs the templates' `MODULE`
      (`UserLibraryTemplate.cmake:167`); `:105-106` forced `.dylib`/`.dll` though the only loader is
      POSIX `dlopen` (`cl_SourceFunction.cpp:224`); `:28` hardcoded a developer's home path.
      Found by Claude, extended by Codex + Grok, 2026-08-29.
      **The C++14 and hardcoded-path halves died with that directory** in the 2026-08-30 migration
      (`devlog/dl20260831_examples_usermat_migration.md`); the wrong-include-root half was carried
      into all four replacement decks unchanged, which is what made it worth fixing rather than
      closing.
      **Fixed 2026-08-31 by Claude** in `examples/{tape_quench_usermat,block3d,undulator2d,disk_pulse}/src/CMakeLists.txt`:
      both branches now derive a `BELFEM_HEADER_ROOT` (`<tree>/src` or `<prefix>/include/belfem`)
      and iterate the same `BELFEM_MODULE_DIRS`, and the `<belfem_user_api>` forwarder is generated
      whenever the flat umbrella is absent rather than only for a source tree. Gate: real `cmake` +
      `make` for all four decks against a staged prefix carrying headers only, and against the
      source tree — 8/8 produce their `.so`. These decks stay `SHARED`, deliberately: the deck's
      `file :` key names the library exactly as it appears on disk and `PREFIX ""` / `SUFFIX ".so"`
      already pin that name on both platforms.
- [x] **D3 — Neither template, nor any `example_user_*.cpp`, is installed. HIGH.** The glob matches
      `*.hpp` only. Fixed by R5. Found by Grok 2026-08-29.
      **Fixed 2026-08-31 by Claude**: both templates and all three `example_user_*.cpp` install to
      `${CMAKE_INSTALL_DATADIR}/${LIBPREFIX}/templates/`. Gate: the rule was configured and
      `make install DESTDIR=` run in isolation (5/5 files land), and both templates were configured
      and built against a staged headers-only prefix, producing `libmyalloy.so` and `libcustom.so`.
- [x] **D4 — `LICENSE` is not installed. MEDIUM.** Installed headers carry the LBNL/UC block
      referring to "the top-level LICENSE file". *(The umbrella and `belfem_version.hpp` carry no
      block at all — see D1.)* Fixed by R7. Found by Codex and Grok 2026-08-29.
      **Fixed 2026-09-04** by R7.
- [x] **D5 — `UserDefinedMaterial::evaluate_derivative_of_polynomial` declared but not defined.
      CRITICAL.** `nm` on the TU's own object showed the symbol `U`. Link error on the next full
      build. Found by Claude 2026-08-29 (compile + `nm`).
      **Fixed 2026-08-29 by Christian**: inline definition at `cl_Material_UserDefined.hpp:441-452`.
      Verified by Claude — symbol now `W`; TU compiles clean under `-DDEBUG` and `-DNDEBUG`.
- [ ] **D6 — `UserLibraryTemplate.cmake` requires a backend the API forbids. HIGH.** `:55`, `:65-66`,
      `:99-101`, `:175` all serve a coupling `cl_SourceFunction.hpp:43` prohibits. All four shipped
      examples compile against the flat set with no backend define. It survived because **no gate
      has ever compiled a source-function plugin** — `example_user_source.cpp` and
      `example_user_defect.cpp` appear in no CMake file. *(Revision 1 attributed that to D3; wrong
      — it explains D6.)* Fixed by R3 + R4. Found by Christian 2026-08-29.
- [x] ~~**D7 — umbrella has no trailing newline.**~~ **FALSE POSITIVE (retracted 2026-08-29).**
      The file's last byte is `0a`. The "\ No newline" marker in the diff belonged to
      `cl_Material_UserDefined.cpp`. Raised by Claude, refuted by Grok.
- [x] **D8 — Dead `Vector` forward declaration in the API. LOW.** `cl_Material.hpp:85` declares
      `template< typename T > class Vector;` and the file contains **zero** uses of `Vector<` —
      sitting in a header whose stated contract is that linalg never appears. Delete before someone
      "completes" it. Fixed by R2. Found by Grok 2026-08-29.
      **Already gone from the tree on 2026-09-04** (no `class Vector;` in `cl_Material.hpp`);
      removed by Christian, date not recorded.
- [x] **D9 — `cl_SourceFunction.hpp:12` guard is `CL_FUNCTION_HPP`. LOW.** Unprefixed and generic;
      unique today, but the file is about to sit in a flat public include directory. Fixed by R2.
      Found by Grok 2026-08-29.
      **Fixed 2026-09-04**: `BELFEM_CL_SOURCEFUNCTION_HPP`; double-include and standalone syntax
      checks pass under `-DDEBUG` and `-DNDEBUG` with no backend define.
- [x] **D10 — A shipped header documents the wrong compile line. MEDIUM.**
      `cl_Material_UserDefined.hpp:120` shows `g++ … -I/path/to/belfem/include`, wrong against
      today's `include/belfem/<module>` layout and wrong after R1. It is one of the 14 files every
      user receives. Fixed by R6. Found by Grok 2026-08-29.
      **Fixed 2026-09-04** for the module layout that ships today: three `-I` lines plus a pointer
      to the installed `UserMaterialTemplate.cmake`. Goes stale again if R1 ever lands.
