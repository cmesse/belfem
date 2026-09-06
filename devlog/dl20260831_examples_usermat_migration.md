# block3d and undulator2d migrated to the `usermat` plugin layout

**Date:** 2026-08-31
**Purpose:** Record the migration of the last two example decks that still used the
pre-0.9.0 `custom` material subsection and the old `MatData/` / `lib/` plugin
directories, onto the `tape_quench_usermat` layout.

---

## What was wrong

`examples/block3d` still declared its bulk conductor through a `custom { }`
subsection. That spelling was renamed to `usermat` in BELFEM 0.9.0 and is now a
**hard error** with a rename message (`cl_MaterialFactory.cpp`, the
`section_exists( "custom" )` arm) — the deck could not run at all. The audit
recorded on the same day states that "no deck uses it"; `examples/block3d` is an
untracked directory and was not in that sweep's scope.

`examples/undulator2d` did not use `custom` in its materials, but carried the
same generation of plugin scaffolding: a `lib/` directory holding a **stale
copy** of `UserLibraryTemplate.cmake` — C++14, a hardcoded `BELFEM_DIR` pointing
at a collaborator's home directory, the material headers included one by one,
and no `BELFEM_USER_DATA_DIR`. The shipped template has since moved on and the
copies never followed; that is the actual failure here, not a defect in the
template.

## What changed

Both decks now mirror `examples/tape_quench_usermat`:

| | before | after |
|---|---|---|
| block3d | `MatData/` → `libmat.so` | `src/` → `usermat.so`, `data/` |
| undulator2d | `lib/` → `libcustom.so` | `src/` → `usersource.so`, `data/` |

- **`src/CMakeLists.txt`** is the reference file, per-deck: C++17, `BELFEM_DIR`
  defaulting to `$ENV{HOME}/codes/belfem` and overridable on the command line,
  the `<belfem_user_api>` shim generation for a source tree, `PREFIX ""` so the
  library appears on disk exactly as the deck names it, and
  `BELFEM_USER_DATA_DIR` compiled in from `../data`.
- **The plugin sources** now include the `<belfem_user_api>` umbrella instead of
  naming `cl_Material.hpp` / `cl_SourceFunction.hpp` and the standard headers
  individually.
- **`block3d/input.conf`**: `custom { file : MatData/build/libmat.so ; }` →
  `usermat { file : src/build/usermat.so ; }`.
- **`undulator2d/input.conf`**: both `current` boundary conditions point at
  `src/build/usersource.so` instead of `lib/build/libcustom.so`.
- **`data/README.md`** in each deck records that the directory is wired through
  `BELFEM_USER_DATA_DIR` and how to add a table. Neither deck ships one today —
  block3d's bulk material is constants only, and the undulator's transport
  current is an analytic trapezoid.

Two documentation slips in the sources were corrected while they were open: the
undulator's current was documented as "interpolation from data file" (it never
read one), and its init-function comment described the material-plugin signature
`<Name>_init(Material*)` rather than the source-plugin `<Name>_init(SourceFunction*)`
that the file actually defines.

Nothing else moved. `undulator2d/ParamsGeo.dat` stays beside `undulator.geo`,
which `Include`s it relative to its own directory.

## Evidence

**Verified.** Both plugins were configured and built out-of-source against this
tree, and the exported init symbols were read back from the libraries:

```
usermat.so     -> bulk_init                          deck label: bulk
usersource.so  -> MyCurrent_init, MyCurrent_reverse_init
                                                     deck labels: MyCurrent, MyCurrent_reverse
```

Both compiled warning-free under `-Wall -Wextra -Wno-unused-parameter`.

**Not verified:** neither deck was executed. Both need a mesh built from their
`.geo` first, and `block3d` additionally addresses several hundred sidesets that
only the generated mesh carries. The plugin now loads and registers, which is
the failure the `custom` spelling caused; whether either deck converges is a
separate question.

## C++17: where the tree actually stands

Checked after the migration, because the two decks' C++14 raised the question.
The framework standard is **C++17 and is not in doubt**: `CMakeLists.txt:62` and
`config/compiler/config_gcc.cmake:2` both set `CMAKE_CXX_STANDARD 17`, and the
real compile line in a configured tree carries `-std=gnu++17`. Both shipped
plugin templates already agree — `UserMaterialTemplate.cmake:207` and
`UserLibraryTemplate.cmake:230`. Only the two examples' private copies were
behind, and they no longer are.

C++14 was **not** breaking those plugins: the `<belfem_user_api>` umbrella still
compiles clean at `-std=gnu++14`, tested. The framework uses almost nothing
C++17-only — one `if constexpr`, no `std::optional`, no nested-namespace
definitions. So raising the examples is trap-prevention, not a bug fix: the
first C++17 construct to land in a header a plugin consumes would otherwise
surface as a plugin-only compile error with no obvious cause.

One place in the tree still contradicts C++17. `config/compiler/config_icc.cmake:44`
appends `-std=gnu++14` to `BELFEM_CXXFLAGS`, which `finalize_compiler.cmake:7`
folds into `CMAKE_CXX_FLAGS`. In practice it loses — CMake emits the
`CXX_STANDARD` flag *after* `CMAKE_CXX_FLAGS` (visible in the position of
`-std=gnu++17` in any `flags.make` here) and the last `-std` on a GCC or Intel
command line wins — so the Intel build compiles at 17 despite the flag. It is a
dead flag that reads as policy, and it should be deleted rather than left to be
believed. Not touched here: it is a compiler-config change, outside this task.

## Left open

- The `-std=gnu++14` in `config/compiler/config_icc.cmake:44`, above.
- **Template and examples disagree on the artifact convention.** The templates
  build a `MODULE` named `lib<name>.so`, with a comment giving the reason (the
  artifact is only ever `dlopen`ed, and MODULE yields the same name on macOS).
  `tape_quench_usermat` — and therefore the two decks migrated onto it — build a
  `SHARED` library with `PREFIX ""`, so `usermat.so` with no prefix. Both work;
  the deck names whatever is on disk. But a user reading the template and a user
  reading an example are told different things, and the template's stated
  reasoning is the better of the two. Worth reconciling in one direction.
- The templates also pass `BELFEM_${BELFEM_BACKEND}` to the plugin, which the
  per-example CMakeLists do not. The template's own header states the plugin API
  is backend-free, so this is likely belt-and-braces rather than a real
  divergence — but it was not chased down here.
- The `custom` subsection is still absent from the schema's refused-input
  vocabulary. With `block3d` migrated, no deck in the tree uses it, so this is
  now purely a completeness gap in `doc/input_schema.yaml`.
