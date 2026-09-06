# `./cmake` folded into `./config/doc`

**Date:** 2026-09-04
**Purpose:** Pre-release source tidy-up — remove the top-level `cmake/` directory by moving its
five Doxygen helper scripts under `config/`, and rewire every call site.
**Module:** meta (build system, documentation toolchain)

## What moved

The top-level `cmake/` directory held nothing but the `make doc` toolchain: five `-P` scripts, no
modules, no `find_*` package files. It sat beside `config/`, which is where every other piece of
CMake machinery already lives (`config/compiler`, `config/system`, `config/linalg`, `config/io`,
`config/numerics`, `config/scripts`). Christian's call for the release tidy-up: fold it in.

```
cmake/normalize_doxyfile.cmake         →  config/doc/normalize_doxyfile.cmake
cmake/normalize_doxygen_layout.cmake   →  config/doc/normalize_doxygen_layout.cmake
cmake/strip_doxygen_html_header.cmake  →  config/doc/strip_doxygen_html_header.cmake
cmake/patch_doxygen_header.cmake       →  config/doc/patch_doxygen_header.cmake
cmake/report_doxygen_warnings.cmake    →  config/doc/report_doxygen_warnings.cmake
```

Done with `git mv`, so the history of each script follows it. `cmake/` was then `rmdir`ed.

## What was rewired

| Site | Change |
|---|---|
| `CMakeLists.txt:575,586,598,608,612` | the five `-P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/…` invocations inside the `doc` target |
| `CMakeLists.txt:539` | the comment pointing a reader at the header patcher |
| `Doxyfile.in:1455` | same, in the `HTML_HEADER` note |
| `config/doc/patch_doxygen_header.cmake:52` | the script's own path inside its `FATAL_ERROR` recovery text |
| `config/doc/report_doxygen_warnings.cmake:72` | the script's own path inside its baseline-update hint |
| `doc/lessons_learned_evidence.md:283,287` | INC-191 and INC-195 pointers |

## What was checked and deliberately left alone

- **`config/system/find_vtk.cmake:41`** — the `lib/cmake/vtk-${VTK_VERSION}` there is VTK's own
  install layout, not this directory. Untouched.
- **No `CMAKE_MODULE_PATH` anywhere in the tree.** Nothing resolved `cmake/` implicitly; all five
  call sites were explicit `-P` paths, so the rewire is complete by construction rather than by
  hoping a search caught everything.
- **Doxygen `INPUT`** is `doc/`, `src/`, `examples/` — `config/` is not in it, so the moved
  scripts do not enter the generated site and the page tree is unchanged.
- **`scripts/check_doc_claims.py`** scans `config/**` for `fno-exceptions` and every tracked
  `*.cmake` for `add_library` kinds. The five scripts contain neither, so no fact count shifts.
- **Dated `devlog/` entries and `todo/closed/`** still read `cmake/…` in ~25 places. Records are
  kept as written; the same rule that preserved the retired `paperN` aliases applies.

## Verification

Each of the five scripts was executed standalone under `cmake -P` from the new location. All five
parse and run, each failing with its *own* scripted guard message (`normalize_doxyfile: no
configuration at`, `patch_doxygen_header: no header at`, and so on) rather than a parse error —
which is what the move could plausibly have broken. `scripts/check_doc_claims.py` reports 38/38.

**`make doc` was not run.** That is the real gate for the `doc` target and it belongs to Christian;
the standalone parse check confirms the files survived the move, not that the target still
assembles a site. Reviewed, not verified.
