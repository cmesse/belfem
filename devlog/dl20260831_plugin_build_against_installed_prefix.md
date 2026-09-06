# A user plugin can now be built against an installed prefix

**Date:** 2026-08-31
**Purpose:** Record the repair of the two mechanisms that made an installed BELFEM
unable to build a user material or source plugin, and the evidence for both.
**Module:** `CMakeLists.txt`, `examples/*/src/CMakeLists.txt`,
`src/physics/materials/User{Material,Library}Template.cmake`

---

## The question that started it

*"When I have installed belfem but not the sources, we only have the include files in
`$DESTDIR/include/belfem`. Do the user defined material scripts work?"*

Short answer: the headers were sufficient, the build glue was not.

## What was measured first

The install rule was replayed into a staging directory — `install( DIRECTORY src/ …
FILES_MATCHING "*.hpp" )` — giving 628 headers under `include/belfem/<module>/` and no
sources. Against that prefix:

- **The API side already worked.** All three `example_user_*.cpp` compiled clean under
  `-DDEBUG` and `-DNDEBUG`, no backend define; a real `-c` plus `-shared` produced a
  `.so` carrying 62 undefined `belfem::` symbols for the host to resolve at `dlopen`
  (`-rdynamic`, `config/compiler/config_gcc.cmake:104`). Nothing about a plugin needs
  the sources — it includes headers and links nothing.
- **The glue side failed twice.** The four example decks' plugin `CMakeLists.txt` — the
  only plugin build files `make install` ships — had a dead install branch:
  `set(BELFEM_INCLUDE_DIRS ${BELFEM_DIR}/include)`, the wrong root and no module dirs,
  with the `<belfem_user_api>` forwarder generated only when `${BELFEM_DIR}/src/...`
  existed. Result against a real prefix:
  `matlib.cpp:18:10: fatal error: belfem_user_api: No such file or directory`.
  And the two templates that *do* handle an installed prefix correctly were never
  installed — the header glob matches `*.hpp`, so no `.cmake` and no `.cpp` reached the
  prefix. The working recipe was the one artifact the install did not deploy.

A single `-I <prefix>/include/belfem` is not enough either, and that was confirmed
separately: the headers include each other by bare name, so the module `-I` list is
mandatory under today's layout.

## What changed

**The four decks** (`tape_quench_usermat`, `block3d`, `undulator2d`, `disk_pulse`) now
derive a `BELFEM_HEADER_ROOT` — `<tree>/src` or `<prefix>/include/belfem`, probed by a
header rather than a directory name — and iterate one `BELFEM_MODULE_DIRS` list for both
branches. The forwarder block triggers on the absence of the flat umbrella instead of the
presence of a source tree, so it serves the installed case and no-ops the day a prefix
ships the real header.

**`CMakeLists.txt`** installs both templates and all three `example_user_*.cpp` into
`${CMAKE_INSTALL_DATADIR}/${LIBPREFIX}/templates/`.

**Both templates** had their status message corrected: they printed "shim for the source
tree" while configuring an installed prefix.

## How it was verified

Executable gates, all rc=0:

| Gate | Result |
|---|---|
| 4 decks × {staged prefix, source tree}, real `cmake` + `make` | 8/8 produce their `.so` |
| `UserMaterialTemplate.cmake` against the staged prefix | `libmyalloy.so` |
| `UserLibraryTemplate.cmake` against the staged prefix | `libcustom.so` |
| New `install( FILES )` rule, isolated project, `make install DESTDIR=` | 5/5 files land |
| `scripts/check_doc_claims.py` after the `CLAUDE.md` edit | 37/37 claims hold |

The staged prefix carried headers only — no `src/`, no library — and `flags.make` was read
back to confirm the installed branch was the one taken.

**Not gated:** nothing was `dlopen`ed by an installed solver binary. These plugins build
and carry the right undefined symbols; whether the host resolves them at load is the
question O5 in `todo/user_api_header_install.md` raises and this session does not answer.
Darwin is untouched and unmeasured.

## Scope note

This repairs the *existing* module-layout install. It neither depends on nor prejudges
`todo/user_api_header_install.md`, which replaces that layout with a flat 16-file API set;
R1 there will have to rewrite the include-path block in all four decks and both templates.
R5/D2/D3 in that plan are closed by this session, out of their planned order — the
dependency on R4 was ordering convenience, not a requirement. R7 (`LICENSE` is still not
installed) remains open.
