# Shared `libbelfem`, `make install` and RPATH — implemented and gated on Darwin

**Date:** 2026-08-27
**Purpose:** Session record for R1–R8b and R10 of `todo/shared_library_and_install_plan.md`, done
on the Mac first because it was the machine at hand; Linux (R4, R9, R11, R12) follows.
**Module:** `CMakeLists.txt`, `config/`, `src/comm` (one function), `tests/` harness
**Platform:** macOS 15.7.9 x86_64, GCC 16.1.0 (SCLS), CMake 4.3.2, Blaze backend, Debug

## Commits (branch `claude`)

| sha | what |
|---|---|
| `4c62bfe1` | 15 printf format/argument mismatches in `BELFEM_ASSERT`/`BELFEM_ERROR` messages (pre-existing; found while asking what breaks on Apple Silicon) |
| `e7922179` | `SplineLookupTable::alpha_custom(Bezier*, Vector, T)` renamed `alpha_composite` — hid the virtual, failed `-Werror=overloaded-virtual` |
| `3cf66776` | R1 + R8b: `CMAKE_POSITION_INDEPENDENT_CODE`, x86-only `-m64`/`-mtune=native`, `BELFEM_DATA` as a ctest property |
| `5fbcb362` | R2 + R3: OBJECT libraries, one aggregate `belfem` target, `USE_SHARED_LIBS` (default OFF) |
| `a75f4b18` | R5–R8: RPATH, `install()` rules, installed data lookup, Darwin re-sign; O7/O8 recorded |

## What the gates showed

Every claim below ran; nothing here is static review.

- **R1:** full rebuild clean, zero warnings; the Fortran objects in `sparse/` carry PIC relocations
  (`otool -rv` shows GOT/BRANCH entries). `make check` 12/12 with `BELFEM_DATA` unset.
- **R3 static:** `lib/libbelfem.a` (674 MB, Debug) from 23 module object sets; `hphirun` links
  exactly one project library; 12/12.
- **R4-on-Darwin (shared):** first link **failed** — `Undefined symbols: _gComm, _gLog`. Both are
  defined in every `main()` (~50 sites incl. `nonfree/`), the library only has `extern`. Not a
  cycle; the plan's cycle analysis was right and this is a different thing. Interim:
  `-Wl,-undefined,dynamic_lookup` on the dylib, Darwin only. Recorded as **O7** with the real fix
  (define in the library) flagged for sign-off. After that: `libbelfem.0.9.0.dylib` +
  `.0` + unversioned symlinks in `lib/`, install name `@rpath/libbelfem.0.dylib`, `make check`
  **14/14** with gastables/gasmodels on — and the gastables suite **ran 10 tests** instead of
  skipping, which is R8b doing its job.
- **R10 (a), relocation:** `make install DESTDIR=<stage>`, tree moved. `otool -l hphirun` has
  `@loader_path/../lib`; `codesign -vv` valid on both executables and the dylib after the
  install-time rpath rewrite (the `install(CODE)` re-sign works); `banner` runs with
  `DYLD_LIBRARY_PATH` unset. The `2D_Tapestack` run then **failed**: `File sp-ap.hdf5 does not
  exist` — R8's compiled-in `/usr/local/share/belfem` does not exist at the moved location. With
  `BELFEM_DATA` pointed at the moved `share/belfem`, the run completes (3 timesteps,
  `iv_results.csv`). Recorded as **O8**; the R9 gate is now two runs.
- **R10 (b), in place:** reconfigured with a scratch `CMAKE_INSTALL_PREFIX`, `make install`, same
  run with `BELFEM_DATA` **unset**: completes. That is R8 proven.
- **R7 churn:** `touch share/material/zz_churn_probe.tmp; make install` (no cmake) — the file
  appeared in the prefix.
- Installed sizes: `examples/` 20 MB (run output excluded, every `.msh` kept), headers 5.9 MB,
  `lib/` 13 MB.
- `scripts/check_doc_claims.py`: 34/34 after the `CMakeLists.txt` and `CLAUDE.md` edits.

## Things worth knowing that are not in the plan

- `USE_GASMODELS` adds a define to every TU, so toggling it is a full rebuild.
- `ld: warning: duplicate -rpath '/opt/scls/lib'` on every Apple link: one entry is CMake's from
  `BELFEM_RPATH`, the other is injected by the SCLS GCC driver (it also adds its own
  `lib/gcc/<triple>/16.1.0`). Harmless; not ours to fix.
- `bin/` held 7 executables, not 4. Decided same day: `material` and `gas` (the latter under
  `USE_GASMODELS`) ship too; `pipette` and `corctest` were driver stubs and are deleted.
  `BELFEM_INSTALL_EXECUTABLES` (`globals.cmake`) is the list; the install rule lives in
  `Add_Executable.cmake` because `install(TARGETS)` across directories needs CMake 3.13 and the
  tree admits 3.11.
- Codesigning on Apple Silicon: the linker ad-hoc signs everything it emits, so third-party
  libraries and the build tree need nothing. Only post-link rewrites (CMake's install-time rpath
  edit) void a signature, hence the `install(CODE)`. The `POST_BUILD` sign is redundant on arm64
  and kept for the debugger's sake. Gated on x86_64, where an invalid signature would not have
  killed the process — the arm64 failure mode is prevented, not reproduced.
- Things that work here and break on Apple Silicon, from a read-only survey (no arm64 machine):
  `-m64` (fixed, now x86-only), `-mtune=native` on Fortran (same), the varargs mismatches above
  (fixed — Apple arm64 passes variadic args on the stack, so a too-short argument list reads the
  caller's frame), FMA contraction differences in `make check` tolerances (not addressed), and
  the `BELFEM_USE_CLANG` fallback that hardcodes an x86 libgcc path (dead path, left alone).

## Left for the Linux session

R4 (`.so.0.9.0` naming, `ldd` clean, `make reset` removes it), R9 (a) and (b), then R11 (flip
`USE_SHARED_LIBS` default, plugin templates) and R12. Open decisions needing Christian: **O2**
(per-module gate renamed — implemented, sign-off owed), **O7** (move `gComm`/`gLog` into the
library, ~50 files), **O4** (`belfemConfig.cmake`, one `install(EXPORT)` away). The Codex + Grok
audit of §8 has not run.

The `build/` tree was left configured with `USE_SHARED_LIBS=ON`, `USE_GASMODELS=ON`,
`CMAKE_INSTALL_PREFIX=/usr/local` (restored after the scratch-prefix gate); nothing was written
to `/usr/local`.
