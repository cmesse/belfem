# Devlog 2026-05-27 - macOS MUMPS Install Names

**Date:** 2026-05-27
**Topic:** Investigation of `hphirun` failing to load MUMPS dylibs on macOS
**AIs involved:** Codex
**Claude Confidence:** N/A
**Codex Audit Confidence:** high
**Literature References:** N/A

## Summary

Investigated the macOS dyld abort where `hphirun` fails to load `../lib/libdmumps.dylib` even though the MUMPS library is under `$SCLS/lib` (`/opt/scls/lib`). The executable already has `/opt/scls/lib` in its `LC_RPATH`; the original failure came from malformed install names embedded in the SCLS MUMPS dylibs. After those external dylibs were fixed, the already-built `hphirun` still had stale `LC_LOAD_DYLIB` entries and had to be relinked or patched.

## Key Findings

- `cmake-build-debug/bin/hphirun` directly records bad MUMPS load commands: `../lib/libdmumps.dylib`, `libmumps_common.dylib`, and `../../lib/libpord.dylib`.
- The executable also records `LC_RPATH /opt/scls/lib`, so missing rpath is not the primary issue.
- `/opt/scls/lib/libdmumps.dylib` has install ID `../lib/libdmumps.dylib`.
- `/opt/scls/lib/libmumps_common.dylib` has install ID `libmumps_common.dylib`.
- `/opt/scls/lib/libpord.dylib` has install ID `../../lib/libpord.dylib`.
- `config/linalg/config_mumps.cmake` links MUMPS with `-ldmumps`, `-lmumps_common`, and `-lpord`; on macOS the linker copies each dylib's install ID into the dependent binary.
- `CMakeLists.txt` sets `CMAKE_OSX_DEPLOYMENT_TARGET` to the full `sw_vers -productVersion` when the environment variable is unset, explaining the `built for macOS 15.7.7` load-command warning in this build.
- After the external SCLS fix, `otool -D` reports `/opt/scls/lib/libdmumps.dylib`, `/opt/scls/lib/libmumps_common.dylib`, and `/opt/scls/lib/libpord.dylib` as the current install IDs.
- A normal `cmake --build cmake-build-debug --target hphirun -j2` did not relink the executable, so stale load commands remained in the existing binary.
- Patching `cmake-build-debug/bin/hphirun` with `install_name_tool -change ...` and re-signing it replaced the stale relative MUMPS load commands with absolute `/opt/scls/lib/...` entries.
- The patched `hphirun` no longer aborts in dyld. It reaches BELFEM startup and fails later because `/Users/christian/codes/belfem/input.conf` does not exist.

## Changes Made / Proposed

- No source or build-system files were modified.
- Patched generated build artifact `cmake-build-debug/bin/hphirun` and re-signed it ad hoc.
- Proposed permanent package-level fix: keep the corrected install IDs and internal dependency names in the SCLS MUMPS dylibs, then force relink BELFEM executables after the package fix.
- Proposed repo-level workaround, if editing is approved: add a macOS/MUMPS post-link `install_name_tool -change` step for generated executables before codesign in `config/scripts/Add_Executable.cmake` and `config/scripts/Add_Test.cmake`.
- Proposed deployment-target cleanup, if editing is approved: avoid defaulting `CMAKE_OSX_DEPLOYMENT_TARGET` to the full current OS patch version; use an explicit supported minimum or let the user-provided environment/cache value control it.

## Open Questions

- Should BELFEM carry a post-link workaround for stale or malformed MUMPS install names, or is a forced relink after SCLS package repair sufficient?
- What minimum macOS deployment target should BELFEM use for local SCLS builds?

## Files Updated

- cmake-build-debug/bin/hphirun
- devlog/dl20260527_mumps_macos_install_names.md
- devlog/README.md
