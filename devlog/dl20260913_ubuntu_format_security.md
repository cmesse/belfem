# Ubuntu's GCC spec file injects -Wformat-security; two core fixes, one of them a latent bug

**Date:** 2026-09-13
**Purpose:** Record why `-Werror=format-security` appeared on the `ubuntu_patch` branch and not on
Rocky, what was changed, and why the remaining nine `-Wformat` pragma sites need nothing.
**Module:** src/core (`fn_sprint.hpp`, `cl_Logger.hpp`); config/compiler (read only)

## Symptom

Release build on Ubuntu 24.04 (gcc 13.3.0-6ubuntu2~24.04.1), `-Wall -Werror` from
`config/compiler/config_gcc.cmake`:

```
src/core/fn_sprint.hpp:47: error: format not a string literal and no format arguments [-Werror=format-security]
src/core/fn_sprint.hpp:53: error: ...
src/core/cl_Logger.hpp:126: error: ...
```

The same tree builds clean on Rocky.

## Root cause: a distro patch, not a BELFEM regression

Ubuntu patches the GCC driver spec file. Verified on the build host:

```
$ gcc -dumpspecs | grep format-security
%{!Wformat:%{!Wformat=2:%{!Wformat=0:%{!Wall:-Wformat} %{!Wno-format-security:-Wformat-security}}}}
```

Unless the command line carries an explicit `-Wformat`, `-Wformat=2` or `-Wformat=0`, the driver
appends `-Wformat-security` to every compile. A bare `printf(m)` warns with no flags at all.
Upstream GCC does not include `-Wformat-security` in `-Wall` (it needs `-Wformat-security` or
`-Wformat=2` explicitly), so on a stock toolchain `-Wall -Werror` never reaches the diagnostic.
Rocky is RHEL-derived; Red Hat's hardening lives in `redhat-hardened-cc1` and is applied only
through the RPM build macros, so a plain CMake-driven `gcc` gets upstream defaults. Ubuntu side
verified here; Rocky side reviewed from packaging knowledge, not run (one-liner to confirm there:
`gcc -dumpspecs | grep -c format-security` prints 0).

It is *our* `-Werror` that promotes the injected warning to an error, not Ubuntu's.

Not taken: adding a bare `-Wformat` to `BELFEM_CXXFLAGS` silences the whole injection (confirmed
with `gcc -Wformat -c`), and would restore the Rocky blind spot in a way no future reader would
connect to the spec file. Also not taken: `-Wno-format-security` project-wide, for the same reason.

## Fixes (two lines, both approved by Christian)

**`fn_sprint.hpp:26`** — the existing GCC branch ignored only `-Wformat`; on GCC that does not
cover the `-security` sub-option (tested: `ignored "-Wformat"` still errors, `ignored
"-Wformat-security"` compiles). Added the second pragma. `sprint()` is a variadic template that
forwards a runtime format by design; `__attribute__((format(printf, ...)))` cannot be applied to a
template, so suppression is the correct tool. Clang branch already had the right name.

**`cl_Logger.hpp:126`** — `std::fprintf( mStream, tMessage.c_str() )` became
`std::fprintf( mStream, "%s", tMessage.c_str() )`. This one is not a pragma case: the warning was
correct. `tMessage` is the already-formatted output of `sprint()`, and passing it as the format
string parsed it a second time. Any `%` surviving the first expansion (a path, a material or mesh
name, `"50%% done"` which `sprint` correctly collapses to `50%`) was reinterpreted as a
conversion reading a vararg that was never pushed: undefined behaviour on every platform,
including Rocky, which simply never reported it. Ubuntu's spec patch surfaced a latent defect
rather than introducing one.

## Why the other nine `-Wformat` pragma sites need nothing

Every remaining GCC `diagnostic ignored "-Wformat"` block was read:

| Site | Region contents | `-Wformat-security` reachable |
|---|---|---|
| `cl_Cell.hpp`, `cl_AR_Vector/Matrix.hpp`, `cl_BZ_Vector/Matrix.hpp` | `print_matlab()`, literal formats | no |
| `banner.cpp` | ~40 `fprintf(stdout, "literal")` | no |
| `stringtools.hpp`, `cl_SourceFunction.hpp`, `cl_Material.hpp` | no printf-family call inside the block at all | no (pragma is dead weight; left alone) |

The diagnostic fires only on a non-literal format with no further arguments. A repo-wide search
of `src/`, `tests/` and `examples/` for that shape (printf/sprintf first argument, fprintf second,
snprintf third not starting with `"`) finds only the two fixed sites. `sprint()` callers are not
at risk: for a template GCC reports at the definition site, which is inside the suppressed region
(confirmed by instantiating from outside the block in the repro). `nonfree/` is not checked out on
this host and was not searched.

## Not an Ubuntu-vs-Rocky difference (correction made in-session)

PIE was first listed as a further `ubuntu_patch` hazard; withdrawn. RHEL 8+ also ships
`--enable-default-pie`, and `CMakeLists.txt:343` sets `CMAKE_POSITION_INDEPENDENT_CODE ON`
globally, so no BELFEM object can trigger a "recompile with -fPIE" link error on either platform.
The same spec line also injects `-fstack-protector-strong`, `-fstack-clash-protection` and
`-fcf-protection`; those are codegen with no build-breaking consequence.

## Status

Reviewed, not verified: the executable gate is Christian's build, which passed `fn_sprint.hpp` and
`cl_Logger.hpp` after the edits. Working tree: `src/core/fn_sprint.hpp` (+1),
`src/core/cl_Logger.hpp` (+1/-1). Not committed by the AI.
