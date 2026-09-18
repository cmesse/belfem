# Prefer the SCLS PMIx at Link Time, and Fail at cmake When MPI Cannot Link

**Date:** 2026-09-18
**Purpose:** A BELFEM 0.9.1 build on a colleague's machine compiled the whole library and then
failed at every executable link with 38 `undefined reference to PMIx_*` errors out of
`libmpi.so`. Two things are wrong, one of them in BELFEM. First, every BELFEM link line carries
`-Wl,-rpath,/usr/lib64:/opt/scls/<flavor>/lib64:...`: the system compiler's library directory
sits ahead of the SCLS prefix, so ld resolves Open MPI's `libpmix.so.2` dependency from
`/usr/lib64` whenever a distro PMIx is installed there, and that PMIx predates the 4.2 API Open
MPI 5 needs. Reproduced on the reference machine with a stub library: identical 38-symbol set.
Second, nothing at configure time performs an MPI link, because CMake's compiler sanity test
runs with plain `g++` before BELFEM swaps in `mpicxx`, so a broken MPI setup is discovered by
the first executable instead of by cmake. The mechanism: keep implicit system link directories
out of `BELFEM_RPATH`, then link a minimal MPI program in `find_mpi.cmake` with the wrapper the
build will use and stop the configure with the linker's own output when that fails.
**Module:** `config/compiler`, `config/system` (+ `doc/mpi_support.md`)
**AIs involved:** Claude (diagnosis, reproducer, plan), Codex (audit), Grok (third voice)
**Status:** ✅ COMPLETE (2026-09-18), build gate owed. Landed: R0 at both sites (`detect_gcc.cmake`
skip for a `/usr` or `/` prefix; `belfem_prune_system_libdirs()` called from
`finalize_compiler.cmake` and from the probe), the link probe with PMIx hint and `ldd` warning
in `find_mpi.cmake`, the summary line, `doc/mpi_support.md` (swept). Verified by scratch-tree
configures on the mkl flavor: rpath begins with the SCLS prefix, PATH-shim and
`LD_LIBRARY_PATH` negatives produce the fatal message and the warning, `check_doc_claims.py`
38/38. Three jury rounds (`devlog/dl20260918_mpi_link_probe_rpath.md`). Residual: `make` on
both trees (Christian), Darwin and Intel configures unexercised, the post-`project()` compiler
swap recorded as pre-existing.

> **Scope guards:**
> - R0 changes the `-rpath` list of every executable (build and install) by *removing* exactly
>   four directories, when present: `/usr/lib64`, `/usr/lib`, `/lib64`, `/lib` (compared by
>   REALPATH). Those are what the loader searches by default, so no library is found from a
>   different place afterwards, except a system copy that was shadowing a prefix copy, which is
>   the defect. A compiler prefix that is not `/usr` or `/` keeps its libdir in the rpath; that
>   entry is why `detect_gcc.cmake:77-82` and `detect_icc.cmake:83` exist.
> - The probe (R1–R4) adds no build flag and no link input. It is a throwaway `execute_process`.
> - **Behaviour change, stated:** today a failed `mpicxx -E` warns and continues
>   (`find_mpi.cmake:60-67`). A failed link probe stops the configure. A wrapper that cannot
>   link a two-call program cannot link `belfem`; postponing that is the 0.9.1 failure with
>   extra steps. (O3, Christian decides.)
> - No attempt to make an older PMIx work. No `-lpmix` from BELFEM: it would hide an absent
>   package and is wrong for a flavor that uses an external PMIx on purpose.
> - `--disable-new-dtags` (pinning the prefix over `LD_LIBRARY_PATH` at run time) stays OUT
>   of scope: a global loader-policy change; `CMakeLists.txt:267-276` already forbids raw
>   `-Wl,-rpath` on targets for a related reason.
> - `ALLOW_UNTESTED_MPI` semantics untouched.

---

## 1. Current Behaviour and How It Fails

| Failure | Mechanism | Evidence |
|---|---|---|
| System PMIx wins over the SCLS PMIx at link time | `detect_gcc.cmake:70-82` appends the compiler prefix's libdir to `BELFEM_RPATH` first; with the system GCC that is `/usr/lib64`, also listed in `CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES`. `CMakeLists.txt:276` turns the list into `-Wl,-rpath,…` on every target (`build/src/executables/CMakeFiles/belfem.dir/link.txt`: `/usr/lib64` first). GNU ld resolves a shared library's NEEDED entries through `-rpath` directories in order, before `LD_LIBRARY_PATH`, before the library's own RUNPATH, before the default directories (`man ld`, `-rpath-link`, items 2–7) | reproducer: stub `libpmix.so.2` with the 24 classic symbols, `mpicxx -Wl,-rpath,<stub>:/opt/scls/mkl/lib64 hello.cpp` → 38 undefined references, set identical to the colleague's log; reversed order → links, `ldd` shows the SCLS copy |
| Configure passes on a tree that cannot link | `project()` runs CMake's compiler test with `$SCLS/bin/g++` (`CMakeLists.txt:15-21`, `:28`); `mpicxx` replaces it afterwards (`detect_gcc.cmake:147-150`, via `detect_compiler.cmake:28-33` from `CMakeLists.txt:171`). No test program is ever linked against `libmpi.so`; the existing probe is `-E` only (`find_mpi.cmake:53-58`) | colleague's log: `libbelfem.a` built, first executable fails; GNU ld's default `--no-allow-shlib-undefined` applies to executables, not archives |
| The failure is late and unreadable | every executable prints the same 38 lines; nothing names PMIx as a package or says which `libpmix.so.2` was used | colleague's log |

**Diagnosis of the colleague's box.** Open MPI 5.0.10 imports 62 PMIx symbols (`nm -D
--undefined-only` over `libmpi.so` + `libopen-pal.so.80` on the reference machine). The
colleague's linker satisfied 24 (`PMIx_Init`, `PMIx_Get`, `PMIx_Fence`, …) and failed on the 38
that became functions in PMIx 4.2 plus `pmix_framework_names`. So a `libpmix.so.2` older than
4.2 was found, not none (an absent library prints `needed by … not found`). SCLS ships PMIx
5.0.10 as `scls-<flavor>-pmix`, required by `scls-<flavor>-openmpi`, and links Open MPI with
`--with-pmix=%{prefix}` (`libmpi.so` NEEDED `libpmix.so.2`, RUNPATH `<prefix>/lib`; verified on
the reference machine, SCLS recipes are not in this tree). Two candidate causes, both fully
consistent with the log: (a) `scls-gcc-pmix` not installed; (b) a distro PMIx in `/usr/lib64`
(EL9 ships 3.2.x) shadowing an installed SCLS PMIx through the `-rpath` order above. The handout
`tmp/pmix_diagnostics_for_colleague.md` distinguishes them. Confidence high on the mechanism,
open on which of (a)/(b) applies there.

**Bottom line:** BELFEM promotes the system library directory above the SCLS prefix on every
link, which is the opposite of the intended preference; and the configure step never performs
the one link that would have failed. R0 fixes the preference; R1–R4 make the remaining failure
modes (absent package, broken wrapper) legible at cmake time.

## 2. Architecture

**R0 — remove the system libdir, keep every prefix libdir.** The defect is one append:
`detect_gcc.cmake:77-82` records the compiler prefix's libdir so that a non-system compiler's
runtime (`libstdc++`, `libgfortran`, oneAPI's `libsvml`) is found at run time; with the system
GCC the prefix is `/usr` and the append promotes `/usr/lib64` above the SCLS prefix at link
time. Rule, at the append site: skip when the prefix REALPATHs to `/usr` or `/`. Rule, in
`finalize_compiler.cmake` after the dedupe: drop entries whose REALPATH is one of `/usr/lib64`,
`/usr/lib`, `/lib64`, `/lib`, whoever appended them. **Rejected:** removing everything in
`CMAKE_<LANG>_IMPLICIT_LINK_DIRECTORIES` (round 2, Codex and Grok independently): that list
is the compiler driver's `-L` set, which on a prefix compiler includes its own libdir, and on
this host already includes `/usr/lib/gcc/x86_64-redhat-linux/11`; a filter on it would strip
the very directory the append exists for. Rejected: reordering the list (fragile) and dropping
the compiler libdir unconditionally (breaks the prefix-compiler case).

**R1–R4 — extend the existing probe pattern.** A second source file, `MPI_Init` /
`MPI_Finalize`, linked with `execute_process( COMMAND ${CMAKE_CXX_COMPILER} … )` exactly as the
flavor probe invokes the wrapper (`find_mpi.cmake:54-58`). That variable is the bare name
`mpicxx` (`detect_gcc.cmake:149`), the driver every `link.txt` invokes, while CMake's cache
and compiler identification still describe the `g++` seen at `project()`; the swap after
`project()` is a pre-existing design outside this plan, and the probe tests the driver the
links use, which is the claim that matters. The call to `MPI_Init` is
load-bearing: under a toolchain that defaults to `--as-needed`, an empty `main` drops `-lmpi`
and the probe proves nothing. The binary is never executed. Rejected: `try_compile()`, which
spawns a sub-configure and reads the compiler from the cache (`g++`, `CMakeLists.txt:15-18`),
the very hole this plan closes; `execute_process` is the file's own precedent and its output is
the linker's own words, which the message must quote. Probe stays in `find_mpi.cmake` for
fail-fast: no MUMPS / PETSc / HDF5 configure when `mpicxx` cannot link `MPI_Init` (O2).

## 3. Gap Table

| # | State | Needed for | Handled today? | Class | Citation / rationale |
|---|---|---|---|---|---|
| 0 | `BELFEM_RPATH` free of the canonical system libdirs | the SCLS copy of any library winning over a distro copy at link time | no | (c) | source: `detect_gcc.cmake:77-82` only (`detect_icc.cmake:83` appends the oneAPI libdir, never a system dir, and is left alone); consumers `finalize_compiler.cmake:86-88,99-102`, `CMakeLists.txt:276,284`; on this host `/lib64 -> usr/lib64`, so REALPATH both sides |
| 1 | The wrapper can link an MPI executable | catching a broken MPI at configure time | no | (c) | `find_mpi.cmake:41-58` preprocesses only |
| 2 | The probe carries the `-rpath` entries the real link will carry | avoiding a **false failure** when the library lives only in a directory the wrapper's own `-rpath` does not name (lib / lib64 split) | n/a | (c) | pass the whole `BELFEM_RPATH` as collected so far, after the R0 filter: compiler libdir if non-implicit, `$SCLS/lib64` or `lib` (`find_scls.cmake:4-9`), `$MPI_HOME/lib64` or `lib` when set (`find_mpi.cmake:15-18`). TPL directories are appended later (`CMakeLists.txt:191-249`); residual in O2 |
| 3 | The message names the package when the failure is PMIx-shaped | one actionable sentence instead of 38 lines | no | (c) | match `PMIx_` **or** `libpmix.so.2[^\n]*not found` (O4 folded in); name `scls-${BELFEM_SCLS_FLAVOR}-pmix` when the flavor is known (`find_scls_flavor.cmake:27-51`); wording must not assert a version, only "missing or incompatible with this Open MPI" |
| 4 | Which `libpmix.so.2` the loader picks | warning about `LD_LIBRARY_PATH` shadowing (loader order: `LD_LIBRARY_PATH` before `DT_RUNPATH`) | no | (c), Linux only | `ldd` on the probe binary, parse `libpmix.so.2 => <path>`; handle `not found` without REALPATH; compare REALPATHs with a trailing slash on `$SCLS` so `/opt/scls/mkl` does not match `/opt/scls/mkl-old` |
| 5 | Flavors without an SCLS PMIx (external PMIx on purpose, or package absent but R2 passed) | not warning where the warning would be false | n/a | (c) | O1 RESOLVED: warn only when `EXISTS "${SCLSLIBDIR}/libpmix.so.2"` (`find_scls.cmake:6,10`, set before `find_mpi` runs) and `ldd` resolves elsewhere; the `lbl` flavor is the example in the comment, not a name the code checks |
| 6 | Open MPI 4.x with internal PMIx, MPICH under `ALLOW_UNTESTED_MPI` | no `libpmix` line in `ldd` | n/a | (a) | skip row 4 silently |
| 7 | macOS | no `ldd`; `otool -L` shows direct dependencies only | n/a | (a) | row 4 gated on `CMAKE_SYSTEM_NAME STREQUAL "Linux"`; rows 0–3 run everywhere |
| 8 | A `libmpi.so` path for the message | the `ldd … | grep pmix` hint | no | (c) | `BELFEM_MPIHOME` exists only under `MPI_HOME` (`find_mpi.cmake:2-5`); in the Open MPI branch use the first token of `${CMAKE_CXX_COMPILER} --showme:libdirs` (one dir here, may be several); omit the hint when that fails |
| 9 | The configure summary | one pasteable line | no | (c) | `config/summary.cmake:19-24`; print the resolved PMIx path when row 4 produced one (Linux only, so not "every tree") |
| 10 | Documentation | `doc/mpi_support.md:30-49` describes the preprocess probe only | partially | (c) | user-facing → Codex language sweep on the touched section |
| 11 | `scripts/check_doc_claims.py` | CLAUDE.md's own claims about the compiler config | n/a | (c) | it reads `CLAUDE.md` and `doc/coding_philosophy.md` only (`check_doc_claims.py:50-52`); run it after R0 because CLAUDE.md describes the build flags, not as a gate for `mpi_support.md` |

### 3.1 Cross-cutting findings

- **Only `-rpath` is affected.** CMake already keeps implicit directories out of the `-L` list
  (`belfem.dir/link.txt` has no `-L/usr/lib64` although `BELFEM_RPATH` feeds `LINK_DIRECTORIES`
  at `finalize_compiler.cmake:99-102`); the shadowing is a NEEDED-resolution effect of `-rpath`.
- **R0 is the fix; the probe is the messenger.** Without R0 the probe reports case (b)
  correctly at cmake time and every build on such a machine stays broken.
- **The probe must carry the prefix `-rpath` entries or it lies in the false-failure
  direction** (row 2). It cannot carry the TPL entries; they never held a PMIx on any machine
  inspected, but that is an observation, not a proof (O2).
- **No project flags on the probe.** `config_gcc.cmake:64-72` documents `-Werror`; the flavor
  probe passes none and the link probe passes none: wrapper, rpath flags, source, output.
- **Probe artefacts** live under `${CMAKE_BINARY_DIR}/CMakeFiles/`, not the build root; the
  PMIx path is a normal variable, never cached (stale on reconfigure).

## 4. Ordered Steps

- [x] **R0** — Two sites, one rule each. (i) `detect_gcc.cmake:77-82`: skip the append when
  `BELFEM_TPATH` (already a REALPATH, `:70`) is `/usr` or `/`. (ii) `finalize_compiler.cmake`
  after `:86-88`: normalise `BELFEM_RPATH` (strip trailing slashes as `:21-28` does for
  includes, REALPATH each entry for the comparison only) and remove entries equal to
  `/usr/lib64`, `/usr/lib`, `/lib64`, `/lib`; guard every `list( REMOVE_ITEM … )` against an
  empty value list. `detect_icc.cmake:83` untouched. Gate: `belfem.dir/link.txt` rpath begins
  with `/opt/scls/mkl/lib64`; R6.2 shim links after R0, fails before it.
- [x] **R1** — Probe source `${CMAKE_BINARY_DIR}/CMakeFiles/belfem_mpi_link_probe.cpp`, a
  complete translation unit (`#include <mpi.h>`, `int main( int argc, char** argv )` calling
  `MPI_Init( &argc, &argv )` and `MPI_Finalize()`, `return 0`), written by `file( WRITE )`
  after the flavor detection, in the Open MPI and unrecognised branches alike.
- [x] **R2** — `execute_process( COMMAND ${CMAKE_CXX_COMPILER} [-Wl,-rpath,<list>] <src> -o
  <bin> )` with the R0-filtered `BELFEM_RPATH` colon-joined as a *user* argument (the wrapper
  appends its own flags after, `--showme:link`); omit the flag entirely when the list is
  empty. Capture `OUTPUT_VARIABLE` and `ERROR_VARIABLE` and concatenate (ld writes to
  stderr). Non-zero → `FATAL_ERROR` quoting the output; only then run the matcher:
  `PMIx_`, `libpmix\.so\.2[^\n]*not found`, or `cannot find libpmix\.so\.2` → prepend the
  package sentence (row 3) with the `ldd` hint when row 8 yields a path. No project flags.
  (after: R0, R1)
- [x] **R3** — Linux, `$SCLS` set, `EXISTS "${SCLSLIBDIR}/libpmix.so.2"`: `execute_process(
  COMMAND ldd <bin> )`; a failing or absent `ldd` is a silent skip, never fatal (the link
  already succeeded). Three shapes: no `libpmix.so.2` line → skip (row 6); `=> not found` →
  WARNING, no REALPATH; a path → REALPATH both it and `$ENV{SCLS}`, compare with a trailing
  slash on the prefix, WARNING when outside (text: what `ldd` showed, and that the loader
  reads `LD_LIBRARY_PATH` before RUNPATH; no claim about which copy was linked). Set
  `BELFEM_PMIX_LIBRARY` (normal variable, never cached) whenever a path was seen. (after: R2)
- [x] **R4** — `config/summary.cmake`: `PMIx       : <path>` under MPIHOME when
  `BELFEM_PMIX_LIBRARY` is set. (after: R3)
- [x] **R5** — `doc/mpi_support.md`: "At configure time" gains the link probe and the
  behaviour change; "At link time" gains the implicit-directory rule; "Checking what you have"
  gains `ldd $(mpicxx --showme:libdirs)/libmpi.so | grep pmix`. Codex language sweep on the
  touched sections only. (after: R4)
- [x] **R6** — Gates, all scratch-directory, no root, `/opt/scls` untouched:
  1. *Positive:* clean reconfigure on this machine (mkl flavor); summary shows the PMIx line;
     `link.txt` rpath begins with `/opt/scls/mkl/lib64`; `make banner` links.
  2. *Real-ld negative for the probe:* a `mpicxx` shim first on `PATH` that execs the real
     wrapper with `-Wl,-rpath,<stub dir>` prepended (the shim passes `-dumpversion` and
     `-showme:*` through, so `detect_gcc.cmake:109-125` still passes); configure must stop
     with the package sentence. The stub is the one from the reproducer (24 classic symbols).
     A prefix copy is **not** a valid negative: the wrapper bakes `-L`/`-rpath` to the install
     libdir (`mpicxx --showme:link`).
  3. *Absent-library negative:* not reproducible with a real ld on a healthy machine, because
     `libmpi.so`'s own RUNPATH always finds the SCLS copy after every `-rpath` miss. The
     `not found` branch of the matcher is covered by gate 5 only, and the plan says so.
  4. *Row 4 warning:* `LD_LIBRARY_PATH=<stub dir>:$LD_LIBRARY_PATH cmake …` (prepended, not
     replaced) on a healthy tree; the probe links (rpath wins at link time), `ldd` shows the
     stub, WARNING appears.
  5. *Matcher unit:* a shim that forwards `-dumpversion`, `-E`, `--showme*` and `-showme*` to
     the real wrapper (otherwise `detect_gcc.cmake:118-125` and the `-E` probe fail first) and,
     on a real link (`-o` present), prints one canned line and exits 1. Run it twice: with
     `undefined reference to PMIx_Info_create` and with `libpmix.so.2, needed by …, not found`;
     both fatal texts must contain the package sentence.
  (after: R5)
- [x] **R7** — `scripts/check_doc_claims.py`; jury round on the restricted diff
  (`tmp/ai_exchange/review_mpi_link_probe.diff`, own files only, shared checkout). Devlog.
  (after: R6)

### 4.0 Implementation Progress (updated 2026-09-18)

**Implemented, gates run, code jury pending.** R0 at both sites (`detect_gcc.cmake` skip for a
`/usr` or `/` prefix; `belfem_prune_system_libdirs()` in a new `config/scripts/` file, called
from `finalize_compiler.cmake` and from the probe). R1–R4 in `find_mpi.cmake` and
`config/summary.cmake`; one `-Wl,-rpath,<dir>` flag per directory (ld64 takes no colon list).
R5 in `doc/mpi_support.md`, language sweep running. R6 gates 1, 2, 4, 5 green in scratch build
directories on the mkl flavor; gate 3 is unit-only by design. `check_doc_claims.py` 38/38.
Not run: `make` (Christian builds), Darwin, Intel.

## 5. Open Design Questions

- **O1 — Which flavors are exempt from the "PMIx outside `$SCLS`" warning?** RESOLVED
  2026-09-18 → (b), installation-derived: warn only when `$SCLS/lib*/libpmix.so.2` exists.
  Raised by Codex and Grok; the name-based rule made the warning text false whenever the SCLS
  PMIx was absent and the system one was new enough to link.
- **O2 — Should the probe carry the full final rpath list?** Keep it in `find_mpi.cmake` for
  fail-fast (Grok R-a); the TPL directories are a residual, not a proven non-risk.
- **O3 — Fatal or warning when the link probe fails?** RESOLVED 2026-09-18, Christian → fatal.
  The argument, not the vote: a warning postpones the same failure to the first executable.
- **O4 — Match `PMIx_` only, or also `libpmix.so.2 … not found`?** RESOLVED 2026-09-18 →
  both, folded into R2 (Grok C10).
- **O5 — Where does the implicit-directory filter live?** (a) once in
  `finalize_compiler.cmake` (one site, catches every future append, but the probe in
  `find_mpi.cmake` runs earlier and must repeat it); (b) at the append sites in
  `detect_gcc.cmake` / `detect_icc.cmake` (the only source of a system directory today, probe
  needs nothing). RESOLVED 2026-09-18, Christian → (b) plus the defensive remove in (a).
  **Rule narrowed after round 2 (Codex F1, Grok C1–C3, independently):** not "anything in the
  implicit list" but "the compiler prefix is `/usr` or `/`" at (b) and "exact canonical
  system libdirs" at (a). Same sites as decided; reported to Christian as a change.

## 6. Interface Design

| variable | set where | consumed where |
|---|---|---|
| `_BELFEM_MPI_LINK_PROBE_SRC/BIN` | `find_mpi.cmake` | `find_mpi.cmake` only |
| `BELFEM_PMIX_LIBRARY` | `find_mpi.cmake` (R3), not cached | `config/summary.cmake` (R4) |

Message texts (final wording audited with the diff):

```
Could not link an MPI program with <wrapper> ( <flavor> ):
<linker output>
[Open MPI's PMIx dependency is missing or incompatible with this Open MPI.
 Under SCLS, install scls-<flavor>-pmix. Check: ldd <libdir>/libmpi.so | grep pmix]
Background: doc/mpi_support.md
```

```
ldd resolves libpmix.so.2 to <path>, not to the copy in <prefix>. The loader reads
LD_LIBRARY_PATH before RUNPATH, so a run in this environment may load that PMIx.
```

## 7. Definition-of-Done Checklist

- [x] Every gap-table row mapped to a step or an open question (rows 0–11 → R0–R7, O2).
- [x] Each claimed gap backed by a citation or a reproducer, not assumption (round 2).
- [x] Ordered steps with dependencies.
- [x] O3 and O5 decided by Christian, recorded in place (2026-09-18).
- [x] R6 gates 1, 2, 4, 5 run and recorded in the devlog (gate 3 unit-only by design).

## 8. Audit Trail

- Exchange thread: `tmp/ai_exchange/review_mpi_link_probe.md`. Round 1 (plan, jury, Codex
  gpt-5.6-terra/high, Grok grok-4.6/high): Codex raised the probe-is-not-the-same-link point
  and the unreliable libmpi path; Grok raised the O1/R3 inconsistency (also Codex), the O4
  wiring gap, the invalid prefix-copy reproducer, and the fatal-on-preprocess-miss change;
  Grok's ld search-order claim was refuted against `man ld`. The rpath-order defect (R0) was
  found by Claude during citation verification and backed by the stub reproducer before it
  entered the plan. Round 2 (revised plan, jury, both xhigh): Codex and Grok independently
  refuted the implicit-directory rule for R0 (a prefix compiler's libdir is implicit too);
  rule narrowed to the canonical system libdirs at the same two sites. Grok supplied the
  robustness list for R1–R3 and the gate-5 forwarders; Codex the multi-directory
  `--showme:libdirs` and the CMake compiler-swap caveat (out of scope). All citations
  re-verified against the tree.
