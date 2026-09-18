# 2026-09-18 — System libdir on every rpath let a distro PMIx shadow the SCLS one; configure-time MPI link probe

**Date:** 2026-09-18
**Purpose:** Record the diagnosis of a 0.9.1 link failure on a colleague's machine, the defect
it exposed in BELFEM's own rpath assembly, the fix, the new configure-time probe, and the
three-round jury that shaped both.
**Module:** `config/compiler`, `config/scripts`, `config/system`, `doc/mpi_support.md`
**Plan:** `todo/closed/mpi_link_probe_plan.md`

## The report

A 0.9.1 build (`gcc` flavor, `/opt/scls/gcc`) compiled every module and then failed at every
executable link with 38 lines of the form

```
/usr/bin/ld: /opt/scls/gcc/lib64/libmpi.so: undefined reference to `PMIx_Info_create'
```

Christian's question was whether SCLS ships PMIx at all. It does: `scls-<flavor>-pmix`
(PMIx 5.0.10), required by `scls-<flavor>-openmpi`, which is configured with
`--with-pmix=%{prefix}`; `libmpi.so` and `libopen-pal.so.80` carry `NEEDED libpmix.so.2` and
`RUNPATH <prefix>/lib`.

## Diagnosis by symbol set

Open MPI 5.0.10 imports 62 PMIx symbols (`nm -D --undefined-only` over the two libraries).
The colleague's linker satisfied 24 of them, the classic API (`PMIx_Init`, `PMIx_Get`,
`PMIx_Fence`, `PMIx_Spawn`, …), and failed on exactly the 38 that became functions in PMIx 4.2
plus one 5.x internal. A missing library leaves all 62 unresolved and prints
`needed by … not found`. So an *older* `libpmix.so.2` was found, not none. EL9 ships a `pmix`
RPM at 3.2.x in `/usr/lib64`.

## The BELFEM defect, found during citation verification of round 1

Every executable link on this machine carried

```
-Wl,-rpath,/usr/lib64:/opt/scls/mkl/lib64:/opt/intel/oneapi/mkl/latest/lib/intel64:…
```

`detect_gcc.cmake` records the compiler prefix's libdir on `BELFEM_RPATH` first, and with the
system GCC that prefix is `/usr`. GNU ld resolves a shared library's own `DT_NEEDED` entries
through the `-rpath` directories in order, before `LD_LIBRARY_PATH`, before the library's own
RUNPATH, before the default directories (`man ld`, `-rpath-link`, items 2–7). A distro PMIx in
`/usr/lib64` therefore beat an installed SCLS PMIx on every link.

**Reproducer** (scratch directory, `/opt/scls` untouched): a stub `libpmix.so.2` defining only
the 24 classic symbols; `mpicxx -Wl,-rpath,<stub>:/opt/scls/mkl/lib64 hello.cpp` fails with the
identical 38-symbol set (`diff` against the colleague's list empty); reversed order links and
`ldd` shows the SCLS copy; the plain wrapper links. The colleague's box may therefore have
nothing missing at all. A handout (`tmp/pmix_diagnostics_for_colleague.md`, ephemeral)
collects the environment and distinguishes "package absent" from "rpath order".

## What landed (uncommitted, Christian builds and commits)

| file | change |
|---|---|
| `config/compiler/detect_gcc.cmake` | skip the compiler-libdir append when the prefix is `/usr` or `/` |
| `config/scripts/belfem_prune_system_libdirs.cmake` (new) | remove `/usr/lib64`, `/usr/lib`, `/lib64`, `/lib` and the multiarch dir from a list, REALPATH both sides |
| `config/compiler/finalize_compiler.cmake` | call the prune after the dedupe, before `CMAKE_BUILD_RPATH` / `CMAKE_INSTALL_RPATH` are set |
| `config/system/find_mpi.cmake` | link a two-call MPI program with `mpicxx` after the flavor probe, carrying the pruned rpath collected so far (one `-Wl,-rpath,<dir>` per entry, ld64-safe); `FATAL_ERROR` with the linker output; for Open MPI, a PMIx-shaped failure adds the package name (`scls-<flavor>-pmix`) and an `ldd` hint on a verified `libmpi.so`; on Linux under SCLS with a prefix PMIx, `ldd` on the probe warns when `libpmix.so.2` resolves outside the prefix |
| `config/summary.cmake` | `PMIx       : <path>` line |
| `CMakeLists.txt` | include the helper |
| `doc/mpi_support.md` | configure-time probe, the rpath rule, the `ldd` one-liner; Codex language sweep (luna/medium) applied |

Behaviour change, stated: a wrapper that cannot link `MPI_Init`/`MPI_Finalize` now stops the
configure. The earlier "preprocess miss warns and continues" path still warns first.

## Gates (scratch build directories, mkl flavor; `make` not run)

| gate | result |
|---|---|
| clean configure | rc 0, no warnings, rpath begins `/opt/scls/mkl/lib64`, summary shows `PMIx : /opt/scls/mkl/lib/libpmix.so.2` |
| PATH shim injecting the stub dir into real links | rc 1, 38 undefined references, package sentence and `ldd` hint present |
| `LD_LIBRARY_PATH` prepended with the stub | rc 0, WARNING names the stub, summary shows it |
| canned-output shim, `undefined reference` and `not found` forms | rc 1 each, package sentence present |
| `scripts/check_doc_claims.py` | 38/38 |

Not run: `make` on any tree, Darwin, Intel. The build gate is owed by Christian.

## Audit trail (`tmp/ai_exchange/review_mpi_link_probe.md`, three rounds)

- **Round 1, plan (Codex terra/high, Grok 4.6/high).** Codex: the probe is not the same link
  as the target's; no reliable `libmpi.so` path for the message. Codex and Grok: R3 silently
  decided O1. Grok: R2 lacked the `not found` match; the prefix-copy reproducer was invalid
  (the wrapper bakes its own `-rpath`); a Grok claim on ld's search order was refuted against
  `man ld`. The rpath defect above was Claude's, from verifying Codex's citation of
  `link.txt`, reproducer-backed before it entered the plan.
- **Round 2, revised plan (both xhigh).** Codex and Grok independently refuted the R0 rule
  "remove every implicit link directory": a prefix compiler's libdir is implicit too and must
  stay. Rule narrowed to the canonical system libdirs at the same two sites Christian had
  decided (O5). Grok supplied the robustness list (stderr capture, empty rpath list, `EXISTS`
  on `SCLSLIBDIR`, three-way `ldd` parse, shim forwarders).
- **Code round (terra/high, 4.6/high).** Codex and Grok: the multiarch directory escaped the
  prune (package finder searches `lib/<arch>`, HDF5 config appends it). Codex: PMIx hint not
  gated on Open MPI; `ldd` hint unchecked. Grok: comment and doc overclaimed "same
  directories", doc date and Related stale, literals not REALPATHed, `PARENT_SCOPE` quoting,
  EOF newline. All fixed; gates re-run green. Residuals accepted: English-locale matcher loses
  only the hint; the probe carries no project link flags by design.

Decisions by Christian: O3 fatal; O5 sites; 0.9.2 requirement is that BELFEM prefers the SCLS
PMIx whenever it is present.

## Residue

- Build gate on both trees (Christian).
- `CMAKE_CXX_COMPILER` is swapped to `mpicxx` after `project()`, which CMake documents as
  unsupported; the cache still says `g++`. Pre-existing; a toolchain file would be the clean
  fix. Not opened as a todo; recorded here.
- Darwin and Intel configures unexercised.
