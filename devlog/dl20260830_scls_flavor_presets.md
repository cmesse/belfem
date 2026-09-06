# SCLS Flavor Presets, and the `-lcblas` Link Blocker

**Date:** 2026-08-30
**Purpose:** Record the `-lcblas` fix that unblocked the test-suite link on the `debug` flavor, the
two reframes it triggered in `todo/linalg_backend_verification.md`, and the SCLS-flavor preset that
landed as a result.
**Module:** `config/linalg`, `config/system`, `CMakeLists.txt` (build system only)
**AIs involved:** Claude (exploration, fix, plan revisions, implementation)

---

## 1. The blocker

`make check` against `/opt/scls/debug` died at the first executable link:

```
/usr/bin/ld: cannot find -lcblas
```

`config/linalg/config_mkl.cmake` emitted `-llapack -lcblas -lblas` unconditionally in its netlib
branch. Of the three SCLS flavors only `gcc` ships a `libcblas.so`, and there it is a symlink to
`libopenblas.so` — the same file `libblas.so` and `liblapack.so` point at, i.e. a redundant third
alias. `/opt/scls/debug` is reference netlib and ships no CBLAS at all, so the flag had nothing to
resolve against and the link search path (`-L/opt/scls/debug/lib64` only) had nowhere else to look.

**The flag was never providing a symbol on this path.** Measured before removing it: `libbelfem.a`
and the `test_physics` objects carry zero undefined `cblas_` references; Armadillo is built
`ARMA_USE_WRAPPER` against `libblas.so`/`liblapack.so` and calls the Fortran symbols; the debug
flavor's `libblas.so`/`liblapack.so` export no `cblas_` symbol; and the debug `libslate.so.2.0.0`
has no undefined `cblas_` either, reaching BLAS indirectly through `libblaspp`/`liblapackpp`.

This was **already an adjudicated item** — D1 in `todo/linalg_backend_verification.md`, closed by
both auditors on link *and* runtime evidence, with the decision being to drop the flag rather than
conditionally link it. A first attempt at a `find_library` probe was backed out for exactly the
reason the plan gives: `find_library` caches, so it would go stale across a `$SCLS` repoint, and the
`gcc` flavor does not need the alias either. Landed as the decided form: the netlib branch now sets
`"-llapack" "-lblas"`, with a comment recording why CBLAS is not a provider here.

## 2. Two reframes, both Christian's

The plan the fix came out of was a much larger piece of work — probe-based *autodetection* of the
backend, round-1 audited by Codex (`terra`/`high`) and Grok (`grok-4.6`/`high`), both returning
blocking findings. It was cut down twice on the same day.

**First: no autodetection.** The knob selects, the probe verifies, a miss is a configure-time error.
That withdrew the tri-state `AUTO` design (R1) and dissolved **D4, D6, D7, D8, O1, O2 and O5** — most
of round 1's output, including its most serious finding, since D4's hazard lived entirely inside the
three-way `STREQUAL` chain that `AUTO` required. None of it was wasted: every one of those was a real
defect in the design as it then stood, and D4 is part of why the simpler design is attractive.

**Second: read `$SCLS_FLAVOR`.** The flavor already states which toolchain was sourced, so the
configure should read it rather than infer anything. This is what the earlier design was circling
without naming, and it made the plan smaller again — O6 then closed with *no candidate ladder and no
vendor knob*, because the netlib names already resolve to whatever a flavor ships.

**The distinction that keeps this from being the autodetection that was just rejected**, written
into the plan's scope guards so a later session does not read it as a flip-flop: a preset reads a
*declaration* and is a default any `-D` overrides; autodetection searches the *host* and lets the
answer vary by machine. The test is whether it reads a statement or goes looking.

Method note recorded in the plan's audit trail: **both round-1 audits were thorough and both were
scoped to the wrong question.** The brief took the knobs as fixed and asked only how to resolve a
backend behind them. A design round asking "what does the environment already tell us?" before "how
do we detect it?" would have reached the smaller plan first.

## 3. What landed

`config/system/find_scls_flavor.cmake` (new, 62 lines), included at `CMakeLists.txt:74` — inside the
"User Settings" block and **above both `option()` calls**. That placement is the whole trick:
`option()` sets a default only when the cache entry does not exist yet, so a default computed after
the declaration does nothing, and does it silently, because the first configure simply caches the
un-preset value. `find_scls.cmake` runs at `:163`, far below the `option()` calls at `:79` and
`:98`; only the environment read moved up, and `BELFEM_CONFIG_DIR` moved to `:69` to make the
include possible (its later duplicate removed).

| `$SCLS_FLAVOR` | BLAS/LAPACK it provides | `USE_MKL` | `USE_DEBUG` |
|---|---|---|---|
| `debug` | reference netlib, real `libblas.so`/`liblapack.so` | OFF | **ON** |
| `gcc` | OpenBLAS — those names are symlinks to `libopenblas.so` | OFF | OFF |
| `mkl` | none of its own; oneAPI MKL at `$MKLROOT` | **ON** | OFF |
| *(no `$SCLS`)* | whatever the devel packages installed | OFF | OFF |

`$SCLS` set but `SCLS_FLAVOR` unset derives the flavor from the basename (the flavors live at
`/opt/scls/<flavor>`). A flavor outside the three — `cea`, or a private prefix — presets nothing and
says so with `STATUS`, deliberately **not** a `FATAL_ERROR`, since pointing `$SCLS` at a hand-rolled
prefix is legitimate and `find_scls.cmake` already aborts on a prefix with no `lib`. A disagreement
between `$SCLS_FLAVOR` and the basename warns naming both and trusts `$SCLS_FLAVOR` as the more
specific statement. `config/summary.cmake` prints `SCLS flavor:` beneath the existing `SCLS:` line.

`CLAUDE.md` updated in the same session, because **the `USE_DEBUG` default flips ON → OFF for anyone
without SCLS**: a fresh non-SCLS clone now builds Release where it used to build Debug. That is the
one deliberate compatibility break here.

## 4. Verified, and what is not

**Verified** — ten cases against a standalone CMake project in the scratchpad, each a fresh
configure unless stated: the three flavors give the table above; no `$SCLS` gives OFF/OFF silently;
basename derivation works with `SCLS_FLAVOR` unset; the mismatch case warns and follows
`SCLS_FLAVOR`; `cea` and a private prefix both preset nothing without aborting; `-DUSE_MKL=OFF`
beats an `mkl` preset; and a tree configured as `debug` then reconfigured as `mkl` keeps its cached
`USE_MKL=OFF`/`USE_DEBUG=ON` while the flavor line reads `mkl`.

That last case is worth keeping in mind: `option()` cannot overwrite a populated cache, so
**repointing `$SCLS` under an existing build tree does not change its knobs.** A build tree is
effectively per-flavor. It is not supported so much as *caught* — the summary now shows the
inconsistency, and once the probe work lands the mismatch becomes a configure error rather than a
link-stage surprise.

**Not verified.** None of this has run inside BELFEM's own configure — the flavor logic was exercised
standalone, and the integration into `CMakeLists.txt` is reviewed only. One `cmake` per flavor in a
real tree is owed, and builds are Christian's to run. The `-lcblas` removal likewise has not been
seen to link; `config_mkl.cmake` is in `CMakeFiles/Makefile.cmake`, so `make check` re-runs the
configure by itself.

`scripts/check_doc_claims.py` is green at 37/37 after every edit. Note that it can no longer see the
`USE_DEBUG` default at all: `:135` matches `option( NAME "…" ON|OFF )` against a literal, and the
default is now `${BELFEM_DEFAULT_USE_DEBUG}`. `USE_TEST` and `USE_VTK` are unaffected. Pinning a
*conditional* default needs a different check; recorded as a known gap in the plan rather than a
silent one.

## 5. Owed

- One real `cmake` per flavor plus a link, per the plan's R7 matrix.
- The probe work itself — R2 (`belfem_probe_link()`) through R7 — has not started. Its surviving
  round-1 findings are D2, D3, D5, D9, D10 and half of D11.
- **An audit round on the landed CMake.** Comment and register text needs no round; executable code
  does, and CMake is executable. Not run here.

## 6. Files touched

| file | change |
|---|---|
| `config/linalg/config_mkl.cmake` | netlib branch drops `-lcblas` |
| `config/system/find_scls_flavor.cmake` | new — flavor read and presets |
| `CMakeLists.txt` | `BELFEM_CONFIG_DIR` to `:69`, flavor include at `:74`, both `option()` defaults from the presets, duplicate `BELFEM_CONFIG_DIR` removed |
| `config/summary.cmake` | prints the flavor |
| `CLAUDE.md` | per-flavor default table replaces "`USE_DEBUG=ON`, the default" |
| `todo/linalg_backend_verification.md` | renamed from `linalg_backend_autodetect.md`; rewritten across both reframes; R0/R0b ticked |
| `todo/README.md` | entry retitled and rewritten |
