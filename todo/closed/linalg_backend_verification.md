# Build Configuration from the SCLS Flavor, and Configure-Time Verification of the Backend

> **CLOSED 2026-09-04** (scope jury — Codex terra/high + Grok 4.6/high, Christian's ruling): the
> flavor presets (R0/R0b) are the whole deliverable, extended the same day by `USE_PARDISO` in
> parity with `USE_MKL`; a flavor whose name contains `mkl` presets both ON, everything else OFF.
> The configure-time link probe (R2–R5), its helper and R6 are **declined**, and so is any
> flavor-versus-knob refusal (tried and withdrawn the same day): any vendor other than MKL is
> assumed to behave like reference BLAS/LAPACK, so there is nothing to verify beyond the
> declaration, and an explicit `-D` or ccmake edit always wins. Owed gate: one fresh configure per
> flavor, Christian's to run. Record: `tmp/ai_exchange/review_linalg_backend_scope.md`
> (distilled into `devlog/dl20260904_linalg_backend_scope_jury.md`).

**Date:** 2026-08-30
**Purpose:** two halves of one problem, both settled on 2026-08-30.

1. **Preset the build knobs from the SCLS flavor.** `$SCLS_FLAVOR` already says which toolchain the
   user sourced — `debug`, `gcc` or `mkl` — so the configure should read it instead of making the
   user restate it. When `$SCLS` is set, the flavor supplies the *defaults* for `USE_MKL` and
   `USE_DEBUG`. When it is not, BELFEM behaves like any other Linux program: look for installed
   devel packages, `USE_MKL=OFF`, `USE_DEBUG=OFF`, and the user sets what they want.
2. **Verify whatever the knobs ended up saying.** Today both branches of `config_mkl.cmake` emit
   hardcoded `-l` flags that nothing ever checks, so a missing backend surfaces as an
   unresolved-symbol wall at the final link, after a full compile. Make the configure link a probe,
   and `FATAL_ERROR` with a useful message when the selected backend is not there.

**The backend is never guessed.** A preset reads a *declaration* — the flavor the user sourced — and
it is only ever a default that `-DUSE_MKL=…` overrides. Nothing searches the host to decide what
BELFEM should link against, and nothing falls back from one backend to another. The two halves are
complementary: the preset makes the common case need no flags at all, and the probe catches the case
where the preset and the environment have drifted apart.
**Module:** `config/linalg`, `config/system` (build system only; no `src/` change)
**AIs involved:** Claude (exploration + plan), Codex (audit), Grok (third voice)
**Status:** **CLOSED 2026-09-04** — see the banner. Before that: PLAN — round-1 audited by both vendors, amended, then **rescoped by Christian's
decision of 2026-08-30: no autodetection.** The tri-state `AUTO` knob is withdrawn; `USE_MKL` stays
the boolean it is today and the work reduces to *verify what was selected, fail loudly when it is
missing*. That decision dissolves R1, D4, D6, D7, D8, O1 and O5 outright (all were consequences of
AUTO) and narrows O2 and O4 — see §4.1 and §5, where each is struck with its reason rather than
deleted. The probe work (R2) and the error text (R5) are unaffected and remain the core of the task.
**One build-system change has landed ahead of the plan: the D1 `-lcblas` removal at
`config_mkl.cmake` (2026-08-30), pulled forward because it blocked the test-suite link on the
`debug` flavor. Nothing else is modified, and no `src/` file is touched.**
**O6 decided 2026-08-30 (Christian): no vendor knob and no candidate ladder** — under `USE_MKL=OFF`
the configure keeps looking for BLAS/LAPACK exactly as it does today (`-llapack -lblas`) and merely
checks the result. That dissolves O3 and half of D11, and resolves O4 by construction (§5).

**R0 and R0b LANDED 2026-08-30** — the flavor presets are in the tree and verified against a
standalone CMake project across ten cases (see R0); what is still owed is one real `cmake` per
flavor inside BELFEM itself, which is Christian's to run. The probe work (R2–R7) has not started.

**Reframed the same day (Christian) with the SCLS-flavor presets, now R0 and §6.** This is what O4
was circling without naming: `$SCLS_FLAVOR` is the declaration the configure should read, so the
flavor sets the *defaults* for `USE_MKL` **and** `USE_DEBUG`, and a tree with no `$SCLS` behaves like
any ordinary Linux build. It makes the plan smaller, not larger — most of what the earlier design
tried to infer is simply stated by the environment. **Two consequences need Christian's eye before
R0 lands, both in §6.2: the `USE_DEBUG` default flips from ON to OFF for anyone without SCLS (a
user-visible change, and `CLAUDE.md` documents the current default), and the presets must be
computed before `CMakeLists.txt:71`/`:90`, which is 68 lines earlier than where `find_scls.cmake`
runs today.** No open design questions otherwise; implementation order is R0 → R2 → R3/R4 → R5 → R6
→ R7.

> **Scope guards:**
> - **No autodetection. The knob selects, the probe verifies, a miss is a configure-time error**
>   (decided 2026-08-30, Christian). `USE_MKL` true means MKL or `FATAL_ERROR`; false means a
>   BLAS/LAPACK distribution or `FATAL_ERROR`. Neither value may fall back to the other backend.
>   This supersedes the round-1 design, in which an unset knob meant "prefer MKL, fall back" — the
>   probe is a *check on a decision already made*, not a way of making it.
> - **A flavor preset is not autodetection, and the distinction is the load-bearing one in this
>   plan.** `$SCLS_FLAVOR=mkl` is the user stating which toolchain they sourced; presetting
>   `USE_MKL=ON` from it reads that statement. Autodetection would have searched the *host* — asking
>   what happens to be installed — and let the answer vary with the machine. Concretely, the preset
>   is deterministic given the environment, is only ever a **default** that any explicit `-D`
>   overrides, and never chooses between two backends that are both present. If a future session
>   cannot tell whether something is a preset or a detection, the test is: *does it read a
>   declaration, or does it go looking?*
> - **All three `/opt/scls` build flavors must configure.** `gcc` and `debug` under `USE_MKL=OFF`,
>   `mkl` under `USE_MKL=ON`, selected by pointing `$SCLS` at the flavor. (`/opt/scls/cea` is *not*
>   a flavor — it holds CEA thermodynamic input decks and no `lib`/`include`, so `find_scls.cmake:4-20`
>   fails there before backend detection ever runs. Corrected after Codex round 1.) `$SCLS/lib64` is therefore the
>   *first* search root, ahead of system paths — mixing an SCLS-built MPI with a system BLAS is the
>   toolchain mixing `find_mpi.cmake:9-11` already warns about.
> - **Build system only.** No file under `src/` is touched. `BELFEM_MKL` / `BELFEM_NETLIB` keep
>   their present meaning and remain the only compile-time backend defines.
> - **Presets fire only when `$SCLS` is set.** With no `$SCLS` there is nothing to read and nothing
>   to preset: `USE_MKL=OFF`, `USE_DEBUG=OFF`, devel packages from the default paths, and the user
>   sets the knobs. That path must stay the plain-Linux path a stranger cloning the repository
>   expects.
> - **The probe gates the linear-algebra subsequence of the link line, nothing more.** It does not
>   and must not claim to validate the full link (PETSc, MUMPS, STRUMPACK, HDF5 all come later and
>   are configured after this point). Failure of any of those remains a link-stage failure.
> - **Not in scope:** replacing the rank-based link assembly in `finalize_compiler.cmake`; adopting
>   CMake's own `FindBLAS`/`FindLAPACK`; touching `USE_MKL_64BIT_API` semantics.
> - **Compatibility promise, and the one place it is deliberately broken.** `-DUSE_MKL=ON` and
>   `-DUSE_MKL=OFF` keep exactly today's meaning, an existing tree's cached values survive
>   untouched (a preset sets a default and cannot overwrite a populated cache), and a backend that
>   cannot be linked is now reported at configure time instead of at the final link. **The
>   exception, stated plainly because it is user-visible: `USE_DEBUG` defaults to OFF without SCLS,
>   where it defaults to ON today** (`CMakeLists.txt:90`). A fresh non-SCLS clone therefore builds
>   Release where it used to build Debug. Decided 2026-08-30 (Christian); see §6.2 for what else
>   that obliges us to change.

---

## 1. Current Behaviour and How It Fails

`USE_MKL` defaults **OFF** (`CMakeLists.txt:71`). Both branches of
`config/linalg/config_mkl.cmake` emit link flags that are never checked against the filesystem or
the linker.

The backend is supplied by the SCLS flavor `$SCLS` points at, since `find_scls.cmake:7` puts
`$SCLS/lib64` on the link path via `link_directories()`. Measured on this host, `USE_MKL=OFF` today:

| `$SCLS` flavor | result | reason |
|---|---|---|
| `/opt/scls/gcc` | links | ships `libcblas.so`; `libblas.so`/`liblapack.so` are symlinks to `libopenblas-r0.3.33.so` |
| `/opt/scls/debug` | **fails at the final link** | `ld: cannot find -lcblas`. Dropping `-lcblas` from the same command links clean — that flag is the *only* failure |
| `/opt/scls/mkl` (current `$SCLS`) | **fails** | the flavor ships no BLAS/LAPACK at all; MKL is the only backend available there |

| Failure | Mechanism | Evidence |
|---|---|---|
| The debug SCLS flavor cannot be built with `USE_MKL=OFF` | The branch hardcodes `-lcblas`, which exists in exactly one of the four flavors | `config_mkl.cmake:74`; `ls /opt/scls/*/lib64/*cblas*` returns only the `gcc` flavor. Link probes above |
| Nothing tells the user which flavor is usable | Neither branch checks any of the `-l` flags it emits | `config_mkl.cmake:41-63`, `:74` |
| A trimmed oneAPI install configures clean and dies at the final link | `find_mkl.cmake` checks only that `$MKLROOT` is an existing **directory**; the five-to-six `-l` flags built from it are never checked | `config/system/find_mkl.cmake:12`; libs assembled at `config_mkl.cmake:41-63`. ScaLAPACK and the BLACS layer are separately deselectable in the Intel installer, so this is the common failure, not an exotic one |
| The MKL library directory is chosen by existence, not by content | `IS_DIRECTORY "${BELFEM_MKLROOT}/lib/intel64"` picks `intel64` whenever the directory exists, even if the libraries live only in `lib` | `config_mkl.cmake:12-17` |
| `-lscalapack` under the netlib branch is unchecked | Emitted at rank 5 whenever `USE_MPI`; resolves on this host only because `find_scls.cmake` calls `link_directories($SCLS/lib64)` | `config_mkl.cmake:70-72`, `config/system/find_scls.cmake:7` |
| The wrong BLACS flavour is emitted on a non-Open-MPI build | `-lmkl_blacs_openmpi_${SUFFIX}` is unconditional | `config_mkl.cmake:62`. Already flagged in `CLAUDE.md`; a link probe catches it as a side effect |

**Bottom line:** the configure never asks the linker whether the backend it just chose actually
resolves, so the answer arrives after a full compile as a wall of undefined references — and which
flavors happen to answer "yes" today is an accident of which `-l` names each one ships.

## 2. Architecture: Why a Link Probe, Not a File Search

The gate must be the *same question the final link asks*: do `dgemm_`, `dgetrf_` (and, under MPI,
the ScaLAPACK/BLACS entry points) resolve against this candidate link line, in this order, with
these linker flags? A `find_library` existence check answers a weaker question and gets two cases
wrong that we care about:

- a `liblapack.so` that exists but is missing the ScaLAPACK layer;
- an MKL tree whose `lib/` holds the sequential libraries but not `mkl_blacs_openmpi_lp64`.

**Admitted gap — the probe cannot detect an integer-width mismatch.** `dgemm_` and `dgetrf_` have
identical symbol names under LP64 and ILP64, so a probe that links them proves nothing about
`USE_MKL_64BIT_API` agreement between BELFEM and the backend. An earlier draft of this section
implied otherwise; that was wrong (Codex round 1). Interface-width validation is out of scope and
stays a link-or-runtime failure.

`try_compile` with the assembled link line answers the real question. This also follows the
precedent already set in this tree: `config/system/find_mpi.cmake:41-58` decides the MPI flavour by
preprocessing `mpi.h` rather than by guessing from paths, precisely because a path is not evidence.

**Rejected alternative — CMake's `FindBLAS`/`FindLAPACK`.** They would replace the hand-assembled
MKL block wholesale, including the `--no-as-needed` / interface-layer / threading-layer ordering
that `config_mkl.cmake:38-63` documents as following the Intel link advisor, and they know nothing
about the rank-based assembly in `finalize_compiler.cmake`. The cost is a rewrite of a block that
works; the benefit is detection we can get in ~40 lines. Deferred, not dismissed.

**Rejected alternative — linking a versioned soname directly.** `find_library` will not find
`/usr/lib64/libopenblas.so.0`, and it should not: linking a build against a runtime package's
versioned soname produces a binary that breaks on the next package update. The correct outcome on
a host without `openblas-devel` is a clear error telling the user to install it. This is a
deliberate strictness, not a gap.

## 3. Gap Table

| # | State | Needed for | Handled today? | Class | Citation / rationale |
|---|---|---|---|---|---|
| 1 | Is MKL present and complete? | verifying `USE_MKL=ON` | no — directory existence only | (c) | `find_mkl.cmake:12` |
| 2 | Which MKL libdir holds the libraries | link line | by directory existence | (c) | `config_mkl.cmake:12-17` |
| 3 | Is a BLAS+LAPACK distribution present? | verifying `USE_MKL=OFF` | **not at all** | (c) | `config_mkl.cmake:74` |
| ~~4~~ | ~~Which system BLAS? (OpenBLAS / FlexiBLAS / netlib)~~ | ~~link line~~ | hardcoded netlib pair | **not a gap** | Decided 2026-08-30 (O6): the netlib *names* are the interface and already resolve to whatever a flavor ships — OpenBLAS on `gcc`, reference netlib on `debug`. Nothing to choose. See R4 |
| 5 | Is ScaLAPACK present under MPI? | link line | no | (c) | `config_mkl.cmake:62,71` |
| ~~6~~ | ~~Distinguishing "user chose OFF" from "nobody chose"~~ | ~~AUTO semantics~~ | impossible — `option()` is boolean | **not a gap** | The distinction is only needed by autodetection, withdrawn 2026-08-30. Unset means OFF, as it always has. See §6 |
| 7 | Is `-lcblas` needed at all? | netlib link line | emitted unconditionally; present in 1 of 4 flavors | (c) → **D1** | measured blocker on the `debug` flavor; see §3.1 |
| 8 | What the configure reports | diagnosis | prints only "MKL" or "BLAS and LAPACK" | (c) | `config/summary.cmake:28-33` |
| 9 | Which SCLS flavor is in use | defaulting `USE_MKL`/`USE_DEBUG` | **not read at all** — `$SCLS_FLAVOR` is never referenced in any CMake file | (c) → **R0** | `find_scls.cmake` reads `$SCLS` only, and not until `CMakeLists.txt:158` |
| 10 | `USE_DEBUG` matching the flavor's own libraries | build/TPL consistency | no — defaults ON everywhere, including against release-built flavors | (c) → **R0** | `CMakeLists.txt:90`; see §6.2a for the compatibility break this fixes |

### 3.1 Cross-cutting findings

- **`-lcblas` is the confirmed blocker (D1), not a tidy-up — but not for the reason first given.**
  It is emitted unconditionally (`config_mkl.cmake:74`) yet exists in only one build flavor, and
  removing it is sufficient to make the `debug` flavor link (measured, §1).

  The draft claimed "no cblas consumer exists in the tree". **Codex refuted that (round 1), and the
  refutation is confirmed by independent measurement:** `/opt/scls/mkl/lib64/libslate.so.2.0.0`
  carries four undefined CBLAS symbols — `cblas_{s,d,c,z}gemm_batch` — and SLATE is linked
  unconditionally at rank 8 whenever `USE_STRUMPACK` is on, which is the default
  (`config/linalg/config_strumpack.cmake:12-38`, `CMakeLists.txt:99`). Blaze and Armadillo are
  indeed not consumers (`blaze_config.hpp:61-62` sets `BLAZE_BLAS_MODE 0` under `BELFEM_NETLIB`;
  `config_matrix.cmake:15`), but they were never the whole population.

  Dropping the flag is nonetheless safe in every configuration that exists, because `-lcblas` is
  never the provider of those symbols:

  | flavor | `libcblas.so` | SLATE undefined `cblas_` | who provides them |
  |---|---|---|---|
  | `gcc` | symlink → `libopenblas.so`, i.e. the same file `libblas.so` and `liblapack.so` point at | 0 | nobody needs to; `-lcblas` is a redundant third alias |
  | `debug` | absent | 0 | nobody needs to; `libblas.so` there is reference netlib and exports **no** `cblas_` symbol |
  | `mkl` | absent | 4 | `libmkl_intel_lp64` (verified `cblas_dgemm_batch`, `cblas_zgemm_batch` defined) |

  OpenBLAS 0.3.33 also defines all four (verified against `/opt/scls/gcc/lib64/libopenblas.so`), so
  a system-BLAS fallback covers SLATE too. Confidence: high — all rows are `nm -D` measurements, not
  inference.

- ~~**D3 (from Codex round 1, accepted) — the probe must cover SLATE's CBLAS symbols.**~~ **DISSOLVED 2026-09-04 with the probe.** BELFEM never calls `cblas_*gemm_batch`; only the `mkl` flavor's SLATE build has the four undefined, and MKL itself provides them. Applied to R4 it would have refused the `debug` flavor (Grok, scope jury). The residual
  risk Codex identified is real and survives the table above: a *reference-netlib-only* fallback
  (`libblas`/`liblapack` with no OpenBLAS, which is exactly what the `debug` flavor ships) combined
  with a SLATE built to need `cblas_*gemm_batch` would probe clean on `dgemm_`/`dgetrf_` and then
  fail at the final link. Keeping `-lcblas` does **not** fix that — netlib CBLAS has no batched GEMM
  either. The fix is to widen the probe: **when `USE_STRUMPACK` is on, the required-symbol set gains
  `cblas_dgemm_batch`.** That converts the one dangerous combination into a configure-time error,
  which is the entire point of this work.
- **The probe's link order must mirror `finalize_compiler.cmake`.** Ranks are emitted 10→0
  (`config/compiler/finalize_compiler.cmake:43-47`), so ScaLAPACK (rank 5) precedes the BLAS block
  (rank 1), which precedes the runtime libraries (rank 0). A probe that reorders these can produce a
  false negative on static archives. The probe reproduces this subsequence and nothing else.
- **D2 (self-found, measured 2026-08-30) — `try_compile` does NOT inherit `link_directories()`.**
  The probe cannot rely on the `link_directories($SCLS/lib64)` that `find_scls.cmake:7` sets for the
  real build; it must pass every search path explicitly as `LINK_OPTIONS "-L…"`. Measured in a
  standalone CMake project: an identical probe links `-lopenblas` from `/opt/scls/gcc/lib64` with an
  explicit `-L` and fails without it, `link_directories()` notwithstanding. Had this gone unnoticed,
  every fallback-BLAS probe would have false-failed and the configure would have rejected a working
  toolchain. Confidence: high — direct measurement. **R2 and R4 must therefore thread the search
  path list into the probe, not assume ambient linker state.**
- **Detection must not run before `detect_compiler.cmake`.** `CMAKE_CXX_COMPILER` is already
  `mpicxx` at the point `config_mkl.cmake` is included (`CMakeLists.txt:155,161,164`), which is what
  makes an MPI-aware probe possible. The new code stays at that inclusion point.

## 4. Ordered Steps

- [x] **R0 — read the SCLS flavor before the knobs are declared, and preset them.** New
      `config/system/find_scls_flavor.cmake`, `include()`d from `CMakeLists.txt` **above** the
      "User Settings" block — i.e. before `:71` and `:90`, not at `:158` where `find_scls.cmake`
      runs (§6.2b). It does four things and nothing else:
      1. Read `$SCLS`. If unset, set `BELFEM_SCLS_FLAVOR` to empty and return — every default keeps
         its plain-Linux value and no message is printed. This is the *majority* path for anyone
         outside this group and must stay silent and boring.
      2. Read `$SCLS_FLAVOR`. If unset but `$SCLS` is set, **derive it from the basename of
         `$SCLS`** — the flavors live at `/opt/scls/<flavor>`, so the basename is the flavor by
         construction, and this keeps a shell that exports only `SCLS` working.
      3. Accept it only if it is one of `debug`, `gcc`, `mkl`. Anything else — `cea`, or a `$SCLS`
         pointing at some private prefix — sets no preset and emits `message( STATUS )` saying the
         prefix is used as a search root but its flavor is unrecognized, so the knobs keep their
         defaults. **Not a `FATAL_ERROR`:** pointing `$SCLS` at a hand-rolled prefix is legitimate,
         and `find_scls.cmake:4-20` already fatals on a prefix with no `lib`.
      4. Compute `BELFEM_DEFAULT_USE_MKL` and `BELFEM_DEFAULT_USE_DEBUG` per the §6.1 table, and
         pass them as the third argument of the two `option()` calls.

      **Mismatch between `$SCLS_FLAVOR` and the basename of `$SCLS`** (e.g. `SCLS_FLAVOR=mkl` with
      `SCLS=/opt/scls/debug`) is a shell that got half-updated, and silently trusting either one
      would preset the wrong backend. Emit `message( WARNING )` naming both, then **trust
      `$SCLS_FLAVOR`**, which is the more specific declaration — and note the probe catches the
      consequence in any case. Keep the warning: this is precisely the state a half-sourced module
      leaves behind.

      Also print the flavor and the two presets it produced through R6, so the configure log answers
      "why is this tree Release?" without anyone having to reconstruct it.
      (after: nothing — R0 is independent of the probe work and can land first)
- [x] **R0b — the `CLAUDE.md` companion edit.** The moment R0 lands, "Debug builds (`USE_DEBUG=ON`,
      the default)" in `CLAUDE.md` is false. Restate it as: `USE_DEBUG` defaults ON under
      `SCLS_FLAVOR=debug` and OFF otherwise. `scripts/check_doc_claims.py` does *not* pin this
      default (it pins `USE_TEST`), so nothing mechanical will catch the staleness — which is an
      argument for extending the checker to `USE_DEBUG` while touching it. (after: R0)

      **As built:** the flag list in `CLAUDE.md` now states the per-flavor default table, that the
      presets are defaults only, and that a cached value is never overwritten. **The checker was
      not extended** — `check_doc_claims.py:135` matches `option( NAME "…" ON|OFF )` against a
      literal, so `USE_DEBUG`'s default is no longer visible to it at all (`USE_TEST` and `USE_VTK`
      still are, and the run is green at 37/37). Pinning a *conditional* default needs a different
      check from the one that file performs; left as a known gap rather than a silent one.
- [x] ~~**R1 — Tri-state `USE_MKL`.**~~ **WITHDRAWN 2026-08-30 (Christian): no autodetection.**
      `USE_MKL` stays the boolean `option( USE_MKL … OFF )` at `CMakeLists.txt:71`, and every
      `if( USE_MKL )` site stays as written — with no `AUTO` value there is no truthy-string hazard,
      so nothing needs converting to a resolved boolean. This step's whole content was machinery for
      a third state that no longer exists. What survives from it is one line of work, folded into R5:
      the stale comment at `config_mkl.cmake:90`. Withdrawing R1 also removes the round's most
      serious finding (D4), which existed only inside the three-way string chain.
- [x] ~~**R2 — `belfem_probe_link()` helper**~~ **DECLINED 2026-09-04 (scope jury): no link probe.** in `config/scripts/`, next to `belfem_find_package.cmake`,
      and `include()`d from `CMakeLists.txt` *before* `find_mkl.cmake` / `find_blas.cmake` use it.
      Given a list of `-L` paths, a list of link items, and a list of symbol names, it writes a
      translation unit declaring and calling each symbol, runs `try_compile`, and returns success
      plus the linker output. Constraints, all of them consequences of the probe running before
      `finalize_compiler.cmake` (Codex round 1, confirmed):
      - **Use the CMake 3.11 `try_compile` signature** — `CMAKE_FLAGS -DCMAKE_EXE_LINKER_FLAGS=…`
        and/or full library paths, **not** `LINK_OPTIONS`, which is 3.25+ (D5).
      - **Forward `CMAKE_C_COMPILER`/`CMAKE_CXX_COMPILER`** so the nested project uses `mpicc`/`mpicxx`
        rather than plain `g++` (D10).
      - **Nothing ambient may be assumed.** Pass `-L` explicitly for both `$SCLS/lib64` and
        `${BELFEM_MKL_LIBDIR}` (D2 + D9), and pass the OpenMP flags, `BELFEM_FORTRANLIBS`, and the
        rank-0 runtime (`-lpthread -lm -ldl`) explicitly too — `finalize_compiler.cmake:7` only applies
        OpenMP to `CMAKE_CXX_FLAGS` *after* this point, `CMAKE_C_FLAGS` is never populated from
        `BELFEM_CFLAGS`, `BELFEM_FORTRANLIBS` is appended per target in
        `Add_BelfemLibrary.cmake:35-38`, and `-L${SCLSLIBDIR}` reaches
        `CMAKE_EXE_LINKER_FLAGS` only at `CMakeLists.txt:321-323`.
      - **Probe in C++, not C**, so the driver matches the production link.
      - **Results must not stick across a `$SCLS` repoint.** `find_library()` caches, and
        `try_compile` caches its result variable; both must be keyed on the candidate paths/options
        or explicitly unset, or repointing `$SCLS` will silently reuse the previous flavor's verdict
        (Codex round 1). This is what makes O1's "re-resolve on reconfigure" actually true.
      (after: nothing)
- [x] ~~**R3 — MKL verification**~~ **DECLINED 2026-09-04 with R2.** in `find_mkl.cmake`, reached only when `USE_MKL` is true: resolve
      `BELFEM_MKLROOT`, pick `BELFEM_MKL_LIBDIR` by *library content* rather than directory
      existence, then probe the assembled `BELFEM_MKL_LIBS` for the required-symbol set: `dgemm_`,
      `dgetrf_`; under MPI additionally `blacs_gridinit_` and `pdgetrf_`; under `USE_STRUMPACK`
      additionally `cblas_dgemm_batch` (D3). A failed probe is fatal — the user asked for MKL and
      MKL is not linkable — but the `FATAL_ERROR` and its text belong to R5, so this step reports
      through `BELFEM_MKL_USABLE` and the captured linker output rather than aborting itself.
      **`find_mkl.cmake:12-13` stays** (a missing `$MKLROOT` under `USE_MKL=ON` really is an error,
      which is what dissolved D6); improve its message to name `$MKLROOT`, the path tried, and
      `-DUSE_MKL=OFF` as the way out. The existing Apple refusal at `config_mkl.cmake:2-4` also
      stays as written, since `USE_MKL` is only ever truthy when the user asked for it (D7
      dissolved). (after: R2)
- [x] ~~**R4 — probe the BLAS/LAPACK link line that the netlib branch already emits.**~~ **DECLINED 2026-09-04 with R2.** **Decided
      2026-08-30 (Christian, O6): no candidate ladder and no vendor knob.** Under `USE_MKL=OFF` the
      configure keeps looking for BLAS/LAPACK exactly as it does today — `-llapack -lblas`, in that
      order, resolved by the linker's own search of `-L$SCLS/lib64` and then the default paths — and
      the only thing this step adds is that the result is now *checked* instead of assumed.
      `BELFEM_BLAS_LIBS` therefore keeps the value it has after D1; no `find_blas.cmake` is created
      and `BLAS_DIR` is not introduced.

      **Why no ladder is needed — measured 2026-08-30, and this is the fact that makes O6 cheap:**
      every flavor that supports `USE_MKL=OFF` exposes the netlib *names*, whatever implements them.

      | flavor | `libblas.so` / `liblapack.so` | what the pair actually links |
      |---|---|---|
      | `/opt/scls/gcc` | both present, **symlinks to `libopenblas.so`** | OpenBLAS, reached by its netlib alias |
      | `/opt/scls/debug` | both present, real files | reference netlib |
      | `/opt/scls/mkl` | **neither** (ships `libscalapack.so.2.2.3` and the `blaspp`/`lapackpp` wrappers only) | nothing — `USE_MKL=OFF` there must error, and now will |

      A four-rung ladder would have spent its complexity discovering that `openblas` and
      `lapack;blas` name the same file on `gcc`. Probing the pair covers both flavors directly.

      Two things do carry over. **Under MPI, probe ScaLAPACK together with the BLAS pair, not
      separately** — the `mkl` prefix ships `libscalapack.so.2.2.3` while shipping no BLAS, so an
      independent resolution could bind mkl-prefix ScaLAPACK to a BLAS from somewhere else (D11,
      second half; the first half dies with the ladder). And the **error text must name the
      packages** a user is expected to install (`blas-devel`/`lapack-devel`, or their OpenBLAS
      equivalents) rather than saying "BLAS not found" — on this host `/usr/lib64` carries no
      `libblas.so`/`liblapack.so` devel symlinks at all, so an unhelpful message here is the likely
      case, not the rare one. Exhausting the search is a `FATAL_ERROR`, never a promotion to MKL.
      (after: R2)
- [x] ~~**R5 — verdict and error text**~~ **DECLINED 2026-09-04: no verdict and no error text — the knobs are defaults the user may override freely; the `-DBLAS_DIR=…` clause below contradicted R4 and is struck with the rest.** in `config_mkl.cmake`. The table is the whole contract, and
      neither row has a fallback edge — that is the point of the 2026-08-30 decision:
      | `USE_MKL` | outcome |
      |---|---|
      | true (`ON`, `1`, `YES`, …) | MKL, or `FATAL_ERROR`. **Never falls back to a system BLAS**, however good the system BLAS is |
      | false (`OFF`, `off`, `0`, `NO`, unset) | a BLAS/LAPACK distribution, or `FATAL_ERROR`. **Never promotes to MKL**, even when MKL is sitting right there |

      CMake's own `option()`/`if()` truthiness handles every spelling in both rows, which is why the
      withdrawal of R1 costs nothing here. Error text in the house style of the MPICH refusal
      (`find_mpi.cmake:74-84`): what was asked for, what was tried, what was found, and the ways out
      (`-DUSE_MKL=ON`/`OFF`, install `openblas-devel` / `flexiblas-devel`, `-DBLAS_DIR=…`). Include
      the probe's captured linker output — that is the line that actually tells the user what is
      missing, and it is the whole reason this task exists. Also fix the stale comment at
      `config_mkl.cmake:90` (the surviving fragment of R1). The PARDISO refusal at `:83-85` keys off
      `USE_MKL` and needs no change (O2/D8 dissolved). (after: R3, R4)
- [x] ~~**R6 — `summary.cmake`**~~ **DECLINED 2026-09-04: the group does not need to know which library won the link.**: print the resolved backend, the libdir, and the library list, so the
      configure log records which backend a tree actually got. (after: R5)
- [x] ~~**R7 — gate.**~~ **REDUCED 2026-09-04 to one fresh configure per flavor (`mkl`, `gcc`, `debug`) checking the presets; Christian's to run.** Configure across the flavor matrix and check each outcome, then build to the
      link stage. `$SCLS` selects the flavor:
      | `$SCLS` | `USE_MKL` | expected |
      |---|---|---|
      | `/opt/scls/mkl` | unset | **presets ON** from `SCLS_FLAVOR=mkl`; resolves to MKL. `USE_DEBUG` presets OFF |
      | `/opt/scls/mkl` | `ON` | resolves to MKL |
      | `/opt/scls/mkl` | `OFF` | resolves to a *system* BLAS if one is installed, else `FATAL_ERROR`. Either way it must **not** silently use MKL. (The draft asserted an unconditional error; that conflicts with R4's system-path search — corrected after Codex round 1. On this host no system devel symlinks exist, so the observed result should be the error) |
      | `/opt/scls/gcc` | `OFF` | resolves to the SCLS OpenBLAS through the netlib names |
      | `/opt/scls/gcc` | unset | **presets OFF**; same resolution. `USE_DEBUG` presets OFF |
      | `/opt/scls/debug` | `OFF` | resolves to reference netlib — this is the case that fails today |
      | `/opt/scls/debug` | unset | **presets OFF**, and **`USE_DEBUG` presets ON** — the one flavor that turns the debug build on |
      | `/opt/scls/debug` | `ON` | `FATAL_ERROR` naming MKL, *not* a silent netlib fallback and *not* an MKL-plus-netlib mixture. `$MKLROOT` is set and valid on this host, so the probe — not the presence of a directory — is what decides |
      | `-DUSE_MKL=off` (lowercase), any flavor | behaves as `OFF` via CMake's own truthiness; no longer a regression test, since D4's chain is gone |
      | **`$SCLS` unset entirely** | unset | `USE_MKL=OFF`, `USE_DEBUG=OFF`, no preset messages; probes `-llapack -lblas` from the default paths. On this host that **errors**, since `/usr/lib64` has no `libblas.so`/`liblapack.so` devel symlinks — and the error text is then the whole user experience, so read it as a stranger would |
      | `$SCLS` set, `SCLS_FLAVOR` unset | unset | flavor derived from the basename; same outcome as the matching row above |
      | `SCLS_FLAVOR=cea`, or `$SCLS` at a private prefix | unset | no preset, `STATUS` message, knobs keep plain defaults; must not abort in `find_scls_flavor.cmake` |
      | `SCLS_FLAVOR=mkl` with `SCLS=/opt/scls/debug` | unset | `WARNING` naming both, presets from `SCLS_FLAVOR`, then the MKL probe decides |
      | existing tree, `$SCLS` repointed `debug` → `mkl` | (cached `OFF`) | preset does **not** fire; the BLAS probe fails and errors. Verifies §6.4's third row |

      The two `AUTO` rows of the round-1 matrix are struck: with autodetection withdrawn there is no
      unset-means-maybe case left to gate. The `debug`+`ON` row replaces the "dangerous AUTO row" and
      tests the same hazard from the deliberate side — a user who asks for MKL inside a netlib-built
      flavor must be told, not quietly accommodated.
      Then `make` far enough to link one executable in a `gcc`-flavor `USE_MKL=OFF` tree and in an
      MKL tree. Until those links run, this work is *reviewed*, not verified (protocol §11).
      Builds are Christian's to run. (after: R6)

## 4.1 Round-1 Defect Tracker (2026-08-30 — Codex `terra`/`high`, Grok `grok-4.6`/`high`, parallel, identical brief)

Every finding below was re-verified against the tree or by measurement before acceptance; auditor
agreement alone was not treated as evidence (protocol §11).

- [x] **D1 — `-lcblas` is dead weight. CLOSED, and the closing evidence is not the one first given.**
  Codex refuted "no cblas consumer exists" (SLATE, via STRUMPACK). Confirmed. The conclusion survives
  on measured grounds (§3.1 table). Grok asked for the one gate it could not run — `readelf -d` on the
  debug `blaspp`, whose `blasppConfig.cmake:23-28` records `-lcblas` although the prefix ships no
  `libcblas`. **Run 2026-08-30: `libblaspp.so.2.0.0` has `NEEDED libblas.so.3 liblapack.so.3` and zero
  `NEEDED` cblas; no debug-flavor DSO carries any undefined `cblas_` symbol.** The recorded string is
  stale metadata, not a runtime dependency. D1 closed on link *and* runtime evidence.
- [x] ~~**D4 — CRITICAL, from Grok. A three-way `STREQUAL` silently violates the override rule.**~~
  **DISSOLVED 2026-08-30 by the no-autodetection decision** — the hazard lived entirely inside the
  three-way chain R1 would have introduced, and there is no longer a third state to chain. Kept on
  the record because it is the reason the tri-state shape was expensive: the finding was correct,
  and the cheapest fix for it turned out to be not needing it. Grok's original text:
  `option()` accepted `off`, `FALSE`, `0`, `NO` case-insensitively; `STREQUAL "OFF"` matches none of
  them, so `-DUSE_MKL=off` falls through the `AUTO`/`ON`/`OFF` chain to the AUTO branch — and on a host
  where MKL is present, AUTO resolves to MKL. **An explicit rejection would have been silently
  overridden**, which is precisely the guarantee this work exists to provide. Measured across nine
  spellings: the naive chain mishandles seven of them; Grok's shape handles all nine.
  ```cmake
  if( USE_MKL STREQUAL "AUTO" )      # string first — AUTO is truthy, so this must precede
      # probe MKL, then the ladder
  elseif( USE_MKL )                  # ON, on, TRUE, 1, YES, and BOOL ON from a legacy cache
      # require MKL
  else()                             # OFF, off, FALSE, 0, NO, and BOOL OFF
      # ladder only, never MKL
  endif()
  ```
- [x] ~~**D5 — from Grok. `try_compile( … LINK_OPTIONS … )` is a CMake 3.25 signature; this tree's floor
  is `cmake_minimum_required(VERSION 3.11)` (`CMakeLists.txt:1`, verified).**~~ **DISSOLVED 2026-09-04 — died with the probe helper.** D2's stated recipe would
  error at configure on any CMake between 3.11 and 3.24. Use the 3.11-safe form —
  `CMAKE_FLAGS -DCMAKE_EXE_LINKER_FLAGS=…` and/or full library paths in `LINK_LIBRARIES`. (The cmake
  on this host is 4.4.2, which is exactly why this would have escaped local testing.)
- [x] ~~**D6 — from Grok. `find_mkl.cmake` is the *detector*, so it cannot be gated on the resolved
  boolean.**~~ **DISSOLVED 2026-08-30.** The circularity was an AUTO artefact: with no AUTO,
  `find_mkl.cmake:2` reads `USE_MKL` — an unambiguous statement of user intent available before any
  probe runs — and the file is only entered when MKL was actually requested. Its `FATAL_ERROR` at
  `:12-13` is therefore **kept, not deleted**: under `USE_MKL=ON` a missing `$MKLROOT` is exactly an
  error. Only its wording is improved (R3). The split D6 asked for still holds in spirit — R3
  detects, R5 phrases the refusal.
- [x] ~~**D7 — from Grok. Apple AUTO would abort the configure.**~~ **DISSOLVED 2026-08-30.** The
  abort required `AUTO` to be truthy on a fresh macOS clone. With the default back to `OFF`,
  `config_mkl.cmake:2-4` fires only when a macOS user explicitly asks for MKL, which is a correct
  refusal and today's behaviour. No change needed. (Related, out of scope: the latent `BELFEM_ACCELLERATE` path at
  `blaze_config.hpp:63-65` *does* want cblas via `vecLib`; no CMake file sets that define today, but
  D1's removal must not be cited as licence to drop cblas on an Accelerate path.)
- [x] ~~**D8 — from Grok. O2 and R1 fight.**~~ **DISSOLVED 2026-08-30.** The conflict needed a
  resolved boolean that could differ from what the user typed. Without AUTO the two are the same
  value, so `USE_PARDISO=ON` under a false `USE_MKL` still hits the existing refusal at
  `config_mkl.cmake:83-85` and no re-keying is required. O2 collapses with it.
- [x] ~~**D9 — from both. The probe's `-L` set must be explicit and must include the MKL libdir.**~~ **DISSOLVED 2026-09-04 — died with the probe helper.**
  Codex and Grok independently noted that `config_mkl.cmake:19-20` appends the MKL libdir to
  `BELFEM_RPATH` only; `finalize_compiler.cmake:74-77` copies that into `LINK_DIRECTORIES` *later*.
  So even the real build has no `-L$MKLROOT/…` at probe time. Pass `-L` for **both** `$SCLS/lib64`
  and `${BELFEM_MKL_LIBDIR}`.
- [x] ~~**D10 — from Grok. The nested project must inherit `mpicxx`.**~~ **DISSOLVED 2026-09-04 — died with the probe helper.** `detect_gcc.cmake:148-150` sets
  the MPI wrappers, so they are available — but a `try_compile` that does not forward them picks
  plain `g++`, and every BLACS/ScaLAPACK probe false-fails.
- [x] ~~**D11 — from Grok. Roots must not be mixed.**~~ **DISSOLVED 2026-09-04 — died with the probe helper.** ~~R4's third rung is load-bearing~~ — the first
  half is moot after O6: with no ladder, `lapack;blas` is not a rung that could lose, it is the only
  thing probed. **The second half stands and is now the whole of D11:** the `mkl` prefix ships
  `libscalapack.so.2.2.3` while shipping no BLAS, so resolving ScaLAPACK and BLAS independently can
  bind mkl-prefix ScaLAPACK to a BLAS from elsewhere. **Probe ScaLAPACK together with the BLAS pair,
  and require them from the same root.**

Converged independently in both audits (strongest signal in the round): `cea` is not a flavor, the
ILP64 overclaim, the R1 call-site miscount, and the R7 `mkl`+`OFF` expectation.

## 5. Open Design Questions

- ~~**O1 — Does `AUTO` belong in the cache as `AUTO`, or should the first configure freeze the
  resolved value?**~~ **DISSOLVED 2026-08-30 — there is no `AUTO`.** A tree's backend is whatever
  the user's `USE_MKL` says, on every reconfigure, and it cannot change without the user changing
  it. The one part worth keeping is the mitigation: **R6 prints the resolved backend on every
  configure anyway**, so a repointed `$SCLS` that silently changes *which* BLAS satisfies the same
  `USE_MKL=OFF` is still visible in the log.
- ~~**O2 — Should `USE_PARDISO=ON` under `AUTO` force MKL rather than error?**~~ **DISSOLVED
  2026-08-30 — the question only existed under AUTO.** `config_mkl.cmake:83-85` keeps erroring with
  "Turn on MKL if you want to use Pardiso", which is now simply correct: the user's `USE_MKL` is the
  only input, so PARDISO+non-MKL is an explicit contradiction and saying so is the right answer.
  D8's re-keying requirement dies with it.
- ~~**O3 — ladder order when a flavor offers several names.**~~ **DISSOLVED 2026-08-30 by O6 —
  there is no ladder to order.** The observation behind it survives and is worth keeping, because it
  is why no ladder is needed: `/opt/scls/gcc/lib64` ships `libopenblas.so` *and*
  `libblas.so`/`liblapack.so` as symlinks to it, so `openblas` and the netlib pair are the same
  library reached two ways. The one real cost of linking through the alias is cosmetic — R6's
  summary line will say "BLAS/LAPACK" where the implementation is OpenBLAS. If that matters later,
  resolve the symlink for the *report* only, never for the link line.
- ~~**O5 — flavor/backend coherence: what should `AUTO` do on a flavor whose own TPLs were built
  against a *different* BLAS?**~~ **DISSOLVED as an `AUTO` question 2026-08-30; its measurements are
  retained because they are the evidence O4 needs.** With no autodetection, nothing can promote the
  `debug` flavor to MKL behind the user's back — a mixed process now requires someone to type
  `-DUSE_MKL=ON` inside a netlib flavor, which R7's `debug`+`ON` row turns into a refusal.
  The measurement stands and is worth keeping, because it is what tells that refusal from a false
  alarm. **On this host (2026-08-30)** `$MKLROOT` is set to `/opt/intel/oneapi/mkl/latest` and that
  directory exists, so a directory-existence test would call MKL usable there — while every BLAS
  consumer inside `/opt/scls/debug` binds netlib:

  | consumer | what it records |
  |---|---|
  | PETSc | `BLASLAPACK_LIB = -Wl,-rpath,/opt/scls/debug/lib -L/opt/scls/debug/lib -llapack -lblas` |
  | `libstrumpack.so`, `libsuperlu.so` | `NEEDED liblapack.so.3`, `NEEDED libblas.so.3` |
  | Armadillo | `ARMA_AUX_LIBS /opt/scls/debug/lib/libblas.so;/opt/scls/debug/lib/liblapack.so` |
  | SLATE | no direct BLAS `NEEDED`; reaches it through `libblaspp.so.2` / `liblapackpp.so.2`, which carry the same two sonames |

  BELFEM's own calls would go to MKL and its TPLs' to netlib: two BLAS implementations in one
  process, each with its own threading layer. **What survives into the new design:** this is why
  `debug`+`ON` must be a `FATAL_ERROR` rather than a shrug, and why R3 probes the MKL *libraries*
  instead of testing whether `$MKLROOT` exists. The half of the question that is still open — may a
  flavor's backend come from outside that flavor at all — is O4, which no longer has a mirror
  image and can be decided on its own.
- [x] **O6 — under `USE_MKL=OFF`, how is the BLAS/LAPACK distribution named? DECIDED 2026-08-30
  (Christian): option (a), and narrower than drafted — no `BLAS_VENDOR` knob, and CMake keeps
  looking for BLAS/LAPACK exactly as it does today.** The link line stays `-llapack -lblas`; only
  the check around it is new. This also dissolves **O3** (there is no ladder, so no ladder order)
  and the first half of **D11**. R4 carries the measurement that makes it safe: on `gcc` those two
  names are symlinks to `libopenblas.so`, on `debug` they are real netlib, and `mkl` has neither —
  so the pair already reaches every backend a flavor offers. The options as they stood:
  - **(a) implicit ladder** — R4 as drafted: try `openblas`, `flexiblas`, `lapack;blas`,
    `openblas;lapack` in order and take the first that links. One knob total. The user says "not
    MKL" and the configure finds whatever BLAS is there, erroring only if none is. Ordering
    questions (O3) stay live.
  - **(b) explicit vendor knob** — e.g. `BLAS_VENDOR` (`openblas` | `flexiblas` | `netlib`,
    defaulting to the first that links, or to a required value). The user names the distribution and
    a miss is an error about *that* distribution, not about BLAS in general. This is the direction
    the "MKL **or** a BLAS/LAPACK distro is selected" phrasing of the decision most literally
    implies, and it makes O3 moot: order stops mattering when the vendor is named.
  Both satisfied the no-substitution rule. (b) was rejected as a knob that buys sharper error text
  at the cost of a second thing to set, document and keep in sync with reality — and the reality is
  that the netlib names already resolve to whatever each flavor ships. R4 recovers the error-text
  sharpness by naming the expected packages directly.

## 6. Interface Design: the Flavor Supplies the Defaults, the Knobs Still Rule

**Decided 2026-08-30 (Christian).** The build does not guess which linear-algebra backend the user
wants. It reads `$SCLS_FLAVOR` if there is one, uses it to *default* the knobs, and then checks that
whatever the knobs say can actually be linked.

### 6.1 The preset table

`/opt/scls` holds three build flavors plus one impostor, verified 2026-08-30 by listing the tree:

| `$SCLS_FLAVOR` | `$SCLS` | BLAS/LAPACK it provides | `USE_MKL` default | `USE_DEBUG` default |
|---|---|---|---|---|
| `debug` | `/opt/scls/debug` | reference netlib — real `libblas.so` / `liblapack.so` | `OFF` | **`ON`** |
| `gcc` | `/opt/scls/gcc` | OpenBLAS — `libblas.so` / `liblapack.so` are symlinks to `libopenblas.so` | `OFF` | `OFF` |
| `mkl` | `/opt/scls/mkl` | none of its own; the flavor's TPLs are built against oneAPI MKL at `$MKLROOT` | **`ON`** | `OFF` |
| *(none — `$SCLS` unset)* | — | whatever `blas-devel` / `lapack-devel` put in the default paths | `OFF` | `OFF` |
| ~~`cea`~~ | `/opt/scls/cea` | **not a flavor** — thermodynamic input decks, no `lib`/`include` | no preset | no preset |

Two facts from that table do real work later. **`debug` and `gcc` are identical to link against** —
both expose the netlib *names*, and only the implementation behind them differs — which is why one
`-llapack -lblas` covers both and no candidate ladder is needed (O6, R4). And **`mkl` is the flavor
whose backend lives outside `$SCLS`**: `/opt/scls/mkl` ships no BLAS at all, so "search `$SCLS`
first" applies to the TPLs and never to MKL itself, which is still resolved through `$MKLROOT`
(R3).

The three flavor directories are structurally identical (`bin doc etc include lib lib64 libexec
sbin share`) with no marker file inside them, so **`$SCLS_FLAVOR` — or, failing that, the basename
of `$SCLS` — is the only signal available.** There is nothing in the tree to sniff.

### 6.2 Two consequences that need deciding with open eyes

**(a) The `USE_DEBUG` default flips for everyone without SCLS.** `CMakeLists.txt:90` is
`option( USE_DEBUG "Compile with debug flags" ON )` today, so a fresh clone builds Debug. Under this
design a fresh clone with no `$SCLS` builds **Release**. That is a defensible default for a stranger
cloning the repository, and it is what Christian asked for — but it is a change in what the project
does out of the box, and it drags one companion edit with it: **`CLAUDE.md` states "Debug builds
(`USE_DEBUG=ON`, the default)" and becomes wrong the moment R0 lands.** The same session must fix
it. (`scripts/check_doc_claims.py` pins the `USE_TEST` default but not this one, so the checker will
*not* catch the staleness — this note is the only guard.)

**(b) The presets must be computed 68 lines earlier than SCLS is read today.** `option( USE_MKL … )`
is `CMakeLists.txt:71` and `option( USE_DEBUG … )` is `:90`, but `find_scls.cmake` is not included
until `:158`. An `option()` whose default is computed after it has already been declared does
nothing at all — and does it silently, since the first configure just caches the un-preset value.
R0 therefore reads `$SCLS`/`$SCLS_FLAVOR` into plain variables *before* the "User Settings" block.
The rest of `find_scls.cmake` (the `link_directories()`, `BELFEM_RPATH` and include-path work) stays
where it is; only the environment read moves up.

### 6.3 What the knobs mean once they are set

Nothing in the `option()` declarations changes except their default expression. All of this task's
remaining value is behind them:

| `USE_MKL` | today | after this task |
|---|---|---|
| true | emits the MKL link block unchecked | probes it; links, or `FATAL_ERROR` naming MKL and what the linker said |
| false | emits `-llapack -lblas` unchecked | resolves a BLAS/LAPACK distribution and probes it; links, or `FATAL_ERROR` naming what was tried |

Why this is the better trade, recorded so the round-1 design is not re-proposed: autodetection buys
convenience on hosts that have several backends, and pays for it with a build whose backend can
change without anyone asking — the failure it invites is *silent and wrong* (two BLAS runtimes in
one process, §5 O5), while the failure this design invites is *loud and immediate* (a configure
error on a host missing a backend). For a framework whose TPL stack is built per flavor against one
specific BLAS, the loud failure is worth more. The preset in §6.1 recovers the convenience without
the cost, because it reads what the user already declared instead of inferring it.

### 6.4 How preset, cache and probe interact

The three mechanisms have to be understood together or R0 will look broken:

| situation | what happens | is that right? |
|---|---|---|
| fresh tree, `$SCLS_FLAVOR=mkl` | `USE_MKL` defaults ON, R3 probes MKL | yes — the point of the feature |
| fresh tree, `$SCLS_FLAVOR=mkl`, `-DUSE_MKL=OFF` | the `-D` wins; R4 probes `-llapack -lblas`, finds no BLAS in `/opt/scls/mkl`, and errors | yes — explicit beats preset, and the refusal is loud |
| **existing** tree, `$SCLS` repointed `debug` → `mkl`, reconfigure | the cached `USE_MKL:BOOL=OFF` survives; **the preset does not fire**, because `option()` cannot overwrite a populated cache | acceptable, *because the probe then fails loudly* — `/opt/scls/mkl` has no `libblas`. The failure mode is a clear error, not a silent mismatch |
| no `$SCLS` at all | `USE_MKL=OFF`, `USE_DEBUG=OFF`, default paths | yes — the plain-Linux path |

The third row is the one worth stating: a build tree is effectively per-flavor, and repointing
`$SCLS` under an existing tree is not supported so much as *caught*. This is the clearest case of
the two halves of the plan covering each other — the preset handles the fresh tree, the probe
handles the stale one. `summary.cmake` (R6) printing the resolved backend and libdir is what makes
the situation legible in the log either way.

### 6.5 Withdrawn: the tri-state `AUTO` knob

~~The round-1 design replaced the `option()` with a cache *string* defaulting to `AUTO`~~ —
withdrawn by the decision above, together with R1, D4, D6, D7, D8, O1 and O5. The measurements
below were run before the withdrawal and are kept for one reason: they document that a cache-type
change from `BOOL` to `STRING` is safe on existing trees, which is the fact a future session would
have to re-establish if this were ever revisited.

Semantics of the withdrawn shape, and why it had that form:

| user passes | cache holds | behaviour |
|---|---|---|
| nothing, fresh tree | `AUTO` (STRING) | probe MKL, then the BLAS ladder |
| `-DUSE_MKL=ON` | `ON` (STRING) | MKL or `FATAL_ERROR` |
| `-DUSE_MKL=OFF` | `OFF` (STRING) | BLAS ladder or `FATAL_ERROR` |
| nothing, **existing** tree configured before this change | `OFF` (BOOL, from the old `option()`) | `set(... CACHE ...)` does not overwrite an existing entry, so the old explicit value survives and is read as OFF |

The last row is the reason for a tri-state string rather than "detect whether the user set the
boolean": on a reconfigure, a cached `USE_MKL` is `DEFINED` whether the user set it or the previous
configure did, so an explicitness test cannot distinguish them. A distinct third value can.

All four rows were measured in a standalone CMake project (scratchpad probe, 2026-08-30), not
assumed:

| probe | result |
|---|---|
| fresh tree, no `-D` | value `AUTO`, cache type `STRING`, `STREQUAL "AUTO"` yes |
| fresh tree, `-DUSE_MKL=OFF` | value `OFF`, type `STRING`, `STREQUAL "OFF"` yes, `if()` reads FALSE |
| tree pre-seeded `USE_MKL:BOOL=OFF` | entry survives as **BOOL**, value `OFF`, `STREQUAL "OFF"` still yes — the compatibility promise holds across the type change |
| `if( USE_MKL )` with value `AUTO` | **TRUE** — confirms the hazard below |

~~`if( USE_MKL )` must not survive anywhere after R1~~ — moot: the truthy-string hazard needed the
string. **Every `if( USE_MKL )` in the tree stays exactly as it is**, and that is now a required
property rather than a risk. The grep that enumerated the sites,
`grep -rn USE_MKL --include=CMakeLists.txt --include=*.cmake .`, is still the check to run — but it
now confirms that *nothing* changed shape.

- **O4 (raised by Codex, round 1) — may the `mkl` SCLS flavor use a non-MKL backend at all?**
  Its TPL stack is built against MKL: SLATE there needs `cblas_*gemm_batch`, and the rest of the
  flavor's libraries were linked expecting MKL's threading and interface layers. A system-OpenBLAS
  fallback under `SCLS=/opt/scls/mkl` would satisfy the probe (OpenBLAS defines those four symbols)
  yet mixes two BLAS implementations in one process — the toolchain mixing `find_mpi.cmake:9-11`
  warns about, in a worse form. Options: (a) allow it, probe-gated, as R4 currently does;
  (b) warn; (c) refuse when `$SCLS` is a flavor whose own TPLs resolve against MKL. Proposed: **(b)**,
  a `message( WARNING )` naming the flavor — refusing would need a reliable way to tell what a flavor
  was built against, which we do not have. **Grok sharpened this into the decisive question of the
  round:** R4's default-path search and R7's `mkl`+`OFF` expectation are in direct conflict. Today
  `mkl`+`OFF` errors only because this host has no `find_library`-visible system BLAS; install
  `openblas-devel` — which the error text itself recommends — and `mkl`+`OFF` will happily configure
  against a system BLAS while MPI, PETSc and STRUMPACK all come from `/opt/scls/mkl`. The same
  conflict produces the `debug`+AUTO row in R7, analysed as **O5** (§5); decide the two together.
  **Pick one:** (i) default system paths are a
  documented escape hatch and the flavor-mixing rows are host-specific, or (ii) when `$SCLS` is set,
  the ladder never falls back to default system paths.

  **Largely resolved as a consequence of O6 (2026-08-30) — flagged rather than silently closed,
  because Christian has not ruled on it directly.** With the ladder gone there is no search of our
  own left to restrict: `-llapack -lblas` is resolved by the linker, searching `-L$SCLS/lib64` and
  then its default paths, which is option (i) by construction. Option (ii) would now mean *the probe
  refusing a link the real linker would accept* — and §2 sets out why that is the one thing the
  probe must never do: a gate that asks a different question from the final link is a gate that
  lies, in whichever direction it errs. So (i) stands, with two consequences worth stating rather
  than discovering later. First, `mkl`+`OFF` on a host that *does* have `blas-devel` installed will
  configure and link against a system BLAS while MPI, PETSc and STRUMPACK come from
  `/opt/scls/mkl`; that is a pre-existing property of the link line and not something this task
  introduces, but it is also not something this task fixes. Second, **R6's summary line is the
  mitigation** — printing the libdir that actually won is what makes such a mixture visible instead
  of silent, which promotes R6 from a nicety to the counterweight for this decision. If the stricter
  behaviour is wanted, the place for it is a separate check ("`$SCLS` is set but the BLAS resolved
  outside it") that *warns*, leaving the probe honest.

## 7. Definition-of-Done Checklist

Items still open on 2026-09-04 belonged to the probe and were struck with it; the owed gate is the
per-flavor configure in the banner.

- [x] ~~Every gap-table row mapped to a step or an open question.~~
- [x] ~~Each claimed gap backed by a citation, not assumption.~~
- [x] ~~Ordered steps with dependencies.~~
- [x] Open questions logged, not decided *by the plan* — O1–O6 were all put to Christian and
      answered on 2026-08-30 rather than assumed away.
- [x] D1 (`-lcblas`) adjudicated by both auditors and closed on link **and** runtime evidence (§4.1).
- [x] ~~D4 regression test in R7~~ — dissolved with the tri-state knob (§4.1).
- [x] ~~D5: no `LINK_OPTIONS` anywhere in the helper (CMake floor is 3.11).~~
- [x] ~~**No backend substitution anywhere in the result.** `grep` the finished CMake for any path~~
      that reaches the MKL block under a false `USE_MKL`, or the BLAS ladder under a true one.
      This is the decision of 2026-08-30 restated as a check, and it is the one that matters most.
- [x] ~~**Every failure path carries the probe's linker output**, not just "not found".~~
- [x] ~~O6 decided before R4 lands~~ — decided 2026-08-30: no ladder, no vendor knob (§5).
- [x] ~~O4 decided before R4 lands~~ — resolved by consequence of O6; the probe mirrors the real
      link line rather than restricting it, and R6's summary is the mitigation (§5).
- [x] ~~**R6 is not optional.** It is the agreed counterweight to O4: a configure that resolves its~~
      BLAS outside `$SCLS` must say so in the log, and it is also where the flavor and its presets
      become visible (R0).
- [x] ~~**A no-`$SCLS` configure is silent about SCLS.** No warnings, no mention of flavors — the~~
      plain-Linux path must not look degraded to someone who has never heard of SCLS.
- [x] **`CLAUDE.md` updated in the same session as R0** (R0b). The `USE_DEBUG` default is stated in
      prose there and nothing mechanical guards it — still true after the edit, and now recorded as
      a known gap in R0b rather than an assumed guard.
- [x] **The presets are defaults, never assignments.** `find_scls_flavor.cmake` writes only
      `BELFEM_DEFAULT_USE_*` and `BELFEM_SCLS_FLAVOR`; it never touches `USE_MKL`/`USE_DEBUG`. Case 9
      (`-DUSE_MKL=OFF` beats an `mkl` preset) and case 10 (a cached value survives a repoint) confirm
      it from the outside.
- [x] ~~`BELFEM_BLAS_LIBS` still reads `-llapack -lblas` when the work is done — if the finished~~
      branch names a library the netlib branch did not name before, O6 was overstepped.
- [x] ~~O5 decided before R5 lands~~, ~~O2/D8 decided before R5 lands~~ — both dissolved (§5).
- [x] ~~R7 configure matrix passes and one executable links.~~

## 8. Audit Trail

- Exchange thread: `tmp/ai_exchange/linalg_backend_autodetect.md` — keeps the original slug, which
  is also this plan's former filename. The plan became `todo/linalg_backend_verification.md` on
  2026-08-30 when autodetection left the scope; the round-1 audits in that thread were run against
  the autodetection design and should be read with §4.1's dissolutions in hand.
- Round 1 (plan audit, 2026-08-30): Codex `gpt-5.6-terra`/`high`, Grok `grok-4.6`/`high`, parallel,
  identical brief.
- **Codex round 1 — returned, reconciled.** Refuted D1's stated rationale (SLATE is a real cblas
  consumer via STRUMPACK) and the claim that the probe detects ILP64 mismatch; refuted "all four
  flavors" (`cea` is a data directory); corrected the R1 call-site count (three consumers, not four)
  and the R7 `mkl`+`OFF` expectation. Confirmed the tri-state cache semantics, the rank ordering, the
  placement, and the override table. Added: validate the cache value explicitly, clear
  `find_library`/`try_compile` caches on repoint, `include()` the helper before use, probe in C++
  with OpenMP/Fortran/runtime libs passed explicitly. Every refutation was re-verified against the
  tree (`nm -D` on the SLATE and BLAS DSOs, `ls` on `/opt/scls/cea`) before being accepted.
- **Grok round 1 — returned, reconciled.** Independently reproduced the flavor table, the cache
  semantics and D2, then found four things Codex did not: **D4** (a three-way `STREQUAL` treats
  `-DUSE_MKL=off` as AUTO and silently promotes a rejected backend — the round's most serious
  finding, since it defeats the plan's own guarantee), **D5** (`LINK_OPTIONS` is a CMake 3.25
  signature against a 3.11 floor — invisible locally, since this host runs cmake 4.4.2), **D6**
  (`find_mkl.cmake` is the detector and cannot be gated on the value it helps resolve; its
  `FATAL_ERROR` must go), and **D7** (Apple AUTO aborts the configure, because `AUTO` is truthy before
  resolution). Also D8, D10, D11 and the `debug`+AUTO gap in R7. Grok flagged one gate it could not
  run — `readelf -d` on the debug `blaspp` — which Claude ran, closing D1.
- Both auditors converged independently on four items: `cea` is not a flavor, the ILP64 overclaim,
  the R1 call-site miscount, and the R7 `mkl`+`OFF` expectation. Convergence is corroboration here
  only because each was separately checked against the tree.
- Claude self-found, before either auditor returned: D2 (`try_compile` ignores `link_directories()`),
  measured in a standalone CMake project.
- **Second reframe, 2026-08-30 (Christian, no audit round) — the SCLS-flavor presets, R0 and §6.**
  Neither auditor proposed this, and in fairness neither was asked to: the brief they were given
  took the knobs as fixed and asked only how to resolve a backend behind them. Reading
  `$SCLS_FLAVOR` removes most of what round 1 was trying to infer, which is worth recording as a
  method note — **both audits were thorough and both were scoped to the wrong question.** A design
  round asking "what does the environment already tell us?" before "how do we detect it?" would
  have reached the smaller plan first. Unaudited so far, and the parts most worth an audit are the
  `option()`-ordering claim (§6.2b) and the preset/cache interaction (§6.4), both of which are
  measurable and neither of which has been measured yet.
- **Rescope, 2026-08-30 (Christian, no audit round).** Autodetection withdrawn: the build verifies
  the selected backend and errors when it is missing, and never chooses between backends. Recorded
  here rather than as a round because it is a *design decision*, not a finding — no auditor argued
  for it and none needs to check it. Its effect on round-1's output is large and worth stating
  plainly: D4, D6, D7, D8, O1, O2 and O5 are dissolved, R1 is withdrawn, and R7's matrix loses its
  two `AUTO` rows. **None of that makes the audits wasted** — every one of those findings was a real
  defect in the design as it then stood, and D4 in particular is part of why the simpler design is
  attractive. The findings that survive untouched are the ones about the *probe* (D2, D3, D5, D9,
  D10, D11) and about `-lcblas` (D1), which is the work that actually remains.
