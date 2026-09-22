# SCLS portability: the source-build route, and Apple Silicon

**Date:** 2026-09-22
**Purpose:** Session record of working off an external documentation patch plan covering
BELFEM's portability story — SCLS as a source-buildable stack rather than a binary
repository, and Apple Silicon's support level
**Module:** `README.md`, `doc/` (`scls_toolchain.md`, `getting_started.md`, `mpi_support.md`,
`README.md`), `CHANGELOG.md`
**AIs involved:** Claude Opus 5 (1M), Codex (`gpt-5.6-luna`/medium on the two guides,
`gpt-5.6-terra`/medium on the two technical pages)
**Verification:** static. `scripts/check_doc_claims.py` 35/38 — the three failures are
pre-existing debt-register row counts, untouched here. Nothing built, nothing run; the
platform claims rest on Christian's rulings, recorded below.

## Summary

A handoff plan arrived asserting that BELFEM's documentation understates SCLS in two ways:
it reads as though the binary repositories are the portability boundary, and it calls Apple
Silicon unverified. The second half of that was a claim about the world, not about the text,
and the tree did not support it — so it was put to Christian before anything was written.

The first half was correct and is now fixed. The second was correct as well, but only
because of evidence that exists outside both repositories.

## The premise that needed a ruling

The plan asked for Apple Silicon to be described as officially supported. Nothing in either
tree said so. SCLS's own `README.md` and `doc/MACOS_BUILD.md` record a first verified GCC
bootstrap on an M2 Pro on 2026-09-10, reaching OpenMPI in the build order, and state plainly
that everything after it — ScaLAPACK, MUMPS, SCOTCH, the SLATE/STRUMPACK stack, PETSc, VTK —
is unbuilt on arm64. The most recent arm64 commit there is `98f39c5`, the same day.

Christian's rulings:

- **The full stack has since been built on arm64**, and BELFEM was configured, built and its
  test suite run against it. That gate is not recorded in either repository; this entry is
  the first written trace of it.
- **Scope: the BELFEM repository only.** SCLS and the public website stay for a session in
  that tree.

The second ruling leaves a known inconsistency, recorded here so it is not rediscovered as a
defect: BELFEM's documentation now names Apple Silicon a supported platform while SCLS's
`README.md` calls the rest of the stack there "still in bring-up", `doc/MACOS_BUILD.md` says
"macOS support is beta", and `web/scls.html.j2:301` still says "no Apple Silicon build has
been run". The last of those is stale in both directions now. BELFEM's `scls_toolchain.md`
points a reader at `doc/MACOS_BUILD.md` for the macOS prerequisites, which are current; it
deliberately does not point there for the validation status, which is not.

## What changed

**The binary repositories are not the portability boundary.** The old wording — "macOS and
other systems build SCLS from source" — read as though a Linux user without a matching
repository had nowhere to go. SCLS's generic Unix builder exists precisely for that case
(older enterprise Linux, a machine you do not administer, a distribution on neither
packaging path), building the same stack from the same recipes into an isolated prefix.
`README.md`, `doc/getting_started.md` and `doc/scls_toolchain.md` now say so, and the
toolchain page's platform table splits the old "macOS and everything else" row into "any
other Linux" and "macOS (Intel and Apple Silicon)".

**Apple Silicon.** `doc/scls_toolchain.md` carried "macOS support is beta … Apple Silicon is
not yet verified". Replaced with a support statement.

**Open MPI is a stack policy, not a missing portability strategy.** `doc/mpi_support.md`
opened by refusing MPICH and Intel MPI without saying what the refusal buys. It now names the
reason: every SCLS flavor is built around Open MPI on every platform SCLS targets, and PETSc,
MUMPS, STRUMPACK and ScaLAPACK are each compiled against that one MPI. Checked while writing:
`recipes/` in the SCLS tree contains exactly one MPI recipe, `openmpi.yaml`.

**One flavor fact worth stating.** `mkl` is Linux-only — `config/linalg/config_mkl.cmake:1-3`
hard-errors on `APPLE`. The `macos` flavor supplies OpenBLAS (`flavors/macos.yaml`,
`linalg: openblas`), and since the flavor preset keys on the flavor *name*, `USE_MKL` stays
`OFF` there without being asked to. Confirmed by Christian mid-session, twice.

## The Codex round

Two sweeps, depth chosen per protocol §9.1 — `gpt-5.6-luna`/medium for `README.md` and
`doc/getting_started.md` (ordinary guides), `gpt-5.6-terra`/medium for `doc/scls_toolchain.md`
and `doc/mpi_support.md` (dense technical pages). Both stamped in
`tmp/ai_exchange/scls_portability_docs.md`. Neither reported a suspected error.

Most of what came back was applied. Three edits were declined, and the reasons are the
interesting part of the round:

- **An inverted causation.** Codex rewrote the flavor note as "The `macos` flavor supplies
  OpenBLAS, so its name does not contain `mkl` and `USE_MKL` remains `OFF` by default." The
  arrow runs the other way: the preset in `find_scls_flavor.cmake` reads the flavor name, and
  the name is not a consequence of the backend. Kept the original direction.
- **An Oxford comma**, added to "EL10, Amazon Linux 2023 and Ubuntu 24.04". The repository
  does not use one.
- **Two pages made verbatim identical.** Codex proposed the same paragraph for `README.md`
  and `doc/getting_started.md`, down to the closing sentence. The two differ on purpose —
  the README is the denser of the two — so the getting-started version keeps its shorter
  "Activation exports …".

## Not done

- The SCLS repository and `web/scls.html.j2` — deferred by ruling, and the inconsistency
  they now hold is described above.
- No executable gate. The arm64 claim in these documents rests on Christian's statement, not
  on anything this session ran or could run from a Linux x86-64 host.
