# SCLS named in the documentation: the dependency stack was never introduced

**Date:** 2026-09-22
**Purpose:** Christian noticed that nothing in the documentation tells a reader — human or AI —
that BELFEM is built against SCLS. The stack appeared in the tree only as a mechanism for
presetting `USE_MKL`, `USE_PARDISO` and `USE_DEBUG`, never as something to install.
**Module:** documentation (`doc/scls_toolchain.md`, `README.md`, `doc/getting_started.md`,
`CLAUDE.md`, `AGENTS.md`)
**AIs involved:** Claude Opus 5 (1M context); Codex `gpt-5.6-luna` / medium for the language
sweep. No jury round: the session is entirely text, which the 2026-08-28 ruling exempts.
**Status:** COMPLETE. New page written and wired into the page tree, five pointers added, sweep
applied, committed. The commit carries this session's hunks only: `CLAUDE.md`, `CHANGELOG.md`,
`doc/getting_started.md` and `devlog/README.md` each also hold another session's uncommitted
work, which stays in the working tree, and `doc/doxygen_nav.dox` is entirely another session's.
**Verification:** static only. `scripts/check_doc_claims.py` 35/38 (the three failures are
pre-existing `todo/debt_register.md` header counts, untouched here);
`scripts/update_doc_index.py --check` clean. Reviewed, not verified.
**Exchange:** `tmp/ai_exchange/scls_toolchain_doc.md` (language sweep).

## What was missing

`$SCLS` is read in four separate places in the configure — the flavor-derived option defaults,
the link/include/rpath directories, the package search order, and a conditional compiler pin —
and is named in `CLAUDE.md`, `README.md`, `doc/getting_started.md`, the CHANGELOG and a dozen
devlogs. In none of them was it ever *introduced*. A reader met `$SCLS_FLAVOR` as a switch that
changes `USE_DEBUG` without being told that SCLS is a thing you install, where it comes from, or
that it is the environment BELFEM is actually built and tested in.

## Two rulings taken before writing

Put to Christian by question form, both answered with the recommended option:

1. **Claim strength: "strongly recommended, de facto required", not "required".** The CMake is
   permissive — `find_scls.cmake` runs only `if( DEFINED ENV{SCLS} )`, and `belfem_find_package`
   falls back to `/usr` and `/usr/local` — so a hard "required" would be contradicted by the code
   a reader can check. The documentation says instead that no other configuration is tested, and
   why an inconsistent stack fails inside a solver rather than at the link stage.
2. **Scope: a dedicated page plus pointers**, mirroring `doc/mpi_support.md`, rather than
   scattering the same paragraph across three files.

## What landed

`doc/scls_toolchain.md` (new, anchor `doc_scls_toolchain`): what SCLS supplies, the
required-or-not question answered honestly, install and activate per platform, the three public
flavors, the four mechanisms that read `$SCLS`, and what building without it demands. Two facts
in it are worth having written down: `scripts/scls_env.sh` is *not* an SCLS script and is easy to
confuse with `activate`, and the Linux binary flavors ship no compiler binary, so the compiler-pin
block in `CMakeLists.txt` is dead on Linux by design and exists for macOS.

Pointers: `README.md` ("Building" now opens on the stack, before the `cmake` block),
`doc/getting_started.md` (new §1.1 ahead of the configure step), `CLAUDE.md` (top of "Build
System"), `AGENTS.md` (routing line), `doc/README.md`, and the regenerated `doc/mainpage.md` /
`doc/doxygen_nav.dox` via a new entry in `update_doc_index.py`'s `PAGE_GROUPS`. CHANGELOG entry
under Unreleased/Added.

## What the sweep caught

Codex (luna/medium, read-only) returned the usual sentence-splitting, most of which was applied,
and one genuine contradiction the prose pass was not asked for: the new page said Open MPI is
"mandatory regardless" while `README.md` and `doc/getting_started.md` kept their older "for
parallel runs" phrasing. `find_mpi.cmake` is gated on `USE_MPI`, which defaults `ON` — so a stock
configure does require Open MPI, and "for parallel runs" was the looser of the two. All three now
say the same thing and name the option. Codex also noticed that CMake was promised in two of the
three package lists and absent from the third; added.

Declined: several rewrites that flattened the direct voice the house style keeps ("Install
SCLS.", "and is nobody's mistake"), and suggestions against two pre-existing paragraphs that this
session did not touch.

## Not done

No build was run, and none is owed — nothing executable changed. The `update_doc_index.py` edit
is one string in a presentation table. Regenerating the page tree also moved a peer's hand-added
`@subpage doc_examples_rutherford` line into its sorted position; same content, generator
ordering.
