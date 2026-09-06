# Repairing `make doc`: Four Version-Skew Breaks in the Doxygen Front-End

**Date:** 2026-08-28
**Purpose:** Record the diagnosis and repair of `make doc`, which aborted before producing a site.
**Module:** build system (`CMakeLists.txt`, `cmake/`, `doc/`)
**AIs involved:** Claude (diagnosis, implementation), Codex + Grok (plan audit and code audit)
**Plan:** `todo/make_doc_repair.md`

## Symptom

`make doc` printed 46 `ignoring unsupported tag` warnings and then died:

```
CMake Error at cmake/patch_doxygen_header.cmake:28 (message):
  patch_doxygen_header: the project-logo markup in the generated header no
  longer matches what this script expects.
make: *** [Makefile:244: doc] Error 2
```

## Root cause

One cause with four faces: **every Doxygen input artifact in the tree was authored against Doxygen
1.18.0, and the only Doxygen installed is 1.9.1** (`/usr/bin/doxygen`, RPM
`doxygen-1.9.1-12.el9_5`). `Doxyfile.in:1` became `# Doxyfile 1.18.0` in commit `bccd5b03`
(2026-08-18); it had been `# Doxyfile 1.8.16` before.

| | Break | Mechanism |
|---|---|---|
| B1 | fresh-tree abort | `doxygen -w html` was invoked with no config argument (`CMakeLists.txt:496-499` as it stood). Doxygen still reads `./Doxyfile` from the working directory and validates `HTML_HEADER` for existence — and `HTML_HEADER` names the very file that command creates. Chicken-and-egg. |
| B2 | the observed abort | the logo-patch needle contained 1.18's `$logosize` marker, which 1.9.1 does not emit. |
| B3 | silent feature loss | 46 tags ignored, including `MATHJAX_VERSION = MathJax_3` — the site was silently on the MathJax **2** CDN. |
| B4 | the run's only `error:` | `doc/DoxygenLayout.xml` is a 1.18 layout dump; 1.9.1 rejects `$SHOW_HEADERFILE` and ~32 of its elements. |

**Which break fires depends on the build tree.** With no `doxygen_header.html`, B1 aborts at command
two. With one present, `-w` succeeds, overwrites it with a 1.9.1 template, and B2 aborts at command
three. Both are real; the reported failure was the B2 path.

The documentation itself was never broken. With B1 and B2 bypassed by hand, 1.9.1 produced 3194
pages with the full hand-written page tree intact.

## Repair

The `doc` target now **derives every Doxygen input from the installed Doxygen** instead of assuming
a version. That principle was already stated in `cmake/patch_doxygen_header.cmake` for the HTML
header; this extends it to the configuration and the layout.

- **B1** — `doxygen -w html` is handed a config copy with `HTML_HEADER` blanked
  (`cmake/strip_doxygen_html_header.cmake`). Generating config-less was rejected: Doxygen's own
  commentary in `Doxyfile.in:1433-1438` says the template depends on configuration settings.
- **B3** — `cmake/normalize_doxyfile.cmake` copies the pristine `configure_file` output and runs
  `doxygen -u` on the **copy**. The copy is not cosmetic: `make doc` does not re-run CMake, so
  normalising in place would leave a permanently downgraded file, and a later Doxygen upgrade would
  then re-add the newer tags at their defaults instead of the values from `Doxyfile.in`. The script
  also summarises the 46 unavoidable `-u` warnings to one line while passing anything else through.
- **B4** — `cmake/normalize_doxygen_layout.cmake` generates the installed version's default layout
  and re-applies BELFEM's delta, **read from `doc/DoxygenLayout.xml`** rather than hardcoded, then
  repoints `LAYOUT_FILE` at the result. The checked-in layout is only ever read.
- **B2** — the logo needle became a regex matching both the unpatched and already-patched forms,
  making the rewrite idempotent by construction.

Content defects cleared along the way: four `$$…$$` blocks in
`thermal_expansion_from_heat_capacity.md` became `\f[…\f]` (Doxygen has no `$$` display math — the
equations had been rendering as **literal LaTeX source** on the site); two duplicate `@param` blocks
on `Bezier::d3ydx3`/`d3xdy3` demoted to plain comments; `::Oxygen` spelled out; three `EF_*::` table
references in `coulomb_gauge_penalty_theory.md` given Doxygen's `%` no-link escape.

## Result (measured, Doxygen 1.9.1)

| | before | after |
|---|---|---|
| `make doc` | **fails** | exit 0 |
| console tag warnings | 46 | 0 (one summary line) |
| log `error:` | 1 | **0** |
| log `warning:` | 200 | **90** |
| warnings outside one file | ~20 | **0** |

All 90 remaining warnings are in `doc/lessons_learned_evidence.md` and are a 1.9.1 markdown-table
limitation, not bad markup — established by three independent checks: backtick parity is even on
every line, a synthetic table with the same constructs emits nothing, and bisection shows the
warnings start ~40 rows into a single 110-row table (line 132 alone: clean). That file is exempt
from prose sweeps, so it is quarantined in the baseline rather than rewritten.

## What the audits caught

Two rounds, Codex and Grok, both read-only. Every finding was re-verified against the tree before
acceptance. The audits materially changed the work, and three findings corrected Claude's own errors:

- **Fabricated identifiers (Grok).** The plan cited `Bezier::dNydxN` / `dNxdyN`, which exist nowhere
  in the tree. They were an artifact of Claude's own `sed 's/[0-9]+/N/g'` while tabulating the
  warning log. Real names: `d3ydx3` / `d3xdy3`.
- **False CMake premise (both).** The plan justified the pristine/working split by claiming
  `configure_file` only rewrites when its input changes. It rewrites on **every** configure —
  confirmed by a scratch project that restored a deliberately mutated output. The split is still
  right, for the narrower reason now recorded.
- **A misleading baseline comment (Grok).** The first re-baselining described the non-`</tt>`
  warnings inaccurately and quietly buried genuine unresolved-reference defects. Those were fixed
  instead of baselined, which is why the number is 90 and not 100.
- **`LAYOUT_FILE` left pointing at the source tree (both)** — without the repoint, the generated
  layout would have been written and ignored.
- **List-splitting in the layout normaliser (Codex).** `string(REGEX MATCHALL)` returns a CMake
  list, so a `;` in a nav-tab title split it into two corrupt fragments (reproduced: two tabs in,
  three items out). Replaced with a consume-loop; the `;` case now round-trips.
- **A too-loose idempotence guard (Codex).** It returned on any `<a href=`, so a stale or wrong URL
  bypassed both the fix and the loud failure. The regex now corrects such a header.
- **Silent nav loss (both).** The no-tabs path warned instead of failing, and the insertion was
  never verified. Both are now `FATAL_ERROR` — losing the custom navigation is exactly what a
  custom layout file exists to prevent.

Grok also established that `includedbygraph` is deliberately *not* carried over: the checked-in
layout hardcodes `visible="yes"`, which was overriding `INCLUDED_BY_GRAPH = NO` in the Doxyfile (set
2026-08-07 over a 362-node graph). Leaving it at the generated default restores the Doxyfile's
intent, and the "Included by graph … too many nodes" warnings are gone.

## Still open

- **O5 — the warning gate cannot fail a build.** `cmake/report_doxygen_warnings.cmake` only calls
  `message(WARNING)`, and `cmake -P` exits 0 on a warning, so a regression past the baseline prints
  and is ignored. Making it fatal is a release-policy decision, left to Christian.
- **O6 — `make doc` writes tracked files.** `update_doc_index.py` runs without `--check`, so
  building the docs can dirty the worktree and needs a writable checkout.
- **O7 — none of this is verified on Doxygen 1.18**, which does not exist on this machine. The
  cross-version behaviour is designed-for, not proven. In particular, whether `doxygen -l` loads
  `./Doxyfile` on 1.18, and whether 1.18's `-w` template is genuinely config-sensitive, are untested.

Per the evidence ladder this work is **verified** for the pipeline itself (the full sequence was
executed from both the header-absent and header-present states and the site inspected) and
**reviewed** for the 1.18 half.
