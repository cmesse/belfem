# Repair `make doc`: Doxygen Front-End Breaks and the 1.18-vs-1.9.1 Version Skew

**Date:** 2026-08-28
**Purpose:** `make doc` aborts before it produces a site. The cause is a chicken-and-egg in the
custom-header generation, compounded by a build environment whose only Doxygen is 1.9.1 while every
Doxygen input artifact in the tree (`Doxyfile.in`, `doc/DoxygenLayout.xml`,
`cmake/patch_doxygen_header.cmake`) was authored against 1.18.0. The fix makes the `doc` target
derive its Doxygen inputs from whatever Doxygen is installed, so the target works on 1.9.1 today and
on 1.18 after an upgrade, without pinning either.
**Module:** build system (`CMakeLists.txt`, `cmake/`, `doc/`)
**AIs involved:** Claude (diagnosis + plan), Codex (audit), Grok (third voice)
**Status:** ✅ **COMPLETE 2026-08-28** for the 1.9.1 half. R1-R7 landed and verified by execution:
`make doc` goes from a hard abort to exit 0, console tag warnings 46 -> 0, log 200 warnings + 1 error
-> 90 + 0, with every remaining warning confined to `doc/lessons_learned_evidence.md` and none
anywhere else in the tree. Two audit rounds (plan, then code) with Codex and Grok; all findings
re-verified against the tree before acceptance and folded in. **Still open: O5** (the warning gate
only `message(WARNING)`s, so it cannot fail a build - a release-policy call), **O6** (`make doc`
writes tracked files), and **O7** (nothing here is verified on Doxygen 1.18; no such binary exists on
this machine). Devlog: `devlog/dl20260828_make_doc_repair.md`.

> **Scope guards:**
> - IN scope: the four mechanisms that stop `make doc` from producing a clean site (§1).
> - OUT of scope: installing or building a newer Doxygen; changing what the site looks like;
>   rewriting `Doxyfile.in` for a specific Doxygen version (rejected by Christian 2026-08-28 in
>   favour of build-time normalisation).
> - OUT of scope: the prose of `doc/lessons_learned_evidence.md`. It is on CLAUDE.md's no-sweep
>   list, and O1 concludes it carries no defect to fix.
> - Compatibility promise kept: the target must not require any particular Doxygen version.
> - **Not promised:** that the result is *verified* on 1.18. Only 1.9.1 exists on this machine, so
>   every 1.18 statement below is an assumption, marked as one. See O7.

---

## 1. Current Behaviour and How It Fails

The `doc` target (`CMakeLists.txt:494-510`) runs five commands sharing one working directory
(`:508`; confirmed by both auditors against the generated recipe in
`cmake-build-debug/CMakeFiles/doc.dir/build.make:69-75`): refresh the markdown page index, generate
the HTML header/footer/stylesheet templates, patch the header to make the logo a link, run Doxygen,
then summarise the problem log. The first command is dropped entirely when Python is absent
(`CMakeLists.txt:477-484`), so "five" is conditional.

All findings below were reproduced on 2026-08-28 against `/usr/bin/doxygen` 1.9.1 (RPM
`doxygen-1.9.1-12.el9_5.x86_64`, the only Doxygen on the machine).

**Which command dies depends on the state of the build tree, and both cases are real:**

- **Fresh tree** (no `doxygen_header.html`): dies at command 2 on **B1**.
- **This worktree today**: `cmake-build-debug/doxygen_header.html` exists — but only because Claude
  created it by hand during this diagnosis. `-w` therefore succeeds, overwrites it with a 1.9.1
  template, and command 3 dies on **B2**. (Found by Grok; the v1 plan's "currently holds no header"
  evidence cell was stale by Claude's own doing.)

| # | Failure | Mechanism | Evidence |
|---|---|---|---|
| **B1** | **Hard abort, exit 1** on a fresh tree | `doxygen -w html` is invoked with no config argument and `WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}` (`CMakeLists.txt:496-499,508`). With no config on the command line Doxygen still reads `./Doxyfile` from the working directory and validates `HTML_HEADER = @CMAKE_BINARY_DIR@/doxygen_header.html` (`Doxyfile.in:1459`) for existence. The command that *creates* the header cannot run because the header does not exist. | `error: tag HTML_HEADER: header file '…/doxygen_header.html' does not exist / Exiting...`, exit 1. The same command in a directory with no `Doxyfile` succeeds. Independently reproduced by Codex in `cmake-build-claude`. |
| **B2** | **Hard abort (deliberate FATAL_ERROR)** — the logo patch cannot find its anchor | `_needle` (`cmake/patch_doxygen_header.cmake:23`) expects `<td id="projectlogo"><img alt="Logo" src="$relpath^$projectlogo"$logosize/></td>`. `$logosize` is a 1.18-ism; 1.9.1 emits the same markup **without** it. The guard at `:25-33` fires as designed. | 1.9.1 header line 29 lacks `$logosize`. Running the script against it yields the scripted FATAL_ERROR. |
| **B3** | **46 console warnings, silent feature loss** — `Doxyfile.in` is a 1.18.0 template | `Doxyfile.in:1` reads `# Doxyfile 1.18.0`; it was `# Doxyfile 1.8.16` until commit `bccd5b03` (2026-08-18). 1.9.1 ignores 46 tags it does not know. | `ignoring unsupported tag 'X'` × 46. Casualties include `MATHJAX_VERSION = MathJax_3` (`:1961`) — the site falls back to the MathJax 2 CDN, confirmed in the generated `index.html:26` — plus `HTML_COLORSTYLE`, `HTML_COPY_CLIPBOARD`, `HTML_CODE_FOLDING`, `PAGE_OUTLINE_PANEL`, `SHOW_HEADERFILE`. **Nuance (Grok):** `MATHJAX_RELPATH` is deliberately blank (`:1997-2005`) precisely so an older Doxygen cannot 404 a v3 path, so the v2 fallback is partly by design. |
| **B4** | **1 Doxygen `error:` and 32 log warnings** — `doc/DoxygenLayout.xml` is also 1.18-authored | The layout uses elements and variables 1.9.1 does not know: `<includes visible="$SHOW_HEADERFILE"/>` (`doc/DoxygenLayout.xml:52`), plus `topics`/`modulelist`/`modulemembers`/`concepts` navindex types and `<properties>` at `:125,:137,:173,:186`. | `error: found unsupported value $SHOW_HEADERFILE …`; 4 navindex warnings + 28 `Unexpected start tag`. **Severity correction (Grok):** this is log noise, **not** functional loss — the three custom nav tabs still render on 1.9.1 (`menudata.js:26-29` carries Getting Started / Modules / Documentation). The v1 plan also over-cited `:77,:100`; those class-level `<properties>` are accepted by 1.9.1 and did not warn. |

### 1.1 What is *not* broken

With B1 and B2 worked around by hand, Doxygen 1.9.1 runs to completion: **exit 0**, 3194 HTML pages,
navigation treeview present, logo link applied, the three custom nav tabs intact, and the whole
hand-written page tree rendered. `HAVE_DOT` is on with a real Graphviz, and `HTML_FOOTER` is blank
(`Doxyfile.in:1469`) so there is no second chicken-and-egg on the footer.

**Caveat (Codex, verified):** `scripts/update_doc_index.py` is invoked *without* `--check`
(`CMakeLists.txt:479-481`) and writes tracked files in place — markdown anchors,
`doc/doxygen_nav.dox`, `doc/groups.dox`, `doc/mainpage.md` (`scripts/update_doc_index.py:440-445`).
So `make doc` requires a writable checkout and can dirty the worktree. Whether that is intended is
**O6**. The generator is otherwise healthy (`26 module(s), 97 page(s)`, no `updated:` lines on the
current tree) and is cwd-safe (`REPO` derived from `__file__`, `:36-38`), but it *can* fail the
target first if a module or page is unfiled (`:53-55,63-64`).

### 1.2 The warning baseline is blown, and the console does not show it

`cmake/report_doxygen_warnings.cmake:20-21` sets both baselines to 0 and records the log as clean as
of 2026-08-18, so "any warning or error is a regression". Under 1.9.1 the log carries **200 warning
lines and 1 error**. Note the trap the script's own comment already anticipates at `:14-16`:
configuration complaints go to **stderr**, not `WARN_LOGFILE`, so B3's 46 tag warnings and these 200
log lines are two disjoint sets. A clean console is not evidence of a clean run, or vice versa.

**Counting correction (Grok, verified):** Doxygen logs the same defect once per pass, so log lines
overcount distinct defects. Baselining is unaffected; a *fix list* built from raw counts is wrong.

| Log lines | Unique | Warning | Source | Class |
|---|---|---|---|---|
| 1 error + 32 | — | layout elements/variables | `doc/DoxygenLayout.xml` | B4 — version skew |
| 80 | 49 | `found </tt> tag without matching <tt>` | all in `doc/lessons_learned_evidence.md` | **O1 — not a source defect** |
| 60 | 15 | `Found unknown command '\alpha' '\frac' …` | all in `src/physics/materials/doc/thermal_expansion_from_heat_capacity.md` | **O2 — genuine content defect** |
| 14 | — | `explicit link request … could not be resolved` | various (incl. `coulomb_gauge_penalty_theory.md` E/C/G ×3) | content defect |
| 6 | — | `Included by graph … too many nodes` | generated | benign threshold |
| 2 | 2 | `multiple @param documentation sections` | `belfem::Bezier::d3ydx3`, `d3xdy3` | genuine content defect |

**Bottom line:** `make doc` is not broken by anything in the documentation content — it is broken
because the target asks Doxygen to validate a file that the same command is about to create (B1),
and because three generated-input artifacts were authored against a Doxygen that is not installed
(B2, B3, B4).

## 2. Architecture: Derive Every Doxygen Input from the Installed Doxygen

Christian's decision (2026-08-28): **auto-normalise at build time**, rather than downgrading
`Doxyfile.in` to 1.9.1 or requiring Doxygen ≥ 1.18.

Doxygen ships the migration tools this needs: `doxygen -u <config>` rewrites a config to the running
version's tag set, and `doxygen -l <layout>` writes that version's default layout. Measured on
2026-08-28, `-u` against the configured Doxyfile: rewrites the header to `# Doxyfile 1.9.1`, drops
all 46 unknown tags, adds the 1.9.1-only tags at their defaults, and **preserves every deliberate
setting**, including the CMake-substituted absolute paths (`PROJECT_NUMBER`, `WARN_LOGFILE`,
`HTML_HEADER`, `GENERATE_TREEVIEW = YES`, `USE_MATHJAX = YES`, `MATHJAX_FORMAT = HTML-CSS`). A
second Doxygen read of the normalised file is then silent.

**Rejected alternatives:** regenerating `Doxyfile.in` for 1.9.1 (re-breaks on the next upgrade and
discards what `bccd5b03` deliberately added); requiring Doxygen ≥ 1.18 (the el9 RPM is 1.9.1, so
`make doc` would be unavailable on a stock RHEL 9 box — the wrong default for a public release).

### 2.1 Why the working Doxyfile must be a copy — corrected after audit

> **CORRECTED 2026-08-28.** The v1 plan justified the pristine/working split by claiming
> `configure_file` "only rewrites when `Doxyfile.in` changes", making the 1.18 content
> "unrecoverable". **Both auditors refuted this, and Claude's own probe confirmed the refutation:**
> a scratch CMake project shows `configure_file` regenerates from the input and overwrites a
> mutated output on **every** configure, input untouched (`PRISTINE` restored after deliberate
> mutation). The CMake docs say the output is rewritten on subsequent runs whenever content differs
> (`configure_file.rst:21-25`). The 1.18 content is never unrecoverable — it lives in `Doxyfile.in`.

The real, narrower risk stands and still forces the split: **`make doc` does not re-run CMake.** An
`add_custom_target` is always out of date but is not a reconfigure, and the Doxygen executable's
version is not a configure dependency. So if `-u` mutates the configure-time output in place, the
downgraded file persists across repeated `make doc` runs — and after an in-place Doxygen upgrade to
1.18 without an intervening `cmake` run, `-u` would upgrade the already-lossy 1.9.1 file and fill
the restored 1.18-only tags with **defaults** rather than the values from `Doxyfile.in`. (That last
step is an assumption about 1.18 behaviour, not verified here.)

**Therefore: never `-u` the configure-time output.** Copy it to a working name first, and normalise
the copy. Grok's simpler form is adopted: leave `configure_file` writing `Doxyfile` as it does
today, and have the `doc` target's first Doxygen-related command copy it to a working file. Same
invariant, no change to `CMakeLists.txt:473`, fewer moving parts than renaming the destination. A
build-time `configure_file` is not an option — `configure_file` is configure-time only.

## 3. Gap Table

| # | Artifact | Needed for | Handled today? | Class | Citation / rationale |
|---|---|---|---|---|---|
| 1 | HTML header template | logo link, treeview, search, MathJax hooks | **no** — B1 blocks creation | (c) explicit | `CMakeLists.txt:496-499` |
| 2 | logo-link patch anchor | linked project logo | **no** — B2, needle is 1.18-only | (c) explicit | `cmake/patch_doxygen_header.cmake:23,36` |
| 3 | patch script's *recovery instructions* | recovering from a future template change | **no** — the command it tells you to run is itself B1 | (c) explicit | `cmake/patch_doxygen_header.cmake:31-32` (Grok) |
| 4 | working (normalised) Doxyfile | the whole run | partly — 46 tags ignored | (a) rebuildable via copy + `doxygen -u` | `Doxyfile.in:1`; §2 measurement |
| 5 | pristine Doxyfile untouched by `-u` | correctness after a Doxygen upgrade with no reconfigure | **no** | (c) explicit | §2.1 |
| 6 | working layout file | navindex, member ordering | **no** — B4 | (c) explicit, see **O3** | `doc/DoxygenLayout.xml`; `Doxyfile.in:894` |
| 7 | `LAYOUT_FILE` pointing at the working layout | making #6 take effect at all | **no** — points into the source tree | (c) explicit | `Doxyfile.in:894` (both auditors) |
| 8 | warning/error baselines | regression gate | stale — 0/0 vs actual 200/1, and the gate does not fail | (c) explicit, see **O5** | `cmake/report_doxygen_warnings.cmake:20-21,35-46` |
| 9 | `$$…$$` display math | rendered equations | **no** — renders as literal text | (c) explicit, **O2 resolved** | `thermal_expansion_from_heat_capacity.md` |
| 10 | duplicate `@param` blocks | clean member docs | **no** — genuine defect | (c) explicit | `cl_Bezier.hpp:233,243` + `cl_Bezier.cpp:308-348` |
| 11 | worktree cleanliness during `make doc` | building from a read-only or clean checkout | **no** — writes tracked files | (b) open, see **O6** | `scripts/update_doc_index.py:440-445` |

### 3.1 Cross-cutting finding

Rows 2, 3, 4 and 6 are one defect wearing four hats: **a Doxygen input artifact was authored against
a Doxygen the tree does not have.** The design answer is uniform — derive every Doxygen input from
the installed Doxygen at build time, never check in a version-bound one — and it is the rule
`cmake/patch_doxygen_header.cmake:6-13` already states for the header. R1–R4 extend that existing
rule to the config and the layout; they introduce no new principle.

## 4. Ordered Steps

- [x] **R1 — Fix B1 by giving `-w html` a sanitised config (variant (b)).** **O4 resolved.**
      Write a working config copy with `HTML_HEADER` blanked, and pass it as the 5th argument to
      `doxygen -w html`. This is Doxygen's own recommended invocation and the only variant that is
      both config-aware and immune to the existence check.
      *Must* blank the tag — passing the real Doxyfile as the 5th argument is still B1.
- [x] **R2 — Make the logo-patch anchor version-tolerant, and fix its recovery text (B2).**
      (after: R1) Replace the literal `_needle` (`cmake/patch_doxygen_header.cmake:23,36`) with a
      match accepting the markup with **and** without `$logosize`, keeping the FATAL_ERROR for a
      genuine template change. Also fix the recovery command at `:31-32`, which as written
      reproduces B1 when run from the build directory (Grok). Update the comment block to document
      the tolerance rather than one version's markup.
- [x] **R3 — Copy the configured Doxyfile to a working name, then `doxygen -u` the copy (B3).**
      Never normalise the configure-time output in place (§2.1). Ordering versus R1 is **not**
      constrained by B1 — see the correction below — but `-u` should precede the main Doxygen run
      so that run is free of tag noise.
- [x] **R4 — Derive a working layout from the installed Doxygen, and retarget `LAYOUT_FILE` (B4).**
      (after: R3) **O3 resolved:** `doxygen -l` into the build tree, then re-apply BELFEM's delta —
      exactly six lines: the three `<tab type="user">` entries (`doc/DoxygenLayout.xml:7-9`, whose
      anchors `doc_getting_started`, `docmap`, `doc_index` all resolve) and the three forced graph
      visibilities (`:53,:54,:265`). Two hard constraints from both auditors: **never mutate
      `doc/DoxygenLayout.xml`**, and **rewrite `LAYOUT_FILE` in the working Doxyfile** — after `-u`
      it still points at the source file (`Doxyfile.in:894`), so without this the working layout has
      no effect.
- [x] **R5 — Fix the genuine content defects.** (no dependency on R4 — Grok correctly flagged the
      v1 ordering as artificial)
      - [x] Convert the 4 `$$…$$` blocks in `thermal_expansion_from_heat_capacity.md` to `\f[ … \f]`
            (**O2**).
      - [x] Remove the duplicate `@param` blocks for `belfem::Bezier::d3ydx3` / `d3xdy3` — the
            header already documents them (`cl_Bezier.hpp:233,243`), the second copy is in
            `cl_Bezier.cpp:308-348`.
      - [x] The unresolved `\ref` warnings (14 log lines, fewer unique).
- [x] **R6 — Re-baseline the regression gate.** (after: R5, and **blocked on O5**)
      Set `BELFEM_DOXYGEN_WARNING_BASELINE` to what actually remains and record which Doxygen
      version produced it — the baseline is version-bound and the file does not say so.
      **`_ERROR_BASELINE` stays 0** (Codex): raising it contradicts the script's own invariant that
      errors mean missing output (`report_doxygen_warnings.cmake:33-39`), and R4 should clear the
      only error. If O1's 49 and O2's residue land in the baseline, say so explicitly in the
      comment rather than letting it happen silently.
- [x] **R7 — Gate.** (after: R6) Christian runs `make doc`. Pass criteria, strengthened after audit:
      exit 0 from a tree with `doxygen_header.html` **deleted** (the B1 case) *and* from a tree
      where it exists (the B2 case); no `ignoring unsupported tag` on the console; log at or below
      the R6 baseline with **0 errors**; the patched header still contains `$treeview` and
      `$mathjax`; `doc/html/index.html` present with the linked logo, the treeview, and the three
      custom nav tabs in `menudata.js`; and a second consecutive `make doc` behaves identically
      (copy-from-pristine idempotence).
      **R7 certifies 1.9.1 only.** It cannot prove the 1.18 half of the compatibility promise (O7).

## 5. Open Design Questions

- **O1 — RESOLVED 2026-08-28 → not a source defect; document, do not edit.** The 80 `</tt>` log
  lines (49 unique) in `doc/lessons_learned_evidence.md` are a 1.9.1 markdown-table limitation, not
  bad markup. Evidence: backtick parity is even on every line in the file (0 odd-count lines); a
  synthetic table reproducing the same constructs (unicode, multiple spans per cell) emits no
  warning; and bisection shows lines 1–131 are clean while the warnings begin ~40 rows into a single
  110-row, 11-column table, needing table context to fire (line 132 alone: 0 warnings). Fixing this
  would mean editing a file CLAUDE.md explicitly exempts from sweeps. Fold into the R6 baseline.
- **O2 — RESOLVED 2026-08-28 → genuine content defect on every version, fix in R5.** Probe: a
  minimal page containing both forms rendered `$$\alpha(T) = \frac{1}{L}\frac{dL}{dT}$$` as
  **literal text** in the HTML (and emitted the unknown-command warnings), while `\f[ … \f]`
  rendered as a proper MathJax block with no warnings. Doxygen's documented display-math syntax is
  `\f[…\f]`. Blast radius is small: `$$` appears in exactly one Doxygen-visible file, 4 occurrences.
  (1.18 behaviour still unverified, but "very likely wrong on every version" — Grok.)
- **O3 — RESOLVED 2026-08-28 → `doxygen -l` + re-apply six lines; see R4.** Diffing the checked-in
  layout against the stock 1.9.1 layout shows it *is* the stock 1.18 default (`version="2.0"`,
  "Generated by doxygen 1.18.0") plus only: 3 `<tab type="user">` entries and 3 forced graph
  visibilities. Everything else is 1.18-vs-1.9.1 default drift. Option (iii) — dropping
  `LAYOUT_FILE` — is **rejected**: it would silently ship a site whose top bar loses Getting
  Started / Modules / Documentation, and R7 as written in v1 would not have caught it (Grok).
- **O4 — RESOLVED 2026-08-28 → variant (b).** Both auditors independently chose it. The deciding
  evidence is in the tree: Doxygen's own commentary at `Doxyfile.in:1433-1438` states the header
  "is dependent on the configuration options used (e.g. the setting `GENERATE_TREEVIEW`)" and
  recommends passing `YourConfigFile`. Claude's measurement that the 1.9.1 template is
  byte-identical with and without a config licenses variant (a) **on 1.9.1 only**; treating it as
  version-independent was an over-read. Variant (c) (`touch`) is rejected — it assumes the check is
  existence-only and can leave a zero-byte header that silently poisons a later run.
- **O5 — OPEN, needs Christian: should `make doc` *fail* when the log exceeds baseline?**
  Codex rates this CRITICAL: `report_doxygen_warnings.cmake:35-46` only calls `message(WARNING)`,
  and `cmake -P` exits 0 on a warning, so the "regression gate" never fails anything. Today a
  baseline breach is invisible in an exit code. Options: (i) leave advisory (status quo, honest
  about being a report not a gate); (ii) `FATAL_ERROR` above baseline, making `make doc` fail;
  (iii) fail on `error:` only, warn on `warning:`. This is a policy call about how strict the
  release build should be, so it is not being decided here. **R6 is blocked on it.**
- **O6 — OPEN, needs Christian: should `make doc` be allowed to write tracked source files?**
  It runs `update_doc_index.py` without `--check` (`CMakeLists.txt:479-481`), so building the docs
  can dirty the worktree and requires a writable checkout. The script already supports `--check`
  (returns 1 if anything is pending). Options: (i) keep as-is — regenerating the nav is the point;
  (ii) use `--check` in the `doc` target and make refreshing an explicit separate step. Interacts
  with packaging and any read-only/CI build.
- **O7 — OPEN, cannot be closed on this machine: is the repair correct on Doxygen 1.18?**
  Every 1.18 claim here is an assumption. R7 certifies 1.9.1 only. Closing this needs a 1.18
  binary; until then the compatibility promise is "designed for", not "verified on".

## 6. Definition-of-Done Checklist

- [x] Every gap-table row mapped to a step or an open question.
- [x] Each claimed failure backed by a reproduction, not an assumption — §1 done 2026-08-28.
- [x] Ordered steps with dependencies; **no R3-before-R1 constraint** (see correction below).
- [x] O1–O4 resolved by measurement; O5/O6 answered by Christian; O7 acknowledged as open.
- [x] `make doc` passes R7 in **both** tree states (header present and header deleted).
- [x] Devlog written and `devlog/README.md` updated.

> **CORRECTED 2026-08-28 (Grok, CRITICAL).** The v1 plan required R3 before R1 "so `-w html` reads
> an already-clean config". That rationale is wrong: unknown tags only *warn*
> (`WARN_AS_ERROR = NO`, `Doxyfile.in:1004`), so normalising tags does nothing about the existence
> check that actually aborts. R1's sanitised config is the B1 fix; R3 is not.
> Grok further worried (medium confidence) that `doxygen -u` might itself validate `HTML_HEADER`,
> which would make R3 a second B1. **Measured and refuted:** `doxygen -u` on a Doxyfile whose
> `HTML_HEADER` names a genuinely absent file exits **0**, while `-w html` on the same file exits
> **1**. Confirmed independently by Codex via a stdin probe. The check is specific to the `-w` path.

## 6.1 Code Audit Round — 2026-08-28 (Codex + Grok, read-only)

Run on the implementation diff after R1-R7 landed. No CRITICAL defect; neither auditor found a
correctness bug in the 1.9.1 pipeline. Every finding below was reproduced before being acted on.
Fixes applied in the same session:

- [x] **D1 (Codex, HIGH) — `string(REGEX MATCHALL)` splits on `;`.** A semicolon in a nav-tab title
      became a CMake list separator and tore the tab into two corrupt fragments. Reproduced: two tabs
      in, three items out. Replaced with an explicit consume-loop; `title="Getting; Started"` now
      round-trips intact. *Fixed 2026-08-28, verified by re-running the normaliser on a doctored
      layout.*
- [x] **D2 (Grok, HIGH) — custom nav tabs could be dropped silently.** The main-page insertion used
      `REGEX REPLACE` and was never checked, so a `doxygen -l` emitting a differently-shaped mainpage
      tag would lose Getting Started / Modules / Documentation and still exit 0. Now: literal
      `string(REPLACE)`, a `FATAL_ERROR` if the anchor is absent, and a post-insertion assertion.
      *Fixed 2026-08-28.*
- [x] **D3 (both, HIGH) — the no-tabs path warned instead of failing.** `cmake -P` exits 0 on
      `message(WARNING)`, so the site could ship with no custom navigation. Now `FATAL_ERROR`.
      *Fixed 2026-08-28.*
- [x] **D4 (Grok, HIGH) — the re-baseline comment was false, and hid real defects.** It described the
      20 non-`</tt>` warnings as unresolved `\ref`/`\includedoc` plus a stray `<module>`; the actual
      composition was 9 coulomb-table lines, 1 `::Oxygen`, and 10 more from
      `lessons_learned_evidence.md`. The fixable ones were **fixed rather than baselined** (`%`
      no-link escape on three `EF_*::` table refs; `HelmholtzModel::Oxygen` spelled out), which is why
      the baseline is 90 and not 100, and why nothing outside one file now warns. *Fixed 2026-08-28.*
- [x] **D5 (Codex, HIGH) — the idempotence guard was too loose.** It returned on any `<a href=`, so a
      stale or wrong URL bypassed both the rewrite and the loud failure. Replaced with one regex
      matching the unpatched *and* patched forms, making the rewrite idempotent by construction and
      correcting a wrong URL. Retested on four cases including a stale link and a structural change.
      *Fixed 2026-08-28.*
- [x] **D6 (Grok, LOW) — tag regexes rejected leading indent.** Nothing guarantees a future
      `doxygen -u` left-aligns assignments; an indented `HTML_HEADER` would have reinstated B1. Both
      regexes now tolerate indent, and both scripts assert the substitution landed. *Fixed 2026-08-28.*
- [x] **D7 (Grok, LOW) — `includedbygraph` deliberately not carried over.** The checked-in layout
      hardcodes `visible="yes"`, overriding `INCLUDED_BY_GRAPH = NO` in the Doxyfile (set 2026-08-07
      over a 362-node graph on `typedefs.hpp`). Leaving it at the generated default restores the
      Doxyfile's intent; the "Included by graph ... too many nodes" warnings are gone. Reasoning
      recorded in the script so it is not "fixed" later by mistake. *Documented 2026-08-28.*
- [ ] **D8 (both, HIGH/LOW) — the gate still cannot fail a build.** `message(WARNING)` only. This is
      **O5**, a policy decision for Christian, deliberately not taken here.

Confirmed correct by both auditors: command ordering (producer before consumer, nothing still reads
the pristine Doxyfile), `configure_file(... COPYONLY)` as the copy primitive and inert as a
reconfigure dependency under `cmake -P`, `(^|\n)TAG` anchoring not matching commented lines, the
Python-absent path, and both content fixes.

## 7. Audit Trail

- Exchange thread: `tmp/ai_exchange/make_doc_repair.md` (ephemeral; distil before sweep).
- **Plan audit round 2026-08-28 — Codex + Grok, both read-only, no worktree breach.** Every finding
  was re-verified against the cited files by Claude before acceptance.
  - **Codex** caught: the false `configure_file` premise (§2.1); that the warning script never fails
    the build (O5); that `LAYOUT_FILE` points into the source tree so R4 was inert (gap row 7); and
    that `make doc` writes tracked files (O6, gap row 11).
  - **Grok** caught: the fabricated Bezier identifiers in R5 — v1 said `dNydxN`/`dNxdyN`, which
    exist nowhere in the tree; they were an artifact of Claude's own `sed 's/[0-9]+/N/g'` when
    tabulating the log, and the real names are `d3ydx3`/`d3xdy3`. Also: that R3-before-R1 was not
    the B1 fix (§6 correction); that `Doxyfile.in:1433-1438` settles O4 against variant (a); that
    B4 is log noise because the custom tabs survive (`menudata.js:26-29`); that B4 over-cited
    `:77,:100`; that the patch script's own recovery command is B1 (gap row 3); and that log lines
    overcount distinct defects (§1.2).
  - **Both independently** chose variant (b) for O4 and flagged the `LAYOUT_FILE` retarget.
  - Claude's pre-audit measurement that `-u` does not validate `HTML_HEADER` resolved Grok's one
    medium-confidence open risk.
