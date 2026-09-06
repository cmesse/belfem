# Doxygen Site: Three Silent Defaults, Three-AI Sweep, Wave 1 Landed

**Date:** 2026-08-18
**Purpose:** Record the three-AI brainstorm sweep over the generated Doxygen documentation, the
root causes it converged on, and the configuration changes applied and verified in this session.
**Module:** meta (documentation toolchain) — `Doxyfile.in`, `cmake/report_doxygen_warnings.cmake`

---

## 1. Trigger

`make doc` produced content the user judged good but structurally hard to navigate. A three-AI
brainstorm sweep was run: Claude (independent pass + empirical verification), Codex and Grok
(blind, parallel, `tmp/ai_exchange/doxygen_doc_ux.md`).

## 2. What was actually wrong

The site was not merely unstructured — it was broken in three places, all caused by **Doxygen
defaults that changed after `Doxyfile.in` was written**. The file is a 2568-line template whose
header still reads `Doxyfile 1.8.16`; the tree now runs Doxygen 1.18.0. Three tags that did not
exist when the template was generated therefore took their modern defaults silently.

| Tag | Default in 1.18.0 | Effect on the BELFEM site |
|---|---|---|
| `IMPLICIT_DIR_DOCS` | `YES` | all 24 curated index READMEs consumed as *directory* documentation |
| `MARKDOWN_ID_STYLE` | `DOXYGEN` | every hand-written GitHub-style TOC in the usage guides dead |
| `HAVE_DOT` (hardcoded `YES`, no Graphviz) | — | 1437 graphs requested, none rendered |

Measured before the fix: **232 `warning:` + 1471 `error:`** lines in
`build/doc/doxygen_warnings.log`.

### 2.1 The severed module tree (`IMPLICIT_DIR_DOCS`)

Since Doxygen 1.9.8 a `README.md` inside a directory is bound to that *directory* rather than
rendered as a page, and its `{#anchor}` is discarded. `scripts/update_doc_index.py` builds the
BELFEM page tree precisely from those anchors, so the whole hand-written tier detached:

- `doc_index.html`, `homology_index.html`, `fem_maxwell_index.html` were never generated
- the `src/homology/doc/README.md` body was reachable only at
  `dir_ef0d80a70cdd471c990285920d5f72e0.html`, titled "Directory Reference"
- all 23 `@subpage <module>_index` lines in the generated `doc/doxygen_nav.dox` were dead
- `doc/mainpage.md:21` — "Find the module that does what you need" → `@ref doc_index` — was dead

The generator script is **not** at fault. Its `{#anchor}` + `.dox`-tree design is sound and works
for every non-README document; a Doxygen default was defeating it.

### 2.2 Dead in-page navigation (`MARKDOWN_ID_STYLE`)

Every usage guide opens with a GitHub-form table of contents
(`src/fem/kernel/doc/dof_manager_usage_guide.md:24-31`). Under the `DOXYGEN` default those slugs
are never minted — 222 unresolved `\ref` warnings, worst in `dof_manager_usage_guide.md` (24),
`graph_usage_guide.md` (22), `iwg_usage_guide.md` (20). The TOCs worked on the repository web view
and were dead on the generated site.

### 2.3 Missing diagrams (`HAVE_DOT`)

`Doxyfile.in` hardcoded `HAVE_DOT = YES` while Graphviz was not installed — CMake's
`find_package(Doxygen)` had reported "missing components: dot". Header pages rendered the caption
"Include dependency graph for cl_Mesh.hpp:" followed by the literal comment `<!-- SVG 0 -->`.
`build/doc/html` held 1437 `.md5` graph stamps and exactly one `.svg` (the Doxygen footer logo).
Christian installed Graphviz 15.1.1 and put it on `PATH` during the session; per his instruction
no discovery logic was added to the build, since `/opt/graphviz/latest/bin` is specific to that
machine.

### 2.4 The guard that could not see it

`cmake/report_doxygen_warnings.cmake` matched `warning:` only, so the 1471 `error:` lines — the
loudest failure in the log — never reached the build output. The script's own comment asserted a
zero-warning baseline while the tree carried 232.

## 3. Three-AI convergence

All three independently found the severed README tree, the dead TOCs, the flat `docmap`, the
audience mix on the landing page, and `examples/` being absent from the site. Divergences that
mattered:

- **Cost.** Codex and Grok both diagnosed the README breakage correctly and both proposed surgery —
  `EXCLUDE_PATTERNS` plus `@includedoc`, or rewriting `update_doc_index.py` — costed M. Neither
  named `IMPLICIT_DIR_DOCS`. The fix is one line. Taking either plan at face value would have meant
  substantial unnecessary script and markdown work.
- **Grok, refuted (2 items).** It listed `HAVE_DOT = YES` and the MathJax/Mermaid stack under
  "leave alone", having read only the first 200 of 1703 log lines and assumed the remainder were
  the same class — it flagged that limit itself. The remainder were the 1471 graph failures. Codex
  separately caught that MathJax 2.7.5 and Mermaid 11 are pulled from CDNs
  (`MATHJAX_RELPATH = https://cdnjs.cloudflare.com/...`, `cdn.jsdelivr.net/npm/mermaid@11`), so a
  local doc build silently requires network access.
- **Grok, uniquely right.** `JAVADOC_AUTOBRIEF = NO` was starving the class list: `cl_Cell.hpp:33-36`
  carries "Cell is a wrapper around the standard vector" while `annotated.html` rendered an empty
  description cell.

## 4. What landed (Wave 1)

| File | Change |
|---|---|
| `Doxyfile.in:359` | `MARKDOWN_ID_STYLE = GITHUB` (new) |
| `Doxyfile.in:937` | `IMPLICIT_DIR_DOCS = NO` (new) |
| `Doxyfile.in:204` | `JAVADOC_AUTOBRIEF` `NO` → `YES` |
| `cmake/report_doxygen_warnings.cmake` | counts `error:` as well as `warning:`, with separate baselines and a distinct message saying output is *missing*, not merely undocumented; warning baseline set to the new floor of 19 |

`HAVE_DOT` was deliberately left at `YES` and no `DOT_PATH` discovery was added — decided
2026-08-18, Christian: Graphviz lives at a machine-specific prefix and belongs on the developer's
`PATH`, not in the repository's build logic.

## 5. Verification

**Verified**, not merely reviewed: `doc/ai_collaboration_protocol.md` §11. Doxygen 1.18.0 was run
against the working tree from the patched `Doxyfile.in` with CMake-style `@VAR@` substitution,
output to a scratch directory so no build tree was touched.

| Metric | Before | After |
|---|---|---|
| `error:` lines | 1471 | **0** |
| `warning:` lines | 232 | **19** |
| unresolved `\ref` | 222 | 1 |
| graph SVGs rendered | 0 | **2687** |
| `*_index.html` module pages | 0 | **27** |
| class rows carrying a description | 26 / 318 | **85 / 318** |
| wall clock | — | 23 s |

Navigation was walked end to end: `index → docmap → mod_fem_kernel → fem_kernel_index` plus all six
of that module's documents now resolve. The DofManager guide's own TOC was checked programmatically —
24 in-page links, 0 dangling, against 24 dangling before.

`JAVADOC_AUTOBRIEF` was A/B tested on its own rather than bundled, because it changes brief
semantics tree-wide: 58 class rows gained a description, **0 rows lost text they previously had**,
and the warning count was unchanged at 19.

## 6. Residual — the 19 warnings

All genuine and small: 8 doc-comment defects in headers (duplicated `@param` at
`src/physics/materials/powerlaws.hpp:125,339`; `@return` documented on void functions at
`src/mesh/cl_Mesh.hpp:909,918` and `src/fem/iwg/cl_IWG.hpp:723,1170`, one of which carries the typo
`collect_nodes_on_wetted_sitdesets`; `@param` without arguments at `cl_FEM_Kernel.hpp:190` and
`cl_FEM_Group.hpp:304`), 6 bare `#autotoc_md` links in module READMEs, 2 relative-path `@ref`s in
`src/fem/interpolation/doc/nedelec_thinshell.md:61,75`, one unresolved
`known-issues-and-critical-bugs` anchor, and one real collision: `src/numerics/spline/doc/README.md:9`
claims the section label `index`, which fights the main page.

One cosmetic consequence of `JAVADOC_AUTOBRIEF`: `src/core/globals.hpp:26` attaches a maintainer
note to `namespace belfem`, so the namespace brief now reads "USER GUIDES:". Changing that block's
`/**` to `/*` fixes it; not applied, as it is a source edit outside the approved scope.

## 7. Wave 2 — navigation (same session)

Approved and landed after Wave 1 reproduced in the user's own build tree (19 warnings, 0 errors,
2687 graphs, 27 module pages, 85/318 class descriptions).

### 7.1 The pages tab: both earlier explanations were wrong

Grok reported that the top navigation bar has no Related Pages tab; Codex and Claude proposed a
`LAYOUT_FILE` tab *reorder*. Building it showed the premise was wrong — Doxygen's default layout
already lists `pages` second. The tab is suppressed because `scripts/update_doc_index.py` makes all
111 pages `@subpage`s of the main page, which leaves Doxygen's page index empty. A reorder would
have changed nothing.

Two variants were built and compared on full runs:

| Variant | Tab bar | Verdict |
|---|---|---|
| A — three `<tab type="user">` entries in `LAYOUT_FILE` | `Main Page \| Getting Started \| Modules \| Documentation \| Namespaces \| Classes \| Files` | **adopted** |
| B — `DISABLE_INDEX = YES`, `FULL_SIDEBAR = YES` | none; no `menudata.js` emitted | rejected |

Variant A wins on the case that matters: every page loads `menudata.js` + `menu.js`, so a developer
who lands on `classbelfem_1_1_cell.html` from search gets a one-click path back to the guides.
Variant B leaves that reader with the sidebar only.

### 7.2 What landed

| File | Change |
|---|---|
| `doc/DoxygenLayout.xml` | new, generated by `doxygen -l` (1.18.0) plus three user tabs |
| `Doxyfile.in:777` | `LAYOUT_FILE` set |
| `Doxyfile.in` `INPUT` | `examples/` added |
| `Doxyfile.in:1011` | `EXAMPLE_PATH = .../examples` |
| `Doxyfile.in:988` | `EXCLUDE_PATTERNS` for `*/main.cpp`, `*/*test.cpp`, `*/test_*.cpp`, `*/poisson.cpp` |
| `scripts/update_doc_index.py` | `MODULE_GROUPS` catalogue table; grouped `docmap`; human module titles; `collect_examples()`; `anchor_for()` extended for the `examples/` tree; coverage guard in `main()` |
| `doc/getting_started.md:70-73` | points at `@ref doc_examples` instead of a repository path |
| `examples/README.md`, `examples/scripts/README.md` | anchors injected by the script |

`docmap` now renders four `<h1 class="doxsection">` groups — Infrastructure, Mathematics, Finite
Elements, Physics and Circuits — each module a link with a one-line blurb, and module pages are
titled `Maxwell` rather than `src/fem/maxwell` (the path moved into the page body). `MODULE_GROUPS`
is the single source for group, title and blurb, and `main()` now exits non-zero if a documented
module has no entry, so adding a module cannot silently drop it off the catalogue.

### 7.3 Verification

Full run from the patched `Doxyfile.in`, CMake-style substitution, scratch output: **19 warnings,
0 errors** — unchanged from Wave 1, so none of this cost anything in log health.

- tab bar: `Main Page | Getting Started | Modules | Documentation | Namespaces | Classes | Files`
- all 12 in-tree drivers gone from the file index; every one had been confirmed to carry its own
  `main()` first. `hphirun`, `hphiTrun`, `electricalCircuit` remain, and `cl_IWG_Poisson.cpp/.hpp`
  — the real physics module, not the driver — survives the `*/poisson.cpp` pattern as intended
- `doc_examples.html` and `doc_examples_scripts.html` render; 11 scenarios and 13 `.conf` decks are
  now reachable from the site

## 8. R7 and R11 — landing-page grouping and the topics tree

**R7 (O2 decided 2026-08-18, Christian: keep the AI-process documents, group them).** The landing
page contents are now emitted in three reader paths — Run a simulation / Extend a module / Project
process — from a `PAGE_GROUPS` table in `scripts/update_doc_index.py`, with a guard that fails the
run if a page under `doc/` or `examples/` is unfiled. Previously the flat alphabetical list put
`ai_collaboration_protocol` and `ai_workflow_best_practices` second and third, ahead of Getting
Started.

**R11, first slice.** `doc/groups.dox` is new: 4 top-level `@defgroup`s and 23 module subgroups,
generated from the same `MODULE_GROUPS` table so the topics tree, the `docmap` catalogue and
`doc/README.md` cannot drift apart. Each group carries its blurb and an `@ref` to the module's
hand-written index, so a group renders usefully even before any class joins it.

Seven entry-point classes were tagged with `@brief` + `@ingroup` + `@see` — none of them had a
documentation block at all before:

| Class | Group | Guide |
|---|---|---|
| `MaxwellFactory`, `IWG_Maxwell`, `MaxwellPostprocessor`, `MaxwellBoundaryConditionFactory` | `grp_fem_maxwell` | Maxwell usage guide / postprocessor recovery theory |
| `Kernel`, `DofManager`, `Controller` | `grp_fem_kernel` | FEM kernel index / DofManager guide / nonlinear controller theory |

### Verification

Full run, scratch output: **19 warnings, 0 errors** — unchanged across R7 and R11, so the topics
tree and the seven new doc blocks introduced no new defects.

- tab bar: `Main Page | Getting Started | Modules | Documentation | Topics | Namespaces | Classes | Files`
- `topics.html` renders all 27 groups in the four-bucket hierarchy
- `group__grp__fem__maxwell.html` lists its four classes, each with the brief from the module
  README's key-class table, plus a collaboration diagram
- the `@see` on `classbelfem_1_1fem_1_1_maxwell_factory.html` resolves to
  "Maxwell Module Usage Guide" — the guide↔API seam both auditors flagged as missing
- class rows carrying a description: 85/318 → **91/318**

Populating the remaining 21 groups is incremental and deliberately not attempted in one pass.

## 9. R11, second slice — infrastructure modules

Seventeen more classes tagged across four modules, taking the topics tree from 2 populated groups
to 6:

| Group | Classes | Note |
|---|---|---|
| `grp_containers` | `Cell`, `Map`, `OrderedMap`, `Set`, `Queue`, `Bitset`, `DynamicBitset`, `ShiftRegister`, `StringList` | `Cell` and `StringList` already had documentation blocks; the tagger merged `@ingroup`/`@see` into them rather than adding a second block |
| `grp_sparse` | `SpMatrix`, `Solver`, `SolverParameters` | |
| `grp_mesh` | `Mesh` | |
| `grp_linalg` | `Vector`, `Matrix` | tagged in all four backend headers; see below |

Briefs were taken from each module README's quick-reference table, so the class list and the guide
now say the same thing.

### 9.1 The linalg backend trap

`src/linalg/cl_Vector.hpp` and `cl_Matrix.hpp` are dispatch headers — they `#include` the Armadillo
or Blaze implementation, and contain no class of their own. The real classes live in
`src/linalg/{armadillo,blaze}/cl_{AR,BZ}_{Vector,Matrix}.hpp`, and both define `belfem::Vector` /
`belfem::Matrix` in the same namespace. Doxygen parses both (no macro expansion) and **merges them
into one entity**.

The first pass gave each a backend-specific brief, and the merge picked Armadillo's — so the class
list read "Column vector (Armadillo backend)" on a machine whose default backend is Blaze. Caught
in verification and corrected: all four briefs are now backend-neutral ("Column vector"), with the
build-time backend choice stated in the detailed text. **Do not reintroduce a backend name in the
brief of a `cl_{AR,BZ}_*` class** — only one of the two can be right for any given build, and the
reader has no way to tell which they are looking at.

`src/linalg` has no other classes at all: everything else is free functions in `fn_*.hpp`. Grouping
those needs file-level `@ingroup` or `@addtogroup`, not the class mechanism used here.

### 9.2 Verification

Full run: **19 warnings, 0 errors** — unchanged across both R11 slices and the linalg correction.
Class rows carrying a description: 91/318 → **104/318** (26/318 before Wave 1). Topics tree:
**6 of 27 groups populated**, 24 classes.

## 10. R10 — MathJax 3 (done, not the way it was framed)

The question asked was whether MathJax 3 is supported on RHEL8. It is the wrong axis: MathJax is
browser-side JavaScript, downloaded by the reader, and never executes on the build host. What
matters is the Doxygen version on the builder — and Christian is building 1.18.0 there, which
settles it.

The naive change would have been unsafe. Verified by probe on 1.18.0:

| `MATHJAX_VERSION` | URL Doxygen emits |
|---|---|
| `MathJax_2` | `https://cdn.jsdelivr.net/npm/mathjax@2/MathJax.js` |
| `MathJax_3` | `https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml.js` |
| `MathJax_4` | `https://cdn.jsdelivr.net/npm/mathjax@4/tex-chtml.js` |

The path *shapes* differ, so `MATHJAX_VERSION` and `MATHJAX_RELPATH` are a matched pair. Pinning a
v3 `MATHJAX_RELPATH` and handing it to a Doxygen too old to understand `MATHJAX_VERSION` yields
`<v3-path>/MathJax.js` — a 404, and all 374 formulas blank with no build error.

**Landed:** `MATHJAX_VERSION = MathJax_3`, and `MATHJAX_RELPATH` returned to empty with a comment
saying it must stay that way. Doxygen then derives the matching CDN path itself, so any machine
still on the distro's 1.8.14 degrades to a working MathJax 2 rather than to 404s. This also removes
the previous pin to `cdnjs.cloudflare.com/.../mathjax/2.7.5/` — one end-of-life release on one CDN,
with that same silent-failure mode if the path is ever retired. `MATHJAX_FORMAT` stays `HTML-CSS`;
Doxygen translates it to `chtml` for v3 on its own.

Verified: 19 warnings, 0 errors (unchanged); the emitted URL is
`https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml.js`; no generated page references cdnjs any
more; a formula-bearing page (`cl_IWG_Timestep.hpp`) still loads it. An unsupported tag was also
confirmed to warn on **stderr only**, never into `WARN_LOGFILE`, so an older Doxygen cannot disturb
the 19-warning baseline.

### 10.1 Found while verifying: the log is not the whole story

Doxygen writes configuration complaints to stderr, not to `WARN_LOGFILE`. The current Doxyfile draws
**14 obsolete-tag warnings** that no one has been seeing:

`CLASS_DIAGRAMS`, `COLS_IN_ALPHA_INDEX`, `DOCBOOK_PROGRAMLISTING`, `DOT_FONTNAME`, `DOT_FONTSIZE`,
`DOT_MULTI_TARGETS`, `DOT_TRANSPARENT`, `FORMULA_TRANSPARENT`, `HTML_TIMESTAMP`,
`LATEX_SOURCE_CODE`, `LATEX_TIMESTAMP`, `OUTPUT_TEXT_DIRECTION`, `RTF_SOURCE_CODE`, `TCL_SUBST`

Two of them are actively configured and silently ignored: `DOT_FONTNAME = Helvetica` and
`DOT_FONTSIZE = 10`. (Corrected after the migration in §11: this had no visible effect. Doxygen's
replacement `DOT_COMMON_ATTR` defaults to exactly `fontname=Helvetica,fontsize=10`, so the graphs
were already rendering as the file intended. The tags were dead, not harmful.) Three more (`CLANG_ASSISTED_PARSING`, `CLANG_OPTIONS`, `CLANG_DATABASE_PATH`)
name an option this Doxygen was not compiled with. `doxygen -u Doxyfile.in` migrates the template
and is the fix; it rewrites comments wholesale, so the diff needs reading rather than trusting.

## 11. R13 — Doxyfile template migration

`doxygen -u` run on the 1.8.16-era template with Doxygen 1.18.0. Note the mechanics: `-u` must be
invoked from the file's own directory (an absolute path silently left the file untouched), and it
writes a `.bak` beside the original.

The migration was reviewed by parsing every `KEY = value` pair before and after — comments ignored,
multi-line continuations joined — rather than by reading a 3100-line textual diff.

| | count | assessment |
|---|---|---|
| keys before / after | 278 / 311 | |
| removed | 17 | all obsolete or naming an option this build lacks |
| **values changed** | **1** | `EXTRA_PACKAGES` list separator normalised; `GENERATE_LATEX = NO`, so inert |
| added | 50 | every one at Doxygen's own default, which was already in force while the tag was absent |

That last row is why the migration is behaviour-neutral by construction: an absent tag and a tag set
to its default are the same thing to Doxygen.

All BELFEM settings survived intact (`IMPLICIT_DIR_DOCS`, `MARKDOWN_ID_STYLE`, `JAVADOC_AUTOBRIEF`,
`LAYOUT_FILE`, `EXAMPLE_PATH`, `EXCLUDE_PATTERNS`, `MATHJAX_VERSION`, `MATHJAX_RELPATH`, `HAVE_DOT`,
`EXTRACT_ALL`), as did all 12 `@VAR@` placeholders. The BELFEM rationale comments did **not** —
`-u` rewrites comment blocks wholesale — so all eight were re-injected above their keys, marked
`# BELFEM:` so a future migration can find them again.

`DOT_FONTNAME` / `DOT_FONTSIZE` were folded into `DOT_COMMON_ATTR = "fontname=Helvetica,fontsize=10"`
and `DOT_EDGE_ATTR`, preserving the intent the old tags expressed.

### 11.1 Verification

- stderr: **14 obsolete-tag warnings and 3 not-compiled-in warnings → 0**
- warning log: 19 warnings, 0 errors, and **byte-identical** to the pre-migration run
- generated tree: 2879 HTML + 2725 SVG, and a full recursive diff against the pre-migration output
  reports **0 differing files**

Byte-identical output across the whole tree is the strongest available evidence that nothing
changed but the configuration file's vocabulary.

One stderr line remains and is unrelated: a Fortran parser complaint about
`src/sparse/splinalg.f90:116`. It is present in the pre-migration runs too.

### 11.2 Newly visible knobs

The migration surfaced tags the old template could not express. None were changed; recording them so
the choice is deliberate rather than forgotten:

- `HTML_COLORSTYLE = AUTO_LIGHT` — `TOGGLE` would give the site a light/dark switch
- `MERMAID_JS_URL = https://cdn.jsdelivr.net/npm/mermaid@11/dist` — the second CDN dependency Codex
  found, now settable; nothing in the tree uses mermaid today
- `PAGE_OUTLINE_PANEL`, `HTML_CODE_FOLDING`, `HTML_COPY_CLIPBOARD` — all default YES already
- `NUM_PROC_THREADS = 1` — 0 would use all cores; `make doc` currently takes ~25-50 s

## 12. R12 — the log reaches zero

**0 warnings, 0 errors.** From 232 warnings + 1471 errors at the start of the session.

### 12.1 The nine documentation warnings, and what they actually were

Only one was what it looked like.

| Reported as | Real cause | Fix |
|---|---|---|
| `multiple use of section label 'index'` | `## Index` in `src/numerics/spline/doc/README.md` mints the GitHub slug `index` under `MARKDOWN_ID_STYLE = GITHUB`, colliding with the main page | heading renamed to `## Documents` |
| 5x `unable to resolve reference to 'autotoc_md'` | not an anchor problem at all: `- [Source code](../)` in five module READMEs — a markdown link to a bare *directory*, which Doxygen cannot resolve | converted to unlinked text naming the path |
| 2x `unable to resolve reference to '../kernel/cl_ThinShellFactory.cpp'` | relative link out of the doc tree to a source file | replaced with the full path in a code span |
| `unable to resolve reference to 'known-issues-and-critical-bugs'` | genuine drift: `iwg_usage_guide.md:35` linked "Known Issues and Critical Bugs" but the heading at :1949 reads "Known Issues and Historical Bug Record" | TOC entry corrected to match |

The `[Source code](../)` case is worth remembering: the warning text names `autotoc_md`, which
points at anchors and sends you looking in entirely the wrong place. The trade-off taken is that
those five READMEs lose a working directory link on the GitHub web view; the generated site has the
Files index for that, and Doxygen has no way to link a directory.

### 12.2 The ten header warnings

- `cl_Mesh.hpp` (2), `cl_IWG.hpp` (2): `@return` on `void` functions — removed. Two mangled words
  in the same comments were fixed while there: "bedimfore" -> "before", "sidesteds" -> "sidesets".
- `cl_FEM_Group.hpp`, `cl_FEM_Kernel.hpp`: `@param` on accessors that take no argument — converted
  to the `@return` they were evidently meant to be.
- `powerlaws.hpp` (4): a red herring. The duplication is not inside that file — the *declarations*
  in `cl_Material.hpp:625,676` repeat the `@param normJ` / `@return` that the definitions already
  carry, and Doxygen merges declaration and definition documentation. The definitions own the
  derivation, the units and the literature references, so the declaration's copies were removed and
  a pointer left in their place.

**Not fixed, flagged instead:** `IWG::collect_nodes_on_wetted_sitdesets` (`cl_IWG.hpp`) has the typo
in the *function name*. Renaming it is a source change with call sites, not a documentation fix, and
was left alone deliberately.

### 12.3 Three modules that had no documentation at all

`src/fem/thermal`, `src/fem/postproc` and `src/visualizer` had no `doc/` directory and appeared
nowhere in the navigation. Each now has a README written from its headers, its CMake guard and its
in-tree consumers, and each closes with an explicit statement of what is *not* yet documented rather
than implying completeness:

- **Thermal** — `ThermalFactory`, `IWG_MaxwellThermal` (deriving from
  `IWG_TransientHeatConduction`), `ThermalBoundaryConditionFactory`, and the `matrices/` element
  kernels for the h- and phi-regions. The two `ThermalFactory` constructors (standalone on a mesh,
  or coupled to an existing magnetic `Kernel`) are the module's whole design. Sole consumer:
  `src/executables/hphiTrun.cpp`. The coupling mechanism itself is *not* documented and needs a
  usage guide.
- **Postprocessing** — `Gradient`, `Surface`, `MeshChecker` plus six `fn_Mesh_*` free functions.
  The README leads with the distinction that matters: this is physics-agnostic and is **not**
  `MaxwellPostprocessor`, which lives in `src/fem/maxwell`.
- **Visualizer** — `VTK_MeshView`, `VTK_Curve`, the `visualize` executable. The README states
  plainly that `USE_VTK` defaults OFF, that nothing else depends on it, and that it is therefore
  the least exercised part of the tree.

The R6 coverage guard did its job: adding the three directories made
`scripts/update_doc_index.py` exit non-zero naming all three as unclassified, before any Doxygen run.

### 12.4 groups.dox is now generated

`doc/groups.dox` was hand-generated in the first R11 slice, which meant the topic tree could drift
from `MODULE_GROUPS`. It is now emitted by `scripts/update_doc_index.py` from that same table, like
`doxygen_nav.dox` — so the catalogue page, the topic tree and `doc/README.md` have one source. The
tree carries 26 modules and 30 groups.

### 12.5 Baseline locked

`cmake/report_doxygen_warnings.cmake`: both baselines set to **0**, with a note that Doxygen writes
configuration complaints to stderr rather than to the log, so a clean log alone is not proof of a
clean run.

## 13. R11a — Finite Elements

R11 split into four sub-tasks by catalogue bucket (R11a-d, see the todo). R11a covers the Finite
Elements group: 15 classes tagged across four modules.

| Module | Classes |
|---|---|
| `fem/iwg` | `IWG`, `IWG_Timestep`, `TimestepMatrices`, `IwgFactory` |
| `fem/interpolation` | `InterpolationFunction`, `InterpolationFunctionFactory`, `IntegrationData`, `EdgeFunction`, `EdgeFunctionFactory` |
| `fem/thermal` | `ThermalFactory`, `IWG_MaxwellThermal`, `ThermalBoundaryConditionFactory` |
| `fem/postproc` | `Gradient`, `Surface`, `MeshChecker` |

Four already carried documentation blocks (`IWG`, `TimestepMatrices`, `InterpolationFunction`,
`EdgeFunction`) and the tagger merged into them rather than adding a second block. The Nedelec
classes point at `@ref fem_interpolation_nedelec` rather than the general interpolation guide, since
that is where their derivation lives.

`grp_fem` was deliberately left empty: `src/fem` has no top-level headers, so the group is an
overview only. That is recorded in the plan so nobody later "fixes" it.

**Verified:** 0 warnings, 0 errors. 10 of 26 module groups populated, 39 classes. Class rows
carrying a description: 104/318 -> **115/318**.

## 14. R11b — Mathematics

Twelve classes across five modules. The homology module carries the most theory documentation in
the tree, so its classes were routed individually rather than all at the module index:

| Class | Points at |
|---|---|
| `CutFactory`, `Homology`, `CutData` | homology usage guide |
| `Cohomology`, `Chain`, `Cochain` | cohomology theory and implementation |
| `SimplicialComplex` | cohomology algorithms |
| `CutProcessor` | thick/thin cuts and conjugate edges |

Plus `Vertex` (`math/graph`), `Tensor` (`math/tensor`), `Quaternion` (`math/quaternion`) and
`Spline` (`numerics/spline`), each to its module's usage guide.

**Verified:** 0 warnings, 0 errors. 15 of 26 module groups populated, 51 classes. Class rows
carrying a description: 115/318 -> **126/318**.

## 15. Prose review policy (2026-08-18, Christian)

Standing instruction from this session: **whenever prose is added or revised, run a Codex language
sweep for readability before treating it as finished.** This generalises the narrower rule in
`CLAUDE.md`, which asks for a Codex prose pass only over `doc/input_file_reference.md`.

The sweep must be told what may not change: technical claims, file paths, code identifiers, CMake
option names, `file:line` citations, the BELFEM header block, and — for the three module READMEs
written in R12 — the closing **Status:** paragraph, which deliberately records what is *not* yet
documented. A smoother sentence that alters a technical claim is a regression, not an improvement.

Applies to reader-facing documentation. Devlog entries and `tmp/ai_exchange/` scratch are records,
not documentation, and do not need it.

### 15.1 First sweep, applied

Run over the three READMEs written in R12. Codex returned 11 edits and reported no
technical-claim problems. Nine were taken verbatim; two were adjusted after reading them against
the files:

- **Thermal, capability list.** Codex asked for a verb-leading first bullet. The other three
  bullets were noun phrases, so taking that edit alone would have broken the parallelism — all four
  were converted instead.
- **Postproc, the `MaxwellPostprocessor` warning.** Codex proposed "`MaxwellPostprocessor` is
  separate." That is plainer but drops the warning function: the sentence exists to head off a
  specific mix-up, and "is separate" does not signal one. Kept as "**This is not
  `MaxwellPostprocessor`.**" — direct, and still a signpost.

The rest were genuine improvements: long sentences split at the em dash, passive constructions made
active ("`Doxyfile.in` excludes them" for "These are excluded"), and "Enable it with `-DUSE_VTK=ON`"
in place of a description of when the module happens to get built. The closing **Status:**
paragraphs came through untouched, which was the main thing to check.

## 16. R14 — the project logo links to belfem.lbl.gov

Feature request: clicking `logo.gif` should go to the project site.

Doxygen has no option for this. `PROJECT_LOGO` is emitted as a bare `<img>` inside
`<td id="projectlogo">`, and the only place to wrap it is a custom `HTML_HEADER`.

### 16.1 Why the header is generated, not checked in

The obvious implementation — `doxygen -w html`, edit one line, commit the file, point
`HTML_HEADER` at it — was rejected after testing the failure mode. **Doxygen does not warn when a
custom header is stale.** Verified on 1.18.0: a header with `$treeview` and `$mathjax` stripped out
built with exit 0, an empty stderr and a clean warning log, and the generated site simply had no
navigation tree and no MathJax.

A header committed against one Doxygen version would therefore lose the navigation tree, the search
box, dark mode or MathJax on a future upgrade, silently — the same class of failure as the three
defaults that started this session, and this time with no warning channel at all to catch it.

So the `doc` target regenerates the header from the **installed** Doxygen on every run and patches
it:

```cmake
COMMAND ${DOXYGEN_EXECUTABLE} -w html ${BELFEM_DOXYGEN_HEADER} ... 
COMMAND ${CMAKE_COMMAND} -DHEADER=... -DURL=... -P cmake/patch_doxygen_header.cmake
COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
```

`cmake/patch_doxygen_header.cmake` wraps the logo in an anchor. It **aborts with FATAL_ERROR** if
the expected markup is not found, rather than writing an unpatched header — an unlinked logo is
exactly the kind of silent regression that would go unnoticed for a year. The error names the
command that re-derives the current markup.

### 16.2 Verification

- logo renders as
  `<a href="https://belfem.lbl.gov" title="BELFEM project site"><img alt="Logo" src="logo.gif"/></a>`
- present on **2893 of 2894** pages, including deep class pages. The one exception is
  `doxygen_crawl.html`, a validator/crawler helper with no title area — correct.
- the custom header dropped nothing: navigation tree, search box and MathJax 3 all still emitted
- 0 warnings, 0 errors
- against the previous build, exactly 2893 files differ — the pages that gained the link, and
  nothing else
- the FATAL_ERROR path was exercised by mangling the markup, and fires as intended

The URL was taken as given from the request; it was not resolved from here.

## 17. R11c — Physics and Circuits

23 classes across five modules: `physics/materials` (4), `physics/database` (2),
`physics/gasmodels` (7), `physics/gastables` (5), `circuit` (5).

Worth recording: **18 of the 23 merged into documentation blocks that already existed.** The physics
modules are the best-documented part of the tree, so the tagger added only `@ingroup` and `@see`
there and left the existing briefs alone — the briefs drafted for those classes were discarded
rather than competing with what the authors had already written. Only `Database`, `Projector` and
the three structural circuit classes needed a block of their own.

Two locator mistakes were caught before tagging: a plain grep for `class Gas` and `class RefGas`
matched **forward declarations** in `cl_GM_EoS.hpp` and `cl_GT_GasData.hpp` rather than the
definitions in `cl_Gas.hpp` and `cl_GT_RefGas.hpp`. The tagger itself rejects forward declarations,
so this would have surfaced as "declaration not found" rather than a wrong tag — but the module
READMEs' own file references are the reliable source and were used instead.

`grp_executables` is deliberately empty, for the same reason as `grp_fem`: `src/executables` holds
three `.cpp` files with `main()` and no classes.

**Verified:** 0 warnings, 0 errors. 20 of 26 module groups populated, 76 classes. Class rows
carrying a description: 126/318 -> **133/318**.

## 18. R11d — infrastructure, and R11 complete

15 classes: `core` (6), `io` (6), `comm` (1), `visualizer` (2).

An error introduced in R12 was caught here: the visualizer README listed **file** names as class
names. The classes are `belfem::vtk::MeshView` and `belfem::vtk::Curve` -- the `VTK_` prefix belongs
to the filenames -- and `BlockActor` / `SideSetActor` were missing from the table entirely. Both
fixed, with a line stating where the prefix actually lives.

**R11 is complete: 24 of 26 module groups populated, 92 classes.** Class rows carrying a
description: **26/318 at the start of the session -> 146/318.**

The two empty groups are correct and are recorded as such in the plan: `src/fem` and
`src/executables` hold no classes -- the first is a parent directory, the second is three `.cpp`
files with `main()`.

## 19. R15 — documenting the linalg public API

`src/linalg` is the largest undocumented public surface in the tree: **zero** `@brief` across all
44 files, roughly **97 declarations**.

### 19.1 The design decision, and a correction

Christian's instinct was to put the documentation in the top-level `src/linalg/fn_*.hpp` files
rather than in the per-backend implementations, on the grounds that a header containing only
`#ifdef`-selected includes is still conceptually "the dot product header".

I had argued that this could only carry a file-level brief, because documenting the function itself
from there would need `@fn`, which I expected to be brittle against the backend return types
(`-> decltype( arma::dot( ... ) )`). **That was wrong.** Tested on 1.18.0: an `@fn` block written
*without* the trailing return type binds cleanly to the backend declaration, with no warnings, and
still binds when both backend headers are parsed in the same run. So the full documentation --
brief, parameters, return -- can live at the public surface, exactly as proposed, and does not have
to be written twice.

Known artifact, accepted rather than fought: the rendered signature keeps whichever backend Doxygen
parsed first, so a Blaze build still displays `-> decltype( arma::dot( ... ) )`. The remedy is to
state the result in words in `@return`, not to try to correct the signature.

### 19.2 R15a — dot, cross, crossmat, norm, trans

The 97 declarations are far fewer distinct *ideas* than that count suggests. The overloads fall into
two kinds:

- **concrete** overloads on `Vector<T>` / `Matrix<T>` -- what a caller writes, and worth an `@fn`
  each;
- **expression-typed** overloads on a generic `ET`, which exist so an unevaluated backend expression
  (`a + b`) can be passed without first being materialised into a temporary.

Documenting every `ET` overload separately would add noise, not information. They are instead
explained once in the file-level block, and only the concrete forms carry an `@fn`. `dot` therefore
needs three blocks rather than 24.

`fn_crossmat.hpp` and `fn_trans.hpp` turned out to be direct implementations rather than dispatch
headers, so their four and one overloads are documented in place.

Reading `crossmat`'s body rather than inferring from its name surfaced behaviour worth recording:
every overload finishes by zeroing result entries whose magnitude is below `BELFEM_EPSILON`
relative to the norm of the result. That is cosmetic dust removal, and a caller relying on the raw
arithmetic would not expect it. It is now documented.

**Verified:** 0 warnings, 0 errors. The five headers appear under the Linear Algebra topic with
their briefs, and `dot` (3 overloads), `cross` and `norm` carry the `@fn` briefs on the namespace
page -- confirming the mechanism binds, not merely that the files are grouped.

### 19.3 The prose sweep caught two wrong technical claims

The Codex readability pass over R15a returned 11 wording edits and, separately, **two flags on
technical claims** -- the "flag, do not fix" instruction earning itself immediately. Both were
verified against the code, and both were errors in documentation written the same session:

**1. `fn_crossmat.hpp` -- "All four overloads finish by zeroing..."** False. Only the two 2D
overloads do the epsilon dust removal; the 3D pair does not.

**2. `fn_dot.hpp` -- `dot( Matrix, Vector )` documented as a matrix-vector product returning a
vector.** False, and the more serious of the two. It forwards to `arma::dot( matrix_data,
vector_data )`, which walks both operands as flat element sequences and returns a **scalar**. The
documentation would have told a reader to expect a vector from a function that returns a number.

Checking those two flags surfaced a third error nobody had flagged: the overloads taking `aScale`
**accumulate** into `aNxA` with `+=`, while the unscaled ones **assign**. A caller following the
original wording ("result, one scalar per column") would have passed an uninitialised accumulator.

And a fourth, found while correcting Codex's own suggested wording: it proposed "three-row matrix",
but there is no assert on `aA.n_rows()` at all. The 2D form reads rows 0-1 and the 3D form rows
0-2; extra rows are ignored. The final wording says exactly that.

All four are now documented correctly. The lesson is not that Codex is good at maths -- it declined
to judge the arithmetic and said so -- but that asking a second reader to *flag rather than fix*
turns a readability pass into a cheap correctness net over prose written from a quick reading of
unfamiliar code.

**Verified after the corrections:** 0 warnings, 0 errors, and the corrected briefs render on the
namespace page.

## 20. R15b — the solvers

### 20.1 A measurement that was wrong

R15 was opened on the claim that `src/linalg` carries "zero `@brief` across all 44 files". That
count covered `src/linalg/*.hpp` and the two backend directories -- it did **not** include
`src/linalg/lapack/`, which had not been looked at. Corrected coverage:

| directory | files with a brief |
|---|---|
| `src/linalg` | 5 of 28 |
| `src/linalg/armadillo` | 2 of 15 |
| `src/linalg/blaze` | 2 of 15 |
| **`src/linalg/lapack`** | **9 of 10** |

The LAPACK layer is the best-documented corner of the module, not the worst. `fn_gesv.hpp` already
carries `@brief`, `@param[in,out]` on every argument, `@return`, and an explanation of when to pass
`AbortOnError = false` to recover from a singular matrix inside an iterative scheme. What it lacked
was any `@ingroup`, so none of it reached the topic tree.

R15b therefore split in two: **wire up** ten already-documented LAPACK headers, and **write**
documentation for the five that had none.

### 20.2 What reading the bodies turned up

`eigen( const Matrix<real> &, Vector<real> & )` computes general, non-symmetric eigenvalues through
`arma::eig_gen`, which returns complex values -- but the BELFEM signature returns a **real** vector.
Any eigenvalue whose imaginary part exceeds `BELFEM_EPSILON` is written as **`BELFEM_QUIET_NAN`**,
silently: no error, no warning, no return code. A caller who cannot rule out complex eigenvalues has
no way to learn this except by reading the implementation. It is now a `@warning` on the function.

`inv2` / `inv3` compute the closed-form adjugate and **return the determinant**, which is deliberate
-- element Jacobians need it anyway. Their singularity check is **relative**: `det^2` against
`eps^2` times the product of the squared row norms, so the tolerance scales with the magnitude of
the matrix. It is a `BELFEM_ASSERT`, so it is compiled out of release builds; a singular matrix in
release yields infinities rather than a diagnostic. Both facts are now documented.

### 20.3 A false alarm worth recording

While reading `fn_inv2.hpp` the singularity assert appeared to be missing the comma before its
message string -- which would not compile in a debug build. It was an artefact of the filter used
to strip comments: the assert continues across lines that **begin with `*`**, the multiplication
operator, and the filter treated those as comment continuations and dropped them. The code is
correct. Verified before reporting anything.

### 20.4 Verification

0 warnings, 0 errors. The Linear Algebra topic now lists **20 headers** with descriptions, up from
5 after R15a: the eight top-level function headers, `fn_inv2` / `fn_inv3`, and all ten of the LAPACK
layer. The `eigen` NaN warning renders on the group page.

Codex readability sweep dispatched; R15b stays open until it is applied.

### 20.5 The R15b sweep — four more technical flags

The second solver-family sweep returned readability edits plus **four technical-claim flags**, all
verified against the code and all correct. Two sweeps have now produced six caught errors; the
"flag, do not fix" instruction is doing more work than the copyediting.

1. **`fn_eigen.hpp` -- the tolerance is not one number.** The documentation said the cutoff is
   `BELFEM_EPSILON`. That is true only of the Armadillo path
   (`fn_AR_eigen.hpp:39`); the Blaze path compares against a **hard-coded `1e-15`**
   (`fn_BZ_eigen.hpp:41`), while `BELFEM_EPSILON` is `10 * DBL_EPSILON`, about `2.22e-15`. The two
   backends therefore apply different thresholds to the same decision. Documented as it is, and
   logged as a code inconsistency to resolve separately.
2. **`fn_gels` is not only for overdetermined systems.** The wrapper handles `m < n` too, through an
   LQ factorization and a minimum-norm solution -- its own parameter documentation says so
   (`fn_gels.hpp:266-270`). The brief now says "least-squares or minimum-norm".
3. **`fn_posv` is symmetric *or* Hermitian.** The real routine is symmetric positive-definite, the
   complex one Hermitian positive-definite, and the pre-existing function documentation already
   said so at `fn_posv.hpp:214`. The new file brief had dropped the complex half.
4. **"yields infinities" was too specific.** A singular matrix in a release build can equally
   produce NaN (zero determinant with zero cofactors) or a very large finite value (nearly
   singular). Reworded.

Verified after the corrections: 0 warnings, 0 errors, corrected briefs rendering.

### 20.6 Open question raised by this work

Christian's proposal, 2026-08-18: since a caller passing a real vector to `eigen` can only ever get
a NaN back, should that be `BELFEM_SIGNALING_NAN` rather than `BELFEM_QUIET_NAN`? Many engineering
problems have no reason to expect complex eigenvalues.

The goal is right; the mechanism would not deliver it. **Nothing in BELFEM unmasks floating-point
exceptions** -- no `feenableexcept`, no `<fenv.h>` -- so a signalling NaN never traps; it degrades
quietly to a quiet NaN on first use. It would rename the problem while looking like a safety net.
BELFEM's own convention supports that reading: two of the three existing `BELFEM_SIGNALING_NAN`
sites are the unreachable `return` *after* a `BELFEM_ERROR` has already fired.

The recommendation instead is the pattern already used two directories away by `gesv` and `posv`:

```cpp
int_t eigen( const Matrix<real> &, Vector<real> &, const bool aAbortOnComplex = true );
```

Abort by default with a `BELFEM_ERROR` naming the offending index and its imaginary part -- always
active, unlike an assert -- with the escape hatch for a caller that legitimately probes a
non-symmetric operator. Two facts make it cheap: `belfem::eigen()` has **no in-tree caller** (the
eigenvalue work goes through ARPACK in `src/sparse`), and if the intended use is symmetric matrices
then the deeper fix is the backend call itself -- `eig_gen` is the general solver, and a symmetric
path returns real eigenvalues by construction, so the case could not arise at all.

Not implemented: this is an API behaviour change, outside a documentation pass, and awaiting
Christian's decision.

## 21. R16a — the eigen tolerance unified

Decided 2026-08-18, Christian: Blaze should use `BELFEM_EPSILON` as well.

`src/linalg/blaze/fn_BZ_eigen.hpp:41` compared the imaginary part against a hard-coded `1e-15`
while `src/linalg/armadillo/fn_AR_eigen.hpp:39` used `BELFEM_EPSILON` (`10 * DBL_EPSILON`, about
`2.22e-15`). The same decision -- is this eigenvalue complex? -- therefore used two different
thresholds depending on which backend the tree was built against.

Changed to `BELFEM_EPSILON`. `typedefs.hpp` was already included, so no new dependency. It was the
only hard-coded epsilon literal anywhere under `src/linalg`.

The `@note` added in R15b, which documented the discrepancy, was rewritten in the same edit -- it
would otherwise have become a description of a state that no longer exists. Documentation that
outlives the behaviour it describes is how this whole session's problems started.

**Two caveats on what "verified" means here.** The Doxygen run is clean (0 warnings, 0 errors) and
the corrected note renders, but that only exercises the documentation. On the compilation side:

- this machine builds with `USE_MATRIX_BLAZE=ON` (`build/CMakeCache.txt:321`), so the edited file
  *is* the active backend path here rather than the dormant one;
- but **nothing in the tree includes `fn_eigen.hpp`**, so a build will not compile it either way.
  The change is a one-token substitution of an already-included macro, which is about as safe as an
  edit gets, but it has not been through a compiler.

## 22. R16b — eigen gets a policy flag, and a symmetric sibling

Decided 2026-08-18, Christian: take the `aAbortOnComplex` design, and add a symmetric wrapper if
both backends offer a symmetric solver. They do -- Armadillo has `arma::eig_sym`, Blaze ships
`blaze::syev` / `heev` (`/opt/scls/include/blaze/math/lapack/`) -- so both landed.

```cpp
int_t eigen    ( const Matrix<real> &, Vector<real> &, const bool aAbortOnComplex = true );
void  eigen_sym( const Matrix<real> &, Vector<real> & );
```

`eigen` now aborts on the first complex eigenvalue with a `BELFEM_ERROR` naming its index and its
imaginary part -- always active, unlike an assert -- and names both escape routes in the message
(`eigen_sym`, or `aAbortOnComplex = false`). With the flag false it fills NaN as before and returns
the count. This is deliberately the same shape as the `AbortOnError` argument on
`belfem::gesv` / `belfem::posv` two directories away.

`eigen_sym` exists because for a symmetric matrix the whole question is void: the eigenvalues are
real by construction, so there is no NaN to guard and no policy to choose. It is the right call for
most of what a finite-element code produces.

### 22.1 Three defects fixed on the way

**Both `eigen` definitions were non-inline free functions in a header.** `void eigen( ... )` with no
`inline`, in a header -- two translation units including it would collide at link time. It has never
bitten because **nothing includes `fn_eigen.hpp`**. Both are now `inline`.

**The two backends read different triangles.** `arma::eig_sym` uses `uplo = 'U'`
(`auxlib_meat.hpp`), and the first draft of the Blaze wrapper passed `'L'` to `syev`. For a genuinely
symmetric matrix that is invisible; for a caller who passes a non-symmetric matrix by mistake the
two backends would have returned **different eigenvalues**. Changed to `'U'`, with a comment saying
why, and the documentation states that the upper triangle is the one read. This is the same class of
divergence as the epsilon mismatch fixed in R16a, caught before it existed rather than after.

**`blaze::syev` overwrites its input**, unlike `arma::eig_sym` which takes a const reference. The
Blaze path takes a working copy so both backends present the same const interface.

### 22.2 Verification, and its limits

Doxygen: 0 warnings, 0 errors. Both functions render with their documentation on the namespace page
and under the Linear Algebra topic, and the signature shows `aAbortOnComplex=true`.

**This code has not been compiled.** Nothing in the tree includes `fn_eigen.hpp`, so `make` will not
reach it, and the Doxygen parse is not a substitute for a compiler -- it is more permissive by
design. This is new logic, not a token substitution: `int_t` return values, a `BELFEM_ERROR` with a
five-argument format string, a `Matrix` copy constructed from `matrix_data()`, and two backend calls
whose signatures were read from the installed headers rather than exercised. If this API is meant to
be live, the cheapest insurance is a unit test that includes the header -- which would also give the
file its first compiler pass in an unknown length of time.

## 23. R15c — the utilities

Measured before assuming, per the note left after R15b: 16 files still ungrouped and undocumented,
two of which -- `fn_to_cell.hpp` and `fn_to_vector.hpp` -- were not on the R15c list at all. They
are now.

All 16 have a `@file` block; `polyval`, `dpolyval`, `sort`, `unique` and `append` also have
function-level documentation, chosen because their contracts are the ones a caller can get wrong
silently.

### 23.1 The convention that matters most

The four polynomial headers share a paragraph stating that **coefficients are in descending power
order** -- `aCoeffs(0)` is the leading coefficient, the last entry is the constant term. Read out of
the Horner loops in `fn_polyval.hpp` and `fn_dpolyval.hpp` rather than assumed: `polyval` seeds the
accumulator with `aCoeffs(0)` and multiplies by `aX` on each step, and `dpolyval` multiplies
`aCoeffs(0)` by `length - 1`. It matches Armadillo and MATLAB, and is the opposite of the ascending
order several other libraries use. Getting it backwards produces a plausible wrong answer rather
than an error, which is exactly why it is worth stating four times.

`dpolyval` also takes the coefficients of the **polynomial**, not of its derivative, and
differentiates analytically -- nothing is a finite difference. Both are now said explicitly.

### 23.2 Contracts that were invisible

- **`unique` sorts as a side effect**, on both backends -- Armadillo through `arma::unique`, Blaze
  through an explicit `std::sort` before `std::unique`. It is not an order-preserving duplicate
  filter, and a caller expecting one would be quietly wrong.
- **`append` requires `aA` and `aB` to be different objects.** The Blaze path asserts it; the
  Armadillo path uses `arma::join_cols` and would silently produce a doubled vector. Documented as
  a `@warning` against relying on either. `aB` is taken by non-const reference but is not modified.
- **`sort` is in place and ascending**, which the name implies but the signature does not.

### 23.3 Verification

0 warnings, 0 errors. The Linear Algebra topic now lists **36 headers** with descriptions -- every
`fn_*.hpp` at the top level plus the ten LAPACK wrappers -- against 5 after R15a. The coefficient
paragraph renders on all four polynomial file pages.

One tidy-up: adding a block above `dpolyval` left it with two stacked comment blocks, the older one
a one-line stub. Removed, and the rest of the group checked for the same pattern with a regex that
distinguishes a genuine stacked pair from the licence banner sitting above a `@file` block.

Codex readability sweep dispatched; R15c stays open until it is applied.

## 24. R17 — and a correction: `fn_eigen.hpp` was never untested

Twice in this session I wrote that **nothing includes `fn_eigen.hpp`**, that it therefore never
reaches a compiler, and that R16b's new code had no build coverage. That was wrong, and the error
was in the search, not the reasoning: I grepped `src/` and never looked at `tests/`.

`tests/linalg/test_LinalgSolvers.cpp:30` includes `fn_eigen.hpp`, the file is listed in
`tests/linalg/CMakeLists.txt:9`, and the suite already contained five `eigen` tests. So the header
does get a compiler pass, and R16b's claim of "no in-tree caller, so the signature is free to
change" was built on the same bad search.

### 24.1 What that meant for R16b

`TEST( Eigen, EigenComplexEigenvaluesReturnNaN )` fed a rotation matrix -- eigenvalues
`cos θ ± i sin θ` -- to `eigen` with the default arguments and asserted that both entries came back
NaN. That is precisely the behaviour R16b changed. **The change broke an existing, deliberately
written test**, and only the absence of a build run hid it.

The test encoded the old contract, so it was updated rather than the code reverted -- the new
contract is the approved one:

| test | purpose |
|---|---|
| `EigenComplexEigenvaluesAbortByDefault` | rewritten from the old NaN test; now `EXPECT_THROW` |
| `EigenComplexEigenvaluesReturnNaNWhenNotAborting` | same matrix with `false`; checks NaN **and** that the return counts 2 |
| `EigenRealEigenvaluesReturnZeroCount` | a real-eigenvalue matrix reports 0 |
| `EigenSymReturnsAscendingRealValues` | `(5±√5)/2`, checked in ascending order |
| `EigenSymDiagonal` | ascending order from an unordered diagonal |
| `EigenSymAgreesWithEigenOnASymmetricMatrix` | the two entry points must not disagree where both are valid |

`BELFEM_ERROR` routes through the same `belfem::assert::belfem_assert` as `BELFEM_ASSERT` and throws
`std::runtime_error`, which the file was already relying on for its non-square test, so the new
abort is testable in the existing style.

The three surviving default-call tests all use real-eigenvalue matrices and are unaffected.

`tests/doc/tests_02_linalg.md` documented the old contract *and* the pre-R16a threshold
("imaginary part > 1e-15"). Both corrected.

### 24.2 The lesson

A negative claim about a codebase is only as good as the directories searched. "Nothing includes
this" needed `src/ tests/ nonfree/`, not `src/`. The failure mode is asymmetric: a false positive
gets checked, a false negative becomes a premise -- and this one was load-bearing for an API change.

**Still not compiled.** The tests are written but `make check` is Christian's to run; brace balance
and style are all that has been verified here.

### 23.4 The R15c sweep — a code bug among the flags

The third sweep returned three technical flags. Two were documentation errors; **one was a bug in
the source**. Running total across the three sweeps: nine caught issues.

**1. A latent type error in `fn_BZ_linspace.hpp`.** The three-argument overload that *returns* a
vector is declared `Vector<T>` but built a `Vector<real>`:

```cpp
template< typename T >
Vector <T>
linspace( const T & aStart, const T & aEnd, const belfem::size_t & aN )
{
    Vector <real> aValues;        // <- not Vector<T>
    ...
    return aValues;
}
```

It has never broken a build for a reason worth noting: **the overload is never instantiated**. Every
caller in the tree -- `cl_TensorMeshFactory.cpp:159,161,171` -- uses the four-argument fill form, and
an uninstantiated template body is only partially checked. The first person to call
`linspace<float>( ... )` would have met it. Fixed to `Vector<T>`; Armadillo has no equivalent
overload, so the bug was Blaze-only.

**2. `append`'s warning overstated the check.** It said the Blaze path asserts that the two vectors
differ. It does -- with `BELFEM_ASSERT`, which is compiled out of release builds, so in release
neither backend checks. The wording now says so, and notes that the const-reference overload carries
the same hazard.

**3. `sum` does not take a matrix.** The brief claimed "vector or matrix". The BELFEM API has
`sum( const Vector<T>& )` only; the matrix-shaped helper is Armadillo-native (`arma::Mat<T>`) and not
part of the wrapper surface.

Codex also independently confirmed the descending coefficient order from three separate sites,
including the Blaze `polyfit` assembly, and confirmed that both `unique` paths sort.

### 23.5 Hidden contracts, now stated

The one-line briefs were hiding preconditions a caller could only learn by reading the source:

- **`polyfit` requires `aX(0) != aX(n-1)`.** Both backends compute
  `tXref = ( n - 1 ) / ( aX(n-1) - aX(0) )` to condition the fit and then scale the coefficients
  back, so equal endpoints are a division by zero, not a diagnostic. Also documented: matching
  lengths, length greater than the degree, and that the coefficient vector is resized to degree + 1.
- **`linspace` needs `aN >= 2`** -- the spacing is the interval over `aN - 1`, and both endpoints
  are included.
- **`r2` returns exactly 1.0** when either the residual or the total variation falls below
  `BELFEM_EPSILON`, rather than dividing two near-zero numbers.
- **`combine`'s output is the last argument**, is resized, and must not alias any input.
- **`reverse` does not reverse in place** -- the public wrapper takes a const reference.
- `to_cell` / `to_vector` copy and share no storage; `max` / `min` require a non-empty input;
  the three polynomial evaluators require a non-empty coefficient vector.

`ddpolyval` gained the function-level block `dpolyval` already had, including the point that it
takes the coefficients of the polynomial itself rather than of either derivative.

Verified after all of it: 0 warnings, 0 errors.

## 25. R17 — compiled, run, green

`make check`, all 12 suites, 100% passed. The `linalg` suite compiled and ran the new work, and the
per-test log confirms the six tests executed rather than being silently absent:

```
Eigen.EigenComplexEigenvaluesAbortByDefault
Eigen.EigenComplexEigenvaluesReturnNaNWhenNotAborting
Eigen.EigenRealEigenvaluesReturnZeroCount
EigenSym.EigenSymReturnsAscendingRealValues
EigenSym.EigenSymDiagonal
EigenSym.EigenSymAgreesWithEigenOnASymmetricMatrix
```

That is the first executable gate over anything in this session's code changes, and it clears the
things that could only have been guesses until now:

- **R16b compiles.** The `int_t` return, the five-argument `BELFEM_ERROR` format string and the
  `Matrix< real > tWork( aMatrix.matrix_data() )` working copy in the Blaze path were all written
  against headers rather than a compiler.
- **The new contract behaves as documented.** `eigen` aborts on a complex eigenvalue by default and
  the abort is catchable as `std::runtime_error`; with `aAbortOnComplex = false` it returns 2 and
  fills NaN; a real-eigenvalue matrix returns 0.
- **`eigen_sym` works on the Blaze path**, including the `syev` working copy and the upper-triangle
  choice, and agrees with `eigen` on a symmetric matrix -- the test that would have caught a
  triangle or ordering mistake.
- **R16a and the `fn_BZ_linspace.hpp` type fix compile**, as do roughly thirty headers' worth of
  documentation comments added across R11 and R15.

Note what the green does *not* cover: this machine builds `USE_MATRIX_BLAZE=ON`, so the Armadillo
paths of `eigen`, `eigen_sym` and `linspace` were not compiled here. They are the mirror images of
what was tested and were written from the same headers, but an Armadillo build is what would prove
them.

**R17 closed. R1-R17 are complete.**

## 26. Where the session ended

| | start | end |
|---|---|---|
| Doxygen problem lines | 232 warnings + **1471 errors** | **0 / 0** |
| graph SVGs rendered | 0 of 1437 requested | **2737** |
| module index pages | 0 (all 24 severed) | **30** |
| API topic groups | none -- no `@defgroup` in the tree | **30**, 24 populated, 92 classes |
| class rows with a description | 26 / 318 | **146 / 318** |
| in-page guide links | 222 dangling | 0 |
| `src/linalg` headers documented | 5 of 44 | **36 under the topic**, every top-level `fn_*.hpp` |
| top navigation | Main Page -> Namespaces -> Classes -> Files | Getting Started, Modules, Documentation, Topics first |

Nine technical errors were caught by the three Codex prose sweeps, of which one was a source bug
(`fn_BZ_linspace.hpp`) and the rest were wrong claims in documentation written the same day. Two
further code defects came out of reading bodies to document them (the `eigen` tolerance mismatch,
the non-inline header functions), and one design change was made and tested (`aAbortOnComplex`,
`eigen_sym`).

## 27. Postscript — the Armadillo gap, and what it cost

The campaign closed with one honest gap: this machine builds Blaze, so the Armadillo paths of
`eigen`, `eigen_sym` and `linspace` had only ever been written against the installed headers. It was
flagged in this devlog, in the plan and in the commit message.

Flagging it was not the same as closing it, and the gap was real: **`eigen_sym` did not compile on
Armadillo.** `arma::eig_sym` writes into an `arma::Col&`, while BELFEM's `Vector<T>` is backed by an
`arma::Mat` (`cl_AR_Vector.hpp:38`), so `vector_data()` returns a `Mat&` and the only route to `Col`
was a user-defined conversion producing an rvalue, which cannot bind to a non-const reference. Fixed
by collecting into a local `arma::Col` and copying out, the way the general `eigen()` in the same
file already handles `arma::cx_vec`.

The lesson is not "the gap was known" -- it was -- but that **`build/compile_commands.json` was
sitting there the whole time**. A `-fsyntax-only` compile with the real flags costs seconds, touches
nothing, and would have caught it before the commit. "The user runs builds" is about not disturbing
a shared build tree; it was never a reason to skip a check that writes nothing.

Both were then done properly: the probe translation unit and the real compile line for
`tests/linalg/test_LinalgSolvers.cpp` compile clean, and a scratch tree configured for Armadillo
builds, links and runs the suite **282/282, all ten eigen tests included**. The Blaze and Armadillo
paths are now both proven.

### 27.1 A linker warning, traced

`ld: warning: ignoring duplicate libraries: '-lgfortran', '-lgomp'` turned out to be three
independent double-adds, each pairing an explicit `-l` with a mechanism that already supplies the
same library: CMake's `CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES` against `BELFEM_FORTRANLIBS`, and
`-fopenmp` against both `BELFEM_OPENMPLIBS` and a third addition in `config_mkl.cmake`.

That third one is worth knowing independently of the warning: **`config/linalg/config_mkl.cmake` is
included unconditionally** (`CMakeLists.txt:140`), so its trailing runtime-libraries block executes
whenever `USE_OPENMP` is on, `USE_MKL` off or not.

Each addition is now conditional on the mechanism that would duplicate it, rather than deleted --
on a toolchain where CMake does not detect the Fortran runtime they are load-bearing. Verified: the
link line carries `-lgfortran` once and no explicit `-lgomp`, and `test_linalg` builds with no
linker warning of any kind.

## 28. Commits

Six, on `betterdoc`: the Doxygen defaults and toolchain; the hand-written navigation tier; the API
topic tree and `src/linalg` documentation; the `eigen` behaviour change with its tests; these
records; and the Armadillo `eigen_sym` fix, followed by the link-line cleanup.


`todo/doxygen_site_structure.md` carries the rest: R10 (MathJax) and R12 (residual warnings,
undocumented modules), plus continuing R11 across the other modules.

### Note on R10, measured 2026-08-18

The tree contains **374** `\f$` / `\f{` formula blocks, all in headers (`cl_IWG_Timestep.hpp:316`,
the `IFG_*` bubble-function headers), so MathJax is load-bearing, not decoration. The site pulls it
from `cdnjs.cloudflare.com/.../mathjax/2.7.5/` with `MATHJAX_FORMAT = HTML-CSS`; `MATHJAX_VERSION`
is absent from `Doxyfile.in` and therefore defaults to `MathJax_2`. Christian's call is that a
network dependency is acceptable because the documentation will be hosted online. The residual risk
is not offline use but the pin: MathJax 2 is end-of-life, and if that cdnjs path is retired the 374
formulas fail silently with no build error.
