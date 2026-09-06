# Doxygen Site Structure and Navigation

**Date:** 2026-08-18
**Purpose:** Make the generated Doxygen site navigable for three readers — a new user who wants to
run a simulation, a new developer who wants to extend a module, and a returning developer looking
something up. Wave 1 (the three broken defaults) has landed; this plan carries the remaining
navigation and content work.
**Module:** meta (documentation toolchain) — `Doxyfile.in`, `scripts/update_doc_index.py`,
`doc/mainpage.md`
**AIs involved:** Claude (sweep + empirical verification + plan), Codex (blind sweep), Grok (blind sweep)
**Status:** ✅ COMPLETE (2026-08-18). R1–R17 all done; see
`devlog/dl20260818_doxygen_site_structure.md`. The Doxygen log is at **0 warnings / 0 errors**, down
from 232 warnings + 1471 errors; 2737 graphs render where none did; the API reference has a 30-group
topic tree where it had none; and `make check` passes 12/12, with six new tests covering the `eigen`
behaviour change. Residual follow-up **closed 2026-08-18**: the Armadillo paths were compiled and run in an
Armadillo-configured tree — 282/282 pass, all ten eigen tests included. `eigen_sym` needed one fix
to get there (`arma::eig_sym` writes into an `arma::Col&`, but `Vector<T>` is backed by an
`arma::Mat`), which is the defect the gap existed to catch.

> **Scope guards:**
> - Graphviz discovery is OUT of scope. `HAVE_DOT = YES` stays; `dot` belongs on the developer's
>   `PATH` (decided 2026-08-18, Christian — `/opt/graphviz/latest/bin` is machine-specific).
> - Replacing `scripts/update_doc_index.py` or adopting a different doc toolchain is OUT of scope.
>   Its `{#anchor}` + `.dox`-tree design is correct; all three AIs agreed.
> - `EXTRACT_ALL = YES` stays. The signal-to-noise problem is an entry-path problem, not a reason
>   to hide API.

---

## Gap table

| # | Gap | Reader hurt | Evidence | Status |
|---|---|---|---|---|
| G1 | Module index READMEs absent from the site | all three | 23 dead `@subpage`, no `*_index.html` | ✅ R1 |
| G2 | Usage-guide TOCs dead | developer | 222 unresolved `\ref` | ✅ R2 |
| G3 | No diagrams | developer | 1471 rename errors, 0 svg | ✅ (Graphviz installed) |
| G4 | Class list has no descriptions | returning dev | `annotated.html` empty `<td class="desc">` | ✅ R3 |
| G5 | Error-blind quality gate | maintainer | `report_doxygen_warnings.cmake` matched `warning:` only | ✅ R4 |
| G6 | No pages tab in the top chrome | all three | `menudata.js`: Main Page / Namespaces / Classes / Files only | ✅ R5 |
| G7 | `docmap` is a flat path-titled list | developer | `docmap.html`, vs the four buckets in `doc/README.md:20-54` | ✅ R6 |
| G8 | Landing page mixes audiences | new user | `doc/mainpage.md:54-63` lists AI-process docs beside Getting Started | R7 |
| G9 | `examples/` invisible to the site | new user | `INPUT` covers only `doc/`+`src/`; `EXAMPLE_PATH` empty; `getting_started.md:70-71` points outside | ✅ R8 |
| G10 | File list polluted by test drivers | returning dev | `petsctest.cpp`, `poisson.cpp`, `corctest.cpp`, 5 `main.cpp` | ✅ R9 |
| G11 | MathJax + Mermaid from CDNs | offline reader | `MATHJAX_RELPATH = https://cdnjs...`, `cdn.jsdelivr.net/npm/mermaid@11` | R10 |
| G12 | No topical API grouping at all | developer | `grep @defgroup src` → **0** | R11 |
| G13 | Guides and API do not cross-link | returning dev | no `@see`/`@ref belfem::` in any `.md` | R11 |
| G14 | thermal / postproc / visualizer undocumented | developer | no `doc/` dir; absent from `doxygen_nav.dox` | R12 |
| G15 | 19 residual doc-comment warnings | maintainer | see devlog §6 | R12 |

---

## Steps

### Wave 1 — the broken defaults (COMPLETE, verified 2026-08-18)

- [x] **R1** `IMPLICIT_DIR_DOCS = NO` in `Doxyfile.in:937`. Restores all 27 `*_index.html` pages and
      the whole `docmap → mod_* → *_index` chain. *(Claude; verified by full run + link walk.)*
- [x] **R2** `MARKDOWN_ID_STYLE = GITHUB` in `Doxyfile.in:359`. 222 unresolved `\ref` → 1.
      *(Claude; verified programmatically — DofManager guide 24 links, 0 dangling.)*
- [x] **R3** `JAVADOC_AUTOBRIEF = YES` in `Doxyfile.in:204`. Class rows with a description
      26/318 → 85/318. *(Claude; A/B tested standalone, 0 regressions.)*
- [x] **R4** `cmake/report_doxygen_warnings.cmake` counts `error:` as well as `warning:`.
      *(Claude; exercised against the new log.)*

### Wave 2 — navigation (open)

- [x] **R5** `doc/DoxygenLayout.xml` checked in, `LAYOUT_FILE` set (`Doxyfile.in:777`). Tab bar is now
      `Main Page | Getting Started | Modules | Documentation | Namespaces | Classes | Files`, present on
      every page including deep class pages. *(Claude; verified — both variants built and compared,
      see O1.)*
- [x] **R6** `docmap` grouped into the four `doc/README.md` buckets with human titles and one-line
      blurbs; module pages retitled (`Maxwell`, not `src/fem/maxwell`) with the path kept in the body.
      `MODULE_GROUPS` in `scripts/update_doc_index.py` is the single table, and `main()` now fails if a
      documented module has no entry, so a new module cannot silently vanish from the catalogue.
      *(Claude; verified — 4 `<h1 class="doxsection">` groups render, 23 modules classified.)*
- [x] **R7** Landing page emitted in three reader paths (Run a simulation / Extend a module /
      Project process) from `PAGE_GROUPS`, with a guard that fails if a page is unfiled.
      *(O2 decided 2026-08-18, Christian: keep the AI-process documents, group them.)*
- [x] **R8** `examples/` added to `INPUT`, `EXAMPLE_PATH` set, `anchor_for()`/`collect_examples()`
      taught about the tree (it has no `doc/` subdirectory), and `doc/getting_started.md:70-73` now
      links `@ref doc_examples` instead of naming a repository path. Pages `doc_examples.html` and
      `doc_examples_scripts.html` render. *(Claude; verified.)*
- [x] **R9** `EXCLUDE_PATTERNS` for the in-tree drivers (`*/main.cpp`, `*/*test.cpp`,
      `*/test_*.cpp`, `*/poisson.cpp`). All 12 candidates were confirmed to carry their own `main()`
      before exclusion. Verified after: the 12 drivers are gone from the file index, `hphirun`,
      `hphiTrun` and `electricalCircuit` remain, and `cl_IWG_Poisson.cpp/.hpp` — the real physics —
      is untouched. *(Claude.)*
- [x] **R10** `MATHJAX_VERSION = MathJax_3`; `MATHJAX_RELPATH` returned to empty and documented as
      needing to stay empty (the version and the path are a matched pair — a pinned v3 path read by
      an older Doxygen 404s and blanks all 374 formulas silently). Removes the end-of-life
      cdnjs 2.7.5 pin. Vendoring was considered and rejected: not needed for an online site, and on
      MathJax 3 it would be one ~1 MB file if that ever changes. *(Christian confirmed Doxygen
      1.18.0 is being built on the RHEL8 server; Claude verified the emitted URL and that an
      unsupported tag warns on stderr only, never into WARN_LOGFILE.)*
- [x] **R13** `doxygen -u` migration of the 1.8.16 template to 1.18.0. Reviewed by parsing every
      `KEY = value` before/after rather than by reading the textual diff: 17 removed (all obsolete),
      1 changed (`EXTRA_PACKAGES` separator, inert with `GENERATE_LATEX = NO`), 50 added — all at
      Doxygen's own defaults, which were already in force. All BELFEM settings and all 12 `@VAR@`
      placeholders survived; the 8 rationale comment blocks did not and were re-injected, now marked
      `# BELFEM:`. Verified: 14+3 stderr warnings → 0, warning log byte-identical, and a full
      recursive diff of the 2879-file output tree reports **0 differing files**. *(Claude.)*
      - Newly visible and deliberately left alone: `HTML_COLORSTYLE = AUTO_LIGHT` (`TOGGLE` gives a
        dark-mode switch), `MERMAID_JS_URL` (the second CDN, unused today),
        `NUM_PROC_THREADS = 1` (0 uses all cores).

### Wave 3 — content (open)

- [x] **R11** Topic tree for the API reference. `doc/groups.dox` is generated by
      `scripts/update_doc_index.py` from `MODULE_GROUPS` (same table as the docmap catalogue), so
      the tree cannot drift: 4 top-level groups + 26 module groups. Classes join with `@ingroup` in
      their own doc block; populating is incremental and split below by catalogue bucket.
      **COMPLETE: 24 of 26 module groups, 92 classes. Class descriptions 26/318 → 146/318.**
      The two empty groups, `grp_fem` and `grp_executables`, are correct — neither directory holds
      any classes. Do not "fix" them.
      - Standing rule: the `cl_{AR,BZ}_*` linalg classes merge into one Doxygen entity, so a
        backend-specific brief is wrong for half of all builds — keep those briefs neutral.
      - `grp_fem` legitimately has no members: `src/fem` has no top-level headers. It renders as an
        overview group. Do not "fix" it.
  - [x] **R11a** Finite Elements — `fem/iwg` (`IWG`, `IWG_Timestep`, `TimestepMatrices`,
        `IwgFactory`), `fem/interpolation` (`InterpolationFunction`,
        `InterpolationFunctionFactory`, `IntegrationData`, `EdgeFunction`, `EdgeFunctionFactory`),
        `fem/thermal` (3), `fem/postproc` (3). Briefs from each module README's key-class table.
        *(Claude; verified — 0 warnings, 0 errors.)*
  - [x] **R11b** Mathematics — `math/graph` (`Vertex`), `math/tensor` (`Tensor`),
        `math/quaternion` (`Quaternion`), `numerics/spline` (`Spline`), `homology` (`CutFactory`,
        `Cohomology`, `Homology`, `SimplicialComplex`, `CutProcessor`, `Chain`, `Cochain`,
        `CutData`). Each homology class points at the document that explains *it* — theory,
        algorithms or thick/thin cuts — rather than all at the module index.
        *(Claude; verified — 0 warnings, 0 errors.)*
  - [x] **R11c** Physics and Circuits — `physics/materials` (4), `physics/database` (2),
        `physics/gasmodels` (7, incl. the concrete `EoS_*` and `Helmholtz` variants),
        `physics/gastables` (5), `circuit` (5). 18 of the 23 merged into documentation blocks that
        already existed, so those keep their own briefs and gained only `@ingroup`/`@see`.
        `grp_executables` stays empty: `src/executables` holds three `.cpp` files with `main()` and
        no classes — like `grp_fem`, it is an overview group. *(Claude; verified — 0/0.)*
  - [x] **R11d** Infrastructure remainder — `core` (6), `io` (6), `comm` (1), `visualizer` (2).
        Also corrected an error introduced in R12: the visualizer README listed *file* names as
        class names; the classes are `belfem::vtk::MeshView` / `Curve`, and `BlockActor` /
        `SideSetActor` were missing from the table entirely. The `linalg` free functions moved out
        to **R15** — they need documentation, not grouping. *(Claude; verified — 0/0.)*
- [x] **R12** All 19 residual warnings closed and the log locked at 0/0. Nine documentation
      warnings (a `## Index` heading colliding with the main page; five `[Source code](../)` links
      to a bare directory, reported misleadingly as `autotoc_md`; two relative links out of the doc
      tree; one genuine TOC/heading drift in `iwg_usage_guide.md`) and ten header warnings
      (`@return` on void, `@param` on no-arg accessors, and four caused by `cl_Material.hpp`
      repeating the `@param`/`@return` that the `powerlaws.hpp` definitions already carry).
      READMEs written for `src/fem/thermal`, `src/fem/postproc` and `src/visualizer`, each naming
      what it does *not* document. `BELFEM_DOXYGEN_WARNING_BASELINE` lowered to 0.
      *(Claude; verified — full run, 0 warnings, 0 errors.)*
      - Flagged, not fixed: `IWG::collect_nodes_on_wetted_sitdesets` has the typo in the **function
        name**. Renaming touches call sites and is a source change, not a doc fix.
      - Trade-off: the five module READMEs lost a `../` directory link on the GitHub web view.
        Doxygen cannot link a directory; the site's Files index covers it.

- [x] **R14** Project logo links to `https://belfem.lbl.gov` (feature request, 2026-08-18).
      Implemented as a build-time-generated `HTML_HEADER` (`doxygen -w html` +
      `cmake/patch_doxygen_header.cmake`), never a checked-in header file: Doxygen does **not** warn
      when a custom header is stale, so a committed one would silently lose the treeview, search or
      MathJax on a version upgrade. The patch script aborts loudly if Doxygen's logo markup changes.
      Verified: link on 2893/2894 pages, treeview + search + MathJax intact, 0 warnings / 0 errors.
      *(Claude.)*
      - If the site ever needs more header customisation, extend the patch script — do not switch
        to a committed header.

- [x] **R15** Document the linalg public API. `src/linalg` is the largest undocumented public
      surface in the tree: **zero** `@brief` across all 44 files, roughly **97 declarations**.
      Decided 2026-08-18 (Christian): the documentation goes in the **top-level dispatch headers**
      (`src/linalg/fn_*.hpp`), not in the per-backend implementations — `fn_dot.hpp` is "the dot
      product header" regardless of which backend it includes, and documenting it once at the
      public surface avoids two copies that can disagree.
      Mechanism, verified 2026-08-18: an `@fn` block in the dispatch header binds to the backend
      declaration cleanly, **including** when both backends are parsed, and with no warnings — the
      trailing `-> decltype(...)` return type does not have to be repeated in the `@fn` signature.
      Add a `@file` block (`@brief` + `@ingroup grp_linalg`) plus one `@fn` per documented overload.
      - Known artifact, accept rather than fight: the rendered signature keeps whichever backend
        Doxygen parsed first, so a Blaze build still shows `-> decltype(arma::dot(...))`. State the
        result in words via `@return`; do not try to correct the signature.
      - Briefs stay backend-neutral, same rule as `Vector` / `Matrix`.
      - Land per function family so each batch is reviewable, and send each batch's prose to Codex
        for a readability sweep (see `CLAUDE.md` §"Prose Gets a Language Sweep").
  - [x] **R15a** `dot` (3 concrete overloads documented; the 21 expression-typed `ET` variants are
        explained once at file level rather than individually), `cross`, `norm` via `@fn` in the
        dispatch headers; `crossmat` (4 overloads) and `trans` documented in place, both being
        direct implementations rather than dispatchers. Verified 0/0, briefs bind on the namespace
        page. Codex readability sweep applied; it flagged **two wrong technical claims** in the
        first draft (the epsilon zeroing is 2D-only, and `dot(Matrix,Vector)` returns a scalar, not
        a vector), and checking them surfaced two more (the `aScale` overloads accumulate rather
        than assign; there is no assert on the row count of `aA`). All corrected, 0/0 after.
        - Recorded from reading the body: every `crossmat` overload zeroes result entries below
          `BELFEM_EPSILON` relative to the result norm. Cosmetic, surprising, now documented.
  - [x] **R15b** Solvers. Split once the tree was actually measured: `src/linalg/lapack/` is
        **9 of 10 files already documented** but had **no `@ingroup`**, so ten headers only needed
        a `@file` block to reach the topic tree; `inv`, `inv2`, `inv3`, `det`, `eigen` needed real
        documentation and got it. The Linear Algebra topic now lists 20 headers, up from 5.
        Verified 0/0. Codex sweep applied; it flagged **four more technical claims**, all correct:
        the `eigen` tolerance differs by backend, `gels` also does minimum-norm, `posv` is symmetric
        *or* Hermitian, and "yields infinities" was too specific for a release-build singular case.
        - Recorded from reading the body: `eigen` writes `BELFEM_QUIET_NAN` for any eigenvalue with
          a non-negligible imaginary part — silently, no error, no return code. Now a `@warning`.
        - `inv2`/`inv3` return the determinant on purpose (element Jacobians need it), and their
          singularity check is relative *and* an assert, so it vanishes in release builds.
        - Note for R15c: the R15 premise "zero briefs across 44 files" was measured without
          `src/linalg/lapack/`. Re-measure before assuming a file is undocumented.
  - [x] **R15c** Utilities, aggregates, polynomials, statistics. Measured first, per the R15b note:
        16 files, including `fn_to_cell` and `fn_to_vector` which were not on the original list.
        All 16 have a `@file` block; `polyval`, `dpolyval`, `sort`, `unique`, `append` also have
        function-level documentation. Linear Algebra topic now lists **36 headers**, up from 5.
        Verified 0/0. Codex sweep applied. It flagged three items, one of which was **a bug in the
        source**: `fn_BZ_linspace.hpp`'s returning overload declared `Vector<T>` but built a
        `Vector<real>` — invisible because the overload is never instantiated. Fixed. The other two
        were doc errors (`append`'s assert compiles out in release; `sum` takes no matrix). Hidden
        preconditions now documented for `polyfit`, `linspace`, `r2`, `combine`, `reverse`,
        `to_cell`/`to_vector`, `max`/`min` and the polynomial evaluators.
        - The four polynomial headers state that coefficients are **descending** (`aCoeffs(0)` is
          the leading coefficient), read out of the Horner loops. Getting it backwards gives a
          plausible wrong answer, not an error.
        - `unique` **sorts as a side effect** on both backends — not an order-preserving filter.
        - `append` requires `aA != &aB`: Blaze asserts, Armadillo would silently double the vector.

- [x] **R16** Two code issues surfaced by documenting `src/linalg`; both are behaviour, not prose.
  - [x] **R16a** `eigen` tolerance unified on `BELFEM_EPSILON` (decided 2026-08-18, Christian).
        `fn_BZ_eigen.hpp:41` was comparing against a hard-coded `1e-15` while
        `fn_AR_eigen.hpp:39` used `BELFEM_EPSILON` (~2.2e-15) — the same decision with two
        thresholds. It was the only hard-coded epsilon literal under `src/linalg`. The R15b `@note`
        documenting the discrepancy was rewritten in the same edit. *(Claude; Doxygen clean, but
        **not compiled** — nothing in the tree includes `fn_eigen.hpp`.)*
  - [x] **R16b** `eigen` given `aAbortOnComplex = true` (BELFEM_ERROR naming the index and its
        imaginary part; false returns the count and fills NaN), matching `gesv`/`posv`. `eigen_sym`
        added because both backends have a symmetric solver (`arma::eig_sym`, `blaze::syev`), which
        removes the question entirely for symmetric matrices. *(decided 2026-08-18, Christian.)*
        - Fixed on the way: both `eigen` definitions were **non-inline in a header** (a link-time
          collision waiting for a second includer); the Blaze wrapper initially read the **lower**
          triangle while `arma::eig_sym` reads the **upper**, which would have made the backends
          disagree on a non-symmetric input; `blaze::syev` overwrites its input, so the Blaze path
          takes a working copy to keep the const interface.
        - Compiled and tested under R17: `make check` 12/12 green, six new tests covering the
          abort, the count and `eigen_sym`.
- [x] **R17** `fn_eigen.hpp` build coverage — **closed by a real test run.** `make check`: 12/12
      suites, 100% passed. The per-test log confirms all six new tests executed
      (`EigenComplexEigenvaluesAbortByDefault`, `…ReturnNaNWhenNotAborting`,
      `EigenRealEigenvaluesReturnZeroCount`, `EigenSymReturnsAscendingRealValues`,
      `EigenSymDiagonal`, `EigenSymAgreesWithEigenOnASymmetricMatrix`). First executable gate over
      R15/R16's code changes. *(run by Christian 2026-08-18.)*
      - The premise had been wrong: the header *is* included by
        `tests/linalg/test_LinalgSolvers.cpp:30` and built via `tests/linalg/CMakeLists.txt:9`. The
        earlier "nothing includes it" came from grepping `src/` only, and R16b consequently broke
        `TEST( Eigen, EigenComplexEigenvaluesReturnNaN )`, which asserted the old silent-NaN
        default. That test was rewritten to expect the abort and five more were added.
      - Not covered: this machine is `USE_MATRIX_BLAZE=ON`, so the **Armadillo** paths of `eigen`,
        `eigen_sym` and `linspace` were not compiled. An Armadillo build would close that gap.
