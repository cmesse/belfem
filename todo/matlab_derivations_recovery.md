# Recover the MATLAB derivation scripts into the documentation tree

**Date:** 2026-09-16
**Purpose:** The pre-reset `tmp/` (now `/home/christian/codes/belfem_transfer/tmp`, mirrored read-only under `tmp/matlab_survey/` for the auditors) holds 56 MATLAB scripts. A few generated tables that the C++ now carries and cannot regenerate on its own (the TET10 Nédélec tables) or verify them symbolically (the PENTA18 second-derivative table); a few more are the only written-down derivations behind a convention the C++ states as a result. One document cites two of them by paths that no longer exist. Keep those, clean them up to the commenting guideline, document them, file them under `src/fem/interpolation/doc/matlab/`, and add asserting drivers so that running them is a gate. Everything else is dropped; the transfer tree keeps the originals.
**Module:** `src/fem/interpolation/doc`
**AIs involved:** Claude (survey via a read-only subagent, plan), Codex `gpt-6-astra`/high + Grok `grok-4.6`/high (jury on the keep set, 2026-09-16), Codex (language sweep of the kept scripts' comments)
**Status:** ✅ COMPLETE (2026-09-16) except O1, which is Christian's call and defaults to drop. Landed: 14 scripts under `src/fem/interpolation/doc/matlab/` (TET10 generator and zero-tests, DefElement cross-check, TET10 Lagrange derivation, four triangle derivation notes, the PENTA18 zero-test, two PENTA6 facet notebooks), two asserting drivers, `compare_tables.py`, five READMEs, `*/doc/matlab/*` excluded from Doxygen, the two dangling citations in `nedelec_derivation.md` repaired, the directory registered in the module doc index. Verified: `check_tet10` and `check_penta18` exit 0 under Octave 9 / `symbolic` 3.1.1 (MATLAB is not installed here), `compare_tables.py` 0 differences over 192 + 108 entries, every runnable note runs clean; the stripped-comment diff against the originals is empty for all 13 unchanged scripts and shows only the preamble for `fragment.m`. Codex language sweep applied (four content errors it flagged fixed: `facet.m` needs the toolbox, `secondordertri.m` prints nothing, a misplaced "column vector" comment, and the DefElement URL is stated as the triangle page). The jury cut the survey's 26 keepers to 14 and withdrew every "exact generator" label but two; 42 scripts stay in the transfer tree. Residual: O1; the `tests/` dated lines and `tests/fem/test_EdgeFunctions.cpp:32`'s `tmp/` path go to comment sweep 2.

> **Scope guards:**
> - No C++ changes. The scripts document and verify tables; they do not alter them. If a script disagrees with a C++ table, that is a finding for a separate plan, not a fix here. (The jury caught R5 breaking this guard; the C++ half of R5 is withdrawn.)
> - One build-documentation change is required and carved out explicitly: `Doxyfile.in` lists `*.m` in `FILE_PATTERNS` with `INPUT = src`, `RECURSIVE = YES`, so the new directories must be added to `EXCLUDE_PATTERNS` before anything is copied (R2a). Without it `make doc` parses MATLAB as source.
> - Scripts are kept as MATLAB; no port to python or Octave. The Symbolic Math Toolbox dependency is stated, not removed.
> - Executable additions are limited to (i) a six-line preamble that makes `fragment.m` runnable, (ii) two new asserting driver scripts (`check_tet10.m`, `check_penta18.m`), and (iii) one python comparison of the embedded MATLAB table copies against the live C++ tables. Every other edit to a kept script is comment-only, checked by diffing against `tmp/matlab_survey/`.
> - MATLAB is not installed on this machine (the `matlab` shell alias points at a directory that does not exist); Octave 9 with the `symbolic` package (sympy-backed) is, so the drivers run there (`octave-cli --eval "pkg load symbolic; check_tet10"`) as the closing gate (R7). A MATLAB run, when one is available, is the same command with `matlab -batch`. Until R7 has run, the kept scripts are **reviewed, not verified**.
> - Comment cleanup follows `doc/commenting_guidelines.md` (license header, no dates, no bylines, no working-record paths). The guideline covers C++ and Fortran; MATLAB's comment character `%` is used by analogy, the guideline is not amended.
> - Runtime verification of the tables already exists and is not replaced: `tests/fem/test_LagrangeInterpolation.cpp` (central-difference second derivatives, `PENTA18` at `:276`), `tests/fem/test_FacetIntegrationPoints.cpp` (every slave orientation, `TET10` at `:215`), `tests/fem/test_EdgeFunctions.cpp` (circulation). What this plan preserves is the symbolic generation and derivation, which no runtime test can regenerate.

---

## 1. Current Behaviour and How It Fails

- `src/fem/interpolation/doc/nedelec_derivation.md:270` (§3.2, the TRI3 curl) cites "the MATLAB fragment `matlab_curl/fragment.m` in the notes repository"; `:443` cites "`tmp/tet10/tet10_generate.m`". Neither path exists in this repository. The TRI3 file the first citation meant is gone from the transfer tree too; the closest survivor is `nedelec/matlab_curl/curl.m`, a scratch derivation of the same operator with the opposite sign convention from `cl_EF_TRI3.cpp:132`. The second citation's file exists and is the generator. Confidence high (lines read; transfer tree listed).
- `cl_EF_TET10.cpp`'s `mG/mH/mU/mV/mW` tables and their fifteen derivative tables were generated by `tet10_generate.m` and are checked against a MATLAB transcription by `tet10_function.m` / `tet10_derivatives.m` (eight of the twelve face rows: the generator emits three candidates per face and the C++ stores all three, the zero-tests compare the two active ones). `cl_IF_PENTA18.hpp::d2NdXi2` (6×18, hand-written) has an exact symbolic zero-test in `penta18.m`; a whitespace-normalised compare of all 108 entries against the header found 0 differences (Claude, pre-jury). None of these generators or symbolic checks is in the tree.
- The old `tmp/` also holds 42 scripts that are duplicates, plots, scratch tests, abandoned alternatives, unfinished derivations, or derivations of a route the C++ does not use (the jury's list, §3).

**Bottom line:** the value is in the TET10 pipeline and `penta18.m`; the rest of the keep set is derivation notes and is labelled as such.

## 2. Architecture: scripts live next to the document that explains the table

A kept script goes under the module `doc/` directory whose prose already explains its table (`src/fem/interpolation/doc/`, `src/mesh/doc/`), in a `matlab/<topic>/` subdirectory with one `README.md` per subdirectory that states, per script: what it generates or verifies, which C++ table or function it pins, what it needs (Symbolic Math Toolbox, another script run first), and how to run it. The prose documents then cite the scripts by their in-tree path. `scripts/` is wrong for them: they are not tooling anyone runs routinely. The directories are excluded from Doxygen (R2a); they are documentation for a reader of the C++, not API.

## 3. Gap Table (the keep set after the jury)

Class: (a) the script reproduces a C++ table exactly, with an asserting driver; (b) derivation note: the only written-down form of a convention or a derivation the C++ states as a result, not a table dump; (c) support file a keeper needs.

| # | Script (old path) | What it is, and what in the C++ it concerns | Class | Target |
|---|---|---|---|---|
| 1 | `tet10/tet10_generate.m` | generator of the `mG/mH/mU/mV/mW` tables in `cl_EF_TET10.cpp::precompute()` (not of the Lagrange `mNxi/mNeta/mNzeta` there, which come from `tet10_lagrange.m`); carries the `lambda = [xi; zeta; eta; tau]` node-map pin; cited at `nedelec_derivation.md:443` | a | `matlab/nedelec_tet10/` |
| 2 | `tet10/tet10_shortcuts.m` | the `xi8`, `eta16`, … shorthand the C++ tables use verbatim; required by 3 and 4 | c | same |
| 3 | `tet10/tet10_function.m` | symbolic zero-test of `mG/mH` (12 rows each) and the 8 active face rows of `mU/mV/mW` against the generator | a | same |
| 4 | `tet10/tet10_derivatives.m` | the same for the fifteen derivative tables | a | same |
| 5 | `tet10/parse_main.m` | prints one table (`mWzeta`) in C++ syntax from the generator; the pattern for the others | c | same |
| 6 | `tet10/defelement.m` | cross-check of the `g,h,u,v,w` construction against the published DefElement expressions (the URL is the triangle page; the body is the tet) | b | same |
| 7 | `tet10/tet10_lagrange.m` | Vandermonde construction of `cl_IF_TET10.hpp::N` in EXODUS node order, the written-down `xi_hat/eta_hat/zeta_hat` pinning | b | same |
| 8 | `nedelec/matlab_curl/curl.m` | scratch derivation of the TRI3 curl operator; the collapse to `C = 2/detJ [s1 s2 s3]` is left commented, and its Whitney sign is opposite to `cl_EF_TRI3.cpp:132`; the `:270` citation's nearest survivor | b | `matlab/nedelec_tri/` |
| 9 | `nedelec/triangle/fragment.m` | the TRI6 edge/face ansatz `theta = a*nabla_j + b*nabla_i` with `a = 4ξ_i(2ξ_i−1)`, `b = 2ξ_j(1−4ξ_i)`; `cl_EF_TRI6.cpp:93ff` carries the expanded polynomials with half these coefficients, so this is the ansatz, not the generator; a loop body that needs a preamble | b | same |
| 10 | `nedelec/matlab_curl/secondordertri.m` | the two-point-per-edge `(a,b)` ansatz behind 9 | b | same |
| 11 | `nedelec/matlab_curl/check_quad.m` | justification of the `0.5*(1±xi)` remap between `[0,1]²` and `[-1,1]²` used at `fn_IF_initialize_integration_points_on_facet.cpp:666-669` | b | same |
| 12 | `penta18.m` | exact symbolic zero-test of `cl_IF_PENTA18.hpp::d2NdXi2` (all 108 entries) | a | `matlab/lagrange/` |
| 13 | `penta/orientation.m` | the PENTA6 (face, slave orientation) parameter maps `eta = f(xi)` behind the slave `intpoints_penta` overload (`fn_IF_initialize_integration_points_on_facet.cpp:635`), as an interactive plot notebook; `elements/penta.m` is its whitespace twin, `penta/surface.m` a truncated earlier version | b | `matlab/facets/` |
| 14 | `penta/facet.m` | which Jacobian rows give the outward normal on each PENTA6 face, as a plot notebook | b | same |

New files: `matlab/nedelec_tet10/check_tet10.m` and `matlab/lagrange/check_penta18.m` (asserting drivers, R3), `matlab/compare_tables.py` (embedded MATLAB copies vs the live C++ tables, R3), one `README.md` per directory (R4).

Dropped (42, all still in the transfer tree): byte-identical duplicates (`tet/tet.m`, `tet10/meshtest/main2.m`, `tet10/interface/init_triangle.m`), the whitespace twin `elements/penta.m` and the truncated `penta/surface.m`; plots and fixtures (`draw_mesh.m`, `testmesh.m`, `orientation_mesh.m`, `set_phi.m`, `make_mesh.m`, `shape_functions.m` whose `Ee2` disagrees with the generator); scratch (`test.m`, an unrelated Kepler fit; `curltest.m`; `penta_orientation.m`, three lines); trivially reproducible algebra (`hex8ts.m`, `pyra/pyra5.m`); abandoned alternatives (`points.m`, `find_shape.m`, `pyra/pyra13.m`, `tet10/nabla_derivatives.m`, `invert_coords.m`); superseded and internally inconsistent (`meshtest/nedelec.m` with the `[3 2 3]` face row, `nedelec_gh.m`, `nedelec_uvw.m`); 2016–2017 generic Legendre/Gauss routines nothing calls (`legendre.m`, `dlegendre.m`, `legendre_zero.m`, `init_gauss.m`); `sideset.m` (no C++ counterpart); the jury's demotions: `pyra/pyra14.m` (stops at the Vandermonde matrix, base at `z = -1` against the C++ base at `zeta = 0`, rational basis against the C++ polynomial one), `elements/hex.m` (`S2b == S2d`, a wrong orientation inventory), `elements/pyra.m` and `elements/tet.m` (plot notebooks that emit no table; the runtime test covers every orientation), `check_linear.m` (integrates over the unit square, not a circulation), `tet/main.m` (`P10/P20/P35` have no consumer in `src/`), `quadcurl_symbolic.m` (bilinear quad; the C++ `mCurv` is the quadratic map's Hessian), `nabla_derivative.m` (uses its cofactors before defining them; the cofactor route is not the C++ route), `tet10/meshtest/init_points.m` and `init_triangle.m` (the magic-table check needs a driver that only the dropped `main.m` provides, tests one of twelve columns, and `tests/fem/test_FacetIntegrationPoints.cpp` checks all of them at runtime); and the three interface T-matrix experiments pending O1.

### 3.1 Cross-cutting findings

- The three TET10 face tables in the old scripts disagree with each other: `meshtest/nedelec.m` `[0 2 3; 3 2 3; 0 3 1; 0 1 2]`, `nedelec_uvw.m` `[0 2 3; 2 1 3; 0 3 1; 0 1 2]`, `tet10_generate.m` `[0 1 3; 1 2 3; 0 3 2; 0 2 1]`. Only the last produced the current C++. Keeping a wrong table next to the right one is worse than none; the same rule dropped `shape_functions.m` and `hex.m`.
- A script that derives a different route than the C++ takes (`nabla_derivative.m`, `quadcurl_symbolic.m`) would mislead a reader about what the C++ does; dropped for that reason, not for being wrong.
- Every date in the kept scripts ("corrected 2026-08-14: eta and zeta exchanged") is history and goes; the pin itself (`lambda = [xi; zeta; eta; tau]`, with the reason that the naive order is left-handed on an EXODUS tet) stays as the contract.
- The zero-tests compare the generator against MATLAB transcriptions of the C++ tables, not against the C++ file. `compare_tables.py` closes that gap textually (as the pre-jury PENTA18 compare did).

## 4. Ordered Steps

- [x] **R1 — Jury on the keep set** (Codex + Grok, blind, 2026-09-16): 33 findings, 31 confirmed by file:line, two refuted; keep set cut to 14; class (a) withdrawn from everything but the TET10 zero-tests and `penta18.m`. Record: `tmp/ai_exchange/review_matlab_recovery.md`.
- [x] **R2a — Doxygen exclusion** (done 2026-09-16) (after: R1): `*/doc/matlab/*` added to `EXCLUDE_PATTERNS` in `Doxyfile.in`.
- [x] **R2 — Copy the keep set** (done 2026-09-16, 14 files) (after: R2a) into `src/fem/interpolation/doc/matlab/{nedelec_tet10,nedelec_tri,lagrange,facets}/` with the original file names (the scripts call each other by name; the citations are repaired in R5). Scripts that call each other sit in the same directory.
- [x] **R3 — Clean up** (done 2026-09-16; stripped-comment diff against `tmp/matlab_survey/`: 13 scripts identical in code, `fragment.m` differs by the preamble and loop only; `compare_tables.py`: 192 TET10 and 108 PENTA18 entries, 0 differences) (after: R2): `%` license header on every script (`scripts/add_license_header.py` extended to `.m`); bylines, dates and `tmp/` paths removed; the `fragment.m` preamble; the 2026-08-14 note rewritten as the node-map contract; a purpose comment at the top of every script saying what it generates, derives or checks and which C++ table or function it concerns, with the class-(b) scripts saying plainly that they are derivation notes and where they differ from the C++ (sign convention in `curl.m`, coefficient factor in `fragment.m`); the two asserting drivers and `compare_tables.py` written. Gate: `diff` of each kept script against its `tmp/matlab_survey/` original shows only comment lines (and the preamble in `fragment.m`).
- [x] **R4 — README per directory** (done 2026-09-16) (after: R3): per script the generates/derives/checks line, the C++ target, the prerequisites (Symbolic Math Toolbox, run order), how to run it, and what output means success; the class-(b) scripts labelled as derivation notes. The four directories registered in `src/fem/interpolation/doc/README.md`, with the note that they are not Doxygen API.
- [x] **R5 — Fix the citations** (done 2026-09-16; `:270` → `curl.m` with its sign caveat, `:443` → `tet10_generate.m`; no C++ edit) (after: R4): `nedelec_derivation.md:270` points at `matlab/nedelec_tri/curl.m` and says what it is (a scratch derivation with the opposite sign convention), `:443` at `matlab/nedelec_tet10/tet10_generate.m`. No C++ edit.
- [x] **R6 — Codex language sweep** (done 2026-09-16, `gpt-5.6-terra`/medium; 50 edits proposed, the substantive ones applied, the accent edits declined to keep the scripts ASCII; four suspected errors, all real, fixed by content) (after: R4) over the scripts' comments and the READMEs, with the no-fact-changes brief; apply and read each edit against the script.
- [x] **R7 — Gate** (passed 2026-09-16 under Octave 9 / `symbolic` 3.1.1 / SymPy 1.14, MATLAB not being installed: `check_tet10` exit 0, "24 value rows and 15 derivative tables match the generator"; `check_penta18` exit 0, "all 108 entries match"; `compare_tables.py` 0 differences; `defelement.m` twenty all-zero blocks; `fragment.m` (after making its preamble portable: Octave rejects `sym(name,[2 1],'real')`), `curl.m`, `check_quad.m` (prints the factor 2), `secondordertri.m`, `tet10_lagrange.m`, `parse_main.m` run without error. Verified, not only reviewed, for everything with a residual; the two facet notebooks are interactive and were not run) (after: R6): `matlab -batch check_tet10` and `matlab -batch check_penta18` from their directories exit 0 (every residual asserted zero); `python3 compare_tables.py` reports 0 differences between the MATLAB transcriptions and the live C++ tables. The class-(b) notes have no gate beyond running without error where they are runnable (`fragment.m` after its preamble, `check_quad.m`, `secondordertri.m`, `tet10_lagrange.m`); `curl.m`, `orientation.m` and `facet.m` are interactive notebooks and are not run headless.

### 4.0 Implementation Progress (updated 2026-09-16)

All steps done 2026-09-16. O1 open.

## 5. Open Design Questions

- **O1 — The interface T-matrix experiments** (`tet10/interface/main.m`, `tet10/meshtest/main.m`, `tet10/meshtest/moments.m`): a coherent, partly failed investigation of the φ↔h face projection; `interface/main.m:91-92` records "the computed T-matrix is wrong", `moments.m` reaches an orthonormalised result checked against `chol`. Nothing in `fem/maxwell/cl_Maxwell_TMatrix.cpp` uses Gram or Cholesky machinery. Default: drop all three (they stay in the transfer tree). If the projection question is still open, `moments.m` is the one to keep, together with `make_mesh.m`, `shape_functions.m`, `draw_mesh.m` and `init_points.m` which it calls. **Christian's call.**
- **O2 — `check_quad.m`**: RESOLVED 2026-09-16 by the jury → keep as a derivation note (the remap is used at `fn_IF_initialize_integration_points_on_facet.cpp:666-669`).
- **O3 — `tet10_lagrange.m`**: RESOLVED 2026-09-16 by the jury → keep as a derivation note (it does invert `V` and writes down the EXODUS pinning).
- **O4 — Runnability of `init_points.m`**: RESOLVED 2026-09-16 → moot, the magic-table cluster is dropped (no driver outside the dropped scripts, one of twelve columns, runtime test covers all).
- **O5 — Names**: RESOLVED 2026-09-16 → original names kept; the four subdirectories carry the topic, and no two kept scripts share a base name.
- **O6 — Doxygen**: RESOLVED 2026-09-16 → `EXCLUDE_PATTERNS` (R2a); the alternative, dropping `*.m` from `FILE_PATTERNS`, would also hide any future Objective-C or MATLAB that someone wants documented and is the larger change.

## 6. Layout and run order

```
src/fem/interpolation/doc/matlab/
  README.md                       what the directory is, that it is not Doxygen API, toolbox requirement
  compare_tables.py               MATLAB transcriptions vs live C++ tables (TET10 tables, PENTA18 d2NdXi2)
  nedelec_tet10/
    README.md
    tet10_generate.m              generator (run first; pins lambda = [xi; zeta; eta; tau])
    tet10_shortcuts.m             shorthand (run by the zero-tests)
    tet10_function.m              zero-test, values
    tet10_derivatives.m           zero-test, derivatives
    check_tet10.m                 driver: runs both zero-tests and asserts every residual zero
    parse_main.m                  prints mWzeta in C++ syntax
    tet10_lagrange.m              Lagrange N in EXODUS order (standalone)
    defelement.m                  DefElement cross-check (standalone)
  nedelec_tri/
    README.md
    curl.m  fragment.m  secondordertri.m  check_quad.m     derivation notes (standalone)
  lagrange/
    README.md
    penta18.m                     zero-test of cl_IF_PENTA18.hpp::d2NdXi2
    check_penta18.m               driver: runs it and asserts
  facets/
    README.md
    orientation.m  facet.m        PENTA6 orientation notebooks (interactive)
```

## 7. Definition-of-Done Checklist

- [x] `Doxyfile.in` excludes the directories (R2a).
- [x] Every kept script has a header, a purpose comment naming its C++ target, no date, no byline, no `tmp/` path; derivation notes say where they differ from the C++.
- [x] Every kept directory has a README that says how to run each script and what output means success.
- [x] No dangling script citation in `src/fem/interpolation/doc/`.
- [x] Codex sweep applied.
- [x] R7 run: both drivers exit 0, `compare_tables.py` reports 0 differences.

## 8. Audit Trail

- Exchange: `tmp/ai_exchange/review_matlab_recovery.md` (jury), `tmp/ai_exchange/matlab_recovery_sweep.md` (language sweep).
- Survey: subagent report distilled into the first §3; every load-bearing claim re-checked by Claude before dispatch (citations at `nedelec_derivation.md:270,443`; the 4×12 table at `fn_IF_initialize_integration_points_on_facet.cpp:898` equals `init_points.m:14-17` up to the `+1` for one-based indexing; all 108 `penta18.m` entries equal `cl_IF_PENTA18.hpp::d2NdXi2`; the three md5 duplicate pairs).
- Jury 2026-09-16: Codex found the TRI6 factor-of-two, the unfinished `pyra14.m`, the broken `nabla_derivative.m`, the assert-less R7, the driverless magic-table check, `hex.m`'s duplicate orientation, the mislabelled TRI3 scripts and the existing runtime verifiers; Grok found the same plus the `:270` section mismatch, the 26/30 count, the fifteen derivative tables, the 8-of-12 face rows, the invented `intpoints_penta_quad_slave` name, the missing `P10/P20/P35` consumer, the `mCurv` mismatch, the guideline misattribution, and the Doxygen `*.m` ingestion. Two Grok findings refuted (a "byte-identical" claim the plan never made; the exchange slug). Every finding re-verified by Claude against the cited lines before this rewrite.
