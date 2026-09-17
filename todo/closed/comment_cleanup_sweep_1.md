# Comment Cleanup Sweep 1

**Date:** 2026-09-16
**Purpose:** Bring the existing source comments in `src/` under `doc/commenting_guidelines.md`: correct the comments that are wrong, remove the ones that carry history, provenance, dead code or narration, triage the TODOs, and trim the Doxygen blocks that restate a name. No executable path changes. The mechanism is a sequence of comment-only commits, one per module, each gated by a comment-only checker, `make check`, and a Codex language sweep of the rewritten comments.
**Module:** `src/` (all modules except the closed cohomology core)
**AIs involved:** Claude (plan, edits), Codex (language sweep of rewritten comments), Grok + Codex (second-tier list after this sweep, see §9)
**Status:** ✅ COMPLETE (2026-09-16). Landed in commits `caf2045` (guideline, tooling, conventions), `05a4eb5` (the sweep and the license headers) and `08c306b` (the MATLAB recovery, O8). Comment-only by `scripts/check_comment_only.sh` over the whole range (266 files); `make check` passed on Armadillo and Blaze (Christian); `make doc` fixed by filing the guideline in `update_doc_index.py`. Census 31,723 → 27,837 comment lines; dated, provenance and working-record lines 128/28/3 → 0 outside string literals and the closed core; commented-out statements 110 → 11 kept by decision; every TODO owned with an exit condition; `@brief` restatements 82 → 47, all read and kept for a unit, an ownership rule or a citation; 121 CLion bylines replaced by the license header and 238 files given one (O6, O7). Rulings: banners stay (O1–O3 earlier), O4 spot-check sufficient, O5–O8 as recorded below. Residual for sweep 2 (§9): missing contracts on `send`/`receive`/`share`, `Cell<T*>` accessors and `Matrix&` arguments; `tests/` dated lines and `tests/fem/test_EdgeFunctions.cpp:32`; the PENTA6TS double include guard (code).

> **Scope guards:**
> - **Comment lines only.** No renames, no reindentation, no include reordering, no code motion. A comment change that reveals a code defect is recorded as a `Dn` row here and fixed in a separate, ordinary commit with its own gate.
> - **Banners are out of scope.** Whether to thin the `//------` population tree-wide is a formatter decision (`doc/commenting_guidelines.md` §6). This sweep removes sub-banners (`// - - -`, `! - - -`) only inside files it touches for another reason, and adds none.
> - **The closed cohomology core is excluded**: `cl_Cohomology`, `cl_Homology`, `cl_SimplicialComplex`, `cl_Chain`, `cl_Cochain`, `fn_Smith` (`.cpp` and `.hpp`). Findings there are listed in §6 for Gregory and are not edited.
> - **`tests/`, `examples/`, `nonfree/`, `doc/`, `devlog/`, `todo/` are not touched** by this sweep. Comments in `tests/` follow the same rules but are a later pass.
> - **Adding missing contracts** (ownership on every `Cell<T*>` accessor, collectiveness on every MPI entry point, column-major on every `Matrix&` argument) is sweep 2 (§9), not this one. Under the current noise the missing facts are not visible.

---

## 1. Current State

Census of `src/` on 2026-09-16 (own measurement, scratchpad script; directives `!$omp`, `!dir$`, `#pragma` and string literals excluded; "noise" is a first-word pattern on comment lines under 45 characters and is an over-approximation that needs hand review):

| module | lines | comment | noise | dead | todo | dated | provenance | `@brief` restating |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `fem/kernel` | 48 366 | 7 381 | 837 | 23 | 9 | 53 | 17 | 0 |
| `mesh` | 52 909 | 4 011 | 565 | 18 | 4 | 6 | 3 | 5 |
| `comm` | 3 679 | 600 | 229 | 0 | 0 | 1 | 0 | 0 |
| `sparse` | 18 070 | 3 152 | 161 | 0 | 1 | 28 | 3 | 0 |
| `homology` (of which closed core) | 14 249 (5 690) | 1 354 (645) | 167 (68) | 11 (1) | 12 (10) | 2 (2) | 1 (1) | 0 |
| `fem/postproc` | 2 653 | 410 | 147 | 5 | 0 | 0 | 0 | 0 |
| `io` | 5 602 | 767 | 140 | 6 | 0 | 0 | 0 | 0 |
| `fem/maxwell` | 9 813 | 1 286 | 115 | 4 | 3 | 6 | 0 | 0 |
| `fem/iwg` | 7 214 | 1 021 | 110 | 3 | 6 | 3 | 0 | 2 |
| `physics/gasmodels` | 16 094 | 1 995 | 90 | 6 | 1 | 0 | 0 | 0 |
| `fem/interpolation` | 21 659 | 1 798 | 80 | 11 | 1 | 2 | 2 | 0 |
| `physics/materials` | 22 621 | 2 095 | 46 | 5 | 0 | 18 | 1 | 75 |
| all other modules | 60 000 | 5 800 | 336 | 18 | 4 | 9 | 1 | 0 |
| **total** | **297 078** | **31 723** | **3 023** | **110** | **41** | **128** | **28** | **82** |

Worst single files by noise: `cl_FEM_DofMgr_DofData.cpp` (251), `commtools.hpp` (221), `cl_FEM_DofManager.cpp` (86), `cl_FEM_DofMgr_FieldData.cpp` (82), `cl_FEM_DofMgr_SolverData.cpp` (81), `cl_IWG.cpp` (75), `cl_SpMatrix.cpp` (65), `cl_Mesh_GmshReader.cpp` (64), `cl_Mesh.cpp` (63), `hdf5_tools.hpp` (63). Worst by history: `cl_FEM_Controller.cpp` (28 dated lines, 14 provenance lines).

**Comments that are wrong today** (found by the Astra review, re-checked against the tree on 2026-09-16 unless marked):

| Where | What it says | What the code does |
|---|---|---|
| `cl_FEM_DofMgr_DofData.cpp`, the loop after `// make list with FACE IDs` | face IDs | the loop counts `tNumCellDofs` and calls `cell_dof_id` |
| `commtools.hpp`, fifteen `// wait until … is complete` lines | a receive completes (5 sites) or a send completes (10 sites) | the `MPI_Wait` is on the request the preceding `MPI_Isend` (resp. `MPI_Irecv`) created; the word is swapped at every one of the fifteen |
| `cl_JcFunction_Database.hpp`, the header of the field-derivative routine | values outside the table window are clamped and their derivative is zero | the `mHaveSelfField && normB < mBmin` branch returns `tC + 2 tB normB`, nonzero (re-checked 2026-09-16 before B1) |
| `cl_FEM_DofMgr_EigenValues.cpp`, the symmetric-MUMPS block | symmetric input must use the lower triangle | MUMPS accepts either triangle for assembled input; the wrapper `cl_SolverMUMPS.cpp` already says so and marks the "lower" wording as a past error |

**Bottom line:** four comments lie, about 130 carry dates, about 30 name a reviewer, about 110 statements are commented out, and roughly a tenth of all comment lines narrate the next statement. The comments that carry a contract or a reason are real and must survive; the review files agree on which they are (§3).

---

## 2. Approach

The rules are `doc/commenting_guidelines.md`; this file only orders the work. Three things the order encodes:

1. **Accuracy before volume.** The four wrong comments are fixed first, in their own commit, so a reviewer judges correctness separately from deletion.
2. **Judgment before mechanics.** The dated and provenance-carrying comments (B2) are the ones most likely to hide a real rule inside a story. They are rewritten by hand, one file at a time, before any pattern-driven deletion touches those files. A useful property: **the date in the comment names the devlog that holds the experiment.** `2026-08-21` in the watchdog comment points at `devlog/dl20260821_*.md`. Before a measurement is dropped from a comment, confirm the devlog of that date records it; if not, add one line there.
3. **Mechanics last, and never blind.** Narration removal (B5) is pattern-assisted but every deletion is read. A line that matches the noise pattern and reads as a warning or a stage heading is rewritten, not removed (`// wait until send is complete` before a buffer reuse is a warning).

Rejected alternative: a scripted tree-wide deletion by pattern (the Fable review's step A for banners, and its noise list applied unattended). Rejected because the Grok review and the Astra review both found the pattern list catches stage headings and warnings, and because banners are out of scope by ruling.

---

## 3. What the Four Reviews Agree On

Read all four (`tmp/commenting/commenting_review_{fable,astra,astra2,grok}.md`) before starting a batch. Their consensus, distilled:

**Keep, do not shorten for length alone:**
- the matrix footprint explanation in `commtools.hpp` (`spacing() * n_cols`, `capacity()` trap) — reduce its second and third copies (matrix `send`, `receive`) to a one-line pointer, keep the first
- the ownership lines in `Mesh::~Mesh()` (`// facets are deleted by sideset`, `// elements are deleted by block`)
- the MUMPS conditioning-slot warning in `cl_FEM_Controller.cpp` (`THIS FUNCTION OWNS SLOT 0 ONLY`) — minus its date
- the Callaway/Debye parameter block in `debye.hpp`, the `ResistivityLaw` enum gloss and the unit annotations in `cl_Material.hpp`, the file-level contract of `cl_Material_Abundance.hpp`
- the thin-shell volume and Jacobian discussion in `cl_FEM_Calculator.hpp`
- the LAPACK contracts in `src/linalg/lapack/`, the solver knobs in `cl_SolverParameters.hpp`, the constants in `core/constants.hpp`
- in `cl_Mesh_PeriodicityFactory.cpp`: original-endpoint matching, the straight / pure-twin / half-cut cases, and why deferred ties must not be overridden — minus phase labels, review credits and dates
- in the DOF manager: the packed connectivity format, `index` versus `my_index`, separate versus combined free/fixed spaces, offset rules, self-links, collective participation
- in the controller: first-iterate timing, reset conditions, initialization order, which solver control applies where

**Remove or move:**
- statement narration (`// delete fields`, `// increment counter`, `// solve system`, `// get my id`)
- the experiment logs inside the controller watchdog, the STRUMPACK ordering notes, `powerlaws.hpp` temperature-derivative history, the `cl_SolverMUMPS.cpp` BLR and memory-cap notes — keep the rule, move the run to the devlog of its date
- the three copies of the MUMPS symmetry note (`cl_SolverMUMPS.cpp`, `cl_IWG.hpp`, `cl_FEM_DofMgr_EigenValues.cpp`) — one full explanation in `src/sparse/doc/sparse_usage_guide.md`, a short local reminder at the wrapper's `BELFEM_ERROR`, one-line pointers at the two callers
- the user-library tutorial repeated in `cl_Material_UserDefined.hpp` and `cl_MaterialFactory.hpp` — into `materials_usage_guide.md`, which already has a section for it
- the three oversized DOF-manager method headers (`SolverData::allocate_matrices`, `SolverData::populate_graph`, `DofData::reorder_dofs`, 51 / 73 / 59 lines) — shared background into `dof_manager_usage_guide.md`, a compact contract stays
- disabled implementations: the old `save_fields` / `load_fields` bodies in `cl_Mesh.cpp`, the diagnostic print in `cl_FEM_DofMgr_SolverData.cpp`, the old print loop in `cl_FEM_Tmatrix.cpp`, the Poisson smoke probe in `cl_CutProcessorManual.cpp` (not closed core), the disabled connectivity calls in `cl_Mesh.cpp`
- `example_user_material.cpp` is a tutorial by design and stays as it is

---

## 4. Ordered Steps

Each batch is one commit (or one per module where a batch spans several), under about 2 000 changed lines, and passes the gates in R0 before it is committed. Christian runs the builds.

- [x] **R0 — Tooling** (before any edit) — done 2026-09-16
  - [x] R0a `scripts/comment_census.py`: the scratchpad census made permanent. Per-module and per-file counts of comment lines, banners, sub-banners, noise candidates, commented-out statements, TODOs, dated lines, provenance lines, working-record pointers, and `@brief` restatements. Excludes directives, string literals, and the six closed-core units (counted separately). Its output table goes into every batch's devlog entry, before and after.
  - [x] R0b `scripts/check_comment_only.sh <commit-range>` (rewritten in B2d after D3; `scripts/strip_comments.py` is its stripper): for every changed `.cpp/.hpp/.h/.c`, strip comments from the old and new revision (`g++ -fpreprocessed -dD -E -P`, which removes comments without expanding macros) and diff; for every `.f90`, strip `!` comments outside string literals with a small awk and diff. Any non-empty diff fails the gate. This replaces the "identical binary" gate the jury refuted (`BELFEM_ERROR` embeds `__LINE__`).
  - [x] R0c `scripts/prompts/comment_sweep.md`, the prompt file for the per-batch Codex language sweep (`gpt-5.6-terra`, `medium`; the dense-document row of protocol §9.1): input is the batch diff, output is a list of rewrites of the *new* comment text only, with the guideline's §10 rules and the "may not change" list from `CLAUDE.md` §"Prose Gets a Language Sweep".

- [x] **B1 — Correct the four wrong comments** (after: R0) — done 2026-09-16, own commit pending. The `cl_JcFunction_Database.hpp` claim re-verified (the `mHaveSelfField && normB < mBmin` branch returns `tC + 2 tB normB`); the header now names the self-field bridge as the exception. The `commtools.hpp` send/receive mislabel turned out to be **fifteen** lines, not one: five `wait until receive is complete` after an `MPI_Isend` and ten `wait until send is complete` after an `MPI_Irecv`; every wait comment in the file was matched to the request it waits on and corrected. `cl_FEM_DofMgr_DofData.cpp` cell-dof heading and the `cl_FEM_DofMgr_EigenValues.cpp` triangle sentence corrected; the dated measurement in the latter is B2d's.

- [ ] **B2 — History, provenance, working-record pointers** (after: B1). Hand rewrite, file by file, largest first:
  - [x] B2a `cl_FEM_Controller.cpp` (28 dated, 14 provenance) and `cl_FEM_Controller.hpp` (5 dated). Highest judgment: the watchdog, the diagnostic arming and capture, the thermal-kernel attachment, the trust-streak member. The guideline's §8 watchdog rewrite is the model. **First half done 2026-09-16** (Christian: "start with the watchdog and diagnostics blocks"): watchdog spare rule and its two header comments, floor-escalation member, conditioning arm/capture/finalize/footer/parse blocks — 17 rewrites, comment-only by the checker; the unlogged 2026-08-21 experiment retained in the session devlog (O3 applied). **Second half done the same day**: the 26 remaining sites, 25 replacements, three unlogged facts retained in the session devlog (INFOG −10 quench, the 19-rejection night, the lagged-reference fix).
  - [x] B2b `sparse/` — done 2026-09-16: 32 sites in eight files (plus `cl_SpMatrix.cpp`, `mumpstools.f90`), the MUMPS symmetry note consolidated into `sparse_usage_guide.md` with pointers from the wrapper, `cl_IWG.hpp` and the eigen caller; two unlogged facts retained in the session devlog.
  - [x] B2c `physics/materials` — done 2026-09-16: 30 sites in 15 files (`powerlaws.hpp` 11, the six decision lines across Chromium/Copper/Silver/WhiteTin/Lead/Nickel, the angle-contract notes, the reviewer asides); five facts retained in the session devlog.
  - [x] B2d `fem/kernel` remainder — done 2026-09-16: 24 sites in eleven files (also `cl_FEM_DofMgr_SolverData.cpp`, `cl_FEM_Postprocessor.cpp`, `cl_ThinShellFactory.hpp`, `core/assert.hpp`, two controller sites a quote character had hidden from the first grep); one fact retained in the session devlog.
  - [x] B2f the modules B2 did not enumerate — done 2026-09-16, 24 sites in 20 files; three facts retained in the session devlog (found by the tree-wide listing on 2026-09-16, about 30 sites): `circuit/cl_ElectricalCircuitFactory.cpp`, `comm/cl_Communicator.{cpp,hpp}`, `core/stringtools.cpp`, `math/graph/graphtools.hpp`, `math/tools/fn_{cardano,ferrari}.hpp` (audit credits naming Claude, Codex and Grok), `fem/interpolation/nedelec/cl_EF_{HEX8,TET10}.cpp`, `fem/iwg/cl_IWG.hpp`, `fem/iwg/cl_IWG_Timestep.hpp` (2), `fem/maxwell/` (`cl_MaxwellBoundaryConditionFactory.cpp`, `cl_MaxwellPostprocessor.cpp`, `cl_MaxwellFactory.cpp` 3, `fn_mesh_config_tag.hpp`, `cl_IWG_Maxwell.cpp`, `matrices/mt_maxwell_h.cpp` with an exchange-thread pointer), `fem/thermal/cl_ThermalBoundaryConditionFactory.cpp`, `numerics/integration/fn_intpoints_gauss_{pyra8,pyra27,pyra64,pyra125,tet35,tet56}.hpp` ("Regenerated 2026-09-03" — six copies of one sentence).
  - [x] B2e `mesh/` — done 2026-09-16: 9 sites; the PENTA6TS CLion byline left for O6.

- [x] **B3 — Commented-out code** (after: B1) — done 2026-09-16; see the devlog for the file list and the age check. The 110 statements plus the block-disabled bodies listed in §3. `git log -S` on anything that looks younger than a month before deleting it. A block that is the naive twin of a tuned kernel becomes a test case, not a comment (`doc/commenting_guidelines.md` §9).

- [x] **B4 — TODO triage** (after: B3) — done 2026-09-16; per-line outcomes in the devlog. 41 lines, of which 10 are in the closed core (§6). For each of the other 31: rewrite as `TODO(owner): … ; remove when …`, or move to a `todo/*.md` file leaving a one-line pointer, or delete as stale. Record each outcome in the devlog. Known clusters: `cl_IWG.hpp` obsolete overloads, `cl_IWG.cpp` "get rid of these matrices", `cl_Maxwell_TMatrix.cpp` "these lines can be deleted", `cl_FEM_DofMgr_FieldData.hpp` "might not need this", `cl_ThermalFactory.cpp` axisymmetry, `cl_FEM_Controller.hpp` "reset omega" beside a duplicate `// end user settings`, `cl_FEM_DofMgr_DofData.cpp` periodic / thin-shell TODOs.

- [ ] **B5 — Narration**, one commit per module, pattern-assisted, every line read (after: B2 for the files B2 touched):
  - [x] B5a `comm/commtools.hpp` (221) — done 2026-09-16, 418 → 81 comment lines, Codex check passed (pointer reworded): the narration is duplicated across every `send`/`receive` specialization; keep a comment in the generic specialization only if it survives the guideline's §2 test, delete the copies; the `wait until send is complete` before a buffer reuse is a warning and is rewritten once at the reuse
  - [x] B5b `fem/kernel/cl_FEM_DofMgr_DofData.cpp` — done 2026-09-16, 691 → 171 comment lines, all 83 sub-banners gone, Codex check applied (one deleted reason restored, four survivors removed)
  - [x] B5c rest of `fem/kernel` — done 2026-09-16 (1,829 → 895 comment lines over seven files); Codex check applied (seven amputated continuations restored — D5 — five survivors removed, four additions reworded): `cl_FEM_DofManager.cpp`, `cl_FEM_DofMgr_SolverData.cpp`, `cl_FEM_DofMgr_FieldData.cpp`, `cl_FEM_Calculator.cpp`, `cl_FEM_Element.cpp`, `cl_ThinShellFactory.cpp`, `cl_FEM_Kernel.cpp` — includes the three long method headers of §3
  - [x] B5d `mesh/` — done 2026-09-16 (964 → 352 comment lines over six files); Codex check applied (two reasons restored, five survivors removed): `cl_Mesh.cpp` (destructor first: it is the ownership template for every owning type), `cl_Mesh_GmshReader.cpp`, `cl_Mesh_OrderConverter.cpp`, `cl_ProtoMesh.cpp`, `cl_FaceFactory.cpp`, `cl_Mesh_VtkWriter.cpp`
  - [◐] B5e `sparse/cl_SpMatrix.cpp`, `io/hdf5_tools.hpp`, `io/cl_HDF5.cpp`, `io/cl_XML.cpp`, `io/cl_Input_Section.cpp` — done 2026-09-16 (643 → 351 comment lines), Codex check passed (two wrong labels removed)
  - [x] B5f `fem/postproc/fn_Mesh_compute_surface_normals.cpp` (143 → 20), `fem/iwg/cl_IWG.cpp` (247 → 121), `fem/maxwell/cl_MaxwellFactory.cpp` (503 → 416), `physics/gasmodels/cl_Gas.cpp` (534 → 330; the entropy-spline-mode explanation kept), `homology/cl_CutFactory.cpp` (222 → 100) — done 2026-09-16, Codex check applied
  - [x] B5g everything else with a noise count above 20 in the census — done 2026-09-16 (eight files; the controller's audit-round labels and the circuit factory's remaining decision IDs went with it), Codex check applied

- [x] **B6 — Doxygen triage** (after: B5), headers only (B6a, B6b and the Codex check done 2026-09-16; B6c open, Christian's):
  - [x] B6a `physics/materials` — done 2026-09-16: `cl_Material.hpp` (29 one-line `@brief Get …` / `Set flag` / `Clear flag` blocks), `cl_MaterialFactory.hpp`, `cl_JcFunction.hpp`, `cl_Material_Metal.hpp`, `cl_Material_YBCO.hpp`; the tutorial relocation of §3; `cl_Material_Abundance.hpp` untouched
  - [x] B6b (done 2026-09-16; `cl_IWG.hpp` and `cl_Mesh.hpp` turned out to carry no `@brief` blocks at all) `fem/iwg/cl_IWG.hpp`, `mesh/cl_Mesh.hpp`, `mesh/cl_Mesh_ConnectivityCalculator.hpp`, and `comm/commtools.hpp` (Codex's B5a check listed 22 Doxygen blocks there that restate their signatures — the utility declarations near the top and the per-overload `send`/`receive`/`broadcast` blocks; the raw-array receive block keeps its preallocation and release-error contract) (the `@brief Destructor (…)` form becomes the parenthetical alone)
  - [ ] B6c run `make doc` once after B6 and spot-check that the generated pages for `Cell`, `Mesh`, `Material` still carry every contract that was kept

- [ ] **B7 — Sub-banners inside touched files** (folded into B2 to B6, no separate commit): `// - - -` and `! - - -` lines in a file being edited go; untouched files keep theirs.

- [ ] **R9 — Closing gate** (after: B6; `make check` passed on Armadillo and Blaze 2026-09-16, `make doc` fixed the same day by filing the page in `update_doc_index.py`; commits and the range gate open): full `make check` on both a debug and a release tree (Christian); `scripts/comment_census.py` final table into the devlog; `scripts/check_comment_only.sh` over the whole range of sweep commits; the Codex language sweep of every batch applied and recorded; this file moved to `todo/closed/` with the Status rewritten as a summary.

### 4.0 Implementation Progress (updated 2026-09-16)

**Implemented and reviewed:** R0a–R0c, B1. Gates that ran: `scripts/check_comment_only.sh` on the B1 working-tree diff (four files, all comment-only); a deliberate `k = 0 → k = 1` in `cl_FEM_DofMgr_DofData.cpp` made it fail (exit 1) and the restore made it pass; the same pair on `mumpstools.f90` (a comment-case change passed, `gMaxNumSolvers = 8 → 9` failed). `scripts/check_doc_claims.py` 38/38 after the `CLAUDE.md` pointer landed. Census unchanged by B1 (it rewrites, it does not remove). No build run; the change is comment-only by the checker, and `make check` is Christian's gate at R9.

**Defects surfaced:** none that are code. D1 (LOW, documentation): the `commtools.hpp` mislabel count in the Astra review (one site) understated the tree (fifteen); recorded above. D5 (MEDIUM, method, found by the Codex check of B5c): deleting narration by exact text match of a short generic line (`system matrix`, `duplicate`, `identity`, `master and slave`, `receive_submesh`, `contains` + `number of elements` + `element id` + `number of dofs` + `dof index`) amputated the continuation lines of seven multi-line comments — two of them payload-layout contracts in `cl_FEM_DofMgr_SolverData.cpp` — leaving the surviving first line truncated mid-sentence. All seven restored 2026-09-16 (the two layouts rewritten as one two-line contract each). Rule for every later narration pass: a candidate line is deleted only if the previous line is not a comment line, or is a comment line that ends a sentence (`.`, `)`, `:`, `!`, `?`) or is a banner; a heuristic scan over all B5 files for the signature (a deleted comment line whose surviving predecessor ends mid-sentence) found no further case — but the Codex check of B5d found one where the deleted line was the *first* line and the survivor the continuation (`cl_ProtoMesh.cpp`), so the guard and the scan must look in both directions. Second lesson from B5f: the guard protected `// orient_terminal_curves_sub()` and `// duplicate_and_relink_facets` in `cl_CutFactory.cpp`, I overrode it as "plain narration" without reading the context, and both were the continuations the guard had seen. **A guard hit is cleared only after reading the lines around it, never by re-running with the guard off.** D4 (LOW, process, found in B3): the first B3 sweep prompt carried an empty diff (`git diff src/ -- ':!…'` is a revision error; the correct form is `git diff -- src/ ':!…'`) and Codex answered anyway, mirroring the prompt's own examples. Fixed the same session by rerunning with the real diff; the per-batch prompt now has to be checked for a non-empty diff block before dispatch (a `grep -c '^[-+]' <diff>` line in the dispatch command). D3 (HIGH, tooling, found in B2d): the first `check_comment_only.sh` stripped C/C++ comments with `g++ -fpreprocessed -dD -E -P`, which rejects a `#` on a macro continuation line (`#aCheck` in `core/assert.hpp`) as an invalid directive; and because the script ran under `set -e`, that g++ failure aborted the run after the previous file with no summary line, so a reader looking only at the last line would have seen `comment-only: src/core/assert.cpp` and taken it for a pass. Neither earlier verdict is affected: `assert.hpp` entered the diff only with B2d, and every earlier run printed its `OK (… files)` summary. Fixed 2026-09-16: comments are now stripped by `scripts/strip_comments.py`, a tokenizer that copies string, character and raw-string literals and preprocessor text through and removes `//`, `/* */` and Fortran `!` comments (keeping `!$` and `!dir$` directives); the checker no longer uses `set -e`, reports a file it cannot strip as a failure, counts the files it checked, and always prints the summary. Re-gated: 40 changed files comment-only; a `#define` value change in `assert.hpp` and a constant change in `mumpstools.f90` each fail it; unit input with `#a` in a macro, `//` inside a string, a raw string and a multi-line block comment strips correctly. D2 (LOW, documentation, found by the Codex sweep of B1): my first B1 wording for the `cl_JcFunction_Database.hpp` header said the self-field tangent is "not zero", but `deval_dB()` returns 0 where `eval()` returns the constant table edge (`tJ0 < tV`); the header now says the tangent is the bridge's own and zero only where the bridge is idle. Fixed 2026-09-16, comment-only gate re-run.

**Codex sweep of B1** (`gpt-5.6-terra`/medium): three edits and one suspected error. The error was real (D2). The eigenvalue-note rewording was applied. The other two edits ask to *delete* the fifteen corrected `wait until …` lines and the corrected `make list with CELL IDs` heading as narration — right by the guideline, but deletion is B5a/B5b, not B1, which corrects only; carried forward to those batches.

---

## 5. Open Design Questions

- **O1 — Banner thinning.** Out of scope by ruling. Left open here only so the census keeps counting banners (about 9 800) for whoever takes the formatter decision.
- **O2 — What the census counts as noise.** The first-word pattern flags `// xi` and stage headings. Options: (a) keep the over-approximation and rely on the hand read, (b) exclude lines that are a single mathematical symbol or that sit directly above a loop opener. Proposal: (a), because (b) hides the mumbling variant of the same line.
- **O3 — Where a rewritten measurement goes when its devlog has no line for it.** Proposal: one line appended to the devlog of the comment's date under a `## Retained from source comments` heading; if no devlog carries that date, the current session's devlog.
- **O4 — Doxygen output.** RULED 2026-09-16: the B6c spot-check is sufficient. `make doc` is live. Does the generated reference need a before/after diff, or is B6c's spot-check enough? Christian's call.
- [x] **O6 — CLion bylines.** RULED 2026-09-16: replace every byline with the license header. DONE the same day: 121 blocks removed by `scripts/add_license_header.py`, gated. 51 files in `src/` open with `// Created by <name> on <date>.`, the IDE template; two of them named `claude` and were removed in B2d. The guideline bans bylines and dates, and the license block already names the developers. Proposal: strip all 51 in one mechanical commit (comment-only, checker-gated). Christian's call, because 49 of them carry his own name.
- [x] **O7 — Files without a license block.** RULED 2026-09-16: every source file carries the header, python files in `examples/` and elsewhere included. DONE the same day: 238 headers added (90 `.cpp`, 46 `.hpp`, 106 `.py`); `.sh`/`.cmake`/`CMakeLists.txt` not asked and not touched. 89 `.cpp`/`.hpp` files in `src/` carry no `Copyright (c)` header (the CLion-byline files among them). Adding the block is a comment-only change but it is not commenting cleanup; it belongs to a release checklist, not to this sweep. Logged so it is not lost.
- [x] **O8 — The TET10 table generator lives in `tmp/`.** RULED and DONE 2026-09-16 → `todo/matlab_derivations_recovery.md`: the generator, its zero-tests and twelve more scripts now live under `src/fem/interpolation/doc/matlab/`, verified under Octave; `cl_EF_TET10.cpp` keeps its "generated symbolically" sentence without a path (no C++ edit). `cl_EF_TET10.cpp` pointed at `tmp/tet10/tet10_function.m` and `tet10_derivatives.m` as the source of truth for its shape-function tables; the pointer is gone (working-record path), but the generator is the only thing that can regenerate the tables and `tmp/` is swept. Proposal: move the two MATLAB files under `scripts/` or `src/fem/interpolation/doc/` and point the comment at them. Christian's call.
- **O5 — Dead-code count.** The Fable review counted 56 commented-out statements on 2026-09-16 with a near-identical regex; this census counts 110. The difference is not resolved; B3 works from the census listing and reports the true number.

---

## 6. Findings in the Closed Cohomology Core (for Gregory, not edited)

`fn_Smith.hpp` and `fn_Smith.cpp` carry eight copies of `// todo: Optimize this to put directly in the while (cont) statement…` (four each) and two further `// todo:` lines; `cl_SimplicialComplex.cpp` has 46 noise candidates and one commented-out statement; `cl_Cohomology.cpp` has two `// todo:` lines. None of these is touched by this sweep. Protocol §7.1 treats comment-only edits there as banned; the list is here so Gregory can decide.

---

## 7. Definition of Done

- [◐] R0a, R0b, R0c exist (2026-09-16); used by B1; every later batch must use them.
- [ ] B1 landed in its own commit; the four rows of §1 are each verified against the tree in the devlog.
- [ ] Every batch's devlog entry carries the before/after census table and the Codex sweep verdict.
- [ ] `scripts/check_comment_only.sh` passes over the full commit range.
- [ ] `make check` green on the release and the debug tree after R9 (Christian).
- [ ] Census after R9: dated lines 0 outside string literals; provenance lines 0; working-record pointers 0; commented-out statements 0 outside labeled explanatory fragments; every TODO owned or moved; `@brief` restatements 0 in the B6 headers.
- [ ] §6 handed to Gregory in a devlog line.
- [x] `CLAUDE.md` carries a pointer to `doc/commenting_guidelines.md` and its one-paragraph summary (landed with B1, 2026-09-16; `check_doc_claims.py` 38/38).

---

## 8. Audit Trail

- Guideline jury: `tmp/ai_exchange/review_commenting_guidelines.md` (Codex `gpt-6-astra`/high, Grok `grok-4.6`/high, 2026-09-16; swept per protocol §10 after the devlog distilled it).
- Source reviews this plan is built from: `tmp/commenting/commenting_review_fable.md` (order of work, census script, noise list), `commenting_review_astra.md` (the reading-cost criterion, the Cell and spline examples), `commenting_review_astra2.md` (the four wrong comments, the ranked targets, the keep list), `commenting_review_grok.md` (P0/P1/P2, the Mesh destructor template, the materials `@brief` pass, the closed-core caveat). Every file:line those reviews cite was pinned to `8465284`; the lines have moved since and this plan cites tokens, not lines.
- Devlog: `devlog/dl20260916_commenting_guidelines.md`.

---

## 9. After This Sweep (steps 4 and 5 of the campaign)

When R9 is done, Grok and Codex each produce an independent second-tier list over the swept tree: what is still narration, what contracts are now visibly missing (ownership on `Cell<T*>` accessors, collectiveness on MPI entry points, column-major on `Matrix&` arguments, status semantics on solver wrappers, scratch-buffer notes), and what the first sweep got wrong. The two lists are reconciled into `todo/comment_cleanup_sweep_2.md`, which adds the missing contracts and repeats the Codex language sweep. That plan is not written until this one closes.

Carried over to sweep 2 from this one: `tests/` was outside this plan's scope and still carries dated lines (`tests/physics/backendfree/test_MaterialBackendFree.cpp:18`, "the 2026-07-02 SplineLookupTable refactor"); `broadcast`/`collect` now state their collectiveness (B6 Codex check) but `send`/`receive`/`share` do not; the double include guard in `cl_Element_PENTA6TS.hpp` (both `CL_ELEMENT_PENTA6TS_HPP` and `BELFEM_CL_ELEMENT_PENTA6TS_HPP`) is code, not comment, and is only noted here.
