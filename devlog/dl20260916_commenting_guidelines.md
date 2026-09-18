# Devlog 2026-09-16 — Commenting guidelines, step 1: consolidation

**Date:** 2026-09-16
**Author:** Claude Code (Fable 5.1)
**Topic:** Consolidate four AI-drafted commenting guidelines and the 2026-09-15 language agreements into `doc/commenting_guidelines.md`
**Status:** steps 1–3 done; plan approved by Christian; R0 (tooling) and B1 (the four wrong comments) done the same day, uncommitted; Codex sweep of B1 pending at the time of writing

## Summary

Christian opened the session with two examples of comments that hurt rather than help: a
seven-line history of a removed setter above `gMaxNumSolvers` in `mumpstools.f90`, and the
dated comments in `cl_SolverMUMPS.cpp`. He had four drafts of a commenting guideline
(Astra, Fable web, Gemini, Grok) in `tmp/commenting/`, a note of the 2026-09-15 language
agreements (`tmp/commenting/language.md`), a *Clean Code* "Comments" chapter extract in `tmp/martin/` (first edition, where it is Ch. 4; the second edition, where it is Ch. 5, entered `literature/books/coding/` the same day and the guideline cites that one — see `dl20260916_clean_code_second_edition.md`),
and the two scientific-software books in `literature/books/coding/`. The five-step plan:

1. consolidate the drafts into `doc/` (this session)
2. jury round on the consolidated document
3. from the four `commenting_review_*.md` findings, plan a first cleanup sweep
4. second sweep, with a Grok + Codex second-tier list
5. repeat

## What was written

`doc/commenting_guidelines.md`, registered in `doc/README.md` under "Guidelines and
Philosophy". Fourteen sections: the one rule; what a comment may say (ten kinds, each with
a greppable anchor in the tree); what it must not say (thirteen anti-patterns with the
census figure for each); placement and length; Doxygen; banners; TODOs; no history, no
provenance, no dates (with two worked rewrites, the `mumpstools.f90` constant and the
controller watchdog); reference implementations; language; extract-vs-comment; comments
edited with the code; a review checklist; relation to the other conventions; sources.

## Decisions taken in the consolidation

Where the four drafts disagreed, the document picks a side and the jury can overturn it:

| Question | Drafts | Chosen |
|---|---|---|
| Banners | Fable: script away ~6 500 lines; Grok/Gemini: house style, leave; Astra: out of scope | House style stays, one per definition, none in bodies, no new sub-banners; tree-wide thinning is a separate formatter decision |
| Reference (unoptimized) code in comments | Fable/Grok: allow a labeled block (Oliveira & Stewart Ch. 12); Astra: tests are the home | Tests are the home; a labeled formula fragment is fine as explanation, an old implementation is not |
| Repeating a warning at two sites | Fable/Grok: once only; Astra: a concise copy at each independently editable site can be justified | Once at the site of the mistake; a one-line pointer at the second site |
| Measurements in a comment | Astra2: benchmark numbers stay where they are the evidence; language.md: no dates | The rule plus the one fact it depends on stays; the deck, date, trace and cascade go to the devlog |
| Doxygen `@brief` on public classes | Grok: one `@brief` per public class; Fable: only when a contract | Class-level `@brief` welcome; member-level only when it states a contract |

Two rules the drafts did not have, added from the tree:

- **Review provenance is banned** (`( Codex+Grok audit )`, `( Grok R1c )`, `Created by claude`):
  31 lines in `src/`. A concurrence is the weakest evidence rung and says nothing about why.
- **Working-record pointers are banned** (`tmp/ai_exchange/…`, `INC-`, `DR-`, `L-NN`): 3 lines
  in `src/`. Extends the existing `doc/`-must-not-cite-devlogs rule to source.

## Census (own measurement, `src/`, `*.cpp` + `*.hpp` unless stated)

| quantity | count |
|---|---|
| source files (`.cpp .hpp .h .c .f90`) | 951 |
| `//` lines | 30 737 |
| banner lines (`//----`, ≥ 8 dashes) | 9 734 |
| sub-banner lines (`// - - -`) | 499 |
| dated comments (`YYYY-MM-DD`, incl. `.f90`) | 128, of which 28 in `cl_FEM_Controller.cpp` |
| review-provenance lines (`codex|grok|claude`) | 31 |
| `tmp/ai_exchange` + `INC-/DR-/L-` refs | 3 |
| `TODO`/`FIXME` (incl. `.f90`) | 44 |
| commented-out statements (regex) | 110 (the Fable review counted 56 with a similar pattern; the difference is not resolved) |
| `@brief` in headers / restating the name | 412 / 88 |
| `// xi` | 24 |

The Fable review's figures (30 689 `//` lines, 9 702 banners) match within the three days
of drift.

## Verification of anchors

Every file and token named in the guideline was re-opened this session. Three of my first
draft's claims were wrong and corrected before the file was registered: the controller
anchor was quoted across a line break (now anchored on `( leaked across retries )`);
`arpacktools.f90` does not use `!>` blocks (only `mumpstools.f90` does); the disabled
`save_fields`/`load_fields` bodies are `/* … */` blocks, not `//` lines, so the line-count
regex does not see them.

One citation not verified against a primary source: the watchdog paper (Chamberlain,
Powell, Lemaréchal, Pedersen 1982) named in the §8 rewrite. Confidence medium; the jury
should check it or the example should drop the citation.

## Jury round 1 (step 2)

Blind `--jury` on the document, Codex `gpt-6-astra`/high (Christian's choice; the tag was added
to the allowlists of `.claude/scripts/ask_codex.sh` and `scripts/cross_review.sh` this session)
and Grok `grok-4.6`/high. Record: `tmp/ai_exchange/review_commenting_guidelines.md`. No auditor
wrote to the tree (`git status` after the round). Of 24 findings, 21 confirmed against the tree,
one partially refuted, one refuted on a token, one routed to Christian.

Both auditors independently: the §8 watchdog exemplar was itself stale (`mWatchdogWindow`
defaults to 30, the "eight" was the experiment) — P0 in a document about stale comments; the
"binaries identical under `-O2`" sweep gate is unsatisfiable because `BELFEM_ERROR` embeds
`__LINE__`; the closed cohomology core (eight of the 44 TODOs are in `fn_Smith`) had no
exclusion; the checklist contradicted the kinds table on rejected alternatives, names and todo
pointers. Codex alone: the deletion test and "reviewer deletes stale comment" could erase an
ownership or MPI contract; `!$omp` directives need a directive exclusion. Grok alone: the
devlog-citation rule was misattributed to `documentation_guidelines.md` (it is `CLAUDE.md`);
two controller error strings carry dates and a regex sweep must skip string literals; the
Chamberlain 1982 citation was unverifiable locally (dropped, per the pre-registered rule);
`// xi` in a Nédélec kernel is a symbol annotation, not mumbling. Neither auditor disputed the
two choices I had flagged as contested (reference implementations in tests; banners kept).

Christian's ruling on the routed question: `Decision: <name>, <date>` lines in source lose both
the name and the date — he takes the blame by default. Recorded in §8.

Revision applied the same session (295 → 310 lines); the resolution entry in the exchange file
lists every change.

## Step 3 — the sweep plan

`todo/comment_cleanup_sweep_1.md`, registered in `todo/README.md`. Built from the four review
files and a per-module census (scratchpad script; directives and string literals excluded; the
closed core counted separately). Order: tooling (permanent census script, a comment-only checker
built on `g++ -fpreprocessed -dD -E -P`, a per-batch Codex prompt); B1 the four wrong comments in
their own commit; B2 dated/provenance/working-record comments rewritten by hand, largest file
first (`cl_FEM_Controller.cpp`: 28 dated, 14 provenance); B3 commented-out code; B4 TODO triage
(31; 10 in the closed core listed for Gregory); B5 narration per module; B6 Doxygen triage; R9
closing gate. Five open questions, banners excluded by ruling, sweep 2 deferred until this closes.
Census totals: 297 078 lines, 31 723 comment lines, 3 023 noise candidates, 110 commented-out
statements, 41 TODOs, 128 dated lines, 28 provenance lines, 82 `@brief` restatements.

## R0 and B1 (plan approved, same day)

Christian approved the plan and asked for R0 and B1. Landed:

- `scripts/add_license_header.py` (new; O6/O7)
- `scripts/comment_census.py` — the scratchpad census made permanent (`--files`, `--closed`, `--root`; directives and string literals excluded; closed core in its own row).
- `scripts/check_comment_only.sh` — strips comments from both revisions of every changed source (`g++ -fpreprocessed -dD -E -P` for C/C++, an awk that respects string literals and keeps `!$` directives for Fortran) and diffs. Controls: a `k = 0 → k = 1` injection in `cl_FEM_DofMgr_DofData.cpp` fails it, restore passes; `gMaxNumSolvers = 8 → 9` in `mumpstools.f90` fails it, a comment-case change passes. My first Fortran negative control was a no-op (sed on the wrong line, nothing changed, "pass" proved nothing) and was redone with the diff line count checked first.
- `scripts/prompts/comment_sweep.md` — the per-batch Codex prompt.
- B1, four files, comment-only by the checker: `cl_FEM_DofMgr_DofData.cpp` cell-dof heading; `cl_JcFunction_Database.hpp` derivative header now names the self-field bridge as the exception to "tangent is zero outside the window" (re-verified: the `mHaveSelfField && normB < mBmin` branch returns `tC + 2 tB normB`); `cl_FEM_DofMgr_EigenValues.cpp` "lower triangle" → "each symmetric pair once, either triangle"; `commtools.hpp` **fifteen** wait comments (5 said receive after an `MPI_Isend`, 10 said send after an `MPI_Irecv`), every wait in the file matched to its request with an awk pass.
- `CLAUDE.md` gained a "Comments" subsection under Coding Standards with the pointer; `check_doc_claims.py` 38/38.

Codex sweep of the B1 diff (`gpt-5.6-terra`/medium): the eigenvalue-note prose applied; two deletion suggestions (the corrected wait lines and the cell-dof heading are narration) deferred to B5, since B1 corrects and does not delete; and one real catch — my "tangent is the bridge's own, not zero" was wrong where `eval()` returns the constant table edge (`tJ0 < tV`, tangent 0). Qualified, gate re-run. A language sweep of four comment lines caught a factual error the jury-style reading of the code had not: worth keeping the per-batch sweep even for small batches.

Reviewed, not verified: no build ran. Nothing committed; the B1 files, the three scripts, the `CLAUDE.md` edit and the document set are all in the working tree.

## Library update: Martin 2025 and Thomas & Hunt 2020

Christian added the full text of *Clean Code* (**second** edition, 2025) and *The Pragmatic
Programmer* (20th anniversary edition, 2020) to `literature/books/coding/`, with navigation files.
Two consequences: the extract I had used in `tmp/martin/` was the first edition's Ch. 4, and the
library's index already states that the guideline cites the second edition — so the guideline was
re-cited to Martin 2025 Ch. 5 with printed-page anchors, and its closing note ("the Martin chapter is
not in the library") corrected. The second edition reverses the first on TODO comments ("TODO means
Don't Do", p. 106); §7 now states BELFEM's owned, exit-conditioned TODO as a deliberate narrow
exception to that position, with the reason (`todo/` is the backlog). The Appendix debate with
Ousterhout is cited in §11 for the one sentence both sides agree on. Thomas & Hunt are cited in
§1 (Topic 9 DRY, Tip 13 "why not how"), §2 (Topic 23, contracts), §12 (Topics 23, 25, contract as
comment then as assertion). `doc/literature_references.md` gained a "Software Craft" section with
all four books (ISBNs from the extracts' front matter, DOIs for the two Cambridge titles), since
the guideline's §4 sends readers there and it listed none of them. `CLAUDE.md` literature routing
gained two rows. Every quotation was re-read in the extracts this session. While this landed, another session had already updated the guideline's Sources entry for Martin to the second edition and swapped the §1 quotation to the 2025 wording; three of my replacements therefore found no match and were redone against the current text (`git status` before the edit, as L-05 asks) — the two edits are compatible and both stand.

## B2a, first half: watchdog and diagnostics blocks

Sixteen comment blocks rewritten in `cl_FEM_Controller.cpp` and `.hpp` (the watchdog's four, the conditioning arm/capture/footer/parse blocks); every date, review credit and `INC-` pointer in them removed, the rule and its one supporting fact kept. Comment-only by `scripts/check_comment_only.sh`. The 2026-08-29/30 history those blocks cited is recorded in `devlog/dl20260830_mumps_error_analysis_split.md` (the slot-0 NaN overwrite, the "scaled matrix" correction, the ICNTL(11) refinement refutation, the 2.96e7 vs 1e4-1e5 measurement, `INC-213`), checked before the dates left the source. The 2026-08-21 watchdog experiment had **no** devlog, so it is retained below (plan O3). Codex sweep of this half (`gpt-5.6-terra`/medium): eight language edits applied (one wording kept: the guard *spares* a cut, it does not prevent cuts); two guideline violations in my own rewrites, both right — "this comment exists because that overwrite has happened" was a past-incident record and became the rule ("do not re-read them here"), and the 2.96e7-vs-1e4 figure was an experiment log already held by `dl20260830_mumps_error_analysis_split.md` and became "can differ by orders of magnitude"; one suspected error, right — |λ|max/|λ|min is a spectral property, not "a 2-norm property" (the induced 2-norm uses singular values and equals the ratio only for a normal matrix, which the same comment says two lines later), corrected. Re-gated comment-only.

### Retained from source comments

- **Watchdog spare rule, measured 2026-08-21 on the coarse `tapestack3d` deck** (formerly in `Controller::watchdog_magnetic`): the AIMD recovery of the line-search relaxation from a backtracked omega of about 0.1 takes 17–22 iterates at beta = 1.1; a no-new-best window of W = 8 cut that recovery mid-climb and the cut cascaded (each halving improved the best residual and was cut again, 2.10 → 1.05 → 0.53 ms). A genuine stall showed omega pinned or shrinking (the hd-floor grinds; one case sat at 383 iterates at 2.10 ms with omega frozen at 0.053), and the spare rule still fires there. A three-AI round the same day (Grok R3) found that every return path must refresh the previous-omega memory after the test, or one spared iterate blinds the detector for the rest of the attempt. The default window today is 30.

## B2a, second half: the remaining 26 sites

The user-settings members in the header (absolute tolerance, update gate, the live absolute
escape), the soft-fail arming, the kappa formatter, the line-search and promotion notes, the
lagged-reference fix, the flat-reject and thermal-freeze counters, the trust streak, the retry
latches, the rejection message, the thermal eigen driver, the PETSc tolerance gate, the chi and
ghost opt-in notes, the penalty log rows, the restart history levels and the BDF restart state.
25 replacements, every reviewer credit, date and plan-step ID (`MIT-3a`, `DQ1`, `EQ1`, `RQ2`,
`R1c`, `SQ5`) removed; the controller now carries **zero** dated or credited comment lines
(from 47). Comment-only by the checker. Devlog coverage checked for each date; the facts below
had none and are retained here.

### Retained from source comments (second half)

- **MUMPS INFOG(1) = −10 at quench, 2026-07-27** (formerly in `Controller::initialize`): the
  singular factorization at a strained iterate is what motivated the soft-fail contract — under
  the controller it is a failed trial that cuts the timestep, not a run killer.
- **Quench night 2026-08-13, 19 rejections** (formerly above the rejection message in
  `reset_timestep`): the rejection print used to be commented out, so every one of the 19
  rejected steps had to be reconstructed after the fact by diffing step headers; that night is why
  the message exists. On its first live print the message named step 683 for an attempt headed
  684, which is where the `+1` on the counter comes from.
- **Lagged-reference fix, 2026-08-10** (formerly in the line-search baseline comment): `mEpsilon0`
  measured under the previous assembly rejected all eight Newton trials flat in omega whenever
  the two assemblies disagreed by more than the acceptance band; the fix replaces it by the
  pre-update residual of the current assembly the moment it is known.

Codex sweep of the second half (`gpt-5.6-terra`/medium): 14 language edits applied, one reworded (the flat-pair exception read ambiguously as written by Codex); the function name in the latch comment taken from the enclosing definition, not from the sweep; two guideline violations in my own rewrites, both right — the PETSc tolerance-gate block still told the tapestack3d story with the observed ~1e-15 (the mechanism stays, the deck and the number go; `dl20260818_reltol_default_strumpack.md` holds them), and "both have changed their default before" was history. No suspected errors. Re-gated comment-only. **B2a closed: the controller carries zero dated, credited or record-pointing comment lines.**

## B2b: `sparse/`

32 sites (the plan said 28; `cl_SpMatrix.cpp` and `mumpstools.f90` are in the module too): `strumpacktools.cpp` (5), `cl_SolverParameters.hpp` (7), `cl_SolverMUMPS.cpp` (8 incl. two duplicated blocks), `cl_SolverParameters.cpp` (2), `cl_SolverWrapper.cpp` (2), `cl_SpMatrix.cpp`, `arpacktools.f90` (2), `mumpstools.f90` (the guideline's own §8 example). Dates, jury credits and "found by review" lines removed; the measurements that a default rests on kept without their dates (the 77 K drift behind `rel_tol = 1e-10`, the 337k-dof GAMG penalty, the 832k-dof deep-tree diagnostic). The three MUMPS symmetry-note copies consolidated: the wrapper keeps the rule and the reason for the error; `cl_IWG.hpp` and `cl_FEM_DofMgr_EigenValues.cpp` keep a short local reason plus a pointer; the measured 1D-Laplacian failure (both SYM values), the MUMPS guide section and the "lower triangle" correction now live once, in `src/sparse/doc/sparse_usage_guide.md` "Exploit Symmetry". Comment-only by the checker; `sparse/` now carries zero dated or credited comment lines.

### Retained from source comments (B2b)

- **BLR default tolerance stalls Newton, 2026-07-06** (formerly in `strumpacktools.cpp`): observed at N = 72k for every rank count; the reason compression is opt-in (`dl20260706_strumpack_blr_default.md` holds the campaign, not this figure).
- **Thread budget, 2026-08-21** (formerly in `cl_SolverWrapper.cpp`): until that date the OpenMP budget used `std::thread::hardware_concurrency()`, the logical count, and recommended it; on the 10-core / 20-thread workstation with ten ranks that recommended 20 threads on 10 cores, the configuration `Allrun` refuses. Grok, same day: an unset `OMP_NUM_THREADS` with several ranks per node is the affinity-mask default and needs no core count to diagnose, so the unknown-budget case must warn, not fall silent.

Codex sweep of B2b (`gpt-5.6-terra`/medium): ten edits, six violations, one suspected error. Applied: the shorter wrapper note (which also removes the `todo/` path Codex flagged as nonlocal — the guide carries it), the IWG pointer, the parameter-flag and thread-budget rewordings, the SpMatrix and `metis nodendp` sentences, the `mumpstools.f90` line (and the guideline's §8 example updated to match the tree). Declined in part, with reasons: the guide keeps the measured 1D-Laplacian failure, because the plan sends the evidence for that rule to the module doc, minus the history sentence about the older comment; the STRUMPACK gate keeps one sentence of mechanism ("at 1e-6 the linear exit test becomes the printed nonlinear residual and the timestep collapses"), because "do not restore the guard" needs its reason; the 77 K drift figures stay in `cl_SolverParameters.hpp` (deck and dof count dropped), because the 1e-10 default rests on them and no devlog holds them. The three other flagged measurements (ZERO_PIVOT run, 832k-dof deck, GAMG 74/700 counts) reduced to their consequence. The suspected error was real: "every BELFEM solve runs SYM = 0" is false — `cl_SolverSUPERLU.cpp` maps both symmetric modes to SuperLU's `SymmetricMode` — corrected to "every MUMPS solve".

## B2c: `physics/materials`

30 sites in 15 files: the six `Decision: Christian Messe, 2026-09-15` / `decided 2026-09-15` lines (Christian's ruling: name and date both go), the YBCO refit history, the `'custom'` rename diagnostic, the embed-guide audit credit, the three angle-contract "since 2026-08-16" notes (the migration warning for plugins compiled against the folded contract dropped; `dl20260816_dr69_unfolded_angle.md` holds it), the `both voices` / `Grok finding 4` / `C5 accepted` credits, and eleven sites in `powerlaws.hpp` including the `tmp/ai_exchange` pointer, the retired-formula history in the T-leg header (rewritten as the reason the closed forms are used) and the "staged out" narrative of the Bézier knot-motion term (rewritten as the finite-difference fact). Comment-only by the checker; `physics/materials` now carries zero dated or credited comment lines. Codex sweep pending.

### Retained from source comments (B2b and B2c)

- **Relative-tolerance drift, 2026-08-13** (formerly in `cl_SolverParameters.hpp`): tapestack3d, 337k thermal dofs, uniform 77 K start: 8.6e-5 K drift in 23 steps at rtol 1e-8, 2.6e-7 K at 1e-10. The figures stay in the header without the deck; the deck is here.
- **GAMG vs ASM, 2026-08-12** (same file): tapestack3d, 337k thermal dofs, rtol 1e-8: GAMG needed ~74 Krylov iterations on the two Picard iterates and ~700 on every Newton iterate of the same step; ASM ~107 on Picard and no degradation on Newton.
- **Deep separator tree, 2026-08-12** (same file): tapestack3d, 832k magnetic dofs, METIS_NodeND path: "supernodal tree was built from etree", 56 levels, STRUMPACK's own stack-overflow warning.
- **MC64 matching, 2026-07-05** (same file): serial `hphirun`, ZERO_PIVOT at Newton iterate 17 without matching.
- **YBCO thermal-conductivity refit, 2026-08-25** (formerly in `cl_Material_YBCO.cpp`): the optical channel's unit defect in `debye.f90` was fixed and alpha( T ) below the split made to follow cp, which brought grueneisen( T ) at 2 K from ~800 to ~2.1 and opened the phonon channel; b, d and Gamma were refitted then.

Codex sweep of B2c (`gpt-5.6-terra`/medium): twelve edits, all applied after verifying the two function names Codex introduced (`djc_eval_dT`, `dn_eval_dT` exist); two violations, both the same 4-33 % knot-motion figure in two places — kept once, in the temperature-derivative header where the "knot motion dominates" rule is stated, and the second copy replaced by a pointer; no suspected errors, and Codex independently re-checked the self-field-bridge sentence from B1 against `deval_dB()` and confirmed it.

## B2d: `fem/kernel` remainder, `core/assert`, the two IDE bylines

24 sites in eleven files: `cl_FEM_Calculator.hpp` (7: the `Grok phase-3 hardening` credit, four audit tags, two dated headings, the retired-formula note), `cl_FEM_Calculator.cpp` (3), `cl_FEM_DofMgr_EigenValues.{hpp,cpp}` (4: the restart-budget incident, the enum-default incident, two rulings), `cl_FEM_DofMgr_SolverData.cpp` (2, one of them a pure "lived here until" note deleted outright), the two controller sites my earlier grep had skipped because their lines contain a quote character, `cl_FEM_Postprocessor.cpp` (the dummy-send incident, now stated as the rule it taught), `cl_ThinShellFactory.hpp` (the physics-position line loses its two names and date), `core/assert.{cpp,hpp}` (the `L-21` pointer and the rejected-proposal date), and the `// Created by claude on 1/11/25.` CLion bylines on both `cl_EF_QUAD4TS` files. Comment-only by the checker; `fem/kernel` and `core` carry zero dated or credited comment lines. Codex sweep pending.

### Retained from source comments (B2d)

- **ARPACK restart budget, 2026-08-28** (formerly in `cl_FEM_DofMgr_EigenValues.cpp`): the product budget was calibrated on a 2026-08-10 measurement at n = 10428 and applied unguarded to a matrix nine times larger, which cut the exterior run from 300 restarts to 95 and broke a converging end; hence "the budget may only ever raise the restart cap".
- **Postprocessor dummy send, 2026-07-19 to 2026-08-27**: the incident is in `devlog/dl20260827_dr100_orphan_dummy_send.md` (RLC np=4 observed, gantry np=10 matching signature); the comment keeps the rule.

**Tooling defect found by this batch (plan D3):** gating B2d, the checker stopped after `core/assert.cpp` with exit 1 and no summary line. Cause: `g++ -fpreprocessed` rejects the `#aCheck` stringizing on a macro continuation line in `assert.hpp`, and `set -e` turned that into a silent abort — the last line printed read like a pass. Earlier verdicts stand (`assert.hpp` was not in any earlier diff, and every earlier run printed its summary). The stripper is now `scripts/strip_comments.py`, a tokenizer for C/C++ and Fortran; the checker fails loudly on a file it cannot strip, counts files, and always prints a summary. Re-gated: 40 files comment-only, both negative controls fail, unit input with `#a` in a macro, `//` inside a string, a raw string and a block comment strips correctly. The lesson is the one from L-09: a gate that can exit without its verdict line fails open.

Codex sweep of B2d (`gpt-5.6-terra`/medium): seven edits applied (the `comm_abort()` declaration pointer kept inside Codex's rewrite); one violation, right — the restart-budget note still carried the nine-times-larger comparison and the outcome, reduced to the invariant and its reason (the figures are in the retained-facts list above); no suspected errors. Re-gated, 40 files comment-only.

Two findings for the plan: 51 files in `src/` still start with a CLion `// Created by <name> on <date>.` byline, and 89 files carry no license block at all; both are outside B2's enumerated scope and are logged as open questions.

## B2e: `mesh/`

Nine sites: `cl_Mesh.cpp` (3: a "found 2026-08-15" on the global-variable ID map, the memdump-count compatibility note now phrased as "older dumps", an "audited" tag) and `cl_Mesh_PeriodicityFactory.cpp` (6: the `match_edges` header loses its plan phase and `todo/` pointer, the two `Grok EDGE-TIES combo III` labels, the `Codex phase-two final audit` credit, and the two dated corc cap-corner references). The topology rules the reviews asked to preserve — original-endpoint matching, straight / pure-mixed / half-cut profiles, the deferred pure tie and why encounter order must not decide — are untouched. Every date checked: `dl20260815_bc_globals.md`, `dl20260825_cap_corner_defect.md` and `dl20260825_halfcut_tie_fix.md` hold the incidents. The remaining match in `mesh/` is the CLion byline on `cl_Element_PENTA6TS.hpp`, which carries Christian's name and waits on O6. Comment-only by the checker. Codex sweep pending.

Codex sweep of B2e (`gpt-5.6-terra`/medium): seven edits, all applied (one restructured so the sentence that follows it still reads; "dof" kept as house vocabulary); no violations, no suspected errors.

## B2f: the modules B2 had not enumerated

Twenty-four sites in twenty files, found by the tree-wide listing: `circuit/cl_ElectricalCircuitFactory.cpp`, `comm/cl_Communicator.hpp`, `core/stringtools.cpp` (the tesla-exponent incident rewritten as the rule), `math/graph/graphtools.hpp` (the METIS heap corruption keeps its call-site mechanism, loses the reproduction date and the pipeline arrow), `math/tools/fn_ferrari.hpp` and `fn_cardano.hpp` (the three-AI verification credits go; "validated numerically against known quartics and random-coefficient sweeps" stays), `fem/interpolation/nedelec/cl_EF_{HEX8,TET10}.cpp`, `fem/iwg/cl_IWG.hpp` and `cl_IWG_Timestep.hpp`, six files in `fem/maxwell/` (including the exchange-thread pointer in `mt_maxwell_h.cpp` and the "5 x 24908 duplicate edge dofs" measurement, which `dl20260901_penalty_opt_in.md` holds), `fem/thermal/cl_ThermalBoundaryConditionFactory.cpp`, and the six integration headers whose "Regenerated 2026-09-03" sentence (four pyramid rules) and "UNUSED since 2026-09-03" (two tet rules) are now the bare facts; `dl20260903_intpoints_exactness.md` holds the regeneration. Two sites deliberately left: a date inside a user-facing error string in `cl_MaxwellFactory.cpp` (program text, not a comment) and the CLion byline on `cl_Communicator.cpp` (O6). Comment-only by the checker (64 files). **`src/` now carries zero dated, credited or working-record-pointing comment lines outside string literals, the closed core and the CLion bylines.** Codex sweep pending.

### Retained from source comments (B2f)

- **Tesla exponents, fixed 2026-08-30** (formerly in `stringtools.cpp` and `cl_MaxwellBoundaryConditionFactory.cpp`): until then every entry of the tesla family carried the VOLT exponents (the "V" block copied with only the string changed), so `check_unit` accepted a voltage where a flux density was required and refused a correct V·s/m²; the scale factors were never wrong. The by-name whitelist in the Maxwell boundary-condition factory was audited under that defect and is kept for that reason.
- **Cardano and Ferrari, 2026-07-23**: both kernels were audited by Claude, Codex and Grok on that day; the resolvent and factorization formulas were derived independently by two of them and agreed; Grok found the q/s Inf/NaN path on an ill-conditioned cubic that the biquadratic fallback now covers.
- **TET10 tables**: the shape-function tables were regenerated on 2026-08-14 after the eta/zeta exchange (`dl20260814_tet10_table_fix.md`); the generator, `tmp/tet10/tet10_function.m` and `tet10_derivatives.m`, lives outside the repository — logged as plan O8.

Codex sweep of B2f (`gpt-5.6-terra`/medium): fifteen edits, all applied (including a shorter version of the whole IWG symmetry-default comment, of which only the last sentence had been mine); two violations, both right — "( same hardening as the Maxwell twin )" was nonlocal provenance and went, and "validated numerically against known quartics and random-coefficient sweeps" in `fn_ferrari.hpp` was validation provenance, dropped once `tests/math/test_Polynomials.cpp` and `test_MathTools.cpp` were confirmed to exercise both kernels (the test is the record, not the comment); no suspected errors. Re-gated, 64 files comment-only. **B2 closed.**

## B3: commented-out code

About 270 lines removed from 46 files: every commented-out statement the census pattern found outside the closed core (110 lines minus the ten that are formulas or prose in statement shape — `// J = [ a b; c d ];`, `// curl( E_k ) = …`, the thermal `dFdX` formula, two controller sentences, the Magnesia fit note and the tutorial line in `example_user_material.cpp` — which stay), plus lines the pattern cannot see (the `//for` loops in the `Mesh` destructor, whose ownership comments stay; `//mesh::ConnectivityCalculator …`), plus ten `/* … */`-disabled bodies: the old `save_fields` and `load_fields` in `cl_Mesh.cpp` with their header declarations, the `dIdV` stub in `cl_FEMTwoTerminals.cpp`, the "eye hack" in `cl_CutFactory.cpp`, two blocks in `cl_BeltedTree.cpp`, the index-restore loop in `cl_ThinShellFactory.cpp`, the print loop in `cl_FEM_Tmatrix.cpp`, the hanging-dof debug print in `cl_FEM_DofMgr_SolverData.cpp`, and the Poisson probe in `cl_CutProcessorManual.cpp` (which announced itself as "safe to be deleted"). Orphaned labels above deleted calls (`// for debugging`, `// Shunn and Ham`, `// make sure that IDs are unique`, `// reset element`) went with them. Three reasons survived as comments where the disabled call carried one: the rejected Wassiljeva mixing rule in `cl_Gas.cpp` now named as the rejected alternative with its VDI citation, the `write_time()` omission in the VTK writer stated with its ParaView reason, and the naive `A( idx( i ), idx( j ) ) += …` form in `SolverData` relabelled as the equivalence it illustrates. No reference implementation worth a test was found: every disabled body was a debug print, a probe, or a superseded path. **Age check:** every removed line was searched with `git log -S` against the tree's history; none was touched by any commit after the 2026-09-04 re-init root, so nothing younger than twelve days was deleted (the pre-re-init history is not in this clone, so the "younger than a month" bar could only be checked to the re-init). Comment-only by the checker, 105 files. Census: commented-out statements 110 → 11 (ten kept by decision, one in the closed core). Codex sweep pending.

Codex sweep of B3 (`gpt-5.6-terra`/medium), two rounds. The first round's prompt carried an **empty diff**: `git diff src/ -- ':!…'` puts the path before the `--` and git rejects it as a revision; the wrapper's stdout began with `fatal: bad revision 'src/'`, and Codex still returned three language edits (on the additions the prompt had named in prose, which it read from the tree) and a "deleted reasons" list that contained exactly the two deletions the prompt had named. A verdict that mirrors the prompt is the tell; the check had not run. Rerun with `git diff -- src/ ':!…'` (475 removed lines in the diff): **pass** — no deleted reason, contract, warning or citation; the six removed `mG/mH( 8, k )` rows of `cl_EF_TRI6.cpp` are not a reference twin (the element allocates eight rows, re-enabling them would write out of bounds); the two retired `hphirun` toggles and the stale "original facet" note are dead, not contracts. One piece of narration debris noted for B5: the dangling `// if this is a facet on a tape,` line in `cl_FEM_Element.cpp`. The three language edits applied. Re-gated, 105 files comment-only. **B3 closed.**

## B4: TODO triage

Every `todo:` / `fixme:` remark outside the closed core, one outcome each:

| site | outcome |
|---|---|
| `core/typedefs.hpp` "move to constants" | deleted, stale wish |
| `cl_FEM_DofMgr_FieldData.hpp` "might not need this" | owned: `initialize_linear_projection_lists()` has no caller (grep: definition and declaration only); delete with its definition |
| `cl_IWG_TransientHeatConduction.cpp` "replace by cp and lambda" plus a dead `//aM +=` line B3's pattern had missed | owned: unit mass coefficient, use the material's cp; dead line deleted |
| `cl_FEM_Controller.hpp` "add flag to reset omega" and the duplicate `// end user settings` | both deleted |
| `cl_Mesh_Field.cpp` "remove field type" | owned: one live enum value, remove once vector fields are ruled out |
| `st_ProtoMesh.hpp` two "still missing" | owned (CurveData is still absent from the struct; the table accessor too) |
| `cl_Maxwell_FieldList.cpp` "add wave BCs" | owned |
| `cl_IWG.hpp` two "normal … todo: delete" | one owned line: `mNormal2D`/`mNormal3D` have no reader |
| `cl_IWG.hpp` two "old function, obsolete soon" | owned, with the two remaining callers named |
| `cl_FEM_DofMgr_SolverData.cpp` "fixme: check this" | deleted (narration plus mumbling) |
| `cl_ThermalFactory.cpp` "axisymmetry could also be a case" | owned: `AxSymmX`/`AxSymmY` exist in `ModelDimensionality` and are not selectable here |
| `cl_IWG.cpp` "get rid of these matrices" | deleted, orphaned by B3 |
| `cl_FEM_DofMgr_DofData.cpp` two periodic / thin-shell todos | one owned line |
| `cl_Maxwell_TMatrix.cpp` "these lines can be deleted" | owned, and corrected: `P` and `Q` feed the assertions below, so they are not deletable but should be guarded |
| `cl_IWG_Maxwell.cpp` "move into calculator" | owned |
| `test_facets.cpp` facet points location | owned |
| `fn_ODE_RK45.cpp` "why is this zero ?" | owned: unexplained, find out when |
| `cl_Gas.cpp` "use inline multiplication" | deleted, unclear wish |

Working-record pointers in ordinary comments, handled under the same rule: pointers to three **closed** plans (`periodic_thin_cut_continuity_fix`, `periodic_cap_cut_emission`, `anderson_picard_acceleration_plan`) dropped, the local reason kept; two **dangling** pointers (`timestep_collapse_mitigation_design.md`, `conditioning_diagnostic_backends.md` — neither exists anywhere in `todo/`) dropped along with the `MIT-3a` step label; the ngspice headers' pointer to the closed parser plan and its eleven `O5`/`O8`/`O9` decision IDs removed; the side-edge-fusing pointer kept, corrected to `todo/deferred/` and stripped of `O1`; the quadratic-shell pointer kept and turned into the owned TODO form, since it names a deferred plan for unimplemented work. Nothing was resolved; only the form changed. Comment-only by the checker, 121 files. Codex sweep pending.

**For Gregory (closed core, not edited):** eight identical `// todo: Optimize this to put directly in the while (cont) statement…` lines in `fn_Smith.hpp` (245, 260, 295, 312) and `fn_Smith.cpp` (96, 113, 167, 198), and two pointers in `cl_Cohomology.cpp` (356: `todo/thin_cut_nonunit_rectification_implementation.md`; 646: `todo/cut_pocket_removal_rules.md`) that the rule would turn into owned TODOs or drop.

Codex sweep of B4 (`gpt-5.6-terra`/medium): seven edits applied (one was a real defect of mine — a lone `*.` left in the ngspice header by the pointer removal); of the fifteen `TODO(cm)` lines Codex judged four still apologies (a topic without a completion condition, a dependent "retire with the one above", two "is not considered" statements), all four rewritten to name the job and the exit; one violation, right — the side-edge-fusing comment pointed at a `todo/` plan without being a TODO, now an owned TODO ("finish or remove the fuse"). No suspected errors. Re-gated, 121 files comment-only. **B4 closed.**

## B5a: `comm/commtools.hpp`

The narration is a template repeated across every `send` / `receive` / `broadcast` specialization: 24 × `get my id`, 20 × `wait until send is complete`, 19 × `send data`, 18 × `get the communication type`, 17 × `offset in data container`, 16 × `compute the chunks for this message`, 15 × `tidy up memory`, 11 × the two-line tag note, and so on — 331 comment lines removed by exact text (a delete list of the unique narration texts, applied file-wide), 418 → 81. Everything that carries a fact stayed: the matrix footprint explanation (once, in `broadcast( Matrix )`; the `send( Matrix )` copy reduced to a two-line pointer; the `receive( Matrix )` note kept because it states the receiver's own capacity guarantee), the `capacity()` warnings, the unchunked-broadcast int limit, the complex-type exclusion from `MPI_MAX`, the one-rank allreduce identity, the section headings. The `+ 1` tag convention that the deleted "For safety, the tag … is incremented" pair explained eleven times is documented once, and more precisely, at the `comm_tag()` definition in `commtools.cpp` (base tag for scalar and size traffic, base + 1 for payloads, and the matching-order consequence) — I had first placed it at the header declaration; Codex corrected the location. The twenty "wait until send is complete" lines are narration of `MPI_Wait`: Codex walked every `MPI_Isend` site and confirmed that no source buffer is touched or freed before its wait, so none of them was a warning. Comment-only by the checker. Codex check: conditional pass — no deleted reason, contract or citation; the added pointer reworded to a full sentence that names the rationale ("footprint note" was a coined label); and a list of 22 Doxygen blocks in the file that restate their signatures, which is B6's scope and is recorded there.

## B5b: `fem/kernel/cl_FEM_DofMgr_DofData.cpp`

The worst file of the census: 691 comment lines and 84 sub-banners over 4170 lines, most of it a per-entity template repeated five times (node, edge, face, cell, lambda): `grab the corresponding bitset` ×14, `allocate memory` ×12, `wait for other procs` ×10, `reset counter(s)` ×14, `loop over …` in a dozen forms, `make list with X IDs` ×5 (including the one B1 had corrected). 433 narration lines deleted by exact text and all 83 sub-banners removed (the `Identify X-Based DOFs` headings they framed stay as plain headings); 691 → 175 comment lines, 4170 → 3617 file lines. Kept: the `Step 1` … `Step 6` stage headings of each entity routine, the `reorder_dofs` stages, the four LINE2/LINE3 hanging-dof cases, the index-space rules (`[free_0…free_n, fixed_0…fixed_m]`, `my_index` in SEPARATE mode, the counter that is deliberately not reset between free and fixed dofs), the cascade and flatten reasons, the Maxwell-specific orphaned and abstract node flags, the order dependencies (`save abstract dofs before we disconnect the mesh`, `we void all indices to provoke an error`). Four rewrites: two typos (`crate`, `Crete`), one history phrasing (`no longer participate`), one warning reworded from a shout to a reason. Comment-only by the checker. Codex check: one deleted reason — the garbled `so we need to identify the over their mesh bases` carried the consequence that hanging candidates are first identified through their mesh bases because the dofs are not yet linked; restored as one line after the surviving premise; four narration survivors it listed removed; the four rewrites passed. **B5b closed.**

## B5c: the rest of `fem/kernel`

Done so far: `cl_FEM_DofManager.cpp` (236 → 84 comment lines) and `cl_FEM_DofMgr_SolverData.cpp` (350 → 197; two plan step IDs `O2`/`O3` that B2 had missed removed; one dead initializer line and one dead `#include` that B3's pattern could not see deleted), plus the common narration set (`wait for other procs`, `allocate memory`, `loop over all …`, `initialize counter`, `get dimension of element Jacobian`, `add contribution to system matrix`, …) applied to all seven files. The two-phase-solve and residual-semantics comments in `SolverData` were left untouched: they are reasons and contracts. Then the five others: `cl_FEM_DofMgr_FieldData.cpp` (192 → 39), `cl_FEM_Kernel.cpp` (125 → 79; the thin-shell exchange contract, the ownership rules and the share/receive collectiveness warning untouched), `cl_FEM_Element.cpp` (152 → 52; two bare "need to do something fancy here for #Facedof" / "here is where we would link #celldof" apologies that carried no `todo` word became owned TODOs), `cl_FEM_Calculator.cpp` (423 → 315; three dead `// uint … ;` declarations that B3's pattern could not see deleted; the weight-normalization reasoning and every dispatch/peer/aura comment kept; one typo), `cl_ThinShellFactory.cpp` (251 → 153; the `#BEGIN/#END SIDE CONNECTOR MOD` position markers removed; a note addressed `@AI:` rewritten as the warning it contained). Seven files, 1,829 → 895 comment lines. Comment-only by the checker. Codex check: **"do not sign off unchanged"** — six deleted fragments had been continuation lines of multi-line comments that my exact-text delete list matched as generic narration (`system matrix`, `master and slave`, `receive_submesh`, `duplicate`, `identity`, and the two five-line payload-layout lists headed `contains` in `SolverData`, which are the positional contracts `compute_dof_dof_connectivity()` reads). All restored, the two layouts as one two-line contract each; five surviving narration lines it listed removed (`auto set blocks…`, two `end loop over all procs`, `bottom layer` / `top layer`); the two owned TODOs in `cl_FEM_Element.cpp` given exit conditions; the sideset warning made concrete; the "plated walls use copper" sentence qualified to what the check enforces (any constant mu); one narration line I had merely typo-fixed deleted. Recorded as plan D5 with the rule that protects continuation lines in every later pass; a heuristic scan of all fifteen B5 files found no further amputation. Re-gated, 126 files comment-only.

## B5d: `mesh/`

Six files, two passes each (a shared narration set, then the file's own singletons): `cl_Mesh.cpp` (282 → 168; every `(owned objects)` / `(references, just pointer storage)` member annotation, the finalize idempotence reasoning, the checksum contract and the facet-recollection assumptions kept; two dead `//Timer` / `//proc_t` lines deleted), `cl_Mesh_GmshReader.cpp` (191 → 50; the gmsh format notes, the "8 and 9 are swapped" node-order fact and the dangling-vertex reasoning kept), `cl_Mesh_OrderConverter.cpp` (160 → 32; almost pure narration, the create-nodes-on-unique/shared-facets stage headings kept), `cl_ProtoMesh.cpp` (133 → 58; the twin-edge and node-pair-key reasoning untouched), `cl_FaceFactory.cpp` (111 → 39; the key128 limits and the face_key_2d contract kept; one typo), `cl_Mesh_VtkWriter.cpp` (87 → 12). Comment-only by the checker. Codex check: two deleted reasons — one continuation amputation of the D5 kind in `cl_ProtoMesh.cpp` (`we need to count the number of facets` / `over the sidesets, because …`, now one sentence) that the heuristic scan had missed because the deleted line came *before* the survivor, and `keep only flagged entities` in the Gmsh reader, which the two named functions' behavior depended on — both restored; five survivors removed, one of them (`backup original node indices` above a loop that counts flagged nodes) inaccurate as well as narrative. The D5 guard now also has to look at the *next* line, not only the previous one. **B5d closed.**

## B5e: `sparse/cl_SpMatrix.cpp`, `io/`

`cl_SpMatrix.cpp` (285 → 161; the eleven `// move …` labels of the move constructor, the save/load narration and the base-conversion narration went; every linked-structure, indexing-base, sentinel-slot and MKL-path reason stayed), `io/hdf5_tools.hpp` (196 → 82; eight `// test if dataset exists`, eight `// check for error`, eight `// check datatype`, one dead `//std::string tLabel` line; the on-disk type-width notes, the transpose layout, the UTF-8 cset reason and the H5Literate contract kept), `io/cl_HDF5.cpp` (68 → 41; fourteen `// call interface`; the navigation-stack reason kept), `io/cl_XML.cpp` (20 → 14), `io/cl_Input_Section.cpp` (74 → 55; the brace-range invariant, the thin-shell `@` rule and the SI-unit assumption kept). The continuation guard from D5 reported one protected line in `hdf5_tools.hpp`; on inspection the only candidate still present was the dead `//std::string tLabel` line, which the guard had no reason to keep (its predecessor is a brace) — I had first written here that the protected line was a `// matrix dimensions` continuation, which was a guess and wrong; the dead line is now deleted and the guard's count is not trusted as a report of what it kept. Comment-only by the checker, 128 files. Codex check: no deleted reason and no truncation; two surviving `// File creation property list identifier` labels in `cl_HDF5.cpp` were wrong (they sat above the `flags` argument of `H5Fopen`) and are removed. **B5e closed.**

## B5f and B5g: the remaining narration

B5f: `fn_Mesh_compute_surface_normals.cpp` (143 → 20), `cl_IWG.cpp` (247 → 121), `cl_MaxwellFactory.cpp` (503 → 416; almost everything left is reasoning about the .bfm cache, the interface pairing and the domain-type rules), `cl_Gas.cpp` (534 → 330; the NASA RP-1311, VDI and Rist citations, the Newton and bisection reasoning and the entropy-scale notes kept; five dead lines B3's pattern could not see deleted; three typos), `cl_CutFactory.cpp` (222 → 100; the `#BEGIN/#END SIDE CONNECTOR MOD` markers and the orphaned `// for debugging` label went). B5g, the tail of the census: `cl_FEM_Controller.cpp` (1501 → 1444: 39 narration lines, three dead `//tOmega *=` lines, and a second class of provenance B2's pattern had missed — audit-round labels without an AI name, `( round-3 C5 )`, `( round-3 N1: … )`, `R3b tripwire:`, `( ts34 trace, dl20260727 )`, `( Garber R9 trace, dl20260806 )`; the same form removed from `cl_FEM_Controller.hpp`, `cl_SolverSTRUMPACK.cpp` and a devlog pointer from `cl_SpMatrix.hpp`), `cl_ElectricalCircuitFactory.cpp` (101 → 71; seven more netlist decision IDs `O5`/`O6`/`O7`/`6b` that B4 had not reached), `cl_EdgeFactory.cpp` (75 → 23), `cl_FEM_DofMgr_BlockData.cpp` (83 → 51; the aura reasoning intact), `cl_CutProcessor.cpp` (113 → 46; the quotient self-cancel reasoning intact), `cl_FEM_SideSet.cpp` (88 → 49; the condensation contract intact), `cl_CurveFactory.cpp` (175 → 127; the Step headings, the node-enumeration diagrams and the traversal reasoning kept), `cl_Mesh_ExodusWriter.cpp` (73 → 43). The D5 guard, now two-directional, protected candidates in six places; each was inspected and cleared by hand. Comment-only by the checker, 133 files. Codex checks: B5g clean (no deleted reason, both directions checked; five survivors it listed removed, one of them — "Each group becomes one boundary condition" — misleading beside a check that refuses more than one). B5f: three lost reasons, all restored — two of them were exactly the continuation lines the D5 guard had protected and I had cleared "as plain narration" without reading their context (`// orient_terminal_curves_sub()`, `// duplicate_and_relink_facets`); the third, the ownership/level snapshot in `cl_MaxwellFactory.cpp`, was weakly worded ("just in case") but stated why two arrays are copied and restored, and is now said properly. Two survivors rewritten or removed. Recorded in the plan under D5: a guard hit is cleared by reading, never by re-running with the guard off. **B5 closed.**

## B6: Doxygen triage

Every `/** */` block in the plan's headers read whole and judged by §5 of the guideline: a block stays only if it states a contract the signature does not. `cl_Material.hpp`: 86 blocks, 27 deleted (`Constructor`, `Get the material type`, `Clear flag`, `Check if …`, the seven protected setters' one-liners, the enum and typedef briefs) and four reduced to the one-line contract they contained (the destructor's ownership, the empty-number return, the RRR definition, the molar-mass unit); every block carrying units on a parameter or a return stayed (they are the contract). `cl_MaterialFactory.hpp` and `cl_JcFunction.hpp`: the default constructor/destructor and `Check if …` blocks deleted, the dependency-declaration and derivative-contract blocks kept. `cl_Material_Metal.hpp` and `cl_Material_YBCO.hpp`: destructor blocks reduced to their ownership sentence, everything else kept (units, formulas, citations). `cl_Mesh_ConnectivityCalculator.hpp`: the constructor and `Compute node-to-element connectivity` deleted; the destructor block became the guideline's own example, `Nothing to release: every container is owned by the mesh.`. `comm/commtools.hpp`: the 27 per-overload `\brief Sends a … / \param aData The … to send` blocks deleted (Codex's B5a list); the tag-ordering tripwire, the allreduce reasons, the raw-array preallocation contract and the offsets contract kept — note for sweep 2 that none of the send/receive/broadcast overloads states its collectiveness, which is the contract they should carry. A tree-wide pass over the 53 blocks the census pattern still flagged found four bare ones (`Destructor` ×2, `Check if a matrix has been populated`, `Get spline for a property`) and 49 that carry units, signatures, citations or ownership. Census: `@brief` restatements by pattern 82 → 48, all remaining ones read and kept on purpose. Comment-only by the checker, 142 files. B6c (`make doc` spot-check) is Christian's.

Codex check (`gpt-5.6-terra`/medium, deleted contracts / reductions / survivors): every non-MPI deletion and all seven reductions clean. One real catch: the deleted `commtools.hpp` briefs were the only lines saying that `comm_barrier` and the `broadcast` overloads involve every rank, and §5 keeps collectiveness. Codex proposed "Collective: every rank must call" for the `collect` overloads too, which is wrong: `collect` is the receiving half of `distribute`/`collect` and calling it on every rank without senders deadlocks (`src/comm/doc/comm_usage_guide.md`, the rank-guard note), so the old "Collects … from all processes" brief had been misleading as well. Applied: a one-line contract on the barrier and the six `broadcast` overloads ("Collective: every rank calls it, and aRoot's data lands on all ranks"), on four `collect` overloads ("Receives only, one slot per rank. Not collective: every other rank must send, through distribute() or send()"), and on the flat-buffer `collect( T *, aOffsets )` ("Rank 0 only, receives only: rank p's block lands at aOffsets( p ) …", its loop starts at p = 1). That closes the sweep-2 note above for `broadcast`/`collect`; `send`/`receive`/`share` still carry no statement, but those are point-to-point and the usage guide is explicit. Two survivors Codex flagged as signature restatements were replaced by the contract each function actually has: `SplineLookupTable::set_spline` (takes ownership; the table deletes it; nullptr keeps the stored spline and refreshes only the dispatch) and `Material::create_spline` (NaN end derivative selects the parabolic end condition, a finite one is imposed as the tangent; requires the property and `T_max`). Gate after the fix: `check_comment_only: OK (HEAD, 141 files)`. B6 closed.

## Final census before the R9 gate

`scripts/comment_census.py` over `src/` after B6, closed core in its own row:

| scope | lines | comment | banner | sub-banner | narration cand. | dead | todo | dated | prov | record | `@brief`? |
|---|---|---|---|---|---|---|---|---|---|---|---|
| start of sweep | 292,287 | 31,723 | — | — | 3,023 | 110 | 33 | 128 | 28 | 3 | 82 |
| after B6 | 292,287 | 27,837 | 9,789 | 454 | 1,288 | 11 | 31 | 4 | 1 | 0 | 47 |
| of which closed core | 5,690 | 645 | 286 | 0 | 68 | 1 | 10 | 2 | 1 | 0 | 0 |

The four dated lines are two CLion bylines (`cl_Element_PENTA6TS.hpp`, `cl_Communicator.cpp`; O6, Christian's call) and two in the closed core; the one provenance line is in the closed core. Outside those two exclusions every dated, provenance and working-record count is zero. The eleven dead statements are the ten formula-shaped lines kept by decision in B3 plus one in the closed core. Every remaining TODO carries an owner and an exit condition. The 47 `@brief` restatements by pattern were all read; each carries a unit, an ownership rule or a citation that the pattern cannot see. Sub-banners were removed only in touched files (B7 as planned). Gate: `check_comment_only: OK (HEAD, 141 files)`. Still owed for R9: `make check` on a debug and a release tree and the `make doc` spot-check, both Christian's; the range gate re-run over the sweep commits once they exist.

## R9: closing gate

Christian ran `make check` on both backends (Armadillo and Blaze): passes. That is the executable gate for the whole sweep; every batch above is now verified, not only reviewed, to the extent that the suite covers the touched files (the MPI coverage caveat of `make check` still applies outside `tests/comm`). `make doc` failed on the new page: `update_doc_index.py` refuses a `doc/` page that is not filed in its `PAGE_GROUPS` table, so a new document cannot land unfiled. `doc_commenting_guidelines` is now filed under "Extend a module" next to the coding philosophy, and `doc/mainpage.md` was regenerated by the script. Christian's rulings on the open questions: O4, the B6c spot-check is enough, no before/after diff of the generated reference; O6, every CLion byline is replaced by the license header; O7, every source file carries the header, including the python files shipped in `examples/` and elsewhere; O8, the MATLAB scripts in the pre-reset `tmp/` are surveyed, the ones worth keeping are cleaned up, documented and moved into the documentation tree.

## O6 and O7: license headers

One mechanical pass with the new `scripts/add_license_header.py` (dry run, then `--apply`): 121 CLion byline blocks removed (the `// Created by <name> on <date>.` line with the bare `//` lines around it; `msh2exo.cpp` carried two, `cl_Element_PENTA6TS.hpp` had one between its two include guards, which the file still has), and the license block added to 238 files that lacked it: 90 `.cpp`, 46 `.hpp`, 106 `.py`. The C and Fortran forms are the existing block verbatim; the python form is the same text behind `#`, placed after a shebang and an encoding line. Every python file in `examples/`, `python/`, `scripts/` and `share/` was checked for third-party attribution first: none (the three hits on "licence" or "written by" are prose inside docstrings). Gates: `check_comment_only: OK (HEAD, 266 files)` for the C++, and the tool compares the `ast` of every python file before and after and refuses on a difference. Not touched: `.sh`, `.cmake` and `CMakeLists.txt` (not asked; 43 + 58 + 5 files), and two python files that already emit an invalid-escape `SyntaxWarning` on their own (`examples/tape_quench_usermat/plot_iv.py`, `scripts/resolve_doc_cites.py`), which is code, not comments. One dated line outside the sweep's scope surfaced on the way: `tests/physics/backendfree/test_MaterialBackendFree.cpp:18` cites "the 2026-07-02 SplineLookupTable refactor"; `tests/` was never in this plan and goes to sweep 2.

## O8: the MATLAB derivation scripts

Christian's ruling: survey the 56 MATLAB scripts in the pre-reset `tmp/` (now `~/codes/belfem_transfer/tmp`), keep what is worth preserving, clean it up, document it, move it into the documentation tree; the keep decision to the jury, the comment language to Codex. Plan: `todo/matlab_derivations_recovery.md`.

**Survey** (read-only subagent, every file read): 22 proposed keepers, mostly labelled "exact generator of a C++ table". Before dispatch I checked the load-bearing claims: both dangling citations in `nedelec_derivation.md` (`:270`, `:443`) are real; the 4×12 magic table in `init_points.m` equals `fn_IF_initialize_integration_points_on_facet.cpp:898`; all 108 entries of `penta18.m` equal `cl_IF_PENTA18.hpp::d2NdXi2` (whitespace-normalised compare, 0 differences); three duplicate pairs by md5.

**Jury** (Codex `gpt-6-astra`/high, Grok `grok-4.6`/high, blind): 33 findings, 31 confirmed against the cited lines, two refuted (a "byte-identical" claim the plan never made; a slug remark). The confirmed ones cut the keep set from 26 files to 14 and withdrew the "exact generator" label from everything but the TET10 zero-tests and `penta18.m`: the TRI6 `fragment.m` carries twice the C++ coefficients and the `:270` citation is the TRI3 curl section, not TRI6; `pyra14.m` stops at its Vandermonde matrix and uses a pyramid with the base at `z = -1` against the C++ base at `zeta = 0`; `nabla_derivative.m` uses its cofactors before defining them and derives a route the C++ (`mCurv` from the map's Hessian) does not take; the magic-table check needs a driver that only a dropped script provides and tests one of twelve columns; `hex.m` lists the same slave ordering twice; `check_linear.m` integrates over the unit square, not a circulation; `tet/main.m`'s `P10/P20/P35` have no consumer in `src/`; `Doxyfile.in` lists `*.m` in `FILE_PATTERNS` with `INPUT = src`, so anything copied under `src/` would be parsed as source; R5 as written edited C++ against the plan's own scope guard; the runtime verifiers (`test_LagrangeInterpolation.cpp`, `test_FacetIntegrationPoints.cpp`, `test_EdgeFunctions.cpp`) already exist, so what the recovery preserves is the symbolic construction, not verification. Record: `tmp/ai_exchange/review_matlab_recovery.md`.

**Kept (14)**, under `src/fem/interpolation/doc/matlab/`: `nedelec_tet10/` (`tet10_generate`, `tet10_shortcuts`, `tet10_function`, `tet10_derivatives`, `parse_main`, `defelement`, `tet10_lagrange`), `nedelec_tri/` (`curl`, `fragment`, `secondordertri`, `check_quad`), `lagrange/` (`penta18`), `facets/` (`orientation`, `facet`). Dropped 39 plus the three interface T-matrix experiments (O1, Christian's call, default drop); all stay in the transfer tree.

**Cleanup** (comment-only, checked by a stripped-comment diff against the originals: 13 scripts identical in code, `fragment.m` differs by the preamble and loop that make it runnable): `%` license header via `scripts/add_license_header.py` (extended to `.m`), a purpose comment on every script naming what it generates, derives or checks and which C++ table it concerns, the class-(b) notes saying where they differ from the C++ (sign convention in `curl.m`, coefficient factor in `fragment.m`), the three "corrected 2026-08-14" lines rewritten as the node-map contract, six German comments in `tet10_lagrange.m` translated, one typo. New: `check_tet10.m` and `check_penta18.m` (asserting drivers), `compare_tables.py` (sympy; the MATLAB transcriptions against the live C++: 192 TET10 entries, 108 PENTA18 entries, 0 differences), five READMEs, `*/doc/matlab/*` in `EXCLUDE_PATTERNS`, the directory registered in `src/fem/interpolation/doc/README.md`, the two citations in `nedelec_derivation.md` repaired (`:270` now names `curl.m` with its sign caveat).

**Gate:** MATLAB turned out not to be installed (the shell alias points at a missing directory). Octave 9 with the sympy-backed `symbolic` package is; the drivers and the runnable notes were run there (`octave-cli --eval "pkg load symbolic; check_tet10"` from the script's directory). `check_tet10`: exit 0, every residual zero, 24 value rows and 15 derivative tables. `check_penta18`: exit 0, all 108 entries. `defelement.m`: twenty all-zero blocks. `parse_main.m`, `tet10_lagrange.m`, `curl.m`, `check_quad.m` (prints the factor 2), `secondordertri.m` run clean. `fragment.m` failed once on my own preamble (`sym( name, [ 2 1 ], 'real' )` is MATLAB-only; Octave answers "MatrixSymbol do not support assumptions"), rewritten with four scalar symbols, then clean. Octave warns that `curl.m` shadows its built-in `curl`; harmless as a script run from its directory, noted in the README. The two facet notebooks are interactive plots and were not run. With that, the TET10 and PENTA18 tables are verified against their generators, not only reviewed; sympy's deprecation notice about non-Expr matrix entries comes from the `symbolic` package itself, not from the scripts.

**Codex sweep** (`gpt-5.6-terra`/medium) over the comments and READMEs: fifty language edits proposed; the substantive ones applied, the accent edits (Nedelec → Nédélec in `%` comments) declined to keep the scripts ASCII. Four suspected errors, all real and fixed by content, not wording: the facets README said no toolbox is needed while `facet.m` calls `syms`; the triangle README said all four notes print a result while `secondordertri.m` ends every line with a semicolon; a translated comment in `tet10_lagrange.m` called row vectors a column vector (the German original was misplaced); and the DefElement URL in `defelement.m` is the triangle page while the body is the tetrahedron, which the header comment now states. Final stripped-comment diff against the originals: 13 of 13 identical in code, `fragment.m` differs by its preamble only. `update_doc_index.py --check`: 0 files out of date.

Files: `src/fem/interpolation/doc/matlab/` (new: 14 scripts, `check_tet10.m`, `check_penta18.m`, `compare_tables.py`, five `README.md`), `src/fem/interpolation/doc/README.md` (section), `src/fem/interpolation/doc/nedelec_derivation.md` (two citations), `Doxyfile.in` (`EXCLUDE_PATTERNS`), `scripts/add_license_header.py` (`.m`), `more/gmsh/meshtool.m` (header only), `todo/matlab_derivations_recovery.md` (new), `todo/README.md`, `tmp/ai_exchange/review_matlab_recovery.md` and `matlab_recovery_sweep.md` (ephemeral), `tmp/matlab_survey/` (read-only mirror for the auditors, ephemeral).


## Committed and closed

`caf2045` (guideline, tooling, conventions, the three same-day devlogs), `05a4eb5` (the sweep and the license headers), `08c306b` (the MATLAB recovery). Range gate after the commits: `check_comment_only: OK (caf2045^..HEAD, 266 files)`. Both plans moved to `todo/closed/` with their Status lines rewritten as summaries; O1 of the MATLAB plan resolved by Christian accepting the default (dropped). Next: sweep 2, the Grok + Codex second-tier list over the swept tree (`comment_cleanup_sweep_1.md` §9).

## Sweep 2: the second-tier lists

Step 4 of the campaign. The same brief to Codex (`gpt-6-astra`/high) and Grok (`grok-4.6`/high), blind, plus a Claude subagent on the ten entry-point headers: what contracts are missing, what is still wrong, what sweep 1 lost. Grok could not run the census or the diff (no shell) and said so; Codex read the diff for the files it opened. The three lists agree on every top item: `send`/`receive`/`share`/`distribute` carry no collectiveness statement while `broadcast`/`collect` now do; `Cell<T*>` non-owning and `Cell( aReserve )` are stated nowhere in the container; the mesh accessors and adopting setters say nothing about ownership; the `Kernel` brief claims to own a mesh it borrows (`mOwnMesh` is never set); three constructors communicate; `Solver::solve` swallows the soft-fail verdict; the wrapper's `get_*` say "if supported" and abort. Sweep 1's own errors, verified in the diff: the `comm_tag + 1` protocol reason deleted at eleven sites, the wait-warning deleted rather than restated once, my B4 TODO rewrite in `cl_IWG.hpp:982` naming a caller that calls a different function, and a commented-out `set_pardiso` B3 missed. Code findings (not comments), verified: `fn_Mesh_integrate_scalar_over_sidesets.cpp:167-171` broadcasts the master rank instead of the value, so non-master ranks return zero; `IWG::N` declared and never defined; wrapper `mX`/`mY` never assigned; `IWG::collect_node_coords` writes `nDim + 1` columns. Plan: `todo/comment_cleanup_sweep_2.md`, pending approval; nothing edited. Record: `tmp/ai_exchange/comment_sweep_2_tier_list.md`.

## Sweep 2 execution (2026-09-17)

Christian's ruling before the start: comments only in this session, behavior-altering changes go to `todo/code_findings_from_comment_sweeps.md` for a code session with the jury loop, obvious fixes allowed. Batches C1–C8 written in one sitting, each gated comment-only and each handed to Codex (`gpt-5.6-terra`/high) with a new brief: for every added sentence, is it true of the code as written, quote the line. That check earned its cost every time. What it corrected: `comm_check` does not abort "on the failing rank only", `BELFEM_ERROR` ends the whole run through `MPI_Abort`; `share()` is sender-only, not root-only, and calling it everywhere hangs only for want of MPI buffering, not by guarantee; `receive( Cell )` fills with 0 rather than "value-initializes"; `Cell::data()` is invalidated by reallocation, not by every size change; `shrink_to_fit()` is a request; a hand offset from `data()` is legitimate with `spacing()`, so the linalg contract says that instead of "never"; `mAbstractNodes` is not a view onto `mNodes` because the cut factory fills it directly (F8 filed); `scale_mesh` zeroes a 2D mesh's third coordinate; `time_step` is not 1-based; the kernel's chained-worker path allocates a placeholder mesh nobody deletes (F9 filed); `collect_fields` barriers inside `FieldData`; `set_solver` synchronizes through the Solver constructor, not itself; `rhs_vector` is the collected global vector on the master; `create_field` links without checking adoption; `Postprocessor::initialize` is a no-op once initialized; `omega` is not clamped; a calculator scratch reference survives, its contents do not; `allocate_work_matrices` has no once-only guard; the calculator's `mM` is not sized where `mK` is; a `Gradient` buffer is sized at construction. Also caught: my helper had placed nineteen contract blocks between a return type and its declaration; fixed by a pass that moves any such block above the return type. Every correction was re-read against the cited line before it was applied.

Declined after reading: the `@brief Set a user-defined function …` blocks in `cl_Material.hpp` and `cl_Material_UserDefined.hpp` (Grok's candidates). Each carries the signature the plugin must provide, the caller, or the library path resolution, which is the contract, so they stay.

C8 was done by a guarded pass over twenty files: a comment line was deleted only if it matched the census narration pattern, stood alone (no comment line before or after it, the D5 rule made mechanical), was at most 60 characters, and carried no reason marker (`so that`, `because`, `must`, `note`, `unless`, …); eight lines survived on that marker, 224 went; Codex checked all 267 deletions against both neighbours and found no lost reason and no cut continuation, and listed 34 survivors that went the same way. C9, the language sweep of the added text in two halves: about seventy edits applied; it also caught one contract error of mine (`time_step() const` described as a writable reference) and two wrong old lines ("new beta" over a bisection midpoint in `cl_GT_RefGas.cpp`, "a nonzero LHS" over a length test in `cl_SolverPETSC.cpp`). Comment-only by the checker throughout (`check_comment_only: OK (HEAD, 54 files)`). Census after sweep 2: 27,837 → 27,281 comment lines; narration candidates 1,288 → about 1,030; the two remaining dated lines and the provenance line are in the closed core. Not committed: the closing gate is Christian's `make check` on both backends.

## Code findings, tier 1 and the two deep scans (2026-09-17, autopilot)

Christian's ruling: sweep 2 committed after `make check` on both trees; F1–F7 in one tier through the plan + audit → code + audit loop, F7 with comment-only write access; F8 and F9 each get two independent deep scans (Codex `gpt-6-astra`/high, Grok `grok-4.6`/high); builds and Valgrind allowed in the sandbox.

**Plan jury** on the tier-1 edits (`tmp/ai_exchange/review_code_findings_tier1.md`): F4's proposed loop bound was right but its column contract wrong (a sideset work matrix is sized `nDim + 1` on purpose); the "MPI only in tests/comm" claim was false (`sparsempi` runs at 2 and 4 ranks); the header of F1 had to change in the same diff; and Grok pointed out that Christian's comment-only grant for F7 is not the named authorization protocol §7.1 requires from Gregory. I had already applied F7 by then; reverted, the closed core is untouched, the list stays in the todo for Gregory.

**Code jury** on the diff (F1 `broadcast( tValue, aMasterRank )` with the barrier dropped and the local renamed; F2 the undefined `IWG::N` deleted; F3 `mX`/`mY`/`x()`/`y()` deleted; F4 `<` bound plus self-sizing and an assert; F5 the outer include guard): code approved by both; both refuted my F4 comment, which claimed `link_to_group()` sets the dimension. The live writer is the derived IWG's constructor (the problem dimension, 2 or 3), Maxwell never sets it, and a 3D sideset IWG gets three columns. Header, assert message and guide line rewritten to say so.

**F8 refuted** by both scans, independently and with the same evidence: `CutProcessor::collect_duplicates()` appends every abstract node to `Mesh::nodes()` (`cl_CutProcessor.cpp:1230`) before the factory moves the pointer list, so `~Mesh` deletes them once; the sweep-2 comment that produced the finding described the unused setter, not the live adoption, and is corrected. Both warned that routing the factory through `set_abstract_nodes()` on top of the append would double-free. Both found the real leak instead, F10: `CutSet::~CutSet` deletes `mBitset` but not `mNodeBitset` (`cl_CutSet.cpp:28-38`); one line added.

**F9 confirmed** by both scans with the same one-line fix: on the chained worker branch of `distribute_mesh()` the child allocates an empty placeholder mesh and `mOwnMesh` is never set, so it leaks once per worker per coupled run; `mOwnMesh = true ;` right after the allocation, and nowhere else (on rank 0 it would delete the file mesh). Codex found three more leaks on the coupled setup path, filed as F11–F13 for a code session: the thermal equation is handed to `create_field()` without `add_equation()`, the thermal boundary-condition factory is never deleted, and the magnetic worker's original empty mesh from `MaxwellFactory` has no owner.

Gates, all run in the sandbox on two cores: debug tree built, `make check` passed (fast 8, mpi 4), Codex audit of the two one-line fixes (both correct), Valgrind on `dipole` for one timestep (`cl_CutSet.cpp:31` in no loss record; the abstract nodes only indirectly lost, so the mesh owns them, as the scans said), Valgrind on `tapestack3d` at two ranks for one timestep (`distribute_mesh():698` in no loss record on either rank). The gates also found what nobody had filed: `Element::relink_dofs` leaks one malloc per consolidated element (F14), one sideset calculator leaks (F15), 2.9 MB is indirectly lost on the master from a root that is not a mesh allocation (F16), the worker's empty mesh from `MaxwellFactory` is confirmed lost (F13), and 120 B at `distribute_mesh():652` on the master (F17). All parked with the others for a code session. Tier 1 committed (`70364c1`); afterwards the release tree was rebuilt and `make check` passed there too.

## Code findings, tier 2 (2026-09-18)

F11–F17 through the same loop. Plan jury: five of seven approved as written; F13 (the factory-created meshes) blocked, because deleting the mesh in `~MaxwellFactory`'s body would run while the factory's own `shared_ptr` still kept the kernel alive, so the kernel's dof managers would die against a freed mesh. Rewritten as the reset-then-delete sequence Grok proposed; the code jury blocked that too, both auditors, because it still made the mesh's lifetime depend on the factory holding the last kernel handle. Third form, Codex-audited correct on all seven points: the kernel claims the mesh at construction exactly as it claims the parameters (`Kernel::claim_mesh_ownership()`), the factory keeps a read-only pointer and deletes the mesh only if it never built a kernel, and the mesh dies with the last kernel handle, after every borrower in the drivers. The jury also added: `mOwnThermalEquation = false` right after `add_equation()` (an exception in `create_field()` would otherwise double-delete), deleted copy and move on both factories, unconditional `free()` of both dof arrays in `Element` (a `malloc( 0 )` result would have escaped the count guard), and two comments corrected. Gates in the sandbox: `make check` debug and release, Valgrind `dipole` serial and `tapestack3d` two ranks, one timestep each: every fixed site absent from every loss record on every rank; the dipole run now reports nothing but Open MPI's own initialization blocks. Rank 0 of the two-rank run exposed two more: the master's node dof containers from `DofData::connect_dofs_to_mesh()` (a million bytes, F18) and 2,211 reader-created elements lost after distribution (F19). Committed by explicit path; a peer session's CMake/PMIx work in the same checkout was left untouched.

## Code findings, tier 3 (2026-09-18)

F18 and F19, the two rank-0 losses the tier-2 gate exposed. Roots read first: F19 turned out not to be a distribution effect at all. In a 3D Gmsh mesh the reader files the line elements into `mBoundaryEdges`, nothing else holds them, and `~Mesh` never deleted the list; the serial 2D deck was clean only because its line elements are facets. F18 was the master's hanging-dof path: `connect_dofs_to_mesh()` fills every basis' dof container over the whole dof list, the list is then split into hanging and non-hanging, and the disconnect walked the non-hanging half. The plan jury (Codex `gpt-6-astra`/high, Grok 4.6/high, on Christian's instruction to use astra for a behavior-affecting case) approved F19 as planned and widened F18: Grok traced the second half of the mechanism, two kernels share each basis' container, the thermal connect mallocs over the magnetic leftover without freeing and zeroes the count, and both `reset_dof_container()` and `~Basis` were keyed on that count, so the replacement was never freed either. A walk over the hanging list alone would have rested on the factories' call order, the failure mode the F13 jury rejected. Landed: the walk, plus the container hardened: pointer initialized to null, free before replace, frees keyed on the pointer. The code jury approved the executable changes and rejected one comment (the `memory()` note claimed the line elements were counted with the block elements; they are not), replaced by real accounting; Codex's allocation check added. Gates: `make check` debug and release, Valgrind on both decks: every fixed site absent on every rank, the master's definite losses down from 2.6 MB to 27 KB, the dipole run still nothing but Open MPI's own blocks. Two more leaks surfaced and are filed: `SideSet::initialize_lookup_tables` (F20, on every rank of every deck, now the largest BELFEM-attributed loss) and zero-length `malloc` in `Basis::set_sources` with the same count-keyed free (F21).

## Code findings, tier 4 (2026-09-18)

F20, the sideset lookup tables, on Christian's instruction with a new standing rule: astra (`gpt-6-astra`) on the Codex seat whenever the code is mesh, kernel, dof manager or their kin. Root read first: the leak is a double build, not a missing delete. The mesh-backed `SideSet` constructor builds its master and slave tables, then hands the same order to its calculator, and `Calculator::set_integration_order()` ends by calling the group's initializer again; the initializer `set_size`d the two cells and `new`ed over the live pointers, and only the enrichment cell deleted first. The destructor deleted the second generation. A factory that sets an order later builds a third one on some sidesets, and `maxwell::TMatrix` builds two at different orders. The plan jury (astra/high, Grok 4.6/high) approved the edit, one private `delete_lookup_tables()` called first in the initializer and from the destructor, and rejected four sentences of my rationale: "nobody keeps a pointer into the old generation" (the calculator does, from its last element link; the safety is that every rebuild in the tree is in setup, before the first link), "`TMatrix` calls the initializer once" (7, then 4 through the setter), "enrichment is refused at `:735`" (that line is a BFM-save abort; the flag is simply never set), and the line numbers. The code jury approved the executable change and rejected my comments: an element relink refreshes only the pointer cache, the edge-function precompute needs `link( Group * )`, and the `.cpp` comment was narration. Rewritten, then a Codex truth check qualified two more sentences (the relink is a sufficient recipe, not a necessity; the setter's caveat applies to a sideset with edge functions). The contract now sits on the public initializer and on `Calculator::set_integration_order()`.

Gates: `make check` debug and release, 17 of 17 on both (the release tree reconfigured under the peer's CMake edits, as in tier 3). Valgrind on both decks: no record cites `cl_FEM_SideSet.cpp` on any rank; the dipole's definite losses fell by exactly the finding's 12 blocks, each `tapestack3d` rank by 28 blocks. Then a gate finding of my own: the deck copies had been resuming from a `memdump.hdf5` since tier 2, the banner says "WARM RESTART, next step 2" and nothing follows it, so the tier-2, tier-3 and first tier-4 Valgrind runs covered setup, restart load and teardown, not assembly. Every site fixed in those tiers lives on setup or teardown, so their verdicts stand, but the `iv_results.csv` identity I had planned as a numerical gate was vacuous. Fresh runs with the memdump moved aside: the 3D deck's rank losses equal the restart run's (20,257 B and 20,302 B, no `cl_FEM_SideSet.cpp` record, no invalid read or write on either rank); its `iv_results.csv` differs from the tier-1 fresh run only in the two quantities of order 1e-9 to 1e-10 (`U_1` at 7.5e-8 relative, `I_2` at 1.2e-3 relative), and two native runs of the same binary differ from each other at the same level, so that deck's two-rank MUMPS solve is roundoff-nondeterministic and byte identity is no gate for it; `I_1` (0.20185037026146 A) agrees to 15 digits in all five runs. The dipole's `iv_results.csv` is byte-identical to the tier-1 fresh run, which the restarts never overwrote; no invalid read or write anywhere, the one uninitialized-value context is MUMPS' own. Gap recorded: both decks are linear, so `TMatrix`'s 7-then-4 rebuild, the one caller that changes the order, ran in no gate.

Filed from the two rounds, not touched: F22 (the setter precomputes the edge functions before it rebuilds the tables; latent while every caller passes the order the tables already have), F23 (`~Block` never frees its enrichment tables; unreachable while enrichment stays off), F24 (`SideSet::mSideSetIntegrationData` is a dead member). F21 stays filed with Codex's correction: a pointer-keyed free alone recovers only the last pair, because the guard lets a zero-count setter run again.

## Not done

- ~~Codex language sweep~~ — done (`gpt-5.6-terra`/medium, the dense-document row): 13 edits, all
  applied; three flagged sentences, two of them real errors (`fn_Smith` is one closed unit with two
  files, not two units; the TODO placeholder sentence contradicted itself) and one unclear
  ("removes nothing that lies"), all rewritten.
- **`CLAUDE.md` pointer** and the one-paragraph summary the Fable review asked for: deferred
  to after the jury for the same reason.
- No source file was touched. Reviewed, not verified: nothing here has an executable gate.

## Files updated

- `doc/commenting_guidelines.md` (new; revised after the jury; Codex-swept)
- `todo/comment_cleanup_sweep_1.md` (new), `todo/README.md` (entry)
- `.claude/scripts/ask_codex.sh`, `scripts/cross_review.sh` (`gpt-6-astra` allowlisted)
- `tmp/ai_exchange/review_commenting_guidelines.md` (ephemeral jury record)
- `doc/README.md` (entry)
- `devlog/dl20260916_commenting_guidelines.md` (this file), `devlog/README.md` (entry)
