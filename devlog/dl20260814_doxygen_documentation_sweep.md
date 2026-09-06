# Devlog 2026-08-14 — Three-AI Sweep of the Doxygen Documentation

**Date:** 2026-08-14
**Topic:** Jury-mode cross review of the generated doxygen documentation (`cmake-build-debug/doc`) and its markdown sources, against three questions: currency vs the code, gaps vs open-source best practice, and the plan.
**AIs involved:** Claude (pre-registration + verification), Grok (substantive audit), Codex (audited the charter — see round shape)
**Claude Confidence:** high on the verified facts; severity labels are judgment
**Grok Audit Confidence:** high on its P0/P1 file:line facts
**Literature References:** N/A (documentation currency), except that ghost-penalty prose was re-checked against Burman & Zunino 2006 as cited in the module doc
**Verification:** two probes (a doxygen coverage run with `WARN_IF_UNDOCUMENTED=YES` into the scratchpad; a mechanical citation scan of all 90 md files) plus source trace on every auditor citation. `make doc` was **not** run; no repository or build-tree file was modified. Reviewed, not verified, for anything requiring a doc build.

## Summary

The documentation splits cleanly in two. The **contract layer is current and unusually strong** for a research code — the `input.conf` reference and schema, the executables guide, the literature list with DOIs, the FVM not-built disclosure, plus a machine-checked input contract (`belfem-conf`) and a doxygen warning baseline of 0. The **module usage guides are the weak half**, and the weakness is not the citation rot I set out to measure: it is prose that describes a superseded system. The Maxwell module documentation still teaches an **XML input file** while the executables have taken `input.conf` for a long time.

The other half of the answer is packaging. Internally this project is over-documented; publicly it is close to invisible. The repository front door is a 16-line copyright notice, `examples/` never reaches the generated site, there is no citation metadata, the version is 0.9.0, and there are zero `install()` rules.

## Key Findings

**Round shape (a defect I own).** `cross_review.sh` treats a `[path]` argument as *the artifact under review*, so Codex audited the charter file instead of the documentation. The round is therefore one substantive audit (Grok) plus a charter audit (Codex) plus my pre-registration. Future non-diff sweeps must carry the brief in the prompt, not in a file passed as the subject.

**The correction to my own pre-registration.** I claimed the citation rot was line drift rather than documentation of vanished code, and that no document was substantively wrong. Grok refuted the general form of that: `maxwell_usage_guide.md` has *clean* citation statistics and is the worst document in the tree. **Citation health and prose health are uncorrelated**, so a mechanical citation checker — the instrument I proposed — would never have found the XML defect. That reshapes the plan away from citation tooling and toward reading the oldest guides.

**The adjudication that mattered.** Grok reported that `CMakeLists.txt:81` sets `USE_TEST` ON while `CLAUDE.md` and `coding_philosophy.md` say OFF, and called the docs false. Verification refuted it as written: the ON is an **uncommitted working-tree edit**; `git show HEAD:CMakeLists.txt` line 81 is OFF, so the documents match the committed tree. What survives is real and latent — `scripts/check_doc_claims.py:224-226` only requires the docs to say OFF *if* CMake says OFF, so the gate goes silent exactly when the documents would become wrong.

**Confirmed P1 prose defects** (each verified against the source): Maxwell guides teach XML input (`maxwell_usage_guide.md:224-229`, `README.md:123,384`) vs `hphirun.cpp:57`; two files still carry "CONFIRMED BUG" / "CRITICAL BUG" banners for fixed code (`dof_manager_usage_guide.md:50`, `hanging_dofs_static_condensation.md:990` — the code now does exactly what those docs prescribe as the fix); `src/mesh/doc/README.md:93-94` lists `HDF5Reader`/`HDF5Writer`, which do not exist (live type `BfmFile`); `src/sparse/doc/README.md:49-53` documents `-DBELFEM_STRUMPACK=ON` when the option is `USE_STRUMPACK`; `examples/README.md:63-64` names `circuit.geo` and `sidecoating.geo` when both directories ship `tapestack3d.geo`, breaking the documented first run.

**Measured, not asserted.** 28 of 507 `file:line` citations are provably broken (8 absent files, 20 past EOF), concentrated in historical change narrative that belongs in devlogs — `iwg_usage_guide.md` is 11/13 bad almost entirely inside its "Known Issues and Critical Bugs" section. A doxygen run with `WARN_IF_UNDOCUMENTED=YES` reports **223 of 327 classes (68%) and 8,206 members undocumented**; the shipped config sets `WARN_IF_UNDOCUMENTED=NO` with `EXTRACT_ALL=YES`, so the 0-warning baseline measures malformed documentation only and must not be read as coverage.

## Changes Made / Proposed

None. This was a review round; per protocol no fixes were applied. The proposed plan is recorded in the exchange file and summarized to Christian, ordered by "what makes a stranger fail in the first hour" (front door, first run, false claims) ahead of completeness work.

## Open Questions

- **`USE_TEST` intent** — is the working-tree ON deliberate (heading for a committed default change) or leftover from the 2026-08-13 suite run? If it is ever committed, `check_doc_claims.py` must become bidirectional first.
- **Publishing internal process docs** — `ai_collaboration_protocol.md`, `ai_workflow_best_practices.md`, `notes.md`, `fvm_summary_old.md`, `fvm_work_in_progress.md` are all in the public doxygen output. Christian's call, not a reviewer vote.
- **Version and tag** — CMake says 0.9.0 against a 1.0 target; CMake, the generated banner, `PROJECT_NUMBER` and the README need to agree.

## Files Updated

- tmp/ai_exchange/review_doxygen_doc_sweep.md (pre-registration, both audits, verification, reconciliation)
- devlog/dl20260814_doxygen_documentation_sweep.md (this file)
- devlog/README.md (index)
