# Devlog 2026-08-14 — Documentation Sweep Fixes + Release Packaging + FVM Move

**Date:** 2026-08-14
**Topic:** Implementation of the doxygen-sweep plan (same-day jury round, `dl20260814_doxygen_documentation_sweep.md`) under Christian's rulings: `USE_TEST=ON` is the release default, and FVM moves to `./nonfree`.
**AIs involved:** Claude (implementation); Codex + Grok jury sweep of the result (separate exchange thread)
**Claude Confidence:** high — every corrected claim was verified against the source before the doc was changed
**Literature References:** Messe et al. 2023 (DOI 10.1088/1361-6668/acf7f9, now in CITATION.cff and README)
**Verification:** executable gates on the changed tree — `check_doc_claims.py` 34/34; `update_doc_index.py --check` clean (23 modules, 85 pages); `belfem-conf drift` clean; full doxygen run (real config, scratchpad output) **0 warnings**; citation scanner 28 broken → 1 (the collaboration protocol's own illustrative example). `make doc` and `make check` in the real build tree remain Christian's runs.

## Summary

Three work streams in one session. **(1) FVM left the open tree** for `nonfree/fvm` (ruling: unfinishable before release); build hooks, CLAUDE.md, the doxygen nav and the one citation into `fvmtest.cpp` were cleaned up, and the claim checker retired its build-status probe automatically. **(2) The false documentation found by the sweep was corrected against the source**: the Maxwell guides' XML input fiction is replaced by real `input.conf` decks and the real `hphirun` driver pattern; the invented factory/postprocessor APIs (`formulation()`, `create_postprocessor()`, `process_recovery()`, `set_current()`, `kernel->set_solver()`) are gone in favor of the verified interfaces; the four fixed-bug banners ("CONFIRMED BUG" accessor, "CRITICAL BUG line 3515" weight, `is_maxwell` gap, factory switch) are rewritten as verified-historical records; the mesh I/O table names `BfmFile` instead of the nonexistent `HDF5Reader/Writer`; the sparse README's `-DBELFEM_*` flags became the real `USE_*` options with their actual defaults; `examples/README.md` names the real `tapestack3d.geo`; stale line citations in the theory docs were converted to searchable anchors. **(3) Release packaging**: root `README.md` is now a real front door (description, build, first run, `make doc`, citation, license — LBNL notice kept below), `CITATION.cff` added, `PROJECT_NUMBER = @PROJECT_VERSION@`, and a new published `doc/getting_started.md` becomes the user entry point, linked first from the mainpage.

## Key Findings (made during implementation)

- **The maxwell usage guide's fabrication went deeper than the sweep saw**: beyond XML, its §4.2/§9/§10 taught five nonexistent API families. Every replacement snippet was checked against `cl_MaxwellFactory.hpp`, `cl_MaxwellPostprocessor.hpp` and `hphirun.cpp` before writing.
- **IWG guide "Issues 3/5" were already fixed in source** (`is_maxwell` includes `MaxwellThermal`; the factory default is an explicit error), while **Issue 4 is genuinely still open** (hardcoded `cp=1.0`, `lambda=111.0` in `cl_IWG_TransientHeatConduction.cpp` — affects only the standalone heat solver). The guide now says exactly that.
- **Bycatch: four pre-existing doxygen warnings the 0-warning gate never saw.** Markdown links of the form `[README.md](README.md)` in four theory docs resolve against doxygen's *working directory*, not the doc's directory — an unresolvable `\ref` to the repo-root README, present at HEAD (verified via a clean worktree run) and invisible to the morning session's mirror. Fixed with explicit `@ref <module>_index` targets; the full-config run is now genuinely 0 warnings.
- `check_doc_claims.py`'s USE_TEST probe is now **bidirectional** (requires the docs to state whichever default CMake sets, and forbids the stale opposite); the one-way form would have gone silent on exactly this session's ON flip. CLAUDE.md, `coding_philosophy.md` (both sites incl. the corrected-against-source ledger) updated to `USE_TEST=ON` per the ruling.

## Changes Made

- Moved: `src/fvm/` → `nonfree/fvm/` (open-repo `git rm`; nonfree devlog `dl20260814_fvm_module_arrival.md`); `src/CMakeLists.txt`, `nonfree/CMakeLists.txt` adjusted.
- New: `README.md` (rewritten), `CITATION.cff`, `doc/getting_started.md`.
- Docs corrected: `src/fem/maxwell/doc/{maxwell_usage_guide,README,contact_impedance_theory,ghost_penalty_stabilization,postprocessor_recovery_theory}.md`, `src/fem/kernel/doc/{dof_manager_usage_guide,README,hanging_dofs_static_condensation}.md`, `src/fem/iwg/doc/{iwg_usage_guide,README}.md`, `src/mesh/doc/README.md`, `src/sparse/doc/README.md`, `src/comm/doc/README.md`, `src/io/doc/io_usage_guide.md`, `src/fem/interpolation/doc/nedelec.md`, `src/homology/doc/{cohomology_theory_and_implementation,thick_thin_cuts_and_conjugate_edges}.md`, `examples/README.md`, `doc/{README,mainpage,coding_philosophy}.md`, `CLAUDE.md`, `todo/instruction_doc_currency.md`.
- Tooling: `scripts/check_doc_claims.py` bidirectional USE_TEST probe; `Doxyfile.in` `PROJECT_NUMBER`; `doc/doxygen_nav.dox` regenerated.

## Jury round on the implementation (same session)

Codex + Grok reviewed the working-tree diff (`tmp/ai_exchange/review_doc_sweep_fixes.md`). **My "corrected against the tree" claim was refuted as a general statement** — the fabrications survived beyond the sections I had rewritten (the README's second XML diagram, `kernel->set_solver`, `process_recovery()/process_filter()`, the 26-variant `h_metal`/`h_hts` kernel family, three references to a nonexistent `material_models.md`, `-DBELFEM_*` flags in five more files, CN/Galerkin advertised as live) and my own §4.2 loop omitted the per-step `set_currents` (a copied driver would run at zero current). All confirmed findings in the doc mandate were fixed in-session and the `-DBELFEM_` class swept to zero; gates re-run green (34/34, nav clean, drift clean, doxygen 0 warnings, 0 past-EOF citations). Grok's residual-friend warning suspicion was refuted by the doxygen probe. In-round disclosure: my bidirectional-probe forbid regex missed backticked stale forms; found by self-verification, fixed, probed 8/8.

**Escalated to Christian; his rulings applied same session:** `test_homology` was missing from BOTH `check`/`check-fast` dependency lists (`CMakeLists.txt:300-321`) while its suite is labeled `fast` — **fixed on his instruction**: `homology` added to both foreach lists (target name `test_homology` confirmed via `Add_Test.cmake:20` + `TESTNAME homology`); `Nickel::create_alpha()` never being called is **work in progress by design** (his ruling — leave it); trailing whitespace in `cl_EF_TET4.cpp:191-193` (trips `git diff --check`) remains with his D2 work.

## Open Questions

- CITATION.cff's `preferred-citation` lists only Christian as author — the SUST co-author list was not verifiable offline; complete it from the DOI record.
- Version: CMake still says 0.9.0 against the 1.0 target; bump at tag time (README/banner/doxygen all inherit).
- The five internal/process pages (AI protocol, workflow notes, `circuit/doc/notes.md`) remain published in the doxygen output — still Christian's call.
- Owed executable gates: a real `make doc` (baseline 0 must hold there too) and `make check` on the reconfigured tree.

## Files Updated

See "Changes Made"; jury-round record in `tmp/ai_exchange/review_doc_sweep_fixes.md` (swept after distillation).
