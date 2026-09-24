# Regression cross-review of devel against v0.9.2

**Date:** 2026-09-23
**Purpose:** Overnight read-only check that the 89 commits since the 0.9.2 tag (37k insertions, 244 files, 42 of them under `src/`) introduced no severe regression: does anything abort a deck that ran, change a result silently, or open an ownership or MPI hazard? Four jury slices over the integrated `git diff v0.9.2..HEAD`, each pre-registered, audited by both seats, verified citation by citation and reconciled
**Module:** `src/fem/kernel`, `src/fem/maxwell`, `src/mesh`, `src/physics`, `src/sparse`, `src/io`, `src/containers`, `src/comm`, `src/core`, `src/homology/cl_CutData`, `CMakeLists.txt`; `tests/` read for wiring only
**AIs involved:** Claude Fable 5.1 (pre-registration, verification, reconciliation); Codex `gpt-6-astra`/high (kernel, mesh), Codex `gpt-5.6-terra`/high (physics, infra); Grok `grok-4.7`/high (all four). The first Codex dispatch of every slice failed on the ChatGPT usage limit; Christian raised the credits and the four legs were re-run at 20:56 with identical prompts, so each exchange file carries one "Audit FAILED" Codex entry from the original round beside the real one
**Status:** COMPLETE — no P0; two P1s, both release-record items (results that changed on purpose and are not in `CHANGELOG.md`); five P2s, all pre-existing defects that the new guards expose or leave uncovered; eight adjudication rows for Christian listed below. No source edited
**Verification:** reviewed, not verified by this session. Executable evidence on the reviewed state, all Christian's: `make check` 17/17 at 19:24 on `build/` and `cmake-build-debug/` (both release, MKL, Armadillo; no debug-assert run exists for HEAD); DR-188's old-vs-new A/B on all ten shipped 3-D decks (232 curves identical). Run here, read-only: `check_doc_claims.py` 38/38, `check_wrapper_policy.py` clean, `check_practice_drift.py --strict` exit 0, `check_citations.py` red (one false positive, below). Codex ran a NumPy probe of the recovery block inverse (errors < 1.4e-15) and a loop-level probe of the periodic matcher; neither is a BELFEM run
**Exchange:** `tmp/ai_exchange/review_devel092_kernel.md`, `review_devel092_mesh.md`, `review_devel092_physics.md`, `review_devel092_infra.md` (each: frozen pre-registration, both auditor entries, per-seat verification, reconciliation table); the four `review_devel092_*.diff` slices are the subjects

## Summary

The delta is healthy. Every behavioral change since 0.9.2 falls into one of three classes, and
none is a regression by the house rule (`doc/lessons_learned.md` L-18: wrong-answer state → hard
fail; retryable → status):

1. **Results that changed on purpose, towards the correct value, and are pinned by a test.** The
   prism side-face surface increment (0.25 → 0.5 on `normal_penta` cases 0–2; 0.9.2 integrated
   half the area of every PENTA quad side face, volume prisms included), the per-point background
   field increment on non-affine facets, the QUAD9/HEX20 bottom-node tables, and the jc/n
   continuation below 10 mT for every deck with a database material.
2. **Aborts that replace a silent wrong answer or undefined behavior.** Non-unit cut coefficients
   (was: zero weight), an incomplete periodic node match in release (was: ASSERT compiled out),
   a non-positive or non-finite relative tolerance, a closed thin-shell terminal curve (was: an
   out-of-range index later), a sideset facet no selected block owns (was: a write to
   `mFaces[ UINT_MAX ]`), a plugin that declared its jc/n arguments in the wrong order (was:
   silently transposed), a table whose lowest-field face rises (was: loaded and clamped).
3. **Paths that used to crash and now run.** Constant-n power-law overloads (null `mNFunction`),
   an unselected thin-shell layer block in the serial linker, a mixed thin-shell selection whose
   facet fill wrote past `set_size`, a quadratic 2-D curve with no sideset B.

The one class a 0.9.2 user will notice without reading the code is class 1, and `CHANGELOG.md
[Unreleased]` lists none of it — that is the P1.

## Key findings

**P1 — release record (both seats and Claude, all four slices agree)**
- `CHANGELOG.md [Unreleased]` names only the rutherford example and the SCLS page. Missing: the
  2× prism side-face integrals, the recovery rewrite (near-singular patches no longer abort),
  the jc/n continuation below 10 mT, the new refusals of class 2, the `$PWD` → cwd change. Also
  `CMakeLists.txt` says 0.9.3 while `CITATION.cff` says 0.9.2 (the 0.9.2 release bumped both
  together, `dl20260918_release_092.md:18-19`) — whether devel carries the next version between
  releases is Christian's ruling (kernel/infra tables).

**P2 — pre-existing defects the new guards expose or miss (single-raiser rows need Christian)**
- `cl_Mesh_PeriodicityFactory.cpp:709-738`: the match loop has no `break`; a source within
  `BELFEM_MESH_EPSILON` of two targets is counted twice, `tCount == n` passes with another source
  unmatched, and in release the over-match writes `aTargetNodes( n )` before the new
  `BELFEM_ERROR` runs. Raised independently by Codex (with a loop probe) and Grok. Fix shape:
  `break` on the hit, count distinct partners.
- `cl_CurveFactory.cpp`: the `sort_end_nodes` refusal (`:1181`) leaves every mesh node at
  `gNoIndex` and walk flags dirty when errors throw; `intersect`/`from_2d_sidesets` restore only
  on success (`:667`, `:775`), and `thin_shell_side_curves`' raw `Curve` (`:204`) is outside
  `tReleaseScratch`. Release aborts; debug and tests see the dirty mesh. Both seats. Fix shape:
  the scope guard the side-curve path already has, in all three builders.
- `cl_FEM_DofMgr_BlockData.cpp:479-481`: the serial Maxwell side-connector loop still indexes
  `mBlockMap` unguarded for an unselected shell (the layer-block loop above it is now guarded;
  the parallel path guards both). Codex. Side connectors are an experimental path.
- `cl_JcFunction_Database.hpp:190-199`: the face-slope certificate is demanded of clamp-only jc
  tables (no `/source`), which never continue below Bmin; a user-built table that loaded in
  0.9.2 is refused for a continuation it will not take. Grok. Shipped tables pass.
- `cl_CutData.cpp`: a caller that catches the new non-unit refusal and calls
  `collect_coefficients` again overwrites the two bitset pointers (pre-existing double-call
  hazard, no production catcher). Grok.

**Adjudication rows (design or policy, not settled by vote)**
1. CutData non-unit coefficient: abort (now) vs warn — Christian/Gregory.
2. Version 0.9.3 in devel before release: intended practice or a slip?
3. Periodic matcher fix (P2 above) — approve the `break` + distinct-count shape?
4. Curve-builder scope guard (P2 above) — approve?
5. Serial side-connector guard (P2 above) — approve, or leave with the experimental path?
6. Slope certificate only where the label continues (P2 above) — approve the reorder?
7. `<= 1e-10` accepted under a message that says "certified not to rise": Codex P2, Claude P3.
8. `Postprocessor` implicitly copyable with two owning maps (Grok, P3): delete copy/move as
   `Kernel` and `Basis` got?

**Corrections to my own pre-registration, accepted from the auditors**
- My one candidate defect (a one-row `/source` curve reading out of bounds in
  `finalize_self_field`) was refuted by both seats: `reduce_self_field_rows` returns before
  exporting fewer than two bins (`cl_JcFunction_Database.hpp:471-472`). Withdrawn.
- The changed prism faces are `normal_penta` cases 0–2, not 1–3, and the direct area test is
  `test_InterfaceOrientation.cpp:2852`; the doubling reaches volume prisms, not only thin shells.
- 0.9.2's null dereference of sideset B hit quadratic 2-D curves only (short-circuit `&&`).
- The DofMgr facet fix upgrades from "abort → runs" to "out-of-bounds write → runs".

**Already ruled, surfaced again by Grok as a refutation:** the closed-terminal refusal in
`MaxwellFactory` makes the `is_closed() && tNumTips == 0` disjunct at `cl_CutFactory.cpp:1245`
unreachable — recorded and ruled on 2026-09-22 (keep the refusal, dead branch noted for
Gregory). Not new.

**P3 and wording (in the tables, no ruling needed):** absolute SVD floor absent from
`postprocessor_recovery_theory.md` §3.2; "exact derivative of eval" briefs vs the admitted
tangent jump across β = 2; the Bernstein bound header vs its scaled return and the n-floor message
printing a bound as an attained value; the 77.5 K calibration scaling the raw linear curve while
runtime uses the projected Hermite; `core_usage_guide.md:1282-1291` stale after the globals
comment; `(int)` narrowing in the MUMPS unknown-bit message under `BELFEM_INT64`;
`NodeFrame` fields unprefixed; a `key_t` map miss prints "Key unknown".

## Tests and tooling

- Every new test file is wired into a suite CMakeLists; no pre-existing assertion was removed
  from a modified test; the `DISABLED_` count rose 7 → 16, all nine new pins documenting known
  register rows (lead cp, YBCO λ, Ni/Fe λ, Kim derivative, power-law tangent, pyramid auto-order),
  none silencing an old test.
- `check_citations.py` exits 1 on HEAD: it reads "Amazon Linux 2023" in `doc/getting_started.md:27`
  and `doc/scls_toolchain.md:63` as an author-year citation. Checker false positive, red gate.
- Gaps in the executable evidence: no debug-assert `make check` exists for HEAD (both trees are
  release); the 2-D curve path (`from_2d_sidesets`) is covered by unit tests but by no deck A/B;
  which shipped decks integrate on prism side faces is not enumerated; the Riva/Denis/Messe
  equation numbers in the jc comment were checked by neither seat.

## Process notes

- Two juries ran with only the Grok seat until the credit top-up; the Codex re-runs used the same
  priming and subject as `cross_review.sh` builds, appended under the same slug and stamped by
  the wrapper, so each slice's record has both seats at the intended depth.
- A `run_in_background` Bash task was not killed at its 600 s cap (a Grok leg ran 17 min inside
  one); the later legs were launched with `setsid nohup` anyway.
- Coordinated with the parallel debt-register session (belfem-de, Christian's
  `session_01PGXN8y…`): it held all `src/` and `tests/` edits until the literal "REVIEW DONE" at
  21:2x, received the per-file list above, and owns the trees from here. `git status` after every
  auditor round showed only that session's `todo/` files.

## Files updated

- `devlog/dl20260923_regression_review_devel_vs_092.md`, `devlog/README.md` (this entry)
- `tmp/ai_exchange/review_devel092_{kernel,mesh,physics,infra}.md` and `.diff`, the saved Codex
  prompts and the two stream logs (ephemeral)

No source, test, register or todo file was changed by this session.
