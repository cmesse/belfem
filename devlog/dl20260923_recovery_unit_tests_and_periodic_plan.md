# Recovery unit tests through a real kernel, a round-off floor in the block inverse, and the periodic-fields plan

**Date:** 2026-09-23
**Purpose:** R9 of `todo/thinshell_nodal_recovery.md`: tests that run the node-local patch recovery
on code-built meshes through a real `Kernel`/`DofManager`, each pinned by a mutation. On the way,
the block inversion was extracted into a testable function, and the round-2 jury found a round-off
defect in it (present since the recovery fix of 2026-09-22), now fixed. Separately: a plan, audited,
for making the recovered fields periodic.
**Module:** `src/fem/kernel` (`cl_FEM_Postprocessor`, new `fn_FEM_recovery_inverse.hpp`), `tests/fem`,
`src/fem/maxwell/doc`
**AIs involved:** Claude (plans, implementation, gates). R9 plan jury: Codex `gpt-6-astra`/high, Grok
`grok-4.7`/high (80 turns). R9 code round 1: Codex astra/high, Grok 4.7/high (80 turns). R9 code
round 2: Codex astra/**xhigh**, Grok 4.7/high — failed at 80 turns, rerun alone at **150** turns.
Periodic plan jury: Codex astra/high, Grok 4.7/high (80 turns). Explore subagent for the periodicity
trace.
**Status:** R9 implemented; periodic plan audited, not implemented. Not yet committed at the time of
writing (Christian runs `make check-fast` in the shared trees first).
**Verification:** scratch worktree at HEAD `42b5e42` + the recovery change + this work, own build
directories (release and debug, Blaze/MKL, `-Werror`); the Armadillo tree syntax-gated only.
11/11 new tests; full `test_fem` 210 release / 212 debug; mutations M1–M5 each red on their own
tests only (table below). The shared trees are not built by me.
**Exchange:** `tmp/ai_exchange/review_recovery_unit_test.md` (+ `_r2.diff`),
`tmp/ai_exchange/review_periodic_fields.md`.

## Summary

`tests/fem/test_PostprocessorRecovery.cpp` (11 tests) runs the real recovery: a test
postprocessor overrides `compute( k )` with a manufactured field at the integration point. Fixtures:
a two-layer PENTA6TS fan on a cylinder (R = 1, t = 2e-3, sagitta over the interior patch 7.6 t), a
2 × 2 × 2 TET4/TET10 cube patch 1e3 from the origin (optionally split by co-located duplicates on
one plane, the cut pattern), and a 1 × 1 × 1e-6 flat TET4. Three direct tests call the extracted
`recovery_block_inverse`.

| mutation | red, and nothing else |
|---|---|
| M1 shell mode off | `ThicknessJumpSurvivesCurvature` |
| M2 centering removed (absolute basis, evaluated at the node) | `LinearFieldFarFromOrigin`, `QuadraticFieldFarFromOrigin` |
| M3 full-matrix truncated pseudo-inverse | `BlockInverseKeepsCoupledConstant`, `BlockInverseDropsSchurRoundOff` |
| M4 siblings left out of frame, Gram and crediting | `SiblingsShareTheUnionFit` |
| M5 round-off floor removed | `BlockInverseDropsSchurRoundOff` |

## What the fixtures needed (both not obvious)

- `Mesh( 3, 0, false )` — the pattern every hand-made fixture uses so `MeshChecker` can run before
  edges exist — also skips `connect_nodes_to_elements()` in `finalize()`. The recovery walks
  node → element; the first run segfaulted in `compute_node_matrices`. Each fixture now calls
  `mesh::ConnectivityCalculator::connect_nodes_to_elements()` after the checker. Not a production
  defect (real meshes are built with connectivity; medium-high, not traced exhaustively).
- Shell volumes are measured by `element id − first id` of the block (`cl_FEM_Kernel.cpp:1026`) and
  facets linked by container position (`cl_FEM_DofMgr_BlockData.cpp:469-472`): wedges inserted in
  tape-facet order with consecutive ids per block (Grok, plan round). The debug assert on the
  ordinal passed.

## What the juries changed

**Plan round (Codex, Grok; all verified):** the parent's R9 also named a quadratic case and the
sibling union (T7, T8 added); M2 must remove the centering consistently, or a red T4 shows a broken
evaluation shortcut rather than conditioning; the degenerate case must couple the constant to the
dropped direction or M3 cannot fail; the calculator is allocated only in `DofManager::initialize`;
F2's κ is (L/h)² ≈ 4e12, not ⁴.

**Code round 1:** Grok's headline — node owners unset, nothing recovered, field stays 0 — was
**refuted by execution** (T2 expects 3.7 everywhere and passes; mutations turn tests red). Confirmed:
T8 passed under side-local fits (a linear field is exact on each half), T6 could not discriminate
M3, the serial-only skip was missing, T7 checked corners only. And my own measurement: **M3 stayed
green through every valid element**, because a node lies in its own patch, so the constant couples
to a dropped direction only as much as the element is flat. Christian: extract the inversion and
test it directly.

**Code round 2 (Codex xhigh; Grok rerun at 150 turns agreeing independently):** the Schur
complement of an exactly singular gradient block keeps round-off (−2.2e-19 in Codex's case,
v = (1, 0.1, 0, 0), w = 0.1), and the cut σ_max/1e10 is relative to that residual, so the noise
was inverted: a constant 3.7 came back as **4.0**. Present in the recovery since 2026-09-22 (the
extraction moved it unchanged); it needs a vanishing Schur block, which no rutherford patch
approached (κ ≤ 12), so the 2026-09-22 A/B stands. Fix: the cut is also floored at
64 ε · max diag(Gᵣᵣ), the scale of the subtracted terms (a healthy patch still inverts exactly to
1.4e-17); `BlockInverseDropsSchurRoundOff` pins it (M5). T8 now compares both siblings with the
unsplit patch. Grok's nits applied: the `aCondition` doc no longer overclaims, a debug assert on the
truncation limit. Left open: Grok's unmeasured doubt (~40 %) that 64 ulps covers round-off summed
over hundreds of rank-1 updates.

## Also this session

- The per-rank line "recovery on proc …: patches, largest condition number" was a probe for
  Christian's κ question (answered: ≤ 15 on 8 ranks); removed with its member at his request. The
  truncation warning stays. `postprocessor_recovery_theory.md` §3.2 updated to match.
- Periodic fields (`todo/periodic_postprocessor_fields.md`): rank 0 combines each master/slave node
  pair after the gather, vectors through R = T_tᵀT_s rebuilt from the plane nodes, **plain mean**
  (Christian). The jury corrected the gate deck (the twisted one is `corc_garber_effect`, not
  `corc_periodic_bc`, which I had taken over from the trace), restricted R's validation to surface
  pairs (layer nodes sit off the plane), made a failed transform a setup error, added the one-sided
  cut duplicates (no periodic pair of their own) and the rank-0 membership rule, and required
  real-pipeline and MPI tests. Medium effort; not implemented.
- Grok's turn cap: two jury legs today needed more than 80 turns; one rerun at 150 delivered.

## Open

- T5 (2D QUAD4TS shell) not built; G5 of the parent plan open.
- `make check-fast` in `build/` and `cmake-build-debug/` (Christian) before the commit.
- The scratch worktree `scratchpad/ab/src` is still registered with git.
