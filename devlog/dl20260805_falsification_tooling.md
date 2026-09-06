# Devlog 2026-08-05 — Falsification Tooling

**Date:** 2026-08-05
**Topic:** Interface/orientation regression battery + `make check-fast`, evidence
hierarchy in the protocol, campaign pages + debt register. Plan:
`todo/falsification_tooling.md` (all O-decisions resolved by Christian).
**AIs involved:** Claude (implementation, builds, tests — newly sanctioned)
**Claude Confidence:** high on everything executable below; battery physics values
are proposals pending Christian's sign-off
**Literature References:** Messe et al. 2023 (paper1) §2.7 (controller context only)
**Verification:** focused regression — `make check-fast` green on BOTH backends
(`cmake-build-arma-test`, `cmake-build-blaze-test`, working tree @ `a24656cb` + this
session's files) incl. from `env -i`; mutation check red/green (below)

## Summary

The falsification side is now tooling. `tests/fem/test_InterfaceOrientation.cpp` locks
the sign/orientation/master-slave/indexing bug class behind a minutes-scale gate
(`make check-fast`: 4 modules, <1 s test time after build, both backends, bare-env
capable). The evidence hierarchy ("reviewed" ≠ "verified", 7-level ladder, review stop
condition, same-session campaign/register rule) is in the protocol §11, the
`/cross-review` reconciliation table, and the devlog template. Six campaign pages and a
45+3-row debt register seed the compression layer (everything `[seeded — confirm]`).

## The battery (10 active tests, 2 documented scaffolds)

Fixtures (`tests/fem/support/cl_TS_TestStack.hpp`): code-built micro-meshes under a
minimal REAL kernel chain (Mesh → KernelParameters → Kernel → DofManagerBase → Block →
`fem::Element` aura ctor + `set_facet`) — O1(b) worked with zero production-code
changes; the designed thickness data path is the one under test. Edges are fixture-owned
(MeshChecker rejects pre-Kernel edges) and built canonical from `get_nodes_of_edge`.

| Test | What it locks |
|---|---|
| Quad4TsUnitCirculation | ∫edge_j E_k·dl = s_k δ_jk, both edges along +tape, rotated variant |
| Quad4TsEdgeDirectionSign | flipped mesh edge flips exactly its own dof |
| Quad4TsFaceActivity | bottom column vanishes at η=+1 and vice versa |
| Quad4TsStokesConsistency | ∮E_k·dl = C_k·Area — catches ONE-SIDED E/C flips |
| Quad4TsInterlayerContinuity | shared-edge tangential trace equal from above/below (the greg2 mechanism) |
| Quad4TsAmpereTelescope | N=4 stack: uniform trace ⇒ I_l=0 per layer; ΣI_l = outer-trace difference |
| Penta6TsUnitCirculation / Hex8TsUnitCirculation | 6×6 / 8×8 circulation matrices = identity |
| Penta6TsEdgeFlipSign | flip diagonal isolation, 3-D |
| TopBottomNodeAlignment | `get_top_nodes(k)` vertically above `get_bottom_nodes(k)` (43474c9f class), 2-D rotated + both prisms |
| DISABLED_GhostElementContract | h_ghost 12-dof contract — needs Calculator fixture (DR-46) |
| DISABLED_Hex8TbUnitCirculation | blocked: NO Lagrange factory case for HEX8TB (cl_IF_InterpolationFunctionFactory.cpp:215) — first executable evidence of the hex8tb_phase2 R4 gap (DR-47) |

## Mutation check (acceptance)

```
sed -i '166,167s/mS\[ 1 \]/-mS[ 1 ]/' src/fem/interpolation/nedelec/cl_EF_QUAD4TS.cpp
make -C cmake-build-arma-test test_fem && ./test/test_fem --gtest_filter='InterfaceOrientation.*'
  → 4 FAILED: UnitCirculation, EdgeDirectionSign, StokesConsistency, InterlayerContinuity
git checkout src/fem/interpolation/nedelec/cl_EF_QUAD4TS.cpp && rebuild
  → 10 PASSED
```
The C-only telescope stays green under an E-only mutation — by design the Stokes test
carries the E↔C cross-check. The battery catches the bug class it was built for.

## EXPECTED values for Christian's sign-off (acceptance 3)

Flagged `// EXPECTED: pending Christian sign-off` in the test file; on sign-off each
flag is replaced by a one-line provenance comment (O6).

1. **Unit circulation convention (QUAD4TS):** ∫edge_j E_k·dl = s_k δ_jk with BOTH
   bottom and top circulating along +tape (the 4a42d982 convention). Derivation:
   E_k = s_k f_k(η) ∇ξ, ∇ξ = t̂/L, f_0 = (1−η)/2, f_1 = (1+η)/2.
2. **Per-layer Ampère (QUAD4TS):** I_l = ∫C·q dA = h_l − h_{l+1} (bottom minus top
   trace) under curl orientation (∇η×∇ξ)_z; uniform trace ⇒ I_l = 0 per layer;
   telescope ΣI_l = h_0 − h_N. Derivation: C = (s_0, −s_1)/(tL), Area = tL
   (rotation-invariant: (∇η×∇ξ)_z = −2/(tL) for any tape angle).
3. **3-D circulation matrices:** PENTA6TS (6 dofs) and HEX8TS (8 dofs) circulation
   matrix over canonical (`get_nodes_of_edge`) edges = identity, dof column k ↔
   element edge k. Empirically confirmed against the current implementation; the
   convention itself needs the sign-off.
4. **HEX8TB circulation = identity over the four longitudinal edges** — NOT yet
   executable (DR-47); value proposed from the h_t stream-function design.

Structural expectations (δ-shape, face activity, Stokes identity, flip isolation,
node-tie alignment) are geometry/consistency facts, not flagged.

## check-fast

`TESTLABELS` support in `Add_Test.cmake`; fast label on containers/linalg/mesh/fem
(math EXCLUDED — see DR-48); root target `check-fast` = build those + `ctest -L fast`.
Measured: 0.93 s (Armadillo) / 0.50 s (Blaze) test time; green from `env -i` through
`scripts/scls_env.sh`. New build trees `cmake-build-arma-test`, `cmake-build-blaze-test`
(Christian's `cmake-build-debug` untouched; its `USE_TEST` stays OFF).

## Defects found by building the battery

- **DR-47:** HEX8TB has no Lagrange interpolation-factory case — kernel chain over a
  HEX8TB block throws; = hex8tb_phase2 R4, now with a reproducer.
- **DR-48:** `tests/math/test_Quaternion.cpp` targets a removed borrowed-buffer
  Quaternion API — `test_math` does not compile at all (latent since USE_TEST=OFF);
  needs reconcile-or-prune (Christian).

## D2 / D3

- Protocol §11 (evidence hierarchy, stop condition, same-session rule), devlog template
  `**Verification:**` field, `/cross-review` evidence column + stop condition — diffs
  shown in-session for approval before commit.
- `devlog/campaigns/`: side_connector_wall_element, 2d_thinshell_validation,
  controller_anderson, bfm_persistence, gasmodels_migration, release_1.0 (periodic +
  double-corc folded in as dormant). ALL claims `[seeded — confirm]`.
- `todo/debt_register.md`: DR-01…DR-48, blocking-1.0 column proposed; seeded in one
  pass from the 56 Open-items sections of the June+ devlogs (timebox held).

## Addendum (same day): DR-48 closed — math tests reconciled

Christian asked for the math-test fix. Three distinct stale layers plus one real
source defect, all resolved; `test_math` = 368/368 on BOTH backends, restored to the
fast set (`check-fast`: 5 modules, Armadillo 0.70 s / Blaze 0.65 s, `env -i` pass).

1. **test_Quaternion.cpp** — eight tests targeted the removed external-buffer mode
   (the header documents the removal, cl_Quaternion.hpp:29-33). §1.2 is now a
   value-type-contract section (`is_trivially_copyable` static_assert + memcpy
   round-trip = the MPI contract; per-object inline storage); the move test now
   asserts rule-of-zero semantics (the old test expected the source pointer nulled);
   conj/inv/binary "returns owning" tests rewritten value-mode; the three
   `*WithExternalBuffer` tests are preserved in an `#if 0` reference block.
2. **test_Tensor.cpp** — `TensorDebug` expected 3x3x3x3-only addition to throw on
   equal 2x2x2x2 operands; the operators were generalized to any MATCHING shapes
   (per-dimension asserts). Tests now check mismatch-throws + matching-legal.
3. **test_TensorKernels.cpp** — the three `*VsReferenceLoop` tests indexed
   `Matrix::data()` linearly: correct under Armadillo, WRONG under Blaze (padded
   columns — the CLAUDE.md landmine). Fixed with a `dense33` accessor-staging helper;
   values unchanged under Armadillo.
4. **Source defect (Blaze-latent):** `fn_ddot.hpp:47` and
   `fn_kelvin_christoffel.hpp:39` called the all-pointer kernel signatures that only
   the Armadillo backend headers provide; under Blaze nothing compiled against them
   (only these tests exercised the wrappers). Both now call the canonical mixed
   signature (tensor data pointer + Matrix/Vector refs) that `cl_Tensor.hpp` already
   uses under both backends.
5. **test_Spline.cpp** — `make_x` moved inside the `BELFEM_SUITESPARSE` guard its
   only consumers live in (the unused-function -Werror under SUITESPARSE=OFF). New
   finding DR-49: the entire spline test section is compiled out in the current
   configurations — spline fast-gate coverage is zero.

**Verification:** focused regression — `make check-fast` green on both backends
including from `env -i`; `test_math` 368/368.

## Open items

- Christian: correction pass on campaign pages + register (acceptance 4); EXPECTED
  sign-off block above (acceptance 3); D2 diffs approval; then commit slicing.
- DR-46/47/48 as registered; greg2 smoke (DR-15) unchanged.

## Files Updated

tests/fem/support/cl_TS_TestStack.hpp, tests/fem/test_InterfaceOrientation.cpp,
tests/fem/CMakeLists.txt, tests/{containers,linalg,mesh}/CMakeLists.txt (labels),
config/scripts/Add_Test.cmake, CMakeLists.txt (check-fast),
doc/ai_collaboration_protocol.md, .claude/commands/cross-review.md,
devlog/campaigns/*.md (6 new), todo/debt_register.md, todo/falsification_tooling.md,
.gitignore (cmake-build-*)
