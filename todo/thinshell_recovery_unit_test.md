# Unit Tests for the Node-Local Patch Recovery (R9 of thinshell_nodal_recovery)

**Date:** 2026-09-23
**Purpose:** Pin the node-local recovery of `Postprocessor` with tests that run the real fit on
code-built meshes through a real `Kernel` and `DofManager`: a curved PENTA6TS shell must keep a
through-thickness jump, a far-from-origin volume patch must reproduce a linear field, and the
inversion must reproduce a constant on a degenerate patch. Mechanism: a test subclass of
`Postprocessor` overrides `compute( k )` with a manufactured field evaluated at the integration
point, runs `initialize()` and `run()`, and reads the recovered nodal field.
**Module:** `tests/fem` (+ possibly `src/fem/kernel`, see O1)
**AIs involved:** Claude (plan, implementation), Codex `gpt-6-astra` + Grok `grok-4.7` (juries)
**Status:** IN PROGRESS — implemented; round-2 code jury: Codex astra/xhigh (2 findings, both fixed
with approval: a round-off floor in `recovery_block_inverse`, and T8 checked against the unsplit
patch), Grok rerun at 150 turns pending. Gate (scratch worktree, Blaze/MKL, `-Werror`; Armadillo
syntax-gated): 11/11 new tests; full `test_fem` 210 release / 212 debug; mutations each red on their
own tests only — M1 → T1, M2 → T4/T7, M3 → both coupled-constant tests, M4 → T8, M5 (floor removed)
→ `BlockInverseDropsSchurRoundOff`. T5 (2D) not built. Shared trees: Christian's `make check-fast`.

> **Scope guards:**
> - Tests only, plus at most the O1 extraction; the recovery itself is not changed.
> - Serial only (the fixtures never partition); the MPI behavior is covered by G4 of the parent plan.
> - Christian's ruling (2026-09-22): a DofManager-based thin-shell fixture, not free functions
>   alone. O1 asks whether the one piece a fixture cannot reach (the degenerate-patch inversion)
>   may be extracted.

---

## 1. What Is Untested Today

No test runs `Postprocessor::recover_fields` on a thin shell; the only evidence for the shell
basis is the rutherford deck (parent plan G1–G4, clean A/B 2026-09-22). `tests/fem/support/
cl_TS_TestStack.hpp` builds a `DofManagerBase` without an IWG, which `Postprocessor::set_field`
cannot take (`iwg()`, block calculators). `tests/fem/test_ThinShellSelection.cpp` (`TwoShellStacks`)
builds exactly what is needed — PENTA6TS layer blocks with tape facets and `ThinShell` records,
edges, a `Kernel`, an `IWG_Maxwell`, and `create_field()` → `DofManager` — and is the template.

| Behavior | Covered by | Gap |
|---|---|---|
| shell basis keeps the jump on a curved patch | rutherford deck only | T1 |
| constant reproduced (intercept exact) | replica in Python only | T2, T6 |
| in-plane linear field on a curved shell | none | T3 |
| centering: linear field far from the origin | none | T4 |
| 2D shell basis (QUAD4TS) | none | T5 |
| truncated inversion on a degenerate patch | Python replica only | T6 (needs O1) |

## 2. Fixtures

**F1 `CurvedShellFan`** (3D): a strip of the cylinder r = R around the z axis, 5 angular stations
over 40° × 3 z rows, split into TRI3 mid-surface facets; one PENTA6TS layer per facet, bottom
nodes at r = R − t/2, top at r = R + t/2 (exact radial extrusion), t = 2e-3 R (the rutherford
ratio 1.6 µm / 0.9 mm). One layer block, domain `ThinShell`, `set_thickness( t )`; one tape
sideset holding the mid-surface facets and a `ThinShell` record, following `TwoShellStacks`
(ghost sideset empty, single layer). Edges created per element before the `Kernel`, as there.
IWG: `IWG_Maxwell` HPhi 3D with the one block as `ThinShell`. The interior node of the strip has a
patch of six wedges whose sagitta (~R·(20°)²/8 ≈ 0.015 R) is ~8 × t: the Cartesian basis fails
there by construction.

**F2 `FarTetPatch`** (3D): 2 × 2 × 2 cubes of edge h = 1e-3, each split into six TET4, translated
to (1e3, −2e3, 5e2): the old absolute order-1 basis had κ ≈ (L/h)² ≈ 4e12 there (Grok: the power is 2, not 4), above the 1e10 truncation cut. One block, domain `Air`,
IWG_Maxwell HPhi 3D.

**F3 `CurvedShell2D`** (2D, optional, O2): an arc of radius R, 5 QUAD4TS elements, thickness t,
IWG_Maxwell HPhi 2D with the block as `ThinShell`.

**Test postprocessor** `ManufacturedRecovery : public fem::Postprocessor`: constructor sets the
block ids and one target field `"f"`; `update_element_dofs()` does nothing; `compute( k )`
evaluates a `std::function< real( const Vector<real> & x, real zeta ) >` at the physical point
`trans( node coords ) * Nvec( k )` and at `zeta = integration()->points()( dim-1, k )`.

## 3. Ordered Steps

- [x] **R1** `tests/fem/test_PostprocessorRecovery.cpp` with F1, F2, the test postprocessor, and
  registration in the `USE_MAXWELL` block of `tests/fem/CMakeLists.txt` (one line).
- [x] **R2** T1 (F1) `ThicknessJumpSurvivesCurvature`: f = ½(1+ζ)·t (0 on the bottom face, t on
  the top face). Every bottom-face node recovers 0 and every top-face node t, to 1e-12·t. A
  companion assertion records why the test discriminates: the interior patch's sagitta exceeds 5 t.
- [x] **R3** T2 (F1) `ConstantReproduced`: f = 3.7 → 3.7 at every node to 1e-13 relative.
- [x] **R4** T3 (F1) `AxialLinearFieldReproduced`: f = 2 + 5 z (z is tangent to every node's
  tangent plane on a cylinder) → exact to 1e-12 at every node.
- [x] **R5** T4 (F2) `LinearFieldFarFromOrigin`: f = 1 + 2(x−x₀) − 3(y−y₀) + 4(z−z₀) with
  (x₀, y₀, z₀) the translation → exact to 1e-10 at every node, interior and boundary.
- [ ] **R6** T5 (F3, optional per O2) 2D analogue of T1. *Not built this round.*
- [x] **R7** T6 per O1 — *and T9 (2026-09-23): the inversion extracted into `fn_FEM_recovery_inverse.hpp` (Christian: "extract + unit test"), tested directly on a healthy, a coupled rank-1 and a coplanar Gram matrix; M3 now turns the coupled case red.* Earlier note: *Built through the fixture (B3): a 1 × 1 × 1e-6 TET4 — passes, truncates every fit, keeps the constant and the in-plane gradient; but does NOT discriminate M3 (see R8).*
- [ ] **R8** Gate: `test_fem` green in both trees (`build/` Armadillo, `cmake-build-debug/` Blaze);
  named mutations turn their own test red, built and run in the scratch A/B worktree (Christian
  runs the shared trees): M1 force `mShellLocal = false` → T1 red, T2/T3 green; M2 drop the
  centering (`mOffset` = absolute coordinates) → T4 red; M3 route the constant through the
  truncation (old full pseudo-inverse) → T6 red.
  *Result 2026-09-23:* M1 ✔ (T1 red only), M2 ✔ (T4, T7, T8 red), **M3 ✘ (all green):** with valid elements the node lies in its own patch, so the dropped direction's component along the constant scales with the flatness — the coupled case needs samples on a plane that misses the node (a one-point rule), which Maxwell's calculators never use.


### 3.1 Amendments from the plan jury (2026-09-23, `tmp/ai_exchange/review_recovery_unit_test.md`)

- [x] **B1** (Codex P1, Grok) The parent's R9 also named a quadratic Cartesian case and the sibling
  union. Add **T7** `QuadraticFieldFarFromOrigin`: F2 built from TET10 (mesh order 2 → `poly2_3d`),
  a full quadratic field exact at every node; and **T8** `SiblingsShareTheUnionFit`: F2 with one
  internal plane of co-located duplicate nodes (`add_duplicate` / `set_original`, the cut pattern),
  a linear field exact and identical on both siblings. If the Kernel refuses the hand-made
  duplicates, T8 is deferred with the reason and R9 stays open for it.
- [x] **B2** (Codex) M2 must remove the centering *consistently* (absolute coordinates in the fit and
  the absolute polynomial evaluated at the node in phase 2), so a red T4 shows the conditioning
  failure, not a broken evaluation shortcut.
- [◐] **B3** (Codex, Grok) T6 through the fixture first: a positive-volume TET4 of 1 × 1 × 1e-6
  (Codex probe: gradient-block κ ≈ 2e12, one direction dropped) at a corner node, whose samples lie
  away from the node so the constant is coupled to the gradient. Extraction (O1) only if the Kernel
  refuses it.
- [x] **B4** (Codex, Grok) T6 asserts: constant exact (the coupled case, where the full pseudo-inverse
  returns 1.85 for 3.7), and the truncation count ≥ 1. The count is private (`mNumTruncated`): a
  const accessor `number_of_truncated_patches()` on `Postprocessor` is the one source addition.
- [x] **B5** (Grok) F1 ordering contract: wedges inserted into the layer block in tape-facet order
  AND with consecutive element ids starting at `elements()(0)->id()` — `Kernel` measures shell
  volumes by `id − first id` (`cl_FEM_Kernel.cpp:1026`), `link_thin_shell_facets_serial` pairs by
  container position (`cl_FEM_DofMgr_BlockData.cpp:469-472`); a mismatch is a debug-only assert.
- [x] **B6** (Grok) Setup sequence as `TwoShellStacks::create_field`: `create_field`, `create_fields`,
  `set_solver`, `DofManager::initialize` (allocates the calculators) before the postprocessor's
  `initialize()`. F1 needs edges (HPhi `edge_h` on ThinShell) created after `MeshChecker`; F2 needs
  none. The ThinShell record takes a `nullptr` ghost sideset; nothing empty is registered on the IWG.
- [x] **B7** (Grok) The test postprocessor returns a member `Vector` of length 1 and uses a member
  coordinate `Matrix`, both sized once.
- [x] **B8** (Grok) M1 flips shell mode where it is decided (`select_recovery_mode`), before the
  node matrices are built.
- [x] **B9** (Grok) F1 pinned: R = 1, Δz = 0.1, t = 2e-3, 5 stations over 40°, 3 z rows, the TRI3
  diagonal fixed so the interior node has valence 6; tolerances relative to each field's scale.
- [x] **B10** (Grok) Winding: mid-surface triangles counterclockwise seen from outside (+r), TET4 with
  positive volume; the Kernel aborts on negative pipette volumes.
- [x] **B11** (Codex) Test comments cite the contracts they pin: `postprocessor_recovery_theory.md`
  §2.1.1 (centering), §2.1.2 (shell coordinates, BELFEM's extension), §3.2 (block inversion), and
  Zienkiewicz & Zhu 1992 for the fit.
- [x] **B12** (Grok) The suite stays single-rank (`test_fem` has no `TESTRANKS`).
- [x] **B13** *met 2026-09-23 via O1:* M3 needs the coupled T6 case to turn red (B3/B4); with an uncoupled rank-1 Gram both
  implementations keep the constant.

## 4. Open Questions

- **O1** **RESOLVED 2026-09-23 → extract + unit test** (Christian), after the fixture route (B3) ran but could not discriminate M3. Original text: the degenerate patch cannot be built with volume elements that the Kernel accepts. Extract
  the inversion (`invert_node_matrix` body) into a header-only `fn_FEM_recovery_inverse.hpp`
  (`src/fem/kernel`), called by `Postprocessor` and tested directly: healthy SPD matrix = exact
  inverse, rank-1 patch keeps a constant, coplanar patch keeps the in-plane coefficients. That is a
  source change (small, behavior-neutral) and needs its own jury pass inside this plan. Or leave
  T6 out.
- **O2** F3 (2D) doubles the fixture work; T5 could follow once F1 is green.
- **O3** Tolerances are relative to the manufactured field's scale; 1e-12 assumes κ ≤ 1e2 on F1
  (measured 16 on rutherford) — the test prints κ via the existing per-rank line.

## 5. Audit Trail

- Parent: `todo/thinshell_nodal_recovery.md` (A13, R9); juries
  `tmp/ai_exchange/review_thinshell_recovery_plan.md`.
- This plan: `tmp/ai_exchange/review_recovery_unit_test.md`.
