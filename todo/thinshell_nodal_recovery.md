# Thin-Shell Nodal Recovery: Shell-Local Patch Fit and a Centered, Rank-Aware Basis

**Date:** 2026-09-22
**Purpose:** Make the recovered nodal `H`, `B`, `J`, `J/Jc` of thin-shell layers correct on curved
shells. The patch fit of `Postprocessor` gets a node-local basis: centered and scaled Cartesian for
every postprocessor, and a curvilinear shell basis `(s1, s2, η)` for linear thin-shell layer blocks,
where `η` follows the element's own through-thickness coordinate. A degenerate patch falls back to
a lower order instead of aborting.
**Module:** `src/fem/kernel` (`cl_FEM_Postprocessor`), `src/fem/maxwell` (`cl_MaxwellPostprocessor`)
**AIs involved:** Claude (probe, plan), Codex `gpt-6-astra` + Grok `grok-4.7` (plan and code juries)
**Status:** PLAN — audited 2026-09-22 (Codex astra/high + Grok 4.7/high), amendments A1–A13 folded in (§4.1); approved with source edits 2026-09-22 (Christian; fallback = truncated pseudo-inverse, κ_max = 1e10); IN PROGRESS — gates G1–G4 PASSED on the clean A/B (2026-09-22 night, §4 R8); G5 and R9 open; R1–R6 implemented in `cl_FEM_Postprocessor.{hpp,cpp}`, syntax-gated only (both backends); code jury done (Codex astra/high + Grok 4.7/high), fixes F1–F7 applied and syntax-gated; R7 docs next, R8 gates (Christian), R9 → follow-up fixture plan (A13).

> **Scope guards (rulings, Christian, 2026-09-22):**
> - Shared layer-interface nodes stay shared (last writer wins); documented in
>   `src/fem/maxwell/doc/thinshell_postprocessor_node_sharing.md` §2.5. No `NodeDuplicates`.
> - The buffer (φ region) keeps its routing through the thin-shell conductor postprocessor; the
>   probe (§1.3) found its displayed tangential field consistent with −∇φ on the flats.
> - Conditioning: center + scale and a rank-aware fallback; no QR/SVD unless κ stays large
>   after the change (measured offline, §1.4: it does not).
> - Out of scope: the solve (it does not read the recovered fields), `HEX8TB` side-connector walls
>   (keep the centered Cartesian path), higher-order thin shells (`QUAD9TS`, `PENTA18TS`: they keep the
>   centered Cartesian path; the Maxwell order gate at `cl_MaxwellFactory.cpp:1576` runs only in the
>   cohomology branch, `:718-721`, so it is no guarantee), `HEX8TS` (no creator in
>   `ThinShellFactory`; keeps the Cartesian path), the MPI aura skip for thin shells
>   (`cl_FEM_Postprocessor.cpp:259`).

---

## 1. Current Behavior and How It Fails

Record of the diagnosis: `tmp/ai_exchange/review_thinshell_nodal_recovery.md` (jury 3/3,
Codex astra/high + Grok 4.7/high).

The recovery fits one polynomial per node in **global absolute** coordinates, `{1, x, y, z}` at
order 1 (`cl_FEM_Postprocessor.hpp:370-376`), volume-weighted over the integration points of the
selected elements touching the node (`cl_FEM_Postprocessor.cpp:817-848`, Gram at `:843`), inverts
the Gram by `posv` against the identity (`:709-712`), accumulates `b = Σ w p y` per element
(`:1046-1064`) and evaluates `p(x_n)ᵀ G⁻¹ b` at the node (`:1142-1148`). The element-side
polynomial values are cached per element because the basis does not depend on the node
(`:853-915`, `check_element_coeffs_done` `:921-943`).

| Failure | Mechanism | Evidence |
|---|---|---|
| Wrong nodal H/B/J inside curved thin layers | a linear Cartesian polynomial spans the same space in any affine frame, so it cannot follow a through-thickness profile whose direction rotates across the patch; both faces of a 1.6 µm layer get the patch mean | rutherford step 100: bend columns carry 0.5 % of the ybco jump; probe §1.4: node error 52 % of the jump |
| Jump displayed in the neighbor layers | shared interface nodes, the ybco pass writes last (`cl_FEM_Postprocessor.cpp:167`, `cl_MaxwellFactory.cpp:2528-2544`) | out of scope by ruling; disappears once the ybco faces carry the right values |
| Origin- and unit-dependent round-off | uncentered, unscaled monomials; κ(G) ≈ 4e12 on flat ybco patches | probe §1.4 |
| Run abort on a degenerate patch | `posv` with `AbortOnError = true` (`fn_posv.hpp:313-316` via `:712`) | source trace |

### 1.3 Buffer probe (2026-09-22, `tmp/ai_exchange/rutherford_figs/` + scratchpad `bufprobe.py`)

−∇φ per buffer wedge (plane fit of nodal φ, no wedge straddles a cut: max φ spread 4.4 A) against
the nodal H on the buffer nodes, tape frame, mT, means over wedges:

| region | −∇φ along / across | nodal H along / across | mean \|diff\| along / across |
|---|---|---|---|
| flats (2520) | 2.98 / 5.65 | 2.99 / 5.55 | 0.08 / 0.23 |
| bends (720) | 6.01 / 10.79 | 5.24 / 7.24 | 0.79 / 3.55 |

On the flats the displayed buffer field agrees with the φ solution; the bend gap is the §1 defect
(its bottom-face nodes carry the ybco patch mean). Limit: the buffer's displayed faces are written
by the ybco pass and by a fit dominated by 50 µm of hastelloy (0.2 : 50 volume weight), so this
does not isolate the buffer element's own `E q`. R3 does NOT make the conductor pass per layer: its
patches still span adjacent layers and stay volume-weighted, so the buffer's top face remains a
hastelloy-dominated fit and its bottom face the ybco pass's value. G1 therefore checks the buffer
per face against −∇φ, not per wedge.

### 1.4 κ and jump-retention probe (offline, scratchpad `kappa.py`, ybco patches of step 100)

Gram of the volume-weighted order-1 fit over the 11-point wedge rule
(`fn_intpoints_gauss_penta11.hpp`), node patches entirely on a flat (2628) or entirely on a bend
(604); manufactured field `f = η`, the exact node value is 0:

| basis | κ median flat / bend | node error of `f = η`, fraction of the layer jump, median (max) flat / bend |
|---|---|---|
| absolute `{1,x,y,z}` (today) | 4.1e12 / 3.1e9 | 0.000 (0.003) / 0.522 (0.555) |
| centered + scaled Cartesian | 4.4e5 / 3.3e2 | identical to absolute (affine invariance) |
| shell-local `(s1, s2, η)` | 16 / 16 (max 41) | 9e-16 / 6e-16 |

**Bottom line:** centering fixes the conditioning, not the bend; the curvilinear thickness
coordinate fixes both.

## 2. Architecture: Why a Node-Local Basis in the Existing Recovery

The fit stays the weighted least-squares patch recovery of Zienkiewicz & Zhu 1992 (Part 1, the local
polynomial fit); only the coordinates the
polynomial is written in change, and they become node-local. That costs the per-element cache of
polynomial values, because the basis now depends on the node, and nothing else. Rejected:

- **Per-face traces** (evaluate `E q` at ζ = ±1 and average per face): exact per element but no
  in-plane smoothing, and a second code path next to the fit.
- **Straight local frame `(s1, s2, d)`, d = distance along the node normal:** affinely equivalent
  to `{x, y, z}`, therefore exactly as wrong on the bend (§1.4 row 2 is the proof).
- **QR/SVD:** κ ≤ 41 after R3 and ≤ 6.5e5 for the Cartesian path; Cholesky is adequate.

## 3. Gap Table

| # | Piece | Needed for | Handled today? | Class | Citation / rationale |
|---|---|---|---|---|---|
| 1 | Node-local centered, scaled Cartesian basis | conditioning, origin independence, every postprocessor | no | (c) | `cl_FEM_Postprocessor.hpp:370-376` and the 2D/higher-order siblings |
| 2 | Curvilinear thickness coordinate η per (node, element, int point) | the bend defect | no | (c) | `η = (1+ζ)/2 · t_e` if the node is on the element's bottom face, `−(1−ζ)/2 · t_e` if on its top face; ζ = last row of `integration()->points()` |
| 3 | Bottom/top face membership of a node in a TS element | row 2 | implicit | (a) | from the element's bottom/top facet node lists, not the local index (A5); factory order today: PENTA6TS 0–2 / 3–5 (`cl_ThinShellFactory.cpp:1555-1561`), QUAD4TS 0–1 / 2–3 (`:1464-1467`); `MeshChecker::swap_quad4` would mix them (`cl_MeshChecker.hpp:207-211`) |
| 4 | Element thickness `t_e` | row 2 | no | (a) | the block's prescribed thickness, as the evaluator uses it (`cl_EF_PENTA6TS.cpp:122-123`, `cl_EF_QUAD4TS.cpp:98-101`), geometric only when it is NaN (A1) |
| 5 | In-plane `s1, s2` | row 2 | no | (a) | projection of `x_k − x_n` on the tangent plane of the node normal = mean of the patch elements' bottom→top directions |
| 6 | Scales `h` and `τ` | κ | no | (a) | `h` = largest distance from the node to a patch integration point (Euclidean, the thickness part is (t/2R)² of it), `τ` = max \|η\| |
| 7 | Rank-aware fallback | robustness | no (abort) | (c) | `posv( …, false )`, κ test on the factor; truncated eigen solve (A3) |
| 8 | Node evaluation `p(x_n) = e₀` | simplification | no | (a) | follows from centering; removes the node-coordinate poly in phase 2 (`:1133-1148`) |
| 9 | Patches spanning two layers (conductor pass, metal/metal interface nodes) | row 2 | — | (a) | η is signed and continuous through the node: negative in the layer below, positive above |
| 10 | Cohomology-cut sibling nodes in the patch (`:703-707`, `:1066-1115`) | row 1 | yes | (b) → O2 | siblings share coordinates; the basis must be built from the ORIGINAL node's frame for all siblings |

### 3.1 Cross-cutting

- **The basis now depends on the node.** `compute_node_matrices` and phase 1 of `recover_fields`
  both evaluate the polynomial per (node, element, k); the element cache
  (`compute_element_coeffs( Element*, … )`, `check_element_coeffs_done`) goes away. Phase 1 changes
  from "one `tB` per element, added to each of its nodes" to "one `tB` per (element, node)".
  Save-time cost grows by the node count per element (6 for PENTA6); no hot path.
- **Flats and non-shell blocks must not change beyond round-off.** On a flat patch η is affine in
  the Cartesian coordinates, and centering is affine; both span the same space as today. This is
  gate G2/G3.

## 4. Ordered Steps

- [x] **R1** Basis API in `Postprocessor`: `compute_local_poly( node, element, k )` filling
  `mPoly`; node frame (origin, `h`, and for shell mode normal, tangents, `τ`) computed once per node
  into member scratch (sized at `initialize`, no per-call `Vector`/`Matrix` temporaries).
  Cartesian mode: `(x_k − x_n)/h` into the existing `poly*_2d/3d`.
- [x] **R2** Shell mode selection per postprocessor: on when every selected block's element type is
  `PENTA6TS` or `QUAD4TS` (reuse the `tIsThinShell` test style of `:244-249`; `HEX8TB`, `HEX8TS`
  stay Cartesian). Mixed selections → Cartesian mode (after: R1).
- [x] **R3** Shell basis: `[1, s1/h, s2/h, η/τ]` (3D), `[1, s/h, η/τ]` (2D); η per gap rows 2–4;
  order 1 only, `BELFEM_ERROR` at setup if a shell-mode postprocessor sees `mOrder > 1` (after: R1, R2).
- [x] **R4** `compute_node_matrices`: per node, accumulate `G` from the node-local basis over the
  patch (originals + selected siblings, as today); drop the element-coefficient cache
  (after: R1–R3).
- [x] **R5** `recover_fields` phase 1: per element and per owned patch node (incl. the sibling
  crediting of `:1066-1115`), `tB_n += w · p_n(x_k) · Y_k`; phase 2: value = row 0 of `G⁻¹ b`
  (after: R4).
- [x] **R6** Rank-aware fallback per A2/A3 (truncated eigen solve, κ_max = 1e10).
- [x] **R7** Docs: `src/fem/maxwell/doc/postprocessor_recovery_theory.md` (node-local basis, shell
  mode, fallback; replace "typically well-conditioned" with the measured κ), and
  `thinshell_postprocessor_node_sharing.md` §2.5 (fixed, how) — docs land without a round, Codex
  language sweep.
- [ ] **R8** Gates (Christian runs builds and decks):
  - **G0** compile: `make belfem_kernel belfem_maxwell` (GCC `-Werror`), then `make check-fast`,
    including the new unit test of R9.
    *[◐] 2026-09-22: `build/` (Armadillo) compiled and linked with the change, test suite green (Christian); `cmake-build-debug/` (Blaze) not yet rebuilt; the R9 test does not exist yet, and no existing test runs `recover_fields` on a thin shell.*
  - ***Clean A/B, 2026-09-22 23:05–23:45 (Claude, scratch worktree at HEAD 42b5e42, own build dir, same configure as `cmake-build-debug`; A = HEAD, B = HEAD + the two Postprocessor files; plus B twice for the noise floor and B at `-np 4`).*** Solver not bitwise reproducible run to run (B vs B: I₂ 5e-10, U₂ 3e-5, air max 1.85e-6), so every gate is read against that floor.
    - [x] **G1 PASS.** A: bend jump in silver/buffer (6.43/6.00 mT), YBCO 0.03, 720/720 columns wrong. B: YBCO 12.45 mT, 720/720 bend + 2520/2520 flat columns right, silver/buffer ≤ 0.09 mT; zero truncations; κ ≤ 12 (air), 3.2 (shells).
    - [x] **G2 PASS** where the patch is a true flat: 2460 nodes clear of the flat/bend junction, H 2.6e-7, J 1.6e-7 of the peak (floor 8e-7). The 168 nodes whose patch reaches within 0.3 mm of the junction differ up to 2.5e-3 (H), 9e-5 (J): the layer there is extruded along averaged node normals and tilts up to 1.25° from y, so it is not a prism and the Cartesian fit carries the bend defect in small; B is the correct one. The "flat by x" criterion of the gate was too loose.
    - [x] **G3 PASS vs floor:** air off the shell A→B max 2.2e-6 / median 5e-9, B→B max 1.85e-6 / median 7e-9.
    - [x] **G4 PASS off the partition boundaries:** `-np 4` vs serial max 5.4e-6 (99 % 6.4e-7); all 320 nodes above 1e-3 are thin-shell nodes on partition boundaries — the thin-shell aura skip, unchanged and out of scope (not re-measured for A).
    - [ ] G5 (a QUAD4TS deck) open.
  - *G1 run 2026-09-22 22:00 (Claude, scratch copy of the deck, `cmake-build-debug` binary with the final code): **PASS** — YBCO holds 99.9 % of the across-tape jump in 720/720 bend columns (was 0.5 %) and 2520/2520 flat columns (114 were wrong before); silver/buffer jumps ≤ 0.09 mT; zero truncated patches; largest κ: air 12, shell conductors 3.2, YBCO 3.2 (was ~4e12). G2 **inconclusive**: median change 1e-7 but 34 flat YBCO nodes change by up to 0.022 Jc (isolated one-element footprints); the old fit on those patches is exact to 50 digits, so the difference is in the element input, and the baseline binary (17:38) predates the peer session's `src/mesh` edits (U₂ differs by 4e-5, `.bfm` byte-identical). A clean A/B needs a baseline from the same tree minus this change. G3: air median 1.8e-8, max 5.8e-6 (gate 1e-6 missed, same confound).*
  - **G1** rutherford (fresh run; delete `*.bfm`, `memdump.hdf5`, `iv_results.csv` first):
    `bcols.py` step 100 — ybco carries ≥ 95 % of the across-tape jump in ≥ 95 % of bend columns;
    `bprof.py` bend stations — B jump within 5 % of μ0 J d from J/Jc; buffer tangential H per face
    against −∇φ (bottom face ≙ ybco top value, top face ≙ hastelloy fit), reported, not gated;
    **zero fallbacks** counted in the shell passes (A3); side-connector walls: none in this deck.
  - **G2** flats: nodal H/B/J on ybco nodes whose patch is entirely flat agree with the baseline
    (`cmake-build-debug/rutherford/baseline_pre_recovery/rutherford.e-s.00100`, Christian 2026-09-22) to 1e-3 of the local field magnitude. The baseline carries κ ≈ 4e12 round-off, so a
    tighter tolerance tests the baseline, not the change (A6).
  - **G3** non-shell fields: air nodes more than one element away from the shell agree with the
    baseline to 1e-6 relative (air patches are not micron-thin; report the max); `I_1`, `U_1` in `iv_results.csv` identical (the solve is untouched).
  - **G4** MPI: rutherford at `mpirun -np 4` vs serial, nodal fields within 1e-8 relative on nodes
    off the partition boundaries (the thin-shell aura skip is unchanged and out of scope).
  - **G5** one deck with a 2D thin shell (QUAD4TS) against its baseline, same tolerance rule as G2.


### 4.1 Amendments from the plan jury (2026-09-22, `tmp/ai_exchange/review_thinshell_recovery_plan.md`)

- [x] **A1** (Codex P1) η uses the block's prescribed thickness `t_b` (`Block::thickness()`), the
  same value the edge evaluator blends with; geometric separation only when it is NaN, and only for
  the normal direction. η = (1+ζ)/2 · t_b below-face node, −(1−ζ)/2 · t_b top-face node.
- [x] **A2** (Codex P1, Grok) `posv` destroys the matrix: keep a copy of `G` before the call for
  the fallback (A3). Guard `h`, `τ` and `G(0,0) = Σ w` for zero/non-finite before any division: a
  zero weight sum or a non-finite coordinate is a defect → `BELFEM_ERROR` (setup, runs once); a
  zero scale (degenerate patch) is left unscaled and caught by A3.
- [x] **A3** (Codex P1, Grok; **ruled 2026-09-22, Christian; amended the same evening after the code jury**) ~~Cholesky stays the default~~ → **always** one path: the constant coefficient solved exactly, c₀ = (b₀ − G₀ᵣ cᵣ)/G₀₀, and the gradient block S = Gᵣᵣ − Gᵣ₀G₀ᵣ/G₀₀ inverted by a truncated SVD (`gesvd`; eigen decomposition of an SPD block), singular values below σ_max/κ_max dropped. A healthy patch truncates nothing and equals the exact inverse; a degenerate one loses only the undetermined direction. κ is exact per node; the run reports the count of truncated nodes and the largest κ. Cost: one small SVD per node at setup (4×4 for Maxwell; ≤ 35×35 for order-4 Cartesian recovery, ~0.1 ms each, split over the ranks). Reason: the Cholesky diagonal is only a lower bound on κ (Codex code round, probe: 7.8e9 reported for κ = 1.6e10), and the truncated pseudo-inverse of the full G loses the intercept (constant 1 recovered as 0.5). After a
  successful factorization, `κ_est = (max L_ii / min L_ii)²` from the factor's diagonal. When
  `info > 0` or `κ_est > κ_max = 1e10`, the node falls back to a **truncated eigen solve** of the
  symmetric fit matrix (`eigen_sym`, equal to its SVD): eigenvalues below `λ_max / κ_max` are
  dropped, the pseudo-inverse keeps every well-determined direction and reduces to the weighted
  mean when only the constant survives. `info < 0` is an illegal argument → `BELFEM_ERROR`. Count
  fallbacks per postprocessor, one warning line per run; G1 requires zero in the shell passes.
- [x] **A4** (Codex, Grok) The node frame (origin, `h`, normal, tangents, `τ`) is built once per node
  from the same element union the Gram uses (original + selected siblings), stored per node next to
  `mNodeMatrices`, and reused in phase 1. Face membership is looked up through the node the element
  actually references (sibling or original).
- [x] **A5** (Grok) Membership from the element's bottom/top facet node lists (the TS element's
  facet accessor), asserted consistent with ζ = ∓1 of the reference element.
- [x] **A6** (Claude, Grok) Gate tolerances as in G2/G3/G5 above; exactness is gated by R9.
- [x] **A7** (Codex) The shell and Cartesian bases get a derivation note in R7: ZZ 1992 fit,
  written in a node-local frame; the curvilinear η is BELFEM's extension, labeled as such, with its
  reproduction property (exact for fields linear in (s1, s2, η) on one layer; a kink between two
  layers of one pass is fit by one line).
- [x] **A8** (Grok) Shell mode = the postprocessor's (single) element type is `PENTA6TS` or
  `QUAD4TS`; do not reuse `tIsThinShell` (it includes `HEX8TS`, `HEX8TB`). The shell basis is order 1
  with 4 (3D) / 3 (2D) coefficients regardless of the mesh-wide `mOrder`; no order error.
- [x] **A9** (Grok) Tangents: `t1 = normalize( e_k × n )` with `e_k` the Cartesian axis least
  aligned with `n`, `t2 = n × t1` — defined for every unit normal.
- [x] **A10** (Grok) Phase 1 reuses one member `tB`-sized scratch per (element, node); no per-pair
  allocation; the per-node accumulators stay the existing map.
- [x] **A11** (Grok) Side-connector walls copy tape H/B after the ybco pass
  (`cl_MaxwellFactory.cpp:2464-2471`); they inherit the corrected values. Named in R7.
- [x] **A12** (Grok) Naming: `a`/`t` prefixes in `src/fem/kernel` (no math exemption there).
- [x] **A13** (Codex) **R9** — *done 2026-09-23 (11 tests, M1–M5; see `todo/thinshell_recovery_unit_test.md`); moved to `todo/thinshell_recovery_unit_test.md` (2026-09-23); was blocked 2026-09-22, routed by Christian to a DofManager-based thin-shell fixture:* `TS_TestStack2D` builds a `DofManagerBase` without IWG, while `Postprocessor` needs a `DofManager` (`iwg()`, block calculators, `_Volumes`); needs a DofManager fixture or the basis/inversion as free functions (plan amendment for Christian). — new unit test `tests/fem/test_PostprocessorRecovery.cpp` on the
  existing `tests/fem/support/cl_TS_TestStack.hpp` fixture where it fits, else a direct fit-kernel
  test: (a) a linear field on a translated and scaled TET patch is reproduced to 1e-12; (b) a field
  linear in η on a curved PENTA6TS stack is reproduced at both faces to 1e-12; (c) a degenerate
  patch takes the eigen fallback, reproduces the well-determined directions and never aborts; (d) a quadratic field is
  reproduced by the order-2 Cartesian basis; (e) sibling nodes use the union frame.

## 5. Open Design Questions

- **O1** Scale for η: `τ = max |η|` in the patch (one-sided patch: the layer thickness; two-sided:
  the thicker neighbor). Alternative: the node's own layer thickness. Proposal: `max |η|` — κ 16 in
  the probe.
- **O2** ~~Siblings of cut nodes: original's frame~~ **RESOLVED 2026-09-22 → A4** (frame from the
  sibling union; siblings are co-located, `add_duplicate` has no periodic caller).
- **O3** Tape edges: a node on the shell's side edge has a one-sided patch in-plane; κ stays
  bounded (probe max 41 includes them). No special case proposed.
- **O4** Should the conductor pass exclude Buffer-type blocks and let the buffer faces be written by
  their metal neighbors only? Not proposed now (ruling: keep routing); G1 decides whether it is
  needed.

## 6. Interface

No public API change. `Postprocessor` gains private scratch for the node frame and a
`mShellLocal` flag; `compute_element_coeffs( Element*, … )` and `check_element_coeffs_done` are
removed; `compute_poly( const Vector<real>& )` is kept for the Cartesian evaluation of a local
coordinate. `MaxwellPostprocessor` unchanged unless R2 needs the block element types from it.

## 7. Definition of Done

- [ ] Every gap-table row mapped to a step or an open question.
- [ ] G0–G5 passed and recorded in the devlog with the numbers.
- [ ] Docs updated (R7) and swept.
- [ ] `thinshell_postprocessor_node_sharing.md` §2.5 states the fix.

## 8. Audit Trail

- Diagnosis jury: `tmp/ai_exchange/review_thinshell_nodal_recovery.md`.
- Plan and code jury: `tmp/ai_exchange/review_thinshell_recovery_plan.md` (this plan's slug).
