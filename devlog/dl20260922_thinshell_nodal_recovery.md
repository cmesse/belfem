# Thin-shell nodal recovery: the current sheet displayed in the wrong layer on curved tapes

**Date:** 2026-09-22
**Purpose:** From a question about the Rutherford deck's terminals to a postprocessor defect: on the
bends of the five-tape `rutherford` example the Exodus nodal B put the YBCO current sheet into the
silver and buffer layers. Diagnosis jury, plan jury and code jury (Codex `gpt-6-astra`/high, Grok
`grok-4.7`/high), a node-local recovery basis with curvilinear shell coordinates, always-SVD
inversion, and the gates handed to Christian.
**Module:** `src/fem/kernel` (`cl_FEM_Postprocessor`), `src/fem/maxwell` (docs)
**AIs involved:** Claude (analysis, plan, implementation), Codex astra (three juries + two language
sweeps), Grok 4.7 (three juries; two legs lost to the turn cap and one to a script edited mid-run)
**Status:** implemented; `build/` (Armadillo) compiles, links and passes the suite (Christian);
**verified on the rutherford deck by a clean A/B overnight** (§6b): G1–G4 pass against the
run-to-run noise floor. Open: G5 (a QUAD4TS deck), R9 (a DofManager-fixture unit test).

## 1. The deck was already right

`cmake-build-debug/rutherford/input.conf` is a symlink to `examples/rutherford/input.conf`, and that
deck already carries Gregory's single-crossing convention (`input curves : [1] ; output curves : [2]`,
30 A), switched on 2026-09-20 by an earlier session together with the README, the generator and its
e2e test. Nothing to fix. What was open: Friday's three-voice round had called the field solution
"false" — that was the 5× current of the grouped-terminal deck, not the distribution. Re-run of
Friday's analysis on today's corrected-deck output (151 steps, debug tree, 17:38 build): every
crossing carries I₁ (29.5–30.0 A at 30 A), Ampère around the cable 4.998 × I₁, and the across-tape
J/Jc asymmetry is antisymmetric about the cable center on both faces and zero on the bends — the
cable-scale critical state of a 12 mm strip, ±0.31 at |x| > 3 mm, the same numbers Friday's Q1 gave.
Distribution physical (medium-high). Open: the YBCO plane integral runs 0.3–2.3 % below Ampère and
the low crossing moves with the cut plane (probably the script's interpolation, unchecked).

## 2. The B field in the tapes

Christian: "does not look physical". Decomposing nodal B into normal / along-tape / across-tape in
the unrolled (s, z) plane and chaining each wedge column through the seven layer blocks by shared
faces (`tmp/ai_exchange/rutherford_figs/bcols.py`):

| region | columns | layer holding the across-tape B jump |
|---|---|---|
| flats | 2406 / 2520 | YBCO, equal to μ0·J·d to ~1 % at every station across the tape |
| bends (r = 0.9 mm) | 719 / 720 | silver (+6.4 mT) and buffer (+6.0 mT); YBCO 0.5 % on average, 8 columns > 25 % |
| the other 114 flat columns | patches reach the bend, 20 of them lie fully on the flat | as the bends |

The total through the stack is right on the bend (13.4 vs 13.1 mT) because the per-layer jumps
telescope to outer minus inner copper face; only the placement is wrong. Along-tape B on the flats
differs by ~0.5 mT between the two YBCO faces, unexplained. J was essentially right: it is `C q` per
element, constant through the layer, so averaging through the thickness does not hurt it — but only
in the YBCO block; the neighbors' nodes carry YBCO's J (the known shared-node leak).

**Cause (jury 3/3, record `tmp/ai_exchange/review_thinshell_nodal_recovery.md`):** the recovery fits
`{1, x, y, z}` in absolute coordinates over the elements touching a node
(`cl_FEM_Postprocessor.hpp:370-376`, `.cpp` `:843`, `:1148`). A linear polynomial spans the same
space in any affine frame, so on a bend, where the sagitta of the patch is 30–55 × the 1.6 µm layer,
it cannot follow the through-thickness jump and returns the patch mean at both faces. The
superconductor pass is created after the normal-conductor pass (`cl_MaxwellFactory.cpp:2528-2544`)
and `synch_target_fields` assigns, so the mean overwrites the correct side-local values on the
shared interface nodes: the jump appears next door. Consumers: Exodus and the printed `J/Jc_max`
only; assembly, Joule heating, `dotQ` read `C q` / `E q` at integration points. Also found: the
magnesia buffer (a φ region) is post-processed through the thin-shell conductor path (Grok, routing
confirmed); the probe shows its displayed tangential field matching −∇φ on the flats to 0.2 mT.
`thinshell_postprocessor_node_sharing.md` said H and B were uncontaminated — corrected (§2.5).

**Grok tooling:** two jury legs died at the 30-turn cap (the allowlist fix Christian made does not
touch that; `GROK_MAX_TURNS=80` does) and one because `ask_grok.sh` was edited at 20:52 while the
leg that started 20:47 was executing it (`line 512: syntax error`; bash reads scripts
incrementally). All three reruns delivered.

## 3. Rulings (quiz)

Fix route: shell-local fit. Shared interface nodes: keep, documented. Buffer: probe first (done,
routing kept). Conditioning: center + scale, rank-aware fallback, "let's see what happens with κ".
Then, after the code jury: **always** the intercept-preserving truncated SVD, no Cholesky-first
(κ-diagonal bound cannot enforce the limit; cost ≤ 0.1 ms per node at order 4, split over ranks,
accepted). Baseline output to `cmake-build-debug/rutherford/baseline_pre_recovery/`. Dead
element-flag loops in `recover_fields`: removed. R9 unit test: a DofManager-based thin-shell fixture,
follow-up.

## 4. What changed (`todo/thinshell_nodal_recovery.md`, R1–R6 + A1–A12)

- Every postprocessor fits in node-local coordinates, `p( (x − x₀)/h )`: same space as before,
  origin-independent, κ 4e12 → 4e5 on the YBCO flats (offline probe `kappa.py`).
- `PENTA6TS` / `QUAD4TS` blocks fit in `(s₁, s₂, η)`: `η = ½(1+ζ)t` from a bottom-face node,
  `−½(1−ζ)t` from a top-face node, `t` the block's prescribed thickness as the edge functions use
  it, ζ the last row of the integration rule; face membership from the element's bottom facet, one
  table per block checked against every element. κ 16 (max 41); the manufactured jump `f = η` is
  recovered to 1e-15 on the bends where the Cartesian basis lost 52 %. A straight local frame
  `(s₁, s₂, d)` would NOT have worked — affinely equivalent to `{x, y, z}` (probe row 2 is the proof).
- Inversion: the constant through the Schur complement of the patch weight, the gradient block by a
  truncated SVD (`gesvd`; singular values below σ_max/1e10 dropped); a constant field is always
  reproduced; each rank prints its patch count, largest κ and the number of truncated patches.
- Per-point finiteness, collapsed-element guard, integration-rule layout as setup errors; maps own
  frames and inverses from creation and are cleared on re-initialize; the per-element polynomial
  cache is gone (the basis depends on the node), phase 1 accumulates one right-hand side per
  (element, target node).
- Docs: `postprocessor_recovery_theory.md` §2.1.1, §2.1.2, §3.2, §4.3, §12;
  `thinshell_postprocessor_node_sharing.md` §2.5. Both swept (luna/medium); the first sweep caught my
  own slip (the H jump is J·d, μ0·J·d is the B jump).

## 5. Juries

Plan round: Codex 7 findings, Grok 7 refutations + 9 risks, all verified against source, 13
amendments (thickness convention, weight sum saved before the factorization, κ bound, sibling
frames from the union, gate tolerances vs an ill-conditioned baseline, target names, side-connector
walls inherit the fix, tangent construction, citation form). My own P-2 (non-co-located siblings)
refuted by my trace. Code round: Codex 4 (intercept lost in the fallback — reproduced, constant 1
→ 0.5; κ lower bound — reproduced, 7.8e9 reported for 1.56e10; NaN through `std::max`; citation),
Grok 4 more (collapsed-element division and wrong error text, retry leak, ζ-row check debug-only,
`h` Euclidean not in-plane), all confirmed; F1–F7 applied with Christian's approval; Python replica
of the block inverse: exact to 6e-17 on a healthy patch. Record:
`tmp/ai_exchange/review_thinshell_recovery_plan.md` (+ `.diff`).

## 6. Gates (Christian runs them)

G0 `make belfem_kernel belfem_maxwell` + `make check-fast` — `build/` done and green; the Blaze tree
not rebuilt yet. G1 fresh `rutherford` (delete `*.bfm`, `memdump.hdf5`, `iv_results.csv`): YBCO
holds ≥ 95 % of the bend jump, μ0·J·d within 5 % on the bends, buffer faces vs −∇φ reported, zero
truncated patches in the shell passes. G2/G3 flats and air vs the baseline (1e-3 of the field / 1e-6
relative; the baseline carries 4e12-κ round-off). G4 `-np 4` vs serial off the partition boundaries.
G5 a QUAD4TS deck vs its baseline. No unit test runs the recovery on a thin shell (R9 follow-up).

## 6a. G1 result (22:00, Claude, scratch copy, final code)

**Pass.** YBCO carries 99.9 % of the across-tape jump in all 720 bend columns and all 2520 flat
columns; the silver/buffer jumps are gone (≤ 0.09 mT); zero truncated patches; largest κ per pass:
air 12, shell conductors 3.2, YBCO 3.2 (the old basis had ~4e12). J/Jc_max 0.7097, as before.

**G2 inconclusive, G3 missed narrowly — both confounded.** Median change on flat YBCO nodes 1e-7,
but 34 flat nodes move by up to 0.022 Jc, each an isolated one-element footprint; on those patches
the old double-precision fit equals a 50-digit reference, so the change is in the element input,
not the fit. The baseline binary (17:38) predates the peer session's `src/mesh` edits; `I₁` is
identical, `U₂` differs by 4e-5, the `.bfm` (mesh, cuts, hanging facets) is byte-identical. Air
nodes off the shell: median 1.8e-8, max 5.8e-6 against a 1e-6 gate. A clean A/B needs a baseline
built from the same tree without this change.

## 6b. Clean A/B overnight (Christian: "option 2", then "make this an overnight run")

Scratch worktree at HEAD `42b5e42`, its own build directory with `cmake-build-debug`'s configure,
`-j12`: A = HEAD, B = HEAD + the two Postprocessor files (incremental rebuild, one object),
B once more for the noise floor, B at `mpirun -np 4`. The peer session's working-tree edits are in
neither, so the change under test is the only difference. The solver is not bitwise reproducible:
B vs B gives I₂ 5e-10, U₂ 3e-5, air max 1.85e-6 — every gate is read against that.

| gate | A → B | floor (B → B) | verdict |
|---|---|---|---|
| G1 bend jump in YBCO | 0.03 → 12.45 mT; silver/buffer 6.43/6.00 → ≤ 0.09 mT; 720/720 | — | pass |
| G2 flat YBCO, patch clear of the junction (2460) | H 2.6e-7, J 1.6e-7 | H 8e-7, J 1.3e-7 | pass |
| G2 within 0.3 mm of the flat/bend junction (168) | H 2.5e-3, J 9e-5 | — | expected: the layer tilts up to 1.25° there |
| G3 air off the shell | max 2.2e-6, median 5e-9 | max 1.85e-6, median 7e-9 | pass |
| G4 `-np 4` vs serial, off partition boundaries | max 5.4e-6 (99 % 6.4e-7) | — | pass |

Every truncation count is zero; largest κ 12 (air), 3.2 (shells), per rank also at 4 ranks. The
junction nodes are the interesting row: the layer nodes there are extruded along node normals
averaged with the bend facets, so the "flat" stack is a slightly skewed prism, and the Cartesian
fit carries a small version of the bend defect that the shell basis removes — the G2 gate's
"flat by |x|" criterion was too loose, not the code. The 320 G4 outliers are all thin-shell nodes
on partition boundaries (the thin-shell aura skip, unchanged, out of scope, not re-measured for A).
The earlier G2 alarm against the 17:38 baseline (34 nodes at 0.022 Jc) is explained by that
baseline predating both this and the peer session's code and by the same solver noise; it is
superseded by this table. Scripts: `tmp/ai_exchange/rutherford_figs/{bcols,gcompare}.py`; runs and
log in the session scratchpad `ab/` (worktree left in place for Christian).

## 7. Also this session

`examples/rutherford/` is still untracked; it rides the next commit. A peer session's edits under
`src/mesh`, `src/fem/maxwell/matrices`, `tests/` (18:23–18:38) are not part of this change.
