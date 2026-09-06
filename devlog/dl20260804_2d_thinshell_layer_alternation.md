# 2-D Thin-Shell Multi-Layer Current Alternation — Root Cause

**Date:** 2026-08-04
**Purpose:** Root-cause investigation of Gregory's observation that in a 2-D
thin-shell tape with N>1 layers the transport current direction alternates
layer-to-layer (J/Jc ≈ −1 / +1 / −1 / +1), while N=1 looks physically correct.
**Module:** fem/interpolation (nedelec), fem/kernel (ThinShellFactory), fem/maxwell
**Session type:** read-only investigation (no source edits). Three blind
parallel analyses: Claude (primary), Codex, Grok — all converged.

## Problem

Model `cmake-build-debug/greg2/`: single HTS tape (2-D cross-section, h-φ
thin-shell, sidesets 5,6), sine transport current 160 A (=Ic) at 50 Hz,
power-law ybco (jc 4e10, n 25). With `layers : tape { ybco : 0.25 mum ; } ×4`
the JJCz output alternates sign per layer through the thickness (visualized in
`greg2/image.png`, range ≈ −0.87…+0.83). With one 1 µm layer the penetration
is symmetric and field-continuous with the air. 3-D thin shells do not show
this. Gregory suspected a sign error in the 2-D thin-shell curl operator.

## Root cause (confidence: high — 3/3 independent analyses agree)

`EF_QUAD4TS` breaks tangential-field continuity across shared inter-layer
edges:

1. **False sign convention in the basis.** `EF_QUAD4TS::E()` negates the
   top-edge dof: `mE(:,1) = -mS[1] * tangent * f_top`
   (`cl_EF_QUAD4TS.cpp:186-187`). The comment at :185 claims
   "`mS[1] = -mS[0]` set in `link()`", but `link()` (:73) only copies
   `edge_directions()` from the element. The ThinShellFactory builds both
   carried edges of a QUAD4TS in the *same* tangent direction
   (`cl_Element_QUAD4TS.hpp:122-128`, edge 1 = {n3,n2} parallel to edge 0 =
   {n0,n1}; clones preserve node order, `cl_ThinShellFactory.cpp:1832`), so
   in practice `mS[0] == mS[1]` and the minus is a real, uncompensated flip.
2. **Adjacent layers share the interface edge.** For identical-material
   layers (`hasDuplicates == false`), block ℓ's top edges and block ℓ+1's
   bottom edges are the *same* Edge objects
   (`cl_ThinShellFactory.cpp:1909-1926`). One shared dof is therefore
   interpreted with `+` by the layer above and `−` by the layer below: the
   through-thickness interpolation is anti-continuous. The represented
   physical h_t — and with it J = curl h — flips sign at every interface.
   This is exactly the observed alternation.
3. **Curl parameters patched to match the broken basis.**
   `EF_QUAD4TS::precompute()` sets both curl parameters to `+0.5`
   ("For symmetric current", `cl_EF_QUAD4TS.cpp:163-165`) instead of the true
   thickness derivatives (−0.5, +0.5). The per-layer current operator becomes
   ∝ (h_bot + h_top) instead of ∝ (h_top − h_bot)/t. Alves et al. 2022
   (paper8) expects the difference form (negative off-diagonal stiffness).
   Consequence: layer currents no longer telescope to the tangential-field
   jump across the whole shell (Ampère); interior dofs count twice with the
   same sign.
4. **The defect is in the solve, not only the output.** The same E/C
   operators feed assembly (`mt_maxwell_h.cpp:127` — M += EᵀE, K += CᵀC;
   `compute_j = C*q` for the power law) and the JJCz postprocessing
   (`cl_MaxwellPostprocessor.cpp:604-613`, `compute_superconductor_ts`).
5. **Why N=1 looks fine.** With one layer no interior edge is shared; the
   flip is an internal dof redefinition with no physical consequence.
6. **3-D differential.** `EF_PENTA6TS::E()` (:285-308) and `EF_HEX8TS` apply
   plain `+mS[k]` to both bottom and top face dofs, with the sign pattern in
   the reference basis. The `-mS` hack exists only in the 2-D file. Commit
   650babb6 already recorded "issues with the curl operator" for 2-D shells.

## Additional latent defects found (same file)

- `cl_EF_QUAD4TS.cpp:127` — the curl coefficient uses `mInvJ(1,0)` (= ∇ξ_y)
  where (∇η×∇ξ)_z requires ∇η_x (= `mInvJ(0,1)`). The two coincide only for
  axis-aligned tapes (Ty = 0 — true for greg2, so latent here, wrong for
  inclined tapes). Confidence high.
- Scaling: with curl parameters ±0.5 *and* `mCoeffs *= 0.5` the η-derivative
  half-factor is applied twice; |C| comes out half the consistent curl of E.
  Possibly why |J/Jc|max ≈ 0.87 rather than ≈ 1 at I = Ic. Confidence
  medium — settle with a patch test during the fix.
- Normal convention mismatch (Codex): `ThinShellFactory::process_nodes_line2`
  uses (dy,−dx) (`cl_ThinShellFactory.cpp:893`) while `EF_QUAD4TS::link()`
  uses (−dy,dx) (:57-60). Flips the global J sign only; cannot alternate.

## Fix (applied 2026-08-04, approved by Christian)

Additional finding during implementation: `Calculator::dV_ts` returns the
edge function's `mDetJ` as the per-point volume increment, and the QUAD
Gauss rules sum to 4 (e.g. `fn_intpoints_gauss_quad12.hpp`) while the PENTA
rules sum to 1 (`fn_intpoints_gauss_penta8.hpp`: 6×1/9 + 2×1/6). The 3-D
convention (validated) is therefore `mDetJ` = physical volume with
unit-sum weights; the old 2-D `mDetJ = thickness·length` over-integrated
the shell block by 4×. Notably the old half-magnitude curl made
K += CᵀC·wdV come out right (¼ × 4), while M += EᵀE·wdV stayed 4× too
large and `compute_j = C·q` under-reported |J| by 2× — which feeds the
power-law ρ(|J|) at n = 25. These compensating errors were removed as one
consistent set.

Changes in `cl_EF_QUAD4TS.{cpp,hpp}`:
1. `E()`: top dof now `+mS[1] * tangent * f_top` (3-D convention; restores
   tangential continuity across shared inter-layer edges).
2. `precompute()`: curl parameters are the true thickness-weight
   derivatives (−0.5, +0.5).
3. `link()`: curl coefficient computed directly as
   `(∇η×∇ξ)_z = mNablaEta[0]·mNablaXi[1] − mNablaEta[1]·mNablaXi[0]`, with
   ∇η derived from the actual bottom→top center vector of the element
   (nodes 2,3 minus 0,1) instead of the rotated facet normal — fixes both
   the transposed-index defect and the stacking-side sign assumption
   (`process_nodes_line2` stacks along (dy,−dx), the EF assumed (−dy,dx)).
   The `mJ2`/`inv2` construction became unnecessary and was removed.
4. `link()`: `mDetJ = 0.25·thickness·length` (physical area / reference
   area, matching the quad weight sum); `mSumW = 4.0` for the record.
5. Consistency hand-check: `C·q` now returns exactly
   J = (c_bottom − c_top)/t for circulations c set on the two edges, and
   the per-layer currents telescope to the outer-trace difference (Ampère).
   Consumers audited: `dV_ts`/`dV_hex` are the only `det_J()` readers;
   `nedelec_data_linear_h` carries no signs; `h_ghost` and the postproc
   build purely on E/C, so they follow the fix; Pipette measures QUAD4TS
   volume geometrically.

Syntax-checked against the real build flags (`/usr/bin/c++ -fsyntax-only`
with the interpolation target's flags.make). Build and smoke runs are with
Christian.

## h_ghost interaction analysis (requested by Christian, 2026-08-04)

Question: does the EF_QUAD4TS fix break `h_ghost` for multilayered 2-D or
3-D shells? Answer: no — 3-D is bit-identical, and in 2-D the fix makes
`h_ghost` consistent with its own design for the first time. Trace:

**Where ghosts exist.** Ghost facets are inter-layer facets created only at
layer boundaries whose materials differ (`hasDuplicates`,
`cl_ThinShellFactory.cpp:2139`); master = lower block's top face
(`top_facet_index`, QUAD→2), slave = upper block's bottom face (QUAD→0,
orientation 1 = forward). At such boundaries the interface edges are
duplicated (`link_elements_with_edges` uses `EdgeDuplicates` for the upper
block's bottom), so `h_ghost` is the ONLY inter-layer coupling there. For
identical-material stacks (greg2: 4× ybco) no ghost facets exist at all —
edges are shared, continuity is strong, `h_ghost` never runs.

**Linking safety.** For sideset groups the group-level `mEdgeFunction` is
never created (BLOCK-only, `cl_FEM_Calculator.cpp:541-547`); the ghost
calculator links `mEdgeFunctionMaster/Slave` with the master/slave BLOCK
element wrappers (`cl_FEM_Calculator.cpp:1358/1369`), which are QUAD4TS
mesh elements (4 nodes — the new `node(2)/node(3)` reads are in range) and
carry the tape facet pointer on every layer (`cl_FEM_DofMgr_BlockData.cpp:484`
serial, `:571` parallel), so the tangent lookup in `link()` works unchanged.

**Sign consistency.** The master EF is precomputed at facet-2 points
(η=+1 exactly → bottom columns identically zero), the slave at facet-0
points (η=−1 → top columns zero). `h_ghost` builds
`Dm·q_m = Em_top·(q_top−q_bot)/hm` and `Ds·q_s = Es_bot·(q_bot−q_top)/hs`
(`mt_maxwell_h.cpp:369-385`) — the one-sided OUTWARD normal-derivative
fluxes. With the fixed basis these equal the represented fields' actual
through-thickness derivatives, and the flux pair (−J_m, +J_s)·T̂ matches
the natural boundary terms ρ(curl h)×n with opposite outward normals —
i.e. the Ern & Guermond structure that the 3-D verification covered
(dl20260319_ghost_thinshell.md) now holds in 2-D with the same plain
`+mS` interface-column convention as the validated `EF_PENTA6TS`. Under
the OLD basis the α-jump term (pure Galerkin in E, basis-agnostic) and the
hand-built Dm flux term DISAGREED: the represented ∂η-field was
∝ −s(q_bot+q_top) while Dm computed ∝ ±s(q_top−q_bot). So 2-D `h_ghost`
was internally inconsistent before the fix and is consistent after it; it
had simply never been exercised (2-D dissimilar-material multilayer was
blocked by the alternation defect anyway).

**dV/dS.** `h_ghost` integrates with the facet `dS` (LINE2 length × 0.5,
weights sum 2), not with the edge-function `det_J` — the mDetJ = area/4
change does not touch it.

**Latent note (pre-existing, zero effect today).** The master facet-2
point sequence is index-reversed (`intpoints_quad` master overload, case 2)
while the ghost slave uses facet 0 forward; whether point k co-locates on
both sides depends on the 1-D Gauss sequence ordering. At linear order this
is provably irrelevant: EF_QUAD4TS `E` has no in-plane variation (∇ξ and
the η=±1 thickness weights are constant per facet), so Em(k)/Es(k) are
identical for all k and only Σw·dS enters. It MUST be revisited before
QUAD9TS ghosts (the kernel currently asserts linear order).

**3-D.** Untouched: only `cl_EF_QUAD4TS.{cpp,hpp}` changed; `h_ghost`,
`EF_PENTA6TS`, `EF_HEX8TS`, meshtools facet indices, and the integration
tables are all as before. 3-D multilayer behavior is bit-identical.

**Recommended ghost smoke case (later, not greg2):** a 2-layer 2-D stack
with dissimilar materials (e.g. ybco + silver) to exercise the 2-D ghost
path for the first time, cross-checked against the equivalent extruded 3-D
stack.

## Post-fix smoke run: second defect found and fixed (same day)

Gregory's re-run of greg2 (4 layers) showed a NEW pathology: Jz odd
(antisymmetric) through the thickness, +Jz_max at one face to −Jz_max at
the other. Quantitative probe of `hphi_results.e-s.00051` (t = 4.9 ms,
I = 159.9 A peak, scipy/netcdf):

- Layer currents −7.8 / −3.2 / +3.3 / +7.9 A — net 0.24 A instead of 160 A.
- The AIR-side solution is CORRECT: the tape-line nodes are properly
  duplicated below/above (201 positions, only the right tip shared), the
  φ jump across the tape ramps 160 → 0 from the left tip (the
  auto-computed cut terminates there; its terminal pair carries the full
  160 A jump — cuts terminate on side-curves, correct) to the shared
  right tip, and Hx is antisymmetric (±3×10⁴ A/m at ±1.8 mm) — a
  Norris-like sheet carrying the full 160 A as seen from the air.
- Net shell current ≈ 0 with an odd J profile algebraically requires
  c_top = c_bottom: the top-surface tie delivers MINUS the circulation.

Root cause (the "D2 top-node pairing cross-wire" already suspected in
todo/2d_thinshell_todo.md): in `hang_thinshell_edges_on_nodes_top`
(cl_MaxwellFactory.cpp:1522) the air nodes come via `to_master_orientation`
(LINE2 = swap, fn_to_master_orientation.cpp:32-37) → position order (0,1),
but `mesh::get_top_nodes` (QUAD branch) returned the facet-2 traversal
{node2, node3} = positions (1,0) — node 3 sits above node 0. The
positional `set_index` pairing therefore swapped each top edge's two φ
sources, and since the T-matrix applies coefficients (+1,−1) in edge-node
order (cl_FEM_DofMgr_DofData.cpp:3632-3641), the top tie became
φ(n1)−φ(n0) = −circulation. With the antisymmetric self-field
(Δφ_above = −Δφ_below) this forces c_top = c_bottom → zero net shell
current; the 160 A rides as a free current sheet in the tie mismatch and
the shell merely screens the ambient parallel field → the observed odd
profile.

Historical note: under the OLD EF_QUAD4TS basis (−mS[1] top dof) this tie
swap and the basis flip CANCELED at the outer surface — that is why the
1-layer case looked right before the basis fix, and why the basis bug
manifested only between layers (alternation). Two independent sign errors,
one masking the other at the boundary.

**Fix:** `meshtools.cpp get_top_nodes()`, QUAD branch — swap the two
corner nodes after `get_nodes_of_facet(2)` so the returned order is
position-aligned with the master orientation, as the hang functions
require. This repairs both `hang_thinshell_edges_on_nodes_top` (air
above) and `hang_thinshell_edges_on_edges_top` (conductor above — its
±1 weight comes from node-index matching, cl_MaxwellFactory.cpp:1729-1738).
`get_top_nodes` has no other consumers; PENTA/HEX branches (3D) untouched.

Codex audit (2026-08-04 16:12): confirmed on all claims, verdict safe for
the smoke rerun. One latent pre-existing item it surfaced: the 2-D
conductor-above path `hang_thinshell_edges_on_edges_top` cannot run today
anyway — the `to_master_orientation(Facet*, Cell<Edge*>&)` edge overload
only handles list sizes 3 and 4, not the single-edge LINE facet
(fn_to_master_orientation.cpp:594). If a 2-D shell ever gets a conductor
block on its top side, that overload needs a size-1 case first. QUAD9TS
corner-swap alignment checked consistent (LINE3 keeps the midnode slot),
but quadratic 2-D shells remain unimplemented upstream regardless.

Expected after re-run: net shell current = I(t); layer currents all the
same sign, roughly uniform through thickness; edge-penetration profile
along the width unchanged (the air solution was already right).

**Validated (Gregory, same day):** re-run looks correct — uniform
same-sign current through the thickness, and the global sign is physical
(field runs clockwise / B in −z sense, current in −z, consistent with
Ampère). Remaining checks: N=1 vs N=4 at equal total thickness, Norris
comparison.

## Validation still open

1. N=1 re-run (field continuity with air at BOTH surfaces — the old top
   trace was sign-flipped relative to the circulation convention, which a
   symmetric single tape can visually mask).
2. N=4 re-run: alternation must vanish; N=1 vs N=4 at equal total
   thickness should agree.
3. |J/Jc|max at I = Ic should now approach ~1 (the 2× under-reported J is
   gone), then compare against Norris.
4. Expectation setting: the 1-layer results will change QUANTITATIVELY,
   not just cosmetically — the old kernels evaluated the power-law ρ at
   |C·q| = |J|/2, i.e. at (1/2)^(n-1) ≈ 6e-8 of the true resistivity for
   n = 25, so the old runs had essentially no resistive dissipation at
   operating current. "Different from before" is expected and correct;
   the reference is the analytical solution, not the old output.

## Attribution

Phenomenon reported by Gregory Giard. Blind parallel audits by Codex and Grok
(exchange thread `tmp/ai_exchange/greg2_layer_alternation_*.md`, ephemeral);
Codex contributed the mS[0]==mS[1] topology proof and the literature
cross-check, Grok the refutation of stale todo claims (ghost-dof sizing,
default edge directions) and the N=2 discriminator.
