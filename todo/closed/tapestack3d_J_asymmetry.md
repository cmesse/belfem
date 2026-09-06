# tapestack3d: |J| Asymmetry in the Outer Copper Layers, Overnight Investigation

**Date:** 2026-08-15 (overnight, read-only)
**Purpose:** Localize the source of the jagged, z-periodic |J| streaks seen in one outer
copper layer of the tapestack3d example, while the other outer copper layer shows smooth
z-invariant bands.
**Status:** investigation complete; verdict below. No source, input, or run state was
modified. The analysis used the exodus files already on disk (read-only), the run log, and
static source reading; Codex and Grok audited independently from opposite ends.

---

## 1. Setup (as read from the inputs and the run)

- **Deck:** `cmake-build-debug/tapestack3d/input.conf` + `tapestack3d.geo`.
  8 tapes, 4 mm wide, 0.1 mm pitch, 10 mm long (periodic in z), soldered stack
  (blocks 2:8 = Pb38Sn62), cylindrical air domain r = 50 mm. Sigmoid transport current,
  amplitude 2 kA, period 10 s; at the analyzed save (t = 5.95 s) I ≈ 1.70 kA. [The
  deck comment's Ic = 512 A is stale: the sp-ap jc database gives ≈ 3·10¹⁰ at the
  operating field, so true stack Ic ≈ 1.5–1.8 kA and the run sits at ≈ 1.0–1.2 · Ic —
  the current-sharing knee, not deep saturation; see §3b-iii.]
- **Layer stack per tape** (7 layers, 8 node planes 0..7, stacked along +y after the
  deck's `-5,-6` flip makes all eight tapes uniform):
  copper 20 µm (L0, planes 0-1) / silver 2 µm (L1) / YBCO 1.6 µm (L2) /
  magnesia buffer 0.2 µm (L3) / hastelloy 50 µm (L4) / silver 1.8 µm (L5) /
  copper 20 µm (L6, planes 6-7). Edge coating (side-connector walls): **on**.
- **Elements:** TET4 air/solder volumes, PENTA6TS thin-shell layers, HEX8TB coating
  walls. Linear order throughout; `cl_Maxwell_TMatrix` (TRI6/TET10 path) is **inactive**
  here. The linear hanging-edge condensation is the live constraint path.
- **Run:** `hphiTrun` (coupled h-φ/thermal), 8 MPI ranks, Blaze, STRUMPACK (magnetic),
  git `5aae789b` + uncommitted changes (which contain the 2026-08-14 `EF_TET4`/`EF_TET10`
  fixes, later committed as `42eb4fc0`). Tolerance 1e-11, Newton, BDF5, Anderson on.
- **Exodus mapping:** the layer blocks are `tape_70_copper` (L0) … `tape_76_copper` (L6),
  each spanning **all eight tapes**; coating blocks `coating_78..108_copper`; J is a
  nodal field (SPR-recovered), element fields are off in this run.
- **Which PNG is which:** by content, the **noisy** screenshot (`yplus.png`) shows the
  outermost −y face of the stack = **tape 1, layer L0** (`tape_70_copper`, node plane 0);
  the **smooth** one (`yminus.png`) shows the outermost +y face = **tape 8, layer L6**
  (`tape_76_copper`, plane 7). (The filename↔camera mapping appears swapped relative to
  the axis triads, but the content mapping above was confirmed by re-rendering the exodus
  data: the plane-0 face reproduces the jagged fingers, the plane-7 face the smooth bands.)

## 2. The symptom, restated after measurement

Re-rendering |J| per layer/plane from `hphi_results.e-s.00119` (t = 5.95 s) and computing
a z-roughness metric (mean |second difference| of |J| along z within 0.05 mm x-bins)
changes the framing of the symptom substantially:

- The noise is a **layer-index effect, not a tape-side or air-facing effect**.
  L0 is noisy on *every* tape it belongs to (tape 8's L0 faces solder and is exactly as
  noisy as tape 1's air-facing L0); L6's outer plane is smooth on every tape.
- Relative roughness matrix at t = 5.95 s (percent of local mean |J|):

  | tape | plane 0 (L0 outer) | plane 1 (L0/Ag) | plane 6 (Ag/L6) | plane 7 (L6 outer) |
  |---|---|---|---|---|
  | 1 | **3.8** | **4.2** | 2.3 | 0.5 |
  | 2 | 0.7 | 0.7 | 0.5 | 0.5 |
  | 3 | 0.8 | 3.3 | 2.6 | 0.6 |
  | 4 | 0.9 | 2.6 | 2.3 | 0.6 |
  | 5 | 1.2 | 3.7 | 2.3 | 0.6 |
  | 6 | 1.8 | **4.7** | 2.5 | 0.7 |
  | 7 | 1.6 | 2.2 | 2.3 | 0.6 |
  | 8 | **3.8** | **4.6** | 3.0 | 0.5 |

- The noisy planes are the **YBCO-side copper/silver system** (planes 0, 1) and, more
  weakly, the Ag/Cu interface plane 6; **plane 7 is clean everywhere**. The outermost
  tapes (1, 8) are worst; **tape 2 is anomalously clean everywhere** (open question).
- The pattern is a **standing, slowly growing mode**, not per-step churn: the detrended plane-0
  pattern correlates at 0.996–0.998 between consecutive saves (50 ms of simulated time,
  dozens of solver steps apart), with amplitude growing roughly in proportion to the mean
  |J| (relative roughness ≈ constant-to-slowly-growing through saturation).
- Time evolution: negligible before t ≈ 3 s, comparable on both layers at t = 5.0 s
  (I ≈ 1.0 kA ≈ 2 · Ic; note the deck comment's current table assumes the older 1024 A
  amplitude — with `amplitude : 2 kA` the Ic crossing is already at t ≈ 4.4 s), then the
  L6 side relaxes to 0.5 % while the L0 side stays/grows at ~4 %.

## 3. Verdict and classification

**(e) "other", with the premise itself refuted.** No defect was found in the
formulation (a), the element interpolation (b), the T-matrix/condensation (c), or the
postprocessor's layer/side handling (d) that could produce this picture. Two real code
defects were found incidentally (§5), but both are demonstrably inert for this run.

The finding breaks into three parts:

**3a. The symmetry premise is wrong (Grok's key contribution, code-confirmed).**
Magnesia carries no `rho` (`cl_Material_Magnesia.cpp`, `MaterialType::NonMetal`), so
`ThinShellFactory::create_buffers` (`src/fem/kernel/cl_ThinShellFactory.cpp:1991-2049`)
reclassifies the buffer layer to `DomainType::Buffer` and hangs its edges on their own
nodes (a φ region), and `create_ghost_facets` (`:2054`, Buffer test at the interface
loop) **skips** any interface involving Buffer. Each tape is therefore two electrically
separate conducting halves: (Cu L0 | Ag L1 | YBCO) below and (hastelloy | Ag L5 |
Cu L6) above, joined only through the edge-coating wall at the rims and, across
tapes, through the solder gaps. The stack is additionally asymmetric in materials (2 vs 1.8 µm
silver; 50 µm hastelloy on one side only). "The two outer coppers should look alike" was
never a valid expectation; only their *mean* current density agrees (both ≈ 1.2·10⁵ A/m²,
consistent with E/ρ at the measured terminal voltage, sanity-checked against
`iv_results.csv`).

**3b. The noise is in the solution, localized where the physics is stiff — and it is a
property of the CONVERGED solution, not of under-convergence.** [CORRECTED 2026-08-15
morning: the overnight draft claimed the 1e-11 tolerance was never reached; that was a
misreading of the log's dB convention. BELFEM prints `10·log10(eps)`
(`cl_FEM_Controller.cpp:1469`, also the `:1157` comment), so −110 dB IS 1e-11, and the
deck's own "−156 dB observed ↔ 2e-16" comment confirms it. Measured over the whole log:
568 of 569 step-final magnetic residuals are ≤ 1e-11 (−110…−156.5 dB; single outlier
4.9·10⁻¹¹); the −106…−108 dB values are mid-step, which is why those iterations
continue. The anti-checkerboarding criterion (deck comment; Messe et al. 2023 §4,
eq. 11) IS enforced.]

The corrected picture: the noisy entities are exactly the conductors that current-share
with the saturating YBCO through the 2 µm silver (planes 0/1, worst at the
high-self-field outer tapes); the clean face is separated from the YBCO by the
insulating buffer + 50 µm hastelloy. The standing z-mesh-scale mode survives *inside*
the 1e-11-converged solution. Newton also shows a stall signature (bit-identical
residual over 4–5 iterations while the relaxation changes, e.g. steps 300, 496, 498, 504
in `out.txt`) before Picard closes the step — at ~2·10⁻¹¹ that is grinding near the
attainable floor, an efficiency wrinkle rather than a correctness defect.

**3b-ii. Cause closed by the morning probe (2026-08-15): the copper is an E-field
microscope of the YBCO's discrete flux-front raggedness.** Direct measurement on
`s.00119`, all eight tapes, copper plane 0 vs the same tape's YBCO nodal J/Jc field:

  | tape | YBCO J/Jc min..max (spread) | YBCO z-ripple % | Cu pl0 ripple % | corr | slope |
  |---|---|---|---|---|---|
  | 1 | 1.23..1.65 (0.41) | 0.59 | 3.81 | +0.67 | 5.1 |
  | 2 | 1.09..1.12 (0.03) | 0.06 | 0.71 | +0.68 | 7.9 |
  | 3 | 1.10..1.12 (0.02) | 0.06 | 0.78 | +0.52 | 5.3 |
  | 4 | 1.07..1.12 (0.05) | 0.10 | 0.85 | −0.05 | −0.2 |
  | 5 | 0.93..1.13 (0.20) | 0.32 | 1.19 | +0.08 | 0.1 |
  | 6 | 0.92..1.13 (0.21) | 0.26 | 1.76 | +0.19 | 0.5 |
  | 7 | 0.93..1.10 (0.17) | 0.25 | 1.62 | +0.48 | 2.0 |
  | 8 | 1.07..1.51 (0.44) | 0.71 | 3.77 | +0.28 | 1.4 |

The mechanism, each link now evidenced:

1. The YBCO sheets are overcritical (J/Jc up to 1.65) with a mesh-scale z-raggedness
   of 0.06–0.71 % — the discrete footprint of the sharp constitutive law's front on
   the 0.05–0.25 mm unstructured surface mesh. Uniformly loaded tapes (2–4, spread
   ≤ 0.05) have essentially no raggedness; tapes with strong in-plane J/Jc structure
   (1, 8; edge concentration from the stack self-field) have the most.
2. The resistivity branch converts J/Jc ripple into E ripple amplified by its local
   logarithmic slope. From the terminal V–I (`iv_results.csv`: E ≈ 3.4·10⁻⁴ V/m at
   J/Jc ≈ 1.4), the piecewise law's effective exponent here is n_eff ≈ 3.6–7.
3. The severed YBCO-half stabilizer reads E directly (J = E/ρ): predicted copper
   ripple ≈ n_eff × YBCO ripple. Measured point-wise slope: 5.1 on tape 1 — matching
   n_eff — with spatial correlation up to 0.70 in the edge bands; and the aggregate
   amplitude tracking holds across all eight tapes, resolving the tape-2 anomaly
   (uniform J/Jc → nothing to imprint → clean copper).
4. The buffered half never sees this E structure: plane 7 clean everywhere.

The point-wise correlation weakens toward the upper tapes (0.28 on tape 8) while the
amplitude tracking stays intact; the likely benign reason is that the *written* J/Jc
uses the postprocessor's jc(|B|, β) at recovered fields, which decouples point-wise
from the assembly's E where the field-angle structure is strong. Aggregate evidence is
the load-bearing part. Standing-and-growing follows naturally: the front pattern is
pinned to the mesh and its amplitude scales with the transport current.

**3b-iii. Round-3 verification (Codex + Grok adversarial + literature, 2026-08-15
afternoon) — mechanism CONFIRMED in corrected form.** Both auditors' challenges were
run to ground with data; the final reconciled account, every number measured:

1. **The run is at ≈ its true critical current, not 3.3 · Ic.** The `sp-ap.hdf5`
   database (`n`/`jc` tables over T × log10 B × angle, log-stored) gives
   jc(77 K, 100–200 mT) ≈ 2.8–3.6·10¹⁰ A/m² — the deck comment's jc = 10¹⁰ (and its
   Ic = 512 A table) is stale. True stack Ic ≈ 1.5–1.8 kA; at t = 5.95 s,
   I = 1.70 kA ≈ 1.0–1.2 · Ic_true.
2. **Local law slope n_eff ≈ 17–22** (database n at the actual |B| ≈ 100–200 mT;
   Codex's correction — my earlier 5.7 secant was a proxy artifact, see 4).
3. **True overcriticality ≈ 1.05–1.15**: with n ≈ 19, (1.10)¹⁹ ≈ 6.7 = the measured
   E/ec from the terminal V–I (E = 6.8·10⁻⁴ V/m at 10 mm periodic length). Closed.
4. **The written J/Jc field is a biased proxy** (mean ≈ 1.4, i.e. ≈ +30 %): the
   postprocessor evaluates jc at recovered B/β nodal fields, and assembly uses
   |n·b| where the postprocessor uses a signed acos (Codex). This bias and its
   texture explain the diluted regression slope (5.1) and the imperfect pointwise
   correlations Grok attacked (W2 resolved).
5. **The discriminating signature triple** (Grok's decisive checks, executed):
   copper noise correlates +0.59 with J/Jc, **−0.19 with the YBCO's own |J|**
   (current sharing — no artifact produces a negative J↔J correlation), and
   **−0.04 with |B|** (common-driver alternative dead). Copper J is computed
   jc-free (C·q), so the JJC correlation cannot be a shared-postproc artifact.
   The clean face's weak +0.44 correlation is rim-weighted (rim +0.39 vs interior
   +0.22) at 7× smaller amplitude — the rim-fed remnant.
6. **Chain, final:** mesh-scale texture (~0.2 % true) in the local overcriticality —
   seeded by the discrete sampling of the anisotropic jc(B, β) and the current
   front on the 0.05–0.25 mm surface mesh — is amplified by n_eff ≈ 19 into a ~4 %
   E texture, which the metallically bonded stabilizer displays directly
   (0.2 % × 19 ≈ 3.8 % = the measured fingers). The severed hastelloy-half face
   shows the 7×-attenuated rim remnant (0.5 %).

Literature placement (Grok's sweep + my pass, curated citations): this is NOT
Messe et al. 2023 §4 checkerboarding (`messe2023.txt:555-567`; that is the
0↔±2jc iteration artifact — ours survives at 568/569 steps ≤ 1e-11); the closest
published analog is Messe et al. 2023 §5's coarse-mesh "staircase" with integrated
losses still agreeing (`messe2023.txt:651-686`), and the supported mitigation is
a-priori graded tip refinement (`messe2023.txt:703-714`). Dular et al. 2021's
oscillations are mixed-space/inf-sup phenomena (different lever); its resistivity
regularization remark (`dular2021.txt:1109-1111`) is context, not a fix here. The
piecewise law is already C1-smooth in log-log (verified, `powerlaws.hpp:335-453`);
do NOT soften it for this — with n_eff ≈ 19 a material change would change physics.

Consequences: the picture is a magnified display of sub-resolution structure at the
current-sharing point; layer-mean currents and integrated quantities are unaffected
(Messe "staircase" precedent). The mesh lever (tapeResolutionTip) attacks the seed
texture; the ×19 amplification is a material property and stays. Ghost-penalty
calibration (A3) demoted to secondary. Report-only by-catch: (a) the assembly/
postproc angle-convention difference (jc(θ) table symmetric to ~1 % → benign for
jc; n(θ) asymmetric ≤ 9 % → minor; unify eventually); (b) plotted JJC overstates
overcriticality ≈ +30 % — read JJCz pictures accordingly; (c) the deck's jc/Ic
comment block is stale vs the sp-ap database and should be rewritten.

**3c. The picture-layer caveats are real but do not touch the outer faces.** The known
interface-node SPR leak (shared `Node*` between layer blocks,
`thinshell_postprocessor_node_sharing.md`) makes nodal J on planes 2–5 unusable: the
YBCO's ~2.5·10¹⁰ A/m² is visibly painted onto the silver-L1/buffer faces (the "buffer"
block renders at ~1.3·10¹⁰ A/m² although magnesia is an insulator). Interface planes 1
and 6 mix two materials' J in one nodal value by construction. Neither mechanism can
write the outer planes 0 and 7, which are single-block nodes.

## 4. Evidence chain, and what was ruled out

Everything below is static source trace plus read-only data forensics on the exodus
files (evidence levels 4–5 of the protocol ladder; nothing here is "verified" by an
executable gate).

**Ruled out, with the checks that ruled them out:**

- **(b) `EF_PENTA6TS` basis/curl.** Hand-derived ∇×(Whitney × (1∓τ)/2) for all six dofs
  and matched them column-by-column against `mCoeffs`/`mCurlPars`
  (`cl_EF_PENTA6TS.cpp:144-259`): `E` and `C` are mutually consistent and analytically
  correct; the top-dof signs follow the *fixed* QUAD4TS convention (`+mS[k]`, no
  hard-coded negation). Codex independently reached the same conclusion. Note the new
  edge-function test battery (`tests/fem/test_EdgeFunctions.cpp`) does **not** yet cover
  the TS elements — worth adding, but no defect found here.
- **The 2D bug pattern (Aug 2026, EF_QUAD4TS + `get_top_nodes`) does not recur in 3D.**
  `get_top_nodes` PENTA branch returns facet-4 nodes {3,4,5}, positionally aligned with
  the bottom {0,1,2} (`cl_Element_PENTA6TS.hpp:100-107`) — no corner swap needed, unlike
  the QUAD branch. `to_master_orientation` TRI3 (`fn_to_master_orientation.cpp:45-81`)
  was checked against `Facet::compute_orientation` (`cl_Facet.cpp:104-130`) for all three
  transposition cases: correct. The edge-on-edge ties
  (`MaxwellFactory::hang_thinshell_edges_on_edges_bottom/top`,
  `cl_MaxwellFactory.cpp:1767-1949`) are self-checking at runtime (node-identity match or
  hard abort), and the run did not abort. The edge-on-node weights
  (`DofData::create_dofwise_t_matrices_master`,
  `cl_FEM_DofMgr_DofData.cpp:3570-3643`) use the same +1/−1 convention for both
  surfaces; the top and bottom paths are symmetric in their sign logic.
- **(c) T-matrix condensation.** `cl_Maxwell_TMatrix` is hard-wired to TRI6/TET10 and
  gated on `max_element_order() == 2` (`cl_MaxwellFactory.cpp:1591`) — inactive on this
  linear mesh (Codex). The active linear condensation recovers hanging dofs as weighted
  source sums after every accepted update (`cl_FEM_DofMgr_SolverData.cpp:2641-2657`).
- **(d) Postprocessor layer/side handling.** One `ThinShellConductor` SPR instance
  handles all non-jc layers with the identical `compute_conductor` kernel
  (`cl_MaxwellPostprocessor.cpp:493-507`); no first/last-layer branch exists anywhere in
  the J path (Grok, confirmed). Block↔material assignment is 1:1; the exodus block labels
  come from that assignment. `copy_seam_fields` copies T/H/B only — J is excluded by
  design (`:529-599`).
- **(e) MPI.** Decisive data point: **every tape's layer elements are owned by exactly
  one rank** (read from the exodus `ElementOwner` field), so no partition boundary
  crosses any tape face — aura-skip truncation, ownership races, and rank-dependent
  patches are all off the table for the tape faces.
- **Wall/SPR cross-writes.** Wall nodes are separate exodus nodes; the wall's
  SideConnector instance claims only wall nodes; tape rim nodes are recovered from tape
  elements only (side-local claim, `cl_FEM_Postprocessor.cpp:265-282`).

**Supporting the verdict:**

- Layer-resolved forensics (§2): noise follows YBCO adjacency, not layer-stack boundary
  roles (master vs slave tie, air vs solder tie — every combination of those appears on
  both a clean and a noisy plane somewhere in the matrix).
- Standing-mode correlation 0.996–0.998 across saves with ~constant relative amplitude:
  a persistent discrete mode of the converged-to-floor solution, not random per-step
  iteration junk and not a postprocessing artifact (a picture artifact could not grow
  smoothly in amplitude while staying spatially frozen as the physics advances).
- Log forensics [corrected, see §3b]: from t ≈ 3.5 s onward the iteration terminates at
  −108…−117 dB — i.e. *at* the 1e-11 tolerance under the 10·log10 convention (the −156 dB
  at t = 0.01 s is the direct-solver roundoff of a then-linear problem); Newton stall
  signature near the floor; Δt collapse at the sigmoid's steepest point (the deck's own
  NOTE predicts trouble exactly there).

## 5. Real defects found incidentally (one confirmed and fixed; one RETRACTED)

1. **`Kernel::compute_element_volumes` thin-shell branch, misindexed SPR weights**
   (`src/fem/kernel/cl_FEM_Kernel.cpp:944-950`; found by Grok, mechanism confirmed by
   reading, then proven inert for this mesh by data). `tSurfaces` is indexed by shell
   facet order, but `f++` advances only for *owned* elements — on any rank whose owned
   elements are not a prefix of the block, every owned element after the first gap reads
   a *different facet's* area. `_Volumes` feeds the SPR weights
   (`cl_FEM_Postprocessor.cpp:825-839, 967-1029`). Wrong-but-positive weights leave
   linear fields exact and corrupt recovery only where the field has curvature — a
   perfect noise-camouflage landmine. **Inert here** only because each tape is
   single-rank-owned, tapes are contiguous in block order, and the eight tape meshes are
   identical translated copies (facet-area sequences match exactly — verified
   numerically), so the misindexed lookup happens to return the correct value.
   **Fix (2 lines):** increment `f` for every element, not only owned ones — i.e. move
   `f++` out of the ownership branch (or index `tSurfaces` by the loop position).
   Expected effect on this run: none (hence "inert"), but it must land before any
   non-uniform tape mesh or finer partition is run.
2. ~~**Tolerance not enforced / accepted-at-floor semantics.**~~ **RETRACTED
   2026-08-15 morning** — an artifact of misreading the log's dB scale as 20·log10; the
   code prints 10·log10 (`cl_FEM_Controller.cpp:1469`), so the −108…−117 dB exits are
   1.6·10⁻¹¹…2·10⁻¹² and the 1e-11 tolerance IS enforced (see the correction in §3b).
   What survives of this item is only the efficiency observation: Newton burning 4–5
   iterations at a bit-identical residual (~2·10⁻¹¹, near the attainable floor) before
   Picard closes the step — worth a look someday, not a correctness defect.

Also worth recording: the module docs in this area have drifted (Grok's audit lists the
stale citations: recovery-theory "Mode 1" line ranges, `compute_conductor` line numbers,
the node-sharing doc's `src/mesh/` factory path, and
`side_connector_wall_element.md`'s "WIP stop" claim vs the built-and-running walls).

## 6. Proposed actions (nothing urgent for the running job)

There is **no minimal code fix that would change these two screenshots**, because the
screenshots are showing (i) a physically asymmetric stack rendered correctly and (ii) a
standing oscillation of the 1e-11-converged solution on the YBCO-side interfaces.
Concretely:

- [x] **A1: land the `_Volumes` indexing fix** (`cl_FEM_Kernel.cpp:944-950`, move the
      `f` increment out of the ownership test). Two lines; no effect on this run's
      results, prevents silent SPR corruption on any future non-uniform shell mesh.
      *(Applied 2026-08-15 with Christian's approval, plus a size-precondition
      `BELFEM_ASSERT` tying block element count to shell facet count; syntax-checked
      with the build tree's flags. **Both auditors approved the fix**: Codex confirmed
      the order invariant for all four `create_elements_on_blocks_*` variants, the
      order-2 layer/block arithmetic, and that no other consumer relied on the old
      compressed index; Grok confirmed completeness and recommended landing it.
      Executable gate = Christian's rebuild + `make check-fast`.)*
- [x] **A1b: `compute_element_volumes` MPI buffer sizing** (found by Codex's audit of
      A1, verified by reading, applied same session): rank 0 sized every remote rank's
      not-my-volumes return buffer with **its own** `tNumNotMyElements`
      (`cl_FEM_Kernel.cpp`, redistribution loop) while filling it with rank p's count —
      buffer overrun whenever rank 0 owns more elements than rank p. Fixed to
      `tProcVolumes.set_size( tProcIDs.length() )`, which matches both the fill loop
      and the receiver's expectation. Syntax-checked together with A1.
- ~~**A2: probe the acceptance path.**~~ **Obsolete** — its premise was the retracted
      §5 item 2; the acceptance path enforces the tolerance as designed. The surviving
      slice (Newton's bit-identical-residual grind near the floor) is folded into the
      retraction note in §5.
- [ ] **A3: ghost-penalty calibration check for the Ag|YBCO interface in deep
      saturation** (`h_ghost`, `mt_maxwell_h.cpp:324-491`; hardcoded `eta = 4.0`,
      `k_reg = 1e-3 Ω`). The standing interface mode lives exactly on the interfaces this
      term controls. A cheap 2D or single-tape sweep of `eta` (4 → 10 → 20) at fixed
      saturation would show whether the mode amplitude responds. Physics call —
      Prof. Sirous / Christian territory, not an AI edit.
- [ ] **A4: extend the edge-function circulation battery to the TS family**
      (QUAD4TS/PENTA6TS/HEX8TB) so the class of defect this investigation *looked for*
      is permanently gated. My hand-check of PENTA6TS (§4) can seed the expected values.
- [ ] **A5: doc refresh pass** over the four stale citations listed in §5.
- [ ] **A7: JJC display bias — diagnose, then fix (postproc-only, pre-release).**
      The written J/Jc (≈ 1.3–1.4) overstates the V–I-implied true overcriticality
      (≈ 1.10) by ~30 %; JJC is display-only (sole writer: `cl_MaxwellPostprocessor`;
      no functional consumer), so this does not affect solutions — but it is a
      headline output for HTS users and should be closed before the August release.
      Step 1, diagnostic probe (scratchpad-probe pattern against the prebuilt libs):
      evaluate jc through BELFEM's own `JcFunction_Database` interpolation at a few
      nodes' recovered (B, β, T) from `s.00119` and compare with the written
      denominator |J|/JJC at the same nodes. Open wrinkle it must resolve: the
      implied jc_used ≈ 1.8·10¹⁰ vs my external nearest-node table read
      ≈ 2.8–3.6·10¹⁰ at the recovered median B — attribution (my read vs recovered
      inputs vs evaluation path) is undecided. Step 2, the fix — only after step 1:
      if the inputs or evaluation path carry the bias, route
      `compute_superconductor_ts`'s remaining field evaluation through the
      assembly's own MaxwellData helpers (the `compute_side_connector` precedent),
      so display and physics share one evaluation path; rewrite the deck's stale
      jc/Ic comment block in the same pass. **The angle-convention part is DONE
      (2026-08-15, Christian's approval): `compute_superconductor_ts` now calls
      `Calculator::bn_angle` instead of its own signed `acos`
      (`cl_MaxwellPostprocessor.cpp`, the β line) — display and assembly share the
      folded |n·b| convention and its small-field guard. Display-only; expected
      JJC change ≤ ~1 % (jc table symmetric); the ~+30 % question stays open until
      the step-1 probe runs.**
- [ ] **A6 (audit by-catch, report-only):** three pre-existing residual risks flagged
      during the fix audit, none introduced or cured by A1/A1b: (i) on non-root ranks
      the distributed sub-mesh's positional facet↔element pairing rests on the
      distributor's non-stable `std::sort` by block index (Grok) — the robust fix is
      identity pairing (element ↔ parent facet), and the new assert catches count
      mismatches but not order mismatches; it is also compiled out in release; (ii)
      the thin-shell counter's `else` branch counts all non-owned elements while the
      pipette branch and the pack loop guard on `owner < mCommSize` (Grok); (iii)
      BFM/proto import paths populate `ThinShell::blocks()` from saved metadata
      without a facet-count check, so malformed imports would trip the new assert —
      which is what it is for (Codex).

If the physics outcome (3a) is itself unexpected — i.e. the buffer was *not* meant to be
insulating, or the halves were meant to be ghost-coupled — then the deck (or
`create_buffers`) is where to intervene, e.g. giving the buffer a large-but-finite rho.
That is a modeling decision, not a bug fix.

## 7. Verification plan for the morning (cheap, in order)

- [ ] **V1 (5 min, no rerun):** view the current last save with the *element*-level
      check: in ParaView, apply a per-block threshold and compare plane-7 vs plane-0
      surfaces of `tape_70_copper`/`tape_76_copper` — confirms the §2 matrix visually.
      (My renders: scratchpad `j_faces.png`; regeneration script embedded in the devlog.)
- [ ] **V2 (one short rerun, 1 rank, coarse save):** run a few saturated steps with
      `mpirun -n 1`. Prediction: **noise unchanged** (MPI/SPR ruled out); if the noise
      vanishes, I am wrong and the `_Volumes`/aura analysis must be revisited.
- [ ] **V3 (one short rerun):** same steps with `edge coating : off`. Prediction (from
      3a): the hastelloy-half blocks (L4/L5/L6) lose their transport current
      (fed only through the walls within a tape; solder still feeds them across tapes —
      so expect a *reduction*, clearest on tape 8 whose top half has no solder above).
      Confirms/refutes the two-halves electrical topology directly.
- [ ] **V4 (no rerun):** dump the edge-on-node constraint rows with the existing
      `BELFEM_PROBE_FUSED_ROWS=1` env probe (`cl_FEM_DofMgr_DofData.cpp:3645-3698`) on a
      1-rank restart and spot-check sign coherence on tape 1 plane 0 vs tape 8 plane 7 —
      closes the last static-only gap in the tie audit.
- [ ] **V5 (discretization experiment, longer; premise updated with the §3b
      correction):** the tolerance is already met, so more iterations will not help.
      Instead, rerun the saturation window with ONE knob changed at a time and watch the
      plane-0/1 roughness: (a) a tighter `maximum timestep` cap through the knee, (b) a
      finer z-resolution at the tape edges (`tapeResolutionTip`), or (c) the `h_ghost`
      `eta` sweep of A3. Whichever knob moves the roughness identifies the mode's
      controlling discretization.
- [ ] **V6:** "fixed" for the original question means: the noisy face is understood, not
      repaired — acceptance = V2 unchanged + V3 confirming the half-split + V5 naming the
      controlling knob. If no V5 knob responds *and* V2/V3 hold, escalate A3 (ghost
      penalty) to a real calibration campaign.

## 8. Codex and Grok reconciliation

Both audits ran blind (Codex from the formulation end, Grok from the postprocessor end;
neither saw my hypotheses; Grok was instructed not to read Codex's entry). Full text in
`tmp/ai_exchange/tapestack3d_j_asymmetry.md` (ephemeral; conclusions lifted here).

**Agreement (all three):** no EF_PENTA6TS sign/pairing defect (the 2D bug does not
recur); T-matrix inactive on this mesh; no first/last-layer branch in the postprocessor
J path; `copy_seam_fields` cannot touch J; the node-sharing leak is real but cannot
paint the outer faces.

**Divergence and resolution:**
- *Codex* ranked "nodal SPR/visualization artifact" first and "side-connector rim
  interaction" second. The data forensics (standing correlated mode, single-rank tape
  ownership, clean linear-regime early steps) demoted the pure-picture hypothesis: the
  mode lives in the solution. Codex's instinct was half right — the *inner*-plane
  pictures are indeed artifacts (the known leak), but the outer-face fingers are not.
- *Grok* refuted the symmetry premise (buffer split) — this survived code confirmation
  and became the backbone of the verdict (3a) — and flagged the `_Volumes` indexing bug,
  which I confirmed by reading and then proved inert for this particular mesh/partition
  by reconstructing the corrupted weights from `ElementOwner`. Grok's top hypothesis
  ("expected physics of the U-jacket-fed substrate half") is close to the final verdict
  but underweighted the log/exodus evidence (Grok read neither); after the §3b
  correction its "converged-solution feature" framing turned out closest to the mark.
- *Claude (this report)* adds what neither auditor had: the plane/tape noise matrix, the
  standing-mode correlation, the tape-ownership fact, the residual log forensics, and
  the hand verification of the PENTA6TS curl — and contributed the one confirmed error
  of the round (the dB misreading, caught and corrected the next morning; see §3b).

**Confidence:**
- Premise refuted / two-halves topology (3a): **high** — code-confirmed at three sites,
  consistent with the measured per-layer means.
- Noise is solution-side, not picture-side (3b, first half): **high** — standing
  correlated growing mode + single-rank tapes + clean plane 7 from the same SPR pass.
- Noise mechanism (3b-ii, as closed by the morning probe): YBCO discrete-front
  raggedness imprinted into the stabilizer at ×n_eff: **high (~85 %)** — three
  independent lines (8-tape amplitude tracking incl. the tape-2 resolution, point-wise
  correlation up to 0.70, slope ≈ n_eff from the independent V–I estimate). The
  residual 15 %: the point-wise decorrelation on the upper tapes is explained but not
  proven, and the mesh-refinement prediction (V5b) has not been run. The earlier ~70 %
  "unreached tolerance" mechanism was retracted outright — a calibration lesson.
- `_Volumes` defect real but inert here: **high** (mechanism by reading; inertness by
  numerical reconstruction from the run's own ownership data).

## 9. Open questions

1. ~~Why is tape 2 clean in every plane?~~ **Resolved by the morning probe (§3b-ii):**
   tapes 2–4 carry a spatially uniform J/Jc ≈ 1.1 (spread ≤ 0.05), so there is no
   front structure to imprint; the copper ripple tracks the per-tape YBCO J/Jc spread
   across all eight tapes. (Residual sub-question: the converged J/Jc distribution is
   2–4-uniform vs 5–7-structured — plausibly physical, since each tape's YBCO sits
   off-center toward −y, making the sheet array asymmetric about the stack midplane.)
2. ~~The exact exit mechanism of the magnetic loop at the −110 dB floor.~~ Resolved
   by the §3b correction: the loop exits because the tolerance is met. Remaining slice:
   the Newton bit-identical-residual grind near the floor (efficiency only).
3. Whether the buffer split is intended modeling (export the answer into
   `doc/input_file_reference.md` — a deck author reading `buffer : 0.2 mum` has no way
   to know that a rho-less material electrically severs the tape).
4. Hastelloy plane-5 nodal J (~7.9·10³ A/m²) is ~30× above E/ρ_hastelloy; most likely
   this is the mixed-plane write (silver side wins at shared nodes), but this pass did
   not check the hastelloy current budget end-to-end.
