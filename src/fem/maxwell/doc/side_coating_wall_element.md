# Side Coatings: the HEX8TB Wall Element {#fem_maxwell_side_coating_wall_element}

**Date:** 2026-07-31 (updated 2026-08-27: renamed from `side_connector_wall_element.md`,
§6.1 thermal ruling added)
**Purpose:** Reference for the side-coating campaign, covering the physics of HTS
tape side conductivity, the failure of earlier constructions, the collapsed-wall theory,
the HEX8TB wall element and its edge function, and the current implementation status.
Code identifiers retain the older "side connector" terminology
(`mCreateSideConnectors`, `h_side_connector`, `side_connector_blocks()`), while the
deck key and domain types use the coating terminology (`edge coating`,
`LeftCoating`/`RightCoating`).
**Module:** `fem/maxwell` (h-φ thin shells), `fem/kernel` (ThinShellFactory),
`fem/interpolation/nedelec` (EF_HEX8TB), `mesh` (element template)
**Status:** implemented end-to-end for the uncoated variant, including mesh
construction, the edge function, IWG + factory wiring, field recovery, and
postprocessing. The `edge coating : on` path runs through assembly in both serial and
multi-rank configurations. Verification gates and the coated variant remain open.
Excluding the walls from the thermal problem is a deliberate modeling decision, not
pending work (§6.1). See "Current status" at the end.

---

## 1. The physical problem

BELFEM models HTS coated-conductor tapes with the h-φ formulation: conductors carry the
magnetic field intensity **h** on curl-conforming Nédélec edge dofs, non-conducting regions
carry a scalar potential φ on nodes. Tapes are collapsed thin shells with one material
layer per virtual level: `Cu / Ag / Hastelloy / buffer / YBCO / Ag / Cu`. The **buffer**
(an oxide barrier required by manufacturing) is essentially insulating and is modeled
φ-only, splitting the stack into two electrically isolated conducting halves.

Real tapes are slit from wide strips and then surround-electroplated: copper deposits onto
the exposed layer cross-sections at the slit edge, forming a U-shaped jacket — top plate,
**edge wall**, bottom plate — one connected conductor. The physical question is:
**during a quench, can current escape around the insulating buffer through this slim edge
passage to the other side of the tape?** This governs current sharing in CORC cables:
tape-to-tape transfer passes through the substrate side and therefore *in series* through
the edge wall. It is not well understood experimentally, which is precisely why we
simulate it.

Two facts shape the model:

- **Injection into the jacket is distributed, not edge-localized.** During a quench, the
  current exits the REBCO vertically through the broad silver interface above it (measured
  Ag/YBCO contact resistances R_ct ≈ 10⁻¹²–10⁻¹⁰ Ω·m², see `contact_impedance_theory.md`
  §10.1) and reaches the wall through the Cu top plate. There is no 20 µm current-crowding
  field structure for a mesh to resolve.
- **The wall's internal physics is 1D ohmic and interface-controlled.** The pristine metal
  contributes ≈ 10⁻⁸ Ω·m of edge; the genuine unknowns are the plating-interface qualities
  (slit-cut Hastelloy oxide, REBCO edge, slitting damage), plausibly 10⁻⁶–10⁻⁴ Ω·m. These
  are *parameters*, not fields. The condensed series resistance of the wall path is

  ```
  r′ = ρ_Cu·h_path/t_w + ( R_ct,top-edge + R_ct,bottom-edge ) / t_w        [Ω·m]
  ```

  r′ is the calibration dial; the collapsed representation loses nothing that was ever
  known about the wall.

The limits behave correctly. r′ → ∞ recovers today's insulated-edge model; r′ → 0 is a
perfect short, and the tape then *behaves* as the single conductor of the pre-buffer
model. This is a limit behavior of the coupled two-conductor model, not a topology change:
the halves stay topologically separate, and the partition remains a solution unknown.

## 2. What failed before, and what the analysis established

Three geometric attempts failed for the same reason: tape width ≈ 4 mm with ~400 µm edge
elements against a ≈ 20 µm plated wall (20:1 length separation):

1. **Meshing the copper strip explicitly** destroys conditioning (20 µm slivers beside
   400 µm bulk).
2. **Extra edge dofs / enrichment** cannot help: a dof whose only job is to represent a
   fraction-η feature (η = t_w/L ≈ 1/20, pinned by physics) is O(η)-degenerate with the
   dofs already present, and the wall runs along the *entire* edge, so the near-null
   directions assemble into a connected, unpivotable subspace along the side curve.
   Ghost-penalty/Nitsche stabilization presupposes a healthy same-phase neighbor, which
   this geometry does not offer — every perimeter element is cut identically. Ad-hoc
   H(curl) enrichment also falls out of the discrete de Rham complex (spurious modes; Monk
   2003, de Rham diagram; Arnold 2018).
3. **The first HEX8TS "wrap" element** (removed June 2026) shorted top and bottom plates
   through shared and *free* edges. It was unphysical: the wrap's binormal-H at the
   shell↔connector fold remained a free parameter, which the solver drove to ~0 in low-ρ
   copper, producing a ~2-order binormal-H discontinuity and suppressed side-curve J/J_c.

The analysis that followed established the correct viewpoint: **the wall's entire
information content is one number per unit edge length — its current.** In the
h-formulation a resistor is `R·I(h)·I(v)` where the current I(·) is a linear functional of
existing dofs (a trace jump or a circulation); the wall term is a curve integral

```
a_wall(h, v) = ∫_Γedge  r′(y) · [[h_∥]] · [[v_∥]]  dy
```

coupling the tangential traces on either side of the wall. No neighboring tape element's
curl is evaluated, and no sub-wall field is resolved. The wall element below carries its
own (trivially representable) curl, but it never asks the adjacent mesh to resolve the 20
µm scale. This "couple, don't enrich" principle governs everything below.

## 3. The current design: a degenerate wall element with fused traces

The revived side connector packages the collapsed-wall physics as a **degenerate volume
element** rather than a bare curve kernel. The element serves four purposes:

1. it uses the standard element/IWG assembly path (no special rank-one kernel);
2. it produces the wall's **current density vector j** in the cross-section for
   postprocessing and, later, Joule heating (`P′ = r′·|I′|²`);
3. it is the natural home of the **coated variant**, where one genuinely new internal dof
   per station is required (below);
4. it gives the wall a geometric representation in the output mesh.

The old wrap's failure mode is *designed away*: the element has **no binormal or normal
edge dofs at all** — only longitudinal traces, all of which are either shared with the
shell stack or fused/hanging. Nothing at the fold is free.

### 3.1 Fusing: two algorithms, selected per side-curve station

The do-nothing condition at an unmodeled lateral wall is n×e = 0 — a zero-impedance
lateral redistribution path (the r′ → 0 limit), *not* insulation. The physical slit edge
instead imposes, per station along the side curve:

- **Uncoated edge** (no plating, or modeling the bare slit): per-layer J·n = 0 plus
  tangential-H continuity with the single-valued air trace is an **essential** condition —
  all side traces fuse to the air trace `g = φ₀ − φ₁`: `h_0 = h_1 = … = h_N = g`. Pairwise
  fusing without the g anchor would float a common mode; all-traces-to-g is the correct
  fuse. Zero new unknowns.
- **Coated edge** (plated wall): the interior side traces fuse to **one new internal edge
  dof `h_in`** per station — the only new unknown — and `h_in` couples to g exclusively
  through the wall element's r′ stiffness. `h_in` must be a real dof owned by an active
  element and must **never** be tied or hung on g: that would short-circuit the wall
  impedance and is the one genuine overconditioning risk of the design. A coated station
  must never ship without its wall element (atomicity).

The uncoated fuse is the first implementation target; the coated variant is designed but
not implemented.

### 3.2 Mesh construction (implemented)

When side connectors are enabled (the factory flag `mCreateSideConnectors`, default
off), `ThinShellFactory::create()` builds one connector block per side curve
(`DomainType::LeftCoating` / `RightCoating`). Per curve it:

- computes the tangent/binormal frame of the curve (`compute_binomial_vectors`), anchoring
  the tape-side sign on a terminal curve — the search fails loudly if no anchor exists,
  which also catches closed-loop tapes (currently unsupported);
- resolves the station-to-edge mapping (`compute_side_edge_indices`) with an order-aware
  stride walk;
- builds a `SideLayer` per layer level: edge-driven node collection, robust against node
  duplication because temporary and layer edges reference original-normalized nodes, and
  the station-to-edge lookup keys are original-normalized as well (a side-curve station
  that is itself a cut duplicate resolves correctly);
  outer nodes and outer edges are created as **hanging** entities whose sources are
  mid-plane curve *nodes* (outer edges hang on their mid-plane edge's end nodes, not on
  the edge entity) — the T-matrix weights are computed later by the dof manager
  (`DofData::create_dofwise_t_matrices_master`);
- assembles one **HEX8TB** element per station and layer gap, edge-driven with a
  per-station orientation flag, so the node-container ordering is not load-bearing;
- builds one **recovery facet** per wall element into a hidden sideset (one sideset per
  side curve): a QUAD4 pairing the shell layer-block element (master, on its lateral
  face) with the wall element (slave, on its inner face). The facet's nodes are linked
  from the master face (canonical order, normal pointing shell → wall), the slave
  orientation is computed by node identity at construction, and the facet id is always
  the wall element's id plus one. These facets are the wall's only channel to the
  neighbor fields: the thickness lookup (§4.1) and the h_b/h_n recovery (§5) both route
  through them.

Current constraints are enforced by loud guards: first-order side curves only, and
connectors are built for **all** layer levels (selective layers are unsupported). The
former work-in-progress stop at the end of the factory path has been removed.
Construction now proceeds through assembly in both serial and multi-rank
configurations (§6).

Because the outer entities own no dofs, **periodicity is inherited through their
sources** — no partner pairing or backup registration is needed, and cross-block periodic
partners (twisted-helix meshes) cannot be corrupted by the connector construction.

### 3.3 The HEX8TB element template

`ElementTemplate<8,8,4,6,1>` (`src/mesh/cl_Element_HEX8TB.hpp`, enum value 125):

- **8 nodes in plain HEX8 ordering** (deliberate; the edge numbering scheme is where it
  differs from HEX8TS);
- **4 edge dofs** on the longitudinal edges (0,1), (3,2), (4,5), (7,6), chosen so that
  **all four point in +curve direction**;
- faces (0,1,2,3) and (4,5,6,7) lie on the two adjacent layer levels;
  `get_edges_of_facet` supports only these two facets;
- the MeshCheckers **error out** on a negative-volume HEX8TB instead of node-swapping — a
  swap would silently corrupt the slot-tied edge assignments. A negative volume means the
  factory built a left-handed element and must be fixed there.

The wall's two thin dimensions are physical model data, stored as block thicknesses; the
mesh node positions of the outer sheet are partly a **drawing proxy** (chosen for
visualization), a fact that drives the edge-function metric below.

## 4. The edge function EF_HEX8TB (implemented)

Reference frame: ξ along the side curve (tangent t), η across the wall width (binormal b),
ζ across the layer gap (normal n). The four reference traces are the ξ-directed subset of
the standard first-family hex Nédélec basis,

```
F_0 = ⅛(1−η)(1−ζ)     F_1 = ⅛(1+η)(1−ζ)     F_2 = ⅛(1−η)(1+ζ)     F_3 = ⅛(1+η)(1+ζ)
```

on the corners (η,ζ) = (−1,−1), (+1,−1), (−1,+1), (+1,+1) matching the four dof edges.

### 4.1 Metric: exact cuboid (3-AI consensus, unanimous)

The element uses **directions from geometry, magnitudes from the blocks**:

- unit tangent t from the two cross-section face midpoints, length L = |mid{1,2,6,5} −
  mid{0,3,4,7}|;
- unit normal n from the layer-face midpoint difference, orthogonalized against t;
- unit binormal b = n × t (right-handed triad t × b = n);
- width w (η) and thickness d (ζ) taken **exactly** from the mesh blocks — never from the
  node positions. w comes from the connector block's own thickness; d is recovered
  through the **recovery facet** (below): the facet always carries the wall element's id
  plus one, and its master is exactly the layer block the wall spans, so
  `mesh->facet(id+1)->master()->block_id()` yields the unambiguous gap thickness. The
  element's physical tag stays free for the material machinery.

Rationale, in order of importance: (i) the node offsets of the outer sheet are a drawing
proxy, so a node-based metric would silently couple the assembled wall conductance to the
visualization; (ii) the block scalars are exact, giving a Jacobian that is constant per
element with zero per-point noise — differencing global coordinates for µm-scale spans
costs relative accuracy, and a Gram-matrix path would carry a condition number of order
(L/w)²; (iii) it is the same "exact collapsed dimension" principle EF_HEX8TS already
applies to its thin direction, extended to two directions.

The cuboid discards in-plane shear/taper of the proxy hexahedron and chord-vs-arc
length. These effects are of the same order as the mesh faceting error everywhere else,
and the physical wall is a rectangular strip by definition.

### 4.2 Geometry Jacobian and Nabla

With BELFEM's row convention J(i,:) = ∂x/∂ξᵢ:

```
J = [ (L/2)·tᵀ ]                Nabla = J⁻¹ = [ (2/L)t | (2/w)b | (2/d)n ]
    [ (w/2)·bᵀ ]
    [ (d/2)·nᵀ ]                detJ = |detJ| = L·w·d/8 = V/8   (constant, > 0)
```

Because the triad is orthonormal, the inverse is analytic: no Gram matrix, no LU.
Everything is computed once in `link()`; `update_nabla()` is a no-op that exists only
because the Calculator calls it before reading `det_J()`. Since the edge function defines
its own right-handed frame, the sign of a Jacobian built naively from node coordinates is
irrelevant to the operators; the MeshChecker guard still catches genuinely inverted
proxies at the factory level.

### 4.3 E and C operators

All four dofs share the same physical direction ∇ξ = (2/L)t:

```
E(:,k) = s_k · F_k(η,ζ) · (2/L) · t                      (h = E ψ = h_t·t only)

C(:,k) = s_k · [ (4/(L·d)) ∂F_k/∂ζ · b − (4/(L·w)) ∂F_k/∂η · n ]     (j = C ψ)
```

so in the tape frame

```
j_t = 0        j_b = ∂h_t/∂s_n        j_n = −∂h_t/∂s_b .
```

Unit circulation holds: on its own edge, F_k = ½; multiplied by (2/L) and integrated over
the physical edge length L, it gives ∫_edge E_k·dl = s_k δ_jk.

**The clean mental model: h_t is the stream function of the cross-flow.** Since t is
constant, j = ∇×(h_t t) = ∇h_t × t — the current flows along level lines of h_t in the
(b,n) cross-section and is divergence-free by construction. The commutation path is a
corner-turning flow: it enters binormally from layer j, climbs normally across the gap,
and exits binormally into layer j+1. j_b and j_n are the same current metered through
different control surfaces. The total current transferred between two boundary points of
the cross-section equals the h_t difference between them — which is exactly the r′ stencil
on the edge traces. The wall stiffness ∫ρ j·j dV must therefore reduce, after
cross-section integration, to r′·(Δh_t)²-type terms; matching the element's effective
resistivity to the physical r′ is fixed at IWG wiring time.

**j_t ≡ 0 is deliberate and correct.** The basis is F(η,ζ)∇ξ, whose curl lies entirely in
the (b,n) plane. Longitudinal transport current belongs to the thin shell; a wall j_t
would double-count it. The IWG must never ask the wall element for tangential current.

### 4.4 Parametrization consistency across element families

The hex family parametrizes edges over [−1,1]; the simplex family (TRI3/TET4 Whitney
functions λ_i∇λ_j − λ_j∇λ_i) uses barycentric coordinates over [0,1]. This causes no
inconsistency: the edge dof is the **circulation** ∫_edge E_k·dl — its unit is ampere, not
ampere per meter — and every BELFEM edge function is normalized to unit circulation, so
the reference parametrization cancels against the ∇ξ factors. Moreover, at lowest order
the tangential *trace* along an edge is constant (1/edge-length) in both families, so a
shared edge is conforming pointwise, not merely in the integral.

### 4.5 Geometry interpolation: HEX8 Lagrange on top of the cuboid metric

Decision (Christian, 2026-08-08): the wall element uses the standard HEX8 Lagrange
functions for nodal interpolation — the interpolation factory maps HEX8TB to the same
`InterpolationFunctionTemplate<HEX, LAGRANGE, 3, 8>` as HEX8 and HEX8TS (the nodes
follow the plain HEX8 convention, so shape index k ↔ node k holds) — while ALL spatial
derivatives and integration weights come from the exact-cuboid metric of §4.1/§4.2.
No new interpolation class is needed.

The two layers meet without touching: plain N(ξ) evaluation is Jacobian-free, and the
assembly kernel consumes only the edge-function operators E, C and dV = V/8 — the
Calculator's dV dispatch prefers the edge function's determinant whenever one is
attached, so a determinant from node coordinates is never taken on this path.

**The committed variational crime, and why it is the lesser of two evils.** On a curved
tape the physical wall is slightly warped and twisted; the cuboid metric ignores this,
while the Lagrange interpolation sees the true node coordinates. The alternative — the
exact isoparametric Jacobian of the proxy hexahedron — was considered and rejected:

- It is not a drop-in substitution. The basis E_k = s_k·F_k·∇ξ and the properties that
  make the element correct (unit circulation, j_t ≡ 0, divergence-free cross-flow) hold
  in the cuboid frame. Honoring an exact Jacobian means the covariant (Piola) transform
  of the basis on a warped degenerate hexahedron — a redesign of the element, not a
  parameter change.
- The transform re-imports the conditioning problem of §4.1: at aspect ratios
  w/L ~ 10⁻²…10⁻⁴, per-point inversion amplifies coordinate noise by O((L/w)²) — the
  artifact class this element exists to avoid — and the positivity of det J would again
  depend on node placement, which the MeshChecker deliberately refuses to repair by
  swapping.
- The payoff would be geometric fidelity the model cannot use. The wall represents a
  lumped transfer resistance whose effective ρ is calibrated against the physical r′
  (§4.3); metric-level constants are absorbed by that calibration.
- The committed error is controlled: second order in w/L and d/L, first order in the
  twist per element along the curve. Through the width no refinement ever occurs — the
  element IS the lumped model — while along the curve, mesh refinement drives the
  per-element twist to zero. The crime therefore does not obstruct convergence in the
  one direction that is actually refined.

Two standing cautions: (i) the argument rests on the sliver aspect ratio — HEX8TB must
not be reused for elements that are not thin in two directions; (ii) any future consumer
needing spatial derivatives of nodal fields on a wall element must take them through the
edge-function Nabla, never through the inverse of the node-coordinate Jacobian — dNdX on
connector groups is undefined by design.

Precedent: the PENTA6TS/QUAD4TS thin shells already split metric from interpolation the
same way — the edge function caches the analytic volume element, and the Lagrange
functions interpolate nodal data only.

## 5. Field recovery for the material law (implemented 2026-08-05, with the wall kernel)

The recovery below shipped together with the `h_side_connector` kernel — the wall
assembles M(µ₀) + K(ρ CᵀC) with the material law evaluated on the recovered
`h = h_t·t + h_b·b + h_n·n`, routed through the master-block/sideset hop described
at the end of this section.

The E operator provides only h_t, but the HTS material law ρ(‖h‖, angle) needs the full
field vector. The thin shell already solves the analogous problem for its normal component
(`compute_hn`, `src/fem/kernel/cl_FEM_Calculator.hpp`: average the master/slave volume
traces — ∇φ on a φ-region, the conductor's own `E * q` on an h-conductor — and
project onto the facet normal). The wall needs **two** recovered components; the agreed
design is:

- **h_n (layer direction):** reuse the shell recovery. Average the recovered h_n of the
  **two layers bounding this wall gap** (never a stack average), project on n. Constant
  per station.
- **h_b (binormal):** the adjacent shells' solved **in-plane field**, projected on b at
  the edge station, is the primary donor. It is tangential to the shell and therefore the
  physically consistent source, and it is smoother than air ∇φ evaluated at the tape-edge
  corner singularity. Air −∇φ serves as fallback/diagnostic only.
- Across the wall cross-section, blend the four per-slot donor values with the bilinear
  corner weights N_k = 2F_k — consistent with the lowest-order element.

Two hard rules, distilled from the old wrap failure:

1. **Never add binormal or normal edge dofs to "improve" recovery** — that is exactly the
   free-parameter pathology that killed the HEX8TS wrap.
2. **Never feed recovered h_b/h_n into the curl/stiffness operator** — they enter the
   material law evaluation only: `h_mat = h_t·t + h_b·b + h_n·n`.

Implementation note: the wall element's faces have no volume neighbors
(the proxy sits inside the air region and its outer nodes hang on mid-plane sources), so
recovery must be routed through the **side-curve facets and calculators of the adjacent
shells**, not through face-neighbor lookups on the HEX8TB itself.

### 5.1 Postprocessing convention: T, H, B are tape-seam copies; J is native

The wall's **output** fields follow a copy convention, not the recovery above
(which serves the assembly-side material law only):

- **T, H and B on wall nodes are the tape seam values**, copied station-ordered
  from the adjacent layer elements (wall node → `original()` → layer rim node)
  onto both wall faces by `MaxwellPostprocessor::copy_seam_fields()`. Physical
  justification: the strip is thin pure metal (µᵣ = 1), so across its few-µm
  width the field equals the seam field to sub-mT accuracy, and B = µ₀H adds no
  independent information. Evaluating the edge-function recovery at nodal
  parametric points of these extreme-aspect slivers is numerically fragile and
  is deliberately not used for output.
- **J = C·q is computed from the wall's own edge dofs** — the edge current
  transfer is the wall's payload and is never copied.
- The copy runs master-only in the postprocessor (after the tape instances have
  recovered their fields): assembly-side field writes are rank-local, and the
  save-time field gather collects dof-carrying entities only, so values written
  during assembly would be lost on parallel runs.

## 6. Current status

Implemented (as of 2026-08-24 — the `edge coating : on` path runs end-to-end
through assembly, serial and multi-rank):

- Mesh construction: side-curve frames, station/edge resolution, SideLayer, hanging outer
  entities, connector blocks with HEX8TB elements; ownership chain verified (mesh owns
  blocks, blocks own elements, ThinShell holds non-owning references). The former
  work-in-progress stop is removed.
- Element template + mesh-level helper coverage (meshtools, enum/to_string, both
  MeshCheckers with the loud no-swap policy, CurvedElementChecker, vtktools).
- Edge function `EF_HEX8TB` (`src/fem/interpolation/nedelec/cl_EF_HEX8TB.{hpp,cpp}`):
  exact-cuboid metric, constant Jacobian/Nabla, E and C per §4; factory case in
  `cl_EdgeFunctionFactory.cpp`; `fn_num_nedelec_dofs.hpp` returns 4.
- **IWG + factory wiring:** the `h_side_connector` kernel assembles M(µ₀) + K(ρ CᵀC)
  with the §5 field recovery and the clamped ρ(T,‖B‖,β) material law. The wiring also
  covers block dispatch on the connector domain types, FieldList dof entries
  (connectors share the conductor edge-h table), group activation, and deck control
  via the `edge coating` / `edge coating width` input keys.
- **Recovery-facet wiring on both parallel variants:** the wall element finds its
  recovery facet through the id+1 invariant (`mesh->facet( element->id() + 1 )`)
  and its master element via `set_master`; the multi-rank path ships the connector
  id list rank-0 → all and re-derives facet + master locally, with always-active
  errors on the owned-element failure paths. The coating block types travel by an
  explicit broadcast, so non-root ranks dispatch the walls correctly.
- **Geometry interpolation:** HEX8 Lagrange reused via case labels; the exact-cuboid
  metric keeps its own derivatives (§4.5).
- **Postprocessor:** connector fields in the output lists; T/H/B seam-copy and
  native J per §5.1 (`copy_seam_fields`, master-only).
- **`.bfm` round trip:** HEX8TB elements/blocks incl. domain type, recovery facets,
  hidden-sideset flags, hanging outer entities with sources and weights, edge slots.
  The mesh-configuration stamp covers the connector settings, so changing them
  rebuilds a cached mesh automatically.
- **Numeric exercise:** coupled magnetic+thermal production runs with
  `edge coating : on` (16-wall tape-stack case) construct the connector geometry
  and run real physics through it.

Pending:

- **Verification gates (§4.3 calibration included):** numeric handedness/orientation
  checks on both connector signs; a single-tape open-curve case comparing wall
  current against analytic r′; multilayer trace probes; a twisted-helix periodic
  regression; and a cut-terminates-on-side-curve station test.
- **Format guard:** the `.bfm` loader does not yet check a format version against
  its own capability — a stale reader given a newer file fails downstream rather
  than at the header.
- **`Element_HEX8TB::get_edges_of_facet`** supports only facets 4/5; harmless until
  facet-based boundary conditions touch connector blocks.
- **Parallel deferreds:** full aura facet shipping (assembly runs owned elements
  only, so aura copies without a shipped facet are skipped by design for now).
- **Coated variant** (`h_in`) and **closed-loop tapes** (no terminal anchor for the
  binormal sign; needs a facet-based anchor and a periodic tangent stencil).

Not pending — ruled out of scope:

- **Wall participation in the thermal problem**, including conduction, heat capacity,
  and the Joule source. Its exclusion is a deliberate and justified modeling choice,
  not missing work. Section 6.1 records the reasoning and the conditions that require
  this decision to be revisited.

### 6.1 Ruling (2026-08-27): exclude the walls from the thermal problem

The coating domain types do not appear in `src/fem/thermal/`. The thermal factory
selects only `Conductor | ThinShell | Ferro | Buffer` blocks, and the dof manager skips
the connector blocks outside the Maxwell problem. This exclusion is intentional. The
walls are electrically necessary because the φ-only buffer separates the tape into two
conducting halves and the wall provides the only inter-half current path. They are
thermally redundant, however, and the one omitted term is smaller than the model's
deliberate thermal approximations. In detail, for the reference geometry (4 mm tape, ≈96 µm stack, 20 µm
copper walls):

1. **Rim conduction is redundant.** Inter-half heat crosses the broad faces. Although
   the buffer blocks current, it conducts heat, and its face conductance per unit tape
   length (λ_MgO·w/t_buffer) exceeds that of the copper rim path (2·λ_Cu·t_w/h) by two
   to three orders of magnitude. This face path is assembled: adjacent layer wafers
   share their interface node set, and buffer blocks run the φ thermal kernel with the
   magnesia conductivity law. The stack also lumps through-thickness on
   microsecond timescales, far below any quench clock.
2. **Wall heat capacity is a ~1% correction** (two 20 µm walls on a 4 mm tape).
3. **ρ(T) feedback already exists without thermal wall dofs.** The wall's material law
   evaluates ρ_Cu(T, ‖B‖, β) at a temperature bilinearly interpolated from the four
   seam supports (`MaxwellData::compute_T_side_connector`); the wall is thermally
   slaved to the tape rim it touches, so an independent wall temperature would add a
   state the physics does not have.
4. **The omitted Joule source is smaller than the error introduced by the adiabatic
   boundary assumption.** Wall dissipation is assembled in the magnetic stiffness but
   consumed by no thermal kernel. With the implemented pure-copper wall
   (r′ ≈ 10⁻⁸ Ω·m), the loss is localized at the transfer front. It is of order 0.1 W
   per front, compared with ~kW/m of stabilizer dissipation; it accounts for under 1%
   of the deposited heat and is zero wherever current sharing is developed. Meanwhile,
   the thermal problem is adiabatic. No convection boundary condition exists (the
   thermal BC factory offers Dirichlet and Neumann only), so the model retains all
   deposited heat. A real LN₂-wetted stack, by contrast, has a thermal time constant
   of order 0.1–1 s against the bath. Over the one-to-tens-of-seconds transients run
   by these decks, the neglected cooling is of the same order as the total dissipated
   energy. It is roughly two orders of magnitude larger than the neglected wall share
   and is conservative where the wall omission is unconservative. Adding a ~1%
   unconservative correction beneath a deliberate ~100% conservative assumption would
   provide no practical benefit.

**Known accounting consequence:** magnetic-side and thermal-side loss totals differ
by the wall share. Consequently, comparing Joule power from the two sides will reveal
a percent-scale gap in `edge coating : on` runs.

**Revisit triggers: both conditions must hold before this ruling is reconsidered:**

1. a cooling boundary condition (convection/Robin to the coolant) is implemented,
   shrinking the conservative margin the wall term currently hides under, **and**
2. r′ is calibrated into the edge-interface band (10⁻⁶–10⁻⁴ Ω·m of edge, §1), which
   makes the wall loss locally first-order at transfer fronts — the regime of the CORC
   tape-to-tape current-sharing studies this element was built for, where the transfer
   current is in series with the wall and an edge hot spot feeds back into jc.

If both conditions are met, the remedy is to deposit `P′ = r′·|I′|²` onto the
existing seam/tape thermal nodes. Even then, explicitly meshing the walls for the
thermal problem remains out of scope: items 1–3 above are geometry-independent, so
thermal wall elements would add cost and dofs while representing physics already
captured by the face path and seam temperature.

## 7. References

**Code:** `src/fem/kernel/cl_ThinShellFactory.{hpp,cpp}` (construction),
`src/mesh/cl_Element_HEX8TB.hpp` (template),
`src/fem/interpolation/nedelec/cl_EF_HEX8TB.{hpp,cpp}` (edge function),
`src/fem/interpolation/nedelec/cl_EF_HEX8TS.cpp` (thin-shell reference implementation),
`src/fem/kernel/cl_FEM_Calculator.hpp` (`compute_hn` recovery precedent).

**Module docs:** `contact_impedance_theory.md` (collapsed-interface kernel family, measured
R_ct), `buffer_cut_topology.md` (cut topology of the split stack),
`ghost_penalty_stabilization.md` (the normal-direction stabilizer, and why it cannot help
in-plane), `thin_shell_virtual_domains.md` (virtual domains and edge allocation).

**Literature:** Messe et al. 2023 — BELFEM h-φ architecture, static condensation;
Alves et al. 2022a — 3D thin-shell TSA term split, lateral entity treatment;
Schnaubelt et al. 2023 — collapsed contact layer, free current coefficient;
Dular et al. 2021 — lateral no-current as essential constraint; Monk 2003 —
H(curl) conformity and first-family hex Nédélec spaces; Arnold 2018 — structure
preservation (why ad-hoc enrichment fails); Bathe 2016 §8.2.6 — round-off/conditioning
budget.
