# Non-tight cohomology generators: three ways to impose the current, and the cut-free formulation in BELFEM

**Date:** 2026-09-02
**Purpose:** The corc_solder decks produce a cohomology generator that is a valid unit cocycle but not
a tight one-sheet surface; the thin-cut factory then destroys it (guard added 2026-09-02, so the run
now aborts instead of solving without a jump). This plan explains the three ways out — tighten the
representative, drop the thin cut in favour of the cocycle as an edge basis function (Pellikka et al.
2013 eq. 4.5), or make the thin cut multi-sheet — and works out what the second one would look like in
BELFEM, since that construction has never been written down for this code base.
**Module:** `src/homology` (CutProcessor, cleaner), `src/fem/maxwell`, `src/fem/kernel`
**AIs involved:** Claude (analysis, plan), Codex gpt-5.6-terra/high (code inventory), Grok grok-4.6/high (formulation)
**Status:** PLAN — drafted 2026-09-02, not scheduled. **Reviewed 2026-09-03 by Gregory Giard, who owns
`src/homology`: he raises two further motivations the draft missed (integer coefficients, element shapes,
§2.1) and rules that BELFEM should eventually carry the thick cut but not now — the majority of problems
work, and the priority is to know the limits rather than lift them (agreed earlier with Frederic). No code.
The cohomology core (`cl_Cohomology`, `cl_Homology`, `cl_SimplicialComplex`, `cl_Chain`, `cl_Cochain`,
`fn_Smith`) is closed to AI edits; every option below lives outside it.

> **Scope guards:**
> - Nothing here changes how generators are computed (Smith form, pairing, `clean_spfa`); only how a
>   given generator is turned into a current constraint.
> - Option 2 must reproduce the present thin-cut solution on a deck where the thin cut is healthy
>   (tapestack3d) to solver tolerance before it replaces anything.
> - Periodic decks are in scope (the cochain already carries periodic mates); 2D is out of scope for
>   the first implementation.

---

## 1. Current behaviour and how it fails

| step | code | what it does |
|---|---|---|
| generator | `Cohomology` (core) | a unit 1-cocycle z on the air edges, δz = 0 on every face (`Cohomology::check()`, `cl_Cohomology.cpp:905`), paired ±1 with the suggested homology cycle |
| thick → thin | `CutProcessor::collect_facets` 3D branch `cl_CutProcessor.cpp:411-599` | per crossed tet: 3 flagged edges sharing a vertex = vertex cap, 4 = 2\|2 diagonal (`CutData::determine_cut_case_3d`); emit the face opposite the apex for **positive** caps only (:429-434); self-cancel (:478-488); drop φ-boundary faces (:489-497); dangling-face loop (:525-570); **new:** abort if the loop removed > ¼ or all faces |
| duplicates | `CutSet::create_duplicates` `cl_CutSet.cpp:42-82`; relink `cl_CutProcessor.cpp:1172-1255` | every node of a thin-cut face gets a duplicate whose sources are the original (weight 1) and the cut's abstract node (weight 1): φ' = φ + I |
| condensation | `DofData::create_dofwise_t_matrices_master` `cl_FEM_DofMgr_DofData.cpp:3478-3647` | hanging dofs replaced by their sources, element matrices become TᵀKT; the abstract node's dof I reaches the element only through these sources |
| interface | `MaxwellFactory::create_hanging_edges_and_facets` `cl_MaxwellFactory.cpp:1612-1713`, T-matrix :3630-3647 | conductor / thin-shell edge dof on the air interface hangs on the two air nodes: h_e = φ_{n0} − φ_{n1} |

**Why it fails on corc.** The positive-case faces close into one sheet exactly when the cochain,
restricted to its band of crossed tets, is the coboundary of a consistent side labelling (then they
are ∂U for U = tets touching the plus side, and 2\|2 tets leave no hole). The corc generator hugs the
conductor over the full circumference (4541 of 14447 edges within 0.3 mm of the tape surface, 10604
conflicts in the band 2-colouring): no side labelling exists, the emitted set has holes, the loop eats
it (7k → 242 faces), Ampère fails. Measured 2026-09-02 (`devlog/dl20260902_corc_strip_diagonal_and_d2_gate.md`
§5). The checker in the core is right: tightness is not a cohomological property, so it cannot test it.

**Bottom line:** BELFEM's current-imposition path silently requires a tight representative, and the
cleaner is not allowed to produce one (Tier-A pocket removal refuses any component touching a
sideset node, `cl_Cohomology.cpp:668-682`, :843-846, rule O1). Either the representative is tightened,
or the imposition stops needing tightness.

## 2. The three options and the recommendation

| option | idea | what changes | risk | effort |
|---|---|---|---|---|
| 1 tighten | replace z by z + δg with g integer on the nodes, minimising weighted support (min-cost tension), then the existing thin cut | new global optimisation in the cleaner (outside the core: a post-pass on the generator); factory unchanged | minimality alone does not prove the positive-case faces close (medium); needs area-like weights | M |
| **2 cut-free** | keep φ single-valued, add I·ψ with ψ = Σ_e z_e w_e (Whitney interpolant of the cocycle) to h in the crossed air elements | no duplicates, no thin cut; an extra dof I in crossed air elements; augmented φ kernel; interface edge relation gains + I z_e; postprocessing. **Also dissolves the non-unit-coefficient restriction and the tet-only restriction — §2.1** | a global dof in volume elements has no carrier today; fill-in ∝ band size; shells crossed on both sides | M |
| 3 multi-sheet thin cut | one duplicate per lift level, φ_λ = φ + λ I, relink by level | CutSet/relink rewrite | on corc the lift is not single-valued on the band (a loop around the cable lives inside the sheath), so it does not even exist there; with condensation it *is* option 2 | L |

**Recommendation.** Option 2 is the correct spine: it makes the imposed current independent of the
representative, which is the property the corc decks lack. Option 1 is worth doing anyway because a
tight z shrinks option 2's coupling column and keeps the thin-cut path usable for the decks where it
works; it is Phase 4 of `todo/thin_cut_nonunit_rectification_implementation.md`. Option 3 is rejected
(Grok Q5, Claude): after static condensation it is option 2 with more bookkeeping, and without
condensation it needs a lift that the corc band does not admit.

### 2.1 Two further reasons for option 2 (Gregory Giard, 2026-09-03)

The draft argued option 2 only from the corc collapse. Gregory adds two restrictions that are properties
of the thin-cut construction itself, not of any deck, and that option 2 removes for free. Both verified
in the tree on 2026-09-03.

**(a) Integer coefficients.** The thin cut can only represent a jump of ±1 per crossed element: `CutData`
holds the cochain as two ±1 bitsets, and `determine_cut_case_3d` (`cl_CutData.cpp:443-663`) classifies
one of seven ± patterns, aborting on anything else. Hence the whole rectification layer: `clean_spfa`
must find a unit representative or the run stops with *"generator N has no unit-coefficient
representative: no thin cut exists on this mesh … the mesh is too coarse along this loop; refine it
there and rerun"* (`cl_Cohomology.cpp:534-543`). Under option 2 nothing anywhere requires |z_e| ≤ 1:
ψ = Σ_e z_e w_e is curl-free because δz = 0 whatever the integers are, and ∮ψ = ⟨z,γ⟩ is the pairing
either way (§3.2). A generator with |z_e| ≥ 2 assembles exactly like a unit one — only the column
Σ_e z_e w_e changes. That retires the mesh-refinement demand, the SPFA feasibility abort, and the open
work in `todo/thin_cut_nonunit_rectification_implementation.md`. Gregory: *"it will definitely fix the
coefficient problem, since the thick cut can handle them directly in the equations."*

**(b) Element shapes.** The cohomology engine is shape-agnostic — it works on the simplicial complex of
any mesh, and Gregory has run pyramids. The thin-cut conversion is not: `CutProcessor` derives **one**
element type for the whole mesh from dimension and max order (`cl_CutProcessor.cpp:33-35`, again at
`:1166-1169`), and every dispatch has exactly four cases, TRI3 / TRI6 / TET4 / TET10, with
`BELFEM_ERROR( false, "This should not happen" )` as default (`:683`, `:802`, `:1220`). The cut cases
themselves are the 6-edge / 4-face tetrahedral pattern table (`cl_CutData.cpp:443-663`). So a pyramid,
prism or hex in the air region is not merely unsupported: because the type is global rather than
per-element, a mixed mesh would send such an element through `flip_node_bitsets_tet4`, which reads local
edges 0–5 of an element that has 8, 9 or 12 — most likely landing on the "inadmissible cut edge pattern"
abort, but by accident rather than by a clean check. Option 2 has no cut cases and no per-shape
bookkeeping: it needs only the Whitney edge functions of the element. Its reach is then bounded by the
Calculator's Nédélec dispatch, which today covers TRI3, TRI6, TET4, TET10, QUAD4, HEX8 and the thin-shell
types but **not** PYRA5 (`cl_FEM_Calculator.cpp:1100-1140`) — so option 2 would extend cuts to hexahedral
and prismatic air meshes, and pyramids would need one more edge-function entry (O6).

**Standing decision (Gregory, agreed with Frederic, relayed 2026-09-03):** BELFEM should eventually
implement the thick cut instead, but not now — the majority of problems work, and while the limits are
known it is a question of priorities. Implementing thin cuts was still worth it: the coefficient problem
had not been hit by anyone before, and the two constructions are now understood well enough to choose
between them. This plan is therefore a design record, not a scheduled work item; what is owed in the
meantime is that the limits stay visible, which is the 2026-09-02 collapse guard (§1) plus O6 below.

## 3. Option 2 explained: from the thin cut to the cohomology basis function

### 3.1 What the thin cut really encodes

In the air h = −∇φ. Around a conductor carrying I, ∮ h·dl = I, so φ cannot be single-valued: it must
jump by I somewhere on every loop around the conductor. The thin cut is that "somewhere": a surface S
with [φ]_S = I. BELFEM realises the jump by duplicating the nodes of S and tying φ' = φ + I through a
hanging dof whose second source is the generator's global unknown I (§1). The physics is in the jump,
not in the surface: any surface S in the same class gives the same h.

### 3.2 The same field without a surface

Write the multivalued φ as φ_single + I·λ, where λ is the integer "sheet number" that increases by one
each time one crosses S. Then h = −∇φ_single − I ∇λ, and ∇λ is a field concentrated on S with zero
curl and unit circulation around the conductor. Now discretise ∇λ directly: on the mesh, the cocycle z
is exactly the edge-wise version of "λ increases by z_e along edge e", and its Whitney interpolant

    ψ = Σ_e z_e w_e        (w_e the Nédélec edge function of edge e, ∫_{e'} w_e·dl = δ_{ee'})

is the discrete ∇λ. It has two properties that need no geometry at all:

- curl ψ = 0 in every element. The curl of a Whitney 1-form is the Whitney 2-form of the coboundary,
  curl(Σ z_e w_e) = W₂(δz), and δz = 0 face by face is the cocycle condition the core already checks.
- ∮_γ ψ·dl = Σ_{e∈γ} z_e = ⟨z, γ⟩ for every mesh loop γ: 1 for the loop the generator was paired with,
  0 for the others (Pellikka et al. 2013 eqs. 4.5, 5.3). The value of I is the pairing, not the shape
  of the support (Grok Q3, high).

So h = −∇φ + Σ_k I_k ψ_k with φ single-valued and one scalar I_k per generator (Pellikka eq. 4.5;
Alves et al. 2022a eq. 2 writes the same h = −∇φ + Σ I_i ψ_i). The thin cut is ψ pushed onto a
surface and condensed into duplicates; ψ itself needs no surface, no duplicates, no consistent
sides. A wandering or conductor-hugging z merely makes ψ's support larger; it changes nothing in the
imposed current. **This is the whole reason option 2 is immune to the corc failure.**

### 3.3 Element level (TET4 in air)

For an air tet K crossed by generator k (at least one edge with z_e ≠ 0; today's `CutData::elements()`,
`cl_CutData.cpp:263-362`, is that list), the element dofs become [φ_1..φ_4, I_k] and

    h|_K = B φ + ψ_k I_k,      ψ_k|_K = E(x) z^k_K      (E: 3×6 Nédélec matrix, z^k_K: the 6 edge coefficients)

Today `maxwell::phi_tet4` assembles M = μ₀ BᵀB dV (`mt_maxwell_phi.cpp:150-169`; general quadrature
:26-46). With the augmented operator B̃ = [B | E z] the same product gives

    M_φφ = ∫ μ BᵀB          (unchanged)
    M_φI = ∫ μ Bᵀ E z       (4×1 per generator)
    M_II = zᵀ (∫ μ EᵀE) z    (1×1)

The Gram is even in the sign of B, which is why BELFEM assembles the φ block without the minus of
h = −∇φ (Grok Q1); the only sign that matters is that z is read on the same mesh-oriented edges the
interface relation uses (§3.4), with no extra flip. That is an orientation test, not a derivation (O1).
Magnetodynamics adds nothing: air has no curl-curl term and the same Gram multiplies ∂t(μh).

### 3.4 Interface edges (conductors and thin shells)

A conductor or thin-shell edge dof on the air interface hangs on the two air nodes with weights ±1,
h_e = φ_{n0} − φ_{n1} (`create_hanging_edges_and_facets`, T-matrix `cl_FEM_DofMgr_DofData.cpp:3630-3647`).
With ψ the tangential trace of the air field on an interface edge in the support is

    h_e = φ_{n0} − φ_{n1} + I_k z^k_e

i.e. the edge dof gets the abstract node as a **third hanging source with weight z_e**. The machinery
exists: `CutSet::create_duplicates` already makes the abstract node a source with weight 1
(`cl_CutSet.cpp:42-82`). Interior conductor edges stay free. Thin-shell edge dofs hang the same way
(`cl_MaxwellFactory.cpp:1742-1757`), so they are covered by the same rule. Open: a support that
reaches a tape from both sides with opposite signs would inject 2I onto the sheet (Grok Q2, O2).

### 3.5 The one missing piece: a global dof inside a volume element

Today I reaches an element only as a hanging source of a duplicated node; no mechanism attaches a
carrier-free global dof to selected volume elements (Codex Q2: circuit unknowns live in the circuit
state, bearings reference an existing node, facet lambda dofs are attached to a `Facet`). The nearest
precedent is `Element::link_lambda_dofs` (`cl_FEM_Element.cpp:1292-1310`), which appends a facet's
lambda dof to the element dof list. Option 2 needs the analogue for volume elements:

- a per-generator "crossed element" list (from `CutData::elements()`), persisted with the cochain;
- `DofData::num_dofs_per_element` (`cl_FEM_DofMgr_DofData.cpp:3447-3463`) is per block; crossed
  elements need +1 per generator, so either the air block's element linking appends generator dofs
  per element (like lambda dofs) or all air elements of a block carry the slot with a zero column;
- the IWG evaluates ψ_k at the integration points from the element's six z coefficients (the
  Calculator already provides E for Nédélec blocks);
- MPI: the distributor must send the cochain to every rank owning a crossed element
  (`cl_Mesh_Distributor.cpp:259-280`, :363-421 today handle abstract nodes and duplicates).

### 3.6 Everything else

- Pins: φ is single-valued, ker ∇ = constants, one pin per air component as today.
- Current condition: I_k is the same abstract-node dof, prescribed (current) or free (voltage);
  nothing changes in `cl_MaxwellBoundaryConditionFactory` except that no `DomainType::Cut` sideset
  exists (its consumer table in `cl_Maxwell_FieldList.cpp:383-388` is dead already, Codex Q1).
- Postprocessing: φ plotted continuous; H = −Bφ + Σ I_k ψ_k on crossed air elements
  (`MaxwellPostprocessor::compute`, `cl_MaxwellPostprocessor.cpp:453-479`) and in `compute_hn`'s
  φ-side trace (`cl_FEM_Calculator.hpp:2626-2743`); B = μH.
- Persistence: the `.bfm` stores abstract nodes, duplicates and hanging sources
  (`cl_Mesh_BfmFile.cpp:328-425`, :1287-1532) but no cohomology record; option 2 needs per generator
  the edge ids, signs, periodic mates (`CutData::collect_edges` :186-228 synthesises the mate with
  the same sign) and the crossed-element list.
- Cost: the I_k row/column couples to every node of the band (corc: ~7.5k of 190k nodes) — the same
  structure as the free-cut lambda dense row BELFEM already handles with direct solvers; a tight z
  (option 1) shrinks it to a sheet's worth.

## 4. Option 1 in more detail

Minimise Σ_e w_e |z_e + g_{head(e)} − g_{tail(e)}| over integer g. The constraint matrix is the graph
incidence matrix, totally unimodular, so the LP relaxation is integral for any positive weights
(classical; Dey, Hirani & Krishnamoorthy 2011 is the same statement for general complexes without
torsion, on the chain side). With w_e = area of the dual face this is the minimum-area representative;
with w_e = 1 the minimum edge count. A closed sheath around the conductor is a coboundary with pairing
0, so it is never the minimiser of the class (refutes Grok's "tight cylinder" counterexample); but
minimality does not by itself prove the positive-case faces close (medium), so the 2026-09-02 guard
stays. Existing pieces (Grok Q4): SPFA solves the L∞ feasibility |c| ≤ 1 (`cl_Cohomology.cpp:445-456`),
greedy sweeps are coordinate descent that fire only at |c| ≥ 2 (:548), Tier-A pocket removal is one
cheap coboundary restricted to structure-free pockets (:668-682). Missing: the global weighted-L1
solve (min-cost flow on the dual, or max-flow between the two lifts of the current support),
area-like weights, and a closedness certificate after it. Lives outside the core as a post-pass.

## 5. Ordered steps (option 2)

- [ ] **R1 — orientation spike (no new code paths).** On tapestack3d (healthy thin cut) read the
      stored z, compute ⟨z, γ⟩ against the imposed I and ∮H from the save, and fix the sign convention
      of §3.3/§3.4 (O1). Gate for everything after.
- [ ] **R2 — persist the generator.** Per generator: edge ids + signs, periodic mates, crossed
      elements; `.bfm` record; distributor transport (§3.5, §3.6).
- [ ] **R3 — element generator dofs.** Append I_k to crossed air elements' dof lists (precedent
      `link_lambda_dofs`); sizing in `num_dofs_per_element`; Calculator access to z_K and E.
- [ ] **R4 — augmented φ kernel.** `phi_tet4` / `phi` with B̃ = [B | E z]; opt-in per generator so
      tight generators may keep the thin cut (O3).
- [ ] **R5 — interface relation.** Third hanging source (abstract node, weight z_e) on interface
      edges in the support; conductor and thin-shell edges (after: R3).
- [ ] **R6 — postprocessing.** H and h_n with the ψ term; φ continuous (after: R4).
- [ ] **R7 — gates.** tapestack3d: option 2 vs thin cut, ∮H = I at every station and |ΔH| within
      solver tolerance; corc_solder_control: winding = I at r = 5 mm along the cell, tape J not
      rim-only; periodic corc; `make check-fast`.

## 6. Open design questions

- **O1** Sign: mesh-oriented z_e with h_e = φ_{n0} − φ_{n1} — one convention or an extra flip? Settle in R1.
- **O2** Thin shells reached by the support from both sides (opposite signs inject 2I on the sheet):
  detect and refuse, or restrict z to one side by construction?
- **O3** Hybrid: keep the thin cut for generators whose positive faces close (no pruning loss) and use
  ψ only for the others, or switch all generators to ψ once R7 passes?
- **O4** Conditioning of M_II on a 14k-edge band; measure on corc before deciding whether option 1
  must precede option 2 in production.
- **O5** MPI: per-element generator dofs on aura elements; who owns I_k.
- **O6** Element coverage: option 2 would need a Nédélec entry for PYRA5 (`cl_FEM_Calculator.cpp:1100-1140`)
  to cover every shape the cohomology engine already handles. Separately, and independent of this plan:
  should the thin-cut path get an explicit "unsupported element type" error, since its type is global and
  a non-tet air element currently reaches the tet bitset routines by accident (§2.1b)?

## 7. Definition of done

- [ ] tapestack3d reproduces the thin-cut solution to solver tolerance with ψ (R7).
- [ ] corc_solder_control winds by I at every station without any thin cut.
- [ ] `.bfm` round trip carries the generator; restart identical.
- [ ] A generator with |z_e| ≥ 2 assembles and imposes the right current (the case the thin cut cannot
      represent at all, §2.1a) — the sharpest single test that option 2 is doing something new.
- [ ] Docs: `src/homology/doc/thick_thin_cuts_and_conjugate_edges.md` gains the ψ section;
      `cohomology_theory_and_implementation.md` §FEM Discretization updated.

## 8. Audit trail

- `tmp/ai_exchange/corc_thin_cut_collapse.md` — the diagnosis round (Codex sol/high, Grok high).
- `tmp/ai_exchange/corc_thin_cut_guard.md` — the guard diff round.
- `tmp/ai_exchange/cut_representative_options.md` — this plan's research round: Codex terra/high
  (code inventory, all cites re-opened), Grok high (formulation; its cylinder counterexample for
  option 1 refuted, its Dey-scope correction taken, its option-3-equals-option-2 argument taken).
- Gregory Giard, 2026-09-03 (§2.1, and the standing decision): the coefficient and element-shape
  motivations, and the ruling that the thick cut is the eventual target but not a current priority.
- Literature: Pellikka et al. 2013 §4–6 (`literature/papers/topology/pellikka2013.txt:369-378`,
  455-462); Alves et al. 2022a eq. 2; Dey et al. 2011; Henrotte & Hameyer 2003 and Kettunen et al.
  1998 are Pellikka's refs [21], [24] for the cut-free lineage — not in the library, not page-verified.
