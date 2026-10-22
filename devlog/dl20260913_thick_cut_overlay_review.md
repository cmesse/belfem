# Devlog 2026-09-13 — Thick-cut overlay plan review

**Date:** 2026-09-13
**Topic:** Independent technical review and comments on `todo/thick_cut_overlay_block.md`
**AIs involved:** Codex
**Codex Audit Confidence:** high on the source traces and algebra; medium on the unresolved buffer extension
**Verification:** static source review on `thickcuts` @ `f81dd8432a1a47e9644c48a7edc7a2c781219ae3`, plus a standalone `python3` integer/rational arithmetic probe described below. No BELFEM build, test executable, MPI run, or solver run. The implementation is not verified.

## Verdict and current scope

The overlay is a credible implementation carrier, and the plan identifies much of the required
mesh, dof, interface, persistence, and recovery work. **Revise R2 and the acceptance gates before
using it as an implementation checklist.** The most consequential gaps are the early thin-only
classifier, the cleanup contract of the core failure exit, and the orientation of coefficients
after edge reconstruction.

Christian clarified the scope during this review: **try thin cuts first; if that fails, use thick
cuts for the entire model. There will be no mixed thin/thick realization.** Thus:

- All generators use thin cuts if the attempt succeeds; all use thick cuts if any generator
  prevents that attempt from completing.
- Per-generator failure information is useful for diagnostics. It does not select individual
  generators for thick assembly.
- Multiple thick generators may still overlap. The full current-current Gram block remains
  necessary.
- The original hybrid skip logic and its effort estimate are superseded. Test access to the
  thick builder can remain useful, but a production `thin | thick | auto` menu is not required
  by this clarified strategy.

No implementation authorization was given or assumed. Proposals involving the closed cohomology
core remain for Gregory/Christian to review and apply under protocol §7.1.

## Findings and comments

### F1 — Decide globally before committing thin-cut topology

**Priority: P1 for the implementation plan. Confidence: high. Applies to R2/R3, O3/O9.**

The current constructor already offers a useful separation: it calls `collect_facets()` before
`collect_nodes()`, `determine_cut_sets()`, node duplication, element relinking, and thin-sideset
creation (`src/homology/cl_CutProcessor.cpp:71-92`). This supports the clarified strategy without
first creating thin cuts and then attempting to undo them.

Make the decision after the thin feasibility attempt and before those topology mutations. On
failure, discard the tentative thin-cut data and build overlays for **every** final generator;
create the abstract-node registry once, preserving its existing ordering. Do not retain a
successful earlier generator's cut when a later generator fails. Tentative mesh flags and
indices also need a defined cleanup path: collection already rewrites edge indices
(`cl_CutProcessor.cpp:163-187`) and element indices (`cl_CutData.cpp:315-321`).

The plan's placement is currently inconsistent: R2 calls the thick builder after constructing
`CutProcessor`, whereas R3 asks for registration before that processor's final `finalize()`.
The constructor finalizes before returning (`cl_CutProcessor.cpp:92`;
`cl_CutFactory.cpp:460-473`). Name the actual construction phase and its finalization obligations.
Also require `copy->node(i) == parent->node(i)` after all relevant relinking, before partitioning
and serialization. `InterfaceProcessor` relinks elements selected through physical master-block
IDs, so a new overlay block does not automatically follow a parent relink
(`src/homology/cl_InterfaceProcessor.cpp:95-141,231-246`). This is an invariant to enforce, not a
claim that every supported air parent is currently relinked.

**Gate:** at least two generators, with a successful early thin candidate and a later failure;
the result has all required thick generators, no committed thin cuts, no thin-attempt duplicate
nodes, and stable abstract-node/current ordering. Interface duplicates created independently
of cuts remain legitimate.

### F2 — Non-unit cocycles must bypass the thin-only local classifier

**Priority: P1. Confidence: high. Applies to R2, O13.**

Replacing the non-tight errors in `collect_facets()` is insufficient. Before reaching them,
`CutData::collect_elements()` invokes `determine_cut_case_2d/3d()` on each candidate
(`src/homology/cl_CutData.cpp:333-354`). The 3D classifier accepts seven unit-cut support patterns
and aborts otherwise (`:480-608`); it then checks coefficients obtained through the ±1 bitsets
(`:612-642`). Those bitsets omit larger coefficients (`:231-258`).

A valid closed integer cochain can fail this classifier. On a TET4, choose local nodal values
`p=(0,1,2,3)` and edge coefficients `z_ab=p_b-p_a`. In BELFEM's local edge order the result is
`[1,1,-2,3,2,1]`: every face circulation is zero, but all six edges are nonzero, giving pattern
63, which the thin classifier rejects. Local exactness is expected on a tetrahedron; this
example is not presented as a global nontrivial cohomology generator.

Separate support/coefficient extraction from thin-case classification. Once the model selects
the thick path, extract every final generator's full integer record without using the thin
classifier. Keep its existing errors for the unit thin path. For the thick path, validate the
actual contract—closedness, orientation, and the intended cycle/terminal pairing—instead of
using admissibility as a thin sheet as the validity test.

**Gate:** this closed six-edge example is accepted by the thick record/kernel path; a deliberately
nonclosed cochain is rejected. The production fallback must reach the builder before any
thin-only classification of an unrectified generator.

### F3 — The owner-applied core exit needs an explicit lifecycle contract

**Priority: P1. Confidence: high. Applies to R2, O13.**

The infeasible branch is more than its final error message. It unflags candidate elements and
edges, flags obstruction edges, changes element levels, and writes `error.exo`
(`src/homology/cl_Cohomology.cpp:496-532`). The comment at `:496-497` explicitly relies on the
process aborting. A return inserted at `:534` would also bypass the normal mesh-flag reset and
progress cleanup (`:629-641`).

Propose a recoverable branch **before** these destructive diagnostics. Define cleanup for flags,
progress state, and temporary data on both exits. With model-wide fallback, ending rectification
for the whole basis is an available policy; there is no requirement to continue finding thin
representatives individually. That choice and any change to pocket removal are the owner's
decision, not an AI patch inferred from the abort site.

There are two cleaning stages: constructor cleaning, then terminal alignment followed by another
`clean()` (`src/homology/cl_CutFactory.cpp:344,360-366`). Terminal alignment creates new linear
combinations, deletes the old generators, and replaces the container
(`src/homology/cl_Cohomology.cpp:1189-1208`). Therefore the overlay must snapshot the **final,
terminal-aligned** generators. An aggregate fallback decision can survive this transformation;
old generator pointers or index-associated per-generator verdicts cannot simply be treated as
the final basis. The enabled policy must be available during construction—a later setter cannot
rescue an abort inside the constructor.

Calling this a known ten-line or trivial change is premature until that API and cleanup contract
are concrete. This does not dispute that non-unit coefficients are a supported mathematical
input to a thick-cut formulation.

### F4 — Canonical storage orientation and live edge orientation are different contracts

**Priority: P1. Confidence: high on the contract gap; no live overlay exists to reproduce it yet.
Applies to R1/R4/R6/R7, O1/O8, §6.**

The plan says `Mesh::compute_edge_directions()` orders edges by original node ID. It actually
compares the element-local endpoint sequence with the **existing edge's** endpoint sequence,
using `original()->id()` to recognize matches and reversals
(`src/mesh/cl_Mesh.cpp:1844-1860`). `EdgeFactory` constructs that existing sequence from increasing
current node indices (`src/mesh/cl_EdgeFactory.cpp:319-336,477-490`). These are not a general
guarantee of increasing original IDs.

Persisting coefficients against an ordered original-ID pair is reasonable. On reconstruction,
explicitly convert from that canonical orientation to the live edge orientation consumed by
`EF_TET4`, which already reads `edge_directions()`
(`src/fem/interpolation/nedelec/cl_EF_TET4.cpp:94-95`). For example, canonical `z=2` on `10→20`
becomes `z=-2` if the rebuilt edge runs `20→10`. Converting canonical-to-live orientation and
letting E handle local-to-live orientation are separate operations, not a duplicate sign flip.
Likewise, expand periodic mates according to their endpoint correspondence; the current thin
path's same-sign copy is a convention to check, not a universal rule about arbitrary ID pairs
(`src/homology/cl_CutData.cpp:249-258`).

**Gate:** reverse endpoint order and permute node storage, including periodic/interface
duplicates and a `.bfm` reload. The reconstructed vector field and signed circulation must stay
unchanged. One healthy deck's overall current sign is insufficient coverage for this contract.

### F5 — The proposed fixture and several acceptance gates do not exercise the promised fallback

**Priority: P1 for fallback coverage. Confidence: high. Applies to R9 and the core-exit test in R2.**

`AnnulusForcedRepair` deliberately adds an exact node coboundary to an already-unit generator,
then expects `clean()` to restore unit coefficients with the pairing preserved
(`tests/homology/test_Cohomology.cpp:255-260,301-340`). It is a **repairable non-unit representative**,
not a class with no unit representative. The fixture uses the planar annulus builder at `:265`.

Keep it for integer-coefficient preservation or direct builder tests, but add a distinct
infeasible difference-constraint fixture to test the core status exit. A local coefficient of
magnitude 2 alone proves neither infeasibility nor that a thin cut cannot represent the same
cohomology class. If the existing cleaner is invoked first, this fixture can become unit again
before the overlay ever sees it.

R9 also primarily compares forced `thick` with `thin`. Under the clarified production strategy,
test the automatic dispatcher itself:

- Healthy model: all thin, no overlays.
- Non-tight failure and SPFA-certified infeasible unit-representation failure: all thick, with
  explicit evidence that the fallback ran.
- Multiple generators: late failure abandons the entire thin attempt; overlapping thick
  generators retain off-diagonal current coupling.
- Invalid cocycle: rejected rather than silently treated as an ordinary fallback.
- Equivalent dispatch result and coefficient identity after reload and under MPI, including
  empty local overlay sets; both prescribed and free current modes.

Direct TET4 builder/kernel tests can exercise the thick implementation independently of this
dispatcher. A test-only forcing hook need not become a production input option.

### F6 — Shell coverage needs a decision for buffer-generated traces

**Priority: P2, with potential formulation impact. Confidence: high on the omitted producer;
medium (~75%) on what extension particular buffer cases require. Applies to R7/R8, O6/O10.**

R7 lists conductor interfaces and the two outer shell-edge hanging helpers. Another producer,
`ThinShellFactory::create_buffers()`, turns nonconducting layers into `DomainType::Buffer` and
hangs their edges on their own nodes (`src/fem/kernel/cl_ThinShellFactory.cpp:2049-2103`). It runs
after layer construction (`:247`). Air and Buffer both dispatch to the scalar-potential kernels
(`src/fem/maxwell/cl_IWG_Maxwell.cpp:286-310`). The nominated `tapestack3d` deck includes a magnesia
buffer (`examples/tapestack3d/input.conf:118-129`).

The plan does not state how a cocycle is extended to these later-created phi layers, or establish
that it vanishes there. Resolve that question before claiming the two outer hanging helpers
complete shell support. Depending on the answer, buffer traces and their volume/recovery terms
need extension, or the relevant case needs an explicit scope restriction and a different first
gate. This review does **not** prescribe that extension or claim every buffer has nonzero support.

### F7 — Correct the support-threshold rationale and the sign inconsistency

**Priority: P2. Confidence: high. Applies to gap row 3, R2, §2 and R6.**

For a closed cochain on a full TET4, the restriction is a nodal difference: `z_ab=p_b-p_a`.
Nonconstant values on four vertices produce at least three nonzero edges; on a TRI3 the minimum
is two. This holds for arbitrary coefficients, not only ±1. Thus the `>= dimension` threshold
at `src/homology/cl_CutData.cpp:294-304` does not lose the support of a valid closed TET4 cochain.
Enumerating nodal values in `{-2,-1,0,1,2}` gave possible nonzero support counts `{3,4,5,6}` for
TET4 and `{2,3}` for TRI3. The general argument is the minimum number of edges between distinct
vertex-value classes in a complete simplex graph.

Using `>=1` as a generic support collector is fine, but describe it as decoupling extraction from
thin-cut assumptions, not as fixing a demonstrated non-unit support loss in v1. Closedness
validation is still necessary.

For the plan's `h=-B phi+Psi I`, the coupling is `-B^T Psi`, as correctly written in R6. The
positive sign in §2 is inconsistent with that convention. The existing phi kernel forms
`mu0 * integral(B^T B)` (`src/fem/maxwell/matrices/mt_maxwell_phi.cpp:34-45`); the reference-tetrahedron
probe below independently checks the cross-term sign.

## Parts of the plan to retain

The basis enrichment is supported by Alves et al. 2022b, Eq. 10, and Pellikka et al. 2013,
Eq. 4.5. Both associate the additional coefficients with cycle integrals through the chosen
dual basis; integer coefficients need not be restricted to ±1.
[Alves et al., author manuscript](https://arxiv.org/pdf/2108.08828),
[Pellikka et al., author preprint](https://www.gmsh.info/doc/preprints/gmsh_homology_preprint.pdf).
The following are implementation conclusions from the source review and algebra, not a
literature endorsement of BELFEM's proposed block carrier:

- Retaining the physical air elements and assembling zero phi-phi contribution on the overlay
  avoids duplicate air energy. Both coupling row and column, and the full overlapping-generator
  Gram matrix, are needed. All-or-nothing selection does not remove overlap.
- Reusing abstract-node dofs fits current assignment. Their ordering is an explicit factory
  contract (`src/fem/maxwell/cl_IWG_Maxwell.cpp:100-119`). Cell dofs are the wrong semantic carrier
  for a current shared across a band.
- Parent-owner inheritance and abstract-dof visibility are separate requirements. The existing
  node-dof proc flags come from geometric element nodes
  (`src/fem/kernel/cl_FEM_DofMgr_DofData.cpp:1173-1216`); shipping abstract mesh nodes alone is insufficient.
- R6's M-channel choice is credible: `Calculator::qold()` and `q()` traverse the actual local dof
  list (`src/fem/kernel/cl_FEM_Calculator.cpp:2894-2915,2995-3000`). Still test fixed-current
  history and free-current transients after adding the slot. Explicitly exclude the new type
  from unwanted `MaxwellData` construction at both `allocate()` and `link_maxwell()`
  (`:1233-1247,1740-1754`); avoiding helper calls inside the new kernel alone does not do that.
- The hidden-block work is real: Exodus selects nonempty blocks without testing hidden status
  (`src/mesh/cl_Mesh_ExodusWriter.cpp:127-142`), and `ProtoMesh` reconstructs a block's type and
  thickness without restoring the hidden bit (`src/mesh/cl_ProtoMesh.cpp:1333-1345`).
- Recovery must run on physical elements and include the correction in the solve-side shell
  trace (`src/fem/kernel/cl_FEM_Calculator.hpp:2698-2701,2742-2746`). Preserve the separate
  prescribed-current and free-current gates.

## Remaining planning comments

R0 remains approval work; this report does not tick it. O3 is resolved by Christian's scope
clarification. Reconcile O2 with mandatory thin-first behavior and distinguish any diagnostic
override from that production policy. Resolve dimensional scope: the guards say TET4/3D, while
the failure table calls the 2D case v1 and O16 leaves it open.

R4 should spell out malformed-record checks (sizes, signed coefficients, missing parents,
generator bounds), the `ProtoMesh` reconstruction path, and how records are available on ranks
that need physical-element recovery or interface sources. Current-owner equality is necessary;
it does not by itself specify all aura/recovery consumers. Use chunked `share`/`receive` or the
appropriate existing wrapper, with matching message order even for empty local sets.

The R6 inert-kernel probe needs a precise observable. Prescribed I remains fixed even if the
volume kernel is disabled, and an intact interface current term may still impose circulation.
Measure the missing volume energy/residual contribution for a kernel-only deletion; distinguish
that from deleting the cocycle correction everywhere.

The 11–14 working-day estimate is historical, not confirmed by this audit. Removing mixed-mode
work should simplify R2, while the core exit, genuine infeasibility fixture, and buffer scope
still need definition. Re-estimate after those are settled; devlog entry counts are not a
measurement of this plan's critical path. R10 remains after implementation and executable gates.

For the protocol's implementation checklist: naming and introduced OpenMP/thread behavior are
N/A at this plan-only stage. Preserve explicit ownership: `Block::~Block()` deletes its elements
(`src/mesh/cl_Block.cpp:34-38`), so parent pointers cannot serve as owned overlay entries. The
record's `int[n][n_edges][|G|]` notation is a layout sketch, not a C++17 variable-length array
declaration. Use owned `Cell` storage for transport records and `Matrix` for the actual algebra,
with setup-time size checks and preallocated kernel scratch. Malformed persisted data needs
release-active errors; expected thin infeasibility needs a status. No C++ implementation,
allocation performance, or MPI execution is signed off by this review.

## Arithmetic probe and limits

Executed with `python3` via stdin on the checkout named in the header, using only `itertools`
and `fractions`. It did not call BELFEM or modify a cohomology implementation. Inputs were the
TET4 edge order from `src/mesh/cl_Element_TET4.hpp:100-140`, nodal values `(0,1,2,3)`, and the
reference tetrahedron with vertices `(0,0,0),(1,0,0),(0,1,0),(0,0,1)`.

| Check | Result |
|---|---|
| Closed local integer cochain | `[1,1,-2,3,2,1]`, zero face circulation, support pattern 63 |
| Nonzero supports for integer nodal values -2 through 2 | TRI3: 2 or 3; TET4: 3, 4, 5, or 6 |
| Whitney interpolation of that nodal difference | `E z = (1,2,3)` at three rational barycentric points |
| Second field from nodal values `(0,2,1,-1)` | `(2,1,-1)` |
| Exact cross-generator mass divided by mu0 | `1/6`, so independent per-generator blocks omit a nonzero term |
| Node-1/current mass entry divided by mu0, zero-based node numbering | `-1/6`, matching R6's sign |
| Canonical coefficient 2 on 10→20, live edge 20→10 | live coefficient -2 |

The constant fields make the two integrals exact: multiply their dot products by the reference
volume `1/6`. These are algebraic discriminators, not a fallback regression. No claim is made
that the earlier thin/thick field comparisons, MPI gates, nonfree model, or restart gates ran.

## Files updated and follow-up

- This report and `devlog/README.md`.
- `todo/thick_cut_overlay_block.md`: current all-or-nothing scope and an unchecked review-comment
  checklist; implementation work remains pending.
- `todo/README.md`: link to this review beside the existing plan.
- `tmp/ai_exchange/thick_cut_plan.md`: short audit handoff; this report is the durable record.

The follow-ups remain attached to the existing unapproved feature plan. This review creates no
new claim of a reachable production-code defect and closes no existing debt-register item.
