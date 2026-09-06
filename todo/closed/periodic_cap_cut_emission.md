# One-Sided Cut Emission on the Periodic Identification (match_edges mismatch)

**Date:** 2026-07-16
**Purpose:** Fix the solve-phase periodic rebuild abort — "number of edges on source
and target periodic surfaces do not match" — root-caused (three-AI, 2026-07-16) to
one-sided emission of thin-cut geometry on the periodic identification: cut faces on
or crossing a cap duplicate nodes/edges on one geometric copy only, so the rebuilt
raw edge objects differ between the two cap collections while node counts fold equal
through `original()`.
**Module:** `src/homology` (CutProcessor / CutData / CutSet), `src/mesh`
(PeriodicityFactory as detector — the assert is CORRECT and must stay)
**AIs involved:** Claude + Codex + Grok (diagnosis); fix pending Christian's direction
**Status:** ✅ FUNCTIONALLY COMPLETE 2026-07-16 (final Codex audit + Christian's field
validation pending). The full periodic seam stack now passes on the double-layer corc:
XOR-once quotient self-cancel (149=149 validated) → facet-mediated pairing →
combo-III direct ties (Grok EDGE-TIES: Whitney I-cancellation; jump stays in node
hanging) → half-cut policy (untied+censused, DOF-inert in air) → symmetric trace-twin
aliases (1:2 realizations along cut boundary curves) → E1 compaction. Final census:
tied 8538 (4306 through duplicates), halfcut 1630, alias 1380, facetMismatch 0,
uncovered 0+0; set_entity_dependencies passes with original-aware orientation. The
run now reaches FIRST JACOBIAN ASSEMBLY and dies at the KNOWN D1 defect of
`todo/maxwell_kernel_collapse_plan.md` (null-material MaxwellData on air blocks) —
next blocker belongs to that plan. Open here: half-cut closure verification IF a
conductor-side case appears (census reports it); audit verdict.

> **Closed 2026-08-09** (todo currentness sweep). Re-verified in tree: the
> facet-mediated `match_edges` rewrite with combo-III ties, the half-cut policy and
> the trace-twin alias asserts are live in `cl_Mesh_PeriodicityFactory.cpp`
> (`:872`, `:928-1148`) and committed; the seam probes that produced the census
> numbers have since been stripped. The D1 blocker named above is also fixed
> (`maxwell_kernel_collapse_plan.md`). The one residual — half-cut closure IF a
> conductor-side case ever appears — is a **watch item**, not work in flight: the
> census assert guards it and will fire loudly if it happens. Kept here as the
> mechanism/formulation record for the seam.

> **Scope guards:**
> - The `match_edges` assert (`cl_Mesh_PeriodicityFactory.cpp:863-868`) is telling the
>   truth — never weaken it or fold edges by original to silence it (that would leave
>   the duplicate's DOFs unconstrained across the seam: helix flipped-J failure mode).
> - Representative choice (sparsity, pockets) modulates the mismatch magnitude but
>   cannot fix it — proven by the 2026-07-16 experiment below. Keep those concerns in
>   their own plans.

## Evidence

1. **Signature:** node counts equal (collect_nodes folds `original()`,
   `cl_Mesh_PeriodicityFactory.cpp:576`), edge counts differ (collect_edges counts raw
   `edge->index()` objects, `:615`, `:628`).
2. **Representative dependence (the controlled experiment):** dense SPFA
   representatives (~4'100 in-plane cap edges/gen) → 12458 vs 12498 (Δ = +40,
   slave-heavy). Sparse greedy representatives (~500–1'500 in-plane) → 10877 vs 10728
   (Δ = −149, master-heavy). The asymmetry scales and flips sign with cap-touching
   cut mass — it is emission machinery, not a fixed mesh feature.
3. **Mechanism chain (Codex, high confidence, file:line):**
   - periodic cap faces survive as thin-cut faces — the boundary-face unflag pass
     covers only non-periodic `mPhiBoundaries` (`cl_CutProcessor.cpp:420-423`,
     `:472-479`);
   - slave/target seam faces are emitted one-sided from `tFace->slave()` after
     `fix_face_slaves()` cleared their master (`cl_CutData.cpp:159-172`;
     `cl_Mesh_PeriodicityFactory.cpp:1266-1307`);
   - one-bit periodic node pairs get UNPAIRED duplicates
     (`CutSet::create_duplicates`, `cl_CutSet.cpp:144-155`);
   - `EdgeFactory::edge_key()` keys on current (not original) node indices
     (`cl_EdgeFactory.cpp:324-337`), so post-relink edge recreation makes distinct
     raw objects on the duplicated side only.
4. **Topology (Grok):** closedness forbids a free cap trace; admissible seam
   configurations are transversal wrapping (historically handled — probe-4b's 76%
   legitimate one-sided pairs), in-plane bands (the suspect class; probe-4b
   sign-coherence history), tangency, quotient pockets (removable, none found on the
   corc caps), and essential crossings.
5. **Known non-causes:** removable pockets (census: zero Tier-A candidates; cap-pure
   boundary sums 466–1099 ≠ observed deltas); the stale `mPairVerdict`/Step-6c
   comments (superseded by Step-9 original-identity keying — do not chase).

## Probe result (2026-07-16, `probe_duplicate_symmetry`)

Christian's matching-tier analysis prompted a discriminating probe (per-original
duplicate multiplicities through the periodic tie; per-original-key raw edge counts):

- **Pre-cut update: perfectly symmetric** (raw = folded = 3197 both sides; 9318 =
  9318 edges; zero mismatches). The asymmetry is entirely cut-created.
- **Solve-phase: 1070 originals with mismatched multiplicities, in BOTH directions**
  (2v3, 1v2, 3v1 …; master 1753 vs slave 1592 duplicates), all on one cap plane in
  the tape-adjacent annulus. Edges: net 149 but sumAbsDiff **1491** — scattered ±1
  asymmetries in both directions, not a single trace.

**Revised diagnosis:** bidirectional per-node differences are what LEGITIMATE
one-sided cut membership produces (probe-4b / dl20260616: period-wrapping
seam→interior one-sided pairs are legitimate, and the NODE periodic machinery was
explicitly extended for them — ONE_SIDED_JUMP, Step-9 keying). "Equal duplicate
counts per original" is therefore NOT a physics invariant. The actual gap: the node
pipeline is duplicate-aware, the EDGE rebuild is not — `match_edges` keys by
original endpoints (collapsing duplicates many-to-one, silently) and asserts strict
raw-count equality, a cut-free-era invariant.

## Candidate fix directions (decision: Christian)

- [x] **(d) — CHOSEN AND IMPLEMENTED (probe-derived): duplicate-aware edge rebuild.** Induce the
  edge pairing from the raw node ties: pair edge (a, b) with (P(a), P(b)) where P is
  the node-level periodic tie including duplicate ties and the one-sided policy;
  edges whose raw endpoints have no tied partner follow their nodes' one-sided
  handling. Replace the count-equality assert with a tie-consistency assert.
  Christian's tier-2/3 ordering scheme (originals first; canonical duplicate
  ordering by cut side, facet index as tiebreaker) lives inside P's construction.
  (a)/(c) below reduce the surface but cannot remove it: period-wrapped conductors
  REQUIRE transversal cap crossings, which duplicate cap nodes regardless.

- [ ] ~~**(a) Keep cuts off the identification** — NOT TAKEN, superseded by (d):~~ treat periodic cap sidesets like
  non-periodic φ boundaries in the cut-face selection (extend the unflag pass at
  `cl_CutProcessor.cpp:472-479`), so in-plane cut faces are pruned and the cut routes
  through the interior. Cheapest if the face-realization freedom allows it everywhere
  (the self-canceling-face removal suggests it does); transversal crossings remain
  and are handled.
- [ ] ~~**(b) Symmetric emission** — NOT TAKEN, superseded by (d):~~ emit in-plane cut geometry on BOTH copies with
  periodic ties (extend `create_duplicates` pairing to the one-bit case with proper
  backups). Physically faithful in the quotient; touches the historically delicate
  duplicate/backup machinery.
- [ ] ~~**(c) Representative-level avoidance** — NOT TAKEN as a fix; survives only as the Phase-4 sparsity goal in `thin_cut_nonunit_rectification_implementation.md`:~~ weighted L1
  min-cost representative with heavy weights on cap edges — cuts avoid the caps
  except transversal crossings. Principled, also delivers guaranteed sparsity;
  largest implementation (new min-cost-circulation solver; complexity table in the
  exchange thread).

~~Recommendation (Claude, medium confidence): (a) first as the minimal intervention if
a probe confirms the in-plane faces are avoidable on corc; (c) as the durable
end-state alongside Phase 4; (b) only if (a) proves topologically impossible for
some class (a cut forced to live in the cap).~~
**Outcome (2026-07-16): the recommendation was overtaken by the probe.** Christian
took (d) — the duplicate-aware edge rebuild — and it is what shipped: the
facet-mediated `match_edges` rewrite with combo-III direct ties, the half-cut
policy, parity-aware trace-twin aliases and E1 compaction
(`cl_Mesh_PeriodicityFactory.cpp`, committed with the seam campaign in `ce6a0e8b`
and neighbors), plus the XOR-once quotient self-cancel in `cl_CutProcessor.cpp`.

## FORMULATION RESOLVED (2026-07-16, Grok F1–F3 + Codex V1–V3, Christian's hypothesis)

**Christian's sign-asymmetry hypothesis: CONFIRMED and sharpened.** The thin cut is
emitted by positive-case elements only (`cl_CutProcessor.cpp:433-435`; doc §
"positive-case elements only"); interior coordination is Face-OBJECT identity
(self-cancel `:461-469`, jump-vs-neighbor across the shared facet), and nothing in
the pipeline consults `Face::periodic()` — so quotient-neighbors across the
identification decide independently. Sharpening (Grok, medium-high): the cut-case
SIGN is computed in the LOCAL edge-orientation frame (`cl_CutData.cpp:615-622` via
`edge_direction`), so a period-consistent cochain can legitimately evaluate `+k` on
element A and `−k` on its image A′ — positive-only emission then fires on exactly one
copy. Not a closedness bug: the positive-case rule meeting a two-copy mesh without a
quotient face object.

**Tie algebra for an in-plane quotient cut face (Grok F1, high confidence):**
- **(i′) one-sided emission + straight orig↔orig ties = CORRECT, minimal.** Emitting
  side's cut element uses `φ′ = φ + I`; neighbor copy keeps originals; quotient
  jump = I exactly once. No dup ties needed.
- (ii′) mirrored emission + crosswise ties + OPPOSITE hanging signs (+I/−I) — also
  correct, but a new condensation convention (all CutSet weights are +1 today,
  `cl_CutSet.cpp:48`); heavy.
- (iii) mirrored + straight ties: **cancels to jump 0** (both copies land on the same
  covering sheet) — silently DESTROYS the essential in-plane cut. Not 2I. Reject.
- mirrored + crosswise + same signs: algebraic contradiction (2I = 0). Reject.

**Recommended implementation of (b): quotient "XOR-once" face-flag reconciliation**,
immediately after the positive-case flag loop and before self-cancel/collect:
for each matched period face pair (F_m, F_s): both flagged → unflag both (quotient
self-cancel, mirrors interior `:461-469`); exactly one flagged → keep a canonical
single emitter (prefer master), clear the other. Straight ties; no new dup classes.
Codex V2 (high): implementable right there — `Face` carries `mPeriodic` /
`set_periodic()` and the face ties exist at CutProcessor time.

**Healing scoped (Grok F2, high):** duplicate-multiplicity equality across the pair
is NOT an invariant — transversal ONE_SIDED_JUMP mismatches are legitimate and must
NOT be healed (blanket healing creates spurious jumps or erases real ones). Christian's
deterministic healing survives only as a narrow repair for mismatches CLASSIFIED as
in-plane emission asymmetry — and XOR-once at the source supersedes it.

**(b) does not retire (d):** after XOR-once, transversal original↔duplicate edge
correspondences remain and still need the facet-mediated pairing in `match_edges`
(the direction-(d) first-contact lesson). (a) — pruning in-plane cap faces — stays
valid where the representative allows isotoping the cut off the identification;
orthogonal, cheapest when true.

**Implementation placement (Codex V2, high):** the reconciliation goes after the
per-CutData positive-case flagging loop and BEFORE the count/collect pass
(`cl_CutProcessor.cpp:425-440`) — doing it after `tFaces` is collected would let a
mirrored/moved flag miss the temporary container. Touch only cap faces with
non-null `Face::periodic()`, process each pair exactly once. Carriers verified:
`Face` and `Facet` both hold `mPeriodic`/`set_periodic()`/`is_periodic()`
(`cl_Face.hpp:40,137-146`, `cl_Facet.hpp:34-35`); ties are set by
`match_facets_and_faces()` during `CutFactory::run()`'s periodicity creation
(before `CutProcessor` is constructed), and `fix_face_slaves()` does NOT clear
`Face::periodic()`.

**Christian's "single-layer corc: coincidence or regression?" — ANSWERED (Codex V3):**
no code invariant ever enforced cap-signature correspondence (high confidence) —
the theory doc itself records that conjugated-facet membership can differ between
periodic partners (`thick_thin_cuts_and_conjugate_edges.md:321-326`). Nothing was
broken on the way; single-layer corc's success was empirical (mild cut/cap
interaction under the old count-only assert), not a stable topological guarantee
(medium confidence on "lucky"). The double-layer case simply has enough cap
crossings to make the missing invariant visible.

## Codex E1–E5 audit of the direction-(d) implementation (2026-07-16)

- **E1 — BLOCKING for the next phase:** one-sided (untied) edges currently remain in
  `Periodicity::master_edges()/slave_edges()`, whose contract downstream is still
  "matched pairs": `set_entity_dependencies()` requires `periodic() != nullptr`
  (`cl_Mesh_Periodicity.cpp:77-90`), BFM save dereferences `periodic()->id()`
  without a null check (`cl_Mesh_BfmFile.cpp:1859-1875`), crosslink/proto assume
  positional pair lists. Lower-blast-radius fix once pairing succeeds: compact the
  containers to tied pairs, count/log one-sided edges separately. Also (physics
  closure, medium-high): a one-sided edge is constrained through the cut relation
  only where the setup paths mark it node-hanging — not guaranteed by matching
  alone; must be verified when the formulation question is settled.
- **E2 PASS** (raw-pointer direction alignment is the correct duplicate-aware test);
  **E3 PASS** (id packing safe on LP64; add a `static_assert(sizeof(luint) >= 8)`);
  **E4 PASS** on the normal reset lifecycle (hardening: also check reciprocity and
  source-set membership); **E5 PASS** (no dependency on the removed reindexing).
- **Addendum after the first-contact run (Codex concurs, high):** raw node-tie
  induction is insufficient as the pairing rule — the needed original↔duplicate
  edge correspondence is inexpressible when `P(duplicate)` is intentionally null.
  Next implementation must first settle the formulation question (facet-mediated
  many-to-one vs symmetric emission vs jump-offset ties), then fix E1's container
  semantics.

## Audit trail

- `tmp/ai_exchange/corc_clean_divergence.md` — Step-0 (H-A..H-D), census, pivot, and
  hybrid-run sections (2026-07-15/16).
- `devlog/dl20260716_census_density_pivot.md`, `devlog/dl20260715_*` (SPFA lineage).
- Related: `todo/cut_pocket_removal_rules.md` (Rule 7), probe-4b devlogs
  (`dl20260615`/`dl20260616`), `todo/closed/periodic_thin_cut_continuity_fix.md`.
