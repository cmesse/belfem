# Implementation Plan: Non-Unit Thick-Cut Rectification (`rectify_to_unit_or_certify`)

> **CLOSED 2026-09-03** (todo/ currentness sweep, round 3): Phases 0–3 are in the tree with tests; Phase 4 (min-cost L1 representative) is carried as option 1 of `deferred/cut_representative_options.md`. Status lines and checkboxes below are as they stood at closure and are not maintained.

**Date:** 2026-06-22 (updated 2026-07-01 — second audit round, QC/null-complex resolutions;
ordered step tracker T1–T9 added)
**Purpose:** Replace the greedy `Cohomology::clean()` with a global feasibility solve that
rectifies a non-unit thick-cut cochain (`|c(e)| ≥ 2`) to unit coefficients when a unit
representative exists, or returns a constructive negative-cycle certificate when none does.
**Module:** src/homology
**Status:** IMPLEMENTED + PIVOTED 2026-07-16. T1-T7 done; T3 verdict = Regime 2
(certificate located the corc throat; refined mesh cuts). **Architecture pivot
(census-driven, Grok+Codex audited):** a pure feasibility θ write-back is unit but
DENSE (corc: 226k support vs 5.4k Smith-form) — `clean_spfa()` is now
**certify-then-greedy**: SPFA proves feasibility (or aborts with the edge
certificate); feasible generators are rectified by the classic greedy sweeps
(`cohomology::rectify_greedy_sweeps`, sparse, multi-trigger caps: plateau ∨ frozen
sweep cap ∨ support-growth fuse); the dense θ write-back survives only as a re-solved
fallback + WARNING. Pocket census/Tier-A pass (`remove_cut_pockets`, rules file) runs
after. **Phase-4 min-cost L1 representative is upgraded from optional to the planned
product end-state** (both auditors; Grok's complexity table in the exchange thread).
Corc now rectifies sparse (23.7k-53.3k) but `match_edges` still aborts with
representative-dependent numbers → genuine one-sided seam-emission defect, tracked in
`todo/closed/periodic_cap_cut_emission.md`. ~~Open: T8 test port, T9 single-corc
regression.~~ **T8 + T9 closed 2026-08-24 (DR-23 struck); Phase 4 is the only
open item, see the 2026-08-25 currency note below.**

**2026-08-09 currentness sweep.** Phase 0/1/2 boxes are now ticked to match reality —
`clean_spfa()`, `rectify_greedy_sweeps()`, `fire_node_coboundary()` and
`remove_cut_pockets()` are all in `cl_Cohomology.{hpp,cpp}` and **committed** (`ce6a0e8b`
and neighbors; the "UNCOMMITTED tree" warnings elsewhere in the campaign are obsolete).
Read Phase 2's first box with the pivot in mind: the dense θ write-back it describes was
implemented, then demoted to a re-solved fallback + WARNING when the census showed it was
unit-but-dense (226k vs 5.4k support); the production path is certify-then-greedy.
The seam defect this plan surfaced at its end — `match_edges` aborting with
representative-dependent numbers — was root-caused and **fixed** under
`closed/periodic_cap_cut_emission.md`, so it is no longer a downstream blocker here.

**2026-08-25 currency note: Phase 3 is CLOSED; only Phase 4 remains.** The
paragraph below is the 2026-08-11 state, kept as written. Since then: the
cohomology-layer suite landed 2026-08-13 (`tests/homology/`, annulus ×
3 algorithms + two 3D periodic-bar fixtures), ran under ctest 2026-08-23
(check-fast 9/9), and the last gap (a layer-tier forced-repair fixture)
closed 2026-08-24 with `Cohomology.AnnulusForcedRepair` (red AND green by
execution). T9 was discharged 2026-08-24 by Christian's live 4-rank periodic
corc run. **DR-23 is STRUCK**; see `debt_register.md` for the full trail.
What is genuinely left is **Phase 4 only** (min-cost L1 representative, the
planned product end-state per both auditors, unscheduled). Related standing
verdict (2026-08-25 jury, `devlog/dl20260825_spfa_double_run_jury.md`): the
double `clean_spfa()` run per `compute_cohomologies()` is required: do not
fold the constructor clean into the factory clean; hygiene residues from that
round are registered as DR-107.

**What is genuinely left: Phase 3 only.** *(historical, 2026-08-11, superseded
by the currency note above)* **T8 is now PART DONE (2026-08-11).** The SPFA
solver half is ported: `tests/math/test_GraphSpfa.cpp`, 16 tests wired into the `math`
suite (`fast` label), covering the degenerate/structural cases, self loops, the cycle
distinctions, 300 randomized instances and two stress cases — with both regression guards
the original scratch harness earned (a zero-weight cycle must read FEASIBLE, and large
feasible instances must survive the spurious `cnt > n` trigger). Verdicts are verified
independently: every constraint is re-tested against the returned theta, and every
certificate is checked to be a closed walk of strictly negative weight.
**Still owed on T8:** the `Cohomology`-level entry points — `clean_spfa()`,
`rectify_greedy_sweeps()`, `fire_node_coboundary()`, `remove_cut_pockets()` — which need a
mesh fixture, exactly as `dl20260715` predicted; and **T9** (end-to-end coarse-CCT + corc
periodic regression). Both remain `debt_register.md` DR-23 and blocking-1.0: the solver now
has coverage, the cohomology layer around it still has **none**. Phase 4 (min-cost L1
representative) was upgraded by both auditors from optional to the planned product
end-state, but is not scheduled.
Thread: `tmp/ai_exchange/corc_clean_divergence.md`.
**Background / theory:** `src/homology/doc/thin_cut_nonunit_rectification.md` (algorithm),
companion `thick_thin_cuts_and_conjugate_edges.md` (cut cases), devlog
`devlog/dl20260612_thin_cut_nonunit_algorithm.md` (3-AI analysis).

---

## Execution Proposal (2026-07-03) — AI-Led Implementation

Christian proposes using a Fable 5 token promo (budget cap ≈ $100; promo ends
2026-07-08) for Claude to implement T1–T9 largely autonomously. Because the rectification
theory is outside Christian's cohomology background, the decision question is whether
Claude can implement it with minimal help. Christian will decide after discussion with
Gregory on 2026-07-06.

### Claude's calibrated confidence

- **T1–T7 (implementation to plan): ~85%.** The deep cohomology reasoning is complete
  and twice-audited (quotient fold, no-mirroring write-back, side-state invariant,
  seven-cut-case enumeration). The remaining core is a difference-constraint feasibility
  solver — queue-based Bellman–Ford with negative-cycle extraction — using standard graph
  algorithms. The plan above is detailed enough to serve as a spec: containers, arc formulas,
  write-back mechanism, and guards are settled.
- **T1–T9 end-to-end within budget: ~65–70%**, conditional on the two preconditions below.
- **Without a self-serve build/run loop: ~40%.** Debugging a solver through pasted logs
  round-trips every iteration through Christian and burns both tokens and time.

### Preconditions (the "minimal help" contract)

- [x] ~~**Dedicated build directory + standing run permission** for Claude~~ — SUPERSEDED:
  Christian granted run permission ad hoc for the corc iteration loop during the
  2026-07-15/16 campaign; the standing rule is unchanged (Christian runs builds). Original
  request: (`make` + the
  diagnostic/test executables), separate from the shared `cmake-build-debug` tree, to close
  the debug loop without collisions.
- [x] **A non-unit reproducer input** — RESOLVED 2026-07-15: the **double-layer corc**
  (`cmake-build-debug/corc`) is heavily non-unit (gen 0 starts 5406/5406 non-unit; gen 6
  stalls greedy `clean()` at 55'511 non-units / 194'288 support, tripping the
  `BELFEM_ERROR` at `cl_Cohomology.cpp:316`). ~~corc does not qualify: its generators are
  unit-coefficient (see Phase 0 caveat)~~ — that caveat applied to the single-layer corc,
  which remains the periodic regression case. Audit thread:
  `tmp/ai_exchange/corc_clean_divergence.md` (Codex + Grok concur, 2026-07-15).
- [x] **Standing edit approval** scoped to the T1–T9 files — granted and used (T1–T7 landed) (`cl_Cohomology.hpp/.cpp` + test
  scaffolding).
- [x] **Christian available at exactly two decision points:** the T3 gate and final sign-off.
  *(T3 gate taken 2026-07-16 — verdict Regime 2, certificate located the corc throat.)*

### Session structure and stop-loss

1. **T1 + T2 first** (solver core + read-only diagnostic) — cheapest steps, highest
   information; ends with the T3 feasible/infeasible verdict reported to Christian.
2. **Stop-loss:** if token burn after T1/T2 projects past budget, stop there — the read-only
   diagnostic is a self-contained deliverable (correct diagnosis instead of a silent hang).
3. **T4–T7** in one push (write-back, guards, messaging, error split), Codex-audited.
4. **T8–T9** (unit tests before end-to-end validation, in that order, to keep runtime
   debugging off the large meshes).

### Expectation setting

A T3 **infeasible** verdict is still a successful outcome. In that case, the deliverable
is a clean `BELFEM_ERROR` with an edge-ID certificate identifying the throat to refine,
replacing today's silent non-termination. The coarse mesh still will not produce cuts until
it is refined; "rectified and running" occurs only when the input is feasible.

---

## Design decisions (settled)

- [x] **Runs on the mesh 1-skeleton directly**, not through a rebuilt `graph::Vertex` object: nodes are
  vertices (`Node::edge(j)`), edges are arcs (`Edge::node(0/1)`), and coefficients come from the
  `Cochain`. The active `math/graph` class is unweighted/undirected with no weighted SSSP, so it
  adds no value here; the archived `fn_Graph_tarjan.*` is a *reference* for cycle bookkeeping only
  (unweighted, combinatorial — not a drop-in for negative-cycle extraction).
- [x] **Solver state lives in external arrays keyed by the vertex index** — `Vertex` is left
  untouched (no new members, `mLevel` kept). This matches the mesh's existing pattern for
  per-entity scratch. The arrays are transient: allocated per solve and freed afterward.
- [x] **Algorithm:** queue-based Bellman–Ford–Moore (SPFA) for difference-constraint feasibility,
  with negative-cycle detection and edge-list extraction. (`θ` is signed — potentials go negative.)
- [x] **Periodic quotient:** fold each slave node onto its `periodic()` master *inline* (no map),
  so periodic partners share one `θ` slot (mirrors the existing slave-flag / `->periodic()`
  redirect in `clean()`). This is required for correct certificates on PBC meshes. Use the typed
  `Node::periodic()` accessor — **no `reinterpret_cast`** (it skips base-subobject offset
  adjustment and is UB-prone here; an implicit upcast or `static_cast` is the only legitimate
  hierarchy cast, and the fold needs none).
- [x] **Keep `clean()`'s signature and call sites**; replace only its body.
- [x] **Realization details (Codex + Grok audited, 2026-06-22):** loop over `mGenerators(1)` with
  per-generator reset; relax only the `original_edges()` domain; scan master∪slave adjacency for
  the quotient; edge-ID certificate; `BELFEM_ERROR` guards; `Cell`/ring-buffer state, int64
  distances. Details in the sections below; provenance under Risks.
- [x] **Shared solver core (2026-07-01):** implement the SPFA + certificate extraction once, as an
  internal `feasibility_solve()`; the Phase 0 read-only diagnostic and the Phase 2 production
  `clean()` body are thin wrappers over it. This eliminates the "diagnostic must mirror
  production semantics" risk: the diagnostic *is* production minus the write-back.
  (Codex + Grok concur, 2026-07-01 audit.)
- [x] **`aPeriodicity == true` is a hard precondition (2026-07-01):** the quotient solve is valid
  only for a complex built with periodicity enabled. `CutFactory` always passes `true`
  (`cl_CutFactory.cpp:430`), so slave seam edges never enter `original_edges()`
  (`cl_SimplicialComplex.cpp:116-131`, `:468-471`); a `false`-built complex on a periodic mesh
  would desync the solver's `rep()` fold from the edge domain. Guard or document this as a
  precondition; do not treat it as an optimization. (Both auditors flagged independently, 2026-07-01.)

## State model (external, index-keyed)

No compaction map: arrays are keyed directly by `node->index()` (the vertex index *is* the
offset). Size them to the node count; slots for nodes outside the cohomology domain stay
untouched. The periodic quotient needs no map either — a slave is canonicalized to its master's
index *inline* via `periodic()` at each access (the pointer is the "map"), so both endpoints of
an arc relax into the master's slot. (Slave slots are dead under the quotient — never read for
output.)

**Containers (BELFEM-fit, per audit):** `Vector`/`Matrix` are linear-algebra-only — use `Cell`
or raw arrays for the scratch. The SPFA worklist must be a **preallocated ring buffer**
(`Cell`/raw + head/tail), *not* `containers/Queue` (it wraps `std::queue`, so `push()` allocates
in the relax loop). Distances are **`int64_t`** (guard against signed overflow on long negative
paths; coefficients write back as `int`, safe since a feasible result has `|c|≤1`).

```text
N          = number of mesh nodes (arrays sized to N, indexed by node->index())
rep(node)  = node->is_flagged() ? node->periodic()->index() : node->index()  // slaves pre-flagged; no map
theta(N)   : Cell<int64_t>   potentials θ (signed, wide)
parentV(N) : Cell<index_t>   predecessor representative (cycle walk-back)
parentE(N) : Cell<index_t>   predecessor mesh-edge id  + a sign/direction bit  (edge-ID certificate)
inQueue(N) : DynamicBitset   SPFA membership
count(N)   : Cell<index_t>   relax counter (neg-cycle trigger when > n_quot)
queue      : Cell<index_t>   preallocated ring buffer of active reps (head/tail), NOT std::queue
```

`rep()` follows `clean()`'s pattern: flag `periodicity()->slave_nodes()` first, then redirect
via `is_flagged()` + `periodic()`. `n_quot` = count of distinct `rep` indices (the quotient
vertex count); `number_of_nodes()` is a safe over-estimate if a tighter bound is inconvenient.
(Assumes `node->index()` is the stable mesh node index at solve time; do not run an
index-mutating graph pass concurrently.)

Constraint arcs are **synthesized on the fly** over the **cohomology edge domain only**
(`mSimplicialComplex->original_edges()`; a null complex is a `BELFEM_ERROR` — the
field-constructor path is dead code, see Risks) — never over raw mesh edges with absent
coefficients read as 0. For edge
`e = (i, j)` with `i = node(0)`, `j = node(1)`, coefficient `c(e)`, and `u = rep(i)`, `v = rep(j)`:
`θ(v) − θ(u) ≤ 1 + c(e)`  (arc u→v, weight `1+c(e)`)  and
`θ(u) − θ(v) ≤ 1 − c(e)`  (arc v→u, weight `1−c(e)`).
`c(e)` from `Cochain::getCoefficient(edge->index())`; the sign convention matches the edge's
intrinsic `node(0)→node(1)` orientation (same convention `clean()` uses — do **not** flip the
weight sign from the coefficient sign).

**Folded self-loop = length-1 certificate.** If `u == v` (an edge whose endpoints fold to the
same rep) and `|c(e)| ≥ 2`, the constraint `0 ≤ 1 − |c|` is violated immediately: this *is* an
infeasibility certificate (the theory's `|⟨c,z⟩| > length(z)` with `length = 1`). Detect and
emit it as a one-edge certificate rather than feeding a degenerate self-loop into SPFA.

---

## Phase 0 — Diagnostic first (QB), no pipeline changes

- [x] Implement the feasibility check (SPFA + negative-cycle detection) as a **read-only
  diagnostic** on a real **coarse-CCT** generator, reporting feasible vs. infeasible. To report
  a certificate loop length, include the minimal parent walk-back; otherwise report only the
  boolean. Production-semantics mirroring comes from the shared `feasibility_solve()` core
  (see Design decisions). Hook it right after generator computation, before `clean()`,
  behind a flag.
- [x] Decide between "feasible-but-greedy-`clean()`-failed" and "genuinely infeasible." Cheap; touches
  nothing. Record the verdict in the devlog before proceeding.
- [x] **Reproducer caveat (2026-07-01, superseded 2026-07-15):** the unit-coefficient
  finding (`todo/closed/periodic_thin_cut_continuity_fix.md`) applied to the
  **single-layer** corc, which stays a periodic regression check. The **double-layer**
  corc is a live non-unit reproducer: generator 6 stalls greedy `clean()`
  (4042 → dip to 1869 → monotone growth to 55'511 non-units, support 29k → 194k).
  Coarse-CCT remains a useful second Regime-2 candidate if gen 6 proves feasible.
- [x] **Do not infer the regime from a hang:** whether greedy `clean()` can cycle on a *feasible*
  input remains unproven (auditor split, 2026-07-01: Grok's simulation of the exact firing
  rule found 0 feasible cycles in ~1700 instances; Codex argues the firing is a monotone
  BF-style relaxation that should terminate). A hang is *consistent with* Regime 2, but does not
  certify it; only the SPFA verdict does.

## Phase 1 — Core solver (new internal function)

- [x] **Loop over generators** `mGenerators(1)`, re-initializing `theta/parentV/parentE/inQueue/
  count` for each (a separate cochain per generator).
- [x] Size the scratch to the node count; index by `node->index()`. Pre-flag periodic
  `slave_nodes()` so `rep()` folds inline (no map).
- [x] Initialize the difference-constraint problem (virtual super-source ⇔ all `θ = 0`, all
  *representative* nodes enqueued).
- [x] SPFA relaxation over the **cohomology edge domain** (`original_edges()` filter). When a rep
  `u` is dequeued, scan the incident edges of **both** the master and (if periodic) its slave
  partner; endpoint folding alone does not cover the full quotient adjacency (`clean()` does this
  second pass, `cl_Cohomology.cpp:269-301`). Canonicalize both endpoints via `rep()`; all reads/
  writes use rep indices only. **Equivalence verified (2026-07-01, three-way):** this rule
  exactly reproduces the quotient coboundary built by `create_complex()`
  (`cl_SimplicialComplex.cpp:144-188`) and the edge set fired by today's `clean()`; sign
  conventions match. `clean()`'s extra `!is_flagged()` skip (`cl_Cohomology.cpp:277`) is redundant given the
  domain filter (slave seam edges are never in `original_edges()`) — the solver does **not**
  need to replicate it.
- [x] Negative-cycle detection (relax-count > `n_quot`; defer Tarjan subtree-disassembly early-exit
  to Phase 4). Also emit the length-1 certificate for a folded self-loop with `|c| ≥ 2`.
- [x] Negative-cycle **extraction** as mesh edge IDs: walk `parentE/parentV` (advance `n_quot`
  steps to enter the cycle, then collect until a rep repeats). Storing `parentE` (edge id + sign)
  is essential — node-only parents lose which oriented arc triggered the relax (parallel edges /
  two arcs per edge). The archived Tarjan is structural reference only.

## Phase 2 — Rectify-or-certify wrapper (replaces `clean()` body)

- [x] **Feasible:** apply `c ← c − dθ` over **all** cohomology-domain edges (not just the current
  support — `dθ` creates new `±1` entries where `c = 0`), using `delta = θ[rep(j)] − θ[rep(i)]`.
  Keep `Cochain` **auxiliary state consistent** (`clean()` updates the boundary cochain manually;
  `addSimplexToCochain` only touches the map) — apply via a `dθ` cochain + `addCochainToCochain`,
  or replicate clean()'s boundary updates; erase zeroed entries. The cochain stays in the same class.
- [x] Add the global post-rectification guards **before `CutData`**, as **`BELFEM_ERROR`** (always
  active — `BELFEM_ASSERT` compiles out in release and this is a production safety gate),
  evaluated on the same `original_edges()` / quotient face domain as the solver: closedness
  (`dc = 0`, per-face oriented sums) and `|c| ≤ 1` everywhere. (Currently absent — the existing
  `check()` only tests *flagged* faces, `cl_Cohomology.cpp:330` — so a surviving non-unit edge
  fails late and undiagnosably.)
- [x] **Infeasible:** fail with a release-active `BELFEM_ERROR` that reports the certificate
  loop as **edge IDs** (node coordinates of the loop are useful secondary output for locating
  the throat in the mesher) and recommends targeted refinement, instead of the current silent
  non-termination / generic "invalid pattern". First-iteration Regime-2 handling remains manual:
  the user refines the reported throat and re-runs; automatic certificate-guided refinement
  stays Phase 4.
- [x] Hand the rectified unit cochain to the existing thin-cut pipeline **unchanged** (justified
  by the seven-cut-case enumeration, §6 of the tech note).

## Phase 3 — Tests & validation

- [ ] Unit: small **feasible** cochain → rectifies to `|c| ≤ 1`, class preserved.
- [ ] Unit: **infeasible** directed 3-cycle with all coefficients 2 → returns the cycle as certificate (the
  case where greedy `clean()` loops forever).
- [ ] Unit: **periodic** case → quotient fold gives correct θ/certificate across the seam.
- [ ] Validation: re-run the coarse-CCT (non-unit) reproducer end-to-end; re-run corc as a
  **periodic regression** check. Its generators are unit-coefficient, so it exercises the
  quotient fold rather than rectification.

## Phase 4 — Optional refinements (defer)

- [ ] Min-cost circulation for the *optimal* representative (weights `w_e` = dual-cell area), if a
  minimal/short cut is ever wanted — the feasibility solve is the `w ≡ 1` constrained case.
- [ ] Certificate-guided local refinement (option B): split along the certificate loop, re-run.
- [ ] Performance: localize the search to the cut support (`c ≠ 0` region) on large meshes;
  consider SPFA→Tarjan early-exit tuning.

---

## Ordered Implementation Steps (T1–T9)

Actionable tracker, in execution order (style per `todo/closed/meshfile_refactor_plan.md` §4). The
phases above carry the rationale and realization details; this section keeps the work
sequence explicit. Source edits start at T1 and require the usual explicit approval per
session.

- [x] **T1 — Scaffold the shared solver core `feasibility_solve()`** *(Phase 1)*.
  DONE 2026-07-15 with one design change (Christian's Q3): the solver lives in
  **`src/math/graph/fn_Graph_spfa.{hpp,cpp}`** as the generic
  `graph::spfa_difference_constraints( nV, tails, heads, weights, theta, cycle )`
  (arc-list API; homology synthesizes the arcs), not as a `Cohomology` member —
  reusable and unit-testable in isolation. Includes: virtual super source ( all θ = 0 ),
  CSR arc buckets, preallocated ring buffer (no `std::queue`), `int64_t` distances,
  negative-self-arc fast path (the folded-self-loop length-1 certificate),
  update-count trigger with **verified** predecessor-graph cycle extraction
  (spurious triggers continue — scratch-test-caught soundness fix), hard op-cap
  `BELFEM_ERROR` failsafe. `Vertex::mLevel` deliberately not used: potentials are
  signed 64-bit, `mLevel` is `index_t`, and the plan keeps `Vertex` untouched.
- [x] **T2 — Phase-0 read-only diagnostic** *(Phase 0)*. SUPERSEDED by Christian's
  2026-07-15 scaffold: `Cohomology::clean()` dispatches through the `mFunClean`
  member-function pointer with `clean_greedy()` / `clean_spfa()` selectable — the
  pointer IS the switch, and the certificate/diagnostic reporting lives in
  `clean_spfa()` itself (per-generator gLog Detailed message; infeasible →
  `BELFEM_ERROR` with edge-ID certificate). The T3 verdict now comes from running
  the double-layer corc with the spfa default (T9).
- [x] **T3 — Decision gate: VERDICT DELIVERED 2026-07-16 — Regime 2 CONFIRMED.**
  Generator 6 of the double-layer corc was **genuinely infeasible**: the certificate
  (plus Christian's nearest-element helper for visual debugging) localized the
  obstruction to an element near the domain edge — too coarse exactly there, which
  nobody would have guessed by inspection. After targeted refinement the mesh
  produces a **valid cut**. The rectify-or-certify pipeline is validated end-to-end
  (helix: rectification path; corc: certificate path). Both prior regimes of the
  greedy stall question are now settled empirically.
- [ ] **T3 — Decision gate on the QB verdict** *(Phase 0)*. Feasible → proceed to T4.
  Infeasible → T4–T6 still land (correct diagnosis instead of a hang), but discuss
  refinement of the reported throat with Christian before further automation.
- [x] **T4 — Replace the `clean()` body (write-back)** *(Phase 2)*. DONE 2026-07-15:
  `clean_spfa()` fires every representative node `u` with weight `−θ(u)` through
  `cohomology::fire_node_coboundary()` — a static helper replicating `clean_greedy()`'s
  firing pattern (master∪slave quotient adjacency, `addSimplexToCochain` sign
  convention, manual boundary bookkeeping) scaled by the weight, so
  `Δc = θ(rep(i)) − θ(rep(j))` on every domain edge and zero entries auto-erase.
  Boundary side-state is maintained by the same scaled replication (`mCoboundary`
  stays stale exactly as under greedy — pre-existing contract). One deliberate
  deviation: the boundary node fold uses `is_flagged()` in BOTH adjacency passes
  (greedy's second pass used `is_periodic()`, which would mis-fold a master endpoint
  of another pair onto its slave).
- [x] **T5 — Global post-rectification guards** *(Phase 2)*. DONE 2026-07-15 (narrowed):
  release-active `BELFEM_ERROR` `|c| ≤ 1` over the **whole generator map** after
  write-back. The explicit `dc = 0` face check was dropped as redundant by
  construction: the write-back only ever adds integer multiples of node coboundaries,
  and `dd = 0`, so closedness is exactly preserved from the input generator
  (upstream's responsibility, as before). Reopen if a non-coboundary write-back path
  is ever added.
- [x] **T6 — Infeasible-branch messaging** *(Phase 2)*. DONE 2026-07-15: release-active
  `BELFEM_ERROR` reporting the negative-cycle certificate as mesh edge IDs (first 32,
  then "..." ) plus loop length and a targeted-refinement recommendation. Node
  coordinates were left out of the message (edge IDs locate the throat in the mesher);
  add them later if locating proves awkward in practice.
- [x] **T7 — Split the `determine_cut_case_3d()` error messages** *(diagnosability)*.
  DONE 2026-07-15 (adapted): `CutData::weight()` only reads sign bitsets and can never
  see a magnitude, so the "non-unit reached CutData" case is representationally
  impossible — the T5 guard in `clean_spfa()` IS that split. What landed instead:
  the inadmissible-pattern `default:` no longer blames mesh coarseness (that case now
  aborts earlier with an edge certificate) and names a generator/closedness defect
  with the pattern code; the two "Invalid cut coefficients" errors now print the
  actual per-edge coefficient values and name a generator orientation defect. Also
  removed a stray block-scope `collect_thin_cut_edges` declaration pasted mid-function
  (`tPattern += 4` branch; committed in `6ee74493`, harmless but accidental).
- [◐] **T8 — Unit tests** *(Phase 3)*. PARTIAL 2026-07-15: a scratch harness
  (`<scratchpad>/test_spfa.cpp`, session-local) covers the solver: (a) feasible chains
  rectify to `|c| ≤ 1`; (b) the infeasible all-2 directed 3-cycle returns a verified
  certificate (the greedy-hang case); plus self-loop (= periodic fold), disconnected
  components, 400 randomized instances (planted potentials / planted cycles), n=2000
  stress, valgrind clean. It caught one real soundness bug (spurious update-count
  trigger) before it ever reached the repo. Solver half ported 2026-08-11
  (`tests/math/test_GraphSpfa.cpp`, 16 tests). **Layer half landed 2026-08-13:**
  `tests/homology/test_Cohomology.cpp` — programmatic annulus (no mesh files), six
  cases (TRI3 and QUAD4 × Pellikka / CCR / PellikkaGeneralized coreduction), driving
  the production sequence (coreduce only, `Homology` SNF, `Cohomology` ctor,
  `updatekGeneratorsFromHomology`, `clean()`), asserting betti numbers, the
  unit-coefficient contract on the cleaned generator (a no-op repair on this
  mesh — see the open item below), and the winding pairing of the generator with
  the inner-ring cycle. Probe-run green (6/6, 45 ms) against the prebuilt libs;
  suite gate `make check-fast` with `USE_TEST=ON` still owed. Two fixture lessons
  recorded in the test header: the mesh MUST be built with
  `aComputeConnectivities = true` or the coreduction silently degenerates, and the
  `reduce_complex*` calls belong to the `mSuggestHomologies == false` branch that
  production never takes. **Periodic quotient landed 2026-08-13, same day
  (stage 2, Christian's go-ahead):** `tests/homology/test_CohomologyPeriodic.cpp`
  — two 3D periodic-bar fixtures with hand-wired periodicity (symmetric
  `set_periodic` links, id-matched edge/face pairs, slave faces re-based to
  slave role as the private `PeriodicityFactory::fix_face_slaves` does),
  matching the production infinite-tape topology: solid torus (betti 1,1,0;
  free axial generator paired ±1 through the fold) and annulus × S¹
  (betti_1 = 2; unimodular pairing matrix against the axial + encircling
  reference cycles). Asserts the CONSTRUCTOR generators:
  `updatekGeneratorsFromHomology`'s 3D branch expects terminal in/out pairs
  from `suggest_Homology`, which raw SNF generators are not — fed anyway it
  yields a det-2 pairing (condition + free cut), per its own Step-4 comment.
  Suite 8/8 green in 7.2 s at probe tier. STILL OPEN: a case that forces an
  actual non-unit repair at the layer tier (all fixture generators come out
  unit already; repair behaviour is covered at the solver tier).
- [ ] **T9 — End-to-end validation** *(Phase 3)*. Run the coarse-CCT reproducer end-to-end
  and corc as a periodic regression. Devlog the outcome; tick the QC risk item on success.

Leave Phase 4 items (optimal representative, certificate-guided refinement, support
localization) deferred and unnumbered; add steps for them only when scheduled.

---

## Risks / open questions

- [ ] **Periodic quotient correctness (QC) — largely resolved by static analysis (2026-07-01),
  runtime verification on a PBC mesh remains pending (Phase 3).** Three-way verification
  resolved: (a) the complex **is already the periodic quotient** — slave entities are unflagged
  before construction (`cl_SimplicialComplex.cpp:116-131`), so slave seam edges carry no cochain
  coefficients and add no duplicate constraints; (b) **no mirroring pass is needed in the
  write-back** — downstream never reads the cochain by slave edge keys; `CutData` propagates
  ±1 bits to `edge->periodic()->index()` itself (`cl_CutData.cpp:230-261`, plus
  `collect_edges` and the `weight()` bitset reads). Residual: that same-sign bit propagation is
  documented correct for **translational** periodic maps only
  (`todo/closed/periodic_thin_cut_continuity_fix.md`, "Translational Periodicity Sign") — a
  pre-existing scope limit, unchanged by this plan; re-check if rotational/anti-periodic
  maps ever activate.
- [x] **Domain scope — resolved by design (2026-07-01):** relaxation must cover the full
  cohomology domain, because smoothness constraints on `c = 0` edges still couple the graph.
  T1 encodes this by enqueueing every representative and relaxing every `original_edges()`
  arc; T4 writes back over that same full domain. Reopen this item if the Phase-4
  support-localization optimization is ever implemented — that change must re-establish
  correctness.
- [x] **Higher-order elements — resolved for supported orders (2026-07-01):** the cochain
  lives only on the corner-edge skeleton. `CutFactory` flags corner nodes before building the
  complex (`cl_CutFactory.cpp:397,422`, filtered by `is_flagged()` in `create_complex`); edge
  keys are formed from the corner endpoints (`edge_key()`, `cl_EdgeFactory.cpp:323`); and for
  quadratic edges `EdgeFactory::grab_nodes()` stores the midside node in slot 2
  (`cl_EdgeFactory.cpp:477-508`), so `node(0)/node(1)` are always corners. A TET10 therefore
  contributes the same six-corner-edge skeleton as a TET4, and the seven-case enumeration
  applies unchanged. Scope: orders 1–2 only — `EdgeFactory` rejects order ≥ 3
  (`cl_EdgeFactory.cpp:272`), so higher orders cannot arise. (Codex caveat, adopted.)
- [ ] **Large-mesh cost:** worst-case SPFA remains `O(VE)`. T1 includes a hard operation cap
  and release-active `BELFEM_ERROR` failsafe; support localization stays deferred to Phase 4.
  Keep this open until T9 measures a representative mesh.
- [ ] **Cochain invariant consistency — narrowed, not closed (2026-07-01; Codex refuted the
  first closure attempt).** Most post-`clean()` consumers read only `getSimplicesMap()`
  (`check()` at `cl_Cohomology.cpp:321`, debug output, `CutData` `collect_edges`/
  `collect_coefficients`). But on the `mSuggestHomologies` path,
  `updatekGeneratorsFromHomology()` runs *after* the constructor's `clean()`
  (`cl_CutFactory.cpp:516`) and recombines the cleaned generators via
  `addCochainToCochain(mGenerators(k)(j), …)` (`cl_Cohomology.cpp:586`) — which **reads the
  source generator's boundary/coboundary side-state** (`cl_Cochain.hpp:294,299`; the
  `!mIsBound`/`!mIsCobound` guards are true for generators, so the reads are live).
  Consequence: T4's write-back **must** keep the side-state coherent — apply the update via a
  `dθ` cochain + `addCochainToCochain` (which maintains it automatically) rather than touching
  the map alone. Close this item at T4 review, once that mechanism is in the code.
- [x] **Numeric width — resolved by design (2026-07-01):** T1 uses `int64_t` distances and
  `index_t` counters; the final write-back fits `int` because a feasible representative has
  `|c| ≤ 1`. Grok's simulation saw transient coefficient growth to `|c| = 7`, confirming that
  wide intermediate storage is required, not optional.
- [x] **Null `mSimplicialComplex` — resolved (2026-07-01):** the field constructor
  `Cohomology(Mesh*, Mesh*)` (`cl_Cohomology.cpp:23-30`) is **dead code**: zero call sites in
  `src/` and `nonfree/`; only the two SimplicialComplex constructors are used
  (`cl_CutFactory.cpp:445,464,479,496`). Today that path would crash anyway:
  `mSimplicialComplex` defaults to `nullptr` and is dereferenced at `clean():199`. Decision: `BELFEM_ERROR` guard on a
  null complex instead of engineering a fallback edge domain. (Codex confirmed, high.)
- [ ] **Greedy termination on feasible inputs — open theory question (2026-07-01):** no proof
  exists in code or docs. Grok's exact-rule simulation found no feasible cycle (~1700 instances,
  transient coefficient growth to |c| = 7 observed); Codex argues from a monotone-relaxation
  termination heuristic. Practical consequence (both agree): never interpret a `clean()` hang as
  a feasibility verdict — Phase 0's SPFA decides. This becomes moot for production once Phase 2
  replaces the greedy body. The transient growth observation independently supports the int64 distance choice.
  **2026-07-15 field observation:** double-layer corc gen 6 stalls the production greedy
  (dip 4042 → 1869, then monotone growth to a 55'511 plateau, support 29k → 194k); regime
  remains UNDECIDED per this item — the stall is a greedy-cascade signature, not a
  feasibility bit (three-way concurrence, `tmp/ai_exchange/corc_clean_divergence.md`).
  Confirmed mechanism: in-sweep cascading — the sweep iterates the live `std::map` while
  `addSimplexToCochain` inserts higher-index offenders that fire in the same pass —
  amplified by the sign-only endpoint rule; no iterator-invalidation or sign bug found.

> **Audit provenance:** this plan was reviewed by Codex + Grok (2026-06-22); all code-grounded
> findings re-verified against `clean()`/`Cochain` before adoption. Thread:
> `tmp/ai_exchange/thin_cut_impl_plan_audit.md`. Key adopted gaps: generator loop, `original_edges`
> domain filter, master+slave quotient adjacency, edge-ID certificate, Cochain aux-state,
> `BELFEM_ERROR` guards, `Cell`/ring-buffer containers, int64 distances, self-loop certificate.
>
> **Second audit round (2026-07-01, thread `tmp/ai_exchange/thin_cut_plan_update_audit.md`):**
> Codex (Q1–Q4) and Grok (Q2–Q3); Claude re-verified all load-bearing citations. Confirmed:
> quotient-adjacency equivalence, no-mirroring write-back, and the dead field-constructor path.
> Adopted: shared solver core, `aPeriodicity == true` precondition, corc-is-unit-coefficient
> reproducer correction, and hang-does-not-certify-regime caveat. Left open: greedy termination on
> feasible inputs (auditor split, moot after Phase 2).

## Key code references

- `Cohomology::clean()` — `src/homology/cl_Cohomology.cpp:195-309` (body to replace; periodic
  slave-flag + `->periodic()` redirect pattern at `:198-241`).
- `Cochain` — `src/homology/cl_Cochain.hpp` (`getCoefficient:165`, `addSimplexToCochain:251`,
  `getSimplicesMap`).
- Adjacency — `Node::number_of_edges()/edge()` (`cl_Vertex.hpp:293,344`), `Edge::node(0/1)`
  (`cl_Vertex.hpp:316`).
- Periodic — `Node::periodic()/is_periodic()` (`cl_Node.hpp:361,373`), `Mesh::periodicity()`.
- `CutData::determine_cut_case_3d()` — `src/homology/cl_CutData.cpp:447-660` (downstream consumer
  of the rectified unit cochain).
