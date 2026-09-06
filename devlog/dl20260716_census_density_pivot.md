# Pocket Census → Density Diagnosis → Certify-Then-Greedy Pivot

**Date:** 2026-07-16
**Purpose:** Autonomous Steps 0–2 (Christian-approved): certify the direction for the
`match_edges` 12458-vs-12498 abort, implement the pocket census + Tier-A pass, and act
on what the census revealed.
**Modules:** src/homology (+ run evidence from `cmake-build-debug/corc`)
**Exchange thread:** `tmp/ai_exchange/corc_clean_divergence.md` (Step-0, census, pivot
sections; Codex + Grok throughout)

## Step 0 — direction check (three-AI)

Claude found a hole in the naive pocket hypothesis: closedness (`dc = 0`) forbids a
one-sided pocket with a free cap trace. Grok confirmed (high) and enumerated the
admissible seam configurations (transversal wrap / in-plane band / tangency / quotient
pocket / essential cut). Codex traced the one-sided emission machinery (H-B): periodic
cap faces survive as thin-cut faces, slave seam faces are emitted one-sided from
`tFace->slave()` (`cl_CutData.cpp:159-172`), one-bit periodic pairs get unpaired
duplicates (`cl_CutSet.cpp:144-155`), and `EdgeFactory::edge_key()` keys on current
(not original) node indices — while `PeriodicityFactory::collect_nodes` folds through
`original()` and `collect_edges` counts raw objects. That is the exact abort signature
(equal node counts, +40 raw edges on one cap). Verdict: census = the right
discriminator; Tier-A flips = safe cleanup but by construction not the fix for a
cap-footprint symptom. One Codex trail discounted per project memory: the unconsumed
`mPairVerdict` (Step 6c) is a documented stale-comment trap, superseded by Step-9
original-identity keying.

## Step 1 — census (implemented + ran on corc)

`Cohomology::remove_cut_pockets( aFireTierA )` in `cl_Cohomology.cpp`: quotient
flood-fill over the zero-coefficient graph (same rep fold and master∪slave adjacency
as the audited firing helper), boundary classification with purity signs, conservative
O1 predicate (any sideset-node contact = structure), cap contact via `is_periodic()`,
other-generator overlap, small-side rule, rebuild-after-every-flip (Rule 5/M1),
per-generator seam counters (in-plane vs transversal support). Census-only logging +
Tier-A firing per the audited rules file.

**Census result (double-layer corc, refined mesh): the real defect is DENSITY, not
pockets.** Per generator: support ≈ 226'300 edges (~60% of the domain) vs ~5'406
original Smith-form support; ~4'100 in-plane cap edges; ~300 pure pockets, all
structure/overlap-contaminated (7 smeared generators cover everything) → zero Tier-A
fires; cap-pure boundary sums 466–1099, no ~40-edge candidate. The abort reproduced
identically (12458 vs 12498).

**Root cause:** the T4 write-back `c' = c − dθ` with a pure *feasibility* θ is unit
but dense — nothing in a feasibility solve drives coefficients to zero. Grok
sharpened it: the BFS-tree warm start *saturates* tree edges at |c'| = 1 (tree arcs
satisfy `θ(v) − θ(u) = 1 + c`), so density is structural, and pocket/component flips
provably cannot recover sparsity from a dense point (essential support untouched;
local exact cancellations only). Consequences: helix's cosmetic pocket scatter, the
epidemic census, and (plausibly, medium confidence) the cap chaos feeding the
one-sided emission machinery — +40 is the net residue of ~4'100 in-plane edges.

## Step 2 pivot — "certify, then greedy, dense fallback" (implemented)

Per generator in `clean_spfa()`:
1. SPFA difference-constraint solve — unchanged. Infeasible → edge-ID certificate
   abort (this WORKED: it located the too-coarse element near the domain edge that
   Christian refined — T3 verdict, Regime 2 confirmed).
2. Feasible → `cohomology::rectify_greedy_sweeps()`: the classic greedy sweeps
   (fires only where |c| ≥ 2 → support stays ~original sparse), now safe because
   feasibility is proven. Multi-trigger bail-out per Grok P2 (plateau detection alone
   is incomplete): exact plateau ∨ sweep cap `max(64, 8·n_nonunit_0)` ∨
   support-growth fuse (4× pre-greedy support), all baselines frozen pre-greedy.
3. Bail-out on a proven-feasible generator → WARNING (never a coarseness error) +
   dense θ write-back after re-solving on the mutated (same-class) coefficients.
   (Deliberate deviation from Grok P4 hygiene: re-solve instead of snapshot/restore —
   greedy mutation invalidates the stored θ; re-solve is simpler on a cold path.)
4. Pocket census + Tier-A pass unchanged after all generators.

Audits: Grok PIVOT — P1 density CONFIRMED (high), P2 hard caps REQUIRED (adopted),
P3 L1 min-cost circulation NOT now (hybrid first; Phase-4 ticket with Grok's
complexity table for successive-shortest-paths at ~4·10⁵ arcs), P4 keep θ fallback.
Codex PIVOT verdict pending at devlog time — findings land in the exchange thread.

## Run outcome (corc, hybrid)

All 14 generator instances (two cohomology passes × 7) rectified **sparse** — no
dense fallback fired. Support dropped 226k → 23.7k–53.3k per generator; in-plane cap
support 4'100 → 508–1'516. **`match_edges` still aborts, but with different numbers:
10877 vs 10728 (Δ = −149, now master-heavy) instead of 12458 vs 12498 (Δ = +40,
slave-heavy); both totals dropped.** The asymmetry scales — and even flips sign —
with the representative's cap-touching cut mass.

**Final verdict on the abort (three-AI convergent):** a genuine one-sided
seam-emission defect in cut construction (Codex H-B chain: periodic cap faces survive
as thin-cut faces `cl_CutProcessor.cpp:472-479` → one-sided slave emission
`cl_CutData.cpp:159-172` → unpaired duplicates `cl_CutSet.cpp:144-155` → raw
duplicate edge objects on one cap, while node collection folds through `original()`).
Rule 7 of the pocket rules applies: representative choice modulates the magnitude but
cannot fix it; the pocket pass correctly fired nothing. Candidate fix directions for
Christian's decision: (a) prune/route cut faces off the identification the way
non-periodic boundary faces already are; (b) symmetric two-sided emission with
periodic ties; (c) Phase-4 weighted L1 representative with heavy cap-edge weights
(cuts avoid the caps except transversal crossings, which historically work). The
density fix stands on its own merits regardless (sparse cuts, meaningful census,
helix cosmetics expected to improve).

## Files touched

- `src/homology/cl_Cohomology.hpp` (+ `remove_cut_pockets` declaration)
- `src/homology/cl_Cohomology.cpp` (`remove_cut_pockets`, `rectify_greedy_sweeps`,
  hybrid restructure of `clean_spfa`, seam counters)
- `src/math/graph/fn_Graph_spfa.cpp` (leftover probe-variable cleanup after
  Christian's probe removal)
- `todo/cut_pocket_removal_rules.md` (already audited), plan + exchange updates

## Addendum — duplicate-symmetry probe (same day, Christian's matching-tier analysis)

Christian proposed the pairing spec (originals first; equal duplicate counts per
original; canonical duplicate ordering) and asked why the node assert passes while
the edge assert fails. Answer: the node assert compares `original()`-folded counts
(blind to duplication); the edge assert is the first raw-object comparison in the
chain. The `probe_duplicate_symmetry` diagnostic (cl_Mesh_PeriodicityFactory.cpp,
before match_edges) then settled the mechanism: pre-cut update perfectly symmetric;
solve-phase has **1070 originals with bidirectional multiplicity mismatches** and
edge sumAbsDiff 1491 (net 149) — scattered ±1 in both directions, i.e. LEGITIMATE
one-sided cut membership (probe-4b lineage), not a one-sided emission bug per se.
Equal-duplicate-counts is not a physics invariant. The true defect: the NODE periodic
pipeline is duplicate-aware (ONE_SIDED_JUMP, Step-9 keying), the EDGE rebuild is not
(match_edges keys by original endpoints — silent many-to-one — and asserts raw count
equality, a cut-free-era invariant that double-layer corc's cap-crossing cuts break
honestly). Fix direction (d) recorded in `todo/periodic_cap_cut_emission.md`: induce
the edge pairing from raw node ties, inherit the one-sided policy, assert
tie-consistency instead of count equality; Christian's ordering scheme becomes part
of the tie construction.
