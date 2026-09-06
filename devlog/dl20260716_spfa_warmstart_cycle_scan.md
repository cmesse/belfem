# SPFA Warm Start + Amortized Cycle Scan (corc gen-6 churn fix)

**Date:** 2026-07-16
**Purpose:** Fix the SPFA churn on double-layer corc generator 6 (6.39M dequeues, no
verdict) reported by Christian's `#probe`; upgrades to `fn_Graph_spfa.cpp` only.
**Module:** src/math/graph
**Exchange thread:** `tmp/ai_exchange/corc_clean_divergence.md`

## Field results (Christian)

- **helix (`cmake-build-debug/helix`): clean_spfa works end-to-end.** The isolated,
  seemingly random cut patches in space are cosmetic — the B-field is physical, and
  static condensation means no unnecessary DOFs are generated. SPFA is a net win over
  greedy regardless of the corc outcome.
- **Double-layer corc generator 6: SPFA churned** — probe showed 6.39M dequeues with
  the queue in a ~37k steady state, no convergence and no certificate.

## Diagnosis

Two compounding solver weaknesses (no soundness bug — both audits' feasibility
verdicts stand):

1. **Detection latency.** The per-vertex trigger (`count > tNumActive`) needs one
   vertex to be updated ~10⁵ times; at ~64 updates/vertex average, a parent-graph
   cycle could exist for hours before anyone looks for it.
2. **Cold-start convergence.** From all θ = 0, integer relaxation must walk the
   potentials down across the whole coefficient-scale excess range — on feasible
   instances with large |c| that is legitimately enormous.

## Fixes (`fn_Graph_spfa.cpp`)

1. **BFS spanning-forest warm start.** Along tree arcs, θ(head) = θ(tail) + w — both
   arcs of each difference-constraint pair are satisfied at init (slack 0 and 2), so
   only off-tree arcs carry violations and relax work scales with the actual excess,
   not the potential range. Tree arcs enter the predecessor graph with the tight
   invariant, so certificates can route through them. Initial queue = tails of
   violated arcs only (ring buffer doubles as the BFS queue first).
2. **Amortized periodic parent-graph cycle scan.** Every `4·tNumActive + 1024`
   updates, run the full cycle extraction from the last-updated vertex — O(1)
   amortized per update. A parent cycle is now caught within ~one scan interval of
   forming. The per-vertex trigger remains as a dormant backstop; the op-cap remains
   the loud last resort.
3. Christian's `#probe` gated to every 65536th dequeue (an `endl`-flush per dequeue
   was itself a tax).

## Verification

- `-Wall -Werror` clean with real `flags.make` flags.
- Scratch suite grown to 22 checks; new cases mimic the corc regime: a 5000-link
  chain with |c| ≈ 10⁴ (cold start would need ~5·10⁷ relaxations — now instant), a
  3000-vertex/30k-edge feasible graph with potential range ±10⁶, and a 3-cycle buried
  in a large-coefficient 30k-edge graph (exercises the periodic scan path). ALL PASS
  in 61 ms total; valgrind clean.
- Codex delta audit (D1 warm-start invariant, D2 init-queue completeness, D3 ring
  capacity, D4 termination) — verdict lands in the exchange thread.

## Expectation for the corc rerun

Either generator 6 rectifies (feasible; write-back + `|c| ≤ 1` gate as on helix), or
the periodic scan produces the edge-ID certificate within roughly one scan interval
(~4·10⁵ updates — seconds, not the 6.4M+ dequeues burned before). If the run instead
hits the op-cap `BELFEM_ERROR`, that is a solver defect report, not a mesh verdict —
bring the probe output back to the exchange thread.

## Files touched

- `src/math/graph/fn_Graph_spfa.cpp` (warm start, periodic scan, probe gate, iostream include)
- scratch `test_spfa.cpp` (+4 checks)
