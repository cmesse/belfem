# Devlog 2026-09-03 — Integration Rules: Exactness Probe, `tet10` Weights, TET Dispatch, Pyramid Tables

**Date:** 2026-09-03
**Topic:** Every Gauss integration table and every `intpoints()` dispatch entry measured for polynomial
exactness; the `gauss_tet10` weights restored from Shunn & Ham 2012; the TET dispatch corrected for
orders 7 and 9; the four multi-point pyramid tables regenerated; an exactness regression test added
**AIs involved:** Claude (probe, patch, verification), Codex `gpt-5.6-terra`/high and Grok
`grok-4.6`/high (jury on the patch before it was applied)
**Claude Confidence:** high — every claim below rests on a numeric probe or an executed test, not on reading
**Literature References:** Shunn & Ham 2012 (`tmp/papers/tet.pdf`, Table 1 and Appendix F); Witherden &
Vincent 2015 and Keast 1986 as cited in the table headers; Felippa 2004 (DOI in the pyramid headers,
construction only — the paper was not on disk)
**Verification:** standalone gates, not `make check`. Strict `-fsyntax-only -std=gnu++17 -Wall -Werror
-pedantic-errors` on every new or changed source: pass. The new test linked against today's
`cmake-build-debug/lib/libbelfem.a` **fails** on TET orders 7 and 9 and PYRA orders 2–9 and nothing else;
linked with the patched `fn_intpoints.cpp` and tables it passes 8/8. `make check-fast` in the shared tree
is Christian's and is owed.

## Summary

The documentation sweep's second round had noticed that the Shunn & Ham tetrahedron tables label "order"
one higher than the degree they integrate exactly. Christian asked for the weights of `tet10` to be fixed
from the paper and for a jury-audited exactness measurement of every table. A compiled probe dumped all
57 tables and all 147 `(geometry, order)` dispatch results, and a monomial test measured each one against
closed-form moments on the reference domains the Lagrange functions use.

Three defects, one of them live in every default pyramid run:

| what | before | after |
|---|---|---|
| `gauss_tet10` weights | the barycentric coordinates copied into `aWeights` (sum 2.5); unreachable | paper's Appendix F weights ÷ 6 (the `tet20` convention); exact through degree 3 |
| TET dispatch | order 7 → `tet35` (degree 6), order 9 → `tet56` (degree 8) | orders 7–8 → `tet46` (degree 8), orders 9–10 → `tet81` (degree 10) |
| `gauss_pyra8/27/64/125` | coordinates scrambled: ∫z = 0.745 vs 1/3, ⟨x⟩ ≠ 0, duplicate points, points outside the pyramid | conical-product rules (Gauss–Legendre × Gauss–Jacobi(2,0)), exact through degree 3/5/7/9, weights positive |

Everything else measured clean: every LINE, QUAD, HEX, TRI and PENTA dispatch entry, and every other TET
entry, returns a rule whose exactness degree is at least the requested order.

## Key Findings

**"Order" meant two things in one directory.** Shunn & Ham's Table 1 numbers their rules by the exponent of
the leading error term — 10 points δ⁴, 20 δ⁶, 35 δ⁷, 56 δ⁹ — which is the exactness degree plus one. The
Witherden & Vincent, Keast and Vioreanu tables, the generic tensor-product rules, and every other geometry's
dispatch use the exactness degree. `fn_intpoints.cpp` mixed the two, so a request for order 7 (which
`auto_integration_order()` issues for every quadratic element) received a degree-6 rule. The auditors'
correction to my wording: the invariant the tree honours, and the test now locks, is *exactness degree ≥
requested order*; many rules deliberately over-deliver.

**The pyramid tables were never rules.** Their weights are the conical-product weights to fifteen digits,
which fixed the intent beyond doubt; the coordinates had the z-nodes in the x column and mapped
Gauss–Legendre nodes in the z column. The existing `test_Integration` passes order 0, which the dispatch
maps to the one-point rule — the only correct pyramid rule. Production goes through
`IntegrationData::populate`, which maps order 0 through `auto_integration_order()` to `pyra27` (linear) or
`pyra64` (quadratic): the scrambled ones. Grok found that distinction; it is why the suite was green while
every default pyramid integral was wrong.

**`tet46` over Keast `tet31` for order 7, unanimously.** Keast's degree-7 rule carries a four-point orbit
of −0.0625 weights against a total of 1/6 (I had written "a negative weight"; Grok corrected it). `tet46`
is exact through degree 8 with all weights positive at 46 points, against 35 today: +31 % quadrature cost
on quadratic tetrahedra. That cost is the one open ruling for Christian.

**A parser near-miss worth keeping.** A literal-only parse of the weight arrays made `tet35` and `tet165`
look badly normalised; they close their last weight with `1./6. - sum( aWeights )`. Reading the file
refuted the probe before it became a finding. For `tet10` both auditors independently preferred orbit-equal
17-digit literals over that closure, which would split an orbit by one ulp.

**Jury defects fixed before applying:** a missing `<algorithm>` include (Codex); tolerance headroom for
1331-point sums, now `1e-12 · max(1, |exact|)` — four decades above roundoff and four below the smallest
under-integration seen (Grok); a stale "not sure if the orders are correct" comment on the PYRA dispatch
(Grok); the `tet31` wording; and a direct test of the undispatched `tet10` table (both).

## Changes Made / Proposed

- `src/numerics/integration/fn_intpoints_gauss_tet10.hpp` — weights from the paper ÷ 6; comment states
  the paper's numbering and the degree
- `src/numerics/integration/fn_intpoints.cpp` — TET orders 7–8 → `gauss_tet46`, 9–10 → `gauss_tet81`,
  with the reason in a comment; PYRA dispatch comment rewritten
- `src/numerics/integration/fn_intpoints_gauss_pyra{8,27,64,125}.hpp` — regenerated conical-product rules
  on the reference pyramid of `cl_IF_PYRA5.hpp` (base [−1,1]² at ζ = 0, apex (0,0,1)); generator:
  `scipy.special.roots_legendre` / `roots_jacobi(n, 2, 0)`, weights `w_ξ w_η w_ζ / 8`
- `src/numerics/integration/fn_intpoints_gauss_tet{20,35,56}.hpp` — comments: Shunn & Ham numbering,
  measured degree, "UNUSED since 2026-09-03" for the two the dispatch dropped (tables kept)
- `tests/fem/test_IntegrationExactness.cpp` (new) + `tests/fem/CMakeLists.txt` — for every geometry and
  every order the dispatch provides, every monomial of degree ≤ order against its closed-form moment;
  plus the `tet10` table directly
- `src/fem/interpolation/doc/interpolation_usage_guide.md` — the paragraph denying a polynomial-exactness
  guarantee replaced by the guarantee the test now locks
- Artifacts for re-verification: `tmp/ai_exchange/intpoints/` (probes, dump, analysis, proposed patch)

## Open Questions

- The +31 % quadrature cost on quadratic TET elements (35 → 46 points). Alternative: Keast `tet31`
  (degree 7, 31 points, four negative weights). Christian's ruling.
- Whether the pyramid defect deserves a DR row: it silently corrupted volume and first moments of every
  default pyramid run for as long as the tables have existed. No shipped example uses pyramids; gmsh
  produces them in hex–tet transitions.
- `make check-fast` in the shared tree, both matrix backends. The standalone gates ran on the Blaze,
  `USE_DEBUG=ON` tree only.
- Not measured this round: `IntegrationScheme::LOBATTO` (LINE only) and `GAUSSCLASSIC`; the negative
  weights already present in `tet5`, `tet11`, `tri4` (unchanged, pre-existing).
- `fn_intpoints.cpp` includes `fn_intpoints_gauss_tet35.hpp` twice (lines 59 and 64) — harmless, left.

## Files Updated

- src/numerics/integration/fn_intpoints.cpp, fn_intpoints_gauss_tet10.hpp, fn_intpoints_gauss_tet20.hpp,
  fn_intpoints_gauss_tet35.hpp, fn_intpoints_gauss_tet56.hpp, fn_intpoints_gauss_pyra8.hpp,
  fn_intpoints_gauss_pyra27.hpp, fn_intpoints_gauss_pyra64.hpp, fn_intpoints_gauss_pyra125.hpp
- tests/fem/test_IntegrationExactness.cpp (new), tests/fem/CMakeLists.txt
- src/fem/interpolation/doc/interpolation_usage_guide.md
- tmp/ai_exchange/doxygen_contradiction_sweep.md (pre-registration, evidence, resolution)
- tmp/ai_exchange/intpoints/ (probeA.cpp, probeB.cpp, dumpA.txt, analyze.py, analysis.txt, new/, proposed.patch)
- todo/doxygen_contradiction_sweep.md (D11, D12 closed; pyramid row added)
- devlog/dl20260903_intpoints_exactness.md (this file), devlog/README.md
