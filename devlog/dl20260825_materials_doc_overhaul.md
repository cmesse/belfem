# Materials Documentation Overhaul

**Date:** 2026-08-25
**Topic:** Three-AI overhaul of every document in `src/physics/materials/doc/` — fact inventory,
blind parallel audits (Codex, Grok), rewrite, Codex language sweep — plus one kernel defect the
audit trail surfaced on the way.
**Participants:** Christian, Claude (Fable), Codex, Grok
**Exchange record:** `tmp/ai_exchange/materials_doc_overhaul.md` (fact sheet, both audits,
verification, reconciliation); sweep results in `tmp/ai_exchange/materials_doc_sweep_*.md`

## 1. Why

The seven documents dated from January–June 2026 and predated the cryogenic α branch, the
per-metal Kohler curves, the Wachtman moduli, the c_v Debye inversion and the formula alloys.
Worse than stale: the README listed 18 pure metals, Bi2212, CuCrZr, Kapton and G10 — none of
which exist; the usage guide and contracts sketched a `Metal` class (`mRRR`, `Database*`
members, `populate_rho_lambda_databases()`, a Fortran `debye_cp`) that was invented; the
field-dependent accessors were documented as `( B, angle, T )` throughout (the real order is
`( T, B, beta )`, and every argument is a `real`, so the swap compiles); the two alloy notes
described as a proposal what the code had implemented in February; and the Callaway guide
carried seven "author to confirm" hedges.

## 2. Method

Phase A: Claude wrote a frozen fact sheet from the source. Phase B: Codex and Grok audited all
seven files blind, in parallel, with a shared checklist (contradictions with source, missing
features, proposal-vs-implemented drift, every hedge settled, structure). Phase C: Claude
rewrote from the reconciled findings. Phase D: Codex swept the prose of each rewritten file
with an explicit do-not-change list. Every `file:line` citation in the new set is generated
from a regex anchor by `scripts/resolve_doc_cites.py` (`@@file|pattern@@` placeholders) and
verified with its `--check` mode — line numbers are no longer typed.

Grok's audit corrected the fact sheet on three points (all nine metals have Kohler curves,
not five; `create_cp` has a quintic beam segment and the second Bézier is `mCpBezierMedium`;
`Metal::rho( T, B, beta )` evaluates without a table while `Alloy`'s asserts one). Codex,
auditing a moving target, added the cache-file name (`Copper_RRR100.hdf5`, the class label),
`load_bh_curve` vs `set_bh_curve` (only the former routes H, μ, dμ/dH — a reader following the
old guide would have got an unrouted curve), and the two `JcFunction` arities. All folded in.

## 3. The document set now

| File | Status |
|---|---|
| `README.md` | rewritten: index, real roster (9 metals, Hastelloy, YBCO, Magnesia, formula alloys), shortest correct call, data-file search, source map |
| `materials_usage_guide.md` | rewritten (rev. 2.0): pitfalls, class tree, real API, the `material` CLI, user-defined ABI, dependency routing, the physical models as implemented (Bloch–Grüneisen with n = 4.5 for Fe/Ni, Hust, Kohler + Pippard, Wachtman E + Grüneisen ν, modified Kim as coded), property table generated from `cl_Material.hpp` |
| `materials_contracts_and_invariants.md` | rewritten (rev. 2.0): sixteen contracts, each naming what enforces it or saying nothing does; the unsourced timing table dropped; thread safety downgraded from contract to design intent |
| `alloy_homogenization.md` | new; replaces `alloy_mechanical_mixing.md` and `alloy_transport_mixing.md` (deleted): what is mixed and how, as implemented; **what the homogenizer can and cannot do** (trustworthy for eutectic solders, Cu–Ag, bronze, Al–Cu, cupronickel, Ni–Cr; wrong magnetically for ferritic Fe–Cr; wrong for austenite, Invar, intermetallics) with the Fe71Cr19Ni10-vs-NIST-304 worked example |
| `thermal_expansion_from_heat_capacity.md` | corrected: cp segment names, optional third Bézier, Lead bullet, Debye-from-c_v, Wachtman provenance of the K in the diagnostic, citations |
| `callaway_thermal_conductivity.md` | rewritten: every hedge settled, v_g formula, drifted Fortran citations regenerated, known-defect note |

`doc/doxygen_nav.dox` updated (new alloy page and the thermal-expansion page, which had never
been in the navigation).

## 4. Defect found on the way — `debye.f90` optical channel

`omega_opt = 100·params(20)·h·c` is already an energy (h·c·ν̃ ≈ 1.0×10⁻²⁰ J at 501 cm⁻¹); the
optical branch then forms `U = hbar·omega_opt/(kB·T)` — a second ħ — so U ≈ 8×10⁻³²/T, V → 1,
and K(5) → 0: the optical phonon–electron channel is inert. YBCO activates it (λ_opt = 1.06,
fitted), so that coupling absorbed a dead term. Correct: `U = omega_opt/(kB*T)` (≈ 9.4 at 77 K),
after which YBCO's λ(T) must be refitted. Confirmed independently by both auditors. **Kernel
untouched; documented as a known defect; decision Christian's.**

## 5. Open, for Christian

- The `omega_opt` fix and the YBCO λ refit (§4).
- Thermal expansion of the non-metals: Hastelloy and YBCO use polynomials that are linear in T
  at low temperature (α(20 K) ≈ 2.3×10⁻⁶ and 1.5×10⁻⁶ against physical values an order of
  magnitude lower); Magnesia uses the pre-fix flat-start Bézier. Same defect the metals had;
  matters for pointwise α below ~30 K, not for the 293→4 K contraction. Recommended order:
  Magnesia (drop-in `create_cryo_expansion`), Hastelloy (refit ΔL/L as a Bézier first), YBCO
  (anisotropy is its larger error).
- Whether read-only concurrent material queries should be promoted to a contract (needs a
  test) or the philosophy's "not thread-safe" stands.
- Whether the alloy theory should stay one document (as written) or return to two (Grok's
  preference).
- No test exercises the roster or the alloy path; a `check-fast` test constructing every
  built-in and `Sn60Pb40` / `Fe71Cr19Ni10` with loose bounds against the tabulated values would
  make the documentation's numbers executable.

## 6. Status

Documentation: reviewed by three voices, prose swept by Codex, citations mechanically verified.
No code changed. Evidence level for every claim in the new documents is source trace or a
labelled estimate.
