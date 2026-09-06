# Coulomb Gauge Penalty Campaign — Final Report

**Date:** 2026-08-23
**Purpose:** Closing report of the overnight task "Coulomb Gauge Penalty
for the h-phi Formulation (theory + G-operator on the Nédélec edge
functions)". Records what was derived, what was audited and by whom, the
disagreements and their resolutions, what is verified versus reviewed,
and the recommended next step. Per the brief, wiring the penalty into
`mt_maxwell_h` was OUT of scope; this report ends with a recommendation
only.
**Module:** src/fem/interpolation/nedelec, tests/fem

---

## 1. Executive summary

1. **The G-operator exists, is tested, and is registered in `make
   check`** on 9 of the 10 Nédélec edge-function classes (all but
   LINE3). Final executed gates: `test_fem` 145/145, `make check` 14/14
   (100%).
2. **The theory verdict is negative for the intended use.** An
   element-local gradient penalty `γ ∫ (∇h):(∇δh)` is identically zero
   on the curl null space of TET4 / TRI3 / thin-shell meshes — the null
   space is spanned by element-wise constant fields (cohomology modes
   included), which every element-local ∇h penalty annihilates. On TET4
   the penalty is even exactly proportional to the existing curl term
   (`∫∇w:∇w = ½∫curl·curl`). The conditioning pathology of the quench
   decks therefore cannot be repaired by this penalty on those meshes.
3. **The draft formula was corrected, not assumed.** Christian's draft
   `α[μ²(GᵀG) + μ EᵀE dH/dt]` with `α = χμ/Δt` is dimensionally
   inconsistent; the audited consistent K-channel form is a
   frozen-coefficient penalty with `γ = χ·ρ*/μ²` (Δt-invariant ratio
   against the resistive term, residual-consistent, no dJdx needed).
4. **What would stabilize instead:** the weak (patch-level) divergence
   penalty `Dᵀ·M_L⁻¹·D` assembled from E and nodal gradients (Monk
   §7.4 route), or tree–cotree gauging (Dular et al. 1997). Tabled in
   the theory note §14; neither was implemented (out of scope).
5. **G remains justified** where its value is nonzero and meaningful:
   hex-dominant meshes (HEX8/HEX8TS/HEX8TB), divergence diagnostics,
   and gradient-based postprocessing.
6. **Collateral find:** DR-98 — `EF_HEX8::C()` returned the NEGATED
   curl. Pre-existing, invisible (no test coverage; K = CᵀC is
   sign-invariant). Fixed with Christian's explicit approval, closed
   with a permanent compiled Stokes gate; both auditors confirmed no
   consumer carried a compensating sign.

## 2. Deliverables

| Artifact | Content |
|---|---|
| `src/fem/maxwell/doc/coulomb_gauge_penalty_theory.md` | The theory: consistent penalty form, the simplex/thin-shell no-go, per-family divergence table, what the measured condition numbers are, and the ranked alternatives |
| `src/fem/interpolation/nedelec/cl_EF_EdgeFunction.hpp` | Pure virtual `G()`, protected `mGrad`, the single layout contract (`G(i+d·j,e) = ∂(w_e)_j/∂x_i`, explicit curl-tie rows) |
| Nine `cl_EF_*.{hpp,cpp}` | Per-element G implementations (structure table §3) |
| `tests/fem/test_EdgeFunctions.cpp`, `tests/fem/test_InterfaceOrientation.cpp` | 23 new G tests across the two batteries, plus the fixtures that made them possible |
| `tests/fem/support/cl_EF_TestVolume.hpp`, `cl_TS_TestStack.hpp` | Midside-offset parameter (first genuinely curved quadratic test geometry in BELFEM) and the TS_TestPrism corner override |
| `todo/debt_register.md` | DR-98 opened and closed |
| Devlogs | dl20260822_coulomb_gauge_stepA, dl20260822_coulomb_gauge_stepBC, dl20260823_coulomb_gauge_stepC_elements, dl20260823_hex8_g_and_dr98, dl20260823_tri6_tet10_g, dl20260823_hex8tb_hex8ts_g |

## 3. Element coverage

| Element | G structure | Point-dependent | Key gates |
|---|---|---|---|
| TET4 | link-time constant, s(A⊗B − B⊗A) | no | FD (exact), tie, Frobenius `∫∇w:∇w = ½∫c·c`, gradient mode |
| TRI3 | link-time constant (2D analogue) | no | same battery, 2D rows |
| QUAD4TS | link-time constant, s·F′·∇η⊗∇ξ | no | η-FD from eval columns, fixture-side ∇η, trace 0 |
| PENTA6TS | s[F·W ∓ ½∇τ⊗w], shared W | yes (τ) | FD (exact), tie at τ≠0, corrected kernel (c₀+c₁+c₂ = 0), Ampère negative control |
| HEX8 | term1 + inverse-map Hessian H⁽ᵃ⁾ = −A·Cm·Aᵀ | yes | FD order ratio 4.000 on distorted hex, Stokes (DR-98 gate), pure-H control G·q = H⁽⁰⁾ |
| TRI6 | term1 + Q-folded constant-curvature channel | yes (curved) | curved FD ratio 4.000, NaN-pinned straight skip |
| TET10 | term1 (combine_functions) + Q-folded channel | yes (curved) | curved FD ratio 4.00001, stride/zero-diag breaks |
| HEX8TB | s(∇F)⊗∇ξ on the imposed orthonormal cuboid frame | yes (η,ζ) | FD 4.1e-12, per-dof Stokes (first curl-sign gate), div ≡ 0 |
| HEX8TS | term1 through the class's per-point thin-shell nablas | yes | FD 7.1e-12 on rotated rectangle, tie on the general quad, in-plane kernel |
| LINE3 | **hard-fail stub** (deliberate) | — | contract decision with Christian; recommendation: permanent stub (theory note §8 — the ambient gradient is not defined from a 1D manifold element alone) |

The HEX8TS scope is a documented three-case contract (exact on affine
quads; omits the symmetric Hessian on planar non-parallelograms, where C
is still the exact curl; shares the thin-shell convention on warped
quads). The tie to the implemented C holds on every geometry by
construction.

## 4. Method and audit trail

Every element ran the pipeline C.1 (math, pre-registered on the
exchange) → blind two-vendor audit → C.2 (plan) → blind audit → C.3/C.4
(implementation + battery + deliberate breaks) → blind audit →
reconciliation. All audits were independent Codex + Grok dispatches;
none was skipped or shortened; the final pair folded C.2 into the C.1
round with both auditors' explicit consent. The C.3 Grok audit of the
final pair was dispatched manually by Christian from a self-contained
embedded-source prompt; the report is recorded verbatim in the exchange.

**Headline catches by round** (each would have shipped a defect or a
vacuous gate):

| Catch | Finder | Resolution |
|---|---|---|
| Draft penalty dimensionally inconsistent | Claude derivation, both auditors confirmed | γ = χρ*/μ² form |
| mRhoMin = 1e-16 claim false (actual 0.0) | Grok (L-08 violation) | corrected, self-annotated |
| PENTA6TS "constant field" kernel vector was the Ampère mode | both auditors | correct kernel c₀+c₁+c₂ = 0; wrong vector kept as negative control |
| Transposed-Jacobian wording in plan v2 | Grok | bound to non-transposed J(m,d) with named columns |
| HEX8 linear-φ "positive control" vacuous (isoparametric exactness) | both auditors | replaced by φ=xy-on-box + pure-H gate |
| Q = mJ·e against the member mJ (never written by TET10, NaN on TRI6-straight) | Grok | per-class storage reads + structural skip |
| NaN·0 = NaN poisoning the straight-path Hessian skip | audit round | structural skip + NaN pin as a break gate |
| TS_TestPrism "rectangular" fixture is a general quad | Codex + Grok + Claude prep, independently | rotated-rectangle override; FD/trace moved there |
| HEX8TB Stokes with q=(1,1,0,0) telescopes to 0 = 0 | Grok | per-dof Stokes |
| HEX8TB G-must-move on ξ-shifted points is blind | Grok | ζ-shifted comparison |
| Missing HEX8TS top-edge flip (top-only mS drop invisible) | Codex + Grok | new test; falsified by a targeted break that only it catches |
| "HEX8TS has no production path" | Grok | withdrawn — factory- and Calculator-live |

**Disagreement handling:** no unresolved disagreements remain. Every
refutation above was settled by computation or source citation, never by
majority; the two cases where I disagreed initially (kernel vector,
linear-φ control) I conceded after reproducing the auditors' numbers.

## 5. Verified vs. reviewed (evidence ladder, protocol §11)

**Verified (executed gates with quoted output):**
- All 23 G tests green; `test_fem` 145/145; `make check` 14/14 — run
  after every source edit of the campaign.
- FD agreement numbers: exact-identity fixtures at roundoff (max
  deviations 4.1e-12 – 7.1e-12, step 1e-4); non-affine fixtures by
  order-of-accuracy ratio (4.000 / 4.00001, window [3,5], coarse step
  1e-2 against the truncation floor).
- Fifteen deliberate-break gates across the campaign, each red with its
  pre-registered signature, every restore grep-verified (standing
  procedure after one scripted restore silently failed).
- DR-98 fix: Stokes both-sides-integrated gate pins the absolute curl
  sign at compiled level.

**Reviewed (static analysis, auditor agreement — NOT verified):**
- The theory note's claims about the quench decks themselves (the 1e17
  condition number was never reproduced here — its provenance is an open
  question, L-04).
- The HEX8TS docstring's skew-parallelogram and warped-quad scope
  statements (mathematically audited, not exercised by a fixture).
- The stabilization alternatives (§13): derived and literature-grounded,
  not implemented.
- Recorded low-severity test residuals: TET10 mT ∈ {1,2,3}
  face-orientation cases (inherited, untested for E/C too); HEX8TB
  fixture axis-aligned; HEX8TS kernel node pairs hardcoded.

## 6. Open questions (Christian)

1. **LINE3 contract** — permanent stub (recommended) or a defined 1D
   convention?
2. **1e17 provenance** — which deck, which solver report, and is it the
   raw or equilibrated condition number?
3. **Penalty target regime** — is the goal solver robustness (then the
   weak-divergence route matters) or gauge cleanliness of postprocessed
   fields (then G-diagnostics may already suffice)?
4. **ρ\* averaging fork** for γ = χρ*/μ² (element-local vs. block
   constant).
5. **Quench-deck mesh families** — if hex-dominant decks exist, the
   pointwise penalty is not vacuous there and could be trialed.

## 7. Recommended next step for mt_maxwell_h (recommendation only)

Do **not** wire the pointwise `γ μ² GᵀG` penalty into the simplex /
thin-shell quench decks — the null space it must control is invisible to
it (§1.2), so it would add cost and a tuning knob without touching the
conditioning. The stabilization candidates worth a plan of their own,
in order:

1. **Weak divergence penalty** (patch-level `Dᵀ·M_L⁻¹·D` from E and
   nodal gradients): acts on exactly the offending modes, fits the
   existing K-channel (auto-Δt-scaled, residual-consistent,
   frozen-coefficient → no dJdx), and needs no new element machinery.
2. **Tree–cotree gauging** (Dular et al. 1997): removes the null space
   instead of penalizing it; larger structural change (dof elimination
   in the DofManager).

If a penalty is nevertheless trialed on hex-dominant decks, use the
K-channel form with `γ = χ·ρ*/μ²`, opt-in behind a new `input.conf` key
— which triggers the two-artifact input-contract rule
(`doc/input_file_reference.md` + `doc/input_schema.yaml` in the same
session).

Deferred, needs explicit clearance: the `Calculator::G()` passthrough
(the assembly-side access point; out of the approved file scope of this
campaign).
