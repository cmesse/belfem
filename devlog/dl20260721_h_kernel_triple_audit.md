# 3-AI Audit: h_picard / h_newton_mu0 / h_newton_mu + MaxwellData H-first Refactor

**Date:** 2026-07-21
**Topic:** Audit of Christian's uncommitted kernel triple and MaxwellData refactor
(unwired; call-logic audit) for `maxwell_kernel_collapse_plan.md` R5/§6.0
**AIs involved:** Claude (primary), Codex (secondary), Grok (tertiary — first run
lost to narration, retry with explicit final-message rule succeeded)
**Claude Confidence:** high (all merged verdicts citation-checked)
**Codex Audit Confidence:** high
**Grok Audit Confidence:** high (~95%) on F1/F2/F5/HTS-tangent
**Thread:** `tmp/ai_exchange/maxwell_h_kernels_audit.md`

## Summary

Christian implemented three collapsed Maxwell kernels — `h_picard`,
`h_newton_mu0`, `h_newton_mu` (`mt_maxwell_h.cpp:22-130`) — and refactored
`calculator::MaxwellData` to an H-first design (compute_h per domain, µ family
µ0/const/µ(|H|) with dµ/dH, B = µ·H; new cache slots H/normH/mu/dmudh). Not
wired to dispatch. Verdict: **do not merge as-is** — one compile break, one
Jacobian-structure defect, one debug-assert landmine; several smaller fixes.
The HTS tangent and the memoization design are confirmed correct.

## Confirmed defects (3-AI merged, worst first)

1. **Build break + µ0-unit trap (F1, unanimous):** free function `compute_bn`
   renamed `compute_hn` (`cl_FEM_Calculator.hpp:2349`) with 11 legacy call
   sites left (`mt_maxwell_h.cpp:233,371,603,766,925,1103,1275,1454,1633,1812,2354`).
   Old function returned µ0-scaled B-normal, new returns H-normal — mechanical
   rename would silently change units in `b = bt + bn` kernels.
2. **`h_newton_mu` dMdx term structurally deficient (F2, unanimous ~95%):**
   writes only `dMdx_times_h += EᵀE·(dmudh·wdV)`. Under the residual
   `R = (αM + ΔtK)q − M·qhist − Δt·f` and
   `dJdx += α·dMdx_times_x − dMdx_times_h + Δt·dKdx_times_x`
   (`cl_TimestepMatrices.cpp:131-155`), the established isotropized convention
   (`phi_ferro`, `mt_maxwell_phi.cpp:76-91`) needs BOTH blocks with field-
   magnitude contractions: `dMdx_times_x += EᵀE·(dµdH·|H₁|)`,
   `dMdx_times_h += EᵀE·(dµdH·|Hₕ|)` with `Hₕ` from `collect_qhist()`. As
   written: wrong units (µ/H), no companion block, no history contraction.
3. **`MaxwellDataValue::n` never marked current (F6/G1, high):**
   `compute_h_ts_edge/node` call `compute_hn` but drop the legacy
   `set(n, aIndex)` bookkeeping; TS HTS ρ/β paths assert `is_current(n)`
   (`cl_FEM_Calculator.hpp:2796-2802`) → debug assert once wired.
4. **Buffer → `compute_h_ts_node` routing (F4, high hazard):** Buffer blocks
   are nodal-φ (Air doftable, `cl_Maxwell_FieldList.cpp:279-285`); the TS-node
   path runs `get_normal_calculator` facet machinery → assert/null-deref on any
   standalone `buffer { blocks: }` deck (legal grammar; none shipped). Correct
   routing: `compute_h_bulk_node`.
5. **`Ctj` workspace misuse (F5, unanimous):** kernels deep-copy
   (`Matrix<real> Ctj = aCalc->matrix("Ctj")`, `mt_maxwell_h.cpp:59,98`) —
   hidden allocation per element — and the workspace is created `d×1`
   (spatial dim, `cl_IWG_Maxwell.cpp:610`) instead of `e×1` (nedelec dofs).
   Fix: `create_matrix("Ctj", e, 1)` + bind by reference.
6. **Misleading comment (G6, low):** `// per-term mu0 (T2)` on the M-term of
   the µ(H) variant (`mt_maxwell_h.cpp:74,115`).

## Open decision (Christian)

**F3 — edge-field sign:** refactor uses `H = +E·q`; legacy assembly uses
`−E·q` (~25 kernels + thermal mirrors) while the postprocessor already uses
`+E·q` (`cl_MaxwellPostprocessor.cpp:463-466`). The tree has carried a
convention split; the refactor sides with the postprocessor. β is mostly
sign-invariant (`abs(dot)` in bj/bn_angle, `cl_FEM_Calculator.hpp:2250-2266`),
so the impact is signed b/h vectors and non-orthogonal TS mixing — but §4.1
bit-parity vs legacy baselines breaks either way. Decision needed: adopt `+E·q`
globally (kernels + thermal + docs) or restore the legacy minus in MaxwellData.

## Confirmed correct

- HTS consistent tangent `dKdx_times_x += Ctj·(dρ/d|j| / |j|)·Ctjᵀ` — exact
  rank-1 form, chain-rule factor correct, `compute_drhodj` is d/d|J|
  (`powerlaws.hpp:1004-1045`). No K-history block needed (residual has Δt·K·q
  only).
- Memoization coherence incl. the deliberate dmudh-before-mu call order
  (`compute_dmu_material` writes both `mMu`/`mdMudH`, marks µ current).
- `h_picard` and `h_newton_mu0` M/K structure; `norm_j` self-computing guards.
- MatrixFlags already set for all Maxwell IWGs (`cl_IWG_Maxwell.cpp:78-84`) —
  Claude's initial "flags must be enabled" concern was overstated (Grok G4).

## Fixes applied (same day, Christian approved F1/F2/F4/F5/F6)

- **F1:** legacy-compatible `compute_bn` wrapper restored (`bn = mu0 * compute_hn`,
  `cl_FEM_Calculator.hpp` after `compute_hn`) — 11 legacy call sites compile
  unchanged; wrapper retires with them at R12.
- **F2:** `h_newton_mu` now mirrors the `phi_ferro` isotropized convention:
  `collect_qhist()` hoisted before the loop, per-point `H1 = |h(k)|` (memoized
  via MaxwellData) and `Hh = |E·qhist|`, dual blocks
  `dMdx_times_x += EᵀE·(dmudh·H1·wdV)` / `dMdx_times_h += EᵀE·(dmudh·Hh·wdV)`.
  `cl_IWG_Timestep.hpp` include added.
- **F4:** Buffer branch removed from the bulk `mFunH` selection — Buffer blocks
  (nodal φ) take `compute_h_bulk_node`; comment records the F4 rationale.
- **F5:** `Ctj` bound by reference in both Newton kernels; workspace created
  `e×1` (nedelec dofs) instead of `d×1`. Plan §6.1 draft updated to the same
  pattern (per-point `Matrix` construction forbidden — Christian's rule).
- **F6/G1:** `set(MaxwellDataValue::n, aIndex)` restored in `compute_h_ts_edge`
  and `compute_h_ts_node`; stale "compute_bn" assert message fixed; misleading
  "per-term mu0" comments corrected. Bonus: corrupted file header comment in
  `mt_maxwell_h.cpp` restored ("christiaMatrices->dKdx_times_x()" → "christian").

## F3 resolved — sign-convention ruling (Christian)

N, B, E, C are physics-agnostic; the rules are **h = −B·φ** (gradient),
**h = +E·q** (edge interpolation), **j = C·q**. Legacy `−E·q` kernels are wrong.
Claude's cancellation sweep initially concluded the error was silent
everywhere; **Codex's adversarial verification refuted the absolute form**
(counterexample verified by Claude): the invariance holds for bulk kernels
(quadratic/norm/abs-protected — `abs(dot)` masks a FULL b-flip), for TS HTS
kernels (`bn_angle` with bn ∥ n correct-signed and bt ⊥ n), and for
ghost/φ/interface/BC/constraint paths — but **TS normal-metal paths are
conditionally exposed**: `h_ts_metal` (`mt_maxwell_h.cpp:249-268`) and thermal
mirrors feed the PARTIALLY flipped `b = bt(wrong) + bn(correct)` into
`bj_angle(b,j)`, and `abs((−bt+bn)·j) ≠ abs((bt+bn)·j)` whenever `j` has a
nonzero normal component (no local guarantee for TS `C·q`). Practical scope:
β-dependent metal rho (Kohler magnetoresistance) on TS layers. Decision
pending (Christian): document the restriction until R12, or flip
`bt = +µ0·E·q` in the TS normal-metal kernels + thermal mirrors now. Ruling
and refined verdict recorded in the plan §6.1 note; full argument in the
thread (`tmp/ai_exchange/maxwell_h_kernels_audit.md`).

## Dispatch collapse + linking fix (same day, follow-up)

Christian rewrote `IWG_Maxwell::link_to_group()`: the whole Conductor material
tree AND the ThinShell case collapse to one algorithm-based pick over the new
triple (Ferro analogously via new `phi_ferro_picard`/`phi_ferro_newton`,
`mt_maxwell_phi.cpp:50,80`). Claude's review found three crossed ternary slots
dispatching `h_newton_calc` — which hardcodes µ0 in its M-term
(`mt_maxwell_h.cpp:2687`) and has no dmudh block — into (a) both constant-µ
Picard slots (wrong M for µ≠µ0) and (b) the field-µ Newton slot (wrong M and
missing tangent, the one case needing `h_newton_mu`). Fixed (approved):
Newton → `is_constant(mu) ? h_newton_mu0 : h_newton_mu`; Picard → `h_picard`
always; `compute_mu` handles µ0-vs-const internally so the µ0 sub-branch is
gone. Verified sound en route: metal `mFundRhodJ = return_zero`
(`cl_FEM_Calculator.cpp:173`) makes the algorithm-based split correct without
the plan's `have(jc)` distinction; `algorithm()` is live per assembly pass so
the hybrid Picard→Newton stage switch re-dispatches. Consequence: `h_calc`,
`h_newton_calc`, and ALL legacy `h_*` material kernels are unreachable from
dispatch — R6-R10 are superseded by one leap; §4.1 verification now compares
the whole rewrite against pre-rewrite baseline `ce6a0e8b`. `mHaveThermal` in
`link_to_group` is now a dead store (cleanup candidate).

## B/β material-derivative machinery (same day, follow-up 2)

Root-caused the Newton pathology's likely fix path: the tangent misses the
∂ρ/∂(|B|,β) channel (same (n−1) amplification as the modeled |j| channel for
anisotropic Jc). Christian drafted the machinery in Material/Metal/Database;
Claude audited (verdict: right architecture — mirrors the dT pattern) and
fixed, approved:
- **B1 (critical):** `evaluate_derivy/derivz` scaled by `inv_element_step(0)`
  (copy-paste from derivx) → now dims 1/2; would have been wrong by grid-step
  ratios across the (T, log₁₀B, angle) axes.
- **B2:** duplicate `Material::drhodT(T,B,β)` inline removed (compile break).
- **B3:** missing `Metal::drhodB/drhodbeta/dlambdadB/dlambdadbeta` override
  bodies added (link break): Kohler-pointer dispatch + WF-quotient
  (−a·b·dρdX/c², only the denominator of λ·ρ/ρ(B,β) is field-dependent).
- **B4:** base drhodB/drhodbeta fallbacks no longer require a Kohler function
  (field-independent metal → 0 is legitimate); assert strings fixed.
- **B5:** duplicate HEX64 `mdFunction3Ddxi` assignment removed.
- **Codex symbolic sweep** (same method as the 2026-07-11 QUAD9 catch):
  quad16deta, hex64deta/dzeta verified correct; **HEX27 N[23] deta/dzeta
  forms SWAPPED** — confirmed by hand against `eval_hex27:500` (breaks
  ΣdNᵢ=0), fixed both lines + removed the two now-unused locals.
  quad4/quad9/hex8 sets double-confirmed correct.

Still open on this track: analytic `drho_powerlaw_dB/_dβ` (the load-bearing
pair for corc), MaxwellData `compute_drhodb/dbeta` dispatchers, and the
E-channel tangent block in `h_newton_mu0`.

## Legacy kernel deletion (same day, follow-up 3 — R12 Maxwell half)

Convergence validated on corc (timestep 42, Δt=15 ms: Picard −64.7 dB in 2,
Newton monotone to −74.4 dB, 4 iterations / 15 s — no switch regression), then
Christian ordered the cleanup. Deleted from `mt_maxwell_h.{hpp,cpp}`
(2724 → 325 lines): all 24 legacy material kernels (`h_metal*`, `h_alloy*`,
`h_hts*`/`h_ts_hts*` incl. defect/piecewise/_t variants), the generic `h`/`h_ts`,
the R5 `h_calc`/`h_newton_calc`, `get_thermal_calculator`, and the dev pragma
blocks. Keep-set: `h_picard`, `h_newton_mu0`, `h_newton_mu`, `h_ghost`,
`save_resistivity`, `get_resistivity`. The legacy `compute_bn` wrapper in
`cl_FEM_Calculator.hpp` lost its last caller and was removed. Straggler sweep
clean (remaining `T_h_*` hits are the live thermal tree — R11's problem).
Side effect: the §6.1 TS-normal-metal sign question is now confined to the
`mt_thermal_h.cpp` mirrors.

**Open caveat, ATTRIBUTION CORRECTED (Christian 2026-07-21):** the current
material law has CONSTANT jc and n, so the ybco powerlaw tangent is already
exact — Claude's "missing powerlaw B/β channel" hypothesis was misdirected.
The one remaining β/B-dependent resistivity in corc is the METAL/ALLOY layers
(copper/silver Kohler, hastelloy table: ρ(T,|B|,β) with `drhodj = return_zero`)
— their field-derivative stiffness is entirely absent from J, and cryogenic
Kohler magnetoresistance is large. Christian's Metal/Database dB/dβ machinery
(audited above) is exactly this tangent; it still lacks consumers
(`compute_drhodb/dbeta` dispatchers + the E-channel block in `h_newton_mu*`).
Note the clean timestep-42 trace ran at Δt=15 ms (~7× weaker per-step
nonlinearity than the pathological Δt=100 ms traces) — not yet proof the
tangent question is closed. Future item: when jc/n become lookup tables,
`drho_powerlaw_dB/_dbeta` (∂ρ_PL/∂jc · ∂jc/∂B chain + parallel-combination
factor) join the powerlaw family.

## T_h_newton — thermal Newton tangent (same day, follow-up 4)

Christian collapsed the thermal dispatch himself (Conductor/ThinShell →
`T_h_picard` (renamed from `T_h_calc`), Ferro/Air/Buffer → `T_phi`, preamble
intact — overtaking the handoff's Mission item 1) and implemented `T_h_newton`.
Claude audit vs the Codex-audited tangent table
(`thermal_matrices_cleanup_and_newton_plan.md` §3), Codex secondary confirmed
all verdicts by independent derivation:
- **N1 fixed:** `dKdx_times_x` was `BᵀB·(dλ/dT·T)`; exact form is the
  mixed-operator outer product `(Bᵀ∇T) ⊗ (dλ/dT·N)` — implemented with a
  preallocated `Btg` (n×1) workspace via a new
  `IWG_MaxwellThermal::create_custom_vectors_and_matrices` override.
- **N2 fixed:** `dfdx()` was an n×1 column with a spurious `T` factor; now the
  n×n quench-feedback matrix `NᵀN·(dρ/dT·|j|²)`.
- **N3:** `MatrixFlag::dFdX` added to the ctor.
- **N4 fixed:** signed history contraction `dot(Nvec, qhist)` replaces
  `norm(N·qhist)`.
- **Codex extras, all fixed:** `T_newton`→`T_h_newton` declaration mismatch;
  missing `inline` on `compute_drhodT_hts` (ODR risk in the shared header);
  `c==a` singularity guard in the ρ_PL reconstruction `b = a·c/(c−a)`
  (piecewise normal regime; limit correctly returns dρ_n/dT).
- Verified good: dM blocks match plan rows 1-2 incl. `collect_qhist`
  (BDF2-5-correct), clamp zeroing coherent, `compute_rho`-before-`drhodT`
  ordering satisfied, M/K/f mirror `T_h_picard`.

**Wired** (Christian's pass-gate met): thermal Conductor/ThinShell dispatch now
picks `algorithm()==NewtonRaphson ? T_h_newton : T_h_picard`. Open: the
`dbdT = 0` jc(T)/n(T) placeholder in `compute_drhodT_hts` (documented in
`todo/powerlaw_jc_n_field_derivatives.md` O3); run verification (Christian).

## Process note

Grok run 1 returned only progress narration (audit lost); retry with an
explicit "final message must BE the audit" rule produced a full report with
derivations. Worth folding into the Grok invocation wrapper.

## Files Audited (uncommitted diff)

- `src/fem/maxwell/matrices/mt_maxwell_h.cpp/.hpp` (new kernel triple)
- `src/fem/kernel/cl_FEM_Calculator.hpp/.cpp` (MaxwellData H-first refactor)
- `src/fem/maxwell/cl_IWG_Maxwell.cpp` (workspace creation, D13 revert context)
- `src/fem/maxwell/cl_MaxwellFactory.cpp` (ThinShell block selection, prior turn)
