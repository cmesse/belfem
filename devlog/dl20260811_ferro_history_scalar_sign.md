# Devlog 2026-08-11 — Ferro Nonlinear-Mass Tangent: Signed History Scalar

**Date:** 2026-08-11
**Topic:** The `dMdX_times_h` producer in both Maxwell Newton kernels contracted the
history with a magnitude where the isotropized tangent requires a signed projection —
a Jacobian sign flip at AC zero-crossings. Two-site fix, no residual change.
**AIs involved:** Claude (implementation); prior independent Codex and Grok audits reached
the same corrected expression; jury round on the diff follows this session
**Claude Confidence:** high (derivation + slot-contract cross-check + thermal-producer precedent)
**Literature References:** Messe et al. 2023 §2.7 (Newton strategy for nonlinear HTS/ferro);
Bathe 2016 §8.4 (consistent tangent vs modified Newton — an inexact tangent alters the
rate, never the converged solution); Hairer & Wanner (1996) II.4 (BDF Jacobian structure)
**Verification:** static source trace + syntax gate — all three touched files pass
`g++ -fsyntax-only` under the build tree's own `flags.make` (Armadillo, Debug, gnu++17).
**Reviewed, not verified**: the executable gate (Newton iteration counts, below) is owed
and is Christian's run.

## Summary

`phi_ferro_newton` (`mt_maxwell_phi.cpp`, `DomainType::Ferro`) and the nonlinear-μ
branch of `h_newton_mu` (`mt_maxwell_h.cpp`, Conductor/ThinShell blocks with
non-constant μ) fill the `dMdX_times_h` slot with `dmudh * norm(hist-field)`.
The slot contract (`cl_TimestepMatrices.hpp:80`: `(dM_ik/dx_j) * h_k`) and the
consumer (`assemble_dJdx`: `dJdx -= dMdX_times_h`, `cl_TimestepMatrices.cpp:142`)
require a **signed** contraction; the thermal producer (`mt_thermal_h.cpp`,
`T0 = dot(Nvec, qhist)`) honors that. The two Maxwell producers did not. Both now
compute the signed projection `(Hvec · Hvec_hist) / H`.

## The Mathematics

The BDF residual carries the mass term `M(q) · (α q − q_hist)`, where
`q_hist = Σ βᵢ q^{n−i}` is the β-weighted history (`IWG_Timestep::collect_qhist()`),
and for a ferromagnetic domain

```
M(q) = ∫ Bᵀ μ(H) B dV ,     Hvec = B q ,   H = |Hvec|
```

(φ-kernel: `Hvec = −B φ`; the sign cancels in every expression below, so the kernel
drops it). Differentiating the mass term w.r.t. `q_j`, with `v = α q − q_hist`:

```
(∂M_ik/∂q_j) v_k = ∫ (dμ/dH) (∂H/∂q_j) (Bᵀ B v)_i dV ,   ∂H/∂q_j = (ĥᵀ B)_j ,  ĥ = Hvec/H
```

— a rank-one, non-symmetric tensor `(Bᵀ B v)(ĥᵀ B)`. BELFEM **deliberately
isotropizes** it (accepted rate trade, unchanged by this fix): the directional
structure is replaced by a scalar multiple of `BᵀB`, the scalar being the projection
of the exact term onto the field direction:

```
dμ/dH · ĥ·(B v) = dμ/dH · ( α H − ĥ·Hvec_hist )
```

`assemble_dJdx` splits this as `α·dMdx_times_x − dMdx_times_h`, so the producers must fill

```
dMdx_times_x :  dμ/dH · H                          (was correct)
dMdx_times_h :  dμ/dH · ( ĥ·Hvec_hist )           (was  dμ/dH · |Hvec_hist| )
```

The magnitude equals the projection only when the history field is parallel to and
co-oriented with the current field. Under reversal (`Hvec_hist = −Hvec`, i.e. an AC
zero-crossing between `t^{n-1}` and the current iterate) the correct scalar is
`−dμ/dH·H` while the magnitude gives `+dμ/dH·H` — a Jacobian sign error precisely
where a ferro AC run transitions. The residual is untouched, so by Bathe §8.4 this is
a modified-Newton defect: the residual's roots — and hence the state any converged run
lands on — are unchanged; what can differ is the iteration path (count, robustness near
zero-crossings, and in principle which runs converge at all).

## The H → 0 Decision (both sites, identical)

`|·|` is not differentiable at 0: `ĥ` is undefined and the subdifferential of `H` is
the whole unit ball, so *no* form of the history scalar is canonical there. Decision:
**take the zero (minimal-norm) subgradient** — `h0 = 0` whenever `h ≤ BELFEM_EPSILON`
(`10·machine-ε`, i.e. only an exact numerical zero; the pre-existing
`|dmudh| > BELFEM_EPSILON` guard is untouched). Consequences, and why this choice:

- The `dMdx_times_x` scalar `dμ/dH·H` vanishes with `H` anyway, so at such a point the
  whole nonlinear-mass tangent contribution drops and the Jacobian degrades to the
  Picard matrix — the correct, symmetric limit, and the same operator the Picard
  branch would use.
- For any `h > BELFEM_EPSILON` the division is safe without further clamping:
  Cauchy–Schwarz gives `|dot(Hvec, Hvec_hist)|/H ≤ |Hvec_hist|`, so the scalar is
  bounded by the history magnitude regardless of how small `h` is — no blow-up mode.
- The old code was also arbitrary at this point (it returned `+|Hvec_hist|`, one
  particular subgradient choice with maximal magnitude and unjustified sign); the new
  choice is documented in the source at both sites so they cannot drift apart.

## Changes Made

- `src/fem/maxwell/matrices/mt_maxwell_phi.cpp` — `phi_ferro_newton`: builds
  `hvec = B·φ` and `hvec_hist = B·φ_hist` in scratch, `h0 = dot(hvec, hvec_hist)/h`
  guarded as above. `+ fn_dot.hpp` include.
- `src/fem/maxwell/matrices/mt_maxwell_h.cpp` — `h_newton_mu` ferro branch:
  `hvec = mx->compute_h(k)` (the field `μ(|h|)` and `dμ/dH` are evaluated at — for
  thin-shell blocks this includes the normal part, keeping the projection direction
  consistent with the linearization point), `hvec_hist = E·qhist` in scratch, same
  `h0`.
- `src/fem/maxwell/cl_IWG_Maxwell.cpp` — `create_custom_vectors_and_matrices`:
  registers the two scratch vectors `"hcur"`, `"hhist"` once per calculator (the
  `"vE"`/`"vC"` pattern), sized `mesh()->number_of_dimensions()` — the row count of
  `B`/`E` on every nonlinear-μ path (jury correction: for bulk *and* thin-shell layer
  element types this equals `mesh::dimension(tType)`; the devlog's first draft claimed
  otherwise). No allocation inside the integration-point loop; the φ site in fact
  loses the former per-point `norm(B*phi)` temporaries.

Deliberately **not** done: the exact rank-one tangent (would make the block
non-symmetric — flagged as a possible future opt-in, not worth the solver cost by
default); any change to the `BELFEM_EPSILON` threshold on `dmudh`; any change to
residual assembly, `collect_qhist()`, or the α-scaling.

## Which Run Confirms It

The residual — and so the state any converged run reaches — is unchanged; identical
physics output proves nothing. The gate is **Newton iteration counts at fixed Δt**:

1. **Treatment arm:** an h-φ deck with a ferromagnetic region (`material::Metal` iron,
   nonlinear B-H) under **AC excitation through at least one zero-crossing** of the
   field in the iron — e.g. the bearing/gauge deck, one full period, fixed Δt, Newton
   (not Picard), BDF1 suffices (the defect is order-independent; BDF2 also covers the
   β-weighting). Compare iterations-per-step, HEAD vs pre-fix (`git stash` A/B):
   expect fewer iterations / no stalls in the steps bracketing the zero-crossings,
   and final fields identical to solver tolerance.
2. **Control arm:** the same deck on a **monotonic ramp** (no reversal): the two forms
   nearly coincide, so iteration counts should be essentially unchanged. This is what
   separates "fixed the sign" from "changed something".

Christian runs the decks; no build or solver launch was made from this session.

## Open Questions

- None on the fix itself. The exact (non-symmetric) rank-one tangent remains a
  documented possible follow-up if ferro AC convergence is still unsatisfying after
  this lands.

## Files Updated

- src/fem/maxwell/matrices/mt_maxwell_phi.cpp
- src/fem/maxwell/matrices/mt_maxwell_h.cpp
- src/fem/maxwell/cl_IWG_Maxwell.cpp
- todo/bdf_nonlinear_mass_verification.md, todo/debt_register.md (DR-64)
- tmp/ai_exchange/ferro_history_scalar.md (ephemeral thread; distilled here)
