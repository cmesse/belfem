# Devlog 2026-04-06 — Nitsche Penalty Rethink

**Date:** 2026-04-06
**Topic:** Read-only review of the ghost/Nitsche penalty for a 1 `mum` HTS layer coupled to 50 `mum` Hastelloy
**AIs involved:** Codex
**Codex Audit Confidence:** high
**Literature References:** Messe et al. 2023 (paper1) Sec. 3, Brenner & Scott 2008 Ch. 10, `todo/thinshell_selective_nitsche_coupling.md`

## Summary

Reviewed whether the current `h_ghost()` penalty is appropriate for a strongly mismatched thin-shell tape stack with `h_HTS = 1e-6 m`, `h_Hast = 50e-6 m`, and `rho_Hast ≈ 123e-8 Ohm m`.

Conclusion: the hard-coded `eta = 4` is probably not the main problem. The dominant effect is the current coefficient model

`alpha = eta * rho_harm * ( 1/hm + 1/hs )`,

with

`rho_harm = rho_m * rho_s / ( rho_m + rho_s )`.

For `rho_HTS << rho_Hast`, this reduces to `alpha ≈ eta * rho_HTS / h_HTS`, so the interface penalty is controlled almost entirely by the HTS side. In that regime, raising `eta` only changes the coupling linearly and does not make the Hastelloy side materially influence the interface stabilization.

Follow-up verification: the cached `element_rho` value is not lagged by a full Newton step or time step. It is written by the block kernels from the current iterate and then read later in the same assembly sweep by the ghost sidesets. Its limitation is that it is an element-mean proxy, not a live interface quadrature coefficient.

## Key Findings

- The ghost kernel uses cached per-element mean resistivities from `element_rho`, not an on-the-fly interface evaluation, in [src/fem/maxwell/matrices/mt_maxwell_h.hpp](/home/christian/codes/belfem/src/fem/maxwell/matrices/mt_maxwell_h.hpp#L27) and [src/fem/maxwell/matrices/mt_maxwell_h.hpp](/home/christian/codes/belfem/src/fem/maxwell/matrices/mt_maxwell_h.hpp#L33).
- `Calculator::q()` reads the current nonlinear iterate directly from the mesh fields in [src/fem/kernel/cl_FEM_Calculator.cpp](/home/christian/codes/belfem/src/fem/kernel/cl_FEM_Calculator.cpp#L2031).
- The conductor block kernels save `element_rho` during block assembly, for example in [src/fem/maxwell/matrices/mt_maxwell_h.cpp](/home/christian/codes/belfem/src/fem/maxwell/matrices/mt_maxwell_h.cpp#L1812), and `DofManager::compute_jacobian_and_rhs()` assembles all blocks before sidesets in [src/fem/kernel/cl_FEM_DofManager.cpp](/home/christian/codes/belfem/src/fem/kernel/cl_FEM_DofManager.cpp#L632) and [src/fem/kernel/cl_FEM_DofManager.cpp](/home/christian/codes/belfem/src/fem/kernel/cl_FEM_DofManager.cpp#L686). So `h_ghost()` sees same-sweep values, not previous-step leftovers.
- `h_ghost()` computes `rho_harm` and `alpha` exactly as `rho_m * rho_s / ( rho_m + rho_s )` and `eta * rho_harm * ( 1/hm + 1/hs )` in [src/fem/maxwell/matrices/mt_maxwell_h.cpp](/home/christian/codes/belfem/src/fem/maxwell/matrices/mt_maxwell_h.cpp#L2078) and [src/fem/maxwell/matrices/mt_maxwell_h.cpp](/home/christian/codes/belfem/src/fem/maxwell/matrices/mt_maxwell_h.cpp#L2087).
- The shell volume operator scales like `rho / h`. This is an inference from `EF_PENTA6TS`: thickness enters the third Jacobian row in [src/fem/interpolation/nedelec/cl_EF_PENTA6TS.cpp](/home/christian/codes/belfem/src/fem/interpolation/nedelec/cl_EF_PENTA6TS.cpp#L132), the inverse Jacobian feeds the curl coefficients, and the element volume is `mThickness * mSurface` in [src/fem/interpolation/nedelec/cl_EF_PENTA6TS.cpp](/home/christian/codes/belfem/src/fem/interpolation/nedelec/cl_EF_PENTA6TS.cpp#L180).
- HTS `rho_piecewise()` is floored at `mRhoMin = 1e-16` in [src/physics/materials/cl_Material.hpp](/home/christian/codes/belfem/src/physics/materials/cl_Material.hpp#L277) and [src/physics/materials/cl_Material.hpp](/home/christian/codes/belfem/src/physics/materials/cl_Material.hpp#L2105), so the ghost penalty can become extremely small but not exactly zero.
- For the stated tape stack and `eta = 4`, `alpha ≈ 4.08e6 * rho_HTS` whenever `rho_HTS << 1.23e-6`. Example: `rho_HTS = 1e-12` gives `alpha ≈ 4.1e-6`, while the Hastelloy volume stiffness scale is roughly `rho_Hast / h_Hast ≈ 2.46e-2`.
- Matching the Hastelloy volume scale would require an effective HTS-side coefficient of roughly `6.0e-9` for `eta = 4`, `2.4e-9` for `eta = 10`, `1.2e-9` for `eta = 20`, or `4.8e-10` for `eta = 50`. This puts the likely size of the `element_rho` proxy error in context: it matters only if the interface-point HTS resistivity is orders of magnitude above the cached mean.
- The internal design note already anticipated this behavior: harmonic weighting deliberately leaves YBCO weakly coupled when its effective resistivity is very small, in [todo/thinshell_selective_nitsche_coupling.md](/home/christian/codes/belfem/todo/thinshell_selective_nitsche_coupling.md#L73) and [todo/thinshell_selective_nitsche_coupling.md](/home/christian/codes/belfem/todo/thinshell_selective_nitsche_coupling.md#L95).

## Changes Made / Proposed

- No source changes made.
- Logged the audit conclusion in `todo/ai_exchange.md`.
- Added this devlog entry.

## Open Questions

- Is the present goal still contrast-robust weak coupling, or do we now want near-strong continuity across the HTS/Hastelloy interface?
- Should `h_ghost()` recover a local interface coefficient on the fly instead of using element-mean `element_rho`?
- If the harmonic coefficient is retained, what `eta` sweep keeps the jump `[\tilde h_t]` acceptably small for the target tape stack?

## Files Updated

- todo/ai_exchange.md
- devlog/dl20260406_nitsche_penalty_rethink.md
- devlog/README.md
