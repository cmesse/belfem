# Devlog 2026-04-07 — High-Resistivity `h_ghost()` Alpha Review

**Date:** 2026-04-07
**Topic:** Read-only review of the current `h_ghost()` penalty scaling for contact and buffer layers with very large resistivity
**AIs involved:** Codex
**Codex Audit Confidence:** high
**Literature References:** Brenner & Scott 2008, Ch. 10; `todo/thinshell_selective_nitsche_coupling.md`

## Summary

Reviewed the current local-worktree `h_ghost()` penalty scaling after the report that a `1e-3 Ohm m` contact layer stalls near `-40 dB` and a `1 Ohm m` buffer layer produces `NaN`s.

Conclusion: the current penalty is now too large on the high-resistivity side. The active code uses
`alpha = eta * ( rho_m / hm + rho_s / hs )`
with no upper cap, plus a non-negligible absolute floor. For micron-scale layers this makes `alpha` jump from order `1` for normal metals to order `1e3 ... 1e6` for contact / buffer layers, while the consistency terms remain weighted by the harmonic mean and therefore stay near the metallic scale. This creates a penalty-dominated interface operator and is a credible explanation for the observed stagnation and `NaN`s.

## Key Findings

- The active local-worktree `h_ghost()` code computes
  - `rho_harm = rho_m * rho_s / ( rho_m + rho_s )`
  - `alpha = eta * ( rho_m / hm + rho_s / hs )`
  in [src/fem/maxwell/matrices/mt_maxwell_h.cpp](/home/christian/codes/belfem/src/fem/maxwell/matrices/mt_maxwell_h.cpp#L1804) and [src/fem/maxwell/matrices/mt_maxwell_h.cpp](/home/christian/codes/belfem/src/fem/maxwell/matrices/mt_maxwell_h.cpp#L1816).
- The same code then applies an absolute floor
  `alpha = max( alpha, 1e-6 * eta / min( hm, hs ) )`
  in [src/fem/maxwell/matrices/mt_maxwell_h.cpp](/home/christian/codes/belfem/src/fem/maxwell/matrices/mt_maxwell_h.cpp#L1819). For `hm = hs = 1e-6 m` and `eta = 4`, this floor is already `4.0`, so it is not negligible.
- The consistency terms still use `rho_harm` through `rwdS` in [src/fem/maxwell/matrices/mt_maxwell_h.cpp](/home/christian/codes/belfem/src/fem/maxwell/matrices/mt_maxwell_h.cpp#L1876). For a metal plus a very resistive contact/buffer, `rho_harm` collapses back to the metal-side resistivity, so the consistency terms stay small while the penalty term grows without bound.
- Scale check for `hm = hs = 1e-6 m`, `eta = 4`, and a metal side `rho = 1e-7 Ohm m`:
  - metal-metal: `alpha = 0.8`
  - `1e-3 Ohm m` contact plus metal: `alpha ≈ 4.0e3`
  - `1 Ohm m` buffer plus metal: `alpha ≈ 4.0e6`
- `git blame` shows the switch to `alpha = eta * ( km + ks )` and the added floor are local uncommitted changes in the current worktree at [src/fem/maxwell/matrices/mt_maxwell_h.cpp](/home/christian/codes/belfem/src/fem/maxwell/matrices/mt_maxwell_h.cpp#L1812).

## Recommendation

- Do not let the physical high resistivity of contact / buffer layers enter the Nitsche penalty unbounded.
- Keep the actual high resistivity in the volume `K` of the affected layer.
- For the Nitsche penalty only, cap the per-side scales `k = rho / h` before forming `alpha`.
- A practical first test is:
  - `k_cap = 1 ... 10 Ohm`
  - `alpha = eta * ( min( rho_m / hm, k_cap ) + min( rho_s / hs, k_cap ) )`
- Revisit or remove the current absolute floor if it was only meant as a temporary safety net.

## Open Questions

- Whether the reported `1e-3 Ohm m` contact quantity is intended as a true bulk layer resistivity or as an interface/contact law that should eventually be modeled separately from the shell volume resistivity.
- Whether the current floor
  `1e-6 * eta / min( h )`
  should be retained at all once an upper cap on `alpha` is introduced.

## Files Updated

- devlog/dl20260407_hghost_alpha_high_rho_review.md
