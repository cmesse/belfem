# HEX8TB wall element: oscillation risk and applicability of ghost / Coulomb stabilization

**Date:** 2026-09-05
**Purpose:** Read-only three-voice analysis of the `edge coating : on` wall element (HEX8TB)
as the thin-shell factory builds it and the Maxwell IWG assembles it. Question from Christian:
can the walls oscillate, and should the Nitsche ghost (`eta`, `k_reg`) or the Coulomb gauge
(`chi`) be switched on with them, and at what values.
**Module:** fem/kernel (ThinShellFactory), fem/interpolation/nedelec (EF_HEX8TB), fem/maxwell
**Voices:** Claude (Fable 5.1), Codex gpt-5.6-terra/high, Grok grok-4.6/high — same brief,
blind, exchange slug `hex8tb_wall_stabilization`. No file modified, nothing compiled or run.
**Status:** reviewed, not verified.

## What the wall actually carries (unanimous, confidence high)

The design reference (`src/fem/maxwell/doc/side_coating_wall_element.md` §3) says the uncoated
wall adds *zero new unknowns*: inner traces shared with the shell, outer traces hanging on the
mid-plane curve nodes (the air trace g = φ₀ − φ₁). The shipped path is different:

| slot | entity | unknown? | evidence |
|---|---|---|---|
| inner bottom / top | `SideLayer::InnerEdges` / `InnerEdgeDuplicates`, one source = the shell's side-curve edge (or that edge's sources if it already hangs) | no | `cl_ThinShellFactory.hpp:269-283` |
| outer bottom / top | `OuterEdges` / `OuterEdgeDuplicates`, **no sources** | **yes**, free `edge_h` | `cl_ThinShellFactory.hpp:285-310` |

The mid-plane-node sources for the outer edges exist only inside `if ( aFuseSideEdges )`, fed
by `mFuseEdgesWhenHavingSideConnectors`, hardcoded `false` (`cl_ThinShellFactory.hpp:355`) and
not deck-settable. The `connect_side_edges( SideLayer*, … )` overload that would attach them
afterwards has no caller; the one call at `cl_ThinShellFactory.cpp:372` is the `Layer*` overload
for the shell's own rim. Codex F1/F2, Grok F1/F6, confirmed by Claude from source after first
asserting the design-doc version — **Claude's first answer was wrong on this point and was
retracted in-session.**

This is not a new defect. The 2026-08-13 ruling (Christian + Prof. Sirous, DR-53 struck
2026-08-26) keeps both fuse flags `false` on physics grounds: the fuse overconstrains the rim,
and the wall-side fuse collapsed Δt to the 100 ns floor on `sidecoatings` at t = 2.0825 ms while
free rims pass the cliff warm and cold. The free outer traces are therefore the *validated*
layout; what is stale is the design document's §3 sentence "nothing at the fold is free".

## Q1 — can the walls oscillate?

- **No singular or checkerboard mode from the wall itself.** The wall kernel assembles
  M = μ₀ EᵀE + K = ρ CᵀC (`mt_maxwell_h.cpp:340-375`). M is a positive definite Gram on four
  linearly independent bilinear traces; K's only null vector is the common mode (all four equal),
  pinned by the inner traces. So the free outer dofs are well-posed. Unanimous.
- **What the free outers change:** the outer face satisfies the wall's natural condition instead
  of tangential-H continuity with air; with the ghost on, the two outer twins at a duplicated
  interface may differ freely (Grok F4). A modeling deviation, not an instability.
- **Grok F2** (outer edges of consecutive stations are not shared, so a station-alternating outer
  pattern is held only by each wall's copper stiffness ρ·d/(L·w)) is correct as a statement about
  weak holding. Claude rates it as noise on the wall J_n output rather than a growing mode: each
  such pattern costs energy in every wall independently.
- **Grok F7** (L/R of the outer subspace ≈ μ₀w²/ρ ~ 1 μs, "can ring in time"): the scaling is
  right, the conclusion is rated **low** by Claude. A mode with τ ≪ Δt under BDF1/BDF2 is in the
  stiffly damped regime; trapezoidal rule would ring, BDF does not. Not measured.
- **Scales** (order of magnitude, w = d default, L = 1 mm): wall through-layer stiffness
  ρ·w/(L·d) ≈ ρ/L, i.e. ~2e-6 Ω for Cu at 77 K; layer through-thickness ρ/d ~ 1e-3 Ω; ghost
  α = eta·k ~ 4e-3…4e-2 Ω. The wall is a ~w/L = 1e-3 perturbation on existing diagonals. The two
  curl factors 4/(L·d) and 4/(L·w) are equal at w = d; a wider coating anisotropizes the two
  cross-flow channels by (w/d)². Not a conditioning hazard (Grok F7, Codex Q1, Claude agree).
- **ρ contrast:** walls are pure metal by the factory gate, beside a YBCO side trace whose own
  stiffness is ~0. The copper wall is then the only stiffness on the YBCO rim's bottom-top
  difference and the only path around the φ-only buffer. Physically intended; a Newton stressor
  at most.
- **Still unverified:** the wall reference lists the handedness/orientation gates on both
  connector signs as pending. An edge-direction slip would show as a station-alternating sign in
  the wall's native J = C·q. That gate, not a penalty, is the thing owed.

## Q2 — ghost and Coulomb on the walls

- **Ghost: cannot touch the wall, and need not.** `h_ghost` runs only on `DomainType::Ghost`
  facets between PENTA6TS layers (`cl_IWG_Maxwell.cpp:426-429`); coating blocks never become
  Ghost. With `eta > 0` the wall picks the correct twins — level-j duplicate at the bottom,
  level-(j+1) original at the top (`cl_ThinShellFactory.cpp:572-580`), the same convention the
  layer element uses (`:1958-1959`) — so the wall spans one block and sits in series with the
  ghost interface above it. It neither bypasses nor short-circuits the Nitsche path (all three).
- **Coulomb: not vacuous, but not a gauge.** The wall kernels do assemble χ·ρ·GᵀG. With
  h = h_t(η,ζ)·t, t constant and ∇h_t ⊥ t: |∇h|² = |∇h_t|² = |∇h_t × t|² = |curl h|² pointwise,
  so GᵀG ≡ CᵀC on HEX8TB and the term is a (1 + χ) rescale of the wall's r′. It leaves the free
  outer traces exactly as free as before (all three, algebraic; no Frobenius gate compiled).
  The HEX8 ∇(xy) kernel argument of the gauge note does not carry over: HEX8TB has no
  binormal/normal edges (Grok F3).
- **Doc findings, not fixed (read-only session):** `coulomb_gauge_penalty_theory.md` §10 says
  the penalty matrix is identically zero on HEX8TB — true for a divergence penalty, false for
  the implemented full-gradient form, which equals the wall stiffness. `side_coating_wall_element.md`
  §3/§3.2 describe the fuse ("outer entities own no dofs") as the shipped state. Grok F8:
  `h_side_connector_newton` carries no ρ(B) tangent for the gauge term while `add_rho_field_tangent`
  covers CᵀC — consistent inconsistency, negligible at any gauge-sized χ.

## Q3 — parameters

None. Keep both blocks absent (`chi` off, no `nitsche ghost penalty` block). If a deck already
uses the ghost for the stack interfaces, `eta : 4` / `k_reg : 1e-3 Ohm` stay for that reason,
not for the walls; raising them cannot reach an outer wall trace. Monitor instead:

1. per-station outer-minus-inner wall trace and the cross-section saddle
   (ψ_in,bot + ψ_out,top) − (ψ_out,bot + ψ_in,top), from the native wall J, never the copied seam H;
2. coating on vs off on the same deck with both penalties absent: Newton iterations, timestep
   cuts, side-curve J/Jc (the old wrap's signature was suppressed rim current);
3. ghost on vs off with coating on: current around the buffer must persist either way;
4. the pending orientation gates and the single-tape analytic-r′ comparison.

Do not use κ(A) as the figure of merit (gauge note §5.3, §11).

## Exchange

`tmp/ai_exchange/hex8tb_wall_stabilization.md` — Codex and Grok entries, ephemeral; this devlog
is the distillation. Grok's tool allowlist was read-only; git status unchanged after the round.
