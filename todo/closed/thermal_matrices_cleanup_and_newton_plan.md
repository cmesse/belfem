# Thermal Matrix Producers: Cleanup + Consistent Newton Tangent

**Date:** 2026-07-09
**Purpose:** Restructure `src/fem/thermal/matrices/mt_thermal_h.cpp` (13 copy-paste
producers, 1314 lines) into single-sourced formulas, and add the missing temperature
derivatives (∂M/∂T, ∂K/∂T, ∂f/∂T) so the thermal Newton solver runs on a consistent
tangent instead of the current frozen-coefficient `M + ΔtK`.
**Module:** `src/fem/thermal` (+ touch points in `src/fem/iwg`, `src/physics/materials`)
**AIs involved:** Claude (analysis + plan), Codex (audit pending), Christian (decisions)
**Status:** BOTH HALVES IMPLEMENTED — verification gate open (re-checked 2026-08-09).
Drafted + Codex-audited 2026-07-09.
**2026-07-21:** the cleanup goal was delivered by the `maxwell_kernel_collapse_plan.md`
R11/R12 thermal collapse, not by this plan's local R2-R3 path. All 13 producers were
deleted and replaced by the single MaxwellData-delegating `T_h_picard`
(`mt_thermal_h.cpp`, 1339 → 55 lines; dispatch = domain type only), resolving D1-D6/D9
below. **Christian's ruling (2026-07-21):** B and β depend on h/j, which are not thermal
dofs. The thermal Newton tangent therefore needs only dcp/dT, dλ/dT, and dρ/dT.
`T_h_picard` is the Picard kernel; the remaining scope is `T_h_newton` (this plan's §3,
minus any B/β channel). R4-R5 rebase onto `T_h_newton`. MaxwellData already carries
`compute_dcpdT`/`compute_dlambdadT`; the dρ/dT dispatcher is the missing piece (O1).

**2026-08-09 currentness sweep — the Newton half is written and wired; what is left is
proof, not code.** In tree at `df3d8f90`:
- `MaxwellData::compute_drhodT` exists with `_bulk` / `_metal` / `_hts` dispatch
  (`cl_FEM_Calculator.hpp:257,398-413,2605-2900`) — O1 is decided and implemented as
  option (c).
- `T_h_newton` (`mt_thermal_h.cpp:54-127`) fills **all four** Newton blocks and is live in
  dispatch: `IWG_MaxwellThermal::link_to_group` picks `T_h_newton` vs `T_h_picard` on
  `algorithm()` for Conductor/ThinShell (`cl_IWG_MaxwellThermal.cpp:66-77`), and the
  `Btg` scratch matrix is created in `create_custom_vectors_and_matrices`.
- The plan's own two scope guards held: no B/β channel appears in `T_h_newton`, and no
  cross-tangent ∂f/∂h term was added.

So D8 → fixed, R4 → done (by a different route than recommended), R5 → half done.
**Remaining, and it is the load-bearing half:** the global finite-difference tangent check
(R5's verification clause), the `tape_hphiTrun` Newton-iteration comparison (R6), and the
Codex audit of the final diff (R7). Also still open: D7/O3 (volume-HTS β = π/2 physics
call) and the D10 comment fix. `mt_thermal_phi.cpp`'s placeholder properties — called out
as a separate small task in the scope guards below — have since been fixed independently:
`T_phi` now reads `density`/`cp`/`lambda` from the group material.

> **Scope guards:**
> - The magnetic Newton tangent (the `dMdX_times_h` partial-history B7 issue, the
>   timestep-24 plateau) is OUT of scope — tracked separately.
> - Cross-tangent terms ∂f/∂h (thermal source vs magnetic dofs) are OUT of scope: the
>   coupled controller freezes j, B per thermal solve; a block-coupled Newton is a
>   different project.
> - `mt_thermal_phi.cpp` (hardcoded placeholder properties) is a separate small task;
>   noted here only so it is not forgotten.

---

## 1. Current state (verified 2026-07-09)

`mt_thermal_h.cpp` holds 13 producers dispatched from
`IWG_MaxwellThermal::link_to_group` (`cl_IWG_MaxwellThermal.cpp:69-186`) on
(DomainType, MaterialType, have_defect, use_piecewise):

| Volume family | Thin-shell family |
|---|---|
| `T_h` (generic, full runtime dispatch) :18 | `T_h_ts` (generic) :148 |
| `T_h_metal` :360 | `T_h_ts_metal` :433 |
| `T_h_alloy` :540 | — (LookupAlloy maps to volume `T_h_alloy`, see D1) |
| `T_h_hts` :598 | `T_h_ts_hts` :668 |
| `T_h_hts_defect` :772 | `T_h_ts_hts_defect` :847 |
| `T_h_hts_piecewise` :956 | `T_h_ts_hts_piecewise` :1026 |
| `T_h_hts_defect_piecewise` :1130 | `T_h_ts_hts_defect_piecewise` :1205 |

Only **two genuine axes** exist:

- **Geometry / field reconstruction:** volume (`b = −µ0·E·q`) vs thin-shell
  (`bt = −µ0·E·q` plus normal-projected `bn` from the master/slave φ gradients,
  plus the facet normal `n` for the HTS anisotropy angle).
- **Material law:** λ, cp, ρ selection — which `T_h` (volume generic) ALREADY
  dispatches completely at runtime from Material flags (`have(jc)`,
  `have_defect()`, `use_piecewise()`, `depends(rho, normB)`, `T_h:96-136`).

Everything else is duplication: the ~25-line preamble (dofmgr → calculator →
element → material → `link`) appears 13×, the M/K assembly lines 13×, the bn
block 6×, the beta guards in three styles. Confidence: high (direct reading).

**Newton status:** `IWG_MaxwellThermal`'s ctor sets `MatrixFlag::dKdX_times_x`,
`dMdX_times_x`, `dMdX_times_h` with a `//todo :: add derivative terms` comment
(`cl_IWG_MaxwellThermal.cpp:29-32`), but **no thermal producer fills any of
them** → `assemble_dJdx` (`cl_TimestepMatrices.cpp:131-155`) adds zeros and the
thermal "Newton" tangent degenerates to `M + Δt·K` (frozen coefficients). It
converges anyway at 77 K because the nonlinearity is mild there; near quench
(the regime this module exists for) the dominant feedback ∂ρ/∂T is invisible to
it. Confidence: high.

## 2. Defect / oddity ledger (cleanup targets)

- [x] **D1 (LOW)** *(resolved 2026-07-21 — dispatch collapsed to `T_h_picard`)* ThinShell + LookupAlloy dispatches to the *volume* producer
  `T_h_alloy` (`cl_IWG_MaxwellThermal.cpp:136-138`). Benign today (alloy law is
  λ(T), ρ(T) — no b needed) but breaks the family pattern; disappears in the
  unification.
- [x] **D2 (LOW)** *(resolved 2026-07-21 — MaxwellData constructor memoizes density through the generic fallback path)* Density handling inconsistent: `T_h` uses
  `have(density) ? density(gTroom) : ref_density()` (`:91`), the ts/hts
  producers call `density(gTroom)` unconditionally (e.g. `mt_thermal_h.cpp:420`,
  `:527`, `:585`, `:654`, `:758`). Centralize one `effective_density()` helper.
- [x] **D3 (LOW)** *(resolved 2026-07-21 — T/ρ clamps are centralized in MaxwellData compute_T/compute_rho; see collapse-plan O2/O3)* The 2026-07-09 safety clamps (`gTmin` on T, `gRhoMin/gRhoMax`
  on ρ, Christian) are applied producer-by-producer; centralize in the T- and
  ρ-helpers so every producer is covered identically.
- [x] **D4 (LOW)** *(resolved 2026-07-21 — single bj_angle/bn_angle helpers inside Calculator/MaxwellData)* Three beta-guard styles coexist (folded+clamped b·j form,
  eps-guarded+clamped b·j form, clamped n·b HTS form). Keep the two physically
  distinct angles (b–j for magnetoresistance, n–b for HTS anisotropy) as two
  helpers matching the magnetic reference (`mt_maxwell_h.cpp:146-148`).
- [x] **D5 (LOW)** *(resolved 2026-07-21 — MaxwellData compute_hn owns the projection; copies deleted)* The bn average+projection block is copied 6×. One helper,
  mirroring `compute_bn` (`mt_maxwell_h.hpp:125-168`). See O2 for where it lives.
- [x] **D6 (HIGH, HPC rule)** *(resolved 2026-07-21 — the allocating defect producers were deleted; MaxwellData compute_x uses preallocated storage)* Hidden per-integration-point allocation in every
  defect producer: `const Matrix< real > Coords = Matrix< real >(aCalc->N(k)*aCalc->X())`
  (e.g. `T_h:100`, `T_h_ts_hts_defect:943`). Allocates a fresh matrix each
  point, each element, each iteration. Replace with a preallocated work vector
  (calculator work storage).
- [ ] **D7 (MEDIUM, physics)** *(π/2 semantics preserved inside MaxwellData bulk-HTS dispatch; O3 still open)* Volume HTS producers compute `beta` from (b, j)
  but then pass `constant::pi/2` to `rho_powerlaw` (`T_h:103-122`, todo comments
  "add angle dependency here"). In the volume there is no tape normal — decide
  the intended angle definition (O3) instead of silently keeping π/2.
- [x] **D8 (HIGH, the Newton item)** ~~Derivative matrices never filled; `dFdX`
  flag not even set.~~ **FIXED — `T_h_newton` fills all four**
  (`mt_thermal_h.cpp:54-127`): `dMdx_times_x` (ρ·dcp/dT·T), `dMdx_times_h`
  (ρ·dcp/dT·T0 against `collect_qhist`, so the history term is exact at every BDF
  order), `dKdx_times_x` (the unsymmetric mixed-operator block
  (Bᵀ∇T)⊗(dλ/dT·N), carried in the preallocated `Btg` workspace — no per-point
  temporary, cf. D3), and `dfdx` (Nᵀ·dρ/dT·|j|²·N, the quench feedback). The
  `dFdX` flag is set in `cl_IWG_MaxwellThermal.cpp:24`. See §3.
- [x] **D9 (LOW)** *(resolved 2026-07-21 — deleted with the producers)* Dead include `fn_sum.hpp`; stray commented-out code blocks.
- [ ] **D10 (LOW, found by Codex 2026-07-09)** Stale comment on the BDF1
  assembly: `A = 3*M + 2*h*K` where the code assembles `M + h*K`. Comment fix only.
  *Still open 2026-08-09; the line has moved to `cl_IWG_Timestep.cpp:798`.*

## 3. The consistent thermal Newton tangent

Literature: Bathe §7.2.2 (incremental equations, nonlinear heat transfer) and
Table 7.2 (the finite element matrices, including the ∂k/∂θ and ∂c/∂θ
contributions); BELFEM hybrid strategy context in Messe et al. 2023 (paper1)
§2.7. The BDF1 residual in BELFEM's assembled form
(`cl_IWG_Timestep.cpp:697-705`) is

```
R(T) = M(T)·(T − T⁰) + Δt·K(T)·T − Δt·f(T),      A = M + Δt·K
```

so the exact tangent is `J = A + dJdx` with the four correction terms mapped
onto the existing `assemble_dJdx` contract
(`dJdx = α·dMdX_times_x − dMdX_times_h + Δt·dKdX_times_x − Δt·dFdX`,
`cl_TimestepMatrices.cpp:131-155`; α = 1 for BDF1):

| Matrix | Element recipe (per int point k) | Physics |
|---|---|---|
| `dMdX_times_x` | `w·Nᵀ·( ρ_m·cp′(T)·T(ξ) )·N·dV`, `T(ξ) = Nvec·q` | heat-capacity change |
| `dMdX_times_h` | `w·Nᵀ·( ρ_m·cp′(T)·T⁰(ξ) )·N·dV`, `T⁰(ξ) = Nvec·qold` | history correction (together: `dM/dT·(T−T⁰)`) |
| `dKdX_times_x` | `w·( Bᵀ·∇T ) ⊗ ( λ′(T)·Nvec )·dV`, `∇T = B·q` | conductivity change; outer-product structure mirrors the magnetic `Ctj⊗(jᵀC)` Newton term (`mt_maxwell_h.cpp:456-460`) |
| `dFdX` | `w·Nᵀ·( ∂ρ/∂T·‖j‖² )·N·dV` | **the quench feedback**: T↑ → Jc↓ → ρ↑ → more heating. Sign handled by `assemble_dJdx` (−Δt·dFdX). j is frozen (from the magnetic iterate) → no ∂j/∂T term, consistent with the coupling scope guard. |

Same pattern as the validated magnetic Newton producers
(`mt_maxwell_phi.cpp:89-90` contracts `dmudH·H1` / `dmudH·H0` — exactly the
x/h split above). Additional wiring: set `MatrixFlag::dFdX` in the
`IWG_MaxwellThermal` ctor next to the three existing flags.

**Caveats to record in code comments:**
- The x/h split is exact for BDF1 only; BDF2-5 need the β-weighted history
  contract (the known B7 producer issue,
  `todo/closed/bdf_nonlinear_mass_verification.md`). Do NOT assume thermal is BDF1:
  `hphiTrun.cpp:81` sets the thermal kernel to the configured controller
  timestep method, so a `method: bdf2` input inherits the B7 approximation
  here too — say so where the term is assembled.
- The clamps create kinks: for T at the `gTmin` floor or ρ at `gRhoMin/gRhoMax`,
  the true local derivative is zero. The derivative helper must evaluate
  consistently with the clamped property (see O1) or Newton gets a tangent for
  a function it is not actually iterating.

**Material derivative API gap:** `Material` exposes magnetic/current
derivatives — `dmudH` (`cl_Material.hpp:764-765`) and `drho_powerlaw_dJ` /
`drho_piecewise_dJ` (`cl_Material.hpp:669-727`, used by the magnetic Newton
producers) — but no temperature derivatives for `cp`, `lambda`, or `rho`.
The dJ family is the natural pattern to mirror for a dT sibling (O1 option c).
Resolution is O1.

## 4. Cleanup design — two options

**Option A (recommended): collapse to 2 producers + inline helpers.**
Keep `T_h` (volume) and `T_h_ts` (thin-shell) only; both use shared inline
helpers in `mt_thermal_h.hpp`:
`grab_maxwell_context()` (preamble), `compute_bn_ts()` (D5),
`guarded_beta_bj()` / `guarded_beta_nb()` (D4), `effective_density()` (D2),
`clamped_T()` (D3), `thermal_rho()` (the material dispatch `T_h` already
contains, plus clamp), `add_newton_terms()` (§3, written once).
`link_to_group` shrinks to the DomainType switch.
*Cost argument:* the dispatch flags are constant per group, so the per-point
branches are perfectly predicted; property table lookups dominate the loop cost
anyway. If profiling ever disagrees, the helpers make per-group resolved
function pointers a local change.
*Con:* diverges from the many-producer style of `mt_maxwell_h.cpp`.

**Option B: one function template, policies, 13 explicit instantiations.**
`template< typename Geometry, typename RhoLaw > void T_h_impl(...)`;
`link_to_group` unchanged; zero runtime branching; formulas still
single-sourced. *Con:* template producers are foreign to the matrices-file
style and harder to step through in a debugger.

Either option writes the Newton terms of §3 exactly once. Decision → O4.

## 5. Plan steps

- [x] ~~**R1** Freeze a golden baseline: run `tape_hphiTrun` for 2 timesteps,
  archive the residual log + one exodus frame. The cleanup (R2-R3) must
  reproduce it (Picard path bit-comparable, Newton path unchanged until R5).~~
  *(superseded — verification ran per `maxwell_kernel_collapse_plan.md` §4.1 against
  baseline commit `ce6a0e8b`; gate passed 2026-07-21)*
- [x] ~~**R2** Introduce the helpers in `mt_thermal_h.hpp`; migrate only `T_h`
  and `T_h_ts` onto them.~~ *(superseded by the MaxwellData collapse — helpers live in
  `calculator::MaxwellData`, not `mt_thermal_h.hpp`)*
- [x] ~~**R3** (after R2, O4 decided) Collapse the 11 specialized producers onto
  the generic pair (A) or the template (B); shrink `link_to_group`.~~
  *(done 2026-07-21 via `T_h_picard` + domain-type-only dispatch; D1-D6/D9 resolved,
  D7 π/2 semantics preserved, D10 still open)*
- [x] **R4** (after O1) **DONE via MaxwellData** — the dT-derivative helpers are
  `compute_dcpdT` / `compute_dlambdadT` / `compute_drhodT` with the bulk/metal/HTS
  dispatch (`cl_FEM_Calculator.hpp:257,398-413,2605-2900`), i.e. O1 option (c)
  (analytic `Material` API), not the finite-difference option (a) that was
  recommended. **The FD unit-check below has NOT been performed** — it is the
  substance of the still-open R6/R7 gate. Original text: Material dT-derivative helper; unit-check the
  derivative helpers against finite differences for each rho overload actually
  used by the producers (plain metal ρ(T), field-dependent ρ(B,β,T), HTS
  power-law and piecewise, defect variants), plus cp and λ.
- [◐] **R5** (after R3, R4) Fill the four Newton matrices per §3, set the
  `dFdX` flag. **Matrices + flag DONE (see D8).** Tangent verification: finite-difference
  the *global* residual on a small case and compare against the assembled Jacobian
  (the only test that catches sign/contract errors — the B7 history lesson).
  **The FD check is the open half of this box** and is the real remaining risk: a
  wrong sign in `dfdx` or a wrong contraction in `dMdx_times_h` degrades the Newton
  rate silently rather than crashing.
- [ ] **R6** Rerun `tape_hphiTrun`; compare thermal Newton iteration counts
  and check `tolerance switch` behavior; devlog the results.
- [ ] **R7** Codex audit of the final diff; distill exchange into this file.

## 6. Open questions

- [x] **O1 — DECIDED (c), implemented.** The dT-derivatives come from the analytic
  `Material` API through `MaxwellData::compute_drhodT` / `_dcpdT` / `_dlambdadT`, with
  per-type dispatch (`_bulk` / `_metal` / `_hts`). Neither (a) nor (b) was taken. The
  validation that (a) was recommended *for* still has to happen — see R5.
  Original text: Material dT-derivatives: (a) centered finite differences in a
  small helper (2 extra table lookups per point, Newton-only, clamp-consistent
  by construction if it calls the same clamped evaluators), (b) analytic spline
  derivatives from the existing spline storage, (c) new virtual API on
  `Material`. Claude recommends (a) first — it validates the term structure
  cheaply — with (b) as the upgrade for the hot HTS path. Christian to decide.
- [x] ~~**O2**~~ *(moot 2026-07-21 — MaxwellData owns the projection)* Where does the shared bn helper live: thermal-local copy in
  `mt_thermal_h.hpp`, or promote the magnetic `compute_bn`
  (`mt_maxwell_h.hpp:125`) into a shared header both modules include?
- [ ] **O3** Volume-HTS anisotropy angle (D7): keep β = π/2 as the documented
  worst-case/isotropic choice, or define β from the local j direction? Needs a
  physics call, not a code call.
- [x] ~~**O4**~~ *(moot 2026-07-21 — the MaxwellData collapse replaced both cleanup shapes; neither A nor B was chosen)* Cleanup shape: Option A (2 producers + helpers, runtime dispatch)
  vs Option B (template + policies). Claude recommends A.

## 7. Exchange files

- `tmp/ai_exchange/thermal_matrices_cleanup.md` — Codex audit of this plan
  (2026-07-09): §3 math confirmed (recipes vs `assemble_dJdx`, signs, x/h
  split vs the magnetic `H1/H0` pattern, `dFdX` sign) at high confidence;
  citation corrections, D10, and the R2/R3 ordering guards folded in above.
  A second round audits the final diff at R7.
