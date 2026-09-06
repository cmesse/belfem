# Devlog 2026-07-13 — Maxwell/Thermal Matrix-Kernel Collapse: Implementation Plan

**Date:** 2026-07-13
**Topic:** Read-only planning session for collapsing the 38 `h_*`/`T_h_*` assembly variants
onto `calculator::MaxwellData`; deliverable is `todo/maxwell_kernel_collapse_plan.md`.
**AIs involved:** Claude/Fable (exploration + plan); Codex (audit round 1, same day); Grok pending
**Claude Confidence:** high (inventory, dispatch shape, trap verdicts — every citation re-read
this session); medium where marked (D1 crash timing, D3 reachability, T9 consumer coverage)
**Codex Audit Confidence:** high — **D1-D5 all confirmed, none refuted** (round 1, 2026-07-13).
Refinements folded into the plan in place: D1 is air-specific (Buffer requires a material,
`cl_FEM_Domain.cpp:36-45`); D2 gains a pre-existing suspect (thermal Conductor dispatch reads
`aGroup->material()`, `cl_IWG_MaxwellThermal.cpp:70-74`); D4 confirmed for the piecewise
overload family too; R3 shadow harness needs an explicit `maxwell()->reset()` between passes;
D3 reachability bounded (no built-in material combines jc with field-dependent lambda).
**Literature References:** Messe et al. 2023 (paper1) §2.5 (equation-object dispatch), §2.6
(material database / Rhyner power law), Eq. 10-13 (Picard→Quasi-Newton, ε<10⁻¹¹,
checkerboarding); Arsenault et al. 2023 (paper3) Section II (magnetodynamic h-φ/T coupling)

## Summary

Produced the required plan (per `todo/plan_template.md`) for replacing the 25 `maxwell::h_*`
variants (`mt_maxwell_h.hpp:40-116`) and 13 `T_h_*` variants (`mt_thermal_h.hpp:21-61`) with
one generic Maxwell kernel + one thermal kernel delegating all material math to
`aCalc->maxwell()`, plus the dispatch-tree shrink (`cl_IWG_Maxwell.cpp:306-397,437-527`;
`cl_IWG_MaxwellThermal.cpp:70-174`). No source touched. Registered in `todo/README.md`.

## Key Findings

**New defects in the current tree (D1-D5 in the plan, found this session):**
- **D1 CRITICAL** — `link_maxwell` is called for every magnetic block with no material/domain
  gate (`cl_MaxwellFactory.cpp:730-737`); air/buffer blocks have null material
  (`cl_MaxwellFactory.cpp:2368-2377`, `cl_FEM_Group.hpp:97`) and the MaxwellData ctor
  dereferences it (`cl_FEM_Calculator.cpp:67,80,85`) → null-deref at first
  `DofManager::initialize` for any h-φ model with air. (Tree flagged "pending build" in the
  helper exchange — plausibly not yet exercised.)
- **D2 CRITICAL** — thermal-side MaxwellData binds `aCalculator->group()->material()`, but
  `ThermalFactory` never assigns materials and `auto_set_materials` is intra-kernel only
  (`cl_FEM_DofManager.cpp:1053-1073`); legacy thermal kernels read the **Maxwell**
  calculator's material (`mt_thermal_h.cpp:32`).
- **D3 HIGH** — one shared β cache slot (`MaxwellDataValue::beta`) is written with the
  bj-angle by `compute_lambda_metal`/`compute_rho_metal` and with the bn-angle by
  `compute_rho_*_ts` (`cl_FEM_Calculator.hpp:2434-2447` vs `:2628-2645`); bulk-HTS relies on
  the reset π/2. First-writer-wins cross-contamination for `have(jc) && depends(lambda,normB)`.
- **D4 HIGH** — UserDefined reduced-dependency overload parity does NOT hold: the helper
  always calls full overloads; the 8-arg defect `rho_powerlaw` asserts full (normB,angleNxB,T)
  jc-dependence (`powerlaws.hpp:229-238`); the 4-arg evaluates ρₙ=rho(T) (`:211`) where the
  legacy 3-arg path uses rho(gTbulk) (`:156`).
- **D5 HIGH** — nothing links the peer kernel's calculator to the paired element after the
  refactor (legacy contract: `get_thermal_calculator`, `mt_maxwell_h.cpp:2474-2492`);
  `Calculator::link(Element*)` resets only its own helper (`cl_FEM_Calculator.cpp:1244`).
  Plan adds `MaxwellData::link_peer()`.

**Brief corrections (scouting facts adjusted after verification):**
- The generic-UserDefined dependency branches exist in `T_h_ts` (`mt_thermal_h.cpp:267-319`)
  only — bulk `T_h` has none.
- Legacy 2D defect `Coords(0,2)` is an **out-of-bounds read** (mX is n×dim,
  `cl_FEM_Calculator.cpp:651-657`) — the helper's mZ=0 (F7) is a bug fix, not drift.
- The vector-aliasing consumer set is narrower than briefed: outside the two matrix files,
  only `MaxwellPostprocessor` touches "b"/"bn" and it overwrites them itself
  (`cl_MaxwellPostprocessor.cpp:637-656`).

**Equivalence verdicts:** all four β conventions reproduce exactly (conditional on D3);
b-sign/mu0 exact; the mu-const/mu(h) b-variants are currently unreachable generalizations
(base `Material` ctor sets constant mu=mu0, `cl_Material.cpp:161`); Newton term FP-identical
via the `drho != 0 && norm_j > eps` guard; `h_alloy`'s skip-b/j optimization falls out of
helper laziness; thin-shell-alloy-uses-bulk-kernel explained (rho is T-only). Divergences
routed to O-items: gTmin/rho clamp ownership (O2/O3), density accessor (O4), MaterialType vs
have()/depends() reconcile (O1), postproc scope (O8, recommended out).

## Changes Made / Proposed

- Created `todo/maxwell_kernel_collapse_plan.md` (38-row inventory, T1-T16 trap verdicts,
  D1-D7 register, R1-R12 steps with shadow-compare verification, O1-O10, target kernel
  bodies, before→after counts). Registered in `todo/README.md`.
- Opened audit thread `tmp/ai_exchange/maxwell_kernel_collapse.md` (D1/D2 refutation is the
  top audit ask; also the shadow-harness legacy-last ordering).
- Cross-links: supersedes the R2-R3 collapse mechanism of
  `thermal_matrices_cleanup_and_newton_plan.md` (Newton work rebases onto the new `T_h`);
  hosts the `rho_lambda_argument_convention.md` decision point (O10).
- Planning phase was read-only. Later same day, with Christian's explicit approval, Claude
  made three source edits to `cl_FEM_Calculator.{hpp,cpp}`: **D11** shared
  `select_link_element_dispatcher()` called from `allocate()` and `link_maxwell`'s rebuild
  branch (fixes stale-memoization on eagerly-initialized thermal calculators); **D10** dead
  `reset_data()` declaration removed; **D5 tail** coupled dispatchers folded to a
  null-or-stale single condition with cached peer + legacy intpoint-match assert restored.
  Not built (user runs builds). D4 proposal ("full-signature policy") filed in plan O1.

## Decision Round (Christian, 2026-07-13, folded into the plan same day)

- **Newton split adopted (plan §6.0):** two versions of each kernel (`h`/`h_newton`,
  `T_h`/`T_h_newton`, dispatch on `have(jc)`); `T_h_newton` is new physics — the
  never-implemented thermal Newton contribution, now buildable with the material API —
  and becomes the rebase target for `thermal_matrices_cleanup_and_newton_plan.md` R4-R5.
- **D1 [◐]** fixed in `allocate()` (Air gate, `cl_FEM_Calculator.cpp:1062`, verified);
  residual: `link_maxwell` rebuild branch (`:1092-1098`) still ungated → R2.
- **D2 clarified:** thermal materials flow from the **Maxwell side** everywhere
  (`mt_thermal_h.cpp:32`; `cl_IWG_MaxwellThermal.cpp:121-128`); thermal groups carry none;
  the Conductor dispatch (`:70-74`) is unexercised (Tape_Quench is thin-shell) and would
  null-deref on a volume-conductor thermal model today.
- **D3 [x]** resolved by policy: (b,n) for HTS/REBCO, (b,j) for metals, never combined in
  one material; action reduced to a convention comment + loud ctor guard.
- **D4 →** policy framing: API unpublished, exact legacy parity not required; documented
  jc-dependency policy + working examples are the acceptance test (O1 reframed).
- **D5** upheld after clarification: the gap is the *peer* kernel's calculator, not the
  assembling one (legacy cross-links per element, `mt_thermal_h.cpp:34`).
- **D6 [x]** self-assign kept deliberately (comment only); **D7 [x]** non-issue
  (once-per-block resize).

## Fix Round (Christian, 2026-07-13 evening, verified by Claude in tree)

- **D1 [x]** — both Air gates in (`cl_FEM_Calculator.cpp:1062` allocate, `:1083` link_maxwell).
- **D2 [◐]** — MaxwellData ctor now binds the Maxwell-kernel block material (`:67`);
  legacy residual: thermal Conductor dispatch (`cl_IWG_MaxwellThermal.cpp:70-74`) still
  reads the null thermal-group material (unexercised; dies at R12).
- **D5 [◐]** — implemented as `mFunLinkElement`/`link_element_maxwell`
  (`cl_FEM_Calculator.cpp:1246,2336-2360`). Claude double-checked the recursion on
  Christian's ask: **no infinite loop** (mElement set before the cross-link; peer's ID guard
  terminates at depth 2). Two new defects found in the process:
  - **D8 CRITICAL** — `mIwgType` never set on block calculators (only the SideSet path calls
    `link(Group*)`, `cl_FEM_SideSet.cpp:88`) → Maxwell-side calcs silently never link the
    thermal peer → stale T in every Maxwell `_t` assembly. Fix: identity-based peer binding.
  - **D9 HIGH** — `peer->element()->id()` null-derefs on the first coupled link
    (`mElement` inits nullptr, accessor unchecked).
  - Also: intpoint-match assert not carried over (re-add); D10 LOW dead `reset_data()` decl.
- **D8 [x] / D9 [x]** — fixed same evening by Christian's three-way dispatcher split
  (`link_element_maxwell` / `_maxwell_thermal` / `_thermal_maxwell`, selected in `allocate()`
  with a `BELFEM_ERROR` fallback, `cl_FEM_Calculator.cpp:1062-1086`) plus null-first peer
  guards. Claude verified: allocate-time `iwg()->type()` is valid; per-element null check is
  cost-noise.
- **D11 CRITICAL (new)** — `link_maxwell`'s rebuild branch doesn't re-select
  `mFunLinkElement`; the thermal kernel initializes eagerly in its factory
  (`cl_ThermalFactory.cpp:242`) before `set_thermal_kernel` (`hphiTrun.cpp:83`), so thermal
  calcs keep `link_element_default` → `MaxwellData::reset()` never fires between elements →
  stale k-indexed memoized values served silently. Fix: shared
  `select_link_element_dispatcher()` called from both `allocate()` and the rebuild branch.

## Decision Round 2 (Christian, 2026-07-14) — all O-items closed

- **O2** confirmed (unified clamp in `compute_T` + `mTClamped` + zero dT-derivatives,
  consistent-tangent analysis grounded in Belytschko §5.9 / Z&T Vol 2 §4.4).
- **O3** resolved with Christian's simpler counter-design: clamp inside `compute_rho` +
  `mRhoClamped`, both flags in MaxwellData, one rho for all consumers (K-path neutrality:
  `Material::mRhoMin = 1e-16` == gRhoMin; `compute_drhodj` = 0 while clamped — expected
  dKdx shadow deviation on deep-subcritical elements, documented).
- **O4** density stays at gTroom — undeformed-mesh physics requirement (transport
  properties are expansion-corrected); doc item → `materials_usage_guide.md`.
- **O5** element_rho mean stays a mesh field, kernel accumulates (unchanged mechanism).
- **O6** confirmed obsolete. **O7** resolved by the `link_element_*` implementation.
- **O8** direction: postprocessors mirror the main kernels' matrix computation; deferred
  past R12. **O9** soft cutover confirmed.
- **O10** resolved AND executed by a parallel Opus session (powerlaw T-bearing overloads →
  `(normJ, T, normB, angleNxB[,x,y,z,t])`; Codex+Grok 0 defects; Claude verified in tree
  before ticking — line counts of both mt_*_h.cpp unchanged, so plan citations hold; §3
  reconciliation note added for the old-order quotes and the ~26-line hpp drift).
- Memory + `todo/README.md` updated for the completed rho/lambda/powerlaw unification
  (doc-pass tail remains in `rho_lambda_argument_convention.md`).

## Open Questions

- ~~O1-O10~~ — all closed by 2026-07-14 (see Decision Rounds above).
- R1 build + smoke run of the wiring fixes — handed to Christian; gates R2.
- The unexercised thermal Conductor dispatch (`cl_IWG_MaxwellThermal.cpp:70-74`, D2
  residual) — dies at R12; would null-deref on a volume-conductor thermal model until then.
- Grok third voice on D1-D11 (optional); R2 helper items then R3 shadow harness.

## Files Updated

- src/fem/kernel/cl_FEM_Calculator.hpp (D10; selector declaration)
- src/fem/kernel/cl_FEM_Calculator.cpp (D11 shared selector; D5 dispatcher fold + assert)
- todo/maxwell_kernel_collapse_plan.md (new)
- todo/README.md (registered)
- tmp/ai_exchange/maxwell_kernel_collapse.md (new, ephemeral)
- devlog/dl20260713_maxwell_kernel_collapse_plan.md (this file)
