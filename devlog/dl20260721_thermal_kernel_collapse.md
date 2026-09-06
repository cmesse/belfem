# Devlog 2026-07-21 — Thermal Kernel Collapse (R11 + R12 thermal half)

**Date:** 2026-07-21
**Topic:** Collapse of the `IWG_MaxwellThermal` dispatch tree onto the single
MaxwellData-delegating kernel `T_h_picard`, deletion of the 13 legacy `T_h_*`
producers, plus two adjacent fixes (Alloy dependency bits, thermal
Buffer/Ferro material sharing). Executes
`todo/handoff_thermal_iwg_collapse_session.md`; mirrors the Maxwell-side
collapse of `devlog/dl20260721_h_kernel_triple_audit.md`.
**AIs involved:** Claude (Fable, primary), Codex (audit)
**Claude Confidence:** high (all claims citation-checked in tree)
**Codex Audit Confidence:** high
**Thread:** `tmp/ai_exchange/thermal_iwg_collapse.md`

## Summary

The thermal half of `maxwell_kernel_collapse_plan.md` R11+R12 is done.
`IWG_MaxwellThermal::link_to_group` no longer consults any material at link
time: Conductor and ThinShell domains dispatch `T_h_picard` unconditionally
(renamed from `T_h_calc`), Buffer/Ferro/Air keep `T_phi`. After Christian's
verification gate passed, all 13 legacy `T_h_*` bodies were deleted:
`mt_thermal_h.cpp` 1339 → 55 lines, `mt_thermal_h.hpp` 77 → 40 lines,
`cl_IWG_MaxwellThermal.cpp` 194 → 91 lines. Straggler sweep clean (zero code
references to any deleted name). The §6.1 sign-convention split is closed —
the last legacy `−E·q` mirrors are gone from the tree.

## Key Findings

- **D2 crash class removed:** the legacy Conductor branch read the thermal
  group's material (`cl_IWG_MaxwellThermal.cpp:72`, null for volume-conductor
  thermal blocks); the ThinShell branch read the Maxwell material
  cross-dofmgr (`:123-129`). Both reads are gone; `T_h_picard` needs no
  link-time material (MaxwellData binds the Maxwell-side material at helper
  construction, `cl_FEM_Calculator.cpp:62-68`).
- **Christian's tangent ruling:** B and β depend on h/j, which are not
  thermal dofs — the thermal Newton tangent needs only dcp/dT, dλ/dT, dρ/dT.
  Hence the collapsed kernel IS the Picard kernel (renamed `T_h_picard`,
  matching the Maxwell `h_picard`/`h_newton_*` naming), and `T_h_newton`
  (owned by `thermal_matrices_cleanup_and_newton_plan.md`) is purely the
  T-derivative companion. Recorded in both plans.
- **D14 (new, found by Codex, fixed):** `material::Alloy` — type `PureMetal`,
  genuine ρ(T, log₁₀B, β) database + WF-quotient λ(T,B,β) overrides
  (`cl_Material_Alloy.hpp:141-184`) — never set normB/angleBxJ dependency
  bits (splines register T-only, `cl_Material_SplineLookupTable.cpp:47`), so
  MaxwellData's property-keyed dispatch routed it T-only on BOTH kernels
  since the collapses. Fix: set the rho/lambda T+normB+angleBxJ bits at the
  end of `Alloy::populate_rho_database()` (`cl_Material_Alloy.cpp:484-499`),
  mirroring Copper — the single funnel where the field database becomes
  valid, so database-less decks keep T-only routing. Not corc-relevant
  (hastelloy is `LookupAlloy`, legacy-thermally T-only `T_h_alloy` anyway).
- **Thermal Buffer/Ferro material gap (reported by Christian, fixed):** the
  thermal kernel never receives materials (`ThermalFactory` assigns none;
  `auto_set_materials` fires only for dofmgr index > 0,
  `cl_FEM_DofManager.cpp:1056`), so `T_phi` null-dereferenced on phi-type
  thermal domains (`mt_thermal_phi.cpp:19-21`). Fix:
  `Controller::set_thermal_kernel` now shares the maxwell-side material
  pointer for thermal blocks lacking one (`cl_FEM_Controller.cpp:1818-1828`;
  non-owning per `Group::set_material`, `cl_FEM_Group.cpp:85-88`).
- **Codex audit deltas already covered by closed decisions:** the T upper
  clamp `[gTmin, Tmax]` (O2) and UserDefined full-signature `jc_eval`/
  `n_eval` routing (O1/D4) are accepted semantic changes, not regressions.
- **Expected numerical shift (documented, not a regression):** legacy
  `T_h_*` carried `b = −µ0·E·q`; MaxwellData implements the corrected
  `h = +E·q` ruling — β-dependent metal properties on TS layers may shift
  slightly (partial-flip effect). Deletion closes the sign split.

## Changes Made

- `src/fem/thermal/cl_IWG_MaxwellThermal.cpp` — dispatch collapsed
  (Conductor+ThinShell → `T_h_picard`), dead includes dropped (194 → 91 lines).
- `src/fem/thermal/matrices/mt_thermal_h.cpp` — 13 legacy producers deleted;
  keep-set `T_h_picard` only (1339 → 55 lines).
- `src/fem/thermal/matrices/mt_thermal_h.hpp` — declarations reduced to
  `T_h_picard` (77 → 40 lines).
- `src/physics/materials/cl_Material_Alloy.cpp` — D14 dependency-bit fix.
- `src/fem/kernel/cl_FEM_Controller.cpp` — thermal material-share loop in
  `set_thermal_kernel`.
- `todo/maxwell_kernel_collapse_plan.md` — R11 [x], R12 thermal-half note,
  §6.1 sign-split closure, D14 entry, Status updated.
- `todo/thermal_matrices_cleanup_and_newton_plan.md` — Status rewritten
  (cleanup half done via the collapse; Newton half = remaining scope);
  D1-D6/D9 ticked, R1-R3/O2/O4 struck as superseded.
- `todo/handoff_thermal_iwg_collapse_session.md` — S1-S7 / T1-T7 ticked.

## Open Questions

- `T_h_newton` (dcp/dT, dλ/dT, dρ/dT terms): recipes in
  `thermal_matrices_cleanup_and_newton_plan.md` §3; the dρ/dT MaxwellData
  dispatcher is the missing piece (its O1).
- O3 (volume-HTS anisotropy angle, β = π/2 semantics preserved) and D10
  (stale BDF1 comment `cl_IWG_Timestep.cpp:706`) remain open in the thermal
  plan.
- R12 final gate: full §4.1 matrix on ≥2 MPI ranks for helix + corc.

## Files Updated

- src/fem/thermal/cl_IWG_MaxwellThermal.cpp
- src/fem/thermal/matrices/mt_thermal_h.cpp
- src/fem/thermal/matrices/mt_thermal_h.hpp
- src/physics/materials/cl_Material_Alloy.cpp
- src/fem/kernel/cl_FEM_Controller.cpp
- todo/maxwell_kernel_collapse_plan.md
- todo/thermal_matrices_cleanup_and_newton_plan.md
- todo/handoff_thermal_iwg_collapse_session.md
