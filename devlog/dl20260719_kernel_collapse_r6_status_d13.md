# Kernel-Collapse R6 Status Assessment + D13 Discovery and Mitigation

**Date:** 2026-07-19
**Purpose:** Status pass over `todo/maxwell_kernel_collapse_plan.md` (half-checked items),
discovery of defect D13, and the approved one-site dispatch revert.
**Module:** `src/fem/maxwell`, `src/fem/thermal`, `src/fem/kernel`
**AIs involved:** Claude/Fable (assessment + fix); Codex not engaged this session

## What was asked

Continue the kernel-collapse plan, starting with a read-only assessment of the
half-checked defects (D2, D4, R2d, T12) and overall plan drift against tree state
`ce6a0e8b`.

## Findings

**R6 landed silently.** Commit `abb71c2e` (2026-07-15 "backup") flipped all four
LookupAlloy dispatch sites to the collapsed kernels — bulk Conductor and ThinShell
`h_calc` in `cl_IWG_Maxwell.cpp`, Conductor and ThinShell `T_h_calc` in
`cl_IWG_MaxwellThermal.cpp` — unrecorded in the plan and unverified per §4.1 (no
baseline commit, no iteration-count parity run).

**New defect D13 (CRITICAL, now in the plan's §3.2 register):** the TS·LookupAlloy
flip dispatches `maxwell::h_calc` on thin-shell groups, but Maxwell-side thin shells
are FEM **SideSets** (`cl_FEM_DofMgr_SideSetData.cpp:95-170`) and no sideset ever
builds a `calculator::MaxwellData`: both construction gates require
`GroupType::BLOCK` (`cl_FEM_Calculator.cpp:675-682`, `:1146-1151`) and the factory
calls `link_maxwell` only on `dofmgr()->blocks()` (`cl_MaxwellFactory.cpp:736`).
First Jacobian assembly on corc (hastelloy tape layer = LookupAlloy) hits the
`maxwell()` accessor assert (`cl_FEM_Calculator.hpp:3023`).

This — not D1 — is the real cause of the 2026-07-16 double-corc stop recorded in
`dl20260716_periodic_seam_phase_two.md`; the D1 air gates were already in tree and
functional. Correction note added to that devlog. Confidence: high on mechanism;
medium-high that the TS branch was corc's exact crash site (backtrace would settle it).

The thermal-side `T_h_calc` flips are NOT affected: thermal shell layers are FEM
blocks on the thermal dofmgr, and `Controller::set_thermal_kernel` relinks
MaxwellData on both kernels' blocks (`cl_FEM_Controller.cpp:1820-1833`) — verified
this session.

**Half-checked verdicts:**

- D2 [◐] stays — legacy thermal-Conductor dispatch residual still live
  (`cl_IWG_MaxwellThermal.cpp:72` reads null `aGroup->material()`); dies at R12.
- D4 [◐] → **[x] closed** — R2c implementation (`Material::jc_eval`/`n_eval`
  dependency routing, `cl_Material.hpp:301-304` + `powerlaws.hpp`) fully in tree
  and Codex-audited.
- R2d [◐] stays — clamp machinery (`mTClamped`/`mRhoClamped`, zeroed consistent-
  tangent derivatives) verified in `cl_FEM_Calculator.hpp/cpp`, but the
  Controller-side converged-at-clamp diagnostic has zero call sites.
- T12 [◐] → **[x] resolved** — stale note; D8/D9 fixed and the intpoint asserts are
  restored in both coupled dispatchers (`cl_FEM_Calculator.cpp:2462-2500`).

## Change made (Christian-approved)

One dispatch site reverted: `cl_IWG_Maxwell.cpp` TS·LookupAlloy case back to legacy
`h_alloy`/`h_alloy_t` routing (with a D13 pointer comment). Bulk-Conductor `h_calc`
and both thermal `T_h_calc` flips stay. The double-corc run should now clear the
`h_calc` assert point.

## Bookkeeping applied

- `todo/maxwell_kernel_collapse_plan.md`: D13 entry added ([◐] — mitigated, sideset
  design residual); D4 ticked; T12 trap row resolved; R6 marked [◐]
  landed-but-unverified with the revert recorded; Status line refreshed.
- `devlog/dl20260716_periodic_seam_phase_two.md`: D1→D13 attribution corrected
  (addendum, original text left in place).

## Open

- §4.1 verification of the three live R6 flips (Christian runs; needs a baseline
  commit + iteration-count parity per the protocol).
- D13 residual design: MaxwellData construction + material binding for thin-shell
  sidesets (the ctor's `dofmgr()->block( group()->id() )` lookup is meaningless for
  a sideset id) — prerequisite for the R9/R10 TS flips and the TS-alloy re-flip.
- Tails unchanged: R2d Controller diagnostic, D2 residual (dies at R12).
