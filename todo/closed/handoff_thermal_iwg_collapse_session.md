# Session Handoff: Collapse cl_IWG_MaxwellThermal Dispatch + Delete Legacy T_h_* Kernels

**Date:** 2026-07-21
**Purpose:** Self-contained brief for a FRESH Fable session to tidy up
`cl_IWG_MaxwellThermal.{hpp,cpp}` and `mt_thermal_h.{hpp,cpp}` the same way the
Maxwell side was collapsed on 2026-07-21: dispatch tree -> collapsed-kernel pick,
then delete the 13 legacy `T_h_*` bodies. Mirrors the Maxwell-side work recorded
in `devlog/dl20260721_h_kernel_triple_audit.md`.
**Module:** `src/fem/thermal` (+ touch points `src/fem/kernel`)
**Status:** ✅ DONE 2026-07-21 — both mission items landed the same day; every trap
(T1–T7) and step (S1–S7) is ticked. Re-verified in tree 2026-08-09:
`cl_IWG_MaxwellThermal.cpp` is 105 lines and dispatches Conductor/ThinShell →
`T_h_picard` / `T_h_newton` (algorithm-selected) and Ferro/Air/Buffer → `T_phi`,
with no material lookup at link time; `mt_thermal_h.cpp` is 129 lines and contains
exactly `T_h_picard` + `T_h_newton` — all 13 legacy `T_h_*` bodies are gone.
The `T_h_newton` wiring that was still "pending its Codex gate" when this brief was
written is live (`cl_IWG_MaxwellThermal.cpp:74`) and is tracked, together with its
verification run, in `thermal_matrices_cleanup_and_newton_plan.md`.
**Closed 2026-08-09** during the todo currentness sweep.

> **Scope guards:**
> - The thermal **Newton tangent** (`T_h_newton`, ∂M/∂T, ∂K/∂T, ∂f/∂T) is OUT of
>   scope — owned by `todo/thermal_matrices_cleanup_and_newton_plan.md` (its
>   O1-O4 decisions are still pending). This session collapses dispatch and
>   deletes dead bodies; it must leave the current frozen-coefficient thermal
>   behavior (`M + ΔtK`) numerically unchanged.
> - The jc/n field-derivative work
>   (`todo/powerlaw_jc_n_field_derivatives.md`) stays in the original chat.
> - `mt_thermal_phi.cpp` remains out of scope. Buffer/Ferro/Air blocks currently
>   dispatch through `T_phi`; preserve that behavior.

---

## 1. Mission

1. Rewrite `IWG_MaxwellThermal::link_to_group()`
   (`src/fem/thermal/cl_IWG_MaxwellThermal.cpp:68-186`) so only
   `DomainType::Conductor` and `DomainType::ThinShell` dispatch the collapsed
   kernel `T_h_calc` unconditionally — no material switch. Keep
   `DomainType::Ferro`, `DomainType::Air`, and `DomainType::Buffer` on `T_phi`
   (`cl_IWG_MaxwellThermal.cpp:177-181`). This removes defect **D2** of
   `todo/maxwell_kernel_collapse_plan.md`: the Conductor branch currently reads
   `aGroup->material()` at `cl_IWG_MaxwellThermal.cpp:72-73`, and the ThinShell
   branch performs a cross-dofmgr material read at `:123-129`.
2. Delete the 13 legacy producers from `mt_thermal_h.{hpp,cpp}`:
   `T_h`, `T_h_ts`, `T_h_metal`, `T_h_ts_metal`, `T_h_alloy`, `T_h_hts`,
   `T_h_ts_hts`, `T_h_hts_defect`, `T_h_ts_hts_defect`,
   `T_h_hts_piecewise`, `T_h_ts_hts_piecewise`,
   `T_h_hts_defect_piecewise`, and `T_h_ts_hts_defect_piecewise`.
   Current definition starts are `mt_thermal_h.cpp:18,147,346,416,520,578,648,
   754,829,938,1008,1112,1187`; declarations are
   `mt_thermal_h.hpp:20-72`. Keep-set: `T_h_calc`
   (`mt_thermal_h.cpp:1297`) plus any helpers it uses.
3. Sweep for stragglers, update both plans' checkboxes, and write the devlog.
   As of the Codex audit, `rg` finds no code callers of the 13 legacy bodies
   outside `cl_IWG_MaxwellThermal.cpp` dispatch plus their declarations and
   definitions; `nonfree/` has no `T_h_*` hits.

This completes the thermal half of R11+R12 of `maxwell_kernel_collapse_plan.md`
(the Maxwell half was done 2026-07-21: `mt_maxwell_h` went 2724 -> 325 lines).

## 2. Current state (verified 2026-07-21)

- `cl_IWG_MaxwellThermal.cpp` (193 lines): timestep flags are set in the ctor
  (`:25-32`, including `dKdX_times_x`, `dMdX_times_x`, and `dMdX_times_h`).
  The `link_to_group()` preamble is split across `IWG::link_to_group(aGroup)`
  at `:56`, local group/calculator binding at `:58-59`, empty-group early return
  at `:60-63`, and `mTimeStepMatrices->initialize(tNumDofs)` at `:66-67`.
  The dispatch tree is `:68-186`; LookupAlloy already dispatches `T_h_calc` on
  Conductor and ThinShell branches (`:80-84`, `:137-141`).
- `mt_thermal_h.cpp` (1339 lines): 13 legacy producers plus `T_h_calc`
  (`:1297`). `T_h_calc` routes material math through the thermal-side
  `calculator::MaxwellData` instance: cp (`:1316`), lambda (`:1320`), rho
  (`:1323`), `|j|` from the Maxwell peer (`:1326`), and density (`:1330`).
- `T_h_calc` does not consult group material inside thermal `link_to_group`.
  It calls `aCalc->maxwell()` at `mt_thermal_h.cpp:1305`; the helper's material
  was bound earlier from the Maxwell peer block in `MaxwellData::MaxwellData`
  (`cl_FEM_Calculator.cpp:62-68`). Peer relinking for thermal assembly is
  `Calculator::link_element_thermal_maxwell` (`cl_FEM_Calculator.cpp:2531-2547`).
- MaxwellData constructor dispatch covers the legacy material families by
  properties: metal/default/HTS rho routing for thin shell (`cl_FEM_Calculator.cpp:147-208`)
  and bulk (`:210-282`), plus lambda routing (`:285-296`). The legacy
  `T_h_ts` UserDefined reduced-overload branches (`mt_thermal_h.cpp:267-319`)
  are intentionally superseded by the O1 full-signature policy recorded in
  `todo/maxwell_kernel_collapse_plan.md:796-813`; do not promise bitwise parity
  for partial-dependency UserDefined materials without that caveat.
- Buffer behavior is current and must be preserved. `ThermalFactory` selects
  Conductor/ThinShell/Ferro/Buffer blocks (`cl_ThermalFactory.cpp:71-99`) and
  activates Buffer blocks (`:172-184`). `create_buffers()` demotes rho-less
  shell layers to `DomainType::Buffer` (`cl_ThinShellFactory.cpp:1730-1736`).
  Thermal dispatch sends `Ferro`, `Air`, and `Buffer` to `T_phi`
  (`cl_IWG_MaxwellThermal.cpp:177-181`), and `T_phi` assembles M/K only, with no
  Joule source (`mt_thermal_phi.cpp:13-37`).

## 3. The Maxwell-side model to mirror

See `devlog/dl20260721_h_kernel_triple_audit.md` sections "Dispatch collapse +
linking fix" and "Legacy kernel deletion" (`:119-138`, `:167-180`). Maxwell
dispatch now picks by solver algorithm and µ-constancy; all material specifics
live in MaxwellData; the legacy bodies were deleted after a straggler sweep.

For the thermal side the pick is simpler: **every Conductor/ThinShell group ->
`T_h_calc`**. There is no algorithm split until the thermal-Newton plan lands.
Leave one comment marking where a future `T_h_newton` pick will go.

## 4. Traps (hard-won this week — read before editing)

- [x] **T1 — do NOT drop the timestep-matrix preamble.** Keep
  `IWG::link_to_group(aGroup)` (`cl_IWG_MaxwellThermal.cpp:56`), the group and
  calculator bindings (`:58-59`), the empty-group early return (`:60-63`), and
  `mTimeStepMatrices->initialize(tNumDofs)` (`:66-67`). The Maxwell rewrite
  lost the initialize call and every run died in `IWG_Timestep::bdf1` with a
  Blaze "Vector sizes do not match" abort on unsized `f()`.
- [x] **T2 — sign ruling side effect (expected, not a regression):** the field
  reconstruction rule is `h = -B·φ` (gradient), `h = +E·q` (edge), `j = C·q`
  (Christian 2026-07-21, recorded in `maxwell_kernel_collapse_plan.md`
  §6.1). The legacy `T_h_*` bodies carry the old `b = -µ0·E·q` / `bt = -µ0·E·q`
  sign (`mt_thermal_h.cpp:63,224,387,488,619,720,795,901,979,1080,1153,1259`)
  and their copied thin-shell normal projection blocks
  (`:205-208`, `:473-476`, `:705-708`, `:886-889`, `:1065-1068`, `:1244-1247`).
  `T_h_calc` gets fields through MaxwellData (`cl_FEM_Calculator.hpp:2620-2651`),
  which implements the corrected convention. Deleting the legacy bodies closes
  the thermal half of the sign split, but TS normal-metal β-dependent layers may
  shift slightly versus old runs. Do not "fix" this back.
- [x] **T3 — D2 is the reason the material switch must go:** volume-conductor
  thermal dispatch must not dereference the thermal group material
  (`cl_IWG_MaxwellThermal.cpp:72-73`), and the ThinShell dispatch must not read
  the Maxwell material just to choose a thermal producer (`:123-129`). `T_h_calc`
  needs no IWG link-time material lookup; MaxwellData already binds the Maxwell
  peer material during helper construction (`cl_FEM_Calculator.cpp:62-68`).
- [x] **T4 — Buffer/demoted-layer behavior:** do not route Buffer blocks to
  `T_h_calc`. `create_buffers()` demotes rho-less shell layers (for example
  magnesia) to `DomainType::Buffer` (`cl_ThinShellFactory.cpp:1730-1736`).
  The thermal kernel includes those blocks (`cl_ThermalFactory.cpp:71-99`,
  `:172-184`) and current dispatch sends them to `T_phi`
  (`cl_IWG_MaxwellThermal.cpp:177-181`). A rho-less layer has no Joule source
  but still conducts heat through `T_phi`'s M/K assembly (`mt_thermal_phi.cpp:13-37`).
- [x] **T5 — `T_h_calc` is R6-era and still needs a gate:** before deleting the
  legacy bodies, Christian should run one coupled deck (hphiTrun-style) with the
  collapsed dispatch and compare temperature fields / iteration counts against
  the pre-rewrite tree (baseline commit `ce6a0e8b` + current uncommitted state;
  ask Christian for the reference run). Delete only after that gate passes.
- [x] **T6 — straggler sweep scope:** after deletion grep `src/ nonfree/` for
  every deleted name. Expected code result: no hits for deleted producers; docs,
  plans, and devlogs may still mention them historically. Remove the declarations
  from `mt_thermal_h.hpp`; trim includes and dev pragmas with the bodies.
- [x] **T7 — concurrent edits:** Christian often edits the same files while a
  session runs. Re-read regions immediately before editing; if an edit anchor
  fails or a file changed on disk, re-read before proceeding.

## 5. Ordered steps for the new session

- [x] **S1** Read `AGENTS.md`, `doc/ai_collaboration_protocol.md`, this file,
  `todo/maxwell_kernel_collapse_plan.md` (Status + R11/R12 + §6.1-6.2),
  `todo/thermal_matrices_cleanup_and_newton_plan.md` (§1 inventory, scope
  guards), and `devlog/dl20260721_h_kernel_triple_audit.md`. *(2026-07-21)*
- [x] **S2** Verify §2's current-state facts against the tree; these files move
  fast. Reconcile any drift before proposing edits. *(2026-07-21 — no drift;
  dispatch `:68-189`, `T_h_calc` at `mt_thermal_h.cpp:1296-1336`,
  Buffer/Ferro/Air → `T_phi` confirmed preserved)*
- [x] **S3** Propose the collapsed `link_to_group` in chat: preserve T1, map
  Conductor/ThinShell to `T_h_calc`, keep Ferro/Air/Buffer on `T_phi`, and keep
  the default error path. Walk it past traps T1-T4 and have Codex audit the
  proposal. Implement only after Christian approves. *(2026-07-21 — implemented and
  approved. Codex audit: `tmp/ai_exchange/thermal_iwg_collapse.md`. Kernel
  renamed `T_h_calc` → `T_h_picard` after Christian ruled that B/β are frozen
  per thermal solve; `T_h_newton` therefore needs only dcp/dT, dλ/dT, and
  dρ/dT. Same-session Codex finding: `material::Alloy` lacked normB/angleBxJ
  dependency bits despite real ρ(T,B,β)/λ(T,B,β) overrides;
  `populate_rho_database()` now sets those bits, restoring field-aware
  routing on both kernels.)*
- [x] **S4** Christian builds and runs the T5 verification gate. *(2026-07-21
  — passed; Christian: "I think we are looking good".)*
- [x] **S5** Delete the 13 legacy bodies, their hpp declarations, and dead
  helpers; run the T6 straggler sweep. *(2026-07-21 — `mt_thermal_h.cpp`
  1339 → 55 lines, `mt_thermal_h.hpp` 77 → 40 lines; keep-set `T_h_picard`
  only; dead includes dropped (globals, fn_norm, fn_sum, DofManager/Kernel/
  Controller); straggler sweep clean, with zero code references to deleted
  names.)*
- [x] **S6** Bookkeeping: tick the thermal halves of R11/R12 in
  `maxwell_kernel_collapse_plan.md`, note the §6.1 sign-split closure, update
  `thermal_matrices_cleanup_and_newton_plan.md`'s Status (single-sourcing goal
  largely met; Newton scope still open), write
  `devlog/dlYYYYMMDD_thermal_kernel_collapse.md`, and update `devlog/README.md`.
  *(2026-07-21 — all done; devlog = `dl20260721_thermal_kernel_collapse.md`)*
- [x] **S7** Report final line counts and any surprises into the fresh exchange
  thread so the original context can be distilled later. *(2026-07-21 — final
  entry in `tmp/ai_exchange/thermal_iwg_collapse.md`; surprises recorded: D14
  Alloy dependency bits, thermal Buffer/Ferro material gap,
  `T_h_calc` → `T_h_picard`.)*

## 6. Conventions

- Read-only until Christian explicitly approves source edits; Christian runs all
  builds and examples.
- Exchange thread: `tmp/ai_exchange/thermal_iwg_collapse.md` (create fresh).
- Confidence labels (high/medium/low) on non-trivial claims.
- Cite file:line for every factual statement about the tree.
