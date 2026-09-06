# Ghost and gauge penalties opt-in; the J speckle traced to the buffer interface

**Date:** 2026-09-01
**Purpose:** Session record — from "a concerning amount of noise in the J field at 15 ms" to a
jury-audited change that makes the Nitsche ghost and the Coulomb gauge opt-in, with the
measurements that motivated it.
**Module:** `src/fem/kernel`, `src/fem/maxwell`, `src/mesh`, `doc/input_*`

## What the measurements said (intrinsic-material `tape_quench`, STRUMPACK, cold)

- **The speckle lives on one node sheet.** Sheet-by-sheet `|J|` on the 15 ms frame: rebco
  3.4e10 with 0 % sign flips; stabilizers 1e3–1e4; the **hastelloy|buffer sheet 1e5 median,
  1.3e8 max, 76 % of nodes reversed**. Thinning the hastelloy next to the buffer from 50 µm to
  1 µm (12-layer stack) made the sheet value **~35× larger** at the matched 3 ms frame (max
  3.5e9 vs 1.0e8) while the hastelloy core read a smooth 15–244 A/m². `J = ΔH / h` of a
  tangential-H mismatch across the rho-less 0.15 µm buffer — a buffer-coupling artifact, not
  Nitsche (no ghost facet exists on a buffer interface, `cl_ThinShellFactory.cpp:2132`) and not
  thick-layer recovery. Remedy is in the buffer formulation / postprocessing, still open.
- **κ is the timestep.** κ/Δt = 2.0–2.4e15 per ms over 30× in Δt (1.2e13 at 0.005 ms, 4.5e14 at
  0.23 ms); the mass term sets the small end of the spectrum. At the deck's 1 ms cap κ ≈ 2e15 ≈
  1/ε — the stabilizer speckle in the 15 ms frame is that floor, not the tolerance.
- **The gauge default was invisible.** `chi = 0` vs `1e-4`: κ identical to two digits at five
  steps. `chi = 0.01`: 2–3.5× lower. `eta = 4e-6` vs `4`: same convergence.

## What landed (plan `todo/penalty_opt_in_plan.md`; two jury rounds)

- `src/fem/kernel/fn_FEM_ghost_switch.hpp` — one reader for three consumers: absent block or
  `eta : 0` = ghost OFF, `eta > 0` = ON; a present block without `eta` is an error.
- `ThinShellFactory( …, aCreateGhostFacets )` fed from the deck by `MaxwellFactory`; OFF creates
  neither duplicate interface dofs nor ghost facets (the old `mUseNitsche`, which had survived as a
  hardcoded `true`).
- `Controller::set_params`: gauge absent → `chi 0` (opt-in); ghost via the helper, the 0 written
  collectively so the IWG's constructor default cannot leak into the log; `k_reg` ignored with a
  notice at `eta : 0`.
- The mesh cache tag carries `thinshell.ghost = on|off` (every pre-existing `.bfm` misses once and
  rebuilds); a named `.bfm` built the other way is refused (a pre-change file counts as ON, the
  only layout the old factory produced); a memdump stores edge/face counts and refuses the other
  layout, older dumps caught by a field-length check.
- Tests `tests/fem/test_GhostSwitch.cpp` (10 cases); schema, reference, module docs; the
  `k_reg : 1e-3 Ohm*m` typo in both `tapestack3d` example blocks fixed (`Ohm`).

**Deck sweep (S3):** single-material or same-label stacks (`undulator2d`, `gantry`,
`tapestack2d_gregory`, `block3d`, `garber`, `corc`) never had duplicates and are unaffected; the
multi-metal decks without a block (`tapestack2d_christian`, `circuit`, `tapestack3d_*`,
`tape_quench_usermat`, `corc2` in the build tree) now run the shared-edge model. No block was added
anywhere (Christian: "I'm not afraid of no-ghost").

**Rulings reversed, on the record:** 2026-08-24 "the ghost is load-bearing and must not switch
off by omission" and 2026-08-27 "chi on by default" — both by Christian, 2026-09-01, on the
evidence above.

## D6 resolved (evening): the `.bfm` reload collapses the DG stack

Git: `mCreateGhostFacets = true` in every commit since 08-26, gate unchanged since 06-18 — the
ghost was ON in the code all week. Caches on disk: the 08-30 usermat cache (fresh, ON) holds
526385 edges = 177673 + 14 × 24908 (9 sheets + 5 duplicates) with a ghost sideset; today's
flag-OFF 12-layer cache holds 501477 = 177673 + 13 × 24908, no ghost. Every run that reloaded the
08-30 cache reported 170941 dofs: `reconstruct_edge_connectivity()` relinks elements to edges by
node pair, the duplicates share their node pairs with the originals and come back as orphans. So
the only genuine DG run on record is usermat run 1 (08-30, cold, MUMPS, 295481 dofs, 84 µs steps,
zero `-9`); every later "ghost" run was shared-edge with an inert ghost sideset. Guard landed
(ghost ON ⇒ cache never reused, named `.bfm` refused); the reload fix is owed in the mesh module.

## Open (was D6 — see above)

On this deck the flag was **inert all day**: hardcoded `false` and hardcoded `true` both gave
170941 free dofs; only the 2026-08-30 usermat run 1 ever showed 295481 = 170941 + 5 × 24908, i.e.
five duplicated interfaces of one edge sheet each. Grok refuted the two candidates I had (plugin
`have(rho)`; edge merge) and named the `.bfm` short-circuit (`create_thinshells()` returns early
on a shell-bearing cache; the old tag ignored the switch). Deciding experiment, no code: with the
new binary, ghost ON + `rm tape.bfm` → dof line; same deck with the cache → dof line; `eta : 0`
→ dof line. R8 (the eta 0 vs 4 gate) means nothing until that is done.

## Also this session

MUMPS memory budget (`dl20260901_mumps_memory_budget.md`): the first capped run's step 1 took one
`-9` → `cap ICNTL(23) = 2485 MB` → success; the conditioning instance got its own 1179 MB; RSS 6.8
GB; a `-20` (reception buffer) at step 2 exposed that `-17`/`-20` need the ladder too (D14, fixed).

Exchange: `tmp/ai_exchange/review_penalty_opt_in.md`.
