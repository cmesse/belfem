# DR-102 closed — the half of the fix that was recorded but never landed

**Date:** 2026-08-28
**Purpose:** finish DR-102 (order-2 face-key defect) and correct a false closure record
**Module:** `src/mesh` (`cl_FaceFactory.cpp`)

## What prompted this

Christian's claim: *"DR-102: replaced `get_nodes_of_facet` with `get_corner_nodes_of_facet`.
Claim: solved."* The claim was checked against the tree rather than against the register.

## Finding

The corner swap in `face_key_3d` (`cl_FaceFactory.cpp:338`) is real and is the DR-102 fix
proper. It matches the row's description exactly, and `git diff` showed it was the *whole*
working diff of the file.

But both the DR-102 status cell and `tmp/ai_exchange/dr102_face_key_order2.md` recorded the
2026-08-24 landing as **"corner swap in `face_key_3d` + `set_size(3|4)` in `face_key_2d`"**.
The second edit was never applied. `face_key_2d` (`:287-328`) still wrote `aWork(0..2)` /
`aWork(0..3)` into a caller-owned scratch it did not size.

That is a false closure record, not a stale one: the register asserted a specific edit that
the tree contradicted.

## Why it mattered more after the swap than before

The hazard has no live trigger — every caller of `create_faces` passes an empty sideset list
(`cl_MaxwellFactory.cpp:1390-1392` passes an explicit `Vector<id_t>()`; `ProtoMesh`,
`CutFactory`, `cl_Mesh_BfmFile`, `tetsplit` take the empty default), so the `face_key_2d`
half of `count_faces` / `find_face_owners` never executes. Grok was right to call it latent.

Codex was right that the swap makes it worse. The scratch is shared between the two key
functions. Before the swap, `get_nodes_of_facet` left it over-sized (TET10 → 6,
PENTA15/18 triangular facet → 6), so the QUAD branch's `aWork(3)` was in-bounds and merely
midside-contaminated. `get_corner_nodes_of_facet` sizes it to exactly 3 for a triangular
facet, so the same write became a genuine out-of-bounds the moment anyone passed a sideset
list. The corner swap inverted the sign of the latent bug, which is exactly why the plan
round had decided to land both edits together.

## What landed

`face_key_2d` now sizes its own scratch: `set_size( 3, nullptr )` in the TRI branch,
`set_size( 4, nullptr )` in the QUAD branch. Seven lines including the comment. The function
no longer inherits whatever `face_key_3d` last left behind.

## Evidence

- **Reviewed, not verified.** `g++ -fsyntax-only` with the mesh module's own `flags.make`
  (`-Wall -Werror -pedantic-errors -std=gnu++17`) is clean. That is the clean-compile rung
  of the closure ladder, not a suite gate.
- **Surviving gate:** `make check-fast` re-run. It confirms no order-1 regression rather than
  exercising the guard — the guard's own trigger has no reachable caller, so it is untestable
  without one being written.
- **Both edits are still uncommitted.**

DR-102 struck under the DR-42/DR-49 exception (Christian's call): the design and code work is
finished and the residue is a single run.

## The transferable lesson

The 2026-08-24 round did everything right — plan, two blind audits, a genuine divergence,
an explicit decision to land both edits — and then recorded an outcome that did not match
what was written to disk. Neither auditor was wrong; the write-up was. A closure record is
a claim about the tree, and the only thing that settles it is the tree.

---

## Follow-on: the landing-claim sweep

Christian's call after the above: run the same check across the register. One
question only, mechanically — **does the tree contain the edit the status cell says
landed?** Not whether the fix is right, not whether the row should close. That keeps
it separate from a currentness sweep and from a verification pass.

Method: all 38 rows split out, status cells filtered for landing verbs, then each
concrete edit claim traced to disk. Rows whose status names no specific edit
(run verdicts, open rows, decisions) are out of scope and marked so.

### Result — 11 checkable landing claims, 1 false

| row | what the status cell claims landed | tree |
|---|---|---|
| DR-37 | `tests/physics/backendfree/` with directory-scoped `remove_definitions`, `check-fast` + `test_physics` edges; R13's dead `BELFEM_CACHE_LOCATIONS` block deleted and replaced by the contract comment | present — `tests/physics/backendfree/CMakeLists.txt:13,26-32`, `CMakeLists.txt:445-446`, `UserMaterialTemplate.cmake:68-78`; `example_user_material.cpp` pulled from `src/physics/materials/` as claimed |
| DR-75 | multi-group hard error after `get_id_groups`, plus reference §11 row, schema `group_count: 1`, usage-guide grammar paragraph | present — `cl_ElectricalCircuitFactory.cpp:567` (cited `:526`, line drift only), `input_file_reference.md:861`, `input_schema.yaml:1480`, `circuit_usage_guide.md:322-329` |
| DR-85 | `bn = dot(bn,n)*n` inserted before `b = mB + bn` in `compute_superconductor_ts` | present — `cl_MaxwellPostprocessor.cpp:866-869` |
| DR-89 | `absolute tolerance` as a full deck key; synchronize widened to 12 uints / 3 reals; STRUMPACK always applies it, PETSc only when deck-stated | present — `cl_SolverParameters.cpp:70-74,137-138,183-189,320-369`, `strumpacktools.cpp:242` unconditional, `cl_SolverPETSC.cpp:540-541` conditional |
| DR-90 | `mNewtonTrustStreak`, x2 recovery gated on `>= 2`, reset **outside** the improvement branch (the Codex audit correction) | present — `cl_FEM_Controller.cpp:1653-1658` (reset outside), `:1691` gate, all six named reset sites |
| DR-94 | `restart` row no longer claims Δt is restored; "raise this knob cautiously" warning; mirrored into the schema | present — `input_file_reference.md:330,334`, `input_schema.yaml:758` |
| DR-95 | O3 conclusion rewritten to "a ramp started small is sufficient" with the T3/T4-REFUTED note; **other half explicitly NOT applied** | both accurate — `todo/closed/dr92_restart_step_ramp.md:243-253` rewritten; `cl_FEM_Controller.hpp:53-56` comment indeed still un-demoted |
| DR-97 | `tThermalBudgetSpent` in `tThermalReset` under the `tUpdateThermal` guard; §4.3 row, schema note, theory doc §1/§4/§6 | present — `cl_FEM_Controller.cpp:1343,1588,1596-1597`, `input_file_reference.md:318`, `input_schema.yaml:637-645`, `nonlinear_controller_theory.md:22,176,278` |
| DR-100 | orphan `[DIAG]` send + barrier deleted from `recover_fields()`'s early return; two-tags-per-rank-pair contract documented on `comm_tag()`; all probes stripped including the committed Waitall decode | present — `cl_FEM_Postprocessor.cpp:954-970` (one barrier, one send), `commtools.cpp:65-76`, `commtools.hpp:1099` decode gone |
| DR-102 | corner swap in `face_key_3d` **+** `set_size(3\|4)` in `face_key_2d` | **FALSE — the second edit was never applied.** Found and landed this session |
| DR-125 | `seed_dof_values()` on both kernels at the end of `reset_timestep()`, and in `reset_thermal()` | present — `cl_FEM_Controller.cpp:2234,2237` (inside `reset_timestep()`, opens `:2185`), `:2423` |

Out of scope, recorded so the count is honest: **DR-92** ("FIXED AND VERIFIED") names
no edit — it is a run-level verdict on a restart anomaly, and the row already says
"cause narrowed but not proven". **DR-120** reads as a landing claim on a verb scan
("landed together so the magnetic/thermal mirror stays honest") but is future tense
inside an `open` row; checked anyway and the row is honest — `cl_FEM_Controller.cpp:3414`
and `:3590` still assign `get_int( "max iterations" )` into `uint` with no validation.

### What the sweep says about the register

One false landing claim in eleven, and it was the one already in hand. The register's
status cells are, on this evidence, reliable about what reached disk — including where
they record that something did **not** land (DR-95's second half, DR-120), which is the
harder kind of honesty and the kind that was verified here rather than assumed.

The one failure has a shape worth naming. DR-102 was the only row where the landing
was **split across two edits decided in the same breath** — one the fix, one a hazard
guard the auditors disagreed about. Every confirmed row above landed as a single
coherent change, or as a change plus documents that the input-contract rule forces into
the same session. Nothing forces the second half of a *divergence resolution* to be
checked, and that is exactly where the write-up outran the disk.

Line-number drift showed up once (DR-75, `:526` → `:567`) and is expected: the register
already treats `file:line` citations as advisory and anchors on searchable tokens.
