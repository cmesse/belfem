# DR-91 residual H-C: seed-mode initialize() ends the post-restore fixed-dof back-write

**Date:** 2026-08-27
**Purpose:** Close DR-91's parked residual H-C — the fixed-dof dof→field back-write during the
full `DofManager::initialize()` that `load_memdump` runs after the fields are restored
**Register:** DR-91 status cell updated (H-C fixed, stale Anderson clause struck); row stays
unstruck until the executable gate runs
**Thread:** `tmp/ai_exchange/dr91_hc_backwrite.md` (pre-registered plan, two audit rounds,
reconciliations)

## The defect

`Controller::load_memdump` restores the mesh fields from the dump and then runs a full
`initialize()` on any manager not yet initialized — on the in-tree drivers that is the
magnetic manager, since `MaxwellFactory` never initializes and both `hphirun` and `hphiTrun`
load the dump before the first timestep (`hphirun.cpp:69→89`, `hphiTrun.cpp:65→96`). Full
`init_dof_values()` writes **dof → field for fixed dofs**
(`cl_FEM_DofMgr_DofData.cpp:2887-2894`), and at that moment a fixed dof holds its
construction default or factory constant — the per-step BC imposition has not run yet. Every
restored field value at a fixed node was silently overwritten.

## Severity: upgraded from latent to live during the plan audit

The DR-91 row had filed H-C as latent ("bites only decks with thermal Dirichlet nodes").
Grok's plan audit refuted the identity claim with a source trace this session confirmed:

- both drivers `fix()` the abstract-node current dofs at `BELFEM_EPS` **before** the dump
  load (`hphirun.cpp:78-89` → `IWG_Maxwell::set_currents`, `cl_IWG_Maxwell.cpp:106-114`),
  while the dump's `phi` slots hold the dump-time currents in amps (written by the
  post-solve writeback, `cl_FEM_DofMgr_SolverData.cpp:2148-2158`);
- the full-init back-write therefore replaced restored amps with ~machine epsilon;
- `shift_fields` runs at `cl_FEM_Controller.cpp:225`, **before** `compute_boundary_conditions`
  at `:274`, and copies the whole clobbered vector into the BDF history that
  `synchronize_history_fields` had just restored.

Blast radius on accepted physics stays unproven: the DR-92 reader table found no assembly
path that reads abstract-node field slots. The history pollution itself is proven by trace.

## The fix

Seed mode threaded to the branch that already existed for DR-91's original fix:

- `DofManager::initialize( const bool aSeedFreeDofsOnly )` — new non-virtual overload
  carrying the old body; the zero-arg `initialize() override` now forwards to
  `initialize( false )` and **must stay** (dropping it would hide the base name and break
  base-pointer dispatch; no default argument on the virtual).
- `init_dofs( aSeedFreeDofsOnly )` passes the flag to
  `mDofData->init_dof_values( aSeedFreeDofsOnly )`.
- `load_memdump` calls `initialize( true )` at both post-restore sites. For the
  factory-initialized thermal manager that is the same early-return no-op as before.

Cold start is untouched: both vendors independently enumerated every `initialize()` caller
(thermal factory, `IWG_Timestep::shift_fields` via base pointer, the lazy
`compute_jacobian`/`compute_rhs`/RHS-helper sites, dev drivers) — all resolve to full mode.
The seed flag is a compile-time literal on an already-collective path, so it is rank-uniform
by construction. Deliberately out of scope, stated in the call-site comment: fixed dof
*values* stay stale until the first `compute_boundary_conditions`, which the drivers run
before assembly.

## The gate: a focused regression, plus its own detection-power proof

New `tests/fem/test_DofSeeding.cpp` (serial, programmatic two-QUAD4 strip, Dirichlet LINE2
sideset with a master link, `IwgType::Poisson`): dofs are fixed and the field
sentinel-filled **before** the manager's first `initialize()`, mirroring the restore order.

- `SeedModePreservesRestoredFields`: after `initialize( true )` the field is intact at every
  node including the Dirichlet ones, free dofs carry the sentinels, the two fixed dofs keep
  their imposed value.
- `FullModeWritesFixedDofsIntoField`: the negative control — full `initialize()` must
  overwrite the two Dirichlet slots. If this stops firing, the suite has lost the ability
  to see H-C and the first test proves nothing.

The fixture's one trap is recorded in-code: the field must be **sized** before the sentinel
fill (`create_fields()` is public; `Mesh::create_field` allocates unsized, and a length-0
`fill()` is a silent no-op that would fake a pass), and `Mesh( 2, 0, false )` must keep its
third argument — connectivity computation would create edges and the Kernel's MeshChecker
refuses pre-existing edges.

## Gate round 1: failed closed, fixture reworked to TRI3

The first execution of `test_fem` aborted both `DofSeeding` tests in fixture construction:
`Calculator::allocate_memory` (`cl_FEM_Calculator.cpp:1107`) threw "Unsupported Element
Type" — its nedelec-data dispatch covers TRI3/TRI6/TET4/TET10/HEX8 and the thin-shell types
but has **no plain QUAD4 case**, and it runs unconditionally for every block Calculator
regardless of physics. The seed-mode code was never reached; this was the fail-closed
construction risk the code audit had named. The fixture now builds four CCW TRI3 on the same
six-node grid; the Dirichlet facet is local edge 2 of triangle 2 (canonical nodes 2→0).
By-catch for Christian to rule on: a plain QUAD4 volume block cannot allocate a Calculator
under any physics — latent only while no deck meshes a 2D block with quads.

## Process and evidence level

Full round both stages: pre-registered plan → Codex + Grok plan audits (no blocker; Grok's
first attempt burned its 30-turn budget and was discarded — rerun at `GROK_MAX_TURNS=80`) →
implementation → Codex + Grok code audits (no blocker). Syntax gate: `g++ -fsyntax-only`
with each target's own `flags.make`, all three TUs clean.

**GATE PASSED — verified by execution (2026-08-27).** Christian ran `test_fem` after the
TRI3 rework: both `DofSeeding` tests green — seed mode preserves the restored field at every
node including the Dirichlet ones, and the full-mode negative control fires, so the suite
provably sees the defect it guards against. With that, DR-91 carried no residual: the row
was struck and moved to `debt_register_closed.md` (live `[P]` count 18 → 17). The QUAD4
Calculator by-catch stays an open question for Christian's ruling, deliberately not filed
as debt. This round's diff is confined to `cl_FEM_DofManager.hpp/.cpp`, the two
`load_memdump` call sites, the new test, and its CMake registration; the unrelated DR-100
probes in the working tree belong to a concurrent session.
