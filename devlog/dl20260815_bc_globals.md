# Boundary-Condition Globals Return to the Mesh

**Date:** 2026-08-15
**Purpose:** Record the campaign that restored per-boundary-condition mesh
global variables (the applied current visible in ParaView, as the legacy
code had), the two pre-existing bugs it uncovered in the globals machinery,
and the audit rounds that reshaped the design twice
**Module:** fem/kernel, fem/maxwell, fem/thermal, circuit, mesh
**Round:** `tmp/ai_exchange/bc_globals_plan.md` + `bc_globals_diff.patch`

## The regression

The legacy branch registered one mesh global per boundary condition
(created at link time, written in `compute(aTime)`), and ParaView displayed
the applied current alongside the fields. The current tree writes NO globals:
the Exodus writer machinery is fully intact, but nothing on the BC path calls
`Mesh::create_global_variable` anymore. Verified on the running study — the
frames carry zero global variables.

## Two pre-existing bugs found by the plan review (Codex), fixed first

Both were dead only because no globals existed, and either would have fired
the moment one did:

1. `Mesh::create_global_variable` inserted every auto-ID variable into the
   ID map under the ARGUMENT `aID` (0), not the assigned ID —
   `DofManager::distribute_globals` would abort on its first lookup of a
   real ID under MPI. Now keyed by `tVariable->id()`.
2. `DofManager::create_globals` called the creator with `(label, id, value)`
   against a `(label, value, id)` signature — compiling because `id_t` and
   `real` interconvert. Fixed.

Plus one latent from the implementation audit: `Mesh::load_globals` now
bumps the auto-ID counter to the highest loaded ID, so a future post-load
creation cannot reuse a dump-occupied ID.

## What shipped

- `PhysicalBoundaryCondition` gains `mLabel` and `update_global()`. The
  helper writes `mValue` into the mesh global; it is called in `impose_bc()`
  BEFORE the early-return switch (function-evaluated types, all ranks,
  rank-deterministic) and in `fix()` (circuit types — rank 0 at
  circuit-solve time). Rank 0 is authoritative for circuit values and is
  what Exodus and the memdump write; worker copies are NOT kept in sync
  (the coupled path's `distribute_globals` runs one step behind the fix and
  the segregated path never calls it) — stated at the helper and in the
  input reference, and no current consumer reads them on workers.
- Naming, injective by construction (reshaped twice by the audits): base =
  section header label (`current : coil1 { }` → `coil1`) else type; a
  TWO-PASS scheme pre-counts base names so a repeated base carries its
  section ordinal on EVERY occurrence (never a bare `current` next to a
  `current_2`); multi-group sections append group ordinals. Thermal
  conditions prefix `thermal_` when unlabeled (both kernels share one
  mesh); circuit terminal pairs use their component `label` else
  `terminalpair_<topology index>`. Bearing publishes nothing (it imposes no
  scalar); Gauge publishes. Any surviving duplicate is a hard setup error.
- Creation happens in `MaxwellFactory::create_magnetic_kernel` and
  `ThermalFactory::create_thermal_kernel` immediately before
  `create_globals()`, on `dofmgr()->mesh()` — the real submesh on workers
  (the factory-ctor stub-mesh trap found in review), on every rank in
  identical order so auto-IDs agree. Circuit-pushed BCs land in the same
  container beforehand and are covered by the same walk. The executables
  load memdumps AFTER the factories, so the exists-guard means "duplicate
  label", never "restart leftover"; a same-deck warm restart reuses the
  same IDs and passes `load_globals`' const-ID check. A deck whose BC ORDER
  changed across a dump refuses the restart — documented limitation.
- Hardening ridden along: the thermal BC factory's silent default
  fall-through (pre-existing) became corrupting with the new labeling loop
  — an unsupported type would have relabeled the previous section's
  condition; it now hard-errors like the Maxwell twin has since 2026-08-11.

## Audit trail

Full two-phase protocol, both voices at each phase. Phase 1 rewrote the
design (six Codex findings, ten Grok sections: the worker stub mesh, the
backwards restart story, the broken auto-ID map, the swapped replication
args, the `label`-key collision, the circuit BCs outside the file list, the
early-return trap, the non-injective naming). Phase 3 passed the core and
failed three more (thermal default fall-through, the false worker-healing
claim, the two-block naming collision), all fixed; a verify pass closes the
round. Legacy's display quirk (zeroing the shown global below 1 mA while
imposing a clamped value) was deliberately not copied — the global always
shows the value actually imposed.

Status: **reviewed, not verified** — 22 syntax gates green across 8 TUs;
nothing built or run. Executable gates for Christian's next rebuild: an
np≥2 run showing the globals in ParaView with the sigmoid ramp, and a
same-deck warm restart round-trip. Input reference §9 updated (no schema
change — no key added; the header-label grammar already existed).

## Same-day follow-up: globals become rank-0-only

Christian, watching the running job stall and distrusting the added
broadcasts on principle ("we just want them written to the mesh, but that
is done by the main proc only"), cut the design down: the globals'
consumers — Exodus writer, memdump, ParaView — all read rank 0's mesh, no
worker code reads a mesh global, so the worker replication built above was
machinery without a customer. (The stall itself was refuted as a broadcast
problem — gdb backtraces on all four ranks showed STRUMPACK iterative
refinement, a separate story — but the simplification stands on its own.)

What changed, all deletion:

- The factory creation walks (Maxwell ~:852, Thermal ~:250) now run under
  a rank-0 guard.
- `DofManager::create_globals()` and `distribute_globals()` deleted
  entirely — declarations, bodies, and every call site (the per-timestep
  pair in `Controller::initialize_timestep` and the two restore-path calls
  in `load_memdump`). The `Cell<string>` label broadcast — the code that
  crashed with "Unknown datatype" the first evening it ever executed — is
  gone with them; the `broadcast( Cell<string> )` overload itself stays
  (live users: material labels).
- `update_global()` keeps its exists-guard, which on workers is now
  structurally false: the same line that was a sync-lag caveat is now the
  design. Header comment and input reference §9 rewritten accordingly.
- The restart ID check was never in the deleted code — it lives in
  `Mesh::load_globals` (rank 0) and survives; the reorder-refuses-restart
  caveat is unchanged.

Net effect: the naming rules, the ParaView visibility, and the warm-restart
persistence are exactly as designed above; what disappeared is every MPI
transfer and the entire class of collectivity hazards that came with it.
Consumer sweep confirmed empty (no worker-side reader in the tree).
Follow-up verify round: `tmp/ai_exchange/bc_globals_rank0.md` — Codex
passed all four questions (collectivity symmetry, warm-restart path,
np==1, comment/doc accuracy) at high confidence; one stale comment in
`cl_Mesh.cpp` fixed as by-catch. Reviewed, not verified: syntax gates
green on all four TUs; the executable gates are unchanged (np≥2 ParaView
check plus warm-restart round-trip on the next rebuild).
