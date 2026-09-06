# Warm-restart regression: the load_fields length guard rejected the BDF history shells

**Date:** 2026-09-02
**Purpose:** Session record — a restart from `memdump.hdf5` aborted in `Mesh::load_fields`
("field 'edge_h0' has 22838 entries, the mesh field has 0") one day after commit 218d9101
added the guard; root cause, the one-condition fix, and the two audit rounds.
**Module:** `src/mesh`

## Symptom

`cmake-build-debug/dipole` (h-φ, MUMPS, BDF, no thermal): run, stop, restart. `load_memdump`
dies at `cl_Mesh.cpp:3227` on the first history field of the dump. The dump is from
2026-08-30, written by a pre-218d9101 binary, and its checksum matches.

## Root cause (confirmed by both auditors, high)

218d9101 made `Mesh::load_fields` require the dumped length to equal the mesh field's
current length before `hdf5::load_vector_from_file` writes it. Three facts break that:

- Edge and face fields are created with length 0 — `cl_Mesh_Field.cpp:52-56`, "size must
  be assigned later". The primary dof fields are sized afterwards by
  `DofManager::create_fields` with the entity multiplicity (`cl_FEM_DofManager.cpp:185-234`,
  called from `cl_MaxwellFactory.cpp:1027`), so `edge_h` passes.
- The BDF history levels `edge_h0..` are created by `IWG_Timestep::create_old_dof_fields`
  (`cl_IWG_Timestep.cpp:302-357`) and nothing sizes them before the load: `create_fields`
  walks only `all_fields`, which deliberately omits the numbered levels. On a cold run they
  get a length at the first `shift_fields` (`Vector` assignment, `:378`). On a warm restart
  the loader itself sized them — `load_vector_from_file` calls `set_size( tLength )` before
  it reads (`hdf5_tools.hpp:672-673`). `Controller::synchronize_history_fields`
  (`cl_FEM_Controller.cpp:5102-5170`) is written on exactly that assumption ("rank 0 holds
  the loaded levels; workers hold the empty shells").
- `load_memdump` loads the fields (`:5368`) before `dofmgr()->initialize( true )` (`:5474`),
  so at guard time `edge_h0` is empty and the dump holds `edge_multiplicity × num_edges`.

The guard's comment — "a dump of a different length would overrun or truncate it silently" —
was false: the loader reallocates. The discretization check the commit was after is its
other half, the `num_edges` / `num_faces` meta check in `Mesh::load_meta` (`:3062-3081`).

## Fix

One condition in `Mesh::load_fields`: an EMPTY target is accepted and sized by the loader; a
SIZED target must still match.

```cpp
BELFEM_ERROR( tTarget.length() == 0
              || tLength == ( hsize_t ) tTarget.length(), ... );
```

Why not exempt the EDGE/FACE entity types, and why not revert the guard: both auditors gave
the same counter-example independently. An order-1 dump loaded by an order-2 deck (or the
reverse) on the same mesh passes the checksum and the `num_edges` check; the sized `edge_h`
guard is the only thing that catches the 2× mismatch. Grok added that an entity-type
exemption would not merely load the wrong length — `create_fields` would then find
`length() != tFieldSize` and `set_size( …, 0.0 )` it (`cl_FEM_DofManager.cpp:248-250`), a
silent wipe of the restored H field. Pre-218d9101 dumps carry no `num_edges`, so for them
the sized-field guard is also the only ghost-switch (`eta`) protection.

Downstream coverage is unchanged: `synchronize_history_fields` verifies every restored
level against its parent's length on every rank before the integrator is told it has
history.

## Audits

Plan round (Codex `gpt-5.6-terra`/high, Grok `grok-4.6`/high; exchange
`tmp/ai_exchange/memdump_history_shell_length_guard.md`): both confirm, both reject the
entity-type exemption with the order-1/order-2 pair. Two of my citations were corrected by
Grok (sizing is `create_fields`, not the DofManager constructor; `initialize` is `:5474`).
Nothing else non-empty and non-authoritative at load time was found in-tree: thermal `T0`
is a node field sized at creation; `element_rho`, `phi`, `Hx..Bz` likewise; circuit state
takes its own group. Latent, not live: a facet (`lambda`) history field would be created at
`number_of_facets` while its parent is `lambda_multiplicity × number_of_facets` — but
`mLambdaDofMultiplicity` is initialised to 0 and never assigned in `src/`. If lambda dofs are
ever enabled, size the history in `create_old_dof_fields`; do not weaken the guard.

Code round (same depth): both accept. Codex: for every non-empty target the original
equality is evaluated unchanged; the only new accept is an empty existing field with a
non-empty dump. Grok: `length()` is `size_t`, `||` short-circuits, nothing else in
`load_fields` changed. Two comment nits taken (node fields are sized by the Field ctor, not
the dof manager; the ctor is why the shell is empty, not where it is sized later). Both
repeat the same caveat, kept as a recorded contract rather than a code change: "empty means
history shell" is a production call-order fact — all three drivers run the factory's
`create_fields` before `load_memdump` — not a `load_fields` invariant; a future caller that
reached it with an unsized primary `edge_h` would have the loader size it and
`create_fields` wipe it afterwards. Grok also notes the pathological 0-edge mesh, where a
primary edge field is legitimately empty, is indistinguishable from a shell.

## Evidence ladder

Syntax-only compile of `cl_Mesh.cpp` with the debug tree's `flags.make`: green. **Not built,
not run** — the gate is Christian's rebuild plus a restart of the dipole deck from the
existing dump, which must pass `load_memdump` and take its first step. No unit test covers
`load_fields` on an empty edge shell; one needs an IWG + DofManager fixture and is flagged,
not written. Reviewed, not verified.
