# Terminal globals unified under `save_IV`; temperature global named by what it is

**Date:** 2026-09-01
**Purpose:** Remove the duplicate `current`/`voltage` mesh globals and name the temperature global `T_max` or `T_fixed`
**Module:** fem/kernel, fem/maxwell

## What was wrong

Christian found a `current` global in the `tape_quench` Exodus output next to
`I_1` and `U_1`. Both carried the same number. Two channels were publishing the
same terminal:

- `cl_MaxwellBoundaryConditionFactory` labelled every value-imposing block with
  its header label, else its type, and `MaxwellFactory::create_magnetic_kernel`
  turned every label into a mesh global updated by `impose_bc`. An unlabeled
  `current { }` block therefore published a global named `current` — the
  literal deck token, which is why no grep for the string found it.
- `Controller::save_IV` publishes `I_k`/`U_k` per abstract dof. Abstract dofs
  are created one per `Current`/`CircuitCurrent` condition, then one per
  `Voltage`/`CircuitVoltage` condition (`impose_voltage_bcs` asserts the
  pairing at `cl_FEM_Controller.cpp:523` and `:537`), so for exactly those four
  types the two channels are the same numbers by construction. `save_IV`
  dominates: it carries the response half too (`U` for a current-driven
  condition), is per group where the block global is per block, and is called
  unconditionally before every Exodus write in `belfem.cpp`.

Separately, the `temperature` global (renamed to `T_max` by Christian this
session) held `gTbulk` frozen at construction, while the console's `T_max` is
`max( T )` live. In a coupled quench run the Exodus `T_max` would have read the
operating temperature for the whole run.

## What changed

**One channel for terminals, deck labels folded into its names.**

- `has_block_global( BoundaryConditionType )` added to
  `en_FEM_BoundaryConditionType`: false for `Bearing` and the four terminal
  types, true otherwise. Both guards in the BC factory (pre-count and naming,
  previously `!= Bearing` twice) and the publication loop in `MaxwellFactory`
  use it, so the three sites cannot drift apart.
- Terminal conditions get the section label on **every** member (empty if the
  section has none), instead of the type fallback on the first member.
- `Controller::create_iv_names()` (new, rank 0, once) names each abstract dof's
  pair: `I_<label>`/`U_<label>` when the driving condition is labelled, with a
  running suffix where a label repeats (`I_coil1_1`, `I_coil1_2` — one block
  with several groups, or two blocks with one label), positional `I_1`/`U_1`
  (zero-padded from ten generators on) otherwise. It walks the conditions in
  the same current-then-voltage order the abstract dofs follow and asserts
  fixed↔current on each pair. Circuit terminal pairs keep the circuit factory's
  label (`component label`, else `terminalpair_<n>`), so they appear as
  `I_terminalpair_2`. Generators beyond the last condition stay positional.
- `save_IV` uses the cached names for the csv header and the globals.

**Temperature global owned by the Controller.** `Controller::save()` creates
(first save; or adopts the one a memdump restored) and refreshes a global named
`T_max` when a thermal kernel is attached and `T_fixed` otherwise, with
`get_Tmax()` — already `max( T )` vs `gTbulk` by the same test. The
`MaxwellFactory` creation block is gone; its reserved-name guard moved into the
publication loop and now covers both names. `T_ad` was proposed and withdrawn:
adiabatic describes the thermal model; a magnetic-only run holds the
temperature, it does not insulate it.

Files: `en_FEM_BoundaryConditionType.{hpp,cpp}`,
`cl_MaxwellBoundaryConditionFactory.cpp`, `cl_MaxwellFactory.cpp`,
`cl_FEM_Controller.{hpp,cpp}`; `doc/input_file_reference.md` §9 (mesh-globals
paragraphs rewritten, the stale `temperature` sentence with them) and the
`initial conditions` row; `doc/input_schema.yaml` (`boundary conditions →
mesh_globals`, `published_as` on the temperature key).

## Evidence

- All four edited translation units pass `g++ -fsyntax-only` with the debug
  tree's own `flags.make` (kernel and maxwell object libraries). Not built, not
  run: **reviewed**, not verified.
- No test references any of these globals (`grep -rn "global_variable\|I_1"
  tests/` — no hits), so `make check` cannot discriminate this change either
  way. The gate is a `tape_quench` run: the Exodus global list should read
  `I_1 U_1 dotQ T_max` and nothing named `current`; a magnetic-only deck should
  show `T_fixed`.
- Breaking for downstream readers of `current`/`voltage`/`temperature`
  globals, accepted: the block globals were one day old for terminals, and a
  pre-change `memdump.hdf5` resurrects the old names beside the new ones
  (documented; delete it once).

## Open

- Codex language sweep over the rewritten §9 paragraphs is owed.
- A jury round on the code is owed per the standing rule for executable
  changes; not run this session.
- `hphirun.cpp`/`hphiTrun.cpp` still call `save_IV`/`save` (unchanged API);
  they are not built.
