# CLI Verbosity Flag and Solver Option Pass-Through

**Date:** 2026-08-12
**Purpose:** GNU-style `-v/--verbose` across all executables; document the existing PETSc/STRUMPACK command-line pass-through
**Modules:** core, executables, physics/gastables, physics/gasmodels, sparse (docs only)

## Context

Request: `hphirun` and `hphiTrun` should accept `--verbose N` / `-v N` to set the
logger info level (bare flag → `InfoLevel::Everything`), and command-line
arguments should reach PETSc and STRUMPACK.

## Finding: the pass-through already existed (no code change)

- **PETSc** — `Communicator::init` passes the full `argc`/`argv` to
  `PetscInitialize`, which loads every argument into the options database;
  `-options_left 0` is set beforehand so PETSc does not warn about
  BELFEM-specific flags. All `-ksp_*`, `-pc_*`, `-snes_*`, `-options_file`, …
  options already worked in every executable. (high confidence — read directly
  from `cl_Communicator.cpp`)
- **STRUMPACK** — `cl_SolverSTRUMPACK` copies `gComm.arguments()` and calls
  `options().set_from_command_line()` after applying the `input.conf`-derived
  settings, on both the serial and the distributed path, so `--sp_*` flags
  already worked and correctly override the input file. The wrapper also
  already enables STRUMPACK's own verbose output at info level ≥ 5.

Both mechanisms are now documented in
`src/sparse/doc/sparse_usage_guide.md` ("Command-Line Pass-Through").

## Changes

1. **`Logger::set_info_level( uint / InfoLevel )`** — new runtime setter
   (`src/core/cl_Logger.{hpp,cpp}`); `mInfoLevel` was previously
   constructor-only, and `gLog` is a global constructed before `main()`.
2. **`belfem::Arguments` base constructor parses verbosity**
   (`src/core/cl_Arguments.{hpp,cpp}`): `-v [N]`, `-vN`, `--verbose [N]`,
   `--verbose=N`; bare flag → `InfoLevel::Everything`. Every executable that
   constructs an `Arguments` (or subclass) gets the flag with identical
   semantics — `material` and the visualizer inherit it with no change.
   Malformed `--verbose=x` is a `BELFEM_ERROR` (setup path). The dead
   commented-out `gComm.set_arguments` block kept its place; the unused
   `commtools.hpp` include was dropped.
3. **Executables wired:** `hphirun`, `hphiTrun`, `electricalCircuit` now
   construct an `Arguments` right after `gComm.init()`, before any factory
   runs, so setup logging honors the level.
4. **GNU-style consistency:** `gastables::Arguments` version flag flipped
   `-v` → `-V` (`--version` unchanged) per the GNU pairing
   `-v`=verbose / `-V`=version. Mildly breaking for anyone typing
   `gastable -v`. Help texts of `gastable` and `gas` now list both flags
   (they listed neither before).
5. **Docs:** core usage guide (`Arguments` section: built-in flags table),
   sparse usage guide (pass-through section as above).
6. **New module doc directory `src/executables/doc/`** — README covering the
   three solver applications and the full command-line interface (BELFEM
   flags, PETSc options, `--sp_*` flags, pitfalls). Registered in
   `doc/README.md` and in the CLAUDE.md module-doc list;
   `check_doc_claims.py` passes (32/32).

## Design notes

- Flag routing follows the pass-through model already in place: every parser
  sees the whole command line and skips what it does not recognize. The
  `-v` tokens land harmlessly in PETSc's options database and are ignored by
  STRUMPACK; conversely `--sp_*`/`-ksp_*` are skipped by `Arguments`.
- The verbosity level is set on **all ranks** (mpirun delivers identical argv).
- No `input.conf` key involved — CLI only, so the input contract
  (`doc/input_file_reference.md` / `doc/input_schema.yaml`) is untouched.

## Verification

Syntax-checked (g++ `-fsyntax-only` with the build tree's `flags.make` flags)
for all six touched TUs: `cl_Arguments.cpp`, `cl_Logger.cpp`, `hphirun.cpp`,
`hphiTrun.cpp`, `electricalCircuit.cpp`, `gastables/main.cpp` +
`cl_GT_Arguments.cpp`, `gasmodels/main.cpp`. **Reviewed, not verified** — no
build or run was performed (user runs builds).
