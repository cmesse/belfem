# Executables Module {#executables_index}

**Date:** 2026-08-12
**Purpose:** Overview of the solver applications and their command-line interface
**Module:** executables

## Overview

This module builds the solver application. `belfem` reads its problem
definition from `input.conf` in the working directory (see
`doc/input_file_reference.md`), and its command line controls the runtime
environment, not the physics. The other listed executables do not read a
deck. `material` and `gas` are property evaluation tools driven by their own
arguments; `db2exo` converts a property database to Exodus; `msh2exo` is a
mesh converter that takes its input file on the command line.

| Executable | Source | Purpose |
|---|---|---|
| `belfem` | `belfem.cpp` | h-ɸ simulation, magnetic-only or coupled h-ɸ/T — the deck decides |
| `material` | `../physics/materials/` | material-property evaluation and table tooling |
| `gas` | `../physics/gasmodels/` | gas-model evaluation |
| `db2exo` | `../physics/database/` | converts a property database to Exodus (needs HDF5 and Exodus) |
| `msh2exo` | `../mesh/` | converts a Gmsh file to an Exodus file of the same name (needs Exodus) |

`belfem` is built here; the rest come from the physics and mesh modules and
are gated: `material` always builds, `gas` needs `USE_GASMODELS=ON`
(default OFF), `db2exo` needs HDF5 and Exodus, and `msh2exo` needs Exodus.
`make install` deploys whichever executables the configuration actually built —
`BELFEM_INSTALL_EXECUTABLES`
(`config/globals.cmake`) is a filter over existing targets, not a build list.

`belfem` unifies `hphirun` and `hphiTrun`. It solves the coupled h-ɸ/T
problem when the deck's `solver` block contains an unlabeled
`linear thermal` or `nonlinear thermal` section. Otherwise, it solves the
magnetic problem only. The selected mode is printed at startup.

A deck is rejected as inconsistent if it declares thermal boundary conditions
(`boundary conditions { thermal { … } }`) without a thermal solver section.
If `initial conditions` gives no temperature, either mode assumes 77 K and
prints a console note.

`hphirun` and `hphiTrun` were the two dedicated executables `belfem` replaces.
They are **retired and no longer built**: their `Add_Executable` blocks are
commented out in `src/executables/CMakeLists.txt`, so a build produces no
`bin/hphirun` and no `bin/hphiTrun`. Their sources remain in the tree and are
still cited as reference drivers, but there is nothing to run — use `belfem`,
which chooses the mode a deck needs rather than requiring you to.
The project library's CMake target is
`belfem_lib` (the archive on disk is still `libbelfem.a`), so `make belfem`
builds this executable.

All are MPI programs:

```bash
mpirun -np <nprocs> --bind-to socket belfem
```

## Command-Line Interface

The command line is shared by three consumers, each of which scans the full
argument list and skips flags it does not recognize. BELFEM, PETSc, and
STRUMPACK options can therefore be mixed freely.

### BELFEM flags

Parsed by the `Arguments` base class (`src/core/cl_Arguments.{hpp,cpp}`),
GNU style; see `src/core/doc/core_usage_guide.md` for details.

| Form | Effect |
|---|---|
| `-v`, `--verbose` | logger info level to `Everything` (5) |
| `-v N`, `-vN`, `--verbose N`, `--verbose=N` | logger info level to N |

Info levels follow `InfoLevel` in `cl_Logger.hpp`:
0 = silent, 1 = minimal, 2 = default, 3 = detailed, 4 = verbose,
5 = everything (including third-party library output).

`belfem` additionally accepts a version flag, handled in `belfem.cpp` rather
than by the `Arguments` base class:

| Form | Effect |
|---|---|
| `-V`, `--version` | print the startup banner (version, build date, git
provenance, enabled parallel features) and exit |

The check runs before the deck is opened, so `belfem --version` works in a
directory that has no `input.conf`. As with every other startup message the
banner is printed by rank 0 only, and all ranks reach `finalize()` together.

### PETSc options

The full command line is handed to `PetscInitialize`
(`src/comm/cl_Communicator.cpp`), so every documented PETSc option works:
`-ksp_type gmres`, `-pc_type hypre`, `-ksp_monitor`,
`-options_file petsc.opts`, `-snes_*`, and so on. Warnings about unused
options are suppressed (`-options_left 0`).

### STRUMPACK options

Flags with the `--sp_` prefix are read by STRUMPACK **after** the
`input.conf`-derived settings are applied, so the command line overrides the
input file (`src/sparse/cl_SolverSTRUMPACK.cpp`). Examples:
`--sp_compression NONE`, `--sp_reordering_method metis`. See
`src/sparse/doc/sparse_usage_guide.md` ("Command-Line Pass-Through") and the
STRUMPACK manual for the full list. STRUMPACK's own progress output is
enabled at info level ≥ 5, i.e. by passing `--verbose`.

### Example

```bash
mpirun -np 4 belfem -v 3 --sp_compression NONE -ksp_monitor
```

## Common Pitfalls

- **`input.conf` is not optional** — the factories read it from the working
  directory; there is no flag to point at a different file.
- **Verbosity is per rank** — every rank parses the same command line, and
  most log messages are emitted by rank 0 only.
- **A typo in a `--sp_` or PETSc flag fails silently** — unknown flags are
  skipped by every parser, so a misspelled option is ignored rather than
  reported. Check the effective solver settings in the log at `-v 3` or
  higher.

## Related Documentation

- `src/core/doc/core_usage_guide.md` — `Arguments` and `Logger` classes
- `src/sparse/doc/sparse_usage_guide.md` — solver configuration and the
  PETSc/STRUMPACK pass-through
- `doc/input_file_reference.md` — the `input.conf` contract
