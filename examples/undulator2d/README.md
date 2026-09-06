# undulator2d {#doc_examples_undulator2d}

**Date:** 2026-09-05
**Purpose:** 2D HTS undulator with a user-defined current source.

---

## What it shows

A two-dimensional undulator with six thin-shell tape groups, an iron yoke, and the
`sp-ap.hdf5` REBCO table at 15 K. `src/current.cpp` defines a trapezoidal current waveform
and builds it into `usersource.so`. Two `current` entries load the same library for
opposite winding directions. `undulator2d.pdf` shows the geometry; `ParamsGeo.dat` holds the
parameters `undulator.geo` reads.

`data/` is on the plugin's include path as `BELFEM_USER_DATA_DIR`. It ships no tables; its
README explains how to read a measured waveform from there instead.

## Files

| file | role |
|---|---|
| `undulator.geo`, `ParamsGeo.dat` | geometry (m); mesh it with gmsh |
| `src/current.cpp`, `src/CMakeLists.txt` | user source plugin |
| `data/` | tables for the plugin (empty) |
| `undulator2d.pdf` | sketch |
| `input.conf` | deck |

## Run

Build the plugin first; the deck loads `src/build/usersource.so`:

```bash
cd src && mkdir -p build && cd build
cmake -DBELFEM_DIR=/path/to/belfem .. && make
cd ../..
```

```bash
gmsh -2 undulator.geo -o undulator.msh
../../build/bin/belfem                 # or mpirun -np 4 ../../build/bin/belfem
```

## Shared launcher

`../scripts/Allrun` generates the mesh, chooses ranks and threads, and can submit to Slurm;
`../scripts/Allclean` removes run output. To use them from this directory, copy the two
wrapper scripts `Allrun` and `Allclean` from `examples/helix/` here and run `./Allrun`
(`./Allrun --dry-run` prints the plan first). If you copy this deck out of the tree, copy
`examples/scripts/` along with it and point the wrappers at its new location. The launcher
does not build plugins; do that step by hand as described above. Per-deck settings go in a
`run.conf` next to `input.conf`; see `../scripts/README.md`.
