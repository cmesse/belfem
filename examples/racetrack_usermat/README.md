# racetrack_usermat {#doc_examples_racetrack_usermat}

**Date:** 2026-09-05
**Purpose:** Three racetrack coils with a user-defined bulk material.

---

## What it shows

Three racetrack coils, each containing a 19-tape thin-shell stack, are driven by
independent ramps of about 28 to 30 kA. `src/matlib.cpp` defines the surrounding bulk
material and builds it into `usermat.so`, which the `materials` section loads with a
`usermat` entry. The tapes use the built-in YBCO power law. Use this deck as a template for
a material that BELFEM does not provide.

`data/` is on the plugin's include path as `BELFEM_USER_DATA_DIR`. It ships no tables; its
README explains how to read measured properties from there.

## Files

| file | role |
|---|---|
| `racetrackAssembly_*.geo` | geometry; mesh it with gmsh |
| `src/matlib.cpp`, `src/CMakeLists.txt` | user material plugin |
| `data/` | tables for the plugin (empty) |
| `input.conf` | deck |

## Run

Build the plugin first; the deck loads `src/build/usermat.so`:

```bash
cd src && mkdir -p build && cd build
cmake -DBELFEM_DIR=/path/to/belfem .. && make
cd ../..
```

```bash
gmsh -3 racetrackAssembly_3coils_19tapes_r2-63.6_r2-58.22_r2-37.5_qcyl.geo \
     -o racetrackAssembly_3coils_19tapes_r2-63.6_r2-58.22_r2-37.5_qcyl.msh
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
