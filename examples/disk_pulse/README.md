# disk_pulse {#doc_examples_disk_pulse}

**Date:** 2026-09-05
**Purpose:** Bulk YBCO disk under a background-field pulse from a user-defined source.

---

## What it shows

This is the `disk_ramp` case with a pulse instead of a ramp. The field rises to 1 T over
100 ms, holds for 100 ms, drops to zero, and the run continues for 300 ms. Screening
currents remain in the disk after the field is removed.

`src/bgpulse.cpp` defines the waveform and is compiled into `bgpulse.so`. The `background`
entry loads it with `type : userdefined`. Note the `units` key there: for a plugin source the
dimension of the returned value comes from that key, never from an amplitude.

## Files

| file | role |
|---|---|
| `disk.geo` | geometry; mesh it with gmsh |
| `src/bgpulse.cpp`, `src/CMakeLists.txt` | user source plugin |
| `input.conf` | deck |

## Run

Build the plugin first; the deck loads `src/build/bgpulse.so`:

```bash
cd src && mkdir -p build && cd build
cmake -DBELFEM_DIR=/path/to/belfem .. && make
cd ../..
```

```bash
gmsh -3 disk.geo -o disk.msh
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
