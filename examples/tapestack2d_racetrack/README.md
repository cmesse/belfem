# tapestack2d_racetrack {#doc_examples_tapestack2d_racetrack}

**Date:** 2026-09-05
**Purpose:** 2D racetrack coil cross section with a table-driven REBCO tape.

---

## What it shows

A racetrack-coil cross section (see `racetrack2d.pdf`) containing a thin-shell tape stack
and a copper bulk conductor. The superconductor uses the `sp-ap.hdf5` table, and the stack
carries a 3040 A, 50 Hz sine. The mesh uses meters. Compare with `tapestack2d_layered`,
which resolves every layer of the tape with built-in materials.

## Files

| file | role |
|---|---|
| `2D_tapestack.geo` | geometry (m); mesh it with gmsh |
| `racetrack2d.pdf` | sketch of the coil the section belongs to |
| `input.conf` | deck |

## Run

```bash
gmsh -2 2D_tapestack.geo -o 2D_tapestack.msh
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
