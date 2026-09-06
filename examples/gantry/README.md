# gantry {#doc_examples_gantry}

**Date:** 2026-09-05
**Purpose:** 2D upper-half model of an HTS gantry dipole with 464 tapes.

---

## What it shows

This is the upper half of the proton-therapy gantry dipole described by Rudeiros Fernandez
et al. (IEEE Trans. Appl. Supercond. 32(6), 2022). It contains 464 Bi-2223 thin-shell tapes
in 8 coils, an iron yoke with a B-H curve, and a midplane symmetry condition. The current
ramps to 340 A over 10 s at 12 K. The header of `input.conf` documents how the gmsh entity
tags map to blocks and sidesets, and why the 464 hand-written cohomology cuts of the old
deck are no longer needed.

`bscco-2223.hdf5` is supplied in `share/material/`. Replace it when data for the actual
conductor is available.

## Files

| file | role |
|---|---|
| `gantry.geo` | geometry; mesh it with gmsh |
| `input.conf` | deck, with a long header describing the model |

## Run

```bash
gmsh -2 gantry.geo -o gantry.msh
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
