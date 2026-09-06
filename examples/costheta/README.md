# costheta {#doc_examples_costheta}

**Date:** 2026-09-05
**Purpose:** 2D cos-theta dipole with an iron yoke.

---

## What it shows

A cos-theta dipole cross section with 68 bulk conductor blocks, air, and an iron yoke using
the `RoxieIron` B-H curve from `bhdata.hdf5`. The current ramps to 8498 A over 10 s. The deck
models one quarter of the magnet and demonstrates the two symmetry planes: the x-axis plane
keeps the current sign and is declared `symmetry`; the y-axis plane inverts it and is
declared `antisymmetry`. The comments in `input.conf` explain why.

## Files

| file | role |
|---|---|
| `costheta.geo` | geometry; mesh it with gmsh |
| `input.conf` | deck |
| `Allrun`, `Allclean` | wrappers for the shared launcher |

## Run

```bash
gmsh -2 costheta.geo -o costheta.msh
../../build/bin/belfem                 # or mpirun -np 4 ../../build/bin/belfem
```

Or `./Allrun`, which does both steps.

## Shared launcher

`../scripts/Allrun` generates the mesh, chooses ranks and threads, and can submit to Slurm;
`../scripts/Allclean` removes run output. To use them from this directory, copy the two
wrapper scripts `Allrun` and `Allclean` from `examples/helix/` here and run `./Allrun`
(`./Allrun --dry-run` prints the plan first). If you copy this deck out of the tree, copy
`examples/scripts/` along with it and point the wrappers at its new location. The launcher
does not build plugins; do that step by hand as described above. Per-deck settings go in a
`run.conf` next to `input.conf`; see `../scripts/README.md`.
