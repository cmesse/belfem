# dipole {#doc_examples_dipole}

**Date:** 2026-09-05
**Purpose:** 2D dipole coil: the simplest deck.

---

## What it shows

A two-dimensional dipole made of two bulk coil blocks in air, driven by a 1 A sine
current of 2 s period with opposite polarity in the two blocks. Nothing else is switched on:
no superconductor, no thin shells, no iron, no coupling. Start here to see a complete deck
and the output it produces before moving to the HTS cases.

## Files

| file | role |
|---|---|
| `dipole.geo` | geometry; mesh it with gmsh |
| `input.conf` | deck |

## Run

```bash
gmsh -2 dipole.geo -o dipole.msh
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
