# tapestack3d {#doc_examples_tapestack3d}

**Date:** 2026-09-05
**Purpose:** 3D stack of eight REBCO tapes, coupled magneto-thermal, periodic.

---

## What it shows

An eight-tape thin-shell stack with solder between the tapes and a periodic boundary
condition. A 2 kA sigmoid drives the coupled h-phi/T problem over 10 s; `belfem` announces
the coupled solve at startup. The superconductor uses `sp-ap.hdf5` with piecewise
resistivity. Commented-out entries show where the Coulomb-gauge penalty (`chi`) and the
Nitsche ghost penalty (`eta`) would go.

Although `tapestack3d.geo` defines `numTapes`, the deck's sidesets, blocks, curves, and
terminals are written explicitly for eight tapes. After changing the geometry, regenerate
the topology with `gmsh -0 tapestack3d.geo` and update the deck.

## Files

| file | role |
|---|---|
| `tapestack3d.geo` | geometry; mesh it with gmsh |
| `input.conf` | deck |

## Run

```bash
gmsh -3 tapestack3d.geo -o tapestack3d.msh
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
