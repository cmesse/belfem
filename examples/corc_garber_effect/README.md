# corc_garber_effect {#doc_examples_corc_garber_effect}

**Date:** 2026-09-05
**Purpose:** One-layer CORC cable with a twisted periodic boundary condition.

---

## What it shows

This deck uses a one-layer, three-tape CORC geometry with a **twisted** periodic condition.
The top plane is the bottom plane rotated by the twist of one period, so each tape continues
into its neighbor and the three tapes form one long conductor. A 90 A, 50 Hz sine drives the
model for 20 ms with Newton and MUMPS. Compare it with `corc_periodic_bc` to see the
difference between translated and twisted periodicity.

CORC(R) is a registered trademark of Advanced Conductor Technologies LLC, Boulder, CO,
used with kind permission to identify the cable technology represented by this example
(https://www.advancedconductor.com/).

## Files

| file | role |
|---|---|
| `corc.geo` | geometry; mesh it with gmsh |
| `input.conf` | deck |
| `Allrun`, `Allclean` | wrappers for the shared launcher |

## Run

```bash
gmsh -3 corc.geo -o corc.msh
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
