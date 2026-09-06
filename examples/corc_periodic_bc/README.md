# corc_periodic_bc {#doc_examples_corc_periodic_bc}

**Date:** 2026-09-05
**Purpose:** One-layer CORC cable with a translated periodic boundary condition.

---

## What it shows

This deck models one layer of three REBCO tapes wound around a round core, as thin shells.
A **translated** periodic condition shifts the top plane along the axis, so each tape
connects to itself. A 160 A, 50 Hz sine drives the cable for 20 ms with Picard and
STRUMPACK. Compare it with `corc_garber_effect`, which uses a twisted periodic condition.

CORC(R) is a registered trademark of Advanced Conductor Technologies LLC, Boulder, CO,
used with kind permission to identify the cable technology represented by this example
(https://www.advancedconductor.com/).

## Files

| file | role |
|---|---|
| `corc.geo` | geometry; mesh it with gmsh |
| `twisted_conductor.pdf` | sketch of the cable |
| `input.conf` | deck |

## Run

```bash
gmsh -3 corc.geo -o corc.msh
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
