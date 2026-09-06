# tapestack2d_layered {#doc_examples_tapestack2d_layered}

**Date:** 2026-09-05
**Purpose:** 2D tape stack resolved layer by layer.

---

## What it shows

An 11-tape two-dimensional stack with each tape layer resolved as a thin shell. Each tape
includes copper, silver, YBCO, buffer, and Hastelloy, using built-in materials. The
superconductor uses the built-in power law. The stack carries a 160 A ramp. Compare with
`tapestack2d_racetrack`, which uses a table-driven superconductor and a copper bulk.

## Files

| file | role |
|---|---|
| `tapestack.geo` | geometry (mm); mesh it with gmsh |
| `input.conf` | deck |

## Run

```bash
gmsh -2 tapestack.geo -o tapestack.msh
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
