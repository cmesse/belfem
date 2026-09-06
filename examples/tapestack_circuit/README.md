# tapestack_circuit {#doc_examples_tapestack_circuit}

**Date:** 2026-09-05
**Purpose:** Two thin-shell tape stacks driven from a SPICE netlist.

---

## What it shows

A three-dimensional two-tape stack whose terminal pairs are driven by a circuit from
`tapestack.cir`. The netlist contains a 500 A, 10 Hz current source in parallel with an
inductor and a capacitor; its `* belfem:` comment lines set the integration order of the
reactive components. The run lasts 150 ms.

## Files

| file | role |
|---|---|
| `tapestack3d.geo` | geometry; mesh it with gmsh |
| `tapestack.cir` | SPICE netlist of the lumped part |
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
