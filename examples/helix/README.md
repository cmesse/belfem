# helix {#doc_examples_helix}

**Date:** 2026-09-05
**Purpose:** Four helical bulk conductors crossing a periodic plane.

---

## What it shows

Four helical copper conductors lie in an air cylinder, one period long, with a translated
periodic boundary condition. The current ramps to 1 A over 10 s at room temperature. This
deck demonstrates periodic boundaries on bulk conductors and serves as the volume-edge
regression case. The `periodic` entry names gmsh **vertices**; the section "Regenerating a
mesh for a periodic deck" in `../README.md` explains what the mesh has to carry.

## Files

| file | role |
|---|---|
| `helix.geo` | geometry; mesh it with gmsh |
| `input.conf` | deck |
| `Allrun`, `Allclean` | wrappers for the shared launcher |

## Run

```bash
gmsh -3 helix.geo -o helix.msh
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
