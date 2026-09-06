# disk_ramp {#doc_examples_disk_ramp}

**Date:** 2026-09-05
**Purpose:** Bulk YBCO disk magnetized by a ramped background field.

---

## What it shows

A bulk YBCO disk in a spherical air domain. A uniform background field along +z ramps to
1 T over one second and is then held. Screening currents build up in the disk; the deck is
the plain three-dimensional bulk HTS case, with the built-in power law and no thin shells,
cuts, or periodicity. `disk_pulse` uses the same geometry with a pulsed field from a
plugin.

## Files

| file | role |
|---|---|
| `disk.geo` | geometry; mesh it with gmsh |
| `input.conf` | deck |

## Run

```bash
gmsh -3 disk.geo -o disk.msh
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
