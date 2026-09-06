# rlc_circuit {#doc_examples_rlc_circuit}

**Date:** 2026-09-05
**Purpose:** Copper inductor coupled to a lumped RLC circuit.

---

## What it shows

A three-dimensional copper inductor is coupled to a lumped circuit through a
`terminal pair`. The circuit contains a 1 V, 50 Hz source, a resistor, the inductor, and a
capacitor. It is defined inline in the `circuit` section, and results are written to
`CircuitResults.txt`. The mesh uses meters. `tapestack_circuit` reads its circuit from a
SPICE netlist instead.

## Files

| file | role |
|---|---|
| `inductor.geo` | geometry (m); mesh it with gmsh |
| `input.conf` | deck, including the circuit |

## Run

```bash
gmsh -3 inductor.geo -o inductor.msh
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
