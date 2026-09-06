# pancake {#doc_examples_pancake}

**Date:** 2026-09-05
**Purpose:** Pancake coil sector with a Python-generated mesh.

---

## What it shows

An eight-tape pancake coil wound 1.75 turns along a Lame spiral, represented by a
30-degree sector of a twelve-coil machine. The sector planes are periodic; the leads exit
radially and end on the outer cylinder, where the stack end faces are the current
terminals. The tapes use `superox.hdf5`, and a 10 kA sigmoid ramps over 10 s.

`python/main.py` generates the spiral, leads, air box, mesh, and matching `topology`
section (`pancake_topology.conf`). Edit the parameters at the top of `main.py` to change
the coil.

## Files

| file | role |
|---|---|
| `python/main.py` | mesh generator; run it from this directory |
| `python/pancake/`, `python/mesh/` | packages the generator uses |
| `python/tests/` | unit tests for the generator |
| `input.conf` | deck |

## Run

This deck ships no `.geo`. The mesh is built by a Python tool, which also prints the matching
`topology` section for `input.conf`:

```bash
python3 python/main.py                 # writes the .msh into the current directory
../../build/bin/belfem                 # or mpirun -np 4 ../../build/bin/belfem
```

Run it from this directory so the mesh lands where `input.conf` expects it. The tool needs
`numpy`, `scipy`, `sympy` and `matplotlib`, and calls the `gmsh` executable.

## Shared launcher

`../scripts/Allrun` generates the mesh, chooses ranks and threads, and can submit to Slurm;
`../scripts/Allclean` removes run output. To use them from this directory, copy the two
wrapper scripts `Allrun` and `Allclean` from `examples/helix/` here and run `./Allrun`
(`./Allrun --dry-run` prints the plan first). If you copy this deck out of the tree, copy
`examples/scripts/` along with it and point the wrappers at its new location. The launcher
does not build plugins; do that step by hand as described above. Per-deck settings go in a
`run.conf` next to `input.conf`; see `../scripts/README.md`.
