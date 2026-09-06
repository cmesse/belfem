# corc_twolayer {#doc_examples_corc_twolayer}

**Date:** 2026-09-05
**Purpose:** Two-layer CORC cable with a Python-generated mesh.

---

## What it shows

This deck models two counter-wound layers of three tapes on a straight core, with solder
between the layers. A twisted periodic condition represents one period of an infinite
cable. A 200 A ramp over 10 s drives the model for 15 s with Newton, Anderson
stabilization, MUMPS, and GMRES. The thermal sections are commented out.

`python/main.py` generates the mesh, writes `corc.msh`, and prints the matching `topology`
section in `corc_topology.conf`. The comments at the top of `main.py` explain the pitch, the
number of turns, and which parameter choices are currently meshable.

CORC(R) is a registered trademark of Advanced Conductor Technologies LLC, Boulder, CO,
used with kind permission to identify the cable technology represented by this example
(https://www.advancedconductor.com/).

## Files

| file | role |
|---|---|
| `python/main.py` | mesh generator; run it from this directory |
| `python/corc/`, `python/mesh/` | packages the generator uses |
| `python/corc/tests/` | unit tests for the generator |
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
