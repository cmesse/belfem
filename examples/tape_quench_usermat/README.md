# tape_quench_usermat {#doc_examples_tape_quench_usermat}

**Date:** 2026-09-05
**Purpose:** Quench of a single tape: coupled magneto-thermal run with user-defined materials, defect and current.

---

## What it shows

One 12.5 mm REBCO tape carries a transport current just above its critical current. A local
defect drives part of the width normal. The coupled h-phi/T run starts at 77 K and follows
the hot spot for 0.3 s. Everything a quench study needs from outside BELFEM comes through
plugins:

- `src/matlib.cpp` defines five tape materials from the measured `cp`, `k`, and `rho`
  tables in `data/`;
- `src/defect.cpp` defines the defect, sized in its comments to give an effective critical
  current of 80 A against 137 A transport;
- `src/current.cpp` reads the waveform from `data/I_vs_t_regular_smooth.txt`.

The material plugin is `usermat.so`; the defect and current source are compiled into
`userdefect.so`. `plot_iv.py` plots the current and voltage traces the run writes to
`iv_results.csv`. The annotated templates in `src/physics/materials/example_user_*.cpp` and
`src/numerics/sources/example_user_source.cpp` point to this deck as their worked example.

## Files

| file | role |
|---|---|
| `tape.geo` | geometry; mesh it with gmsh |
| `src/matlib.cpp` | user material plugin (`usermat.so`) |
| `src/defect.cpp`, `src/current.cpp` | user defect and current source (`userdefect.so`) |
| `src/user_table.hpp` | two-column table reader shared by the plugins |
| `data/*.txt` | measured properties and the current waveform |
| `plot_iv.py` | post-processing |
| `input.conf` | deck |

## Run

Build the plugin first; the deck loads `src/build/usermat.so` and `src/build/userdefect.so`:

```bash
cd src && mkdir -p build && cd build
cmake -DBELFEM_DIR=/path/to/belfem .. && make
cd ../..
```

```bash
gmsh -3 tape.geo -o tape.msh
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
