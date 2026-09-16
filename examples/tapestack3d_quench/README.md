# tapestack3d_quench {#doc_examples_tapestack3d_quench}

**Date:** 2026-09-15
**Purpose:** Quench of an eight-tape REBCO stack: the `tapestack3d` geometry driven by a
user-defined current ramp and a local Ic defect.

---

## What it shows

The geometry, layer stack and periodicity are those of [tapestack3d](@ref
doc_examples_tapestack3d) — eight 4 mm wide thin-shell tapes with solder between them, over a
10 mm length whose end planes are identified. What differs is the drive, and with it the
question the deck asks: instead of a prescribed sigmoid current, a plugin supplies both the
transport current *and* a local defect, so the run follows the current sharing from end to end.

- A linear 200 A/s ramp, as in an Ic test, starts from zero and passes 2 kA at the 10 s end of
  the run.
- At 3 s the topmost tape loses 90 % of its jc over a 1.5 mm disk at mid-length (z = 5 mm, away
  from the periodic end planes). The defect switches on logistically and is effectively fully
  on by 3.5 s.
- The solder diverts that tape's share into the seven healthy tapes, which still have margin at
  that current. The defect alone does not quench the stack.
- The ramp does not stop, so the stack runs out of margin a few seconds later and quenches,
  starting at the defect where the diverted current crosses the resistive solder.

Both plugin functions read their timing from one header, `src/ramps.hpp`, so the current and
the defect cannot drift apart — the comment block at its top is the written form of the story
above. They compile into a single plugin, `userdefect.so`; the deck names that same file in its
`materials { ybco { defect { } } }` block and in its current boundary condition.

The coupled h-phi/T problem starts at 77 K and `belfem` announces the coupled solve at startup.
Compared with `tapestack3d`, the deck also tightens both nonlinear tolerances to `1e-9`
(current sharing between eight tapes is the quantity of interest), switches the conditioning
estimates off, and turns `edge coating` off — the sharing of interest here goes through the
solder, not around the tape edges.

Although `tapestack3d.geo` defines `numTapes`, the deck's sidesets, blocks, curves, and
terminals are written explicitly for eight tapes, and `src/defect.cpp` hard-codes the tape
plane positions it gates on. After changing the geometry, regenerate the topology with
`gmsh -0 tapestack3d.geo`, update the deck, and re-check the constants near the top of
`defect.cpp`.

## Files

| file | role |
|---|---|
| `tapestack3d.geo` | geometry; mesh it with gmsh |
| `src/ramps.hpp` | the time program: current ramp and defect switch-on, shared by both plugins |
| `src/current.cpp` | user current source (`MyCurrent`) |
| `src/defect.cpp` | user Ic defect (`MyDefect`), including the defect geometry |
| `src/CMakeLists.txt` | builds both into `src/build/userdefect.so` |
| `input.conf` | deck |

## Run

Build the plugin first; the deck loads `src/build/userdefect.so`:

```bash
cd src && mkdir -p build && cd build
cmake -DBELFEM_DIR=/path/to/belfem .. && make
cd ../..
```

The plugin needs CMake 4.0 or newer, the same floor BELFEM itself asks for. A distribution's
CMake is often older, so use the one from the toolchain you build BELFEM with rather than
whichever comes first on your `PATH`.

```bash
gmsh -3 tapestack3d.geo -o tapestack3d.msh
../../build/bin/belfem                 # or mpirun -np 4 ../../build/bin/belfem
```

The first run builds `Copper_RRR50.hdf5` and `Silver_RRR10.hdf5` in the run directory; later
runs reuse them. `sp-ap.hdf5` ships in `share/material`.

## Changing the scenario

The timing lives in `src/ramps.hpp`: the ramp rate and its safety cap, and the defect's
switch-on time, transition width and fuzzyness. The defect's *shape* — where it sits, how wide
it is, how deep it goes — is the constant block near the top of `src/defect.cpp`. Rebuild the
plugin after editing either; the deck does not need to change.

Note that the ramp rate and the simulation time are tied: `simulation time : 10 s` at 200 A/s
is what puts the stack past its self-field Ic within the run. Slow the ramp without
lengthening the run and the stack never quenches inside the simulated time.

## Shared launcher

`../scripts/Allrun` generates the mesh, chooses ranks and threads, and can submit to Slurm;
`../scripts/Allclean` removes run output. To use them from this directory, copy the two
wrapper scripts `Allrun` and `Allclean` from `examples/helix/` here and run `./Allrun`
(`./Allrun --dry-run` prints the plan first). If you copy this deck out of the tree, copy
`examples/scripts/` along with it and point the wrappers at its new location. The launcher
does not build plugins; do that step by hand as described above. Per-deck settings go in a
`run.conf` next to `input.conf`; see `../scripts/README.md`.
