# tape_with_contact_and_defect {#doc_examples_tape_with_contact_and_defect}

**Date:** 2026-10-08
**Purpose:** One REBCO tape, coupled magneto-thermal and periodic, with a contact resistance between the REBCO film and the silver cap and a Jc defect that forces current across it.

---

## What it shows

A single 4 mm tape from the `tapestack3d` deck, without solder or neighboring tapes, sits
in a periodic air cylinder that is 10 mm long. A 150 A sigmoid drives the coupled h-phi/T
problem from 77 K. The superconductor uses `sp-ap.hdf5` with `resistivity type : riva`:
the power law runs in parallel with the normal-state channel and remains finite when a
defect or a temperature near T_crit drives the film normal.

The contact is the `rint` layer, which is 0.1 µm thick and sits between the silver cap and
the REBCO film. It comes from a user material, as in `tape_quench_usermat`. The plugin
sets its resistivity to `rho = R_ct / h`. For a thin layer, the resistance to current
crossing it is `rho * h` per unit area, so `rint` imposes the specific contact resistance
`R_ct` on current moving between the REBCO and the stabilizer.

In a uniform periodic tape, current never needs to cross the contact. The deck therefore
adds a permanent Jc defect at mid-length. It is 2 mm long, spans the full width, and has
smoothed ends and a Jc multiplier of 0.35. This represents a local loss of REBCO
thickness, the defect definition of Wozniak et al. 2025 (Sec. III-E). At 77 K, `sp-ap`
gives the sound tape roughly 280 A before the self-field reduction, and the defect section
a critical current of roughly 100 A. Because the defect leaves no superconducting path
around it, the current above that has to leave the REBCO through `rint`, run along the
stabilizer past the defect, and cross the contact again on the other side.

`src/matlib.cpp` also carries, commented out, an alternative defect: a 0.4 mm crack
that enters from one edge, the typical shape of slitting or handling damage. It needs
a mesh resolution of 20 to 50 µm around the crack, which the default `tape.geo` does
not provide.

## Where the contact sits and how large it is

| interface | contact layer? | why |
|---|---|---|
| REBCO / Ag cap | **yes, `rint`** | the interfacial resistance that controls current transfer into the stabilizer |
| REBCO / buffer / Hastelloy | no | the oxide buffer is already an insulator |
| Ag / Cu, Ag / Hastelloy | no | metallic bonds whose resistance is negligible next to the REBCO/Ag interface |
| tape / solder / next tape | not in this deck | applies only to stacks and coils (turn-to-turn) |

Specific contact resistance of the REBCO/Ag interface, as tabulated in
`src/fem/maxwell/doc/contact_impedance_theory.md` §10:

| case | `R_ct` [Ω·m²] | [nΩ·cm²] | `rho` at h = 0.1 µm [Ω·m] |
|---|---|---|---|
| as-supplied commercial tape | 2.5e-12 | 25 | 2.5e-5 |
| **typical (deck default)** | **1e-11** | **100** | **1e-4** |
| degraded or aged interface | 1e-10 | 1000 | 1e-3 |

Values above about 1e-4 Ω·m make the magnetic system much harder to condition.
Turn-to-turn contacts in no-insulation coils are 1e-9 to 1e-7 Ω·m², which this layer model
is not designed for.

## Knobs

All three knobs live in `src/matlib.cpp`; rebuild the plugin after changing them.

| knob | where | note |
|---|---|---|
| `rint_Rct` | top of `matlib.cpp` | the specific contact resistance, Ω·m² |
| `rint_h` | top of `matlib.cpp` | must match the `rint` thickness in `input.conf` |
| defect size and depth | `my_defect()` | the refined box in `tape.geo` covers z = 3 to 7 mm; widen it if the defect grows |

Point ids in `tape.geo` double as mesh node ids, which the `periodic` block in
`input.conf` relies on. Keep them gap-free when editing the geometry.

## Files

| file | role |
|---|---|
| `tape.geo` | geometry; mesh it with gmsh |
| `src/matlib.cpp` | user plugin `usermat.so`: the `rint` material and the `Defect` |
| `src/user_table.hpp` | two-column table reader for the plugin |
| `data/cp_HTS.txt`, `data/k_HTS.txt` | thermal properties of `rint` |
| `input.conf` | deck |

## Run

Build the plugin first; the deck loads `src/build/usermat.so`. The CMakeLists requires
CMake 4.0 or newer:

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

`../scripts/Allrun` generates the mesh, chooses ranks and threads, and can submit to
Slurm. `../scripts/Allclean` removes run output. To use them from this directory, copy the
two wrapper scripts `Allrun` and `Allclean` from `examples/helix/` here and run `./Allrun`
(`./Allrun --dry-run` prints the plan first). If you copy this deck out of the tree, copy
`examples/scripts/` along with it and point the wrappers at its new location. The launcher
does not build plugins; build the plugin by hand as described above. Per-deck settings go
in a `run.conf` next to `input.conf`; see `../scripts/README.md`.
