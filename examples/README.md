# BELFEM Examples {#doc_examples}

**Date:** 2026-08-31
**Purpose:** How to run the shipped example decks, and what each one needs before it will run.

---

## 1. Generate the mesh first

Examples ship the **geometry**, not the mesh — a `.geo`, from which you build the `.msh`
yourself:

```bash
cd examples/helix
gmsh -3 helix.geo -o helix.msh
```

Sixteen of the eighteen decks ship a `.geo` file. `corc_twolayer` and `pancake` build their
meshes with `python/main.py` instead, run from the deck directory. `examples/scripts/Allrun`
handles both cases. Each deck has a `README.md` describing what it shows and how to run it.

Use the `.geo` basename that the deck's `mesh { file : ... ; }` entry names. gmsh 4.x writes
format 4.1, which is what the periodic and cut machinery expects — several decks address
geometry entities (`topology { periodic { source : 2,3,4 ; ... } }`) by the **vertex IDs** that
only a format-4.x file carries.

## 2. Run

```bash
../../build/bin/belfem
```

The `bin/` path follows your build directory — `build/bin/` after the quick start in the
[Getting Started](@ref doc_getting_started) guide, `cmake-build-debug/bin/` in a stock IDE setup.

`belfem` reads the deck and selects the physics automatically. An unlabeled `linear thermal` or
`nonlinear thermal` solver section requests the coupled h-ɸ/T problem. Otherwise the run is
magnetic-only. It announces the choice at startup, so there is nothing to pass on the command
line and nothing to get wrong.

`hphirun` and `hphiTrun` are **retired and no longer built** — their `Add_Executable` blocks are
commented out in `src/executables/CMakeLists.txt`, so no `bin/hphirun` is produced. They were the
magnetic-only and coupled halves of what `belfem` now selects on its own. Their `.cpp` sources are
still in the tree and are still cited as reference drivers, but there is nothing to run.

A deck with a `circuit { }` section is an ordinary `belfem` run, magnetic-only or coupled as
the deck decides. Its section is consumed by `ElectricalCircuitFactory` inside `belfem`.

`examples/scripts/` holds a shared launcher that handles the runs above, including mesh
generation, rank and thread selection, and optional Slurm submission — see
`examples/scripts/README.md`. Per-deck settings go in a `run.conf` beside the deck's
`input.conf`.

The `Allclean` script in the same directory removes generated output: results (`*.exo`,
`*.e-s.*`, `*.csv`), the cached `*.bfm` mesh and `memdump.hdf5`, run logs, circuit output and
Slurm files. It deliberately leaves your `.msh` and any material databases in place. `Allclean
--mesh` also removes a `.msh` — but only where a `.geo` sits beside it, since a deck shipping a
mesh without its geometry cannot regenerate one.

## 3. First run builds material databases

A deck using a metal with an `RRR` value and angle-dependent resistivity (copper, silver)
builds a lookup database on first use and caches it as `<Label>_RRR<value>.hdf5` — for example
`Copper_RRR50.hdf5`. This takes a while; later runs in the same directory reuse it.

Two things worth knowing:

- **The cache is directory-local.** Copying a run directory copies the databases with it;
  starting a fresh directory rebuilds them.
- **A database written by an older BELFEM is detected and rebuilt automatically**, with a
  warning naming the file. Older files stored the ratio as `rrr`; the current format stores
  `RRR`.

## 4. Running in parallel

```bash
mpirun -np 4 ../../cmake-build-debug/bin/belfem
```

Nothing extra is needed — a first parallel run builds any missing material database itself.

## 5. The decks

Each of the eighteen decks has its own `input.conf` and `README.md`. The columns show which
BELFEM features each deck exercises.

| deck | dim | conductors | thermal | periodic | circuit | plugins | mesh from |
|---|---|---|---|---|---|---|---|
| `dipole` | 2D | bulk copper | | | | | `.geo` |
| `costheta` | 2D | bulk + iron (B-H) | | | | | `.geo` |
| `gantry` | 2D | 464 thin-shell tapes + iron | | | | | `.geo` |
| `tapestack2d_layered` | 2D | thin shells, every tape layer resolved | | | | | `.geo` |
| `tapestack2d_racetrack` | 2D | thin shells (table) + copper bulk | | | | | `.geo` |
| `undulator2d` | 2D | 6 thin-shell groups + iron | | | | source | `.geo` |
| `helix` | 3D | 4 bulk copper helices | | translated | | | `.geo` |
| `disk_ramp` | 3D | bulk YBCO | | | | | `.geo` |
| `disk_pulse` | 3D | bulk YBCO | | | | source | `.geo` |
| `rlc_circuit` | 3D | bulk copper inductor | | | inline RLC | | `.geo` |
| `tapestack_circuit` | 3D | 2 thin-shell tapes | | | SPICE netlist | | `.geo` |
| `tapestack3d` | 3D | 8 thin-shell tapes + solder | coupled | yes | | | `.geo` |
| `pancake` | 3D | 8-tape pancake coil | | yes | | | `python/main.py` |
| `corc_periodic_bc` | 3D | 1-layer CORC, 3 tapes | | translated | | | `.geo` |
| `corc_garber_effect` | 3D | 1-layer CORC, 3 tapes | | twisted | | | `.geo` |
| `corc_twolayer` | 3D | 2-layer CORC, 6 tapes | | twisted | | | `python/main.py` |
| `racetrack_usermat` | 3D | 3 racetrack coils, thin shells | | | | material | `.geo` |
| `tape_quench_usermat` | 3D | 1 thin-shell tape | coupled | | | material, defect, source | `.geo` |

**Where to start.** Start with `dipole`, the smallest complete deck. Then try `helix` for
periodic bulk conductors, `disk_ramp` for a bulk superconductor, `tapestack2d_layered` for
thin shells, `corc_periodic_bc` for periodic thin shells, `rlc_circuit` for a lumped circuit,
or `tapestack3d` for the coupled thermal problem.

### Decks that load a plugin

`disk_pulse`, `racetrack_usermat`, `tape_quench_usermat` and `undulator2d` name a shared
object in their deck (`file : src/build/usermat.so ;` and siblings), so the two-step recipe
above is not enough. Build the deck's `src/` first, pointing it at this clone:

```bash
cd examples/<deck>/src
mkdir build && cd build
cmake -DBELFEM_DIR=/path/to/belfem .. && make
```

That writes `src/build/*.so`, which is where the deck looks. Then generate the mesh and run
`belfem` as above. The three plugin kinds, user material, user source and user defect,
each have an annotated template in the source tree: `src/physics/materials/example_user_material.cpp`,
`src/numerics/sources/example_user_source.cpp`, `src/physics/materials/example_user_defect.cpp`.

### Regenerating a mesh for a periodic deck

If you regenerate a mesh for a deck that declares `periodic { source : … ; target : … ; }` or
`bearing { nodes : … ; }`, note that those entries name **vertices**, and BELFEM builds
vertices from **point elements** (gmsh element type 15), taking the id from the entity or
geometry tag. gmsh writes point elements only for points in a physical group, or when saving
all elements — so a mesh produced without them will abort in `MaxwellFactory::create_periodic`
with *"Key not found in map"* (in a release build `BELFEM_ASSERT` is compiled out, so the
readable "invalid vertex id" message in `Mesh::vertex` never appears). This is a property of
how the mesh was written, not of the format version.

### Do not write a `bearing`

Gauge pinning is automatic: BELFEM sets one pin per connected φ component and persists it in
the `.bfm`, so no deck here writes a `bearing` any more. The key still works, but a
hand-written `bearing { nodes : … ; }` **replaces** the automatic pins rather than adding to
them — every φ component the deck does not name is then left floating. It is an expert
function and is not recommended; `doc/input_file_reference.md` has the details.
