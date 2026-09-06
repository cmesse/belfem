# Shared example scripts {#doc_examples_scripts}

**Date:** 2026-08-31
**Purpose:** one maintained launcher and cleaner shared by the decks in
`examples/`

## Usage

From a deck directory:

```bash
./Allrun                # run locally
./Allrun --dry-run      # print the plan, touch nothing
./Allrun --slurm        # submit to Slurm, then run inside the allocation
                        # ( the partition, account and paths default to LBNL
                        #   Lawrencium; edit them in Allrun for another site )
./Allclean              # remove run output ( keeps the mesh )
./Allclean --mesh       # also remove generated meshes
```

The decks that use it (`costheta`, `corc_garber_effect`, `helix`) hold a two-line
wrapper that `cd`s to itself and execs the shared script, so the working
directory is always the deck ( the FEM executables read `input.conf` from
the cwd ). Any other deck can be run the same way with `../scripts/Allrun`.

## Per-deck overrides: `run.conf`

`run.conf` is optional. Defaults live in `Allrun`; it overrides only the
keys it sets, so a deck records its differences and inherits later
improvements. It is parsed as `KEY=value`, not sourced — a shipped
example should not execute arbitrary shell before it has printed what
it intends to do. Unknown keys are reported and ignored.

| key | default | notes |
|---|---|---|
| `EXECUTABLE` | auto | empty = decide from `input.conf`, see below |
| `NTHREADS` | 2 | OpenMP threads per rank |
| `NRANKS` | auto | only with `RANK_POLICY=fixed` |
| `RANK_POLICY` | `even` | `even` \| `power2` \| `fixed` |
| `EXTRA_ARGS` | *(none)* | e.g. `-v` for verbose |
| `MESH_BEFORE_SUBMIT` | true | mesh on the login node, not in the job |
| `MPIRUN`, `GMSH`, `LOGFILE` | | tool names / log path |
| `BELFEM_SLURM_*` | lr4 defaults | partition, account, qos, walltime, nodes, cores-per-node, bin dir |

`BELFEM_SLURM_ACTIVATE` is deliberately **not** deck-settable: that path is
`source`d inside the allocation, so honouring it from a deck would reopen
the arbitrary-shell hole that parsing `run.conf` exists to close. Set it in
the shared script.

Example, for a deck that wants verbose output on four threads per rank:

```
NTHREADS=4
EXTRA_ARGS=-v
```

No shipped deck currently needs a `run.conf`. The three former files set only
`NDIMS`, which no longer exists — see below.

## How the mesh is generated

`Allrun` meshes only when there is something to mesh. If the `.msh` named in
the deck's `mesh { file : ... ; }` block is missing, it looks for the matching
`.geo` beside it and runs gmsh on that.

Not every deck has to have one. `corc_solder` (not shipped in `examples/`) builds
its cable in `python/main.py`, and no `.geo` exists or ever will. So when the
`.geo` is absent the launcher looks
for a `main.py` in the deck or one level below it, and if it finds one it says
which script to run rather than reporting a missing geometry file — a file the
reader would otherwise go hunting for. With neither a `.geo` nor a generator,
it says that too. Either way it stops, because there is no mesh to run on.

## Why there is no dimension setting

`Allrun` always calls `gmsh -3`, and that is not a deck-settable choice.

gmsh meshes each dimension *up to* the number given. For a geometry with no
volume, the 3D pass finds nothing and returns — measured at 31 microseconds on
`2D_tapestack.geo`, with output byte-identical to `gmsh -2` on the two 2D decks
that were meshed both ways (`2D_tapestack.geo` and `costheta.geo`;
`undulator.geo` was not re-meshed). So telling it `-2` bought exactly nothing.

The reverse mistake was not free: `gmsh -2` on a 3-D geometry exits 0 and
writes a `.msh` with every tetrahedron missing — `helix.geo` yields its 26084
triangles and none of its 217223 tets — which passes the "did gmsh produce the
file?" check and reaches the solver as a surface-only mesh.

A knob whose intended setting is a no-op and whose wrong setting fails
silently is worth deleting. So it was (2026-08-31).

## How the executable is chosen

The launcher runs `belfem`, which selects the physics from the deck.
An unlabeled `linear thermal` or `nonlinear thermal` **section header**
in `input.conf` requests the coupled h-ɸ/T problem. Without either
section, the run is magnetic-only. The chosen mode is printed at startup.

`EXECUTABLE=` in `run.conf` still overrides the launcher's choice, but
`hphirun` and `hphiTrun` are no longer built — their `Add_Executable` blocks are
commented out in `src/executables/CMakeLists.txt` — so pointing it at either one
finds no binary.

A `circuit{}` section is not an executable hint either. Decks with one
are ordinary `belfem` runs (magnetic-only or coupled, as the deck decides)
whose section is consumed by `ElectricalCircuitFactory` inside `belfem`.

## Ranks and threads

`NRANKS` defaults to the largest **even** count that fits:
`2 * floor( (physical_cores / NTHREADS) / 2 )`. Assembly is serial per
rank and dominates a step, so ranks are what buy throughput.

Why even, precisely: STRUMPACK distributes the elimination tree by mapping
subtrees onto ranks in proportion to their estimated **work**, not by
bisecting the communicator. An odd rank count forces the top-level split
to be uneven — measured at 5 ranks, the sibling groups ran at 127.7
against 185.4 GFLOP/s with the faster side waiting. An even count does not
*guarantee* balance ( uneven work can still split 3-against-5 on 8 ranks );
it removes the case where imbalance is unavoidable.

`RANK_POLICY=power2` is available for factorization-dominated runs but is
not the default — on a 24-core node it would idle 8 cores to help the
phase that does not scale.

Oversubscription is refused outright. See `doc/parallel_execution.md`.
