# pancake_tape {#doc_examples_pancake_tape}

**Date:** 2026-09-24
**Purpose:** Single-tape pancake coil whose solder is one solid that touches itself.

---

## What it shows

One 4 mm tape wound 200 turns from r = 7 mm to r = 17 mm around the z axis, a turn pitch of
50 µm. The solder between the turns forms one solid that touches itself, and its contact
surface is the tape, which is a thin shell. A solder sliver sits inside the first turn. The inner lead
leaves the spiral along a circle of radius 3 mm tangent to its inner end, wraps 90 degrees into
the bore, bends the hard way in the unrolled plane of that small cylinder until it runs along z,
and climbs to the top cap on a solder bar. The outer lead bends the easy way in the plane and
leaves radially to the domain cylinder on a solder bar. The bars' end faces are the current
terminals; the tape's terminal curves attach to them.

The mesh is structured where it matters: it has one hexahedral cell across the pitch, six Kuhn
tetrahedra per cell, and air meshed by gmsh around the closed shell. BELFEM opens the thin
shell block by block, so the generator cuts the solder into half-turn blocks and the tape into
half-turn sidesets on the same columns. The generator refuses any mesh in which a tape facet does not
separate two different blocks, air sits on the tape's inner side, a block lies on
both sides of one sideset piece, or a node carries 256 or more elements.

## The solder/air sideset

In addition to the tape pieces, terminals, and three outer boundary faces, the generator
writes one more sideset, `SolderAir`: every triangle where the solder or a lead bar touches
the air, that is, the top and bottom faces of the winding, the inner face of the sliver, and
the exposed faces of both bars. The tape facets that face the air are not in it; they
belong to the thin shell. It is the last sideset, so all other ids remain unchanged
(813 with 200 turns, 94 with 20). The deck does not list it: BELFEM types a sideset between a
conductor block and the air block as a conductor/air interface on its own.

BELFEM needs it for the cuts. Before accepting a thin cut, BELFEM prunes the cut sheet's
dangling edges; only edges on an interface, boundary, or periodic sideset are exempt from
that pruning. The cut sheets here live in the air and end on the solder surface. Without a
sideset on that surface, their edges are not exempt, and pruning removes the sheet
until the cut is empty. The generator checks that every `SolderAir` triangle lies
between exactly one solder tet and one air tet, and that `SolderAir` and the air-facing tape
facets together cover the whole conductor/air contact.

## Files

| file | role |
|---|---|
| `python/main.py` | mesh generator; run it from this directory |
| `python/pancake/` | the generator package |
| `input.conf` | deck; the generator writes its `topology` and `boundary conditions` sections |
| `Allrun`, `Allclean` | wrappers around the shared launcher in `../scripts/` |

## Run

This deck has no `.geo` file. Generate the mesh from this directory, then run the deck through
the shared launcher:

```bash
python3 python/main.py                 # writes pancake.msh and pancake_topology.conf
./Allrun                               # ./Allrun --dry-run prints the plan first
./Allclean                             # removes run output, keeps the mesh
```

Run the generator from this directory so `pancake.msh` lands where `input.conf` expects it.
If the mesh is missing, `./Allrun` stops and identifies `python/main.py`; run it as shown above,
not from inside `python/`. The tool requires `numpy` and calls the `gmsh` executable. The
launcher chooses ranks and threads and can submit to Slurm; per-deck settings go in a
`run.conf` next to `input.conf` (see `../scripts/README.md`).

```bash
python3 python/main.py 20              # 20-turn smoke case at the same pitch
```

This builds a 20-turn smoke case at the same 50 µm pitch, r = 7 mm to 8 mm, in a few seconds.
Edit the parameters at the top of `main.py` to change the coil. The deck's physical ids
follow the number of turns. After changing it, or for the smoke case, paste the generated
`pancake_topology.conf` into `input.conf` in place of its `topology` and `boundary
conditions` sections. The shipped `input.conf` holds the 200-turn ids.

The tape layer stack in `input.conf` is a placeholder that fits the 50 µm pitch, and the
solder is the built-in Pb38Sn62 until a user material replaces it.

## Status

This example builds its mesh and passes cut creation, but it does not yet run to a solution.
A 20-turn run currently stops in the kernel while the dofs are set up, on a known defect: an
interface node whose air elements are all taken over by the duplicates of two cut sheets
that end on the same interface edge gets no φ dof, so the hanging interface edges that
depend on it have no dof to take their weights from. A debug build stops at an assertion in
`DofData::create_dofwise_t_matrices_master`; a release build crashes in the same function.
The 200-turn case has not been run past that point either.
