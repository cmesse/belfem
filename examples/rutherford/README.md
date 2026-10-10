# rutherford {#doc_examples_rutherford}

**Date:** 2026-09-18
**Purpose:** REBCO tapes wound around a flat Rutherford-type core, with a Python-generated mesh.

---

## What it shows

This deck models the flat, flexible REBCO cable of Romancew (2026): REBCO tape wrapped
helically around a Rutherford-type copper core. The core is a flat bar with rounded edges, a
stadium in cross-section, and the tapes follow it at a constant lay angle measured from the
cable axis. The deck winds five 4 mm tapes at 30 degrees, the most that fit on the core,
covering 89 % of its surface; the reference study wraps a single tape (`numTapes = 1`). A
translated periodic condition represents a fifth of a pitch of an infinite cable, after which
every tape continues into its neighbor. A 30 A ramp over 10 s drives the model for 15 s
with Newton, Anderson stabilization, MUMPS, and GMRES. The copper core is not modeled: its section is meshed as air, so the deck is about the
tape alone.

`python/main.py` generates the mesh, writes `rutherford.msh`, and prints the matching
`topology` section in `rutherford_topology.conf`. The parameters follow the reference study:
tape width 2, 4 or 6 mm (default 4), lay angle 30 or 45 degrees (default 30), edge radius
0.9 mm (a 1.8 mm thick core). The core width is not documented there; 12 mm is a placeholder.

## Geometry

The cable axis is z. The core cross-section is a stadium of width `cableWidth` along x and
thickness `2 edgeRadius` along y: two flat faces joined by semicircles of radius
`edgeRadius`. Its perimeter is

    P = 2 (cableWidth - 2 edgeRadius) + 2 pi edgeRadius

The extruded core is a developable surface, so a tape at lay angle θ is a straight band in
the unrolled (arc length, z) plane, advancing `tan θ` of perimeter per unit z. One turn
takes the pitch `L = P / tan θ` along z, the tape covers `tapeWidth / cos θ` of perimeter at
any z, and the tape is only bent, never stretched, in this idealization. On the round edges
the tape crosses the curvature obliquely: its centerline sees the bending radius
`edgeRadius / sin² θ` and the twist rate `sin θ cos θ / edgeRadius` (Euler's theorem on the
cylinder), so a 30 degree lay bends the tape half as much as a 45 degree lay. The generator
prints both numbers.

Multiple tapes per layer sit at equal perimeter spacing. Stacked layers are counter-wound
(criss-cross) at the same pitch, which makes their lay angles differ slightly because the
outer stadium has the longer perimeter; `layAngle` is the layer-0 value. The end planes are
pure translations of each other, so `numPitches * numTapes` must be a whole number (a third
of a pitch with three tapes closes with the tapes permuted). The
half-pitch closure by a 180 degree rotation about the axis, which the stadium's symmetry
would allow, is not implemented.

The closure decides how many conductors the model has. Because tape k continues into tape
k + q across the periodic planes, with q the number of tape positions advanced over the
model length, a layer of n tapes is gcd(q, n) galvanically separate conductors. With five
tapes over a fifth of a pitch that is one conductor: a single tape of five times the length,
crossing each end plane five times. Every conductor is one current condition and consumes
one cohomology generator, and a periodic model spends a further generator on the free axial
loop. The condition is defined on **one crossing** of the conductor, the back-plane and
front-plane curves of its first tape (`[1]` and `[2]`), with the tape current as amplitude:
the current is the same on every crossing, so one crossing fixes it, exactly as
`corc_periodic_bc` does. Listing several crossings of one conductor in one bracket group is a
different condition (the group chain links the conductor once per crossing) and is not what
is wanted here.

## Mesh

The tape strips are built from rows at constant z, spaced so the step along the tape equals
`tapeResolution`; the slots across the tape are the same size. A strip quad is therefore a
parallelogram with one side along the perimeter and one along the tape, and the periodic
cut planes are rows of the mesh rather than slices through elements. Split along the short
diagonal, the quads give near-equilateral triangles at 30 degrees. The bare core between the
tape turns is meshed at `gapResolution`. The flat/round junctions cross these elements; the
chord error is a few micrometers at the default resolution. End caps are meshed once in the
plane through the exact ring positions and instantiated at both ends, so the periodic pairs
match node for node. With the defaults, the mesh has about 15 000 nodes and 74 000
tetrahedra and takes twenty seconds to build.

## Files

| file | role |
|---|---|
| `python/main.py` | mesh generator; run it from this directory |
| `python/rutherford/`, `python/mesh/` | packages the generator uses |
| `python/rutherford/tests/` | unit and end-to-end tests for the generator |
| `input.conf` | deck |

## Run

This deck ships no `.geo`. The mesh is built by a Python tool, which also prints the matching
`topology` section for `input.conf`:

```bash
python3 python/main.py                 # writes the .msh into the current directory
../../build/bin/belfem                 # or mpirun -np 4 ../../build/bin/belfem
```

Run it from this directory so the mesh lands where `input.conf` expects it. The tool needs
`numpy` and calls the `gmsh` executable. The tests run from `python/`:

```bash
cd python
python3 rutherford/tests/test_stadium.py
python3 rutherford/tests/test_winding.py
python3 rutherford/tests/test_e2e.py        # calls gmsh, about half a minute
```

## Reference

Romancew, E. (2026). *Development and Characterization of Flat and Flexible REBCO Cables for
High-Field Accelerator Magnets.* M.S. thesis, Northern Illinois University. Tape widths, lay
angles and the bending study; the core cross-section is not given there.
