# PENTA6 facet orientation notebooks

**Date:** 2026-09-16
**Purpose:** The interactive studies behind the PENTA6 slave-face parameter maps and the per-face outward normals in `fn_IF_initialize_integration_points_on_facet.cpp`.
**Module:** `src/fem/interpolation`

Both scripts are plot notebooks meant to be run in the MATLAB desktop, one case at a time; neither emits a C++ table and neither is run headless.

| Script | What it shows | How to use it |
|---|---|---|
| `orientation.m` | the five EXODUS faces of a PENTA6, every slave ordering of each face (`S1a … S5c`), and for each `(face, orientation)` the map `eta = f(xi)` from a point on the master face onto the slave element | set `S` to one slave ordering, keep only the matching `eta` block (each later block overwrites the earlier one), run; the master point `p` (red cross) and the slave point `q` (blue circle) coincide when the map is right |
| `facet.m` | for face `f`, the outward normal at the face center, built from two rows of the Jacobian `dN/dxi * X`; the `if` chain records which rows, with which signs, give the outward normal on each face | set `f`, run; the arrow must point out of the element |

The runtime check of the same conventions is `tests/fem/test_FacetIntegrationPoints.cpp`, which maps quadrature points through every slave orientation on real coordinates. `orientation.m` needs no toolbox; `facet.m` builds its Jacobian from `syms` and needs the Symbolic Math Toolbox or the Octave `symbolic` package.
