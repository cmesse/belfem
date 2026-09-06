# Postprocessing Module Documentation {#fem_postproc_index}

**Module:** `src/fem/postproc`
**Purpose:** Index of documentation for BELFEM's mesh-level postprocessing utilities

---

## Overview

The `postproc` module provides mesh-level utilities for solved meshes. It covers
geometric measures (volume, surface, edge lengths), surface normals, scalar-field
integration over sidesets, gradient recovery, and a ghost-mesh builder for triangulated surfaces (`Surface`).

**This is not `MaxwellPostprocessor`.** That class lives in `src/fem/maxwell` and recovers
the electromagnetic fields (B, H, J, J/Jc) from the h-φ solution. That class is
physics-specific; the geometric utilities here (volume, surface, normals, edge lengths,
sideset integrals, gradient recovery) work on any mesh. `mises_planestress()` is a plane-stress
elasticity routine and `Surface` handles TRI3 sidesets only.

---

## Key Classes

| Class | Files | Role |
|-------|-------|------|
| **`Gradient`** | cl_Gradient.{hpp,cpp} | Gradient recovery; derives from `Postprocessor` |
| **`Surface`** | cl_Surface.{hpp,cpp} | Surface extraction and evaluation |

## Free Functions

| Function | Header | Computes |
|----------|--------|----------|
| `mises_planestress()` | fn_FEM_mises_planestress.hpp | von Mises stress under plane stress — stub: aborts with "missing implementation" (plane-stress C-matrix not populated) |
| `compute_volume()` | fn_Mesh_compute_volume.hpp | Volume of a block |
| `compute_surface()` | fn_Mesh_compute_surface.hpp | Area of a sideset |
| `compute_surface_normals()` | fn_Mesh_compute_surface_normals.hpp | Outward normals on a sideset |
| `compute_edge_lengths()` | fn_Mesh_compute_edge_lengths.hpp | Edge lengths |
| `integrate_scalar_over_sidesets()` | fn_Mesh_integrate_scalar_over_sidesets.hpp | Sideset integral of a scalar field |

Within `src/` the module currently has no compiled consumer; `src/fem/kernel/staticheat.cpp` includes the normals header, but its call is commented out and the driver is not built.

The module also includes two development drivers with their own `main()` functions:
`normaltest.cpp` and `pentatest.cpp` (built only under `USE_EXAMPLES`). `Doxyfile.in` excludes them from the
generated API reference with `EXCLUDE_PATTERNS`; they are not part of the library.

---

## See Also

- **Source code** - header files in `src/fem/postproc/`
- [Maxwell module](../../maxwell/doc/README.md) - `MaxwellPostprocessor`, the field recovery
- [Mesh module](../../../mesh/doc/README.md) - the mesh entities these utilities traverse
- [FEM kernel](../../kernel/doc/README.md) - the `Postprocessor` base class

---

**Status:** this index was written from the module's headers. The individual algorithms —
in particular the gradient recovery — are not yet documented.
