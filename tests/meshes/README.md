# Test Meshes

**Date:** 2026-08-13
**Purpose:** Checked-in mesh fixtures for the test suite.

This directory holds mesh files that tests load from disk. Decided 2026-08-13
(Christian): fixtures live here, colocated with their consumers, **not** under
`share/` — `share/` carries runtime data shipped to users (`share/material`),
and a `USE_TEST=OFF` install should not carry test fixtures.

## Conventions

- **Prefer no file at all.** A mesh built programmatically inside the test
  (see `tests/homology/test_Cohomology.cpp`, `tests/mesh` tensor meshes,
  `tests/fem/support/cl_TS_TestStack.hpp`) is easier to debug and needs no
  asset. Check a mesh in only when the loader path itself is under test or the
  geometry is impractical to build in code.
- **Check in the `.geo` alongside the generated `.msh`.** The `.geo` is the
  source of truth; the `.msh` is committed so running the tests does not
  require gmsh. Regenerate with the gmsh version noted in a comment at the top
  of the `.geo`.
- **Name by what the mesh is**, lowercase_with_underscores:
  `annulus_tri3.msh`, `periodic_strip_quad4.msh`.
- Keep meshes minimal — tens of elements, not thousands; the `fast` test label
  assumes it.
