# Facet Normals and Pipette: Restore Test Coverage

> **CLOSED 2026-09-03** (todo/ currentness sweep, round 3): folded into the test-hardening campaign plan as a single coverage row. That plan was removed from the tree on 2026-09-05, so the gap (no `test_Normals.cpp`, no Pipette test) is real and currently untracked. Status lines and checkboxes below are as they stood at closure and are not maintained.

**Date:** 2026-08-05
**Purpose:** Rebuild the facet-normal and volume-integral (pipette) tests that were
lost when the `geometry::` namespace was retired, using the surviving reference
fixture.
**Module:** `src/fem/kernel`

---

**Status:** OPEN — gap identified during the `tests/old` triage (2026-08-05). The
old tests were deleted rather than revived because they exercise an API that no
longer exists; the reference data they consumed was kept.
**Re-verified 2026-08-09 (currentness sweep): unchanged and still fully valid.**
`tests/fem/` contains `test_AndersonMixing`, `test_FacetIntegrationPoints`,
`test_Integration`, `test_InterfaceOrientation` and `test_LagrangeInterpolation` — there is
still **no** `test_Normals.cpp` and no pipette coverage, so `Calculator::normal_*` and
`cl_Pipette` remain at zero tests. The fixture `tests/fem/support/test_database.hdf5` is
still present. R1 (confirming the outward-normal convention) is the blocking decision and it
is Christian's.

---

## 1. What happened

`tests/old/fem/` held nine data-driven tests: eight `fn_normal_<type>.cpp` files
covering TRI3/TRI6, QUAD4/QUAD9, TET4/TET10 and HEX8/HEX27, plus
`cl_IntegrationData_Interface.cpp`. All of them were written against a
`FEM_geometry.hpp` header that exposed a `belfem::fem::geometry` namespace
(`normal_tri`, `normal_quad`, `normal_tet`, `normal_hex`, `test_pipette`,
`test_intpoints`, `create_test_mesh`). That header and namespace have been
removed from `src/`; the normal computation now lives as member functions on the
calculator (`src/fem/kernel/cl_FEM_Calculator.hpp:1249-1284` —
`normal_tri_straight`, `normal_tri_curved`, `normal_quad_straight`,
`normal_quad_curved`, `normal_hex`).

The tests were therefore not revivable, only rewritable, and were deleted. An
earlier triage (`todo/closed/tests/existing_tests_triage.md`, row for
`test/fem/fn_normal_*.cpp`) had marked them "keep" — that verdict predates the
header removal and no longer holds.

## 2. Current coverage gap

| Code under test | Location | Covered today? |
|---|---|---|
| `Calculator::normal_tri_straight` / `_curved` | `cl_FEM_Calculator.hpp:1249-1257` | no |
| `Calculator::normal_quad_straight` / `_curved` | `cl_FEM_Calculator.hpp:1259-1267` | no |
| `Calculator::normal_hex` | `cl_FEM_Calculator.hpp:1284` | no |
| `Pipette` (element volumes / surfaces) | `src/fem/kernel/cl_Pipette.{hpp,cpp}` | no |

Adjacent and already covered: `tests/fem/test_FacetIntegrationPoints.cpp` maps
facet integration points from volume parametric space through every orientation
permutation, which is what `cl_IntegrationData_Interface.cpp` was for. That one
does **not** need to come back.

## 3. Available fixture

`tests/fem/support/test_database.hdf5` (64 KB, moved out of `tests/old/fem/`)
carries per-element-type reference data under the `normals` group: `Nodes`,
`Normals`, `Surface`/`Surfaces`, `Volume`, and for the 3D types also `Points`
and `Weights`. The data is independent of the deleted API and is still valid as
ground truth. The `interfaces` group in the same file fed `create_test_mesh` and
is not needed unless the interface test is rebuilt.

## 4. Steps

- [ ] **R1** — Confirm the sign/orientation convention of the `Calculator`
      normals matches the stored `Normals` columns (outward positive). If it
      does not, establish which is authoritative before writing expectations —
      this is a physics call, not a test-style call.
- [ ] **R2** — Add `tests/fem/test_Normals.cpp` driven by
      `tests/fem/support/test_database.hdf5`, covering TRI3/TRI6, QUAD4/QUAD9,
      TET4/TET10, HEX8/HEX27 against the `Normals` reference columns.
- [ ] **R3** — Extend to surface and volume checks (the old `test_pipette` role)
      against the stored `Surface`/`Surfaces` and `Volume` entries, exercising
      `Pipette` directly.
- [ ] **R4** — Wire the fixture into `tests/fem/CMakeLists.txt`. The old
      `file(COPY ... DESTINATION /tmp)` pattern with a global `gDatabase`
      pointer should not be reproduced; open the file from a path known to the
      test instead, so the suite stays runnable in parallel.
- [ ] **R5** — Label the test `fast` if it stays under the ~10 s per-test
      budget.

## 5. Out of scope

- Reviving `cl_IntegrationData_Interface.cpp` — superseded by
  `tests/fem/test_FacetIntegrationPoints.cpp`.
- Reviving `tests/old/maxwell/` — those drove `fem::MaxwellJob::initialize_test`
  / `run_test` and per-element `IwgType::MAXWELL_HPHI_*` enum entries, none of
  which exist any more (the enum collapsed to a single `IwgType::Maxwell`,
  `src/fem/iwg/en_IWGs.hpp:73`). An end-to-end h-φ regression test is worth
  having, but it is a new design, not a restoration, and belongs in its own
  plan.
