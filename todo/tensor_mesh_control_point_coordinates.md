# Tensor-Mesh Control-Point Coordinates

**Date:** 2026-10-22
**Purpose:** Decide whether the stored B-spline control-point coordinates of the tensor meshes
should follow the element-parameter convention the basis functions use.
**Module:** `src/mesh`
**AIs involved:** Codex (raised), Claude (trace); decision Christian's
**Status:** OPEN — question for Christian. Cosmetic for the projector, which never reads these
coordinates (`cl_DatabaseProjector.cpp:216-300`).

---

## 1. The Convention (Christian, 2026-10-22)

Each element is ξ ∈ [-1, 1]. Linear B-spline control points coincide with the nodes; quadratic
ones sit at ξ = -2, 0, 2, cubic ones at ξ = -3, -1, 1, 3, i.e. one element width apart and starting
half an element (quadratic) or one element (cubic) before the element. These are the Greville
abscissae, and they agree with `bspline_quad` / `bspline_cub`: for cubic at ξ = -1,
(8·(-3) + 32·(-1) + 8·1) / 48 = -1.

## 2. What the Factories Do

The node grid spans `origin` to `limit = origin + order · nel · step` (`cl_TensorMeshConfig.cpp:68`),
so `step()` is the node step and the element size is h = order · step. The control points are
spaced by `step()`:

| order | convention (ξ) | `create_bsplines_*` start, spacing | stored, in ξ |
|---|---|---|---|
| 2 | -2, 0, 2 | x0 − 0.5·step, step (`cl_TensorMeshFactory.cpp:610-617`) | -1.5, -0.5, 0.5 |
| 3 | -3, -1, 1, 3 | x0 − step, step (`:640-648`, `:710-722`) | -1.67, -1, -0.33, 0.33 |

The T-matrix and the node weights do not depend on these coordinates, and
`tests/mesh/test_TensorMeshBSpline.cpp` computes the Greville positions from the lattice index
rather than reading them.

## 3. Steps

- [ ] **R1** — Christian: is the stored position meant to follow §1? If yes, space by
  `element_step()` and start at x0 − (order − 1)·h / 2 in every `create_bsplines_*`; plan + jury.
- [ ] **R2** — if the coordinates are meant to stay as they are, document what they denote in
  `cl_TensorMeshFactory.hpp`.
