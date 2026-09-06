# DR-153: MeshChecker Orphan Deleted — and It Was Holding a Real Fix Hostage

**Date:** 2026-08-31
**Purpose:** Close DR-153 (orphaned `src/fem/postproc/cl_MeshChecker.{hpp,cpp}`); port the pyramid
dispatch fix the orphan was carrying before deleting it.
**Module:** fem/kernel, fem/postproc

## The row understated itself

DR-153 was filed as a maintenance/review hazard: a second copy of `MeshChecker` that no
`CMakeLists.txt` compiles, which had already misled a vendor in the DR-111 round. Investigating
before closure showed the orphan was worse than a reading trap — **it had absorbed a functional
bug fix.** Commit `d205680a` ("fixes for pyramid", 2026-03-12) changed the PYRA13/PYRA14 dispatch
from `swap_pyra5` to `swap_pyra14` in the postproc copy only. The compiled kernel copy still
dispatched both types to `swap_pyra5`, so on a clockwise gmsh mesh carrying second-order pyramids
the live checker reverses the base corners while leaving the mid-edge nodes bound to the pre-swap
edges — a corrupted element rather than a reoriented one. This made the orphan's capture count
three: INC-185 (doc fix on the dead header), the DR-111 vendor cite, and `d205680a` itself.

## The derivation that confirmed the March fix

Rather than trusting `d205680a` blind, the swap set was derived from first principles, with the
known-good `swap_pyra5` (swap nodes 0,2) as the anchor:

1. Corners 0↔2 with 1, 3, 4 fixed is, in the `cl_IF_PYRA14.hpp` parameter coordinates, the
   reflection **(ξ,η,ζ) → (−η,−ξ,ζ)** — Jacobian determinant −1, orientation-reversing, exactly
   what a negative-volume repair needs. Since PYRA5 works, this reflection is the canonical repair.
2. Pushing all 14 node coordinates through the map: slots 1, 3, 4, 10, 12, 13 are fixed points;
   the swap pairs are **{0↔2, 5↔6, 7↔8, 9↔11}** — precisely the existing `swap_pyra14` body
   (`cl_MeshChecker.hpp:333-340`).
3. Cross-checked against the `cl_Element_PYRA14.hpp` edge topology: every relabeled edge lands on
   an existing edge carrying its own midnode (e.g. edge 0 = slots (0,1,mid 5) maps onto old edge
   1 = (1,2,mid 6); edge 4 = (0,4,mid 9) onto old edge 6 = (2,4,mid 11)).
4. **PYRA13 rides along:** its parameter coordinates are identical to PYRA14 on nodes 0-12
   (checked in `cl_IF_PYRA13.hpp`), and the swap touches no index above 11, so one dispatch entry
   serves both.

The independent derivation and Christian's March fix converge, so the port is a one-token change:
`src/fem/kernel/cl_MeshChecker.cpp`, PYRA13/PYRA14 case, `swap_pyra5` → `swap_pyra14`.

## What landed

- `src/fem/kernel/cl_MeshChecker.cpp`: PYRA13/PYRA14 now dispatch `swap_pyra14` (PYRA5 unchanged).
- `git rm src/fem/postproc/cl_MeshChecker.{hpp,cpp}` — kernel copy ruled canonical; it is a strict
  superset (parallel/edges/faces guards, the INC-516 non-silent flip message,
  `set_mesh_checker_flag()`, the TS element types).
- `src/fem/postproc/doc/README.md`: `MeshChecker` row removed from the class table. The
  `@ingroup grp_fem_postproc` tags existed only on the orphan header, so the class had been
  documenting itself into the postproc doxygen group from a file nothing compiles; it drops out
  of that group with the deletion.
- A latent trap died with the orphan: both headers shared the include guard `CL_MESHCHECKER_HPP`
  while kernel and postproc are both on the include path (`src/fem/kernel/CMakeLists.txt:47-48`,
  `src/fem/maxwell/CMakeLists.txt:25-26`); kernel happened to precede postproc, so the right
  definition won by ordering luck alone.
- Register: DR-153 struck and archived to `todo/debt_register_closed.md`; `[W]` recount 4 → 3;
  `check_doc_claims.py` 37/37.

## Verification status

The DR-153 fix criterion RAN: `grep -rn cl_MeshChecker.cpp --include=CMakeLists.txt src/` returns
exactly `src/fem/kernel/CMakeLists.txt:29`, and the tree holds exactly one MeshChecker pair. All
remaining consumers (`cl_FEM_Kernel.cpp:28`, `cl_MaxwellFactory.cpp:45`) resolve uniquely.

The ported dispatch is **reviewed, not verified**: no CW second-order-pyramid deck exists and
`make check` does not reach the checker's pyramid path — a rebuild plus suite run is a regression
gate only. The port landed on Christian's direct instruction with the in-session parametric
derivation as its review; no vendor round was run. `swap_pyra14`'s body itself is unchanged from
what has sat in both headers since `e5f533ae`.
