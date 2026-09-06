# Cos-theta iron model: B-H curve from hdf5 + coil symmetry planes

**Date:** 2026-07-02
**Purpose:** Wire up the cosine-theta quadrant magnet: iron yoke picks its nonlinear
B-H curve by name from an HDF5 database, and add the two symmetry planes with the
correct current-mirror behaviour.
**Module:** fem/maxwell, physics/materials, mesh

## Context
Modelling a **quadrant** of a cosθ dipole (`cmake-build-debug/costheta.{geo,msh,input.conf}`):
cables as `coil` blocks (1:68), iron yoke as `Ferro` blocks (69:70), air (71:74).
The iron nonlinear μ(H) comes from `bhdata.hdf5` (groups `RoxieIron`, `SAE1010`,
`SAE310@4K`). Symmetry planes carry the current mirror.

## Findings
- `BhCurve(path,label)` already selects a named group (`bnur`/`hnur`/`bsat`/`hsat`)
  from the hdf5 (cl_BhCurve.cpp:41-56) — matches bhdata.hdf5 layout. Loader was fine.
- Gap was in `MaterialFactory` input parsing: only `builtin`/`custom` paths existed;
  `iron` dead-ended at the `//todo:: add Roxie Iron here?` (cl_MaterialFactory.cpp:240).
- **Blocker (found via Codex audit):** base `Material::set_bh_curve()` is a hard
  `BELFEM_ERROR` (cl_Material.cpp:593) and was **non-virtual**; only `Metal` stores a
  curve and overrides the virtual `mu/H/dmudH_bhcurve`. So `new Material(Ferro)` +
  `load_bh_curve` would crash. Fix: instantiate `Metal` and make `set_bh_curve` virtual
  so base `load_bh_curve` dispatches into it.
- Symmetry weak form (verified): `symmetry_phi_2d` penalizes `n×∇φ`
  (mt_maxwell_symmetry.cpp:44) ⇒ H∥=0, field normal ⇒ φ=const. `air/ferro symmetry`
  → `impose_dirichlet(0)` (cl_MaxwellFactory.cpp:626-632) = **current-preserving** plane.
  `*antisymmetry` → sideset deactivated (cl_MaxwellFactory.cpp:604-610) ⇒ natural
  `B·n=0` = **current-inverting** plane.
- `domain_type(string)` had no `"… antisymmetry"` parse (en_DomainType.cpp), and the
  `Domain` ctor switch (cl_FEM_Domain.cpp) didn't classify AntiSymmetry as a sideset —
  both needed for explicit tagging.
- Codex nuance: untagged one-sided Air/Ferro boundaries auto-classify as `*AntiSymmetry`
  in `Topology` (cl_Topology.cpp:221,251) — so the outer boundary defaults to `B·n=0`.

## Changes
1. `cl_MaterialFactory.cpp` — a `curve:` key builds a `material::Iron(NaN,false)` and calls
   `load_bh_curve(new BhSplineCurve(bhfile|"bhdata.hdf5", curve))`. (Christian added the
   `Iron : Metal` class + the `iron`/`ferro` case in `create_material`, and split `BhCurve`
   into an abstract interface + concrete `BhSplineCurve`; the factory builds the concrete so
   base `Material` stays linalg-free.)
2. `cl_Material.hpp` / `.cpp` — `load_bh_curve` now takes a pre-built `const BhCurve*` (base
   no longer constructs a spline object); `set_bh_curve` virtual, `Metal` overrides. Fixed
   `create_bh_curve` to return `BhSplineCurve`. BhCurve doc corrected (HDF5, not ASCII).
3. `en_DomainType.cpp` — parse `air/buffer/ferro/conductor antisymmetry`.
4. `cl_FEM_Domain.cpp` — treat the four AntiSymmetry types as sidesets.
5. `cmake-build-debug/input.conf` — `material: iron` on the Ferro block; `air symmetry`
   {330,331,332} (x-axis, preserving); `air antisymmetry` {337,336,293,335,303,334}
   (y-axis, inverting).
6. `maxwell_usage_guide.md` §7.2 — documented "magnetic wall (anti-symmetry) = B·n=0 =
   natural BC = no term needed", vs the flux-normal symmetry plane (B×n=0) which must be
   enforced. Clarified that only `mt_maxwell_symmetry.cpp` exists because the antisymmetry
   side is natural, and that `air symmetry`'s hard `impose_dirichlet(φ=0)` currently
   shadows the general `symmetry_phi` (B×n=0, φ-free) weak form in the IWG.

## Physics mapping (matches user's rule)
- x-axis mirror preserves current polarity → φ=const → `air/ferro symmetry`.
- y-axis mirror inverts current → natural `B·n=0` → `air/ferro antisymmetry` (deactivated).

## Open / to watch (Codex Q1, medium-high)
Coil transport current + cuts are generated from terminal IDs, independent of the
symmetry BC; there is **no** automatic half/full-current correction. Input current must
be the modeled-conductor current, and no cable should be bisected by a symmetry plane.
Verify after first run. `phi=0` symmetry + bearing at node 273 is redundant-harmless.

## Handoff
Build + run handed to user (Christian builds).
