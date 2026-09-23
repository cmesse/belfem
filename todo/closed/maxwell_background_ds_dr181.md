# DR-181: `background_phi` Reads the Surface Increment at Point 0 for Every Point

**Date:** 2026-09-22
**Purpose:** `background_phi` integrates the background-field sideset kernel with `aCalc->dS()`
(defaulted index 0) inside its loop over integration points, so on any facet whose surface
increment varies across the face — curved second-order facets, HEX/PENTA faces — K and f carry
the point-0 increment at every point. The fix is two tokens (`dS()` → `dS( k )`); the plan's
weight is in proving the mechanism, sweeping for siblings, and building a fixture that is red on
the current code.
**Module:** `src/fem/maxwell` (kernel `matrices/mt_maxwell_background.cpp`); reads
`src/fem/kernel/cl_FEM_Calculator`; gate in `tests/fem`
**AIs involved:** Claude (exploration + plan), Codex (plan/code audit — critical-framework row,
`gpt-6-astra`), Grok (third voice)
**Status:** CLOSED 2026-09-23 — landed in the combined commit with the mesh bundle and DR-188. Plan jury (Grok on its third attempt) and code jury folded in; `BackgroundField.*` ran red on the unfixed kernel (K 17.6 %, f 18.8 %) and green after; DR-181 struck and archived, by-catch DR-183 refuted. Devlog `dl20260922_debt_register_oneliners_dr181_dr188.md`.

> **Scope guards (from the task brief):**
> - IN: the two `dS()` calls in `background_phi`; any other defaulted-index `dS()`/`dV()` inside a
>   loop over integration points found by the sweep (extend-don't-branch — none found, see §1.2);
>   one fixture that is red on the current code.
> - OUT: the Calculator's normal/increment caching design (`mNormalIndex`), the `0.5`/`0.25`
>   chain-rule factors already pinned by DR-178, the narration comments in the kernel, the
>   `Dirichlet → BackgroundField` stamping observed in §3.3 (reported, not touched).
> - No deck, no input key, no documentation contract changes.

---

## 1. Current Behaviour and How It Fails

### 1.1 The mechanism, traced through the Calculator

`background_phi` (`mt_maxwell_background.cpp:43-50`):

```cpp
for ( uint k=0; k<aCalc->num_intpoints(); ++k )
{
    const Matrix< real > & B = aCalc->Bm( k );
    aMatrices->K() += w( k ) * trans( B ) * B * aCalc->dS();
    aMatrices->f()   -= w( k ) * trans( B ) * H * aCalc->dS();
}
```

`dS( const uint aIndex = 0 )` (`cl_FEM_Calculator.hpp:1125`, body `:2079-2089`) forwards to
`mFundS`; for a Cartesian model that is `dS_cartesian` (`:2265-2274`), which **calls the normal
callback at `aIndex` and returns the member `mSurfaceIncrement`**. So `dS( k )` is valid at any
point and needs no prior call — the normal callback is the thing that computes the increment, and
`dS` invokes it. Whether the returned value depends on `k` is decided by the callback's caching
rule (`cl_FEM_Calculator.cpp`):

| master geometry / order | callback | recompute rule | increment per point? |
|---|---|---|---|
| TRI, linear | `normal_tri_straight` `:2245-2289` | once, when `mNormalIndex == BELFEM_UINT_MAX` (`:2247`) | no — constant edge length |
| QUAD, linear (2-D) | `normal_quad_straight` `:2355-2403` | once (`:2357`) | no — constant edge length |
| TET, linear | `normal_tet_straight` `:2466-2517` | once (`:2468`) | no — flat triangle |
| TRI / QUAD, quadratic, `mIsCurved` | `normal_tri_curved` `:2294-2350`, `normal_quad_curved` `:2408-2461` | when `mNormalIndex != aIndex` (`:2300`, `:2414`), from `Jm( aIndex )` | **yes** |
| TET, quadratic, `mIsCurved` | `normal_tet_curved` `:2523-2580` | when `mNormalIndex != aIndex` (`:2529`) | **yes** |
| HEX (any order) | `normal_hex` `:2692-2761` | when `mNormalIndex != aIndex` (`:2694`), from `Jm( aIndex )` | **yes** — a bilinear quad face has a constant Jacobian only if it is a parallelogram |
| PENTA (any) | `normal_penta` `:2585-2681` | when `mNormalIndex != aIndex` (`:2587`) | **yes** on tapered side faces |

Selection is by `mesh::geometry_type( mGroup->master_type() )` and interpolation order
(`:1536-1590`). `link()` resets `mNormalIndex = BELFEM_UINT_MAX` and `mSurfaceIncrement = NaN`
(`:1751`, `:1888`).

The failing trace on a per-point family (high confidence, source-traced):

1. `k = 0`: `Bm( 0 )` → `Jm( 0 )`; `dS()` → `normal_X( 0 )`: `mNormalIndex` is `UINT_MAX` ≠ 0 →
   computes the point-0 normal and increment, sets `mNormalIndex = 0`.
2. `k ≥ 1`: `Bm( k )` → `Jm( k )` (the Jacobian cache moves to `k`); `dS()` → `normal_X( 0 )`:
   `mNormalIndex == 0 == aIndex` → **cache hit, no recompute** → returns the point-0 increment.

So `K = Σ_k w_k B_kᵀ B_k · dS_0` and `f = −Σ_k w_k B_kᵀ H · dS_0` instead of `dS_k`. On the
"once" families rows 1-3 the two are identical to the bit, which is why the flat decks are exact.

| Failure | Mechanism | Evidence |
|---|---|---|
| Wrong K and f on curved / tapered background facets | point-0 increment reused at every k | `mt_maxwell_background.cpp:48-49`; `cl_FEM_Calculator.hpp:2265-2274`; `.cpp:2694`, `:2300`, `:2414`, `:2529`, `:2587` |
| Latent: wrong even on a straight edge under an axisymmetric model | `dS_axsymmx/y` (`hpp:2279-2301`) multiply by the radius interpolated **at `aIndex`**, so `dS()` is the radius at point 0 | not reachable by Maxwell today: `AxSymm` has 0 hits in `cl_MaxwellFactory.cpp` and `cl_IWG_Maxwell.cpp`, 0 axisymmetric example decks (`rg -i axisym examples`) |
| Exact on flat TRI3 / straight LINE2 facets | callbacks rows 1-3 compute once | `:2247`, `:2357`, `:2468` |

**Bottom line:** the surface increment is per-point wherever the callback is per-point, `dS( k )`
is the correct and cheap call (the Jacobian at `k` is already cached by `Bm( k )`,
`Bscalar_master` `hpp:1905-1918`), and the defect is exactly the two defaulted calls.

### 1.2 Sweep: every defaulted-index `dS()` / `dV()` in the tree

Scope: `rg -uu` over `src/` and `tests/` for `->dS( )`, `->dV( )`, `.dS( )`, `.dV( )`, plus the
explicit-index forms `dS( 0 )`, `dV( 0 )` and the sibling pattern `->normal( )`. `nonfree/` is
**absent in this checkout** (`test -d nonfree` → absent) and is therefore **not swept** — a
session with the tree present owes the same `rg -uu` pass before this row is closed.

| # | site | inside an integration-point loop? | verdict |
|---|---|---|---|
| 1 | `mt_maxwell_background.cpp:48` `dS()` | yes (`k`) | **defect — this plan (R1)** |
| 2 | `mt_maxwell_background.cpp:49` `dS()` | yes (`k`) | **defect — this plan (R1)** |
| 3 | `mt_maxwell_phi.cpp:164` `phi_tri3`: `dV()` | no — single `B( 0 )`, closed-form TRI3 mass | correct: linear simplex, constant Jacobian |
| 4 | `mt_maxwell_phi.cpp:176` `phi_tet4`: `dV()` | no — same for TET4 | correct |
| 5 | `mt_maxwell_phi.cpp:204` `phi_tri6_tet10`: `dV()` after the loop | no — the non-curved branch only (`:184`); `dV_tri6_tet10` pins the index to 0 when `!mIsCurved` (`hpp:2108-2118`), so `dV( k )` would return the same value | correct by construction; the curved branch goes through `phi()` |
| 6 | `cl_MaxwellPostprocessor.cpp:551`, `:864` `tCalc->normal()` | no — outside any loop, comment states "we are always having linear elements" | not the pattern; a postprocessing-only linear assumption, noted for the register, not widened here |
| 7 | `dS( 0 )` / `dV( 0 )` explicit | — | 0 hits outside `cl_FEM_Calculator` |

Result: **the sweep widens nothing.** The plan stays at the two tokens.

---

## 2. Architecture: Why the Two-Token Fix Is the Whole Fix

`dS( k )` is what every sibling sideset kernel already does: `mt_maxwell_symmetry.cpp:54, 85,
111, 140`, `mt_maxwell_h.cpp:570` (`real wdS = w( k ) * aCalc->dS( k )`),
`cl_IWG_StaticHeatConduction.cpp:146-191`. The Calculator's contract is that the caller asks by
index; the default exists for single-point closed-form kernels (rows 3-4 above), not for loops.

Rejected alternatives:
- **Changing the default or removing it** from `dS`/`dV` — touches the critical kernel and the
  five legitimate single-point callers; out of proportion for a two-token defect, and the
  default is the right tool for `phi_tri3`/`phi_tet4`.
- **Hoisting `w( k ) * dS( k )` into a `wdS` temporary** as `h_ghost` does — fine, but a second
  `dS( k )` is a cache hit (`mNormalIndex == k`), and the house rule for bug fixes is the value
  in place, not a rename. Offered to the jury as optional, not recommended.

---

## 3. Gap Table and Consumer Impact

| # | State | Needed for | Handled today? | Class | Citation / rationale |
|---|---|---|---|---|---|
| 1 | per-point `dS` in `background_phi` K | correct stiffness on curved/tapered facets | no | (c) R1 | `mt_maxwell_background.cpp:48` |
| 2 | per-point `dS` in `background_phi` f | correct load on curved/tapered facets | no | (c) R1 | `:49` |
| 3 | a fixture red on the current code | the gate | no | (c) R2 | §4 R2 |
| 4 | flat-facet decks unchanged | consumer promise | yes — bit-identical, see §3.2 | (a) | `:2247`, `:2357`, `:2468` |
| 5 | `nonfree/` siblings | closing the row | unknown — tree absent | (b) O3 | §1.2 |

### 3.1 Who reaches `background_phi`

`IWG_Maxwell::link_to_group` binds `mFunMKF = maxwell::background_phi` for
`DomainType::BackgroundField` (`cl_IWG_Maxwell.cpp:421-424`), for any facet dimension. A group
gets that type from:

- a `background` boundary condition: the mesh sideset (`cl_MaxwellFactory.cpp:175-179`) and then
  the dof-manager group with the enrichment integration order (`:865-870`);
- a `backgrounddirichlet` or `dirichlet` boundary condition: the **dof-manager** sideset is
  stamped `BackgroundField` while its nodes are pinned (`:908-912`) — see §3.3. The
  `backgrounddirichlet` arm is **dormant handling**: the name parses
  (`en_FEM_BoundaryConditionType.cpp:103-105`), but `MaxwellBoundaryConditionFactory` refuses it —
  its Background/Neumann/Dirichlet case excludes it and it falls to the unsupported-type error
  (`cl_MaxwellBoundaryConditionFactory.cpp:191-214`). Only `dirichlet` is deck-reachable.

### 3.2 Shipped decks

`rg -n -i "background|dirichlet" examples --glob input.conf` → exactly two hits:
`examples/disk_pulse/input.conf:78` and `examples/disk_ramp/input.conf:77`, both
`background : bgfield { sidesets : 7:14 ; direction : 0, 0, 1 }`. Sidesets 7-14 are the eight
spherical patches of the 200 mm air sphere (`disk.geo:69-80` `Circle(13..24)`, `Surface(7..14)`
`:83-97`) — the **surface is curved**, but the mesh is first order (`disk.geo` has no
`SetOrder`/`Mesh.ElementOrder`; 0 hits for `order`), so the facets are flat TRI3 under TET4
masters → `normal_tet_straight` → the increment is computed once per facet →
**`dS()` and `dS( k )` return the same double, and both decks are bit-identical before and
after the fix — conditional on `disk.msh` having TET4 masters on sidesets 7-14.** Both decks load
`disk.msh` (`disk_ramp/input.conf:8`, `disk_pulse/input.conf:8`), and the first-order claim is
inferred from `disk.geo` alone; no mesh file sits next to the decks. Source-traced for TET4
masters (`cl_FEM_Calculator.cpp:2468`); **not executed** end to end.

Result change on shipped decks: **none, under the TET4-master condition above.** The same decks re-meshed at order 2 would flag every
sphere-facet TET10 curved (`cl_Mesh_CurvedElementChecker.cpp:648` `check_tet10`; midside nodes on
a sphere are off the chord), and then the defect would be live on all 8 sidesets.

### 3.3 Observation outside scope (reported, not fixed)

A `dirichlet` boundary condition stamps `DomainType::BackgroundField` onto the dof-manager
sideset (`cl_MaxwellFactory.cpp:908-912`). The draft of this plan inferred that such a sideset
would be bound to `background_phi` and assemble `∫ BᵀB` and `−∫ BᵀH` "on rows that are fixed
anyway". Both halves fell in the plan jury:

- *"rows that are fixed anyway" is wrong* (astra): the factory pins the sideset **nodes**, but a
  boundary facet links **all master nodes** (`cl_FEM_Element.cpp:578-582`), so a `background_phi`
  contribution would not be confined to pinned rows.
- *the assembly premise is refuted on the controller path* (Grok, source trace, re-verified by
  Claude): the controller runs `compute_boundary_conditions` (`cl_FEM_Controller.cpp:282`) before
  `compute_jacobian_and_rhs` (`:1090`); `impose_bc` → `SideSet::impose_dirichlet`
  (`cl_FEM_PhysicalBoundaryCondition.cpp:262-266`) resets the group to `DomainType::Dirichlet`
  and `GroupActivationMode::GeometryOnly` on every step (`cl_FEM_SideSet.cpp:365-368`);
  `is_active()` is true only for `GeometryAndDofs` (`cl_FEM_Group.hpp:607-609`), and the sideset
  loop in `compute_jacobian_and_rhs` skips inactive sets (`cl_FEM_DofManager.cpp:614-618`). So a
  Dirichlet sideset never reaches `link_to_group` as `BackgroundField` during the timestep solve.
  What remains is a stale domain-type write that the next `impose_dirichlet` overwrites; the window
  between the factory stamp and the first `impose_bc` was not searched beyond the controller (low
  residual).

Register row **DR-183** was filed from the draft's premise before the refutation landed; its
disposition is Christian's ruling (O2). This plan does not edit the register.

---

## 4. Ordered Steps

### R1 — the fix (`src/fem/maxwell/matrices/mt_maxwell_background.cpp:48-49`)

- [x] Replace the two defaulted calls with the loop index. Intended diff:

```diff
                 for ( uint k=0; k<aCalc->num_intpoints(); ++k )
                 {
                     const Matrix< real > & B = aCalc->Bm( k );
-                    aMatrices->K() += w( k ) * trans( B ) * B * aCalc->dS();
-                    aMatrices->f()   -= w( k ) * trans( B ) * H * aCalc->dS();
+                    aMatrices->K() += w( k ) * trans( B ) * B * aCalc->dS( k );
+                    aMatrices->f()   -= w( k ) * trans( B ) * H * aCalc->dS( k );
                 }
```

What could break: nothing on the "once" families (identical value); on per-point families the
second `dS( k )` is a cache hit (`mNormalIndex == k`), and `normal_X( k )` reads `Jm( k )`, which
`Bm( k )` has just cached (`Bscalar_master`, `hpp:1905-1918`) — no extra Jacobian evaluation. The
debug assert `adS >= 0.0` (`hpp:2083`) now fires per point rather than once; a negative increment
at a later point would have been silent before. No comment change is owed — the loop states its
own contract now.

Confidence: **high** that this is the complete fix (source-traced in §1.1 and §1.2).

### R2 — the gate: a background facet whose increment varies (after: nothing; lands with R1)

- [x] New `tests/fem/test_BackgroundField.cpp`, registered under the existing
  `if( USE_MAXWELL )` block of `tests/fem/CMakeLists.txt:44-52` (`list( APPEND SOURCES
  test_BackgroundField.cpp )`); `test_fem` already carries `TESTLABELS fast` (`:2`), so it rides
  `make check-fast`.
- [x] Fixture `BackgroundHexStack` (anonymous namespace in the test file, modeled on
  `TS_TestGhostStack`, `tests/fem/support/cl_TS_TestStack.hpp:631-925`, minus thin shell and
  edges):
  - `Mesh( 3, 0, false )`; one **HEX8** air element, id 1, block 1 (`mesh::Block( 1, 1 )`,
    `DomainType::Air`). Bottom face at z = 0 on the deliberately non-parallelogram corners
    `(0,0) (1.3,0.1) (0.9,1.1) (0,0.8)` (the `TS_TestPrism` footprint, `:342-343`), top face the
    same xy at z = 0.5; counterclockwise from +z, top above bottom → right-handed.
  - A QUAD4 facet element (id 100) on the bottom-face nodes in HEX8 facet-4 order **0, 3, 2, 1**
    (`cl_Element_HEX8.hpp:84-89`), `mesh::Facet`, `set_master( tHex, 4 )` —
    `mesh::bottom_facet_index( HEX8 ) == 4` (`meshtools.cpp:708-709`), the `normal_hex` case-4 arm
    (`cl_FEM_Calculator.cpp:2734-2741`); **no slave**; `set_sideset_id( 10 )`;
    `mesh::SideSet( 10, 1 )` with `DomainType::BackgroundField`; `finalize()`; `MeshChecker`
    before the Kernel (fixture `:829-835`).
  - `KernelParameters`, `Kernel`, `IWG_Maxwell( Formulation::HPhi, ModelDimensionality::ThreeD,
    false, false )`, `set_blocks( {1}, {Air} )`, `set_sidesets( {10}, {BackgroundField} )`
    (`cl_IWG.cpp:1603-1613`), `add_equation`, `create_field`, `create_fields`,
    `set_solver( gDefaultSolver )`, `initialize()`, then `link_to_group( sideset( 10 ) )`, which
    binds `background_phi` (`cl_IWG_Maxwell.cpp:421-424`). The `"h"` vector the kernel reads is
    **not** created by `link_to_group`: it comes from `create_custom_vectors_and_matrices`
    (`cl_IWG_Maxwell.cpp:453-473`), run from `Calculator::allocate_memory`
    (`cl_FEM_Calculator.cpp:1085`) during sideset-calculator setup.
  - **Boundary-condition ownership (plan-jury P1, astra + Grok).** `add_boundary_condition` only
    stores the pointer (`cl_FEM_Kernel.cpp:743-746`), and `Kernel::~Kernel` deletes every
    registered BC (`:158-161`) — the draft's window `:124-150` stopped short of that loop, and its
    stack-allocated `tBC( 0 )` would have been an invalid delete. The fixture therefore allocates
    `new PhysicalBoundaryCondition( 0 )`: `set_type( Background )`, `set_direction( 0, 0, 1 )`,
    `set_domains( {10} )`, `fix( 1.0 )` (`value()` returns `mValue`, `mScale` defaults to 1,
    `cl_FEM_PhysicalBoundaryCondition.cpp:134-145`), `add_boundary_condition`, and keeps only a
    borrowed pointer it never deletes. The split follows the ghost fixture (`:874-923`): the kernel
    deletes the dof manager, the IWG and the BCs; the fixture deletes the kernel, the parameters
    and the mesh (`mOwnMesh` / `mOwnParameters` default false, `cl_FEM_Kernel.hpp:80`, `:97`).
    `Kernel::boundary_condition( id )` is a linear search over the BCs' domains
    (`cl_FEM_Kernel.hpp:397-411`).
  - Link through the **element** path, `tCalc->link( tGroup->elements()( 0 ) )`. For a slave-less
    facet it sets up the master at `cl_FEM_Calculator.cpp:1772-1776` and treats the slave arm
    (`:1798-1809`) as optional; `link( mesh::Facet * )` asserts a slave (`:1848-1849`) and is not
    usable for a boundary facet.
- [x] Part A (`BackgroundField.SurfaceIncrementVariesOnHexFace`) — the fixture discriminates
  (L-09: prove the path can see the defect): `Σ_k w(k) dS(k)` equals the **absolute** shoelace
  area of the bottom quad to 1e-12 relative, **and** `max_k |dS(k) − dS(0)| > 1e-3 · dS(0)` — if
  this ever fails, the face became a parallelogram and the gate is blind. The absolute value is
  required: facet 4 runs 0,3,2,1, the reverse of the counterclockwise corner cycle, so a signed
  shoelace is negative, while `dS` returns a norm. The default rule here is the **8-point
  `gauss_quad8`** (linear element → integration order 4, `fn_intpoints_auto_integration_order.cpp:28-30`;
  GAUSS scheme default, `cl_FEM_KernelParameters.hpp:53`; `fn_intpoints.cpp:161-165`), with
  Σw = 4. On a planar bilinear face `|X_ξ × X_η|` is linear in (ξ, η), so any such rule
  integrates it exactly. The draft's reading of `normal_hex`'s "0.25 × the surface" comment
  (`:2756`) as a constant `dS = area/4` holds only for a parallelogram and is not what Part A pins.
- [x] Part B (`BackgroundField.KernelUsesPerPointSurfaceIncrement`) — the red gate:
  `TimestepMatrices tMat; set_flag( K ); set_flag( F ); initialize( 8 )`;
  `fem::maxwell::background_phi( tCalc, &tMat )`; relink, then a reference assembled in the test
  from the same calculator, `Kref += w(k) · Bm(k)ᵀ Bm(k) · dS(k)`,
  `fref −= w(k) · Bm(k)ᵀ H · dS(k)` with `H = value · direction = (0,0,1)`;
  `EXPECT max|K − Kref| < 1e-12 · max|Kref|`, same for f. Prediction, pre-registered: **red on the
  unfixed kernel, green after R1.** The draft's 1e-2 to 1e-1 deviation band was a guess and is
  not claimed.
- [ ] Both trees (debug asserts the increment sign per point; release exercises the `#else`
  branch of `dS`, `hpp:2085-2087`): `test_fem --gtest_filter='BackgroundField.*'`, red with R1
  reverted and green with it, then `make check-fast`. **Owed — Christian runs builds.**
- [x] Code-round fix pass (E1–E3), test file only, approved by Christian: (E1) both TESTs skip
  before building the fixture when `comm_size() != 1` (`commtools.hpp`, the
  `test_SideSetRebuild.cpp` guard); (E2) every value entering a `std::max` reduction — the
  per-point `dS` of Part A, every `K`/`Kref` and `f`/`fref` entry of Part B — is
  `ASSERT_TRUE( std::isfinite( … ) )` first, since `std::max( a, NaN )` returns `a`; (E3) the
  fixture deletes move construction and move assignment beside the deleted copies. Syntax-gated
  with the `test_fem` flags (`-std=gnu++17 -fsyntax-only`, exit 0) — below the compile rung.

Fixture risks, resolved by source reading (medium-high confidence; only a run confirms):

- (i) **No mesh edges are needed.** `IWG_Maxwell` sets `mUseEdges` for h-φ
  (`cl_IWG_Maxwell.cpp:58-60`), and `FieldList::initialize` pushes the conductor list into the dof
  list regardless of which blocks exist, so `"edge_h"` is in `all_fields()`. But every consumer
  that would touch a mesh edge is gated on the **block's** edge dof types, which are empty for an
  air block (`mFields.Air` is `"phi"` only, `:63`):
  `DofData::count_edge_dofs` returns 0 when no selected block has edge dofs
  (`cl_FEM_DofMgr_DofData.cpp:1080-1094`); `compute_max_ids` guards on `edges_exist()` (`:734`);
  `DofManager::create_fields` sizes `"edge_h"` to `edge_multiplicity · number_of_edges() = 0`
  (`cl_FEM_DofManager.cpp:185`); the boundary facet links `FacetAndMaster`
  (`cl_FEM_Element.cpp:287-293`), where `link_edge_dofs` returns early on an empty type list
  (`:1057-1060`); the block element's `link_dofs` only computes edge directions under
  `has_edges()` (`:172-175`). In the calculator, `"edge_h"` makes `allocate_memory` build master
  edge functions for the HEX8 master (`cl_FEM_Calculator.cpp:1036-1054`), and `link` calls
  `EF_HEX8::link( master )` (`:1827`), which reads node coordinates and the master's
  edge-direction **bitset** (`cl_EF_HEX8.cpp:178-195`, `cl_FEM_Element.hpp:634-647`) — no mesh
  edge. The canonical-edge fallback of the ghost fixture was therefore **not** used.
- (ii) The "Calculator has no plain QUAD4" refusal concerns 2-D QUAD4 **blocks**; here QUAD4 is a
  facet under a HEX8 master, and the sideset calculator dispatches on the master
  (`:1536-1590`; nedelec `HEX8` arm `:1115-1124`).
- (iii) A single-element mesh through `Mesh::finalize()`/`MeshChecker` is the proven
  `TS_TestPrism` route (`:342-381`); the added boundary facet without slave is the ghost
  fixture's `set_master` pattern minus `set_slave`. An air block without a material: no
  `material` reference in `cl_FEM_DofMgr_BlockData.cpp` or `cl_FEM_Block.cpp`.

Why not cheaper: a calculator-only pin (`Σ w dS(k)` on a tapered face) proves the premise but is
green on the current kernel; hacking the DR-178 ghost fixture by moving a node makes a trapezoid
side face of a PENTA6TS, which is the one case Christian ruled unreachable in production, on a
fixture that already carries two documented hacks. The HEX8 boundary facet is the production
shape of the bug (a background sideset on a hex air mesh) and needs no support-header change.

### R3 — close out (after: R1, R2, gate run by Christian)

- [ ] Register row DR-181: status cell → FIXED with SHA, evidence = the two gtest names and both
  trees; note the `nonfree/` sweep owed (O3); if O2 says yes, file the §3.3 row.
- [ ] Devlog `devlog/dl20260922_maxwell_background_ds_dr181.md` + `devlog/README.md` line.
- [ ] Move this file to `todo/closed/` and update `todo/README.md`.

---

## 5. Open Design Questions (not silently decided)

- **O1** — **Ruled (Christian): R2 as written, the kernel-level HEX8 fixture**, with the plan
  jury's corrections. Draft text: Default: the kernel-level HEX8 fixture of R2 (red on the current code).
  Alternative if Christian wants the round cheaper: the calculator-only pin plus
  `make check-fast` and `disk_ramp` as a no-change smoke — honest but not red.
  **Recommended: R2 as written.**
- **O2** — **Premise refuted on the controller path (§3.3); DR-183 disposition is Christian's
  ruling.** Draft text: §3.3: should the `dirichlet → BackgroundField` stamping (`cl_MaxwellFactory.cpp:908-912`)
  be filed as its own register row (ruling, P3), or is it intended and the observation dropped?
  Not a widening of this plan either way.
- **O3** — **Still owed at R3.** `nonfree/` is absent here; the row closes with the sweep owed there, or Christian
  confirms no `dS()`/`dV()` loop exists in that tree.
- **O4** — Not ruled; left as residue per the recommendation. Pin the 2-D curved-edge family (`normal_tri_curved`) separately? The same two tokens
  cover it; a TRI6 fixture needs midside offsets and the curved flag
  (`cl_EF_TestVolume.hpp:27-45` has the machinery). **Recommended: no** — one red fixture on the
  per-point family is the gate; note as residue.

---

## 7. Definition-of-Done Checklist

- [ ] Gap-table rows 1-3 land in R1/R2; row 5 resolved by O3.
- [x] Every claim above carries a `file:line`; the `nonfree/` absence is stated as scope.
- [ ] R1 → R2 → R3 in order; the gate (Part A + Part B) runs in both trees, by Christian.
- [x] O1-O4 answered in place, not silently (O2's DR-183 ruling and O3's sweep remain open, stated).
- [ ] `test_fem --gtest_filter='BackgroundField.*'` red before R1, green after; `make check-fast`
  green; `disk_ramp` or `disk_pulse` unchanged is a claim of bit-identity (§3.2), optional to run.

## 8. Audit Trail

- Exchange thread: `tmp/ai_exchange/review_maxwell_background_ds_dr181.md` (plan round, then the
  code round on an own-files-only diff: `mt_maxwell_background.cpp`, `tests/fem/CMakeLists.txt`,
  `tests/fem/test_BackgroundField.cpp`).
- Codex seat: `gpt-6-astra` (reads the Calculator, a critical-framework file); Grok third voice.
- Every auditor `file:line` re-verified against source before inclusion.
- Plan round (2026-09-22): astra (full audit) and Grok (third attempt, blind) — both flagged the
  stack-allocated BC deleted by `Kernel::~Kernel` (P1, folded into R2); astra flagged the
  "fixed anyway" dismissal, the dormant `backgrounddirichlet` arm and the conditional
  bit-identity (folded into §3.1-§3.3); Grok refuted the Dirichlet-stamp assembly premise by
  source trace and supplied the absolute-area, `"h"`-origin and `:1772-1776` corrections (folded
  into §3.3 and R2). Grok and the planner both named the edge risk; resolved by source reading in
  R2 (i).
- Code round: own-files diff at `tmp/ai_exchange/review_maxwell_background_ds_dr181.diff`.
  astra requested three test-only changes (serial guard, NaN-blind reductions, deleted moves);
  applied as E1–E3 in R2, production fix unchanged.
