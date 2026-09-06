# Thin-Shell Nedelec Elements {#fem_interpolation_nedelec_thinshell}

**Date:** 2026-04-17 (rev. 2026-06-05: side-connector material removed; rev. 2026-07-31: HEX8TB returned as the wall element, notes updated)
**Module:** `src/fem/interpolation/nedelec`
**Purpose:** Document BELFEM's reduced-dimension Nedelec elements used for thin-shell tape models: the surface elements `QUAD4TS` and `PENTA6TS`.

**Historical note:** earlier revisions of this file also documented a side-connector element (`HEX8TS`, and before it an experimental `HEX8TB`) used to model a wrap-around metallic connector on the tape edge. That wrap construction was removed from BELFEM in June 2026 after it was found to be unphysical (it could not represent the binormal-H field at the shell/connector fold and suppressed the side-curve current density). Side connectors returned in July 2026 as a **different** element: the degenerate `HEX8TB` wall element with four longitudinal edge dofs and no binormal/normal dofs — the free-parameter failure mode of the wrap is designed away. See `src/fem/maxwell/doc/side_coating_wall_element.md` for the full theory and status; references to the old `HEX8TS` wrap remain stale.

---

## 1. Overview

BELFEM models REBCO tapes as thin stacks of material layers whose thicknesses are far smaller than the other geometric dimensions. Resolving every layer with regular 3D volume elements would require extreme aspect ratios and a prohibitive number of DOFs. The thin-shell approach collapses each layer to a reduced element and performs the through-thickness behavior in the formulation rather than in the mesh.

BELFEM uses two reduced edge-element types in this context:

| Element | Base geometry | Embedding | DOFs | Use case |
|---|---|---|---|---|
| `QUAD4TS` | `LINE2` | line in 2D | 2 | 2D thin shells |
| `PENTA6TS` | `TRI3` | surface in 3D | 6 | 3D thin shells |

`QUAD4TS` and `PENTA6TS` are the classical thin-shell elements in the Alves sense: one reference direction is the collapsed through-thickness direction.

---

## 2. `QUAD4TS` and `PENTA6TS`

### 2.1 Geometry and Jacobian handling

`QUAD4TS` and `PENTA6TS` are surface elements embedded in a higher-dimensional physical space. Their geometric Jacobian is therefore rectangular, so the implementation uses the standard Gram-matrix / pseudo-inverse construction:

```text
G   = J J^T
J^+ = J^T G^-1
```

The columns of `J^+` are the physical gradients of the in-surface reference coordinates. The through-thickness direction is handled separately with an explicit shell thickness and a shell normal.

This is the pattern used in:

- [`cl_EF_QUAD4TS.cpp`](../nedelec/cl_EF_QUAD4TS.cpp)
- [`cl_EF_PENTA6TS.cpp`](../nedelec/cl_EF_PENTA6TS.cpp)

### 2.2 DOF structure

The reduced shell elements carry only face-edge DOFs:

- `QUAD4TS` has 2 edge DOFs, one on each shell face.
- `PENTA6TS` has 6 edge DOFs, three on the bottom triangular face and three on the top triangular face.

There are no through-thickness edge DOFs. The through-thickness dependence is represented by the shell formulation itself.

### 2.3 Layer interfaces inside the tape stack

Inside a multilayer shell stack, neighboring shell layers may either:

- share edge DOFs directly,
- use duplicated edges and ghost-facet coupling on H-H interfaces, or
- collapse to nodal `phi` behavior for insulating layers such as the buffer.

The mesh-side plumbing for this lives in `src/fem/kernel/cl_ThinShellFactory.cpp`, especially the layer creation, duplicate-edge creation, ghost-facet generation, and buffer handling.

### 2.4 Buffer treatment

For a `phi`-formulation buffer layer, the shell edge DOFs are not meant to represent a conducting layer. The buffer block is converted to nodal behavior and the shell-edge values are rerouted to the corresponding nodal sources in `create_buffers()`. The buffer is therefore locally insulating even if the surrounding conductor is globally connected elsewhere.

---

## 2a. `QUAD4TS` layer elements are left-handed by construction

The 2-D layer elements produced by `ThinShellFactory` have a **negative** scalar Lagrange
determinant. This is a deliberate consequence of the construction, not a mesh defect, and it must
not be "repaired" by renumbering the corners.

`process_nodes_line2` builds each `LINE2` facet normal as `n = ( t_y, -t_x )` — the facet tangent
rotated *clockwise*. Layer offsets ascend, so the top curve always lies at `+h·n`, and
`create_elements_on_blocks_line2` winds the quad `( bottom0, bottom1, top1, top0 )`. Its signed
area is `h·( t × n ) = -h·|t|`, negative for every facet at every orientation — `n` is *derived
from* `t`, so reversing the facet reverses both.

That ordering is the only legal one. Two contracts pin it:

- `Element_QUAD4TS::get_nodes_of_edge` declares edge 0 to be nodes `{0,1}` (the bottom curve) and
  edge 1 to be `{3,2}` (the top curve), while `ThinShellFactory::link_elements_with_edges` assigns
  the bottom layer's mesh `Edge` to local edge 0 *positionally*. So `{0,1}` must stay the bottom
  curve.
- The edge signs `mS` are derived from the element's node order (`Element::compute_edge_directions`
  compares the mesh `Edge`'s node ids against `get_nodes_of_edge`), whereas the Nédélec tangent
  `∇ξ` is taken from the **facet** (`EF_QUAD4TS::link`). Their product is the basis direction, so
  the local edge order must follow the facet edge direction.

Keeping both forces `( bottom0, bottom1, top1, top0 )`, which is clockwise. **No permutation of the
four corners is both counter-clockwise and sign-neutral.** The 3-D sibling has no such conflict:
`process_nodes_tri3` uses the right-hand normal `( B-A ) × ( C-A )`, and the `PENTA6TS` stacking
agrees with it, so those elements are positively oriented.

### What this means for consumers

The map is still a diffeomorphism, so only *measures* are affected, never gradients:

- `N` is reference-space, and `B = J⁻¹·dN/dξ` returns true Cartesian gradients whichever
  orientation the map has (the inverse divides by the signed determinant, so the sign cancels).
- The **volume weight** is the one orientation-sensitive quantity. Change of variables asks for
  `|det J|`, which is why `Calculator::dV_quad4ts` and `Pipette::measure_quad4ts` exist: they are
  QUAD4TS-only wrappers taking the absolute value. That is the exact weight, not a clamp.
- The magnetic solve is unaffected in practice — it reaches those wrappers too, but takes the
  edge-function branch where the value is `0.25·thickness·length`, positive for any thickness that
  reaches assembly, so the absolute value is the identity.

Two traps follow, and both are real:

- `MeshChecker`'s swap table lists `QUAD4TS`. Its `swap_quad4` exchanges nodes 1 and 3, which on a
  thin-shell element turns the *tangential* pair into the *through-thickness* pair — the edge
  function would then read the tape length as its thickness. The checker never sees these elements
  (they are created later, after edges exist), but the table should not be taken as an
  endorsement.
- A negative layer thickness in the deck inverts the stack and flips all of the above. Nothing
  rejects it at parse time; `Kernel::compute_element_volumes` catches it later.

---

## 3. Implementation map

| File | Role |
|---|---|
| [`cl_EF_QUAD4TS.cpp`](../nedelec/cl_EF_QUAD4TS.cpp) | 2D shell edge function |
| [`cl_EF_PENTA6TS.cpp`](../nedelec/cl_EF_PENTA6TS.cpp) | 3D shell edge function |
| `src/fem/kernel/cl_ThinShellFactory.cpp` | Mesh creation for thin-shell layers |

---

## 4. Summary

BELFEM uses two reduced Nedelec element families for tape-shell work:

- `QUAD4TS` for 2D thin shells,
- `PENTA6TS` for 3D thin shells.

Both are the literature-based thin-shell elements derived from the collapsed-shell construction (Alves). References to the old `HEX8TS` **wrap** side-connector construction are stale — only that construction was removed (June 2026); the `HEX8TS` element type and `EF_HEX8TS` themselves remain in the code as thin-shell machinery. `HEX8TB` is current again since July 2026 as a different element (the 4-dof side-connector wall element; see the historical note at the top and `src/fem/maxwell/doc/side_coating_wall_element.md`).

---

## 5. References

- **Alves, B. de Sousa et al. (2022b)**, "A thin-shell H-phi-formulation for superconducting devices." *Supercond. Sci. Technol.* **35**, 024001.
- **Messe, C. et al. (2023)**, "BELFEM: a special-purpose finite-element code for the magnetodynamic modeling of high-temperature superconducting tapes." *Supercond. Sci. Technol.* **36**, 114001.
- **General Nedelec references** — see [nedelec.md](nedelec.md).
