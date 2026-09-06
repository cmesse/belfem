# Thin-shell geometric construction and periodicity propagation {#mesh_thin_shell_geometry_and_periodicity}

**Date:** 2026-06-16
**Purpose:** Describe how `ThinShellFactory` builds the thin-shell virtual domain geometrically (protoshell → surface → extruded layers) and how periodicity links are propagated from the original seam nodes to the extruded layer nodes — including the key subtlety that the extruded periodic nodes need **not** lie on the periodic plane.
**Module:** src/mesh (geometry/periodicity); implemented in `src/fem/kernel/cl_ThinShellFactory.cpp` (namespace `mesh`)
**Status:** Verified against the source by Codex 2026-06-16 (corrections to the P/terminal, normal-computation, edge-propagation, and wrapper-edge claims folded in).

---

## 1. Protoshell

A **protoshell** describes a thin shell before it is given thickness. It carries:
- `mSideSetIDs` — the sidesets whose facets form the shell **surface**, and
- `mThicknesses` — the per-layer thicknesses (from `input.conf`).

`ThinShellFactory::create( Protoshell* )` (`cl_ThinShellFactory.cpp:103`) turns a protoshell into the extruded virtual domain.

## 2. The surface node set S

Let **S** be the set of all nodes on the protoshell surfaces. It is obtained from the facets on the protoshell's sidesets:
`collect_sidesets()` (`:435`) → `collect_facets()` (`:458`) → `collect_nodes()` (`:512`).

## 3. Periodic-plane intersection P and the terminals

Conceptually, for a periodic plane **A**, the relevant nodes are the intersection **P = A ∩ S** — the surface nodes that lie on the periodic plane. **Operationally** the code does not intersect geometrically: it identifies **P** as the subset of collected surface nodes **S** that already carry periodic partners (`is_periodic()`), selected in `create_nodes_on_layers()` (`cl_ThinShellFactory.cpp:1125-1136`). Those partners come from the source/target periodic-plane pairing done earlier on the original mesh.

The **terminal curves** are *separate* Protoshell geometry (`Protoshell::terminal_curves()`), attached from mesh curves touching the thin-shell sidesets (`cl_MaxwellFactory.cpp:2582-2609`). In a periodic-seam case the terminals may *contain* nodes of P, but P is **not defined from** the terminal curves.

*Nuance (cohomology):* for cohomology generation, the thin-shell side nodes are **duplicated** and the slave-side elements relinked (`cl_CutFactory.cpp:1600-1808`), creating a **degenerate, zero-width hole** — the **slit** the protoshell sits in; terminal loops are closed with the duplicate segments (`:2049-2189`). The P that `ThinShellFactory` uses for the periodic links holds only the **original / master-side** nodes (`reset_node_indices()` asserts master nodes are non-duplicates, `cl_ThinShellFactory.cpp:396-400`; `collect_nodes()` returns `mMasterNodes`, `:537-540`).

## 4. Extrusion: S → L

To build the virtual volume, the surface S is **extruded**:
1. Compute the surface **normals** at the nodes via the element-specific `process_nodes_*()` routines selected in `create()` (`cl_ThinShellFactory.cpp:126-146`; definitions at `:555` line2, `:648` line3, `:720` tri3, `:827` tri6). *(`compute_binomial_vectors()` at `:1966` is a side-curve helper and is **not** part of the extrusion path.)*
2. Step **layer-wise along the normal**, using the thicknesses from the protoshell (`compute_distances()`, `:984`), creating the layered thin-shell node set **L** (`create_nodes_on_layers()`, `:1022`). The per-node step is literally `tX += tD * normal` (`:1081-1097`).

To navigate the node sets cleanly, the node indices of S are **temporarily overwritten to be contiguous** (`reset_node_indices()`, `:393`), and the **same** local indices are reused in each layer-wise subset of L.

## 5. Periodicity propagation: P ⊂ S ⇒ Q ⊂ L

Because **P ⊆ S**, the extruded set L contains a corresponding subset **Q** (the extrusions of the P nodes, across all layers). Each layer node `tC ∈ Q` is paired with `tD = tNodes(tA->periodic()->index())` and registered in the periodic backup (`create_nodes_on_layers()`, `cl_ThinShellFactory.cpp:1125-1161`) — so **Q inherits its node periodicity explicitly from P**. Their **edges**, however, are **not** copied directly: periodic edge pairs are re-derived later by `PeriodicityFactory::collect_edges()` / `match_edges()` from the periodic wrapper facets — which is why a conductor wrapper's facet edge containers must be linked before that pass runs (see §6).

**Critical subtlety — Q need not lie on the periodic plane A.** The extrusion steps along the shell **normal** `n`, which is **not necessarily in-plane to A**. So an extruded layer node `q ∈ Q` is periodically linked (it inherits its partner from the corresponding `p ∈ P`) but can sit **off** the plane A.

This is why the periodic seam facets of the thin-shell layers **cannot** be recovered by a purely geometric "all nodes on the plane" selection (`PeriodicityFactory::select_sidesets()`, `cl_Mesh_PeriodicityFactory.cpp:1433`): the layer nodes Q are off-plane by construction. They are periodic **by inherited link**, not by geometric position. Accordingly, the thin-shell layer periodic facets are selected by **inherited flags** (flags 4/5): `flag_periodic_nodes()` (`cl_ThinShellFactory.cpp:2941`) and `flag_layer_nodes()` (`:2967`) mark the P/Q nodes, and `create_periodic_sideset()` (`:2990`) selects facets whose nodes are all flagged — geometry is not consulted.

### 5.1 Assumption — the shell normal never lies in a periodic (or boundary) plane

The extrusion and seam logic **assume** that the thin-shell normal `n` is **never coplanar with a periodic plane A** (nor with any boundary plane). Two consequences follow, both relied upon by the periodic rebuild:

1. **No PENTA face is ever a periodic seam facet.** A PENTA lateral face is QUAD4 and is spanned by `n` (it is the swept image of a surface edge along `n`). If `n` lay *in* A, such a lateral face could become coplanar with A and be selected as a seam facet. Because `n ∉ A`, this never happens. Therefore **every 3D periodic seam facet is a TRI3 backed by a TET4 volume element** (the surrounding bulk/air mesh), never a PENTA face. This is the invariant that makes the edge-sourcing in §6 safe: `TET4::get_edges_of_facet` is defined for all four faces, whereas `PENTA6TS::get_edges_of_facet` only covers the triangular caps (3,4) — the QUAD4 laterals (0–2) deliberately have no through-thickness edges in the reduced TS edge model.

2. **The off-plane layer nodes Q never collapse onto A.** Since extrusion moves strictly along `n ∉ A`, an extruded node `q ∈ Q` is pushed off A by a non-degenerate amount; it cannot accidentally re-land on A and be mis-selected by geometric plane tests.

This is a *modeling* assumption, not a checked precondition. An optional defensive guard could verify that no protoshell-surface facet normal is in-plane with any periodic/boundary plane and error out otherwise; this is **low priority** and not currently implemented.

## 6. Consequence for edges (ties to the periodic rebuild)

The generated periodic **wrapper facets** for the layer seam are created in `create_periodic_sideset()` (`cl_ThinShellFactory.cpp:2990`) as hidden `DomainType::Periodic` facets. It calls only `set_master()` / `set_sideset_id()` — it does **not** allocate or link the wrapper's edge container. So a thin-shell periodic wrapper has **no edges regardless of layer type** (`Facet::edge()` delegates to the wrapper element, `cl_Facet.hpp:405-416`; `has_edges()` is true only after `allocate_edge_container()`, `cl_ElementTemplate.hpp:397-412`).

**This does not cost the conductor layers their periodic edge pairs.** For a **conductor** layer these facets must contribute periodic **edge** pairs (the layer carries H/edge DOFs); for ferro/void/non-conducting-buffer layers they carry none (see `src/fem/maxwell/doc/thin_shell_virtual_domains.md`). `PeriodicityFactory::collect_edges()` (`cl_Mesh_PeriodicityFactory.cpp:600`) gets the conductor case right **without** the wrapper owning edges: it tries the wrapper element first and, when that has no edge container, **reconstructs the edge set from the master volume element** via `master()->get_edges_of_facet( index_on_master() )` (`:624-635`). A bare `has_edges()` skip in its place would silently drop exactly these constraints.

The reconstruction leans directly on the §5.1 invariant. Every 3D periodic seam facet is a TRI3 backed by a TET4, and `TET4::get_edges_of_facet` is defined for all four faces — whereas `PENTA6TS::get_edges_of_facet` covers only the triangular caps (3,4), the QUAD4 laterals (0–2) having no through-thickness edges in the reduced TS edge model. A PENTA lateral turned seam facet would leave the fallback with nothing to source, which is the concrete cost of violating §5.1. The differing edge enumeration and intrinsic direction on the two sides of the seam are harmless: `match_edges()` (`:886`) pairs by endpoint originals and aligns direction by node swap.

---

## See Also

- [periodicity.md](periodicity.md) — periodic master/slave pairing and rebuild
- `src/fem/maxwell/doc/thin_shell_virtual_domains.md` — virtual domains and the edges-iff-conductor rule
- `src/fem/maxwell/doc/thin_shell_facet_orientation.md` — facet master/slave orientation
