# The `.bfm` Mesh File Format {#mesh_bfm_file_format}

**Date:** 2026-07-01
**Purpose:** Definition of the BELFEM `.bfm` (HDF5) enriched-mesh file format — layout, datatypes,
semantics, and the save/load contract.
**Module:** `src/mesh` (`cl_Mesh_BfmFile.{hpp,cpp}`, `cl_ProtoMesh.{hpp,cpp}`)

---

## 1. Purpose and Design Philosophy

A `.bfm` file is an **enrichment cache**, not a general mesh exchange format. Maxwell runs spend
significant time in the cohomology/cut, thin-shell, hanging-entity, and periodicity factories; the
`.bfm` stores the *result* of that enrichment so a rerun (parameter change or restart) can load the
finished mesh and skip the factories entirely.

Three principles govern the format:

1. **The input file remains the source of truth for physics.** The `.bfm` stores geometry,
   topology, and enrichment state only. Regular block→material assignment is rebuilt from
   `input.conf` on every run; physical fields and the time cursor live in a **separate** restart
   file. (The path-taking `Mesh::save_fields` / `load_fields` overloads are commented out at
   `cl_Mesh.hpp:1075-1078`; what remains takes an open `hid_t`, so it is not a
   filename-level restart API.) The one exception: thin-shell *layer*
   materials are stored in `/thinshells`, because the input-domain pass skips `DomainType::ThinShell`.

2. **All cross-references are by entity ID, never by array index.** References resolve through the
   mesh ID maps (`Mesh::node(id)`, `element(id)`, …), the same maps the MPI `ProtoMesh`
   reconstruction uses. IDs are stable across finalize/reorder (see
   `mesh_contracts_and_invariants.md` §4), so the format is robust against reordering — the
   historic index-keyed writer was not.

3. **Writer and reader share one reconstruction engine.** `BfmFile` populates and consumes the
   `proto::*Data` structs of `ProtoMesh` — the same ID-keyed structs the MPI `Distributor` uses to
   send a mesh to worker ranks. The HDF5 file is effectively a serialized `ProtoMesh` stream.

The file is written and read **on the master proc only** (rank 0). Workers receive the mesh
through the normal MPI distribution afterwards.

### The checksum gate

`/meta/checksum` stores the checksum of the **base** (pre-enrichment, post-scaling) source mesh.
On a rerun, `MaxwellFactory::read_mesh` compares the freshly-loaded `.msh` checksum against the
sidecar. If it matches, the `.bfm` is loaded and the enrichment factories are skipped; if it does
not, the `.bfm` is regenerated. Storing the *base* value (not the enriched mesh's) is what makes the
gate meaningful across runs.

---

## 2. Encoding Conventions

| Convention | Meaning |
|---|---|
| `id_t` | `uint32` — entity IDs (nodes, elements, …). `0` and `gNoID` are invalid IDs. |
| `index_t` / `uint` | `uint32` — counts and enum values |
| `real` | `float64` |
| `uchar` | `uint8` — small enums, local indices, orientation codes, bit flags |
| `suint` | `uint16` — `EntityType` codes in `/hanging` |
| `string` | HDF5 variable-length byte string |
| `[N]` | fixed-length dataset, one row per entity |
| `[N][var]` | HDF5 **variable-length (vlen)** dataset (`hdf5::Dataset` helper, `cl_HDF5_Dataset.hpp`): one variable-length row per entity |
| *(opt)* | group/dataset written only when non-empty; the reader guards with `group_exists` / `dataset_exists` |

Coordinate matrices are stored **transposed**: in memory the canonical layout is `(dim × n)`
(column-major, one entity per column); on disk the dataset is `(n × dim)` — one entity per row,
so the file reads naturally in HDFView.

Enum values (`DomainType`, `ElementType`, `EntityType`) are stored as their numeric codes; the
authoritative definitions are `en_DomainType.hpp`, `Mesh_Enums.hpp`.

---

## 3. Group Reference

Groups appear in `BfmFile::save()` order. `(opt)` marks groups and datasets that are written only
when non-empty; the loader guards them with `group_exists` / `dataset_exists`. Counts used below:
`B` blocks, `S` sidesets, `N` nodes, `E` elements, `Fc` facets, `Ed` edges, `F` faces, `Cp`
control points, `H` hanging entities per type, `T` thin shells, `V` vertices, `C` curves.

```
/meta
    dimensions     scalar         uint       (num spatial dimensions, 2 or 3)
    entities       [7]            index_t    (entity counts indexed by EntityType;
                                              NODE/EDGE/FACE/FACET/CONTROLPOINT populated)
    groups         [4]            index_t    ([num_blocks, num_sidesets, num_curves, num_thinshells])
    checksum       scalar         uint64     (base-mesh checksum, see section 1)
    belfem         scalar         string     (optional; writing build's version)
    git            scalar         string     (optional; writing build's commit, "-dirty" suffixed)
    branch         scalar         string     (optional; writing build's branch)
    config         scalar         uint64     (optional; fingerprint of the SETTINGS the mesh was
                                              enriched with — layers, cuts, coating, topology.
                                              The checksum above identifies the base mesh only,
                                              so this is what catches an edited input deck on an
                                              unchanged .msh. Absent on files written before it
                                              existed, which the reader treats as a mismatch)
    config_text    scalar         string     (optional; the canonical text behind `config`, so a
                                              mismatch can name the setting that changed)
/blocks
    ids            [B]            id_t
    labels         [B]            string
    elements       [B]            index_t    (element count per block → positional membership, section 4)
    types          [B]            uint       (ElementType per block — single-type per block;
                                              authoritative, disambiguates e.g. TET20 vs HEX20)
    domains        [B]            uint       (DomainType)
    tags           [B]            uchar      (bit flags: 1 = hidden, 2 = has edges, 4 = has faces)
/sidesets                         (same dataset set as /blocks; domains incl. Ghost = 37)
    ids, labels, elements, types, domains, tags
/nodes
    ids            [N]            id_t
    coords         [N, dim]       real       (saved transposed → N×dim on disk, HDFView-readable)
    abstract       [*]            id_t   (opt) abstract-node IDs (cut jump carriers)
    pinned         [*]            id_t   (opt) auto-pin node IDs (Mesh::autopins())
    orphaned       [*]            id_t   (opt) orphaned-node IDs
    duplicates     [*][var]       id_t   (opt) per original: [original_id, dup_id, dup_id, ...]
/elements
    ids            [E]            id_t       (in block order, section 4)
    physical       [E]            uint   (opt) physical tags
    topology       [E][var]       id_t       (node IDs)
/facets
    ids            [Fc]           id_t       (in sideset order, section 4; every ID appears once)
    physical       [Fc]           uint   (opt)
    topology       [Fc][var]      id_t   (opt) node IDs — only if aSaveTopology (save passes false;
                                              loader re-derives facet nodes from the master element)
    elements       [Fc][var]      id_t       (master[/slave] element IDs)
    indices        [Fc][var]      uchar      (idx_on_master, idx_on_slave, orient_on_slave;
                                              row length encodes the master/slave case)
/edges                            (opt; H(curl) meshes)
    ids            [Ed]           id_t
    topology       [Ed][var]      id_t       (node IDs; element→edge incidence is reconstructed
                                              on load by node-key matching, IDs preserved)
/faces                            (opt; 3D, order ≥ 2)
    ids            [F]            id_t
    elements       [F][var]       id_t       (master[/slave] element IDs)
    indices        [F][var]       uchar      (idx_on_master[, idx_on_slave, orient_on_slave])
/control_points                   (opt; tensor/B-spline meshes)
    ids            [Cp]           id_t
    coords         [Cp, dim]      real       (saved transposed)
    topology       [E][var]       id_t       (control-point IDs per element, one row per element in block order; empty row when the element has none)
/hanging                          (opt; one subgroup per entity type with >0 hanging)
    /nodes | /edges | /faces | /facets | /control_points
        ids        [H]            id_t       (the hanging targets)
        topology   [H][var]       id_t       (source entity IDs)
        types      [H][var]       suint      (EntityType per source)
        weights    [H][var]       real       (T-matrix weights)
/periodic                         (opt; only if has_periodicity)
    planes         [2, 3]         id_t       (row0 master / row1 slave × 3 plane-node IDs)
    nodes          [2, *]         id_t   (opt) per col: (id, periodic()->id) — saved transposed
    edges          [2, *]         id_t   (opt) authoritative partner pairs
    faces          [2, *]         id_t   (opt)
    facets         [2, *]         id_t   (opt)
/thinshells                       (opt)
    ids            [T]            id_t       (sideset ID of each shell's facet sideset)
    ghost          [T]            id_t       (ghost sideset ID; 0 = none)
    layers         [T][var]       id_t       (layer block IDs, incl. buffer layers)
    thicknesses    [T][var]       real       (layer thicknesses, aligned with layers)
    materials      [sum T]        string     (layer material labels, flattened across shells;
                                              de-flattened on load via the layers row lengths)
    coatings       [T][var]       id_t   (opt) side connector block IDs, one per side curve;
                                              written only when some shell has connectors
    seams          [T][var]       id_t   (opt) recovery sideset ID per connector block,
                                              aligned with coatings by position
    widths         [T][var]       real   (opt) connector block width, aligned with coatings
                                              (the only record of it: block thickness is not
                                              part of the group data)
/vertices                         (opt)
    ids            [V]            id_t       (vertex element IDs)
    nodes          [V]            id_t       (the vertex's node ID — may differ from the vertex ID)
/curves                           (opt)
    ids            [C]            id_t
    labels         [C]            string
    sidesets       [C, 2]         id_t       (terminal sidesets a/b; 0 = none; saved transposed)
    types          [C]            uchar      (segment ElementType per curve)
    closed         [C]            uchar      (1 = closed loop)
    segments       [C][var]       id_t       (segment IDs per curve)
    lengths        [C][var]       real       (arclength table per curve)
    topology       [sum C][var]   id_t       (node IDs per segment, one row per segment,
                                              in segments order across all curves)
/circuit                        NOT WRITTEN -- see the note below
    time         scalar   real    circuit clock mTime at dump (O2)
    delta_time   scalar   real    Δt of the last completed step — self-description for prev_x;
                                  overwritten by set_timestep(input Δt) before first use (O2)
    x            [n]      real    solution vector (node voltages + unknown branch currents)
    prev_x       [n]      real
```

> **`/circuit` is not written today.** `BfmFile::save()` finishes after the curve data and
> emits no circuit group (grep for `circuit` in `cl_Mesh_BfmFile.cpp` returns nothing). The layout
> below records the intended block; a reader must not expect to find it in a current file.

Notes:

- A sideset may legitimately have **zero** facets: after the thin-shell factory, the original tape
  sidesets remain on the mesh as empty groups (their facets moved to the aggregate sideset), and
  they are still looked up by ID when the input domain types are re-applied on reload. The loader
  therefore keeps zero-count sidesets (`ProtoMesh::create_sidesets( aKeepEmpty = true )`); the MPI
  distributor path continues to drop them on workers.
- An entity is *hanging* when it carries source entities (static condensation / T-matrix
  constraints — periodic slaves, interface hanging edges, cut duplicates). The loader rebuilds the
  mesh-level source containers from `/hanging`; this is also how periodic DOF constraints are
  restored, since the save runs after `Periodicity::set_entity_dependencies()`.
- `/periodic` pairs are stored by authoritative partner (`entity->periodic()->id()`), not by list
  position, so reload restores true partner links without geometric re-matching.
- The `.bfm` stores the *realized* thin-shell mesh — layer blocks, facets, and ghost sidesets are
  regular entities in the groups above; `/thinshells` only records the `ThinShell` objects that
  tie them together, plus the layer materials.

---

## 4. Positional Membership

Two membership relations are positional rather than ID-keyed, mirroring how the MPI distributor
streams entities:

- **element → block:** `/elements` rows are consumed in `/blocks` order, `/blocks/elements` rows
  at a time;
- **facet → sideset:** `/facets` rows are consumed in `/sidesets` order, `/sidesets/elements`
  rows at a time.

Everything else — facet master/slave, edge/face incidence, duplicates, hanging sources, periodic
pairs, thin-shell blocks, curve segments — is by entity ID.

Container order in the file is the save-time container order of the rank-0 mesh, and the loader
preserves it. This makes the reload deterministic (same DOF enumeration as the run that wrote the
file); it does *not* rely on that order for correctness, since all references are ID-keyed.

---

## 5. The Reload Contract

What the reader (`BfmFile::load`) does:

1. Rebuilds all entities via the `ProtoMesh` builders in save order (nodes → elements → facets →
   edges → faces → control points → duplicates → hanging → periodicity → thin shells → vertices →
   curves).
2. Restores the base checksum (`Mesh::force_checksum`), so the sidecar gate works on the next run.
3. Performs a single `Mesh::finalize()`, which recomputes everything deliberately **not** stored:
   ID maps, node↔element/facet/node connectivity, edge directions (from `original()->id()`
   ordering), facet orientations, sideset/block index stamps.

What is deliberately **not** in the file:

| not stored | why / how it comes back |
|---|---|
| connectivity (node→element, node→facet, node→node, element→element) | recomputed by `Mesh::finalize()` / `ConnectivityCalculator` on load, and at `Kernel::distribute_mesh` time |
| edge directions, curved-element flags | recomputed (deterministic from node IDs / geometry) |
| element→edge/face incidence | reconstructed by node-key matching from `/edges`, `/faces` |
| regular block→material | input file is authoritative (rebuilt every run) |
| physical fields, global variables, time cursor | separate restart file (`Mesh::save_fields` / `load_fields`), checksum-guarded to the mesh identity |
| entity owners (`proc_t`) | assigned by the partitioner at `Kernel` construction on every run |
| the Maxwell topology type map | **synthesized, not stored**: the fresh path builds it *before* enrichment, so a reload must not scan the enriched mesh with `Topology::run()`. `Topology::run_on_enriched_mesh()` reproduces the pre-enrichment view by excluding everything reachable from `mMesh->thin_shells()` and skipping sideset-type re-detection |

Consumers must respect one invariant when extending the format: **anything the fresh path derives
before enrichment must either be stored or synthesized with the pre-enrichment view in mind** —
scanning the enriched mesh would give a different answer.

---

## 6. Versioning and Open Items

- `/meta/format_version` is not yet written. Add it (and a reader gate) before the next
  intentional format break. Deliberate breaks so far (all pre-1.0, no shim): the thin-shell
  `blockids` dataset was renamed `layers` when the side-connector record was added, so an
  older `.bfm` fails at thin-shell load with "Dataset layers does not exist" — delete the
  cached file and it regenerates from the source mesh.
- Global variables (total current, energy, …) live in the fields restart file
  (`Mesh::save_globals` / `load_globals`), not in the `.bfm`.