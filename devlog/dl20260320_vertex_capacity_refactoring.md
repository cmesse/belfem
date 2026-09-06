# Vertex Capacity Refactoring and Latent Bug Fixes

**Date:** 2026-03-20
**Purpose:** Devlog for mesh::Vertex counter/capacity split and bugs surfaced by the fix
**Module:** mesh, fem/interpolation

## Summary

Refactored `mesh::Vertex` to use separate `uint8_t` counter and capacity fields,
fixing a longstanding inconsistency where a single `uint` counter served
double duty as both allocation size and fill position. Also changed
`graph::Vertex::mVertexCounter` from `uint` to `uint16_t`.

## Problem

The original design used one counter per entity type (nodes, edges, faces,
facets, elements) for three incompatible roles:

1. **Counting phase:** accumulate via `increment_X_counter()`
2. **Fill position:** reset to 0 by `allocate_X_container()`, then used by `add_X()` as write index
3. **Bounds/size:** used by `insert_X()` for bounds checking and by `number_of_X()` to report size

Roles 2 and 3 conflicted: `insert_X()` needed the counter at the allocated size,
but `add_X()` needed it at 0. As a result, `insert_X()` always failed its own
assert in debug builds after allocation, and `number_of_X()` returned 0 after
positional fills.

## Solution

Split each counter into two `uint8_t` fields (counter + capacity), fitting in
the same 2 bytes (vs. the old 4-byte `uint`). Net memory savings: 10 bytes per
mesh entity.

- `mXCapacity` — set during allocation, never reset until container is freed.
  Used for bounds checking in accessors and `insert_X()`.
- `mXCounter` — reset to 0 on allocation. `add_X()` increments it sequentially.
  `insert_X()` increments it only when writing to a previously-nullptr slot
  (requires `std::fill` to nullptr on allocation).
- `number_of_X()` returns counter (semantic count of populated entries).
- Accessor asserts check against capacity (prevents buffer overrun).

Overflow asserts added: counter < 255 for `uint8_t` mesh counters,
counter < 65535 for `uint16_t` graph vertex counter.

### Files Modified

- `src/math/graph/cl_Graph_Vertex.hpp` — `mVertexCounter` → `uint16_t`, overflow assert
- `src/mesh/cl_Vertex.hpp` — counter/capacity split, all accessors check capacity
- `src/mesh/cl_Vertex.cpp` — allocate/reset/add/insert rewritten for dual-field design

## Latent Bugs Surfaced

### 1. Missing unflag in `connect_faces_to_faces()` (cl_Mesh_ConnectivityCalculator.cpp)

With the old code, `Face::number_of_edges()` always returned 0 (because
`insert_edge()` never touched the counter). The edge-based neighbor loop in
`connect_faces_to_faces()` never executed, hiding a bug: collected faces were
not unflagged, causing double-collection and Cell out-of-bounds when a face
was reachable via two edges.

**Fix:** Added `tOther->unflag()` in the edge collection loop (line ~953).

### 2. Missing `populate_for_slave_penta()` (cl_IF_IntegrationData.cpp)

The `populate_for_slave()` switch handled TRI, QUAD, TET, HEX but not PENTA.
Never triggered before because face-edge connectivity wasn't working for the
code path that reaches PENTA sidesets.

**Fix:** Implemented `intpoints_penta()` slave overload and
`populate_for_slave_penta()` based on orientation tables derived in
`tmp/penta/orientation.m`. PENTA has mixed facet types:

- Faces 0–2: QUAD (4 orientations each, 12 cases) — maps [-1,1]² → PENTA parent via `(1±x)/2`
- Faces 3–4: TRI (3 orientations each, 6 cases) — permutations of barycentric coords

### Files Modified

- `src/mesh/cl_Mesh_ConnectivityCalculator.cpp` — unflag fix
- `src/fem/interpolation/fn_IF_initialize_integration_points_on_facet.hpp` — slave penta declaration
- `src/fem/interpolation/fn_IF_initialize_integration_points_on_facet.cpp` — slave penta implementation
- `src/fem/interpolation/cl_IF_IntegrationData.hpp` — `populate_for_slave_penta()` declaration
- `src/fem/interpolation/cl_IF_IntegrationData.cpp` — PENTA case + `populate_for_slave_penta()`

## Next Steps

- **Write proper test functions** for the PENTA orientation point mappings before
  continuing integration. The 18 orientation cases (12 quad + 6 tri) need
  verification against the Octave reference (`tmp/penta/orientation.m`).
- Run full test suite once the build passes with the new Vertex layout.
- Watch for further latent bugs surfaced by now-correct face/edge connectivity.

## Status

**In progress.** Core refactoring done. PENTA slave integration implemented
but untested — needs orientation tests before further work.
