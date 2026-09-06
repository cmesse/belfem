# Periodic Input Extension Plan

**Date:** 2026-05-15
**Purpose:** Add explicit input-file syntax for declaring periodic sideset pairs while preserving the existing plane-based periodic definition.

## Current State

The currently implemented syntax is:

```txt
topology
{
    periodic
    {
        source : 1, 2, 3 ;
        target : 5, 6, 7 ;
    }
}
```

`source` and `target` are three geometry/node IDs defining the master and slave planes. `PeriodicityFactory` then selects all mesh sidesets lying on those planes and matches nodes by projected in-plane coordinates.

This works for `cmake-build-debug/input.conf`, where sideset `7` is on the `z = 0` plane and sideset `8` is on the `z = 17.39482682` plane. The limitation is that the input file cannot explicitly say that only sideset `7` should be paired with sideset `8`.

## Recommended User-Facing Syntax

Use the existing `source`/`target` plane points for orientation, and add optional explicit sideset filters:

```txt
topology
{
    periodic
    {
        source : 1, 2, 3 ;
        target : 5, 6, 7 ;

        source sidesets : 7 ;
        target sidesets : 8 ;
    }
}
```

This is the safest first extension because it avoids inferring the in-plane coordinate system from arbitrary sideset node ordering. It also supports multiple sidesets on each side:

```txt
source sidesets : 7, 9 ;
target sidesets : 8, 10 ;
```

The counts must match unless the implementation deliberately treats the IDs as unordered sets on two planes. For the first implementation, require matching counts to catch input mistakes early.

## Optional Later Convenience Syntax

After the explicit-filter path is tested, consider a compact pair syntax:

```txt
periodic
{
    sidesets : 7 @ 8 ;
    source : 1, 2, 3 ;
    target : 5, 6, 7 ;
}
```

This would require a small `input::Section` helper for parsing one or more ID pairs from a key value. Do not make this the first implementation unless we want to touch the input parser now.

## Implementation Plan

### 1. Keep Existing Syntax Fully Backward Compatible

`MaxwellFactory::create_periodic()` must continue to accept the current `source` and `target` keys exactly as implemented today.

Validation:

- `source` exists.
- `target` exists.
- Each contains exactly three IDs.

### 2. Extend `MaxwellFactory::create_periodic()`

Add optional parsing for:

- `source sideset` or `source sidesets`
- `target sideset` or `target sidesets`

Use existing `input::Section::get_ids()`; this does not require an input parser change for the recommended first syntax.

Validation:

- If either source-side or target-side sidesets are provided, both must be provided.
- Source and target sideset counts must match for the first implementation.
- Every listed sideset must exist on `mMesh`.
- Listed sidesets must be geometrically compatible with the corresponding plane. This check can be delegated to `PeriodicityFactory` by verifying that all nodes lie on the plane within `BELFEM_MESH_EPSILON`.

### 3. Extend `mesh::PeriodicityFactory`

Add optional explicit sideset selection state:

```cpp
Vector< id_t > mMasterSideSetIDs ;
Vector< id_t > mSlaveSideSetIDs ;
```

Add public setters:

```cpp
void set_master_side_sets( const Vector< id_t > & aSideSetIDs );
void set_slave_side_sets( const Vector< id_t > & aSideSetIDs );
```

Change `create_tree()` / `select_target_nodes()` or their helper path so that:

- Without explicit sidesets, current plane auto-detection is unchanged.
- With explicit sidesets, only the listed sidesets are used.
- Listed sidesets are still checked against the plane equation.

### 4. Keep Matching Semantics Unchanged

Do not change the node matching algorithm in the first pass:

- Master nodes still build the projected k-d tree.
- Slave nodes still project through the slave transform.
- Edge and face pairs are still derived through `Periodicity::update()`.

This keeps the change isolated to input selection, not periodic matching.

### 5. Consider Whether `input::Section` Needs a New Helper

For the recommended syntax, no new helper is required.

If we choose the compact `sidesets : 7 @ 8 ;` syntax later, add a helper such as:

```cpp
void get_id_pairs( const string & aKey,
                   Vector< id_t > & aLeftIDs,
                   Vector< id_t > & aRightIDs ) const;
```

It should parse comma-separated `A @ B` pairs and reject ranges unless we define clear pair semantics.

### 6. Test Plan

Add or run a serial periodic regression using `cmake-build-debug/input.conf` with:

```txt
source : 1, 2, 3 ;
target : 5, 6, 7 ;
source sidesets : 7 ;
target sidesets : 8 ;
```

Checks:

- `mMesh->has_periodicity()` is true before cut creation.
- `Periodicity::master_side_sets()` contains sideset `7`.
- `Periodicity::slave_side_sets()` contains sideset `8`.
- Original periodic nodes, derived periodic edges/faces, and post-thin-shell duplicate periodic nodes are present.
- DOF-level periodic hanging constraints are created after the Step 2 `DofData` work.

### 7. Documentation

Update:

- `src/mesh/doc/periodicity.md`
- `todo/periodic_bc_fix_plan.md`
- Any Maxwell usage guide section that documents `topology`.

## Open Design Decision

The first implementation should **not** infer planes from sideset IDs alone. That inference is attractive for input brevity, but it risks ambiguous in-plane orientation and wrong node pairing on symmetric or curved meshes. Keep plane points mandatory until a robust orientation rule is designed and tested.
