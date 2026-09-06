# Devlog 2026-05-13 - Periodic Branch Gap Audit

**Date:** 2026-05-13
**Topic:** Second-opinion audit of `periodic` versus `periodic_new` periodic cohomology changes
**AIs involved:** Codex
**Codex Audit Confidence:** high for the branch gap; medium for exact cherry-pick conflict behavior because Git index writes are sandbox-blocked
**Literature References:** N/A

## Scope

Read-only comparison of the active `periodic_new` branch against the older `periodic` branch after Claude reported that periodic cohomology fixes were never merged into the sideconnector line.

## Findings

The correction is valid. `periodic_new` contains mesh-level periodic hooks in `MaxwellFactory` and periodic duplicate propagation in `CutFactory`, but it does not contain the full periodic cohomology/input machinery from `periodic`.

Files with periodic branch changes missing or incomplete on `periodic_new`:

- `src/homology/cl_Cohomology.cpp`: periodic-aware `clean()` logic is missing from `periodic_new`.
- `src/homology/cl_SimplicialComplex.{cpp,hpp}`: `aPeriodicity`, periodic chain/cochain coupling, and `mOriginalEdges` are missing from `periodic_new`.
- `src/homology/cl_CutData.cpp`: periodic cohomology-edge handling is missing from `periodic_new`.
- `src/io/cl_Input_Section.cpp`: `unique()` calls remain active on `periodic_new`, so periodic source/target node association can be reordered.
- `src/fem/maxwell/cl_MaxwellFactory.{cpp,hpp}`: periodic input setup and `PeriodicityFactory` plumbing from `periodic` are missing from `periodic_new`.
- `src/homology/cl_CutFactory.{cpp,hpp}`: the `periodic` branch adds `set_periodicity()` and creates periodicity after edges/faces exist; `periodic_new` has separate sideconnector changes in the same files.
- `src/mesh/en_DomainType.cpp`: `periodic_new` has `DomainType::Periodic` and specific `*Periodic` enum values, but the string parser does not parse `periodic`, `air periodic`, `buffer periodic`, `conductor periodic`, or `ferro periodic`.

Two periodic-branch hunks should not be carried over blindly:

- `MaxwellFactory::create_magnetic_kernel()` contains an unconditional `mMesh->save("mesh.exo")`.
- `CutFactory` has an uncommented `write_debug_cohomology()` call.

## Merge Assessment

The four periodic commits identified by Claude are the relevant commits:

1. `45bf086` - bulk periodic cohomology/input plumbing.
2. `6c5c985` - `Cohomology::clean()` periodic refinement.
3. `b1371f8` - further `clean()` and `SimplicialComplex` refinement.
4. `be97d6b` - removes `unique()` calls from input ID parsing.

The postproc merge commits on `periodic` should not be merged wholesale into `periodic_new`.

`cl_Input_Section.cpp`, `cl_CutData.cpp`, and `cl_CutProcessor.cpp` look low-risk to port manually. `cl_SimplicialComplex.cpp` has one sideconnector-line change in the element selection predicate that must be preserved. `cl_Cohomology.cpp`, `cl_CutFactory.cpp`, `cl_MaxwellFactory.cpp`, and `en_DomainType.cpp` need manual merge decisions.

## Recommendation

Add a Step 0 to the periodic BC plan: manually port the periodic cohomology/input machinery from the four periodic commits into `periodic_new`, preserving sideconnector-line fixes and excluding debug artifacts. Re-audit the merged result before doing the DOF registration and MPI periodic-constraint fixes.

## Files Touched

- `devlog/dl20260513_periodic_bc_audit.md`
- `devlog/dl20260513_periodic_branch_gap_audit.md`
- `devlog/README.md`

