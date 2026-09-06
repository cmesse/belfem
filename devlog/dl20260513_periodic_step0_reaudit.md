# Devlog 2026-05-13 - Periodic Step 0 Re-Audit

**Date:** 2026-05-13
**Topic:** Re-audit of merged periodic cohomology/input Step 0 on `periodic_new`
**AIs involved:** Codex
**Codex Audit Confidence:** high that Step 0 code is present and builds; medium that behavior is correct before a periodic end-to-end test
**Literature References:** N/A

## Scope

Read-only source audit after Step 0 commits were applied:

- `3940461`
- `2ca7a04`
- `cec4fb7`
- `75d35a7`

## Findings

The Step 0 merge delivered the expected periodic cohomology/input pieces:

- `MaxwellFactory` reads `topology/periodic`, creates a `PeriodicityFactory`, passes it to `CutFactory`, and `CutFactory` creates mesh periodicity after edges/faces exist.
- `CutFactory` now builds `SimplicialComplex( mMesh, true )`.
- `SimplicialComplex` has the `aPeriodicity` parameter, periodic chain/cochain coupling, and `mOriginalEdges` tracking.
- `Cohomology::clean()` flags periodic slave entities, redirects coboundary cleanup through master-side entities, and guards against reduced-away edges using `original_edges()`.
- `CutData` mirrors cohomology edges and plus/minus bitsets onto periodic counterparts.
- `CutProcessor` has the null-master guard.
- The five `unique()` calls in `Input_Section::get_ids()` are commented out.
- The sideconnector-line predicate in `SimplicialComplex` is preserved as `is_flagged() && dimension() == tDim`.
- The known periodic-branch debug artifacts were not imported: no unconditional `mMesh->save("mesh.exo")`, and `write_debug_cohomology()` calls remain commented.

## Caveats

- `domain_type()` currently parses only the generic string `"periodic"` to `DomainType::Periodic`. It still does not parse `"air periodic"`, `"buffer periodic"`, `"conductor periodic"`, or `"ferro periodic"` to the richer periodic side-set domain types.
- `Domain::Domain()` does not handle `DomainType::Periodic` as a sideset/block domain. This is acceptable if `topology/periodic` is only a special section consumed by `MaxwellFactory::create_periodic()`, but it does not make periodic sideset sections usable.
- This audit did not run an end-to-end periodic input case, so behavior still needs a regression test after Step 1/Step 2 are in place.

## Verification

- `git diff --check HEAD~4..HEAD` passed.
- `make -C cmake-build-debug hphirun -j2` passed.

## Files Touched

- `devlog/dl20260513_periodic_step0_reaudit.md`
- `devlog/README.md`

