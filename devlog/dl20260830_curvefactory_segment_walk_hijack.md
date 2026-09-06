# CurveFactory segment-walk hijack: side curves scrambled when the start segment is not cell slot 0

**Date:** 2026-08-30
**Module:** mesh (CurveFactory), consumer fem/kernel (ThinShellFactory)
**Exchange:** `tmp/ai_exchange/curvefactory_next_segment.md` (5 auditor entries, swept after distillation)

## Symptom

`examples/tape_quench_usermat`-derived deck with `edge coating : on` (3D, TRI3 shell)
aborts in `ThinShellFactory::compute_side_edge_indices` with a bare
`Key not found in map`. A diagnostic upgrade of that abort (now a named
`BELFEM_ERROR` reporting curve id, station, both nodes with originals/indices/
coordinates, radix and map size — kept, it is an always-on setup-path guard, not a
scratch probe) showed side curve 13, station 1 pairing a node at x = +5.96 mm with
the rim endpoint at x = −6.25 mm: the curve's node list was out of order.

## Mechanism (confirmed by Codex terra/high + Grok 4.6/high, plan round)

`CurveFactory::next( ..., Segment * )` substituted
`aSegments( aSegment->index() )` for empty adjacency slots, intending "the current
segment", to be dropped by the id filter. But `sort_segments` overwrites
`Segment::index()` with its traversal counter while the cell is reordered only by
the final `sort`, so the stand-in resolves by *position* — usually `aSegments(0)`,
which carries the curve-minimum element id and therefore always wins the
smaller-id candidate tiebreak. Empty slots exist exactly at the two degree-1
endpoints of an open curve, so on the first `next()` call the walk teleports to
slot 0 unless the start segment *is* slot 0 — the coincidence class that let
`examples/sidecoating` and `examples/3D_tapestack` pass. On the failing mesh one
rim started at slot 0 and one did not: one side curve sorted, the other scrambled.

## Fix (plan+audit -> code+audit, both vendors, 2 code rounds)

`src/mesh/cl_CurveFactory.cpp/.hpp`:

1. Both candidate loops skip `gNoIndex` slots; no stand-in.
2. The no-candidate branch returns `aSegment` itself (flagged; the caller breaks
   before writing another index — no duplicate walk rank possible).
3. The `aSegment` parameter lost its `const` (sole caller holds a mutable pointer).
4. The two live adjacency fills use `gNoIndex` instead of `BELFEM_UINT_MAX`
   (divergent under `BELFEM_INT64`), routed through `set_size( r, c, T )` — Grok's
   round-1 catch: the 3-arg Matrix *constructor* fill is typed `real` on both
   backends, and a 64-bit `gNoIndex` does not fit double's 53-bit mantissa; the
   constructor-filled first draft of this fix would have been UB exactly on the
   configuration it targeted. The four caller-site `BELFEM_UINT_MAX` constructor
   inits are dead values (`set_size` refills unconditionally on both backends) and
   were left untouched.

For the failing class (open polyline) both auditors rate the fix necessary AND
sufficient. No behaviour change for closed loops, LINE3, n=1, n=2, or open curves
starting at slot 0 (shape-table traced independently by both vendors, twice).

## Not this bug, recorded for later

- `sort_end_nodes` accepts disconnected degree-1/2 unions; an incomplete walk
  leaves duplicate segment indices and a `gNoIndex` OOB write in the permute.
- Closed-loop convention: `orient_segments` leaves the wrapping segment 0→n−1, so
  `collect_nodes`' closed assert fails (ties into O6 in
  `todo/hex8tb_phase2_fem_wiring.md`).
- `next(Node*)` comment says "smaller ID", code returns the larger; dead on open
  paths.
- `create_protocurves` unflags nodes once before the per-curve loop; node sharing
  across proto-curves would poison the second walk.
- Physics, separate question: on this mesh the rim lies on coplanar air/air
  sidesets (11/12), so the coating wall grows into the plane of existing facets.

## Status

Reviewed, not verified. Owed gates: rebuild + rerun of the coating deck
(`cmake-build-debug/tape_quench_usermat_coating`), and a no-regression run of a
previously-working coating deck. `devlog/README.md` index line pending — the file
was hot in a parallel session at write time.
