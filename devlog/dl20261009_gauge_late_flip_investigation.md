# Devlog 2026-10-09 — DR-199 filed: the `gauge` boundary condition fixes its dofs after the free/fixed split

**Date:** 2026-10-09
**Topic:** Christian's instruction via the Server session ("file a dr entry and start an investigation using the jury") on the hazard both code auditors of the autopin fix raised: a `gauge` condition fixes its dofs only at the first timestep. Diagnosis, an executed reproducer, a jury round, the row. No code change
**AIs involved:** Claude Fable 5.1 (reading, reproducer, row); Codex gpt-6-astra/high and Grok grok-4.7/high (blind jury, same brief); Server session (instruction, report)
**Claude Confidence:** high on the mechanism (executed in debug, release and on two ranks); medium on the rank-1 crash site (not traced) and on the release-mode order of the two failure sites (source-traced, not backtraced)
**Literature References:** N/A — a lifetime defect, not a formulation question

## Summary

A `gauge` block is a bearing with a time-dependent value: created on point ids like a bearing,
imposed per timestep through `Bearing::impose_dirichlet`. Unlike a bearing it is pre-fixed by
no factory. Its first imposition therefore arrives after `shift_fields()` has run the one-time
`DofManager::initialize()`, which synchronizes the Dirichlet flags, splits the dofs into free
and fixed, reorders them and allocates the matrices. `Dof::fix` sets a flag and a value only.
The next assembly skips the dof's row and addresses the free×fixed Dirichlet block with the
dof's free index; the fixed-value gather, sized from the frozen count, runs one entry long.

The strip reproducer (two air components, a copper column, a point element at node 11) shows
it: with `gauge { nodes : 11 ; type : constant ; amplitude : 0 ; }` setup completes with the
two automatic pins and the first assembly aborts — debug in `assemble_jacobian`
(`aColIndex out of bounds ( 21 >= 2 )` from `SpMatrix::position`), release with Armadillo's
`Mat::operator(): index out of bounds` from the fixed-value gather after a raw one-past write
into the sparse block, two ranks with a SIGSEGV on rank 1. The same deck with `bearing
{ nodes : 11 }` passes assembly. A `gauge` block without `type` crashes earlier: the source
function's member pointer is never initialized and the global update calls it. A typed gauge
with a unit is refused at deck read, because the Maxwell units table gives gauge no unit
(the contract says dimensionless); the thermal factory's gauge takes K and late-flips the same
way after `ThermalFactory` initializes the field.

Both auditors returned REAL on the late flip and on the untyped crash, found no re-split, and
agreed on the containment: refuse a first free→fixed flip after `is_initialized()` in
`Bearing::impose_dirichlet`, as the Dirichlet-sideset path already does, and refuse an absent
`type` at deck read. The feature-preserving repair — pre-fix gauge points at factory time in
both factories — needs rulings: whether a gauge sets `tHaveBearing` (and so disables the
automatic pins), whether a gauge point becomes a contact of `impose_autopins`, and the thermal
order. No shipped deck uses the key. Predates the autopin fix; 12a8eec did not change its
reach.

## Where the auditors corrected the brief

- The release abort is the gather's Armadillo check, which survives `NDEBUG`
  (`ARMA_CHECK_CONFORMANCE` is the SCLS default); the Dirichlet-block write at
  `assemble_jacobian` is a raw one-past write in release. I had the two sites swapped.
- The key is documented: `doc/input_file_reference.md` says `type` is effectively required
  and gives the Maxwell gauge no dimension and the thermal gauge K. I had written that its
  semantics appear nowhere.
- "Any unit is rejected" overshoots: a dimensionless amplitude passes.
- MPI: Codex read worker bearings as always empty; Grok read `link_bearings_with_dofs` as
  inserting dofs on every rank for a local node. The source and the two-rank crash on rank 1
  support Grok. Recorded as an adjudicated split, not a vote.
- `impose_autopins` names the gauge in its comment and excludes it on purpose; "handled nowhere
  in the factory" was too strong.

## Disagreement left for Christian

Whether this blocks 0.9.3. Codex: block on the containment (medium confidence). Grok: a product
decision; the technical fact is that a documented key is unsafe. Claude: not a blocker for the
shipped decks (loud abort, no deck uses the key), containment recommended before the release
because the key sits in the documented table and the untyped form is a crash on a legal-looking
deck.

## Changes Made

`todo/debt_register.md`: DR-199 filed, `[CODE] [P]`, P2; the preamble's `[P]` count 2 → 3.
This devlog and its index line. Nothing in `src/`, nothing in `doc/`.

## Files Updated

`todo/debt_register.md`, `devlog/dl20261009_gauge_late_flip_investigation.md`,
`devlog/README.md`. Exchange: `tmp/ai_exchange/gauge_late_flip_claim.md`,
`tmp/ai_exchange/review_gauge_late_flip.md`. Reproducer decks in the session scratchpad only.
