# One mesh global per boundary condition, plus the operating temperature

**Date:** 2026-08-26
**Topic:** Christian's request off the back of the DR-108 verification run: the exodus was
carrying one global variable per tape.

## What was wrong

The gantry deck writes its 464 tape currents as a single `current` block with an unbracketed
`input curves : 1:464` list. An unbracketed list expands to one condition per id, and every
condition was publishing its own mesh global, so the exodus held `current_1 … current_464` —
464 entries, all carrying the identical 42.0112 A. The values were equal by construction: a
block has one value function, so the per-group ordinals never encoded a difference. corc
wrote 6 copies of one number for the same reason. Nothing was wrong with the physics; the
output was simply unusable in ParaView.

## What landed (commit 2320d938)

- `cl_MaxwellBoundaryConditionFactory.cpp`: only the first-created member of a block carries
  the label, so a block publishes exactly one global. The rest keep an empty label, which
  `update_global()` already no-ops on, and every condition still receives its function and
  direction — the boundary-condition objects are untouched, since `set_currents` and the
  cohomology need all 464 of them. The section-ordinal suffix stays, so two blocks of the
  same type still publish `current_1` and `current_2`.
- `cl_MaxwellFactory.cpp`: the operating temperature publishes as the global `temperature`
  from `gTbulk`, inside the existing rank-0 block, with a hard error if a boundary condition
  already claims that name and a NaN skip for decks that never set it.
- `doc/input_file_reference.md`: the reference documented the retired per-group rule and is
  corrected, with the temperature global, its distinctness from the nodal field `T`, and the
  warm-restart caveat now written down.

## Gates

| deck | before | after |
|---|---|---|
| gantry (1 block, 464 curves) | 464 × `current_N` | `current`, `temperature` |
| corc (1 block, 6 curves) | 6 × `current_N` | `current`, `temperature` |
| 2D_Undulator (2 userdefined blocks, 168 curves) | 168 × `current_N` | `current_1`, `current_2`, `temperature` |

The Undulator is the interesting one: its two globals hold **+0.0467 and −0.0467 A**, so
genuine per-block variation survives while the per-member duplication collapses. That is the
`userdefined` path the audit flagged as unprovable, exercised rather than argued.
`make check-fast` 9/9, `check_doc_claims.py` 34/34, and a final probe-free rerun with no
memdump confirmed two globals tracking the ramp (3.82 A at the first save, 21.18 A at the
sixth, temperature steady at 77 K).

## The one warning worth repeating

Memdump load re-adds any dumped label the current run has not created, and the ID-mismatch
check never inspects the label *set*. A dump written before this change therefore resurrects
all 464 retired names alongside the new two. **Delete `memdump.hdf5` after any change to how
globals are named.** The pre-fix dump in `cmake-build-debug/gantry` was moved aside to
`memdump_pre_dr108fix_20260826.hdf5.bak` — it also carried wrapped-constraint state from
before the DR-108 fix, so restarting from it would have continued wrong physics.

## Round notes

Plan and code both audited by Codex and Grok. Neither found a defect in the change; both sent
back stale text, all fixed in-round — the reference lead-in still said "every value-imposing
*condition*", the factory comment still described the retired ordinal, and my "pure in those
two arguments" wording was corrected to the real hazard, which is per-instance init side
effects in a user library. Grok also required the temperature global, the `T` distinction and
the restart caveat to be documented, and supplied the two naming-check decks.

The one vendor split — whether a `userdefined` block can produce per-member values — was
settled by reading the construction site: the source function is built inside the per-member
loop, so instances really are per-member (Codex right on the fact), but all are initialised
from the same `(file, label)` of the same block, so divergence needs an impure library (Grok
right on the consequence). Recorded as an accepted assumption rather than defended against.

Round record: `tmp/ai_exchange/bc_globals_one_per_condition.md`.
