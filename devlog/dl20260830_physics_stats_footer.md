# The physics summary was always one save behind, and the asterisk was on the wrong line

**Date:** 2026-08-30
**Purpose:** record why `print_physics_stats()` could not be correct where it was called, the move
into `print_footer()`, and the box-framing handoff that move needed
**Module:** `src/fem/kernel` ( `cl_FEM_Controller` )

## The question

Christian added a physics summary row to the per-timestep box — applied current, J/Jc maximum,
maximum temperature — and marked J/Jc with a trailing asterisk when the value could not be
recomputed, because it is only refreshed by a postprocessor pass. The question was whether the
asterisk was placed correctly, or whether the display lagged a postprocess behind.

It lagged. Worse, the lag fell precisely on the rows that claimed to be fresh.

## Why the predicate could not work there

`JJCx/y/z` are written by `MaxwellPostprocessor` and by nothing else, so the field is only current
immediately after `DofManager::postprocess()`. That call lives in `Controller::finalize()`. The
summary was printed from the certified-exit branch of `iterate_coupled()` / `iterate_magnetic()`,
which sets `mTripExit` and returns — before the caller ever enters `finalize()`.

So at print time the postprocess for the step being announced had not run. `save()` was a correct
predictor of *whether this step would be postprocessed* — `mTime` is already the new time by then,
advanced in `initialize_timestep()` — but a prediction about the future is not a statement about
the field's current contents. On a saved step the code took the "fresh" branch and rescanned the
field left behind by the *previous* saved step, printing a stale number with no asterisk. On an
unsaved step it printed `mLastJJcMax`, which is the same stale number, with one.

Both rows carried the same value. The only difference was that the unmarked one lied about it.
No predicate evaluated at that point in the loop could have been right: between two saves the
field is stale unconditionally.

## The fix: record the fact, not a proxy for it

`finalize()` now stores what actually happened rather than what the save grid implies:

```cpp
mPostProcessed = mKernel->dofmgr()->postprocessors().size() > 0 && aPostProcess ;
if ( mPostProcessed ) { ... postprocess() ... }
```

and the summary is printed from `print_footer()`, which runs at the end of `finalize()` — after
that pass. The un-asterisked number is now live in the strict sense, and the asterisk means the
value was measured at the last save.

A third case appeared once the branch was written honestly: when the deck has no postprocessor at
all, `JJCz` never exists. That cell is now left **empty** rather than printing `0.000000`, which
would read as a measurement of zero rather than an absence.

## A second defect the move fixed for free

In the segregated path the summary sat at the *magnetic* convergence exit, ahead of the thermal
sub-step loop. `T_max` therefore reported the temperature from before the sub-steps advanced it.
Printed from the footer it is sampled after them.

## The framing handoff

Moving the row across a function boundary split a box between two emitters. The row and the
"timestep succeeded" row above it are both three-column ( 20/25/24 ), so the divider between them
must be `┼`; a `┴` followed by a `┬` renders as a doubled rule. But in the segregated case the
thermal sub-step boxes come in between and close the section themselves, so there the summary must
open a fresh one with `┬`.

Rather than have `print_footer()` re-derive the deck topology, the emitters state what they did.
`mBoxSectionOpen` is set by whichever exit left a three-column section open:

- `iterate_coupled()` prints its row and leaves the section open
- `iterate_magnetic()` leaves it open only when `mKernel2 == nullptr`; with a thermal kernel it
  closes as before
- `iterate_thermal()` closes its own sub-step box and clears the flag

`print_physics_stats()` picks its opening divider from the flag and clears it. Nothing prints
between the trip exit and `print_footer()` — `adjust_timestep()` and `compute_conditioning()` are
both silent — so the handoff cannot be interrupted.

## Same defect shape, one row down

The "Time for postprocessing" row was gated on `mLastSave`, a shadow of the sticky `mSave` flag,
while `mPostprocesingTime` is only written where the pass ran. That is the same
proxy-instead-of-fact pattern, and it could print a stale duration. It reads `mPostProcessed` now.

`mLastSave` is consequently unread: still assigned in three places in `adjust_timestep()`, read
nowhere. Left in place pending a decision.

## Gates

Syntax-only compile of the TU with the tree's own flags from `belfem_kernel.dir/flags.make`,
including `-Wall -Werror -pedantic-errors` — clean. Frame widths measured at 76 display columns for
every emitted line, counting the combining asterisk ( U+20F0 ) as zero width. Christian then built
and ran it and confirmed the box renders correctly.

`make check` was not run. The segregated three-emitter path has not been seen in real output.
