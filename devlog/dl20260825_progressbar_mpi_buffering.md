# Progressbar under MPI: the missing flush

**Date:** 2026-08-25
**Purpose:** Record the root cause and fix for the progress bar arriving in one burst under `mpirun`
**Module:** `src/core` (`Progressbar`); symptom reported on the `Cohomology` bar

## Symptom

Christian: under `mpirun` the progress bar (example: `Cohomology::clean_spfa`) shows nothing at all
for the whole run, then the finished 100 % bar appears at once. In a plain terminal the same bar
moves, if jerkily.

## Root cause

`Progressbar` redraws in place with `\r` and never emits a newline, and nothing in the class
flushed the stream. Whether the frames reach the user is therefore decided entirely by the C
library's buffer policy for `stdout`:

| stdout is | glibc policy | buffer | effect on a bar that never writes `\n` |
|---|---|---|---|
| a terminal | line buffered | 1024 B | flushes only when the buffer *fills*, ≈ every 13 frames — the jerky-but-moving case |
| a pipe (this is `mpirun`) | fully buffered | 4096 B | a short bar never fills it; everything lands at exit |

Measured with a standalone probe (40 steps, one frame per step) linking the real `Progressbar`:

```
OLD, stdout on a pty : 5 writes  (1024 B each)     <- coarse, but visible
OLD, stdout on a pipe: 1 write   (3296 B, t=0.79s) <- the entire bar at exit
NEW, stdout on a pipe: 40 writes (76 B each, one per frame)
```

3296 B < 4096 B is the whole story: the bar is simply too small to fill a pipe buffer, so it is
never flushed until `exit()`. Nothing about MPI itself is involved — the rank-0-only guard was
already correct (`CutFactory::run()` gates the cohomology path on `mCommRank == mMesh->master()`),
and there was never a second writer. `mpirun` only changes what `stdout` *is*.

A second defect made the bar far more expensive than intended: `mStep` was never assigned, so the
`tStep > mStep` throttle always compared against 0 and every `step()` call redrew the whole bar.
That is invisible while nothing flushes; with a flush per redraw it would have been one `write()`
per element in the Exodus writer.

## Fix — `src/core/cl_Progressbar.{hpp,cpp}`

- **`flush()`** after every frame, and after `reset()` / `finish()`. Because a frame is built with
  buffered `fprintf` calls and flushed once at the end, the complete ~76-byte frame leaves in a
  single `write()` — the MPI I/O forwarder cannot tear it apart either.
- **`draw( aStep, aProgress )`** extracted from `step()`; it is the one place that assigns
  `mStep = aStep`, which repairs the throttle. Redraws are now bounded by `mWidth` (65) regardless
  of how many times `step()` is called.
- **`finish()`** forces the closing frame only when `mStep < mWidth`, so a bar that already reached
  100 % is not redrawn.
- **Guards:** `mNumSteps > 0` (the old code divided by it unconditionally), the width product is
  computed in `index_t` and clamped to `mWidth`, so an overshooting caller can neither overflow the
  multiplication nor run past the end of the field.
- **Cursor:** `reset()` hides the cursor and only `finish()` gave it back. A destructor now calls an
  idempotent `show_cursor()`, so a bar destroyed without `finish()` no longer leaves the terminal
  blind. This closes the ordinary-teardown half of the leak routed to Christian in the DR-107 round;
  **the `BELFEM_ERROR` half is still open** — `error_abort()` calls `std::abort()`, so no destructor
  runs on the "mesh is too coarse along this loop" exit from `clean_spfa`. Restoring the cursor from
  `assert.cpp` is Christian's call and was not touched here.

## Gates

Executed:

- `-fsyntax-only -std=gnu++17` clean.
- Standalone probe (`Progressbar` + `Logger` + `assert` + `stringtools`, no build tree), stdout on a
  pipe and on a pty via `pty.fork()`, chunk arrival times stamped: old = 1 chunk at exit, new = one
  chunk per frame. This is the measurement quoted above.
- Edge probe: `mNumSteps == 0` (no division), 100 000 steps (65 frames, i.e. bounded by `mWidth`),
  `step( 999 )` on a 10-step bar (clamped, one frame), destruction without `finish()`
  (hide/show cursor balanced 4/4 across all four cases).

Not executed: `make check-fast`, and a real `mpirun hphirun` on a cohomology deck.

## Correction, same session: the app-side flush is necessary but NOT sufficient

The section above originally closed by treating the launcher side as an unmeasured "if". It was
then measured, and the answer is the opposite of what that hedge implied: **under Open MPI 5.0.10
the flush alone does not fix the bar.** An MPI probe (rank 0 draws the bar, `mpirun -np 2`) gives

```
mpirun               probe_mpi_old   3 chunks   ( everything at exit )
mpirun               probe_mpi       3 chunks   ( everything at exit )   <- flush does not win
mpirun --stream-buffering 0   probe_mpi_old  86 chunks
mpirun --stream-buffering 0   probe_mpi      64 chunks
mpirun --stream-buffering 1   either          3 chunks   ( line buffering cannot help a bar with no newline )
```

`strace -f -tt` locates the two halves exactly. The rank writes its 76-byte frames one per
`fflush`, correctly spaced 20 ms apart — the class-side fix works. The single 3043-byte delivery
comes out of the **prterun** process, at job end: Open MPI's I/O forwarder is holding the stream.
So the two fixes are independent and both are needed — the flush gets the frame out of the rank,
`--stream-buffering 0` gets it out of the launcher.

One trap worth recording: the flag is **not** reproducible through the environment.
`mpirun --stream-buffering 0` sets exactly `OMPI_MCA_ompi_stream_buffering=0` in the child and
nothing else (diffed `mpirun -np 1 env` with and without), yet exporting that variable by hand
does not stream — verified forwarded to the rank, 3 chunks in three repeats against ~90 for the
flag. The flag therefore also carries a forwarder-side setting that does not travel in the
environment. Use the command-line option.

## Launcher — `examples/scripts/Allrun`

`--stream-buffering 0` added to the `mpirun` invocation, next to the existing `stdbuf -oL`
(which addresses the *other* half: line-oriented output; it can do nothing for a bar that never
writes a newline). Capability-probed in the style the script already uses for `stdbuf`, with a
NOTE on the failure path. The probe is `"$MPIRUN" --stream-buffering 0 --version`, which
validates the option and exits **without spawning ranks** — ~4 ms, and safe to run inside an
allocation (an unsupported flag exits 213, a supported one 0). Not added to the `srun --mpi=pmix`
path: there is no equivalent there.

### Pre-existing bug found in the same block, fixed

```bash
status=${PIPESTATUS[0]}
tee_status=${PIPESTATUS[1]}     # <- always unbound
```

Bash rewrites `PIPESTATUS` after **every** command, and a plain assignment is a command. Reading
`[0]` into `status` therefore leaves a one-element `PIPESTATUS` behind, so `${PIPESTATUS[1]}`
trips `set -u` and kills the script *at that line* — taking the `exit "$status"` below it with it.
Net effect at HEAD: `./Allrun` exited 1 on a perfectly successful run and 1 on a failed one, i.e.
the solver's status never propagated at all and the careful comment above those lines described
something the code could not do. Fixed by copying the whole array in one command
(`pipe_status=( "${PIPESTATUS[@]}" )`).

Verified with a stub deck and stub launchers: solver exit 0 → `Allrun` exit 0; solver exit 42 →
`hphirun exited with status 42` and `Allrun` exit 42; HEAD → exit 1 in both cases. Launcher
accepting the flag → `--stream-buffering 0 -np 2 stdbuf -oL hphirun`; launcher rejecting it →
flag omitted plus the NOTE. `--dry-run` unchanged and still probes nothing.

`doc/parallel_execution.md` updated: the flag added to the recommended command line, and a
`--stream-buffering 0` note beside the existing `stdbuf -oL` one explaining why they are not the
same knob.

## Documentation

`src/core/doc/core_usage_guide.md`: the `Progressbar` implementation block updated (private
`draw`/`flush`/`show_cursor`, `mWidth`, bounded frame count), a new "Output buffering" note under
the class, a `Progressbar` row in the MPI behavioural-differences table, and the MPI best-practice
entry extended. One stale pitfall corrected while there: "calling `step()` after `finish()` causes
assertion failure in debug builds" was never true — there is no assertion in the class; it is now
(and was always) a silent no-op.
