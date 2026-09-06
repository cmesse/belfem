# DR-96: BELFEM_DUMP_SYSTEM_STEP off-by-one and atoi hardening

**Date:** 2026-08-28
**Purpose:** Fix the env-gated system-dump hook so the LIVE-vs-RESTORE differential it was
built for actually fires when driven from the WARM RESTART banner
**Module:** fem/kernel

## What was wrong (DR-96, all three findings re-confirmed against the working tree)

1. **Off-by-one against the banner.** The WARM RESTART banner printed the PRE-increment
   `mRunningTimeStep`, while the timestep loop increments (`cl_FEM_Controller.cpp:214`, magnetic
   path `:346`) before any assembly, and the `BELFEM_DUMP_SYSTEM_STEP` filter compared against
   the POST-increment value. Live-log evidence from the DR-92 campaign: banner "timestep 224",
   first step header "BDF5 225". Copying the banner's number targeted a step whose assemblies
   had already happened — the differential silently dumped nothing.
2. **Bare `( uint ) std::atoi`.** Empty or non-numeric text folded to 0 (dumps only at step 0,
   silently); a negative wrapped to a huge `uint` (never matches, silently); overflow is UB.
3. **Textual duplication.** The 8-line block existed twice (`iterate_coupled` + magnetic path),
   each with its own function-local `static int tDumpCount`, making the 4-dump budget per-path.

## The fix (one pass, probe exemption — env-gated diagnostic instrumentation, zero production impact)

**Direction call on the off-by-one:** the filter's numbering was the *correct* dialect — it
matches the step headers ("BDF5 225") the LIVE arm's operator copies from. Teaching the filter
to accept banner numbering would have broken live-log targeting the other way. So the banner
was moved to the filter's dialect instead: it now prints the step the first assembly will carry,

```
 resuming at t = 3100.0000 ms, next step 225, delta t = 50.0000 ms
```

with a comment stating why (`the loop increments mRunningTimeStep before assembling`), so the
value is safe to copy into `BELFEM_DUMP_SYSTEM_STEP` verbatim.

**Parsing:** digits-only validation (`strspn` over `"0123456789"` covering the full string,
empty rejected) followed by `strtoul`, behind one `BELFEM_ERROR` — per the philosophy doc's
tier rule (missing environment / bad input is always-active, and the check runs once per
assembly, not per element). `strtoul` alone is insufficient: it silently accepts and wraps a
leading `-`, which is exactly failure mode 2.

**Deduplication:** both copies folded into a new private
`Controller::dump_system_if_requested()` (declared beside `impose_voltage_bcs()`, which shares
the same "used by both iterate paths" role), called immediately after
`compute_jacobian_and_rhs()` on both paths. The `static int tDumpCount` now lives once, so the
4-dump budget is global across paths, and a third caller costs one line.

## Evidence

- Full-flag `-fsyntax-only` compile of the TU with the build tree's own
  `flags.make` (defines + includes + `-Wall -Werror -std=gnu++17`): clean.
- **Reviewed, not verified.** Nothing was built or run. The runtime gate is owed and recorded
  in the register row: on a warm restart, `BELFEM_DUMP_SYSTEM=1
  BELFEM_DUMP_SYSTEM_STEP=<banner "next step" value>` must produce `sysdump_0.hdf5` at that
  step, and a mistyped value (`abc`) must abort with the new `BELFEM_ERROR`.

## Gate run (same session, later): both halves PASSED — row struck

Christian directed the rebuild into the dedicated sandbox `cmake-build-claude/` (Debug,
`USE_TEST=ON`). Full incremental rebuild clean; `make check-fast` 10/10 (11.1 s).

Before the gate ran, a parallel session extended the helper — the gate therefore verified the
*merged* form, not my original: per-field dump budgets as Controller **members**
(`mDumpCountMagnetic` / `mDumpCountThermal`; one shared counter let whichever field assembled
first starve the other, and a function-local static is wrong the moment a second Controller
exists), a field tag in the filename (`sysdump_magnetic_0.hdf5`), a nullptr guard, and the
thermal assembly sites (`iterate_coupled` thermal half + `iterate_thermal`) instrumented too —
four call sites now. The DR-96 core (banner numbering, loud validation) is intact inside it.

Setup: scratchpad copy of `examples/circuit` (9372 nodes, 3D thin-shell tapestack, MUMPS,
shipped `tapestack3d.bfm` so no cut generation), warm-restarting from the example's shipped
`memdump.hdf5` (t = 1.5125 ms, dumped `running_timestep` 7, dt 0.1 ms), `simulation time`
shortened to 1.8 ms in the copy. The shipped example directory was not touched. Serial run,
`mpirun -np 1`, `OMP_NUM_THREADS=1`.

- **Numbering**: banner printed `resuming at t = 1.5125 ms, next step 8`; the run's first
  step header is `BDF1 8`. Banner and step headers agree — the original defect (banner 224,
  first header 225) is gone.
- **Positive half**: `BELFEM_DUMP_SYSTEM=1 BELFEM_DUMP_SYSTEM_STEP=8` produced
  `sysdump_magnetic_0.hdf5` and `_1.hdf5` — exactly step 8's two assemblies. Steps 9 and 10
  then ran to completion and added **nothing** (budget of 4 not exhausted), so the filter hits
  the target step and does not leak past it. Exit 0.
- **Abort half**: `BELFEM_DUMP_SYSTEM_STEP=abc` aborted at the first assembly with
  `invalid BELFEM_DUMP_SYSTEM_STEP 'abc' - expect a plain nonnegative integer`
  (backtrace through `Controller::dump_system_if_requested` ← `iterate_coupled`), exit on
  SIGABRT — loud, as designed.

**Verified** in the evidence-ladder sense: an executable gate ran and passed. Register row
struck; `[P]` count 25 → 24, recounted mechanically. Logs kept in the session scratchpad
`dr96_gate/` (`gate_mistyped.log`, `gate_positive.log`).

## Process notes

- Christian ruled: fix now under the probe exemption (no plan+audit → code+audit round);
  the row is marked **fixed, gate owed** rather than struck, following the register's
  strike-after-gate precedent (DR-92, DR-115).
- Probe-removal did NOT apply: DR-92's O2/O6 (Codex's solver-entry-state theory) are still
  open, and this hook is the named decisive instrument for the live-vs-restored
  assembled-system comparison.
- Files touched: `src/fem/kernel/cl_FEM_Controller.cpp` (helper definition, two call sites,
  banner), `src/fem/kernel/cl_FEM_Controller.hpp` (declaration), `todo/debt_register.md`
  (DR-96 status + gate). Uncommitted; the controller also carries parallel sessions'
  uncommitted work from this same day.
