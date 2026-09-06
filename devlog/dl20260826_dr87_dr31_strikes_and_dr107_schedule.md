# DR-87, DR-31 and DR-53 Struck by Ruling; DR-107's Gate Scheduled for 02:00

**Date:** 2026-08-26
**Purpose:** Record two strikes taken on Christian's ruling, the findings that produced the DR-87
ruling, and the scheduled autopilot job for DR-107.
**Module:** cross-cutting (todo/, fem/kernel, sparse)

## What was asked

Run DR-87 then DR-107 so both could be struck, on autopilot at 02:00. Preparing DR-87's run turned
up three findings that made the run the wrong move, and Christian ruled: **strike 87 on the
existing evidence**, and **strike DR-31** with it.

## Why DR-87 should not have been run — the finding that mattered

The obvious objections were resources: **32 GiB available against a 47 GiB peak** (ParaView holding
6.0 GiB, CLion 5.8 GiB) — the exact shape that killed DR-77 at 51.6 GiB with swap exhaustion — and
**disk at 97%**, 37 GiB free of 945 GB.

The finding that actually decided it is different and would have wasted the night silently:

> **The deck is no longer in the configuration DR-87 was measured in.** The row's measurement was
> **STRUMPACK magnetic + PETSc thermal**. `cmake-build-debug/tapestack3d/input.conf:15-18` now reads
> **`linear magnetic { library : mumps }`** — switched during the 2026-08-24 four-run
> BDF1/BDF5 × mumps/strumpack × gauged/ungauged comparison that closed DR-02.

A run tonight would have measured **MUMPS** allocator behaviour, not the STRUMPACK regime the
`malloc_trim` mitigation was gated against, and would have walked into DR-106 (MUMPS −9 workspace
exhaustion answered by a timestep cut, still open). It would have produced hours of numbers that
could not support the strike they were run for. **A gate re-run against a drifted deck is not the
same gate** — worth adding to the standing traps if it recurs.

Two smaller points also argued against: DR-87's gate had **already PASSED and been verified by
execution**, and what remained was a **watch with no defined pass threshold** — the row says the
risk is a creep past ~53 GiB "over hours" but never says how many flat steps are enough.

## The two strikes

Both are recorded **DR-42/49-style exceptions**: design work finished, residue is not a decision.
Applied per the register's own rules — ID and description struck, **status column left unstruck**
because that is where the closure evidence lives, exception recorded rather than applied silently,
live triage tag dropped (retired rows do not carry one), and the row **filed into the archive** in
the same action.

**DR-87.** The gate's evidence stands: 13.7 GiB released at a step boundary against ~3 GiB pre-trim,
trough 46.04 → 31.86 GiB, in-step PEAK flat at 47.19 / 47.24 / 47.38 GiB against pre-trim peaks that
climbed 50.4 → 57.8 GiB over hours. **The watch survives as a standing caution, not as debt:** peak
past ~53 GiB means the ratchet is back, next rungs `mallopt(M_MMAP_THRESHOLD)` or jemalloc. PEAK is
the metric; the trough is sampling-limited at 60 s. It rides the next long STRUMPACK campaign. The
row's warning stays live — **do not "fix" a symptom here with `krylov method : preonly`**, because
the factorization is degraded to ~4 digits and the outer GMRES is what reaches 1e-11.

**DR-31** is the cleanest exception the register has had: its gate is not merely unrun but
**unrunnable by design**. The fix is written, committed (`0671b29a`), tree-verified, and covers
`edge_h` and `face_h` through one entity-type-agnostic path. Reaching it needs an order-2 mesh on
≥ 2 ranks, and the 2026-08-24 chain found and fixed **DR-102** and **DR-103** on the way before
terminating at an always-active `BELFEM_ERROR( max_element_order() == 1, "Not implemented for higher
order" )` in `create_hanging_edges_and_facets` (`cl_MaxwellFactory.cpp:1399`). Every in-tree 3D h-φ
deck computes cohomologies, so no order-2 deck can reach the Maxwell solver. Order-2 is not a 1.0
feature and the framework refuses it loudly, so the row's own ask resolves itself. **Struck ≠
verified applies with full force: that code has never executed.**

## DR-53 — struck on the same principle, with a check that mattered

> Christian: *"DR-53 is stricken. Changing the default setting should be considered experimental,
> not for user audience."*

This settles the **default question** the row had carried since 2026-08-09 — *keep the flag `false`
versus remove it altogether* — as **keep it, classified experimental, and do not expose it.**

The check that made the strike safe rather than merely convenient: **the fuse flags are not
deck-settable.** `mFuseEdges` and `mFuseEdgesWhenHavingSideConnectors` are hardcoded `false` at
`cl_ThinShellFactory.hpp:353-354`; there is **no matching string literal anywhere in `src/`** and no
entry in `doc/input_file_reference.md` or `doc/input_schema.yaml`. So the two-artifact input-contract
rule **does not fire here**, and — the substantive point — no user can reach the ~2e-6 residual floor
or the 100 ns Δt collapse. Getting there requires editing and rebuilding the source, which is exactly
what makes it an experiment rather than a supported configuration. That is what closes it as *debt*.

Two things deliberately survive the strike:

- **The standing constraint.** Neither flag may be re-defaulted to `true`. The reproducer is
  `sidecoatings` at t = 2.0825 ms; fuse-off passes that cliff warm and cold.
- **The physics half of O1** — graded two-branch authority versus the HEX8TB wall element as the
  proper physics — is **not** settled by this ruling. It stays in
  `todo/side_edge_fusing_cut_aware_plan.md` as a *research question, not debt*. The flag site already
  records the position (fusing is the mathematically cleaner continuity statement, yet measured runs
  converge slower for no better result — it appears to overconstrain the rim, like weakly enforcing a
  B·n = 0 the formulation already satisfies), so nobody re-runs the experiment blind.

The plan's tracker was updated in the same turn per standing rule: O1's default half struck and
marked resolved, and the living Status line rewritten to say the debt is closed and only the research
question remains.

Register integrity after all three: **40 live + 70 archived = 110 IDs, no gaps, no overlap, no struck
rows in the live file, every live row tagged.**

## DR-107 scheduled

`at` **job 19, Thu 2026-08-27 02:00** → `~/belfem_autopilot/dr107_gate.sh`. Chosen over CronCreate,
which is session-only and fires only while the REPL is idle — too fragile for a 6.8-hour wait.

The script guards before it acts: refuses if a compile or solver is already running in the shared
tree, and refuses below 10 GiB free disk. It records HEAD and the full uncommitted bundle under
gate, builds with `-j16`, runs `make check`, extracts the homology lines and the ctest totals, and
appends a verdict to `~/belfem_autopilot/dr107_SUMMARY.txt`. **If the build fails it aborts before
testing and says so** — the bundle is large and has never been compiled as a whole.

The build is also the regression gate for that whole bundle: the homology cleanup (−1201/+135),
DR-108's `uint8_t` fix, DR-30, DR-101, DR-102, DR-103.

## Backup taken before anything

`cmake-build-debug/tapestack3d/memdump.hdf5` is the **sole surviving dump** (t = 7.1 s, step 3349,
218 MB). Copied to `memdump_t7100_step3349_KEEP.hdf5` before any run was contemplated. DR-89 and
DR-90 both lost their gates to deleted dumps; this one should not be swept.

## Owed

- DR-107's verdict, from `~/belfem_autopilot/dr107_SUMMARY.txt`, and the strike if homology reads 9/9.
- DR-95's other half — the `cl_FEM_Controller.hpp:53-56` comment demotion.
- Nothing was built or run this session; the 02:00 job is scheduled, not executed.
