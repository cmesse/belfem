# DR-129: Kernel::mController initialization — plan, two audit rounds, fix landed, gates run

**Date:** 2026-08-30
**Purpose:** Session record — DR-129 (the register's only P1) taken from row to code-audited fix
in one day: plan + xhigh plan audit + implementation + high code audit, unanimous approve, then
ALL THREE run gates green the same evening (debug suite, release suite, coupled thermal smoke).
Ready to strike.
**Module:** `src/fem/kernel`, `src/fem/maxwell`, `src/fem/thermal`, `tests/fem`

## What happened

1. **Plan** (`todo/dr129_kernel_controller_init.md`): re-derived the row from the tree. Two
   corrections to the row before any auditor saw it: the consumer lives in `src/fem/maxwell/`,
   not `src/fem/iwg/` as filed, and `has_controller()` did not exist anywhere in `src/` — the
   row's fix shape assumed it did. Reader census: **seven** `Kernel::controller()` call
   expressions, not the one the row named.

2. **Plan audit, both vendors at xhigh** (safety-boundary tier — subject is lifetime; depth
   chosen by Christian). The round earned its cost:
   - **Codex found the blocker the row never knew about:** `cl_ThermalFactory.cpp:138` built the
     thermal kernel as `make_shared<Kernel>(Kernel(mKernelParameters))`. User-declared dtor + no
     copy/move members → C++17 suppresses the implicit move, the implicit **copy** ctor runs, and
     it copies the indeterminate `mController` AND `mMyCommIndex` — on **every coupled run**
     (`belfem.cpp:157`, `hphiTrun.cpp:80`). Claude's calibration trace: on the live path the
     temporary owns nothing (`KernelParameters(Kernel*)` selects the non-owning
     `distribute_mesh` branch), so no double-free — formal UB only. But the *unused* mesh-based
     `ThermalFactory` ctor would have made the temporary own its comm tables in parallel:
     copy + destroy = use-after-free. Same implicit-copy family as DR-141.
   - **Both vendors independently refuted two plan claims:** `have_thermal()` had **zero
     consumers** (the `:266` write was a dead store — the "silent wrong physics" severity story
     died), and the battery order-shuffle gate was wrong (post-fix the member is
     deterministically null; no in-tree shuffle target exists).
   - **Unanimous O1 recommendation:** (c) — `BELFEM_ERROR` in `controller()` plus the predicate.
     Deciding reason: `MaxwellData` binds `const real & mTime` to `controller()->time()` in a
     member-initializer list and *cannot* take a local guard; the accessor is the only
     release-loud site.
   - Grok's by-catch: `IWG_Maxwell::mFunMKF` also lacked `= nullptr`; and the thermal
     `initialize()`-before-`set_thermal_kernel` window is a recorded future hazard (do NOT move
     `MaxwellData` construction earlier).
   - Claude's own miss, caught in reconciliation: the caller grep pattern `TypeName(` missed
     stack declarations (`ThermalFactory tTFactory(...)`) — match the type name, not the paren.

3. **Christian's rulings:** O1 = (c); R1 scope = all three members; implement now.

4. **Landed** (7 files): `mController = nullptr`, `mMyCommIndex = 0`, `mFunMKF = nullptr`;
   `has_controller() const` (inline, `is_master()` precedent); `BELFEM_ERROR` in `controller()`;
   both-arm guarded `mHaveThermal` assignment through a `tKernel` local; ThermalFactory
   constructs in place; DR-46 fixture `link_to_group` re-enabled; `GhostElementContract` gains
   `ASSERT_FALSE( tStack.iwg()->have_thermal() )` — `have_thermal()`'s first consumer ever.
   All four TUs pass `g++ -fsyntax-only` with the build tree's `flags.make` defines (needed
   `-isystem /opt/scls/gcc/include` + oneapi MKL includes on top — flags.make alone does not
   carry the TPL include roots).

5. **Code audit, both vendors at high:** APPROVE / APPROVE, zero logic defects, `git status`
   clean after the Grok round. Both walked the re-enabled fixture path line by line
   (`mTimeStepMatrices` allocated in the `IWG_Timestep` base ctor; Ghost FieldList/dof arms
   wired; nothing after the switch touches the controller). Two Grok comment nits applied.
   Recorded residuals, deliberate: implicit `Kernel` copy ctor still exists (Rule-of-5 hygiene,
   another day), no `EXPECT_THROW` on the bare-`controller()` abort, empty-group return still
   skips the `mHaveThermal` write (unobservable).

## Gates

- **Debug `make check` — RAN GREEN 2026-08-30 (Christian's build), same session.**
  `GhostElementContract [ OK ] (3 ms)` in the 15:55-rebuilt `test_fem`
  (`LastTest.log:17041-17092`), suite green; the new assertion string is present in the binary,
  so the pass is the new test. The re-enabled controller-less `link_to_group` is **VERIFIED** at
  the named-suite rung.

## The release gate, and a false start on it

The release sub-gate RAN the same evening and PASSED: `build/` reconfigured to `USE_DEBUG=OFF`
(`CMAKE_BUILD_TYPE=Release`, `-O2 -DNDEBUG`), `test_fem` relinked 17:13, `make check` 17/17 suites
green, `InterfaceOrientation.GhostElementContract [ OK ] (3 ms)` at
`build/Testing/Temporary/LastTest.log:13062-13113`. This is the decisive run — under `NDEBUG` the
accessor's old `BELFEM_ASSERT` does not exist, so a pre-fix binary had nothing between it and the
dereference of the indeterminate member.

The mode was confirmed **from the linked binary, not from the cache**: an ASSERT-only string
(`"Don't know what to partition"`) has zero hits in `test_fem`, proving `BELFEM_ASSERT` really is
compiled out; while `"Controller is not set for this kernel"` (R4's `BELFEM_ERROR`) and
`have_thermal` (R5's contract assertion) are both present, so O1(c)'s always-active guard survives
into release exactly as intended.

That check earned its keep. A first attempt at this gate, 35 minutes earlier, was reported as
release and was not: `build/CMakeCache.txt` still read `USE_DEBUG:BOOL=ON` and its `flags.make`
still read `-Og -g -DDEBUG`. `USE_DEBUG` is a plain `option( ... ON )` (`CMakeLists.txt:90`) with
no forcing, so it takes the **cached** value unless `-D` appears on the `cmake` invocation itself —
and a bare `make` that triggers a reconfigure re-reads the old mode silently. Nothing in the 17/17
output distinguishes the two builds. **Confirm build mode from `flags.make`, or from an ASSERT-only
string in the binary; never from the tree's name or from what the configure was meant to say.**

## The coupled smoke, and why it was closed early

Third sub-gate RAN the same evening: `tape_quench_usermat` under `prterun -np 8`, from a release
build, warm-restarting at t = 11.0000 ms / step 203. Genuinely coupled — the run announces
"Thermal solver section found in the input file: solving the coupled h-ɸ/T problem" and prints both
Magnetic and Thermal Picard blocks every step. It carried the rebuilt `ThermalFactory` construction
site at startup and ran **483 converged coupled timesteps** (step 203 to step 686) at
~8e-11 magnetic / ~3e-13 thermal, with **zero DR-129-class failures**.

Closed at that point rather than waiting hours for the deck to finish, and the reason is worth
stating because it generalizes: **the R8 site executes once, at startup.** It was behind the run
inside the first minute. Every subsequent timestep exercises the deck's physics, not this row's
defect — so a longer run buys no evidence for DR-129, only for whatever else the deck is being run
to answer. Match a gate's duration to when its subject actually executes.

The mode was again confirmed from compiler artifacts: the run's tree is *named* `cmake-build-debug`
but was reconfigured `USE_DEBUG=OFF` (`-O2 -DNDEBUG` in `flags.make`, `bin/belfem` relinked 17:18),
and the solver process started 17:19:57 — after the link, so the running binary is the new one.

By-catch, deliberately unfiled and belonging to no DR row here: 8 `MUMPS soft fail: error code -20`
and 2 `PARPACK failed to compute the eigenvalue conditioning estimate`, all in the
eigenvalue/conditioning instrumentation path, every one recovered (26 further timesteps succeeded
after the last). Handled soft-fail paths, not aborts.

## Closed

DR-129 was STRUCK and archived to `todo/debt_register_closed.md` on Christian's ruling the same
evening, and the plan file moved to `todo/closed/`. Filed 2026-08-28 as by-catch of the DR-46
fixture work, closed 2026-08-30: plan, xhigh plan audit, implementation, high code audit, and all
three run gates inside one day. **The strike clears the debt register's only P1 row.** Nothing
executable is owed.

## Lessons

- The xhigh plan round found a production-path defect (D3) that three careful readings of the
  *accessor* call graph could not, because the read was not a call: **enumerate writers and
  copies, not just named reads**, when the subject is an uninitialized member.
- A grep for `TypeName(` does not find `TypeName tVar(...)` declarations.
- `have_thermal()` going from zero consumers to one (as a test assertion) is the cheapest kind
  of contract pin: it converts a dead store into an observable.
- An append-mode log is not a run. `out.txt` here accumulates every run of the deck, and the first
  count taken off it ("523 steps") silently spanned two of them — a fifth run had warm-restarted
  from an earlier memdump while the count was being read, so its steps were added to the fourth
  run's and its lower step numbers made the range look wrong in the other direction. Segment such a
  log on its start banners before counting anything, and pin the run to a process start time and a
  binary mtime. Corrected here to 483 steps (203 to 686) for the run that used the gate binary.
- Match a gate's duration to when its subject executes. A once-at-startup fix is fully exercised
  in the first minute of a multi-hour run; waiting for the run to end confuses "the deck finished"
  with "the defect was tested".
- A build tree's *name* and its cache's *intent* are not evidence of its mode. When a gate's whole
  meaning depends on the mode (here: "does this pass with the assert compiled out?"), verify the
  mode from an artifact the compiler produced — `flags.make`, or the presence/absence of a
  tier-specific string in the linked binary.

**Session:** Claude (Opus 5 → Fable 5) with Codex `gpt-5.6-terra` and Grok `grok-4.6`;
rulings by Christian. Exchange: `tmp/ai_exchange/dr129_kernel_controller_init.md` (distilled
here and in the plan; GC-eligible per §10).
