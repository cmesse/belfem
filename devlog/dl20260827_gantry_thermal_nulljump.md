# Gantry np=8 SEGV: null thermal kernel dispatch; `belfem` unified executable

**Date:** 2026-08-27 (evening session)
**Topic:** gantry `hphiTrun` first-timestep segfault — root cause, and the `belfem`
executable that makes the failure class unrepresentable
**Exchange:** `tmp/ai_exchange/belfem_unified_executable.md` (pre-registration, Codex +
Grok plan audits, reconciliation, code audits), diff in
`belfem_unified_executable_code.diff`

## Symptom and the false trail

`mpirun -np 8 hphiTrun` on `cmake-build-debug/gantry` (2D, 464 thin-shell tapes, iron
yoke, first run at 12 K) segfaulted on every launch in the first iteration of the first
timestep — PETSc's handler caught SEGV on ranks 4 and 0 and aborted before a core could
be written. Christian's debug script appeared to work, suggesting a race. It was not:
the debug attempt ran **hphirun**, the crash was in **hphiTrun** — two different
executables, and the deck (no thermal solver sections, no thermal BCs) was only ever
valid for the first one. Build-tree forensics were a second distraction worth recording:
the two binaries were linked 65 minutes apart (18:15 vs 19:17) across an actively edited
tree, so "which code is in which binary" had to be reconstructed from object-file
timestamps before any diff-based suspicion made sense. None of the day's DR fixes was
involved.

## Root cause (verified by execution — valgrind memcheck, np=8)

All ranks: *jump to address 0x0* through the function pointer `mFunMKF` in
`IWG_MaxwellThermal::compute_mkf()` (`src/fem/thermal/cl_IWG_MaxwellThermal.cpp`),
called from the first thermal element assembly. The pointer was declared without an
initializer and assigned only inside `link_to_group()`'s domain-type switch. The
observed 0x0 was luck — uninitialized is UB, not nullptr. Which group reached
`compute_mkf()` without a kernel remains **open**: the fem Block constructor keeps
`number_of_elements()` and the element Cell consistent (Grok's refutation of the
early-return theory), so the trigger is not the zero-element guard alone. The hardening
below makes the next occurrence self-identifying.

## Ruling (Christian)

- One unified executable `belfem`: the deck decides the physics. An unlabeled
  `linear thermal` or `nonlinear thermal` section under `solver` selects the coupled
  h-ɸ/T problem; without one, magnetic-only. A deck with no thermal setup must never
  reach thermal assembly.
- `initial conditions/temperature` is the material-law reference, not a thermal-solve
  request; absent, **77 K is assumed in both modes** with a console note (extends the
  existing `hphirun` behaviour).
- `hphirun`/`hphiTrun` stay untouched through the release; deprecation is post-release.

## What landed

- **`src/executables/belfem.cpp`** (replaces the scaffold): `hphirun`-shaped main +
  optional thermal. Collective deck discriminator via `InputFile::section_exists` (a
  rank-split decision would deadlock the factory collectives); startup line names the
  selected mode; 77 K default applied on all ranks between the MaxwellFactory
  constructor and `create_magnetic_kernel()` so the temperature mesh global is
  published; thermal BCs without a thermal solver section are refused as inconsistent;
  `set_thermal_kernel()` is only ever called with a real kernel; the segregated loop is
  reachable only in thermal mode. Registered in `src/executables/CMakeLists.txt`.
- **Pointer hardening** in `cl_IWG_MaxwellThermal.{hpp,cpp}`: `mFunMKF = nullptr` at
  declaration; the dispatch switch runs before the zero-element early return (so the
  IWG-wide pointer never carries the previous group's kernel) with `default:` resetting
  to nullptr; an always-active `BELFEM_ERROR` fires for any *assembled* group without a
  kernel, naming group id and domain type; a hot-path `BELFEM_ASSERT` backstops
  `compute_mkf()`. Empty groups of unhandled types keep skipping silently, as before.
- **Input contract + docs, same session:** `doc/input_file_reference.md` (§4.1 mode
  selector, §10 temperature default) and `doc/input_schema.yaml` (belfem_mode_selector
  notes, initial-conditions note); CLAUDE.md executable inventory (the claims script was
  already red because the scaffold existed unlisted — now 34/34);
  `src/executables/doc/README.md`, `examples/README.md`,
  `src/fem/maxwell/doc/README.md`.
- **Launcher switched over** (Christian's follow-up ruling — `hphirun`/`hphiTrun`
  retire soon): `examples/scripts/Allrun` now always launches `belfem` and its
  in-script thermal discriminator is retired (the executable does that job); the old
  binaries stay reachable via `EXECUTABLE=` in `run.conf` until removed.
  `examples/scripts/README.md` updated to match.

## Audit trail

Three-vendor round on the plan (Codex: accept with three corrections; Grok: accept with
five requirements; all citations spot-verified, two by execution). Key audit
contributions baked into the code: the always-active post-switch error (asserts compile
out in release), the `hphirun`-based main shape (`hphiTrun`'s segregated branch calls
`initialize_thermal()` unguarded; `set_thermal_kernel` dereferences its argument), the
empty-group whitelist in the moved switch, and the 77 K publish-order constraint. Code
audits dispatched on the diff (verdicts land in the exchange).

## Status / owed gates

Both edited TUs are syntax-clean under the build tree's `flags.make` flags;
`check_doc_claims.py` 34/34. **G1 verified by execution (2026-08-28, Christian's
build + run): `belfem` builds, selects magnetic-only on the gantry deck, and runs
the formerly-crashing np=8 configuration.** Remaining gates:

- [x] Christian: `make` (`make belfem` builds the executable, `bin/belfem`; the
  library target was renamed `belfem_lib` to free the name — archive still
  `libbelfem.a`) — done 2026-08-28
- [x] G1 — gantry deck, `belfem`, np=8: selects magnetic-only, runs timesteps (the
  formerly-crashing configuration) — green 2026-08-28
- [ ] G2 — a coupled deck (tapestack3d / tape_quench): selects coupled mode, matches an
  `hphiTrun` reference over the first timesteps
- [ ] G3 — temperature-stripped deck: prints the 77 K notice and runs (77 K is
  over-critical for the gantry tape — the gate is "runs", not "matches 12 K physics")
- [ ] G4 — `make check`
- [ ] Codex language sweep over the touched user-facing prose

**Release note:** after the hardening, a hand-launched `hphiTrun` on a magnetic-only
deck no longer segfaults — it starts an expensive, unwanted coupled solve.
`belfem` (or `hphirun`) is the right launcher for such decks; `hphiTrun`'s always-on
thermal factory is deliberately unchanged during the freeze.
