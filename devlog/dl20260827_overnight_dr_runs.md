# Overnight DR autopilot: DR-107 gate rerun, Batch B, DR-105 deck reconstruction

**Date:** 2026-08-27 (evening session; results land overnight)
**Purpose:** Discharge the register's executable residue while the machine is free —
Christian's go: "do these runs on auto pilot", 8 MPI procs approved.

## Starting state

- The 02:00 `at` job for DR-107's gate **aborted itself** — its guard found a solver running
  (`~/belfem_autopilot/dr107_SUMMARY.txt`). The guard worked; the gate never ran.
- The homology bundle it was meant to gate **was committed anyway** (`d8e5fa20`/`dcca23cc` — tree
  clean, `cl_SideSetFactory.*` gone), so the `make check` homology 9/9 gate is now owed as a
  post-commit regression check. Struck is not verified; this run is the verification.
- `hphiTrun` in `cmake-build-debug` was linked 02:15, the matmerge merge landed 02:59 — the shipped
  binary predated HEAD's material headers (relevant below).

## What was launched

`~/belfem_autopilot/night_20260827.sh`, background, sequenced:

1. **Phase 1 — DR-107 gate:** `dr107_gate.sh` verbatim (build at HEAD + `make check`, homology
   lines extracted). Also the first build of HEAD after the matmerge merge.
2. **Phase 2 — Batch B** (`run_gate_batches.md` §2): `examples/RLC_Circuit/Allrun` and
   `examples/circuit/Allrun`, 2 h timeout each. RLC at its default rank count doubles as the
   post-merge re-confirmation of the committed DR-100 fix (`7a211d59`).
3. **Phase 3 — DR-105 plugin rebuild**, **3.5 — 5-min smoke probe**, **4 — the long run**:
   reconstructed tape_quench deck, `mpirun -np 8 hphiTrun` via the shared `examples/scripts/Allrun`
   (run.conf: NRANKS=8 fixed, NTHREADS=2), disk watchdog kills the solver below 12 GiB free.

## DR-105 deck reconstruction (the actual work of the session)

DR-105 was `[RUN-BLOCKED]`: the tuned deck died with `build/`. Reconstructed at
`cmake-build-claude/tape_quench_dr105` from three recorded sources — the §6 lesson
("a gate whose recipe is recorded survives the deletion of its artifacts") did its job:

1. **Substrate:** `tmp/examples/Tape_Quench/CustomMat` (January lineage, byte-identical matlib.cpp
   per dl20260823).
2. **The three dl20260823 plugin repairs re-applied:** `buffer` → `custom_buffer` (four sites),
   `set_constant(T_crit, 92.5)` in `hts_init`, jc/n reshaped to `MatFunc3` and registered through
   the 3-dependency `(normB, angleNxB, T)` overload (the 1-dependency overload never reaches
   `mJcFunction`).
3. **The T5 recipe from the DR-105 register row:** strumpack magnetic / petsc thermal, `anderson
   depth : 3`, Picard with `tolerance switch : 1e-4`, chi off (no gauge block), `initial timestep :
   0.002 ms`, `maximum timestep : 0.02 ms`. Everything else stays at the reference deck's values —
   the lost `input.conf` carried Christian's settings and only the register-recorded subset is
   reproducible; deltas from the lost original are possible.
4. **Defect restored to the measured form:** the January `defect.cpp` had the hard 1e4× step active
   (the deterministic limit cycle); the tuned matrix ran the smooth 100× Gaussian, which was present
   but commented out. Swapped.

**Plugin build repairs** (template CMakeLists were stale against today's tree): C++17 (gnu
extensions), include list extended (linalg/blaze, numerics/bezier, MKL includes), and — decisive —
the tree's backend/feature macros mirrored from `material.dir/flags.make` (`-DBELFEM_BLAZE` etc.).
Without them the header chain does not even parse; with a *stale binary* they still crash:

**Measured ABI-drift crash, worth keeping:** the pre-flight smoke (fresh plugin, 02:15 binary)
aborted inside `hts_init` — the first *virtual* call (`set_user_defined_function`, jc, MatFunc3)
landed on `evaluate_polynomial`'s vtable slot ("return polynomial for jc of hts is not defined",
`cl_Material_UserDefined.hpp:390`). Disassembly pinned it: `hts_init+0x9c` is the return address of
the first `callq *%r10`. The plugin was built against HEAD headers, the binary against pre-merge
ones — the exact silent-mismatch mechanism the plugin-ABI memory warns about, this time loud only
because the drifted slot happened to hold an asserting function. Phase 1 rebuilds the binary at
HEAD; phase 3.5's smoke probe re-checks before the long run is allowed to start.

## Outputs to keep (analysis obligation — announced per standing rule)

Until the morning session has read them (keep through **2026-08-29** at least):

- `~/belfem_autopilot/night_SUMMARY.txt`, `night_*.log`, `dr107_*.log` — verdicts
- `cmake-build-claude/tape_quench_dr105/` — `belfem.log`, `allrun_console.log`, exodus outputs,
  `memdump.hdf5` if written; the deck itself is the DR-105 gate artifact now
- `examples/RLC_Circuit/belfem.log`, `examples/circuit/belfem.log`

## Mid-night event: the first launch died at 2% — and that is the headline

Phase 1's plain `make` failed in 40 s: **the DR-37 backend-free gate fired**
(`tests/physics/backendfree/test_MaterialBackendFree.cpp:38`, its own `#error`). Root cause: the
matmerge merge added `#include "cl_Bezier.hpp"` to `cl_Material.hpp:23`, which pulls
`cl_Vector.hpp` — the exact linalg leak the gate (landed the same day, "reviewed, not verified")
exists to forbid. Two same-day sessions collided; the gate worked on its first real `make`.
Consequence: **HEAD does not build under the default configuration** — filed as **DR-114 (P1,
`[CODE]` `[P]`)**. The night script had also fallen through to Batch B on the stale 02:15
pre-merge binary (its only build check was "hphiTrun exists"); those partial Batch B results are
void. Reworked and relaunched: phase 1 is now `make -k` (backendfree is the only tolerated
failure, any other is reported) + freshness check on both solver binaries + direct `ctest`
(the `check` target's command would be skipped on the failed dependency). The DR-107 verdict
tomorrow must therefore say "ctest, with the backendfree dep broken" — not a clean `make check`.

## Overnight results as they landed (appended live)

- **Phase 1 GREEN:** `make -k` at HEAD — the only failing targets are the two backendfree-gate TUs
  (every error line traces to `op_VectorPlus.hpp` through the gate's include); both solvers
  relinked fresh; **ctest 14/14 passed, homology suite Passed** (5.4 s). DR-107's gate is
  discharged in substance, with the caveat that it ran as `make -k` + direct `ctest` because of
  DR-114, not as a clean `make check`.
- **RLC_Circuit GREEN at np=4:** full transient, 754 full-matrix collects (the previously-crashing
  phase), `iv_results.csv` written. DR-100's fix re-confirmed on merged HEAD.
- **NEW FINDING — `examples/circuit` aborts at np=4 (exit 134, ~100 s in):** bounds assert
  `aIndex < this->length()` at `cl_BZ_Vector.hpp:349` via
  `ElectricalCircuit::shift_back()` ← `update_components()` ← `Controller::reset_timestep()` ←
  `iterate_coupled()` — the circuit history shift overruns a `Vector` on the first timestep
  REJECTION. Debug-only catch; a release build would read out of bounds silently. Serial probe
  (5 min, timeout) never reached a rejection and stayed healthy, so serial-vs-parallel is
  undiscriminated; the trigger is the reset path, not steady iteration. Evidence:
  `~/belfem_autopilot/circuit_np4_crash_belfem.log` (np=4 crash, backed up before the probe
  overwrote the deck log) and `circuit_serial_probe.log`. Filed as **DR-115 (P2)** the same
  night — candidate mechanism neighborhood: DR-39's circuit residue / restart_circuit work.
- **DR-105 smoke probe PASSED** on the fresh binary (the vtable mismatch was binary staleness, as
  diagnosed). First phase-4 launch was refused by the launcher — 8 ranks × 2 threads exceeds the
  10 physical cores — relaunched standalone with `NTHREADS=1` (`dr105_phase4.sh`).

- **DR-105 np=8 result (04:05–04:20):** the launcher first refused 8×2 threads on 10 physical
  cores (relaunched at 8×1); the run then made a real measurement — **263 accepted steps to
  t = 0.375 ms in ~15 min, then a hard wall**: Δt collapsed to the 1 µs floor at step 264 and
  the 20 floor retries exhausted (`cl_FEM_Controller.cpp:2267`, the dl20260823 abort). T5 rode
  through the 0.77 ms wall *serially*, so np=8 walls earlier than the recipe's measurement.
  Confounds: rank count, debug binary, and the unrecorded remainder of Christian's deck settings.
  Evidence: `belfem_np8_wall0375ms.log` + `np8_run_artifacts/` in the deck dir.
- **Phase 5 launched (04:2x):** serial discriminator — identical deck, 1 rank, no wall limit,
  cold start (np=8 memdump moved aside). If serial also walls near 0.375 ms, the reconstruction
  differs from the lost deck (tolerances are the suspects); if it rides through 0.77 ms, the wall
  is a parallel effect. One variable changed; no physics keys touched on autopilot.

## Morning after: DR-114 fixed by Christian, struck same day

Christian's fix, landed the next morning: `alpha_custom` + `create_low_temperature_alpha` moved
with their Bezier dependency into a new `SplineLookupTable` class (`Metal` now derives from it);
`cl_Bezier.hpp` and `fn_polyval.hpp` (which includes `armadillo.hpp` — the second, easily-missed
violation) left `cl_Material.hpp`, restoring its pre-merge include set. The 11 metal call sites
were requalified `Material::alpha_custom` → `SplineLookupTable::alpha_custom` (the explicit
qualifier is load-bearing: each metal's 1-arg override hides the inherited 3-arg overload).
Gate verified by execution: `make test_material_backendfree` green, then a plain `make` to 100%
exit 0 — the first clean full default-config build since the matmerge merge. **DR-114 struck**
(Christian's ruling). The DR-37 gate's value is now demonstrated by its first catch, one day
after landing. **DR-107 struck the same morning** on Christian's own `make check` confirmation —
its gate discharged by three independent green runs (overnight ctest, morning make check,
Christian's make check). Batch A is closed in `run_gate_batches.md`.

Separately, the DR-105 serial discriminator answered its question and then some. It ran 9.8 h to
**t = 9.14 ms** at a steady ~1.0 ms/h — 24× past the np=8 wall, 9.6× past T5 — settling the
discriminator: **the np=8 wall is a parallel effect, not a reconstruction defect.** It then
aborted on territory no run of this deck had ever reached: at the quench front,
`rho_piecewise`'s `n > 1.0` assert (powerlaws.hpp:763) fires inside the 90 K < T ≤ 92.5 K
window — the plugin's `hts_n` polynomial falls back to 1.0 at 90 K, but the dl20260823 repair
set `T_crit = 92.5` (the in-tree YBCO builtin value), so the superconducting branch stays
active past the fit's domain. Deck-consistent fix: `T_crit = 90.0`. The guard is debug-only;
release would divide by zero via `1/(n-1)` in the flux-flow branch. T5 never saw this because
no tuning run got past ~1 ms.

## DR-82 gate discharged (same morning, Christian's prompt)

Christian proposed a size-check `BELFEM_ERROR` for DR-82; inspection showed his own 2026-08-16
fix already landed the stronger form (width conversion via `datatype<T>()` memory type + an
always-active class guard, `check_read_datatype<T>()`). What was owed was the verification gate.
Discharged with a scratchpad probe against the prebuilt archives (io/core/comm + hdf5/mpi/petsc,
plus `gComm`/`gLog` definitions the probe must supply itself): (a) u64 dataset → `uint32_t`
scalar, value intact — the formerly-corrupting case; (b) float dataset → integer read, refused
with the intended message. **Rungs 1+2 verified by execution.** On the recommendation that a verified-fixed P1 corruption
headline misstates release risk during freeze week, **Christian struck DR-82** and the hardening
residual (rung 3: fixed-width file types on create; plus the writers' LE-copy BE-host note) was
spun off as **DR-116 (P3, `[CODE]` `[P]`)** — the DR-108→DR-110 spin-off pattern.

## Register maintenance (same session)

- DR-105 retagged `[RUN-BLOCKED]` → `[RUN]`; gate cell now names the reconstructed deck.
- `run_gate_batches.md`: status update at top; DR-105 struck from §6.
- Morning session owes: read `night_SUMMARY.txt`, judge phase by phase, amend DR-107 / DR-100 /
  DR-105 rows with the verdicts. **None of tonight's runs is a strike by itself.**
