# DR-127: Certified Exit for the Nonlinear Loops

**Date:** 2026-08-29
**Purpose:** Land the DR-127 fix — the nonlinear loops now exit only on a state whose residual
was actually measured — via the full plan+audit→code+audit ceremony, overnight per Christian's
approval.
**Module:** fem/kernel (Controller, DofManager, SolverData), executables

## The defect (recap)

Always-accept Picard commits x_{k+1} while the printed ε certifies x_k: the residual is
evaluated on the entry state of the trip, the update is committed unconditionally, and the loop
exits on a tolerance test of the *previous* state. A 2-iterate exit certifies x₁ and
time-steps an unchecked x₂. The lagged post-update residual is provably (1−ω)·r — vacuous.
Filed 2026-08-28 by the jjc-noise jury round (see `todo/tapestack3d_jjc_noise_2125ms.md`);
design decision in the source, not a bug, but its contractivity assumption fails on
near-null systems (demonstrated on the then-singular DR-126 deck).

## The fix (plan v4, Fork A: fold the driver loop into the Controller)

Three audit rounds on the plan (v1/v2/v3 each rejected with substance, all integrated); frozen
spec `tmp/ai_exchange/dr127_certified_exit_plan_v4.md` = v3 + N1–N16 normative + N17–N19
deconfliction with the parallel eigen/ARPACK session.

**Stage A — two-phase solve seam (SolverData/DofManager).**
`solve()` intercepts iterative single-RHS problems and splits into:

- `compute_residual()` — loads, Dirichlet multiply, ‖b‖, r = Ax − b, one-shot per assembly
  (`mResidualReady` guard, invalidated by `compute_jacobian_and_rhs`), no barrier. Newton
  saves b into the new `mRhsOriginal` before overwriting with r; Picard NaNs
  `mPreUpdateResidual`. Retained scratch (`mFixedValuesScratch`, `mCollectedFields`) — no
  per-call allocation.
- `solve_from_residual()` — all-ranks barrier pairing, per-algorithm master solve + update,
  Newton restores b for the post-update recompute, Picard `mRhsBackup` = r. The old
  monolithic Iterative case is excised to a `BELFEM_ERROR` stub.

**Stage B — Controller.** `iterate_coupled()` restructured to head/body/tail:

- HEAD: fresh assembly + `compute_residual()` of the *committed* state → certificate.
  Body solves run only while uncertified or under `min iterations`; the first body trial
  consumes the head assembly (`tFirstTrial`), so a solving trip costs exactly what it does
  today. Thermal head sits at the Gauss-Seidel position and doubles as the co-certificate
  when it lands on a no-mag-body trip — an exit costs two assemblies only when both fields
  are active.
- CERTIFIED EXIT: the trip returns (`mTripExit`) only when the magnetic head is certified
  and the thermal side is certified, absent, or flat-stall latched (the stall accept prints
  an explicit "ACCEPTED AT STALL" uncertified label).
- In-loop machinery (handoff, line search, AIMD, watchdogs, history, prints) keeps today's
  residual variables and semantics untouched — the head certificate is an *additional*
  measurement that gates only the exit. AIMD pairs against the previous trip's residual
  (`tEpsilonPrev` capture), not the head — for Picard the body result equals the head, and
  pairing against it would freeze the omega adaptation.
- Tail bookkeeping gated on `tSolveM` (strikes, stagnation, watchdog, `++mIteration`);
  thermal budget/reset checks deliberately ungated. Dead `mEpsilonFirst` family removed.
- New public wrappers `solve_coupled()/solve_magnetic()/solve_thermal()` own the loop;
  `iterate_*` is now private and the `run_*` predicates were deleted outright (see the
  audit section). Segregated loops got the same head + certified early-return.
  `solve_thermal()` returns on `mResetThermal` so the driver's outer time2 loop re-inits
  at the shrunken sub-step — the contract the code audit had to correct, see below.

**Stage C — drivers.** All 7 `while(run_*()) iterate_*();` loops in `belfem.cpp`,
`hphirun.cpp`, `hphiTrun.cpp` collapse to `tControl->solve_*();`; the reset()-retry
structure is preserved.

**Stage D — records.** `nonlinear_controller_theory.md` (certified-exit section),
`maxwell_usage_guide.md`, `maxwell/doc/README.md`, `circuit_usage_guide.md` snippets;
input contract: `doc/input_file_reference.md` §4.2 + `doc/input_schema.yaml` — since this
change `min iterations` counts BODY solves; `min iterations : 0` permits a zero-solve
predictor accept; the thermal min now also binds on the coupled path (declared, D6).

## Semantics changes declared

1. The exit residual is now always the committed state's, measured from a fresh assembly.
2. `min iterations` counts body solves; 0 allows a zero-solve exit.
3. Thermal `min iterations` joins the coupled path.
4. A predictor that already meets tolerance can exit without solving (when min = 0).

## Evidence status

**Verified.** Every touched TU passed a GCC syntax gate under the real build
flags (`compile_commands.json`); no stray callers of the privatized API in `src/` or `tests/`.
The code-audit jury round on the full diff is recorded in
`tmp/ai_exchange/review_dr127_certified_exit.md` (pre-registration with 4 declared plan
deviations + 8 self-flagged risks). Christian rebuilt and ran the test suite on the morning
of 2026-08-29: **all tests pass** — the build and suite gates are green. Still owed as
behavioral gates: the tapestack3d rerun (accepted-state residual window at t ≈ 2.1 s must
stay smooth) and a segregated hphirun/hphiTrun deck smoke exercising the corrected
`solve_thermal()` retry path.

## Code-audit round (same night)

Jury on the full working-tree diff, blind, both vendors: Codex "reject pending
correction", Grok "do not land as-is" — **independently converging on the same two
defects**, both confirmed against the tree and fixed the same night:

1. **P0 — `solve_thermal()` dropped `mResetThermal`**: `reset_thermal()` sets only that
   flag; the smaller sub-step exists only after the driver's next `initialize_thermal()`
   (which recomputes `mDeltaTime2` and clears the flag). My wrapper looped on
   `mReset || mTripExit` and would have retried the failed sub-step at the stale Δt₂
   forever. Fix: return on `mResetThermal` too — the driver's outer time2 loop re-inits,
   restoring the old `run_thermal()` contract exactly.
2. **P1 — `tMagneticDoomed` ungated on skip trips**: a certified magnetic head at
   `mIteration == max iterations` froze the thermal side of a healthy step and tripped
   the deadlock `BELFEM_ERROR`. Fix: the same `tSolveM` gate that round-3 R6 gave
   `tMagneticReset`.

Also applied from the round: dead `Vector` allocation in `SolverData::residual()`
removed (was 2× per trip on the hot path); `iterate_magnetic` ω clamp-before-set order
aligned with the other paths; the uncertified `run_coupled/run_magnetic/run_thermal`
API **deleted** (not just privatized — zero callers, and a future internal call would
have silently reintroduced the DR-127 bug); write-only `mEpsilonFirst*` family removed;
stale `file:line` cites and one doc drift fixed. Grok's flag that the coupled thermal
`min iterations` is a live behavior change is correct and *declared* (plan v4 D6, both
input-contract artifacts). All TUs re-passed the syntax gates after the fixes.
Full reconciliation table: `tmp/ai_exchange/review_dr127_certified_exit.md`.

## Run gate: tapestack3d rerun (PASSED, 2026-08-29)

The deck that filed the bug — `cmake-build-debug/tapestack3d`, sp-ap table with
`resistivity type : piecewise`, the same input.conf as the defect run — rerun on the
fixed build, 137 timesteps to t = 3.3 s.

**The 2.125 s frame is no longer special.** Across the whole 2.00–2.40 s window,
max|B| and max j/jc grow by exactly ×1.047 per frame, every frame, tracking the
transport-current ramp:

| t [s] | max\|B\| [T] | max j/jc | step ratio |
|---|---|---|---|
| 2.100 | 3.2181e-03 | 4.5098e-02 | ×1.047 |
| **2.125** | **3.3690e-03** | **4.7212e-02** | **×1.047** |
| 2.150 | 3.5269e-03 | 4.9426e-02 | ×1.047 |
| 2.175 | 3.6922e-03 | 5.1742e-02 | ×1.047 |

Defect signature for comparison: max|B| ×3 **in one step** at this frame, tape-1
transverse current 5.8×/7.0× per frame, I₁ ×40. Now the largest frame-to-frame I₁
ratio anywhere in the run is ×1.115 (at t = 2.25 s), and I₁ grows ×8.0 smoothly over
1.17 s while I₀ ramps 9.5 → 17 A. The scatter measure `std(j/jc)` over active nodes
tracks the same ×1.047, so the growth is the flux front, not noise; the active-node
count climbs monotonically 19306 → 19486.

**Controller behavior on the certified exit:** every step takes exactly 2 Picard
bodies and exits on a head certificate at 2.6e-15…5.0e-15 (the printed `certified
exit: magnetic …` line is the committed state, which is the entire point of DR-127).
Zero timestep cuts, zero resets, zero divergence strikes, zero stalls across 137
steps. Conditioning flat at 1.26e9 (it was 1e17–1e21 before DR-126).

**What this run does and does not prove.** It is a fully-coupled run on a
well-conditioned system, so it verifies that the certified exit is correct, cheap and
non-disruptive — it does *not* independently prove DR-127 would rescue a near-null
system, because DR-126 removed the singularity that made the always-accept
contractivity assumption fail. The two fixes cannot be attributed apart from this
single run, and that is fine: the deck is smooth and every exit is now measured.
`solve_thermal()` is untouched by a fully-coupled deck, so the audit-corrected
`mResetThermal` return still owes a segregated smoke that reaches a thermal sub-step
cut.

**Segregated and warm-restart paths exercised (2026-08-29, `out.txt` / `out2.txt`).** The
deck was rerun with `coupling : segregated`: 69 steps, each taking 2 magnetic bodies to a
certified exit and 2 thermal bodies to a certified sub-step exit, with zero resets or
warnings — so `solve_magnetic()`, `solve_thermal()`, both segregated heads and both
certified early-returns all execute. A subsequent warm restart from `memdump.hdf5` resumed
at step 70 with the BDF history restored (Δt reopening at 5 ms and ramping back to 25 ms)
and ran clean, so the certified exit composes with `load_memdump`. Neither run hit a
thermal sub-step cut on its own, so the audit-corrected `mResetThermal` return in
`solve_thermal()` was forced with a temporary probe (Christian's call: no new DR row for one
unexecuted branch). The probe discarded the first successful thermal sub-step and called
`reset_thermal()` in its place, once. Result (`out3.txt`, warm restart at step 99): the probe
line printed, `solve_thermal()` returned to the driver, `initialize_thermal()` re-ran with the
doubled coupling factor, and the step completed as **two** sub-steps at the halved Δt₂
(9.307e-11 → 1.650e-11, 1.176e-11) before advancing to step 100. **No hang** — which was the
failure mode the original wrapper would have produced. Probe removed the same session; `grep`
confirms zero residue. This was DR-127's last gate.

Reading note for these logs: a printed `magnetic : 2.220e-16` (equivalently −156.54 dB) is
the `BELFEM_EPS` clamp in `SolverData::residual()` (`cl_FEM_DofMgr_SolverData.cpp:2577-2582`),
applied when the residual norm underflows machine epsilon — a floor, not a measurement. The
fully-coupled run's 3.9e-15 / 5.0e-15 values are genuine.

**DR-128 is live in this very run.** max|B| stays between 2.7 mT and 5.6 mT across the
window, entirely below the sp-ap table's B-axis floor at 10 mT — so jc and n are
B-independent (clamped) for the whole deck. The noise is gone with the clamp still in
place; DR-128 remains open on its own merits.

## Files

- `src/fem/kernel/cl_FEM_DofMgr_SolverData.{hpp,cpp}` — two-phase seam
- `src/fem/kernel/cl_FEM_DofManager.{hpp,cpp}` — forwarders + invalidation
- `src/fem/kernel/cl_FEM_Controller.{hpp,cpp}` — head/body/tail + solve_* wrappers
- `src/executables/belfem.cpp`, `hphirun.cpp`, `hphiTrun.cpp` — driver loops
- `src/fem/kernel/doc/nonlinear_controller_theory.md`, `src/fem/maxwell/doc/*`,
  `src/circuit/doc/circuit_usage_guide.md` — docs
- `doc/input_file_reference.md`, `doc/input_schema.yaml` — input contract
- `todo/debt_register.md` — DR-127 row
