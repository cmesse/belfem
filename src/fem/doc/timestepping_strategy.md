# Timestepping Strategy {#fem_timestepping_strategy}

**Date:** 2026-08-09
**Purpose:** How BELFEM chooses the timestep: time integration schemes, the PID step-size
controller and its legacy fallback, failure and floor policy, warm-start behavior, and the
literature behind the design.
**Module:** `src/fem` (`kernel/cl_FEM_Controller.{hpp,cpp}`, `iwg/cl_IWG_Timestep.{hpp,cpp}`)

Companion document: `src/fem/kernel/doc/nonlinear_controller_theory.md` covers the *inner*
loop (Picard/Newton hybrid, relaxation adaptation, stagnation guards, watchdogs). This
document covers the *outer* loop: how Δt grows, shrinks, is bounded, and is restored across
restarts.

---

## 1. Two Nested Controllers

Every timestep runs two feedback loops:

1. **Inner (per iteration):** the nonlinear solver adapts the relaxation ω after each
   iterate (Messe et al. 2023, Eq. 14) and switches between Picard and Newton
   (`cl_FEM_Controller.cpp`, the loop inside `solve_coupled`). Its outcome — the iteration count for an
   accepted step, or failure — is the *cost signal* for the outer loop.
2. **Outer (per step):** `adjust_timestep()` sizes the next Δt from that cost signal;
   `reset_timestep()` handles rejected steps. This document describes the outer loop.

The inner loop treats small residual changes as noise: when an iterate stays within the
**neutral band** (`mOmegaNoiseBand`, 5% ≈ 0.2 dB) of its predecessor, the controller holds
ω and the divergence counter instead of decaying them. Without that band, "improvement" at
a residual noise floor is a coin flip, and the strict-decrease rule collapses ω
geometrically to its minimum, freezing the iteration. The stagnation band (`stall
tolerance`, default 0.2 dB) is set to the same noise scale so the Newton→Picard demotion can
fire on a flat residual.

## 2. Time Integration Schemes

The default and validated baseline is **BDF1** (implicit Euler; Messe et al. 2023,
Eq. 9). BDF2–BDF5 are available with variable-step coefficients and a startup ramp
(`cl_IWG_Timestep.cpp`, `compute_bdf_coefficients`): during the first steps of a BDF-p run,
the scheme uses the highest order supported by the populated history. A warm restart from a
memdump that carries the BDF state (written since 2026-08-15) resumes at the earned order.
Since 2026-08-30 that state is mandatory: an older dump without it is refused by name rather
than re-anchoring the ramp from BDF1, because every dump the current writer produces carries it.

Higher BDF orders interact with step-size *changes*: the admissible step-ratio bounds shrink
with order (Grigorieff 1983; Hairer & Wanner). The controller therefore clamps growth by the
**active** order (`adjust_timestep`):

| active order | max growth per step |
|---|---|
| 1–2 | 1.5 |
| 3 | 1.4 |
| ≥ 4 | 1.2 |

Shrinking (×0.5) is always benign. Guidance: HTS flux-front transients have limited temporal
regularity, so orders above 2 rarely justify their tighter ratio bounds and higher
perturbation sensitivity. Prefer `bdf1` or `bdf2`; treat `bdf5` as experimental. BDF5 is
also only A(α)-stable (α ≈ 51.8°).

## 3. Step-Size Control: PID (default) and Legacy

### 3.1 The signal

Both modes control on the **iteration count** of the last accepted step, expressed as the
cost error

```
e = N_iterations / N_target        ( "target iterations" in the deck )
```

`e > 1` means the step was too expensive (shrink); `e < 1` means it was too cheap (grow).

### 3.2 Legacy mode

The original rule (`mUseLegacyTimestepControl = true`, internal switch in
`cl_FEM_Controller.hpp`) is

```
φ = sqrt( N_target / N ) = e^(−1/2)
```

In control-theoretic terms, this is a **memoryless I-controller with gain 0.5**. That
controller class is known to overreact to noisy, quantized signals and oscillate between
ceiling growth and rejection cascades (Söderlind 2002). It uses only the last step's count
and forgets rejections immediately.

### 3.3 PID mode (default)

The default mode filters the same signal through a multiplicative (log-space) PID
controller (Valli, Carey & Coutinho 2002):

```
φ = (e₁/e₀)^kP · (1/e₀)^kI · (e₁²/(e₀·e₂))^kD
```

where `e₀`, `e₁`, and `e₂` are the cost errors of the last three accepted steps
(initialized to 1). Defaults (`mCtrlKp, mCtrlKi, mCtrlKd` in `cl_FEM_Controller.hpp`):

| gain | value | rationale |
|---|---|---|
| kP | 0.15 | damps reaction to single-step noise |
| kI | 0.30 | steady pull toward the iteration target; legacy ≡ kI = 0.5 alone |
| kD | 0.0 | derivative action amplifies quantization noise — off by default |

The legacy rule is exactly the (kP, kI, kD) = (0, 0.5, 0) point of this controller, so an A/B
comparison is a gain setting, not a separate formula (note: PID mode also applies the
post-rejection hold below; zero `mPostFailureHoldSteps` for a pure-formula A/B).

After the PID or legacy ratio is computed, both modes share the same order-aware growth
clamp (§2), the [0.5, φ_max] ratio window, the `minimum`/`maximum timestep` clamp, and
save-grid trimming.

**Rejection memory (PID mode):** after a timestep cut, the error history resets to 1 (a
failed attempt's cost says nothing about the next accepted step), and growth is vetoed
(φ ≤ 1) for the next `mPostFailureHoldSteps = 2` passes through `adjust_timestep` — standard
practice in production BDF codes (Brenan, Campbell & Petzold 1996, Ch. 5). Because the first
accepted step after a cut skips `adjust_timestep` entirely (its predecessor's iteration
counter is zero), growth effectively resumes on the fourth accepted step after a cut. This
prevents the controller from repeatedly re-climbing the same cliff and crashing back to the
floor.

## 4. Failure Policy and the Timestep Floor

A rejected timestep is an **expected algorithmic failure**, not an error: it feeds the retry
policy (cut ×0.5, restore the savepoint, re-run) rather than aborting — see the error-tier
discussion in `doc/coding_philosophy.md`. Two terminal conditions bound the retries:

1. **Persistently singular systems:** a soft solver failure (e.g. MUMPS INFOG(1) = −10) cuts
   the timestep; eight consecutive such cuts abort (`mSolverFailCount`).
2. **The timestep floor:** a cut clamps at `minimum timestep` and never halves below it
   (before 2026-08-09, cuts walked Δt below the configured minimum indefinitely). A cut
   demanded while Δt already sits **at** the floor cannot change the timestep — instead of
   aborting on the spot, the controller **widens the iteration budgets**: each floor retry
   doubles `max iterations` and `watchdog window` (both fields), up to
   `mFloorEscalationCap = 4` times the deck values, prints a warning box, and retries at the
   minimum. A slowly converging attempt gets the iterations it needs; the budgets return to
   the deck values with the next accepted step.
3. **The floor backstop:** escalation is not unbounded. After `floor retries` consecutive
   attempts at the minimum (default 20, `0` disables the backstop) the run stops with a
   diagnosis. The rationale is that a residual which responds to neither the timestep nor a
   4× iteration budget will not converge by repeating the attempt, and an unattended job
   would otherwise spend its entire allocation at the floor. Everything up to the last
   accepted step is already on disk, so this is a bounded loss.

A residual floor independent of Δt (solver noise floor, κ(A)·ε_mach, a model inconsistency)
is exactly the case that walks into the backstop: more iterations at a smaller step cannot
remove it. Diagnose it, then either fix the cause or accept it — see below.

The sanctioned way to *accept* a run whose residual sits at a known noise floor above the
relative tolerance is the **absolute tolerance** escape: `absolute tolerance` in the
`nonlinear` (magnetic) and `nonlinear thermal` sections terminates the loop on the absolute
residual ‖Ax−b‖. It is dimensional and disabled (0.0) by default.

To *measure* whether a floor is κ(A)·ε_mach, enable
`solver { compute conditioning : true ; }`. The timestep footer then shows each field's
spectral ratio, `|λ|max/|λ|min`, on its own row. You can also set the key per field in
`linear magnetic { }` / `linear thermal { }`; a per-field setting overrides the shorthand.

Since 2026-08-30 this key controls the **eigenvalue estimate only**. The MUMPS error-analysis
rows are a separate, independent diagnostic behind their own key:

```
solver
{
    compute conditioning : true ;    // |λ|max/|λ|min rows  — ARPACK, one extra solve/step/field
    mumps error analysis : true ;    // MUMPS ADD COND1/COND2 rows — ICNTL(11), no extra solve
}
```

The two report **different quantities and must not be compared**. `compute conditioning`
computes the matrix's spectral ratio alone. That ratio equals κ₂ only for a normal matrix;
the h-φ Jacobian is not normal, hence the label. `mumps error analysis` reports the
Arioli–Demmel–Duff COND1/COND2 pair: a componentwise 1-norm condition estimate for the solved
*system with its actual right-hand side*. On the same tapestack3d thermal operator, the
ratio was 2.96e7 while COND1 was in the 1e4–1e5 range. Both are correct; they answer
different questions.

Two separate things make the numbers differ. It is worth being precise about both: a gap of
several orders of magnitude can look like a bug until its causes are clear.

**They are different quantities.** COND1 is a componentwise estimate that includes the
right-hand side: it answers "how much forward error can this particular solve carry?" The
spectral ratio is a property of the matrix alone and ignores the right-hand side entirely; it
answers "how hard is this problem?" A favorable right-hand side on a badly scaled matrix can
give a small COND1 and a large ratio, and neither is wrong.

**Under Newton they are also computed from different matrices.** The eigenvalue path takes
the assembled system matrix, while a Newton body solves the tangent matrix — the same system
matrix plus the `dJdx` term. So the ratio measures the difficulty of the *problem*, and COND1
measures the conditioning of the *solve that actually ran*. This is by design, not a
discrepancy to reconcile: the ratio is the more useful of the two for deciding whether a
residual floor is conditioning-limited, precisely because it does not move when the tangent
does. Under Picard the two operators coincide.

With both keys on for a MUMPS field, the footer's closing block reads:

```
   │                    Time for eigenvalue analysis    :           412 ms │
   │                    |λ|max/|λ|min, Magnetic         :        1.43e8   │
   │                    MUMPS ADD COND1, Magnetic       :        2.10e4   │
   │                    |λ|max/|λ|min, Thermal          :        1.21e6   │
   │                    MUMPS ADD COND1, Thermal        :        3.40e3   │
   │                    Time for timestep iteration     :      0 min 23 s │
   ├───────────────────────────────────────────────────────────────────────┤
   │                    2026-08-30  21:41:07                               │
   └───────────────────────────────────────────────────────────────────────┘
```

Rows appear only for the diagnostics a field asked for; a `MUMPS ADD COND2` row joins its
COND1 when its term can contribute; it is omitted, not printed as `n/a`, when it cannot —
usually because no matrix row fell into its category, and also in the rare case where rows
did but omega2 is exactly zero.

Both diagnostics sample **once per timestep, at the first iterate**. For the ADD pair, the
controller arms MUMPS's error analysis for that single solve and disarms it immediately.
Before 2026-08-10, leaving it armed ran the analysis for every solve of every iteration,
even though only one value per step was read. The first iterate is also the better sample because it
is assembled at the converged previous step. This lets the series be read as a trend. Nothing in the
iteration scheme consumes either number, which is why one sample per step suffices.

The eigen estimate costs one extra solve per timestep per field, timed in the footer. Once
the ratio exceeds about `tolerance / ε_mach`, it cannot resolve the small end of the
spectrum. A mixed h-φ system reaches ~1e17; there the row prints `n/a` and a one-line
banner explains why. The ADD pair is inexpensive but available only with MUMPS. Elsewhere, the key is a no-op
and the parser warns.

## 5. Warm Starts

Restarting is **opt-out** (since 2026-08-10): by default a rerun resumes from an existing
`memdump.hdf5`, restoring fields, time, and Δt. It announces the resume point with a WARM RESTART banner
(this resume used to happen silently). A deck sets
`timestep { restart : false ; }` to ignore the dump (a log line says so) and start fresh;
the dump is then overwritten at the first save.

The memdump (`save_memdump` / `load_memdump`) persists the **current Δt** (`delta_time` in
the `meta` group) alongside the fields, but **since 2026-08-19 the loader does not adopt it
verbatim: the first step after a warm restart is capped at the deck's `initial timestep`**
Older dumps without the key fall back to the deck value; the restored Δt is clamped
to the `minimum`/`maximum timestep` window, then capped.

**Why a restart may not re-enter at the earned step.** In the measured cases, the first
post-restart Δt determined whether the process entered the pathological solver path, and
reducing Δt afterwards did not recover it. Measured on tapestack3d from one dump: re-entering at the earned 50 ms gives a first preconditioned residual of 1.1e4 with
every linear solve hitting its iteration cap, the step grinds and is rejected. Shrinking
to 25 ms *within the same process* leaves it at 9.9e3. Re-entering the same dump at 5 ms
gives 0.903 in 16 Krylov iterations, converges in two Picard iterates, and that process then
climbs back to 50 ms and stays healthy. A *live* run at the identical state and Δt is also
healthy, so the earned step is not itself too large. What differs is the step the process
started from. Working hypothesis, inferred from measurements and source tracing but not
independently proven: solver-entry state derived from the first Jacobian is reused later in
the process, and because BDF assembles `J = αM + ΔtK`, a large first Δt may put that reused
state in a worse regime.

The cost is the growth back: `initial × 1.2ⁿ` under the order-clamped limiter, ~13 steps on
tapestack3d (5 → 50 ms), and longer on decks whose `initial timestep` is far below their
working step. Those steps are cheap against a rejected grind, but they are not free.

**This deliberately reverses only the Δt half of the original restart-cliff policy.** The Δt half is
reverted; the BDF-history half (`bdf_h`, `bdf_step_count`, `bdf_last_dt`) is kept, and it is
what actually prevents the order cliff, and the two are separable. `bdf_last_dt` is
*not* usable as the cap: every dump is written at a save point, so the last completed step is
a save-point value, not necessarily a safe first step for a new solver process. The BDF history *spacing* is persisted with it (`bdf_h`, `bdf_step_count`,
`bdf_last_dt` — the last being the just-completed step the first post-restart shift pushes
into the history; it is distinct from the dumped next-step `delta_time`). A resumed run
therefore continues at its earned order with the exact variable-step coefficients of the
uninterrupted trajectory. Dumps without these keys are refused (since 2026-08-30): they can
only have been written by a binary predating 2026-08-15, and silently rebuilding the order from
BDF1 was the restart cliff this feature exists to remove. Use `restart : false` to ignore such
a dump and start fresh.

## 6. Deck Keys (quick reference)

See `doc/input_file_reference.md` §4.2–§4.4 for the authoritative list. The outer loop reads
`initial/minimum/maximum timestep`, `simulation time`, `adapt timestep`, `scheme`, `save
every`, and `target iterations`. The PID gains and the legacy switch remain internal members
(no deck keys) until the mode is validated.

## 7. Literature

Not held in `./literature/`; DOIs for independent access:

- Gustafsson, Lundh & Söderlind 1988, *A PI stepsize control for the numerical solution of
  ordinary differential equations*, BIT 28:270–287. doi:10.1007/BF01934091
- Gustafsson 1994, *Control-theoretic techniques for stepsize selection in implicit Runge-Kutta
  methods*, ACM TOMS 20(4):496–517. doi:10.1145/198429.198437
- Söderlind 2002, *Automatic control and adaptive time-stepping*, Numerical Algorithms
  31:281–310. doi:10.1023/A:1021160023092
- Söderlind 2003, *Digital filters in adaptive time-stepping*, ACM TOMS 29(1):1–26.
  doi:10.1145/641876.641877
- Söderlind & Wang 2006, *Adaptive time-stepping and computational stability*, J. Comput.
  Appl. Math. 185:225–243. doi:10.1016/j.cam.2005.03.008
- Valli, Carey & Coutinho 2002, *Control strategies for timestep selection in simulation of
  coupled viscous flow and heat transfer*, Commun. Numer. Meth. Engng 18(2):131–139.
  doi:10.1002/cnm.475
- Brenan, Campbell & Petzold 1996, *Numerical Solution of Initial-Value Problems in
  Differential-Algebraic Equations*, SIAM. doi:10.1137/1.9781611971224
- Messe et al. 2023, Section 4 — the BELFEM baseline: BDF1, Eq. 14 relaxation rule,
  halve-on-failure, grow-when-cheap.
