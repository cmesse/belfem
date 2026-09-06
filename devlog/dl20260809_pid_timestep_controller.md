# PID Timestep Controller + Δt-Collapse P1 Fixes (sidecoatings)

**Date:** 2026-08-09
**Topic:** Root-causing the sidecoatings Δt collapse (three-AI jury), fixing the round-1 P1
findings, and implementing the jury-audited PID step-size controller with legacy fallback.
**Participants:** Christian (rulings, deck A/B, concurrent C7 fix), Claude (analysis, plan,
implementation, verification), Codex + Grok (three blind jury rounds).
**Threads:** `tmp/ai_exchange/review_fem_controller_timestep.md` (round 1),
`review_pid_timestep_draft.md` (round 2), `review_pid_timestep_impl.md` (round 3).
**Plan:** `todo/pid_timestep_controller_plan.md`.
**State:** implemented, NOT compiled (build is Christian's); round-3 jury verdict appended
below when reconciled.

## The failure

`cmake-build-debug/sidecoatings` (3D thin-shell tapestack, edge coating on, fully coupled,
originally bdf5): at timestep 98 the printed Δt is 0.0000 ms against a 100 ns deck minimum.
Log signature: magnetic Picard contracts 5.0e-5 → 1.1e-5 in two iterations, promotes to Newton
at ε_switch = 1e-4, the residual rises to ~2.2e-5 (−46.5 dB) and stays pinned ±0.2 dB for 30
iterations while relax collapses 0.5 → 1.3e-4; thermal is converged at −100 dB throughout;
the watchdog cuts at iteration 32 = 2 + window. The floor is Δt-INDEPENDENT — identical at
every halving — so cutting can never succeed.

## Root-cause chain (round 1, all verified)

1. **C1 (P1, Claude+Codex blind agreement):** `reset_timestep()` halves Δt with no
   `mDeltaTimeMin` clamp and no abort; the theory doc promised the bound. Hence the grind
   below the deck floor.
2. **C2 (P1, routed to Christian):** Newton solves with the consistent tangent J but is judged
   by the lagged-operator residual ε = ‖A(x_old)x_new − b‖/‖b‖ — a hybrid Messe 2023 Eq. 13
   does not cover (the paper's quasi-Newton uses A itself). The Newton entry visibly RAISES
   the lagged ε and never recovers. Same family as the Anderson-ε false-convergence issue.
3. **C3 (P1):** the Eq.-14 AIMD requires strict improvement every iterate (grow ×1.1) and
   halves ω otherwise; at a noise-floor residual improvement is a coin flip → geometric ω
   collapse → frozen state → self-confirming stall. The 10-consecutive-divergence ω-restart
   can never fire at a coin-flip floor.
4. **C4 (P1):** the stagnation flat-band (0.001 dB) sat two orders below the floor's noise, so
   the designed Newton→Picard demotion was dead; the only exit was the watchdog cut.
5. **S4 (root question, Christian):** the −46 dB floor itself. Candidates: κ(A)·ε_mach
   (~1e11 · 1e-16 ≈ 1e-5 — fits), the edge-coating wiring, C2's metric mismatch. Christian's
   bdf1 A/B (this evening) discriminates the scheme hypothesis; the MUMPS conditioning probe
   (`mComputeConditioning`, currently not deck-enableable) would give κ directly.

Also confirmed in round 1 (backlog, not this session's fixes): Grok's rank-0-only circuit BC
payload (P1 for MPI+circuit runs, needs a probe), segregated-path state leaks
(`mNumIterationsDiv` shared magnetic/thermal, `initialize_magnetic` hygiene, unclamped first
`set_omega`), circuit NaN-residual accepted as success, terminal-pair indexing disagreement
between `initialize_timestep` and `compute_circuit_current`, Timer Rule-of-Five. One Grok
claim (log10 of non-positive residual) REFUTED — `residual()` guards zero. C7
(`compute_hanging_dofs` unreachable) was independently root-caused in the parallel viz session
and its fix (call hoisted out of the algorithm switch) rides in the same working tree.

## What was implemented (Christian's execute order)

**P1 package** (`cl_FEM_Controller.{hpp,cpp}`):
- C1: cut path clamps to `minimum timestep` (one attempt runs exactly AT the floor) and a
  further cut there aborts with a diagnostic pointing at the residual floor, the tolerance,
  and the absolute-tolerance escape. Applies to both controller modes.
- C3: neutral band `mOmegaNoiseBand = 0.05` (≈0.2 dB) in all four AIMD tails — a wobble within
  the band holds ω and the divergence counter; genuine regressions still decay.
- C4: `mStallBand` default 0.001 → 0.2 dB, matched to the noise scale so the Newton→Picard
  demotion can fire (deck key `stall tolerance` still overrides).
- C5: magnetic `absolute tolerance` wired live (was parsed-but-dead; Claude+Grok blind
  agreement): `mEpsilonAbs` captured per solve, consumed in `run_coupled`/`run_magnetic`,
  inert at the 0.0 default. This is the sanctioned accept-at-floor escape.
- C6: the progress watchdog re-anchors on the normal Picard→Newton promotion (the escalation
  path already did), so Newton is judged against its own entry residual, not Picard's best.
- Setup validation (round-2 finding): `target iterations` > 0; `initial timestep` present,
  positive, inside the min/max window; min > 0; max ≥ min.

**PID controller** (plan §3, post-audit): log-space multiplicative PID on the iteration-count
error (Valli, Carey & Coutinho 2002), gains (kP, kI, kD) = (0.15, 0.30, 0), legacy rule
≡ pure-I kI = 0.5 behind internal `mUseLegacyTimestepControl`; post-cut growth hold (2
adjusts, SUNDIALS practice) + error-history reset; order-aware growth clamps unchanged.
Memdump now persists `delta_time` in the meta group (lvalue-`hid_t` pattern, read inside the
open group, all-ranks broadcast + clamp) — warm starts resume at the earned Δt instead of the
deck initial (the observed restart cliff). Docs: new `src/fem/doc/timestepping_strategy.md`
(+ README), `doc/input_file_reference.md` synced (absolute tolerance live, stall band default,
floor semantics, initial-timestep validation).

## Jury round 2 (draft audit) highlights

Codex+Grok caught, pre-implementation: a real compile blocker in the memdump sketch
(`active_group()` returns `hid_t` by value, helpers take `hid_t&`), the consecutive-cuts abort
contradicting the coding-philosophy retry tier (feature dropped — floor abort is the sole new
terminal condition), the hold under-count through the `mIteration0` gate (documented), the
degenerate e₀ = 0 sample (floored), and the honest R5 criterion: **fewer cascades or a clean
abort — not "sidecoatings reaches tolerance"** (the floor is a separate disease). Full defect
log in the plan §6.

## Jury round 3 (implementation diff) — closed same evening

Codex + Grok on the controller-only diff; six correctness-tier findings, all verified and
FIXED the same turn: the ε/ε_abs pair mismatch after exhausted backtracking (Codex — the
pre-registered B1 risk, real), the thermal flat-stall gate ignoring the new magnetic absolute
criterion (Codex), the thermal promotion not re-anchoring the thermal watchdog (Grok), memdump
`delta_time` not NaN-hardened (Grok), negative `absolute tolerance` silently disabling the
escape (Grok), and the zero-iteration divide surviving on the legacy branch (Grok — ASSERT
hoisted over both modes). Grok's two policy flags stand as R5 watch items, not code changes:
the PID default and the 0.2 dB stall band change behaviour for every deck without a deck key —
Christian's execute-order ruling, with the internal legacy switch, the (0, 0.5, 0) gain point,
and the `stall tolerance` key as rollback levers. Both auditors independently verified the
warm-start broadcast shape, the floor policy table, the neutral-band structure, and the error
tiers as clean.

## Addendum 2026-08-10 — batch 2 (Christian's revisions, implemented on his order)

Christian committed batch 1 as `1d6ef305` and ruled three revisions plus one cosmetic change,
all implemented the same day (uncommitted, not compiled; jury round 4 on the incremental diff,
thread `review_pid_timestep_impl_b2`):

1. **Floor policy revised — escalate first, backstop second.** The 2026-08-09 floor abort is
   replaced by bounded iteration-budget escalation: each cut demanded at `minimum timestep`
   doubles `max iterations` and `watchdog window` (both fields) up to 4× the deck values,
   prints a warning box, and retries at the floor; `finalize` restores the deck values on the
   next accepted step. Both round-4 auditors then pushed back on the unbounded retry as their
   only P1 (an unattended job spins at the floor forever, each attempt 4× the cost), and
   Christian took their middle ground: a new deck key `timestep { floor retries : n ; }`,
   default 20 and `0` = unlimited, stops the run with a diagnosis after that many consecutive
   floor attempts. Results up to the last accepted step are already on disk, so the loss is
   bounded; `absolute tolerance` remains the accept path for a genuine noise floor.
2. **Warm restart made visible and controllable (opt-OUT, Christian's final ruling).** Root
   cause of the "keeps continuing" surprise: both drivers call `load_memdump("memdump.hdf5")`
   unconditionally (`hphirun.cpp:83`, `hphiTrun.cpp:90`), so ANY rerun silently resumed if a
   dump sat in the directory. New key `timestep { restart : ... ; }`, default TRUE (existing
   crash-resume workflows keep working) — but a resume now prints a WARM RESTART banner with
   resume t / timestep / Δt instead of happening silently; `restart : false` ignores the dump
   with a log notice and starts fresh. (First implemented default-false; Christian flipped it
   to opt-out the same hour.)
3. **`solver { compute conditioning : true ; }`** — the diagnostic member existed but its
   constructor-time MUMPS arming ran before `set_params`, making it unreachable from any deck
   (round-1 finding). This is the κ(A) probe for the −46 dB floor question.
   **Reworked the same day** on Christian's observation that one sample per timestep should
   suffice: the cost was never in the read (`compute_conditioning()` already ran once per step
   from `finalize`) but in the arming — `MumpsErrorAnalysis::Full` was set once and never
   cleared, so MUMPS ran its error analysis on *every* solve (30+ per timestep) while one
   value was read. The controller now arms around the **first iterate only** and disarms as
   soon as it reads, per field (`arm_conditioning_*` / `capture_conditioning_*`, rank-uniform
   since ICNTL must agree across the communicator). First rather than last: last is not
   knowable before the residual test without a re-solve, and the first iterate is assembled at
   the converged previous step, so κ is sampled at a comparable state every timestep and the
   series is trendable. Thermal κ added (`mConditionNumber2`); a field whose solver cannot
   supply κ cheaply prints `n/a`, and no eigenvalue fallback is run for the thermal field.
4. **Controller box header** shows the configured scheme (BDF1…BDF5, Explicit, CN, Galerkin)
   instead of "Timestep"; the footer's closing row now carries κM and κT beside the timing
   (Christian's layout, 22+48 columns — the κ glyph is two bytes in UTF-8, so the cells are
   built from literal box characters and fixed-width fields, never `%-Ns` padding, and the
   condition numbers use a compact `1.43e8` formatter so the column never shifts).

Docs synced same turn. `doc/input_file_reference.md`: new "keys read directly in the `solver`
section" table (`compute conditioning`), new `timestep → restart` row, revised
`minimum timestep` / `absolute tolerance` / `max iterations` / `watchdog window` /
`target iterations` rows, a "not deck-exposed" paragraph naming the internal controller
constants (PID mode and gains, hold, escalation cap, neutral band) so nobody hunts for keys
that do not exist, a "When Δt collapses toward `minimum timestep`" field-guidance block in
§4.5 (the two-cause split and the κ(A) measurement that separates them), a §13 pitfall for
the default-on restart, and a §2 note that `memdump.hdf5` is a fixed name. The §4.2/§4.3/§4.4
`Site` anchors were **re-baselined** — they had drifted by roughly 300 lines across the whole
table (pre-existing, flagged in round 3) and now match `set_params` in the current tree.
`src/fem/doc/timestepping_strategy.md` §4/§5 updated for both policy changes.

## Open / next
- R5 A/B (Christian): legacy vs PID, bdf1 vs bdf5, `edge coating : off` discriminator; the
  conditioning probe needs a one-line enable if κ is wanted.
- C2 residual-semantics ruling (Christian) — mitigated but not resolved by C3/C4.
- Backlog from round 1: circuit MPI payload probe, segregated-path hygiene, terminal-pair
  indexing trace, Timer Rule-of-Five.

## Addendum 2 (2026-08-10 evening) — greg5 diagnosis, line-search fallback, tangent jury

The first κ measurements landed and immediately earned their keep. greg5 (2D thin-shell,
constant jc = 4e10 / n = 25 power law, 160 A ramp, point bearing, bdf1) with Picard: healthy —
κM pinned at 6.4e4 the whole run, 5.7% rejected attempts, Δt small because
`target iterations : 50` is a COST setpoint and ω̄ ≈ 0.2 makes iterations cheap-but-many.
Claude wrongly recommended trying Newton; Christian's out2 run showed the documented
net-current + point-bearing pathology precisely: Picard contracts at ω = 1, promotion at
−44 dB, the line search rejects ALL eight trials to ω = 0.5·2⁻⁷ = 0.00391, the step dies, and
Δt collapses 22 ms → 0.5 µs while κ(A) stays at 6.4e4 (the Picard operator is fine — the
TANGENT is the problem, and κ(A) structurally cannot see it).

**Line-search fallback implemented (Christian's order):** an exhausted backtracking budget on
a PROMOTED Newton tangent now restores the last accepted iterate and latches `mJustPicard` for
the rest of the timestep (half-throttle Picard resume, stall window cleared, watchdog
re-anchored, notice line printed) instead of cutting Δt — a Δt cut cannot fix a tangent whose
direction is unusable at any scale, and promotion re-fired at every step size. The escalated
tangent (`mForceNewton`) keeps the old cut path: there Picard has already broken down. While
ANY Picard latch holds (magnetic mJustPicard, thermal mJustPicard2), the promotion branch is
now skipped entirely — this also fixes a latent flaw in the batch-1 C6 re-anchor, which would
otherwise refresh the watchdog's best-residual clock on every latched iterate and disable the
guard (all four promotion sites guarded). Docs: nonlinear_controller_theory.md §2 (third
latch) + §5 line-search paragraph, input reference §4.5 Newton guidance updated (Picard still
the right setting for this deck class; the fallback makes Newton survivable, not free).

**Newton tangent derivative jury dispatched** (Christian's request, thread
`review_newton_tangent_derivatives`, brief `newton_tangent_brief.md`): matfix's PARTIAL
tangent converges better on greg5 than HEAD's full one — verify every coefficient of the live
chain (4-arg powerlaw dρ/dJ, j-channel dyad, constant-jc field-channel gating, BDF scaling of
dKdx) and adjudicate wrong-coefficient vs correct-but-stiff-tangent. Claude's pre-registered
hand checks: parallel-combination derivative and dyad assembly CORRECT; hypothesis H1 =
coefficients right, failure = (n−1)≈24× stiffening along ĵ meeting the bearing gauge mode
under the lagged-A acceptance metric. Verdict pending at write time.

**Tangent jury verdict (same evening, reconciled):** the round's premise fell — `git show
matfix` proves matfix carries the IDENTICAL active tangent on the greg5 path (dispatch, dyad,
field hook); matfix "worked better" because its moved-baseline rescue is unconditional while
HEAD gates it on a thermal kernel (Garber R9), so magnetic-only runs never rescue. The
unclamped derivative chain is confirmed correct by three independent derivations; the
constant-jc field channels are exactly zero. TWO real mechanisms survive, both to Christian:
(1) CONFIRMED wrong-derivative band — `mRhoMin = 1e-16` floors ρ_PL for all J < 0.874·jc on
greg5's constants while `drho_powerlaw_dJ` differentiates the unfloored law, injecting an
O(n)·ρ spurious rank-1 term at the flux front's subcritical skirt (Claude's pre-registered
"impact ≈ 0" formally RETRACTED — both auditors caught it blind); proposed fix: zero the
derivative when the floor binds. (2) The lagged-A acceptance metric: the ω→0 limit of a trial
is the FRESH entry residual while the reference is the LAGGED one — a gap > 0.3 decade
rejects all eight trials regardless of direction (round-1 C2 again). Discriminating probes
proposed (per-IP J/jc + clamp fraction; smallest-ω trial vs fresh entry residual). Thread
`review_newton_tangent_derivatives`, closed with reconciliation.

## Addendum 3 (2026-08-10 night) — Christian's rulings on the tangent verdict

1. **`mRhoMin = 0`** (`cl_Material.hpp`): the 1e-16 floor is gone, which eliminates the
   confirmed wrong-derivative band identically — with no floor, `drho_powerlaw_dJ`
   differentiates exactly what `rho_powerlaw` evaluates, everywhere. ρ_PL = 0 is IEEE-clean
   in the parallel combination (1/0 = inf → ρ_eff = 0, the physically right superconducting
   limit); the derivative already guards |J| < ε.
2. **Lagged-A metric, first structural fix — the honest line-search reference.** The
   backtracking reference is now the pre-update residual of the CURRENT assembly
   (`SolverData::mPreUpdateResidual`, captured in the Newton branch between the r = Ax−b
   multiply and the solve, broadcast with the other residual scalars, forwarded via
   `DofManager::pre_update_residual()`), replacing the previous iterate's ε which was
   measured under the PREVIOUS assembly. Rationale (verification-tier, from the tangent
   round): the ω→0 limit of any trial IS the pre-update residual, so with the stale
   reference any fresh-vs-lagged assembly gap > 0.3 decade rejected all eight trials flat in
   ω regardless of direction — the greg5 death signature, and the very pattern matfix's
   unconditional moved-baseline rescue papered over. With the honest reference, the search
   can only fail on genuinely explosive directions, which the same-day line-search fallback
   latch then catches; the two changes layer. Picard iterates set the reference to NaN
   (no line search); the controller falls back to mEpsilon0 if the value is not finite.
   Scope guard: the REPORTED per-iterate ε, the AIMD adaptation, promotion thresholds, and
   convergence tests are all untouched — this fixes the acceptance metric only. The full
   C2 decision (what ε should MEAN for Newton iterates) remains open with Christian.

**Batch-3 jury (same night, reconciled):** the core — honest pre-update reference, exhaustion
latch, promotion-suppression guards, opt-out restart, floor escalation, mRhoMin = 0 — verified
clean by both auditors ("ship-quality for the nonlinear-control fixes", Grok). Grok's headline
R1 (κ analysis staying armed through iterate-0 backtracks) REFUTED by structure — iteration 0
is always Picard — but hardened for free (capture moved inside the trial loop). Fixed from the
round: stale-κ carry-over across timesteps (NaN resets at all three initializers — a frozen
thermal step now prints n/a instead of last step's κT), the conditioning_string 7-character
overflow at mantissa-rounding boundaries (Codex+Grok blind agreement; per-precision carry
thresholds), the floor-backstop message (post-increment count, unconditional "4×"), and the
powerlaws floor doxygen (now states the zero default and the raise-it-again contract).
Codex's piecewise note (exact ρ = 0 at J = 0 under the zero floor) confirmed as fact and left
by ruling — same limit the parallel combination reaches, benign in the M-dominated transient;
on record for any future K-only consumer. Third round of DR-02 dump findings (TRUNC silently
overwrites sysdump_*.hdf5 on rerun) appended to Christian's tracker. Thread
`review_controller_metric_batch3`, closed with reconciliation.

**greg5 out3 — reproducer-tier confirmation (2026-08-10 16:36):** same deck as out2 (Newton,
bdf1, target 20), rebuilt with batch 3. The run that died at t = 8.2 s after 629 attempts now
completes the FULL 15 s in 169 attempts — 2 rejections, 1 honest watchdog cut in the early
ramp, Δt riding the 200 ms deck maximum, κM stable at 6.4e4–1.1e5, and the line-search
fallback latch fired ZERO times. Verdict of the discriminating experiment: the lagged-A
acceptance metric was the entire greg5 disease (with the clamp band removed jointly by
mRhoMin = 0) — no gauge-mode explosions remain. The step-9 frame shows the repaired
machinery end-to-end: Newton 4 accepts a trial the stale reference would have rejected
(−47.4 lagged vs −37.5 fresh), demotes, Picard contracts, re-promotes, and Newton finishes to
−70 dB with ω growing 0.39 → 0.94 at Δt = 114 ms. Consequence flagged for later: the §4.5
"use Picard on point-bearing net-current decks" guidance was written against what we now know
was largely a metric artifact — revisit after A/B on more decks, not from one run.

**sidecoatings κ measurement (2026-08-10 16:54, run in progress) — S4 CLOSED at measurement
tier.** κM = 1.6e16 … 3.5e18 (vs greg5's 6.4e4 — eleven to fourteen decades worse). κ·ε_mach
spans ~3.5 … ~780: the magnetic solves sit at or beyond total-digit-loss in the worst
directions, so the original −46 dB residual floor was conditioning-limited, and tolerance
1e-6 is unreachable in principle at these κ. The trace also explains the old limp-along: κ
falls ~1 decade when Δt halves (A = M + ΔtK, cuts buy M-dominance), so the old controller was
trading timestep for conditioning all the way down the cliff; the cuts in this run correlate
with κ spikes (3.5e18 at the step-19 failure → 1.5e17 on the halved retry). The K-side
carries the pathology; prime suspect is the edge-coating wall wiring (greg5 = plain 2D
thin-shell at κ 6.4e4; sidecoatings adds the HEX8TB walls with µm thicknesses — the pending
R5 r′ wall-term theory). Discriminator: `edge coating : off` rerun, one deck key. κT = n/a as
designed (PETSc thermal).
