# R6 coupled restart, and one deck key that moved two ways at once

**Date:** 2026-08-29
**Purpose:** record the DR-124 R6 gate ( conducted, passed ), the Δt-consistency defect it
exposed ( folded into DR-112, plan at `todo/restart_timestep_consistency.md` ), and the
phase-sign convention change ( DR-148 )
**Module:** fem/kernel, numerics/sources, circuit

## DR-124: the gate that was not a formality

R6 asked whether a dumped circuit state survives into the first post-restart solve. It ran
on the `belfem` executable — `hphirun` is being retired — at np=4, against
`examples/RLC_Circuit` shortened to 8 ms. The shortening was forced: the time loop writes
an Exodus mesh *and* an 18 MB memdump every timestep, ~29 s of wall per step against ~2
Picard iterations converging to 2.2e-16. The seam is what R6 tests, not the duration.

**The control came first.** The 4 ms leg and the 8 ms baseline are bitwise identical over
all 27 common rows and all 7 output channels. The run is deterministic at np=4, so anything
after the seam belongs to the restart and not to noise. Without that control the rest of
this entry would be uninterpretable.

**The result: the restore is exact.** With the first post-restart step taken at the dumped
timestep, the resumed run is bitwise identical to the uninterrupted baseline over all 20
compared rows, every channel, zero deviation. DR-124's v2 per-component restore is verified
at coupled level. Row struck.

## What the gate actually caught ( filed as DR-149, then WITHDRAWN )

The first attempt did not look like that at all. With the deck's own `initial timestep` the
first post-restart step advanced the clock by 1.0e-5 s and returned currents **bitwise
identical to the baseline row 2.0e-4 s after the seam** — twenty times the time advance —
while the node voltages correctly tracked the new time ( `V0 = 0.961389 =
sin(2*pi*50*0.00411258)` ). That step rose 144 A where its neighbors rise about 4.6 A.

Sources evaluate at the new time; the reactive components integrate a full old timestep.
The trigger is the deliberate DR-92 cap in `load_memdump`, which holds the first restored
step at `mDeltaTimeInitial`. Grok's caveat on the DR-124 row predicted this interaction
before anyone ran it.

The dt-matched rerun is what separates the two findings: state restore is exact, and the
defect is Δt-consistency on the first step. It is **deliberately not localized** — the
controller does call `set_timestep -> shift -> compute_MNA_matrix` on that step, so the stale
Δt more plausibly arrives through the restored magnetic BDF history feeding the FEM-coupled
terminal pair. Naming a module without evidence would have been a guess.

**It was first filed as DR-149 and withdrawn the same session.** Christian's ruling: the
register is inflating one ID per finding, and a session should extend the row it is already
working rather than branch a new one. The evidence moved into DR-112 — whose loader contract
this genuinely is — and the work into `todo/restart_timestep_consistency.md`, whose R1 is
localization and blocks every other step. The policy is now recorded in the register's own
rules. That is the better call: DR-112 and this defect are two symptoms of one contract, and
splitting them would have meant two rows fixed by one change.

One hypothesis died on measurement: "the physics is permanently ahead of the clock."
Re-comparing against a baseline shifted by the 1.9e-4 s lead makes the fit *worse*
( 452 A vs 137 A ). It is one bad step that then decays, not an accumulating offset.

## DR-148: a key that advanced one waveform and delayed three

Codex found it as by-catch of DR-145 and Grok independently flagged it; both said keep it
out of that row. Probed: at `phase = +pi/2, T = 1` the sine's peak moves from `t = 0.25` to
`t = 0` while the sawtooth's reset moves from `t = 0` to `t = 0.25`. Same deck key, opposite
directions, because the sine reads `phase()` through `sin(omega*t + phase)` and the other
three read `aTime - time_offset()`.

Christian ruled: BELFEM follows the SPICE definition, a positive phase advances. I had
recommended the opposite — delay internally, negate at the ngspice boundary — and was
overruled. The ruling is the better contract: it needs no negation layer, and it leaves the
ngspice passthrough already correct, since SPICE `SIN` advances too.

Landed as three sign flips. Sine, both setters, the ngspice factory and ramp/sigmoid are
untouched; ramp and sigmoid are genuine start delays on the separate `offset` key, which
both auditors confirmed no factory ever feeds `phase` into. The triangle's `+3T/4` is a
shape origin, not a phase term, and stayed.

The proof is one identity: `f_phase(t) == f_zerophase(t + T/4)` now holds at exactly
0.000e+00 for all four waveforms. Before the change it held for the sine alone. That test
is the one that would have caught the original inconsistency, and it is now in the suite.

Because this changes what a deck key *means*, both input-contract artifacts were updated
along with `circuit_usage_guide.md` — the one row in this session where the two-artifact
rule genuinely fired ( DR-117 and DR-145 did not need it; both vendors agreed on all three
calls ).

## What the audits caught that I would not have

- **An ID collision.** A parallel session filed DR-146 and DR-147 while this round was in
  flight. Grok checked the live register and caught it before a duplicate DR-146 was
  created. The register is contended; ids must be re-read at filing time, not at planning
  time.
- **My test table was measured against the wrong tree.** I computed the "today" column as
  advance-without-wrap when the tree was delay+wrap, which turned both of my claimed guards
  into discriminators. The real guard is `phase = 0`, invariant under both the wrap and the
  sign.
- **A counterfactual comment.** The `+1.8` overshoot note described reverting the wrap while
  keeping the new convention, not any tree that ever existed. Corrected, and the
  intermediate tree's actual `-0.2` verified.

## Status

DR-124 struck. DR-148 landed, probe-verified and approved by both vendors at code stage,
`make check` owed. The Δt-consistency defect lives in DR-112's extended scope with a written
plan and a reproducer in hand; DR-149 was withdrawn. `check_doc_claims` 36/36 — note it counts `[P]` rows but not `[F]`, and
the `[F]` prose needed a manual bump to three. Uncommitted; `example_user_source.cpp` is
staged in the shared index by another session and must not be swept into a commit here.
