# BDF Warm Restart, Second Round: the History the Fields Never Carried

**Date:** 2026-08-16 (overnight campaign)
**Purpose:** Record the DR-78 investigation and fix — the warm restart's
first-ever execution with history present crashed; the restore side was
the defect, and the save side was accused, convicted, and acquitted in
one night
**Module:** fem/kernel
**Round:** `tmp/ai_exchange/bdf_history_restore_plan.md` +
`bdf_history_restore_impl.patch`

## The crash

The first warm restart from a history-carrying dump (tapestack3d, five
completed 1 ms steps) died on the first BDF5 assembly — parallel as a
rank-1 SEGV under NDEBUG, serial (Christian's lldb repro) as a bounds
throw in `Calculator::qold(1)`: a length-0 vector behind a non-null qold
table entry. The integrator restore itself was correct; the *fields* it
promised were not there.

## What the investigation established

- **Parallel, proven:** `load_fields` is rank-0-only and the restore
  synch uses `all_fields()`, which never contains the numbered history
  labels (`phi0..4`, `edge_h0..4`, `T0..4`). Workers enter BDF5 with the
  empty shells `init_qold_table` created. The guarding asserts in `qold`
  are compiled out under NDEBUG — hence a raw SEGV, and had the shells
  been sized, a silently wrong right-hand side.
- **Save side: accused, then acquitted (F7, retracted in phase 3).** The
  intermediate finding "nothing collects history before save_memdump"
  was true but its conclusion was wrong: `DofManager::solve()` ends
  every step with `distribute( all_fields )` after the master writes the
  full solution, so rank 0's fields are solution-fresh each step and its
  shift builds the true history. A collect preamble was implemented on
  the false conclusion and REMOVED the same night when Codex's phase-3
  audit (which caught that `collect_fields` skips EDGE/FACE entirely)
  forced the re-trace that refuted F7. The save side was never broken.
- **Serial, eliminated but unpinned:** four hypotheses (mesh-identity
  aliasing, double-creation, missing dof-field history, shift ordering)
  all refuted — Codex statically, the rest by a miniature reproducer
  (scratchpad copy of Tape_Quench) that *refuses to crash*: serial warm
  restart works there, even from a dump with a history level deliberately
  truncated to length zero (the resumed shift heals shallow holes before
  any read). Production's serial mechanism needs an ingredient the mini
  deck lacks; the strongest surviving suspect is parallel-written dump
  state. Two probe attempts were vacuous (resumed time ≥ simulation
  horizon — zero steps assembled) before the harness was corrected;
  recorded as its own lesson.

## The fix (restore side)

- The fix is **restore-side only** (the save side needed nothing — see
  above): `Controller::synchronize_history_fields()`,
  called inside the `tHaveBdf` guards after the state broadcast and
  **before** `restore_history_state`: distributes the numbered levels to
  the workers (their field shells exist — `initialize()` creates them
  first) and then hard-verifies, on every rank, every level the resumed
  step's SHIFT will read — depth `min(count, order−1)`, one deeper than
  the scalar restore's own bound because the shift pushes levels down
  before `collect_qhist` reads (Codex phase-3 catch) — length-matched to
  the parent field, `BELFEM_ERROR` naming field, rank, dump, and
  promised depth. Deeper levels may be legitimately empty on a mid-ramp
  dump and refill through the shift before any read. The label
  convention got a single source: `IWG_Timestep::history_field_labels`
  over `mDofLabels` — the list `shift_fields` actually rotates, which
  `dof_fields()` subtly is not.
- Old dumps without bdf keys: both additions dormant; cold starts
  untouched.

The verification sweep is also the instrument this campaign was missing:
whatever produced production's length-0 level, the next real warm restart
converts it from a SEGV into one line naming the culprit field.

**The fix's own near-miss.** Grok's audit (its CLI recovered the next
morning) caught that the first implementation distributed *all* history
levels while verifying only the ones the resumed step reads — so on a
parallel mid-ramp restart the deeper, never-filled EDGE/FACE levels would
reach `FieldData::distribute`, which walks entity indices with no empty
guard, and reproduce the exact crash the function exists to prevent. The
serial probe could not see it (no worker loop) and neither could the
saturated production dumps (all five levels filled). The distribute list
is now built at the verified depth, so the two sets are identical by
construction. Grok also confirmed the F7 retraction independently —
tracing master-only field writes through the rank-0 scatter, including
the rejected-step and final-step-before-save cases I asked it to
attack — and refuted, correctly, my prose claim that the depth is
"one deeper than the scalar restore's r" (true mid-ramp only; the two
agree once the ramp saturates).

## Process

Plan + Codex audit (Grok's CLI died twice at dispatch — recorded, owes a
phase-3 pass); probes run in the scratchpad against the prebuilt binary,
never touching the shared tree. The production machine also delivered its
own lesson tonight: three OOM kills of the study, the last explained in
part by hand-launched runs missing `OMP_NUM_THREADS=2` (STRUMPACK showed
T=4 in the factorization trace — Allrun exports 2).

Status: **reviewed, round closed with both voices** — Codex phase 3
(three blockers: F7 retraction, depth off-by-one, label single-sourcing)
and Grok's combined pass (the parallel mid-ramp distribute hole, plus
independent confirmation of the F7 retraction and the depth formula),
all resolved same-session. Syntax gates green on both TUs. Executable
gates at the next rebuild: mini-deck A/B on the fixed binary with an
np=2-written dump AND a **parallel mid-ramp restart** (dump written at
step 1-4 of BDF5, restarted on np=2 — the case both probes were blind
to);
`make check-fast`; a production warm restart.
