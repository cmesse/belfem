# Cold Re-Read of the 2026-08-27 Strike Round

**Date:** 2026-08-27
**Purpose:** Christian asked for a cold re-read of the six debt-register rows struck today
(DR-41, DR-69, DR-76, DR-78, DR-84, DR-101) — six strikes in one day, two days before the
design freeze, is exactly the window where a strike taken slightly too fast is cheapest to
make. Three-voice round: Claude pre-registered, Codex and Grok audited blind in parallel.
**Exchange:** `tmp/ai_exchange/strike_cold_reread.md`
**Status:** read-only round; no source, register, or todo file touched.

## Method

Claude verified every checkable source claim in the six devlogs against the tree and
pre-registered a risk ranking (DR-69 first, DR-78 second, rest low) plus what would change
its mind, BEFORE dispatching both auditors with instructions to refute. The audit question
was framed as: does the cited evidence discharge the debt THE ROW described, or a different
proposition?

## Verdicts — unanimous across three voices

| row | grade | one line |
|---|---|---|
| DR-41 | SAFE | scope bookkeeping; fvm state matches the devlog word for word |
| DR-76 | SAFE | all four design claims present in `strumpacktools.cpp`; successors DR-88/89 live |
| DR-84 | SAFE | dispatch + both TS normal-fold paths verified; gate predates strike |
| DR-101 | SAFE | backup/flip/restore in all three PARDISO phases; red-to-green gate on record |
| DR-78 | QUESTIONABLE | parallel face implemented and gated; serial reporting face loud but ungated |
| DR-69 | QUESTIONABLE | evidence proves the control and the provenance, not the fixed table path in a solve |

No strike is UNSAFE: all six are process-correct closes of what Christian ruled. The two
QUESTIONABLE grades are about evidence coverage, not about the rulings or the code.

## The one substantive finding: DR-69's G3 fell out of the ledger

`dl20260816_dr69_unfolded_angle.md` names TWO remaining gates: G1 (constant-jc control) and
G3 ("tapestack3d frame A/B — expect asymmetric jc/ρ changes between the ±x tape halves; the
physics payoff gate"). Today's strike devlog says "both DR-69 gates discharged" — but counts
G1 plus the raw-data provenance check (Gate 2, done 2026-08-16). G3 is neither discharged,
waived, nor mentioned. Under constant jc the table path (`JcFunctionDatabase::eval` /
`wrap_angle` / the [0, π] coverage ctor) is bypassed by construction — `jc_eval` returns the
constant when `mJcFunction == nullptr` — so no executed gate ever showed the fixed table path
sampling the previously hidden pinning lobe inside a solve. Also: no committed test touches
`JcFunctionDatabase` or `wrap_angle` (only ModifiedKim evenness is tested), and the Gate-1
artifact directory `cmake-build-debug/tapestack3d_dr/` was already deleted by the time of
this re-read — the strike's numbers are a prose-only record, same day.

Mitigating, per Grok: table-Jc decks DID run after the 2026-08-16 unfold (the tapestack3d
campaign deck and the 2026-08-23 `2D_Tapestack` DR-84 gate frame both load `sp-ap.hdf5`) —
the fixed path likely executed; it was just never analyzed as G3. Cheapest discharge is
therefore a frame analysis of an existing output, not a new run.

**Decision owed (Christian):** re-open a slim DR-69 residual for G3, record G3 as waived, or
discharge it against an existing frame.

## DR-78, the second QUESTIONABLE

The np≥2 gate covers the pinned parallel mechanism exactly (rank-0-only `load_fields`,
history fields absent from `all_fields()`). The original SERIAL `qold(1)` bounds throw stays
mechanism-unpinned and has no red-to-green record; the per-rank length check
(`cl_FEM_Controller.cpp:4043`) makes the failure loud, not unreachable — and passes the
0==0 edge (empty parent AND empty history), where downstream qold guards are ASSERT-only
under NDEBUG. Acceptable as struck if the serial throw is accepted as unreproduced; the
0==0 edge belongs in DR-112's fix shape.

## Register-currentness by-catch (verified stale by read, NOT fixed — parallel session owns the file tonight)

1. Freeze-list DR-89 wording (`todo/debt_register.md:23`): still "library-default absolute
   tolerance"; abstol 1e-14 is always applied since 2026-08-17 — DR-89's remaining work is
   its run gate.
2. `todo/powerlaw_jc_n_field_derivatives.md:15-16`: β channel still "gated on DR-69's
   signed-angle work"; DR-69 is closed — pointer should move to the DR-07 residual.
3. DR-85's live description may predate the bn projection now in
   `compute_superconductor_ts` (`cl_MaxwellPostprocessor.cpp:864-867`) — medium confidence,
   check before editing the row.

## Confidence

High on every SAFE grade and on the G3-drop finding (direct reads). Medium on "the fixed
table path executed post-fix in the campaign runs" (binary rebuild dates not verifiable).
This round is reviewed-plus-verified-citations, not an executable gate: nothing here was run.
