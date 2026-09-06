# Evening campaign: riva jury, gauge Newton tangent, and the quench-front wall

**Date:** 2026-08-27 (evening/night session, continuous with dl20260827_overnight_dr_runs)
**Purpose:** Record the riva ship/archive jury, the gauge-tangent fix (planned, audited,
implemented, audit-corrected in one evening), and the configuration-elimination ladder that all
of tonight's tapestack3d_coarse runs together built at the quench front.
**Module:** physics/materials, fem/maxwell, fem/kernel (analysis only there)

## 1. riva jury (Christian: "Let's call the jury!")

Question: what can the riva law do that piecewise cannot; ship or archive. Full round in
`tmp/ai_exchange/review_riva_ship_or_archive.md` (pre-registration, blind Codex+Grok, verified
reconciliation). Headline verdicts, both auditors independent: **BELFEM's riva is totalized
Duron — the thesis's distinctive ρ_ηβ model is NOT implemented** (the tree's own
literature_references.md already said so); riva's real capability is numerical totality, which
`powerlaw` lacks at the tangent extremes and piecewise lacks until R8. Auditors split on
disposition (Codex: archive-for-freeze 3/5; Grok: ship-all-three with the §5 quench guidance
retracted, reject both archive options); both reject my pre-registered archive-powerlaw idea
(default-key migration in freeze week) — routed to Christian, undecided at close. The
resistivity_laws.md §5 line "quench studies: use riva" is unsupportable pending the wall
diagnosis and should be caveated before the freeze.

Verification found the riva-wall mechanism candidate both auditors missed: the failing iterate
missed the −90 dB target by 0.6 dB and the magnetic conditioning under riva was 2.7e19 vs
piecewise's 2-11e17 — a conditioning-mediated residual floor (H-κ), single-raiser, probe owed
(riva warm-restart from the 1050 ms memdump; R3 Δt-ladder of the floor plan).

## 2. n-clamp physics question (Christian: "What does physics say?")

Rhyner 1993's family spans n = 1 (ohmic) to n = ∞ (Bean): a clamp n ≥ 1+ε is the physically
correct limiting form, but the resulting ρ = ec/jc is garbage when jc is a sentinel — so the
deck-side fix (T_crit at the fit boundary) remains primary. Made partly moot the same evening by
the cleaned sp-ap/sst-1 tables (n interpolant ≥ 1.01 by Bernstein control net, jc floored) and
by the R-campaign's central n_eval floor: the hazard is confined to user plugins.

## 3. Gauge Newton tangent (plan → audit → code, one evening)

Christian recalled the skipped gauge Newton cross-term; source confirmed it
(`mt_maxwell_h.cpp`: power-law channel consistently linearized, gauge channel not), and
`coulomb_gauge_penalty_theory.md` §7 turned out to already contain the exact formula AND the
deliberate-omission policy. Plan `todo/gauge_newton_tangent_plan.md`; blind jury returned in
~50 min; Christian ordered coding before the verdicts ("while you code…"), audits landed
mid-implementation and corrected it in-place (per-element q() binding, Gq/GtGq scratch shapes —
Grok C7 adopted verbatim). Shipped: consistent j-channel in `h_newton_mu0`/`h_newton_mu`
(covers 2D/3D/thin-shell — same kernels), theory doc §7 updated, both TUs compile -Wall -Werror.
**Audit-forced retractions recorded:** the term is O(χ) against the physical tangent at this
deck's chi = 1e-4 (not "95 % wrong" — that was the action-on-q ratio against the gauge block
alone), and the 4.17 s wall was thermal-first, not this term's signature. Owed: G-FD
directional-derivative probe (the only true verification), `make check-fast` — both blocked
tonight by the in-flight libbelfem CMake refactor (executable/library target `belfem`
collision at configure).

## 4. The configuration-elimination ladder (tapestack3d_coarse, quench front t ≈ 4.2–4.6 s)

| config | outcome |
|---|---|
| piecewise Newton + PETSc thermal + chi (clean) | walls 4.17 s: PETSc DIVERGED_BREAKDOWN (960 its), thermal watchdog cuts; recovered at 0.5 ms, ground on |
| piecewise Picard-only, chi off | antiphase coupled limit cycle (mag/thermal alternating blowups), Δt cut is the correct lever there |
| thermal Newton (earlier trace, 4.6 s) | thermal Newton anti-converges while its Picard descended — promote-destroys-convergence, thermal has no demote latch |
| **gauged Picard + thermal MUMPS (G1, tonight)** | **advances 4.18 → 4.25 s grinding, then walls at 4250.9 ms: thermal CONVERGED (−63 dB, MUMPS clean, zero breakdowns), magnetic Picard non-contractive at Δt down to 2 µs** |

Elimination result: thermal is fixed by MUMPS (DR-77's second measurement — T-field crosscheck
vs the PETSc run still owed before wider trust); the coupling was not the binding constraint in
G1 (thermal converged); **Picard is structurally insufficient at the quench front** (fixed-point
radius ≥ 1 in the current-sharing state, Δt-invariant). The one untested configuration is the
one every eliminated line points at: **Newton magnetic + thermal MUMPS + chi + the consistent
tangent** — buildable as soon as the libbelfem configure collision is resolved.

## 5. Register/ruling state at close

- "Gauging on by default" (Christian's conditional ruling): still pending — G1's wall means no
  configuration has passed yet; the gauge-tangent objection is resolved, the O(χ) finding means
  chi was never the big lever; the ruling's basis will be whichever configuration traverses.
- riva disposition: with Christian (jury split recorded).
- DR-77: thermal-MUMPS-beside-magnetic-MUMPS ran breakdown-free for hours on the coarse deck —
  amendment material once the T-field check runs.
- Follow-ups owed: G-FD probe; O2 (gauge field-channel, HTS+metals); riva H-κ probes;
  thermal demote-latch + coupled-mode discriminator (floor plan); resistivity_laws.md §5 caveat.
