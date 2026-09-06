# Relative Tolerance: the A/B That Ended the Timestep Collapse, and the Default That Now Reaches STRUMPACK

**Date:** 2026-08-18
**Purpose:** Session record of (1) Christian's 1e-8 vs 1e-10 A/B on the tapestack3d
magnetic linear tolerance, which identified the timestep-collapse mechanism, and
(2) the audited code change making BELFEM's class default 1e-10 reach STRUMPACK
always, with the deck key demoted to an expert override.
**Threads:** `tmp/ai_exchange/reltol_ab_prereg.md` (pre-registration, RESULT, both
audit rounds), predecessor `step124_failure_investigation.md`
**Verification:** the A/B itself is executed evidence; the code change is
implemented + reviewed (two-vendor plan audit + two-vendor code audit, both
approving), riding the next rebuild. The running campaign states 1e-10 explicitly
and is unaffected.

## The A/B (Christian's experiment, over Claude's pre-registered objection)

Same restart state (t = 2300 ms dump), one knob: magnetic `relative tolerance`
1e-8 → 1e-10. Steps 165–170:

| | eps8 | eps10 |
|---|---|---|
| nonlinear iterates/step | 5–11 | **2** |
| Δt trajectory | 7.1 → 3.3 ms, then step-175 triple rejection to 0.15 ms | 7.1 → **12.3** ms, zero rejections |
| Krylov exits | 20 maxit-50 exits inside the collapse | none; ~16 its to rel.res ~3e-15 |

Mechanism, order-matched in the logs: **the printed nonlinear residual IS the
linear solve's exit residual.** At 1e-8 GMRES stopped on a ~4e-9 plateau and
Picard "stalled" at −85 dB — the controller then promoted to Newton on exit-test
noise, and the ω sawtooth / Δt collapse followed. At 1e-10 the next Arnoldi
steps drop to ~1e-15, Picard reaches its 1e-11 target in two iterates, and
Newton is never invoked. The collapse class chased across DR-88/DR-89/DR-90 was
substantially a measurement artifact of the linear exit test.

Claude's pre-registered floor model (absolute floor ~4e-9, predictions P1–P3)
was refuted on all three counts; the named reasoning error — generalizing from
20 pathological Newton solves while ordinary solves in the same file reached
1e-13 — is recorded in the thread. The executable gate beat the static analysis.

## The default change (both vendors approving, no blocker)

Discovered while planning: the class default 1e-10 (2026-08-13, thermal drift)
never reached STRUMPACK — the DR-76 have-flag gate (REFINE-era) meant an
unstated deck inherited the **library default 1e-6**, two decades looser than
the 1e-8 that collapsed. Same silent-contract-change pattern as DR-89: the
2026-08-16 switch to PREC_GMRES changed the meaning of "left alone" under a
guard calibrated for REFINE.

Landed (`strumpacktools.cpp`, `cl_SolverParameters.{hpp,cpp}`,
`cl_FEM_Controller.{hpp,cpp}`):

- `set_rel_tol( relative_tolerance() )` **always**, mirroring DR-89's abs_tol
  treatment; lineage comment ends "Do not re-introduce the guard."
- Unmapped Krylov methods (`cg`, `cgs`, `ibcgs`, `tfqmr`) on a STRUMPACK block
  are a **hard error** at setup — they used to fall through silently to the
  library's AUTO = Richardson-REFINE, the one path that could still reach the
  DR-76 stall geometry. Both auditors demanded error over silent remap.
- `mHaveRelativeTolerance` kept as a set-but-unread record with a poison-pill
  comment (positional synchronize payload stays width 12).
- The constant 1e-10 kept; coupling to the nonlinear tolerance rejected by both
  auditors (layering violation, unmeasured 1e-12).
- PETSc headroom gate NOT extended to STRUMPACK: the working default pairing is
  itself lin > nonlin and relies on GMRES overshoot.

Doc sweep (Input Contract, same session): `input_schema.yaml` +
`input_file_reference.md` rows for `relative tolerance` AND `krylov method`,
`solver_memory_and_compression.md`, the S9 comments in both Controller files,
and the tapestack3d deck comments. The retired theory — "STRUMPACK's exact
factorization already delivers ~2e-16, so a looser linear number is legitimate"
— was found in SEVEN places by the code audit (both vendors, convergent list)
after the first pass claimed it retracted; all fixed and re-gated. DR-76's
register row records the gate retirement; DR-88's evidence (gathered at 1e-8)
is flagged for re-measurement; DR-90's collapse traces are artifact-suspect.

## Open

- Unit test for the unmapped-Krylov error path (review-only until added).
- `set_relative_tolerance` lacks the >0/finite validation of its absolute
  sibling (pre-existing).
- `set_from_command_line` runs after the wrapper and can undo rel_tol/maxit
  (pre-existing; expert-paragraph material).
- Re-measure DR-88's first-iterate baseline and the J-asymmetry standing mode
  on eps10-era saves; revisit DR-90's scope.
