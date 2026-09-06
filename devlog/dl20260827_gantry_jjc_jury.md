# Gantry J/Jc = O(1000) adjudicated by jury; deck aligned to the design paper

**Date:** 2026-08-27
**Topic:** gantry example — memdump J/Jc magnitude, modeling vs post-processor; deck alignment
**Exchange:** `tmp/ai_exchange/review_gantry_jjc.md` (pre-registration, Codex audit, Grok
partial, verification, reconciliation), brief in `tmp/ai_exchange/gantry_jjc_brief.md`

## Question

`cmake-build-debug/gantry/memdump.hdf5` (t = 9.8 s, I = 340.6 A, T = 77 K) shows |J/Jc|
median 409, max 7315 on the 19488 conductor nodes; expectation was O(1). Modeling problem or
post-processor problem? Adjudicated with a blind jury round (`cross_review.sh --jury` on an
evidence brief).

## Verdict — unanimous (Claude pre-registration, Codex full audit, Grok partial): modeling

- **The postprocessor reports the model faithfully.** Implied jc = Jz/JJCz per node
  (3.1e4..4.4e7 A/m²) sits exactly inside the loaded table's 77 K manifold
  (`share/material/bscco-2223.hdf5`, log10(jc) storage, 76 K slice = 1.2e4..8.9e7 A/m²).
- **Solve and postproc agree by construction.** `compute_superconductor_ts`
  (`cl_MaxwellPostprocessor.cpp:798`) and the solve-side `compute_rho_piecewise_ts`
  (`cl_FEM_Calculator.hpp:3400`) use the same unfolded `bn_angle` and the same
  `JcFunctionDatabase::eval`. All Codex file:line citations re-checked and confirmed.
- **The solver had already quenched the tape.** `element_rho` ceilings at 1.072e-7 Ω·m =
  the piecewise law's normal-state branch (the global rho clamp defaults are a no-op,
  `cl_Communicator.cpp:55`); |Jz| = I/A_tape uniform, no critical-state profile.
- **Root cause:** 340 A at 77 K in a 1–5 T field through a tape whose loaded table says
  Ic(77 K, sf) = 115 A. The deck's "Ic = 495 A" comment predated the table.

Evidence rung: probe (independent h5py reads by two parties + source trace). Reviewed, not
verified — no solve was re-run. Process notes: Grok hit its 30-turn ceiling before writing a
formal entry; its partial transcript independently confirmed B-reconstruction and angle
parity and is weighted as a partial voice. `git status` clean after the round.

## The legacy-table comparison that explained the old O(1) expectation

Old dataset (`tmp/bscco/bscco.hdf5`, converted 2026-08-27 from the 2022-era data):

1. **No data above 33 K** (T axis 4.2..33 K) — a 77 K run clamped to the 33 K edge.
2. **jc values ~1000× inflated**: at matched 33 K / 1.4 T, old/new = 1291. Old jc × its
   intended 1.98e-6 m² area gives Ic(33 K, sf) = 735 kA (absurd); old jc × 1.98e-9 m²
   gives 930 A at 4.2 K — right next to the new table's 750 A. The legacy conversion
   normalized over an area 1000× too small (0.396 mm thickness slipped to µm scale).

So legacy runs showed J/jc ≈ 0.0004..0.055 at "77 K" for two stacked spurious reasons, which
is where the O(1)-or-below intuition came from.

## Deck aligned to the design paper (edits applied on Christian's request)

Paper: J. L. Rudeiros Fernández et al., "Mechanical and Thermal Analysis of an HTS
Superconducting Magnet for an Achromatic Gantry for Proton Therapy", IEEE Trans. Appl.
Supercond. 32(6), 4401805, 2022, doi:10.1109/TASC.2022.3158366 (local copy
`tmp/gantry.pdf`). Operating point 340 A at **12 K**; 464 turns (4 double-pancake coils × 2
layers × 58 turns — matches the deck); peak |B| in conductor 2.42 T, peak By midplane 2.46 T;
design tape DI-BSCCO HTi-CA, Ic(77 K, sf) ≥ 170 A (193 A measured), Ic(12 K, |Bx|) = 485 A →
the paper's 1.43 margin at 340 A.

`cmake-build-debug/gantry/input.conf` changes:

- header: design-paper citation with DOI added
- `temperature : 77 K` → `12 K` (Table II)
- `amplitude : 347.2 A` → `340 A` (Table II, I_op)
- stale "Ic = 495 A" material comment replaced by the paper's tape data plus the caveat
  that `bscco-2223.hdf5` holds the weaker AMSC 115 A tape (Turrioni et al. 2008), so the
  simulated margin is smaller than the paper's

Predicted with the current AMSC table at 12 K / 340 A: J/jc = 0.50 (sf, best angle) to 1.14
(2.42 T, B ⊥ tape); Ic(12 K, 2.42 T) = 298 A worst-angle vs the paper tape's 455 A. Expect
O(1) values with possible local over-critical spots at perpendicular-field corners until a
DI-BSCCO Ic(B,θ,T) table replaces the AMSC one. 12 K is inside the table's measured range
(4.2–33 K), unlike the 77 K extrapolation corner.

## Open items

- Rerun the gantry deck at 12 K / 340 A; validation anchors from the paper: peak By on the
  midplane 2.46 T, peak |B| in conductor 2.42 T, J/Jc ≈ 0.7 at the load-line point.
- Tape geometry mismatch left in place (mesh has 5 mm tapes at 0.414 mm pitch; paper tape is
  4.6 mm × 0.4 mm insulated) — a remesh decision, not a deck edit.
- P2 (Codex, single-raiser): `JcFunctionDatabase` owns a raw `Database*` without deleted
  copy/move — latent copy-safety, unrelated to this symptom.
- P2: legacy `tmp/bscco/bscco.hdf5` carries the ×1000 normalization error; add a provenance
  note if the file is kept.
