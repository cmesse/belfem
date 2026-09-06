# The 1.5e5 Preconditioned Residual: a Warm-Restart Anomaly, Three Audits, and a Lost Dump

**Date:** 2026-08-19
**Purpose:** Session record of the post-merge restart anomaly ( DR-92 ), its
three-voice sweep, and the discriminator that exonerated the merged binary.
**Threads:** `tmp/ai_exchange/bnorm_anomaly.md`; log record
`cmake-build-debug/tapestack3d/out_3100_anomaly.txt`

## What happened

The first run of the merged binary ( c9361d2d + betterdoc ), warm-restarted
from a t = 3100 ms / Δt = 50 ms dump, showed magnetic STRUMPACK solves with
`GMRES it. 0` around 1.5e5 ( vs O(1) in the previous run at the same step ),
every solve burning its full 50 Krylov iterations to exit at rel 1e-8..1e-5,
and step 195 rejecting twice — after holding a −109.05 dB iterate one dB from
target and losing 28 dB to noisy Newton corrections.

## What the sweep established ( Claude + Codex + Grok, unanimous )

- `GMRES it. 0 res` is the **left-preconditioned** residual ‖M⁻¹b_scaled‖ —
  solution-scale for Picard, correction-scale for Newton — proven at source
  level in STRUMPACK's GMResMPI.cpp / SparseSolverMPIDist.cpp. The raw ‖b‖
  is never printed; the original "70,000× jump in ‖b‖" framing was wrong
  ( Claude's own premise, corrected mid-sweep ).
- Refuted with citations: DR-91 seeding as a cause ( free dofs only; assembly
  reads mesh fields; Claude's Dirichlet-path mechanism specifically false ),
  the merged eigen changes ( no caller on any assembly/solve path ), any
  scale-touching code change ( identical min-pivot line = identical ‖A‖₁ ),
  and a memdump format change. Dump content scan: no NaN, physics-consistent
  norms, and ironically the cleanest φ history of any dump on disk
  ( max|φ| per BDF level = I(t) to four digits ).
- The Δt = 25 ms retry reproduced the anomaly ( kills "Δt = 50-specific
  operator" ) but its it.0 decayed 1.5e5 → 169 as iterates improved —
  consistent with the first restored-state corrections being genuinely huge.

## The discriminator

New binary + `memdump_2350.hdf5` ( Δt = 0.154 ms restored ): **healthy** —
it.0 = 0.0297, 16 Krylov its, Picard-2 at −131 dB, two iterates/step, thermal
entry −104.75 dB ( DR-91 verified a third time ). **The binary is exonerated**;
the anomaly follows the restored state. Leading hypothesis ( DR-92 ):
restart ≠ continuation at FULL step — every previously healthy restore
re-entered at ≤ 10.6 ms; this was the first 50 ms leap off a reconstructed
BDF stencil, at 58 A. Falsifiable for free at the next 50 ms-dump restart;
proposed code fix, if confirmed: cap the first post-restore Δt in
load_memdump.

## The loss, and the lesson

The t = 3100 dump — the only reproducer — was overwritten during test setup:
Claude handed Christian a restart recipe whose first line was the backup
copy, and it was skipped under a mistaken belief about the state's age. The
anomaly is now permanently unreproducible; its record is the preserved log.
New standing rule ( in Claude's memory ): **Claude makes the endangered-file
copies himself, in the same message that proposes any restart — never
delegated.**

## Consequence for the release

Main merge unblocked: full `make check` green on the merged tree, both code
audits clean, DR-91 execution-verified on two independent restores, and the
one anomaly shown to not follow the binary. DR-92 carries the open question.
