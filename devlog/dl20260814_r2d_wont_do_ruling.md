# R2d Converged-at-Clamp Diagnostic: Won't-Do (DR-02 Is Purely a Run Again)

**Date:** 2026-08-14
**Purpose:** Record Christian's ruling against the proposed Controller diagnostic,
and the evidence that decided it.
**Module:** `todo/` only — no source touched.

Claude proposed landing the deferred R2d half of the kernel-collapse plan: sticky
clamp counters in the Calculator plus a once-per-step "converged with N clamped
evaluations" line in the Controller, consuming the `T_clamped()`/`rho_clamped()`
accessors that have had zero consumers since R2d's clamp half landed.

Christian was not convinced and pointed at the live evidence,
`cmake-build-debug/tapestack3d/out.txt` (63k lines, 8 procs, BDF5, the production
quench workload). The printout decided it against the proposal:

1. **The run never clamps.** The field sits near 77 K throughout; the diagnostic
   would have printed nothing in exactly the workload it was meant to serve.
2. **The path into clamping is already voiced.** Runaway T degrades convergence
   before it reaches a ceiling, and `print_thermal_stall_warning` covers that;
   the printout's healthy signature (−90 dB residuals, orderly Δt adaptation,
   clean rejection-retries at steps 343/349) is what the warning's absence means.
3. **A genuinely converged-at-ceiling state is visible in the results**: a
   `T_max` plateau in the exodus T field, which quench analysis always inspects.
   The "silent wrongness" scenario required nobody ever looking at T.

Against two existing observability channels, a third was not worth perturbing a
Controller that has just reached a stable state after the Anderson,
false-convergence and PID campaigns.

Recorded: the plan's R2d box is closed with the strikethrough-and-ruling
convention; the DR-02 row drops its "one live code item" amendment and is
**purely a run gate again**. The accessors stay in the tree — unconsumed,
trivially small, available if a future campaign wants them.
