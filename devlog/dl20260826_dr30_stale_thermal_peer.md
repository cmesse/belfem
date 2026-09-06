# DR-30: stale thermal-peer temperature in superconductor postprocessing

**Date:** 2026-08-26
**Module:** fem/maxwell (postprocessing), fem/kernel (calculator linking — read only)
**Exchange:** `tmp/ai_exchange/dr30_stale_thermal_peer.md` (plan pre-registration, both audits,
reconciliation)

## What was fixed

`Calculator::link_element_maxwell_thermal` keeps the thermal peer on its *previous* element
when the current element is missing from the thermal group's aura (silent-keep, documented in
the code comment at `cl_FEM_Calculator.cpp:3016-3037`). `compute_superconductor` and
`compute_superconductor_ts` guarded only `block_exists`, then read the peer's `N(aK)`/`q()` —
another element's temperature — into `jc(...)` for the displayed J/Jc field.

Fix (Christian-approved, plan+audit → code+audit with Codex and Grok): in both functions, use
the peer's temperature only if the peer's linked element id equals the maxwell calculator's
current element id; otherwise `Temp = gTbulk`, the same contract as the existing
block-absent branch. Minimal diff; `Temp`/`T` names kept per minimal-rename policy.

## The load-bearing audit correction

The plan claimed the guard was defensive-only because the element-field loop is owned-only.
Both auditors independently refuted this: the *primary* caller is nodal field recovery
(`Postprocessor::run` → `recover_fields`), which walks `mMyElementIndices` — built with a
second **aura** pass for non-thin-shell blocks so recovery patches span the full disc at
partition boundaries (`cl_FEM_Postprocessor.cpp:252-264`). Each kernel expands its own aura
from its own owned set, so maxwell conductor aura elements can be absent from the thermal
group: the stale path is reachable today on ≥2 ranks. Claude verified the trace against the
tree. Thin-shell recovery (side-local, no aura pass) and `compute_element_data` (owned-only)
never trigger the guard.

## Evidence

- `g++ -fsyntax-only -std=gnu++17` with the target's `flags.make` flags: PASS.
- Static only — **reviewed, not verified**. Serial runs cannot exercise the fallback; the
  discriminating gate is the owed 2-rank thermal-coupled corc run (register row DR-30).
- Grok read-only check: target-file mtimes predate audit dispatch; tree clean.

## Residuals (kept on the DR-30 row)

1. T is not in `mPostprocessorSourceFields` — postproc correctness still rests on a preceding
   thermal solve+distribute.
2. `gTbulk` on maxwell-only aura cells biases partition-boundary nodal J/Jc when T is far
   from bulk; the physics-complete fix is expanding the thermal aura to cover maxwell's, not a
   postproc guard.
3. Pre-existing interpolation mismatch: postproc `norm(N*q)` (no clamp) vs assembly
   `dot(Nvec,q)` + clamp. Out of scope here.
