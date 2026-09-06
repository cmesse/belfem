# 2-D current-BC sign inversion — fix plan

**Date:** 2026-08-28
**Status:** DONE and gate-verified. Fix B landed, input contract updated in
the same session, both audit rounds closed (Codex + Grok; Grok refuted the
voltage *mechanism* and the durable record was corrected). Gates O1-O4 all
passed on the rebuilt binary with fresh cut generation. Remaining: **O5** debt
row and **O7 / G-V1**, a 2-D voltage-driven run, both owed and neither
blocking. Changes are UNCOMMITTED.
**Exchange:** `tmp/ai_exchange/current_sign_2d.md` (ephemeral)

## Defect

A declared positive `current` amplitude on a 2-D deck produces a physically
*negative* current (into the page). Measured on `cmake-build-debug/gantry.exo`
(all 464 tapes declared +340 A): every tape carries exactly -340.0 A by Jz
integration, confirmed by Ampere loops on H (-19720 A per 58-tape coil).
3-D is correct (corc: current flows input terminal -> output terminal).
Post-processing is internally consistent; the inversion is on the solve side
(cut cochain orientation from the suggested homology).

## Steps

- [x] D1 — Plan audit: **Codex** confirmed the sign chain independently and
      endorsed the repair point. **Grok produced no verdict** — exhausted its
      30-turn budget ($0.27) with no final message; retried on a narrowed code
      audit instead.
- [x] D2 — **G-B1 discriminator — DONE 2026-08-28: bulk is ALSO inverted.**
      Answered without a new run from `cmake-build-claude/costhetatest/
      hphi_results.e-s.00001` (68 bulk `coil` blocks, +169.96 A declared):
      CCW Ampere loops on air-side H give -11414..-11420 A over three nested
      boxes, ratio **-0.988** to the declared +11557 A. Uniform, no mixed
      signs => **Fix B selected**, Fix A excluded.
- [x] R1 — **Fix B applied**: `Homology::reorient_generators` now multiplies
      every suggested 1-generator by `-1`; the dimension branch is gone (3-D
      already had `-1`, so it is a no-op there). Comment rewritten with the
      measured evidence. `src/homology/cl_Homology.cpp`, anchor
      `reorient_generators`. ~~Fix A~~ excluded by G-B1.
- [x] R2 — Input contract: document "2-D positive current = +z (out of plane),
      right-hand rule" in `doc/input_file_reference.md` AND
      `doc/input_schema.yaml` (same session; behavior definition of a key)
- [x] R3a — Code audit (Codex): no blocking defect. Three findings accepted
      and applied (doc wording overstated "only mirrored"; comment trimmed;
      `operator*(-1.0)` → `-1`). One rejected with reason (the "unrelated
      prose" hunks predate this session — earlier prose sweep, not mine).
- [x] R3b — Code audit (Grok): **refuted Claude's mechanism, confirmed the
      conclusion.** The "e_i -> -e_i on an assembled system, both I and V
      columns negate" story is wrong — the system is REASSEMBLED, declared
      scalars enter the new coordinates unchanged, so save_IV's columns stay
      and their geometric meaning flips. Durable record corrected. Also found:
      the missing "+V drives +z" doc statement (added), the stale "terminals
      point inwards" comments in cl_CutFactory.cpp (rewritten), and a
      pre-existing FALSE comment at cl_FEM_Controller.cpp:501 stating the
      current/voltage ordering backwards (rewritten, by-catch)
- [x] O1 — **PASSED**: fresh-mesh costheta on the rebuilt binary, Ampere ratio
      -0.988 -> **+1.0006 / +0.9984 / +0.9918**
- [x] O2 — **PASSED**: gantry fresh mesh, 0.2 s, user's own mixed-polarity
      deck. curves 1:232 declared -340 -> 232/232 tapes at -3.74 A; curves
      233:464 declared +340 -> 232/232 at +3.74 A (ramp value at t=0.1 s);
      **0 sign errors**. Also clears the positional abstract-dof ordering
      contract across 464 separate conditions, untestable until now
- [x] O3 — **PASSED**: corc 3-D no-op verified by execution, not only by
      inspection. Ampere ratio to declared current -0.6851 pre-fix vs
      **-0.6852** post-fix (operating points 186 A vs 0.48 A); direction still
      input terminal (z=0.0377) -> output (z~0)
- [x] O4 — **`make check-fast` 10/10 green** (homology 6.08 s). As predicted
      this guards against collateral damage only; no test asserts a current sign
- [ ] O7 — **G-V1 (owed, Codex's request):** 2-D voltage-driven regression —
      positive declared V drives +z current; saved I*V has the expected
      dissipation sign; current- and voltage-driven runs agree on IV quadrant
- [x] O5 — **DR-133 filed** in `todo/debt_register.md`: `orient_terminal_curves_2D`
      is geometrically vacuous (`tR` reduces identically to `-(facet tangent)`;
      never consults the enclosure side). P2, deliberately NOT fixed in this
      round — repairing it in the same diff would have confounded the
      discriminator that selected the sign fix
- [x] O6 — Devlog `devlog/dl20260828_current_sign_2d.md` + README index

## Gap table

| # | gap | risk | covered by |
|---|-----|------|-----------|
| 1 | bulk 2-D sign unknown | fix at wrong site flips bulk wrongly | D2 (G-B1) |
| 2 | additional hidden sign between cochain and phi-jump | one-line fix lands on wrong layer | D1 audit Q1 |
| 3 | 2-D circuit / save_IV voltage co-flip | inconsistent I-V output | R3 audit |
| 4 | existing 2-D decks change meaning | user surprise | pre-release; sole user requested it; O2 restores gantry deck |
