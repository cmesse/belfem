# Linalg R4.3: the gas and ODE consumers use wrapper expressions

**Date:** 2026-11-08
**Purpose:** sub-round R4.3 of `todo/linalg_wrapper_hardening_plan.md`, which belfem-de handed over and
Christian approved here. Christian was away and asked Claude to run `make check` itself.
**Module:** `src/physics/gasmodels`, `src/numerics/ode`
**AIs involved:** Claude Opus 5.5 (primary). Codex `gpt-6-astra`/high and Grok `grok-4.7`/high on the
code jury. Peer belfem-de (linalg) owns the plan text and Appendix D.
**Status:** COMPLETE.
**Verification:** executed by Claude at Christian's request in two fresh scratch trees (the shared
`cmake-build-*` trees were not touched): Blaze release with `USE_GASMODELS`,
`USE_NONFREE` and `USE_PYTHON` on, `make -j12 check` with no prior `make`, 22/22 tests, 2461 gtest
cases; Armadillo debug with `USE_GASMODELS=ON`, 19/19. Before that, the compile gate (16/16).
**Exchange:** `tmp/ai_exchange/review_linalg_wrapper_hardening.md` (R4.3 entries; diff `…_roundD3.diff`)

## Summary

**Change.** Eight hatch lines became wrapper expressions:
- `cl_Gas.cpp` ×5 (`aX ∓= mFormationTable * mFormationWork`)
- `cl_GM_EoS_Cubic.cpp:115` (`tDeparture += tY( k ) * mDepartureCoefficients( k )`)
- `fn_ODE_RK45.cpp:166` (`aY = yb`)
- `fn_ODE_DOP853.cpp:286` (`aY = ya`)

Each right-hand side is the backend expression the hatch already used; only the left-hand side moved
onto the round-C compound operators or the `Vector` copy assignment.

**Handoff corrections.**
- `cl_Gas.cpp` lives in `gasmodels/`, not `gastables/`.
- RK45 copies `yb`, its 5th-order stage, not `ya`. Appendix D's row was wrong; the code was right.
- My pre-registration claimed the gas lines gained a debug assert. They did not: the right-hand side
  was already the wrapper product.

**Gates.**
- Compile: each file compiled with `g++ -c` under its target's `flags.make`, on Blaze and Armadillo, in
  release and debug, with `-Werror -pedantic-errors`. 16/16, no diagnostics, and `nm -C` showed each
  object built against its own backend.
- Code jury: no defect.
- `make check`: the two trees in the header.

**Test gap noted by Grok, left open.** The compound-assignment value tests cover `M * s`, while
`remix()` uses `s * M`. Both go through the same overload, so the second form is an untested twin.

## Process notes

- **The plan file went in belfem-de's commit.** Within belfem-de's plan this session edited only the R4.3
  checkbox line and an appended outcome paragraph. belfem-de folded R4.3 into the Status line and
  committed the whole file in 6f8059a, so the R4.3 commit carries only the four sources.
- **Shared tree changed between messages.** belfem-de reported `cmake-build-debug` as Armadillo + debug
  at 11:51; by 12:19 it was Blaze + debug again. The gates did not depend on it.
- **The same session's DR-201 work** is recorded in the nonfree devlog, because it documents the nonfree
  test layout.

## Files Updated

`src/physics/gasmodels/cl_Gas.cpp`, `src/physics/gasmodels/cl_GM_EoS_Cubic.cpp`,
`src/numerics/ode/fn_ODE_RK45.cpp`, `src/numerics/ode/fn_ODE_DOP853.cpp`, `devlog/README.md`.
`todo/linalg_wrapper_hardening_plan.md` (R4.3 line and outcome paragraph) went in 6f8059a.
