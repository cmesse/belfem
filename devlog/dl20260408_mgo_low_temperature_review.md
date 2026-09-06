# Devlog 2026-04-08 — MgO Low-Temperature Review

**Date:** 2026-04-08
**Topic:** Read-only investigation of negative `cp`/`lambda` output for the new Magnesia material
**AIs involved:** Codex
**Codex Audit Confidence:** high

## Summary

Reviewed the new `Magnesia` material implementation after the `material mgo` table showed `-0.000` for `cp` and negative `lambda` values at very low temperature. The negative sign does not come from a physically negative closed-form property law. It is introduced by the interpolation layer used for printing. For `cp`, the observed `-0.000` at `0 K` is a signed-zero / roundoff artifact. For `lambda`, the situation is more serious: the cubic spline built on 4 K tabulation points overshoots below zero between `0 K` and `4 K`, and the underlying MgO conductivity model also stores switch points in `log(T)` but compares them directly against `T`.

## Key Findings

- The printout uses the spline-backed accessors, not the raw MgO closed-form functions. `Magnesia::create_cp()` and `Magnesia::create_lambda()` both call `create_spline(...)` in [cl_Material_Magnesia.cpp](/home/christian/codes/belfem/src/physics/materials/cl_Material_Magnesia.cpp#L124) and [cl_Material_Magnesia.cpp](/home/christian/codes/belfem/src/physics/materials/cl_Material_Magnesia.cpp#L191), and `Material::set_spline()` switches evaluation to `cp_spline()` / `lambda_spline()` in [cl_Material.cpp](/home/christian/codes/belfem/src/physics/materials/cl_Material.cpp#L433).
- The spline is generated on a fixed 4 K grid in [cl_Material.cpp](/home/christian/codes/belfem/src/physics/materials/cl_Material.cpp#L588). `Spline::eval()` is an unconstrained cubic polynomial evaluator in [cl_Spline.hpp](/home/christian/codes/belfem/src/numerics/spline/cl_Spline.hpp#L432), so it does not preserve positivity between tabulated points.
- `cp` itself is nonnegative in the raw MgO model. `gamma` is set to `0.0` and `beta` positive in [cl_Material_Magnesia.cpp](/home/christian/codes/belfem/src/physics/materials/cl_Material_Magnesia.cpp#L41), and the low-temperature law is `(gamma + beta*T*T)*T` in [cl_Material_Magnesia.cpp](/home/christian/codes/belfem/src/physics/materials/cl_Material_Magnesia.cpp#L237). The printed `-0.000` at `0 K` is therefore just endpoint roundoff.
- `lambda` is strictly nonnegative in the raw MgO model as written, because every nonzero branch returns `exp(...)` in [cl_Material_Magnesia.cpp](/home/christian/codes/belfem/src/physics/materials/cl_Material_Magnesia.cpp#L258). So the reproduced negative values at `1 K`, `2 K`, and `3 K` must come from the spline interpolation, not from the closed-form conductivity function itself.
- There is an additional MgO conductivity bug: `create_lambda()` computes curvature-switch roots in the `x = log(T)` coordinate, but stores them directly in `mTLambdaSwitch` at [cl_Material_Magnesia.cpp](/home/christian/codes/belfem/src/physics/materials/cl_Material_Magnesia.cpp#L200) and later compares them against `T` in [cl_Material_Magnesia.cpp](/home/christian/codes/belfem/src/physics/materials/cl_Material_Magnesia.cpp#L261). The roots are therefore used with the wrong units. This distorts which branch is sampled into the spline and likely worsens the low-temperature behavior substantially.
- Reproduced with the current binary: `./cmake-build-debug/bin/material mgo -t 0 8 1` prints `lambda = -185.204` at `1 K`, `-229.105` at `2 K`, and `-131.705` at `3 K`, while `4 K` jumps to `106.998`. That pattern is consistent with cubic overshoot in the first spline interval.

## Changes Made / Proposed

- No source edits.
- Wrote this devlog to record the root cause before any fix is attempted.

## Open Questions

- Should MgO `lambda(T)` remain on the custom path instead of being resampled onto the generic 4 K spline?
- If spline interpolation is still desired, should it be positivity-preserving in log-space rather than an unconstrained cubic in temperature?

## Files Updated

- /home/christian/codes/belfem/devlog/dl20260408_mgo_low_temperature_review.md
