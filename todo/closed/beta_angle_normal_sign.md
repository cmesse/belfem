# HTS Beta Angle: Possible Sign Error from Uncontrolled Normal Direction

**Date:** 2026-02-25
**Purpose:** Document potential bug in β angle computation for HTS thin shells
**Module:** fem/maxwell

## The Issue

In `mt_maxwell_h.cpp`, all `h_ts_hts*` functions and `h_ts()` compute the
angle between the magnetic field B and the tape normal n:

```cpp
const Vector< real > & n = tCalc->normal();
beta = norm_b < 1e-6 ? pi/2 : acos( clamp( dot(n, b) / norm_b, -1, 1 ) );
```

The normal `n` comes from `tCalc->normal()`, which is computed from the
master element's face via `mMasterIndex = facet->index_on_master()`. The
sign of this normal depends on which volume element is master — a convention
that is **not controlled for a specific outward direction**.

If n happens to point "downward" (into the conductor) rather than "upward"
(out of the conductor), then `dot(n, b)` flips sign, giving β → π−β.

## Where It Appears

9 functions in `mt_maxwell_h.cpp`:

| Line | Function |
|------|----------|
| 454 | `h_ts_hts` |
| 611 | `h_ts_hts_defect` |
| 768 | `h_ts_hts_piecewise` |
| 925 | `h_ts_hts_defect_piecewise` |
| 1101 | `h_ts_hts_t` |
| 1294 | `h_ts_hts_defect_t` |
| 1487 | `h_ts_hts_piecewise_t` |
| 1680 | `h_ts_hts_defect_piecewise_t` |
| 1988 | `h_ts` (dispatcher) |

## Does β → π−β Matter?

The standard elliptical Kim-type model decomposes the field into:

```
B_perp = B * cos(β)     (perpendicular to tape, along c-axis)
B_par  = B * sin(β)     (parallel to tape, in ab-plane)
```

and computes:

```
Jc = Jc0(T) / (1 + sqrt( (B_par / Bk)^2 + (B_perp / Bc)^2 ))^alpha
```

Since `sin(β) = sin(π−β)` and `cos²(β) = cos²(π−β)`, any model that uses
B_perp² or |B_perp| is invariant under β → π−β. The standard models should
be safe.

**However**, if any material model uses `cos(β)` linearly (not squared or
with abs), the sign matters. This would be unphysical (a tape has no
preferred "up" vs "down"), but could exist as an implementation oversight.

## Recommended Fix

Replace the β computation with a sign-invariant version:

```cpp
// Option A: use abs on the dot product
beta = norm_b < 1e-6 ? pi/2
     : acos( clamp( std::abs( dot(n, b) ) / norm_b, 0.0, 1.0 ) );
```

This restricts β to [0, π/2], which is the physically meaningful range
(angle between field and tape normal, regardless of which side of the tape
the normal points to).

Alternatively, pass the decomposed components directly to the material model:

```cpp
real B_perp = std::abs( dot(n, b) );           // |B · n|
real B_par  = std::sqrt( norm_b*norm_b - B_perp*B_perp );  // |B - (B·n)n|
```

This avoids the angle computation entirely and makes the sign independence
explicit. The material model would take `(B_perp, B_par, T)` instead of
`(|B|, β, T)`.

## Action Items

1. Check `rho_powerlaw()` and `drho_powerlaw_dJ()` — do they use `cos(β)`
   linearly or only through `sin(β)` and `cos²(β)` / `|cos(β)|`?
2. Apply one of the fixes above to all 9 functions
3. The same check applies to `rho_piecewise()` and `drho_piecewise_dJ()`

## Related

- `todo/master_slave_swap_analysis.md` — full analysis of master/slave swap
  impact, which identified this as the only physics-sensitive location
- The `h_ts_metal` functions use `abs(dot(b, j))` (angle between field and
  current, not tape normal) — these are already sign-invariant
