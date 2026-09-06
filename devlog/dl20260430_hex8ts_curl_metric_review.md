# 2026-04-30 HEX8TS Thin-Shell Curl Metric Review

## Context

Read-only Codex review of the current `cl_EF_HEX8TS.*pp` adaptation against `todo/hex8ts_thin_shell_curl_metric.md`.

## Findings

- Correction after user clarification: for these `HEX8TS` elements the actual top-bottom element thickness is the intended thickness scale, not the connector wrap thickness. The earlier connector-thickness objection is withdrawn.
- The current geometry plumbing is now split between `Element::facet()` and `Element::reference()`: `facet()` carries the generated normal side surface used by `EF_HEX8TS::link()`, while `reference()` carries the adjacent shell/binomial facet. This makes the `aElement->facet()` geometry dependency reasonable for the current design.
- Follow-up after the latest changes: the regular-thin-shell fallback now uses the FEM element block thickness, which matches the PENTA/QUAD thin-shell pattern. The connector branch still reads `reference()->master()->block_id()`; this is only correct if that reference facet's master is guaranteed to be the layer block spanned by the connector row. In the current construction the connector map stores the same per-curve shell facet for every through-stack row, so this likely does not identify the row-specific layer thickness.
- The adaptation is still not correct yet. `update_nabla()` now builds the correct 2D Gram matrix as `mJ2 * trans(mJ2)`, but the copied components into `mInvJ` are still permuted. This breaks even simple axis-aligned mappings.
- `mDetJ = 2.0 * norm(mW) * mThickness` is not consistent with the constructed HEX reference mapping. If `mW = cross(dx/dxi, dx/deta)` and the reference thickness coordinate spans `[-1,1]`, the volume determinant should scale as `0.5 * norm(mW) * mThickness`, and the normal gradient column should scale as `2 / mThickness`, not `1 / (2 * mThickness)`.
- Separate regression observed in `mt_maxwell_h.cpp`: `h_tb()` now comments out the curl stiffness `K += C^T C rho dV`. If that was not intentional for debugging, it removes the resistive Maxwell stiffness contribution.
- Follow-up with the MATLAB cross-check: latest code now uses `mDetJ = 0.5 * norm(mW) * mThickness`, which is the correct physical-thickness scaling. The pasted MATLAB script uses `/ hz`; that only agrees for `hz = 1` and should be `* hz` if `hz` is the physical element thickness. The remaining code issue is the manual expansion of `J2' * inv(G)`: rows use wrong inverse-Gram entries and access `mInvG(0,2)` even though `mInvG` is `2 x 2`.

## Confidence

High for the pseudo-inverse component-order and determinant/normal-column scaling issues. Medium for the connector-thickness lookup concern because it depends on the final reference-facet invariant. High that the commented-out `h_tb()` stiffness is a behavioral regression if left in production.

## Files Reviewed

- `todo/hex8ts_thin_shell_curl_metric.md`
- `src/fem/interpolation/nedelec/cl_EF_HEX8TS.cpp`
- `src/fem/interpolation/nedelec/cl_EF_HEX8TS.hpp`
- `src/fem/interpolation/nedelec/cl_EF_PENTA6TS.cpp`
- `src/fem/kernel/cl_FEM_Calculator.cpp`
