# Removal of Side Connectors from BELFEM

**Date:** 2026-06-05
**Author:** Claude (documentation update at user request)
**Branch:** `periodic_new`
**Scope:** Documentation only. Source-code removal of the side-connector construction is tracked separately.

---

## Decision

The **side-connector** construction — the wrap of `HEX8TS` elements added on the lateral edge of a thin-shell tape stack to model the surround-plated copper that shorts the top and bottom conducting layers — has been **removed from BELFEM because it turned out to be unphysical**.

The deciding problem is the one diagnosed at length in the (now deleted) `shell_connector_coupling.md` and `hex8ts_thin_shell_curl_metric.md`: the `HEX8TS` wrap cannot represent the binormal-H field at the PENTA6TS ↔ HEX8TS fold. The connector is not "wetted" (it has no `compute_bn` analogue feeding it transport information), so `H_binormal` in the wrap is left as a free parameter that the solver drives toward zero in low-resistivity copper. In the pure-copper test this produced an ≈ 2-orders-of-magnitude discontinuity in `H_binormal` across the fold and a corresponding suppression of `J/J_c` along the YBCO side curve. Neither the geometric-hanging path nor the `h_penalty` weak-coupling path reached a usable, convergent fix.

Rather than carry a construction that produces unphysical side-curve currents, the decision is to remove it.

## Documentation changes made in this session

**Deleted (files entirely about side connectors):**

- `side_connector_review.md` (repository root) — the review-and-cleanup plan for the wrap.
- `src/fem/maxwell/doc/side_connector_modeling.md` — modeling stance and CutFactory-on-envelope architecture.
- `src/fem/maxwell/doc/shell_connector_coupling.md` — binormal-H discontinuity diagnosis and the two proposed fix paths.
- `todo/geometric_hanging_connector_edges.md` — geometric-hanging design derivation.
- `todo/hex8ts_audit_synthesis_2026-04-17.md` — consolidated HEX8TS audit synthesis.
- `todo/hex8ts_thin_shell_curl_metric.md` — HEX8TS curl-metric report.
- `todo/hex8ts_review_20260417.md` — HEX8TS edge-function review.

**Edited (living reference docs and indexes that mentioned the wrap among other content):**

- `src/fem/interpolation/doc/nedelec_thinshell.md` — retitled to *Thin-Shell Nedelec Elements*; removed the `HEX8TS` section, overview row, and implementation-map rows; left a short historical/removal note. `QUAD4TS` and `PENTA6TS` content unchanged.
- `src/fem/interpolation/doc/nedelec.md` — removed the `HEX8TS` rows from the element and DOF-count tables; reworded the "beam variants" references and the unit-circulation example that pointed at `EF_HEX8TS`.
- `src/fem/interpolation/doc/README.md` — removed the `HEX8TS` bullets and quick-reference row; dropped "and beam" from the `nedelec_thinshell.md` summary.
- `src/fem/interpolation/doc/interpolation_usage_guide.md` — removed the `HEX8TS` DOF-table row and the stale-`det_J()` gotcha's `HEX8TS` example (kept the gotcha as a general note for reduced/curved elements).
- `src/fem/maxwell/doc/README.md` — removed the `side_connector_modeling.md` and `shell_connector_coupling.md` index entries.
- `todo/README.md` — removed the `geometric_hanging_connector_edges.md` index entry.
- `todo/edge_function_quadratic_shells.md` — dropped `HEX8TS` from the linear-thin-shell list.
- `todo/hdf5_writer_repair.md` — changed a side-connector test-case suggestion to a thin-shell one.
- `todo/periodic_bc_fix_plan.md` — removed the now-moot "side-connector periodic handling" future-work item.

**Deliberately left unchanged:**

- All existing `devlog/dl*.md` entries — these are the immutable dated record of the audit history that led to this decision. They still reference `HEX8TS`, side connectors, and the deleted docs; that is correct for a historical log.
- `todo/ai_exchange.md` — append-only AI-to-AI exchange log; treated as historical record.
- Git-lineage references to the `sideconnectors` branch in `todo/periodic_bc_fix_plan.md` and `todo/README.md` — accurate statements about how `periodic_new` was derived.

## Follow-up (not done here)

- Source-code removal of `cl_Element_HEX8TS.*`, `cl_EF_HEX8TS.*`, the `ThinShellFactory` side-connector path (`create_side_connectors` and helpers), the `InterfaceTsConnector` domain type and `h_penalty` kernel, `LeftCoating`/`RightCoating` block dispatch, and the postprocessor connector-block plumbing. This devlog covers documentation only.
