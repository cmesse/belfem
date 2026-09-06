# Devlog 2026-08-09 — SideLayer Viz-Decoupling Regression Round + Δt-Collapse Root Cause

**Date:** 2026-08-09 (evening; follows dl20260809_exodus_empty_blocks_jury.md)
**Topic:** End-of-day jury round on the working-tree diff (ExodusWriter fix + SideLayer visualization decoupling bundle); Δt collapse at t = 2.0825 ms root-caused
**AIs involved:** Claude (pre-registration + verification), Codex + Grok (blind jury)
**Claude Confidence:** high (collapse mechanism, reproducer-backed); medium-high (latent findings)
**Codex Audit Confidence:** medium ("would not sign off as-is")
**Grok Audit Confidence:** high (writer/flags/postproc), medium-high (SideLayer/seam-T)
**Literature References:** none required (infrastructure/constraint wiring)
**Verification:** END-TO-END REPRODUCER for the headline item — Christian's A/B: `mFuseEdgesWhenHavingSideConnectors = false` passes the t≈2.08 ms cliff both warm (.bfm reload) and cold (fresh create); all other findings source-trace tier

## Summary

The frozen jury round targeted the full working-tree diff: Claude's ExodusWriter
consistent-block-list fix plus Christian's SideLayer visualization decoupling
(decoupled inner/outer wall sheets, flag-slot API, seam-T rework,
`MaxwellPostprocessorType::SideConnector`, `ProtoMesh::normalize_edge_hangs`).
Mid-round, the sidecoatings run collapsed its timestep at t = 2.0825 ms with a
state-independent magnetic residual floor of ~2e-6 (line search rejecting the
Newton direction at every relax down to 1e-4, thermal at machine floor) —
diagnosed as a constant constraint inconsistency, and confirmed by Christian's
A/B the same evening.

## Key Findings

- **P0 (RESOLVED by switch, default decision open): Δt collapse = outer-edge tape
  fusing fighting the decoupled sheets.** With the new inner wall sheets hanging
  1:1 on tape edges, the fuse branch (cl_ThinShellFactory.hpp:297-309) still
  hangs the OUTER wall edges on tape temp-edge nodes; at rim stations where the
  cut-composed λ branches differ above/below (90/348 stacks per the 2026-08-06
  authority census) the two constraint families clash — a constant ~2e-6 magnetic
  residual no step length can remove, so the watchdog collapses Δt to the 100 ns
  floor. Fuse-off passes the cliff warm AND cold (reproducer tier). Consistent
  with the side-edge-fusing campaign's standing "ship = free rims" conclusion.
  Default (false vs. removal) routed to Christian / Prof. Sirous.
- **P1 latent: seam-T station permutation.** The reworked seam-T read uses
  recovery-facet nodes in master canonical order while the bilinear wall
  interpolation assumes wall (xi,zeta) station order
  (cl_FEM_Calculator.cpp:615-621); the code itself documents the master/slave
  face-cycle offset as non-constant (cl_ThinShellFactory.cpp:632-637). No
  alignment proof exists; the "ORIGINAL-NORMALIZED" comment (:435-438) is false
  as written. Not the collapse driver (static, consistent). Gate before thermal
  wall production runs: one-element probe comparing facet node ids vs wall
  tTapeFace ids for both connector signs.
- **P2: `normalize_edge_hangs` ignores the wall slot rule on reload** (Codex's
  find): slot reconstruction takes the FIRST twin for HEX8TB slots 2/3
  (cl_ProtoMesh.cpp:1727-1748) while the hang normalizer re-points every 1:1
  edge hang to the LAST twin (:1814-1851). Real contradiction; runtime scope
  narrow (fires only where the top boundary layer carries a duplicate sheet —
  tonight's warm start converged).
- **P2 latent (cluster):** collect_blocks filter inherits the
  capacity-vs-fill `number_of_elements()` semantics (round-1 item);
  SideConnector missing from the postproc element-field branches
  (cl_MaxwellPostprocessor.cpp:314-336 — inert while the factory passes false);
  misleading "not in the gather" comment in the SideLayer d==1 slot resolution
  (construction-sound, Grok+Claude both traced it); outer wall nodes
  intentionally outside duplicate lists — needs a design note.
- **Watch items:** opposite-oriented 1:1 edge twins under normalize (unproven
  either way); NEW — wall viz duplicates (`tDupI`) now populate tape nodes'
  duplicate containers, which dof machinery (cl_FEM_DofMgr_BlockData.cpp:310,
  :343), homology (cl_InterfaceProcessor.cpp:78), and the distributor walk —
  not the collapse driver, but consumers assuming "duplicates = cut sheets" now
  see viz nodes; needs a census probe.
- **Clean:** ExodusWriter unification held under both auditors; flag-slot API
  extension complete (no missed overrides); SideConnector postproc
  wiring/linkage verified end-to-end by Grok; Q1/Q3/Q5/Q6/Q7/Q9 closed.
- **P3:** coding_philosophy flags example uses slot 0 for a function-local dedup
  while its own convention reserves slot 0 for cross-function protocols (3/3);
  progressbar sideset count skew in the writer.

## Changes Made / Proposed

- No source edits this round (Christian flipped `mFuseEdgesWhenHavingSideConnectors`
  himself at 20:47 and verified the A/B).
- Proposed, pending Christian: seam-T wall-ordered read (or order proof + comment
  fix); fuse default ruling; slot-aware hang normalization (or documented
  restriction); doc example slot fix; SideConnector element-field branch when
  element fields ship.

## Open Questions

- Fuse default: `false` in source vs. removing the branch — physics ruling
  (free rims), Christian + Prof. Sirous.
- Seam-T order probe result (gate before thermal-wall production).
- Duplicate-list pollution census; Q8 orientation census.

## Addendum (same night) — seam-T copy root-caused, three voices + numeric probe

Christian reported the wall seam-temperature copy from the adjacent PENTA6TS
elements "does not work"; focused Codex+Grok consult on thread
`tmp/ai_exchange/seam_temperature_copy.md`. Outcome (details + reconciliation
table in the thread):

- **The read-route equivalence holds** — `Calculator::q()` reads the mesh field
  directly (cl_FEM_Calculator.cpp:2815-2822), the thermal solver writes every
  accepted iterate back into the field, and the layer rim nodes ARE thermal dof
  nodes (H1 refuted 3/3). Both of Christian's routes read the same store.
- **P1 (dominant in MPI runs, = hazard P7b confirmed live):** the viz copy is
  rank-local and the save-time field gather collects only dof-carrying entities.
  Numeric probe on e-s.00004 (9 ranks): wall 27 has zero master-owned nodes →
  100% frozen at 77.0; wall 29's owner-0 nodes are 100% copied, all others
  frozen. Proposed fix: move the VIZ copy to the master-side SideConnector
  postprocessor; keep the physics Tseam read in the frame prep.
- **P1 (always on):** facet-order permutation proven 3/3 — `set_master(...,true)`
  copies PENTA canonical face order (face 2 = {0,3,5,2}, up-leg first;
  HEX8TB face 2 = {2,3,7,6} ≠ station {3,2,6,7}); `compute_orientation()` is
  computed at build and never used by the read. Fix: restore the station-ordered
  element route `T( node(tTapeFace[k])->original()->index() )`
  (cl_FEM_Calculator.cpp:462-463).
- **New P1-candidate (Codex side-find, Claude-verified):**
  `compute_hanging_dofs()` is DEAD in the vector solve path — the call at
  SolverData:2299 sits behind an always-aborting `default:` after all case
  breaks; dead since 22432b05 (2025-03-07). Only the multi-RHS branch calls it.
  The HEX8TB walls are the first heavy hanging-EDGE consumers; consequence trace
  (who refreshes hanging-dof field values for wall q() reads?) routed to
  Christian.
- P2: Grok side-find — Picard branch skips the mFieldValues refresh before the
  residual multiply (Newton refreshes at :2223); known Anderson/residual
  territory.

### Fix round (same night, Christian's approval, incl. side-finds)

All three edits landed uncommitted, `g++ -fsyntax-only` green with real
per-target flags; build + serial/MPI rerun is Christian's gate:

- `cl_FEM_Calculator.cpp`: Tseam read restored to the station-ordered
  `original()` element route; rank-local viz writes removed from assembly;
  contract comment corrected; unused `tFacet` local dropped.
- `cl_MaxwellPostprocessor.{hpp,cpp}`: new master-only
  `copy_seam_temperatures()` (guards: rank 0, `field_exists("T")`) iterating
  Left/RightCoating blocks, station-ordered copy onto both wall faces; wired
  as the SideConnector case in `run()`.
- `cl_FEM_DofMgr_SolverData.cpp`: `compute_hanging_dofs` moved out of the
  algorithm switch — restores pre-22432b05 behavior (hanging dofs follow their
  sources after every accepted vector-path update; soft-fail returns still skip
  it). The Picard pre-update residual was NOT changed (standing 2026-08-08
  ruling, Messe 2023 §4); a documenting comment now marks it deliberate.
- Follow-up after Christian's first run (T copy confirmed working): the
  recovered wall H spiked to ~15 kT near the YBCO layer — nodal-parametric
  edge-function evaluation on extreme-aspect slivers, not physics. On
  Christian's ruling the copy was generalized: `copy_seam_fields()` now copies
  T, Hx/Hy/Hz, Bx/By/Bz from the tape rim originals onto both wall faces
  (physical basis: thin µᵣ=1 strip → seam field = wall field to sub-mT;
  B = µ₀H adds nothing); J = C·q stays native. Convention documented in
  `side_connector_wall_element.md` §5.1; the underlying recovery defect is
  hidden, not fixed — future session if computed wall-H is ever needed.

## Files Updated

- tmp/ai_exchange/review_sideconnector_viz_exodus.md (regression round record)
- tmp/ai_exchange/seam_temperature_copy.md (seam-T investigation record)
- devlog/dl20260809_sideconnector_viz_regression_round.md (this file)
- devlog/README.md (index)
