# Debt Register: Scaffolded-Code Rulings — Three Rows Struck, Two Closed as Keep-by-Design, One Documented

**Date:** 2026-08-13
**Purpose:** Record Christian's rulings on DR-52, DR-53, DR-64, DR-66, DR-67, DR-68 and
place the documentation where each ruling says it belongs.
**Module:** `todo/`, `src/fem/kernel`, `src/fem/iwg`, `src/fem/interpolation/doc`

## Struck as done (Christian's instruction)

- **DR-52** (Anderson false-convergence residual): fixed and committed since `4f2c11cd`;
  the last residual, the greg3 A/B run gate, is retired with the row. Its entry in the
  register's "why it stays unstruck" table is struck alongside. The row's blocking-1.0
  flag (the last live YES in the register) is retired with it.
- **DR-64** (ferro history-scalar sign): fixed 2026-08-11; the Newton-iteration-count
  A/B gate is retired with the row.
- **DR-66** (gesvd work-buffer floor): fixed and probe-verified earlier today
  (`dl20260813_dr66_gesvd_workbuffer_floor.md`); struck per the DR-42/49 pattern — the
  `make check-fast` (`USE_TEST=ON`) gate survives in the status column, shared with
  those rows, and has still never run.

## Closed as keep-by-design (dormant capability, not dead code)

- **DR-67** — the sideset branch of `DofManager::compute_jacobian`:
  `compute_jacobian_on_sideset()` is not needed for the maxwell-thermal problem, but
  the capability could solve a future problem where surface matrices matter. Ruling:
  keep. Documented at the flag declaration in `cl_IWG.hpp`, including the one condition
  on any future writer: the sideset loop must gain the `is_active()` guard it alone
  lacks (it is the only `link_to_group( tSideSet )` site without one).
- **DR-68** — the bubble-enrichment machinery: not used, but deleting it would be a
  waste. Ruling: keep as scaffolding; it can be compile-gated (e.g. a `BELFEM_BUBBLE`
  define) if it ever gets in the way — noted, not implemented. Dormancy documented in
  `src/fem/interpolation/doc/README.md` at the `bubble/` entry, including the reason
  the experiment was dropped (bubble was the wrong enrichment space — hierarchical per
  Dular et al. 2021 — for the iron phi enrichment) so the next reader does not
  re-derive it.

## Documented position, row stays open

- **DR-53** — side-edge fusing: Christian and Prof. Sirous believe fusing the side
  edges is the mathematically more correct continuity statement, **yet the code
  converges slower with it and produces no significantly better result**. The working
  suspicion is that the fuse **overconstrains** the problem — analogous to weakly
  enforcing a B·n = 0 that the formulation already fulfills automatically at the
  boundary. Both fuse flags stay `false`, now on physics grounds and not only on the
  2026-08-09 Δt-collapse reproducer. The position is recorded in three places: the
  DR-53 row, the flag-site comment (`cl_ThinShellFactory.hpp`, `mFuseEdges` /
  `mFuseEdgesWhenHavingSideConnectors`), and the O1 section of
  `todo/side_edge_fusing_cut_aware_plan.md`. Remove-vs-keep of the flag itself remains
  the open O1 tail, so the row stays live.

## Source touched (comments only)

`cl_ThinShellFactory.hpp` (fuse-flag comment extended with the physics position) and
`cl_IWG.hpp` (dormant-capability note at `mComputeJacobianOnSideset`). Both TUs
(`cl_ThinShellFactory.cpp`, `cl_IWG.cpp`) syntax-check clean with the tree's own flags
(Blaze + MKL, `-Wall -Werror -pedantic-errors -std=gnu++17`). No behaviour changed.
