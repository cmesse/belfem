# Devlog 2026-10-09 — DR-199 fixed: gauge points are fixed at setup, a late first flip is refused, a condition without `type` is refused at deck read

**Date:** 2026-10-09
**Topic:** The fix for DR-199 on Christian's instruction ("propose a fix and have the jury audit it", then "when the jury returns, you can start implementing"), with his three rulings taken in question form: a pre-fixed gauge point is a contact of the automatic-pin rule; an absent `type` is refused at deck read; the thermal factory is included
**AIs involved:** Claude Fable 5.1 (plan, implementation, gates); Codex gpt-6-astra and Grok grok-4.7 (plan jury at high, code jury at xhigh)
**Claude Confidence:** high on the serial mechanism and the tests (executed); medium on the two-rank gate (see below); thermal statically reviewed, functionally untested
**Literature References:** N/A

## Summary

A `gauge` is a bearing with a time-dependent value. Both factories now fix its points before
the free/fixed split, in the loop that fixes bearings, without disabling the automatic pins;
`impose_autopins()` takes the point's node as a contact, so only the gauge's own φ component
loses its pin; the per-timestep imposition becomes a value update and writes the scaled,
published value (a `userdefined` gauge in a declared unit was imposed unscaled — Grok's
finding in the plan round). `Bearing::impose_dirichlet` refuses a first fix of a point after
the system was initialized, on the rerouted source as well as on the direct dof, in every
build, so no caller can reach the frozen split again. Both boundary-condition factories refuse
a value-carrying condition without `type` when the deck is read; before, the source function's
member pointer stayed uninitialized and the first timestep crashed. Every value-carrying block
of the 21 shipped decks carries a `type` (section scan), and no deck uses `gauge`, so no
example changes: costheta at 1 s is bit-identical to 12a8eec.

## Juries

Plan round (Codex astra/high "revise", Grok 4.7/high "spine holds"): nine items folded in,
all re-read against the tree — printf arguments in the sketches, the Dirichlet test fixture
that the `type` refusal would have rejected, a gate that observes a nonzero time-dependent
value through the real dispatch, the unscaled gauge arm, a side record for the skip message
keyed by vertex id, the stale "worker bearings never hold dofs" sentence, both branches
guarded, thermal left functionally untested and said so, the documentation line numbers.
Code round (astra/xhigh + Grok xhigh): both approve the spine (predicate, pre-fix order,
`mValue`, the contact flag on the plane the search reads, kernel lifetime in the tests). Fixed
after the round: the refusal text overclaimed a history the predicate does not record ("was
free when" → "is free although the system is already initialized"); the ramp test now asserts
the exact 0.5 and 1.5; `<cmath>` included; the unreachable "point does not exist" check in
`impose_autopins` is an assert (the pre-fix loop refuses first); the schema's contact
sentence names the gauge point. And one correction to the record: **worker bearings never
hold dofs** — `connect_dofs_to_mesh` runs only on the master (`cl_FEM_DofMgr_DofData.cpp:3024`
inside `collect_hanging_dofs`' rank-0 branch), so `link_bearings_with_dofs` inserts nothing on
a worker; the plan round had adjudicated the opposite from the strip's rank-1 crash, which
was the mesh. The bearing comment now says master-only, agreeing with the comment above it.
Two findings are REFUTED: the input contract IS updated (reference, schema, changelog sit in
the worktree, excluded from the audited diff by the comment rule). Residuals recorded below.

## Gates

- `tests/fem`: `DirichletPin` 9/9 and `AutoPinGauge` 10/10 (19 tests), release tree. New:
  bearing late first flip refused on the direct dof and on the alias's source; pre-init fix
  then post-init value update accepted on both; the gauge point as contact (node 11 fixed, its
  component unpinned with `( gauge point 11 )`, the other pinned); constant-3 value on the dof
  after `initialize()` + `compute_boundary_conditions( 0.5 )`; ramp values at 0.25 s and
  0.75 s differ; `.bfm` reuse under a gauge deck; untyped `gauge` and `dirichlet` refused at
  construction with `'type' is required`. `make check-fast` 8/8.
- Strip reproducer (189-node variant), serial: the pre-fix binary 1879433 aborts with
  Armadillo's `Mat::operator(): index out of bounds`; the fixed binary completes setup with
  the gauge point as the contact of its component and reaches the step controller (the
  zero-load deck is then rejected there, as its controls are). The 18-node and 189-node strips
  crash on two ranks with or without a gauge, before the first step — a property of these
  hand-written meshes under the partitioner, not of this change; the two-rank gate was moved
  to costheta.
- Two ranks on a real mesh: costheta regenerated with `Mesh.SaveAll = 1` and a physical point
  at the yoke's outer apex (geometry tag 329; the deck's `gauge { nodes : 329 ; type : constant ;
  amplitude : 0 ; }`), 1 s. Fixed binary: serial 6 steps, prescribed dofs 173 → 174 (the
  point), the component's pin already left to the symmetry plane; `mpirun -np 2` the same
  counts and the terminal history to 1.65e-9 of the serial run. The pre-fix binary 1879433
  aborts on the same deck at the first assembly (prescribed 173: the point not fixed). Two
  detours recorded for the next person: costheta's geometry has no physical groups, so a lone
  `Physical Point` makes gmsh export only that point unless `Mesh.SaveAll = 1` precedes it
  (the "Key 1 not found in map" that followed was the reader on a one-element mesh, not a
  reader limitation); and the hand-written strip meshes crash on two ranks with or without a
  gauge, before the first step — not investigated.
- costheta 1 s, fixed binary against 12a8eec: same pin, terminal history and fields identical
  to the last digit (no gauge in the deck; the contact set only grows for gauge decks).

## Changes Made

`src/fem/kernel/cl_FEM_Bearing.cpp` (two refusals, the contract comment),
`src/fem/kernel/cl_FEM_PhysicalBoundaryCondition.cpp` (gauge arm imposes `mValue`),
`src/fem/maxwell/cl_MaxwellBoundaryConditionFactory.cpp` and
`src/fem/thermal/cl_ThermalBoundaryConditionFactory.cpp` (`type` refusal),
`src/fem/maxwell/cl_MaxwellFactory.cpp` (gauge pre-fix; gauge contacts and their message),
`src/fem/thermal/cl_ThermalFactory.cpp` (gauge pre-fix), `tests/fem/test_AutoPinGauge.cpp`,
`tests/fem/test_DirichletPin.cpp`, `doc/input_file_reference.md`, `doc/input_schema.yaml`,
`CHANGELOG.md`, `todo/debt_register.md` (DR-199 status), `todo/dr199_gauge_fix.md` → closed,
`todo/README.md`, `devlog/README.md`, this file.

## Residuals

- Thermal gauge: pre-fixed by the same three lines, no deck-driven thermal test fixture exists;
  statically reviewed, functionally untested. The bearing refusal is its safety net.
- The `mValue` change (scaled value imposed) has no regression that distinguishes it: the only
  non-unit scale is a `userdefined` gauge in a declared unit, which needs a plugin library; the
  constant and ramp tests pass with either implementation (Codex).
- The `dirichlet` arm of `impose_bc` still imposes the unscaled `compute( aTime )` while its
  global publishes `mValue` — the same hole for a `userdefined` Dirichlet in a declared unit,
  pre-existing, not touched (Grok).
- A gauge on a periodic alias fixes the source, but the contact flag sits on the alias's
  original; if that node lies outside the air/buffer/ferro adjacency the component keeps its
  pin beside the fixed source. No periodic gauge fixture (Grok, medium).
- The strip meshes' two-rank crash without any gauge (not investigated; not a row unless
  Christian wants one — hand-written `.msh` 2.2 with line sidesets and no physical names).
- DR-199 stays open until Christian strikes it; the status cell carries the gates.

## Files Updated

See Changes Made. Exchange: `tmp/ai_exchange/review_dr199_gauge_fix.md` (plan and code
rounds, diff), `tmp/ai_exchange/review_gauge_late_flip.md` (the diagnosis).
