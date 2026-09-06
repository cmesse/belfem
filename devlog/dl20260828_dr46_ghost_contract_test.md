# Devlog 2026-08-28 — DR-46: the ghost 12-dof contract test is live

**Date:** 2026-08-28 (overnight autopilot session, ~02:00-04:00)
**Topic:** DR-46 — Calculator-level fixture for the `h_ghost` contract; scaffold enabled, gate run
**AIs involved:** Claude (primary), Codex + Grok (two jury rounds: plan, code)
**Claude Confidence:** high
**Codex Audit Confidence:** high (plan round: approve with corrections; code round: see below)
**Literature References:** Burman & Zunino 2006 (heterogeneous-DG penalty, via the kernel's own
comment block); Rivière/Wheeler/Girault 2001 cited by the ghost doc — see Open Questions
**Verification:** **verified — focused regression**: `make check-fast` green (8/8 suites),
`test_fem` 150/150 with `InterfaceOrientation.GhostElementContract` enabled, on
`cmake-build-debug` at working tree of 2026-08-28 ~03:00; plus 3 shuffled-order private-binary
runs and a standalone probe ladder (evidence level 2 with level-4 probes beneath it).
Valgrind on the new test: **0 errors**; the leak summary above the battery baseline traces to
production `fem::SideSet` teardown paths (`SideSetData::create_sidesets`,
`SideSet::initialize_lookup_tables`, `Group::create_calculator`) — a pre-existing teardown
leak class production never exercises repeatedly; not a fixture ownership defect (report-only)

## Summary

DR-46 ("ghost 12-dof contract test needs a Calculator-level fixture — a layer above anything
the battery builds") is **done**: `TS_TestGhostStack` in `tests/fem/support/cl_TS_TestStack.hpp`
stands up the full production `Kernel → DofManager → IWG_Maxwell(HPhi)` chain over a hand-built
two-prism ghost mesh, and the formerly `DISABLED_GhostElementContract` scaffold in
`test_InterfaceOrientation.cpp` now runs three probe-calibrated assertion groups against the
real `maxwell::h_ghost`. **No production file was touched** — the 2026-08-13 "layer above"
verdict and the 2026-08-27 "one production-side blocker (C10)" assessment were both overturned
by an executable probe: the heavy route through real production constructors needs no seams.

The session followed probe → plan → jury (Codex+Grok, blind) → implement → gate → jury on the
diff. Yesterday's feasibility round (`tmp/ai_exchange/review_dr46_ghost_fixture.md`, never
distilled until now) is folded into this devlog; both of tonight's rounds live in
`tmp/ai_exchange/review_dr46_fixture_plan.md` and `review_dr46_ghost_code.md`.

Jury outcomes: plan round — Codex approve-with-corrections (drop the redundant re-link,
rewrite the scaffold prose; both adopted pre-implementation). Code round — Codex
request-changes with four findings (gTbulk RAII guard for ctor-throw safety, delete the
owning fixture's copy ctor, `ASSERT_GT` before the ratio divisions, `aNorm`→`tNorm`), all
four adopted and all gates re-run green. Grok's headless wrapper died at max turns in both
rounds; its recovered partials contributed one confirmed find (the `ThinShell` ctor hides its
sidesets — documented in the fixture) and four claims refuted by measurement or source read
(teardown inversion, geometric-artifact rho scaling, "no ghost element instantiated",
"MeshChecker lacks PENTA6TS").

## Key Findings

### The fixture route (executable evidence, standalone probe against `libbelfem.so`)

- The lightweight route (light ctors + protected-member pokes) is dead — yesterday's round:
  `const` zero element count in `Group`, no `set_master`/`set_slave` on `fem::Element`,
  `link_edge_dofs` aborts on `DofManagerBase`. The HEAVY route works instead: build the
  post-ThinShellFactory mesh state by hand and let the production chain consume it.
- What the mesh must carry (all found by probe iteration): two PENTA6TS layer blocks
  (thickness set, `DomainType::ThinShell`), shared middle node row with **duplicated**
  interface edges (the `hasDuplicates` state of `cl_ThinShellFactory.cpp:142`), a ghost TRI3
  facet (master = lower top face, slave = upper bottom face, orientation 1), a ghost sideset
  (`DomainType::Ghost`), a tape sideset with the mid-surface reference facet, and a
  `mesh::ThinShell` container — `BlockData::link_thin_shell_facets_serial` reads it to give
  the layer fem elements their facet, which `EF_PENTA6TS::link` dereferences.
- Order matters: `finalize()` → `MeshChecker` (must run before edges exist; sets the flag the
  Kernel ctor checks) → hand-built edges + `finalize_edges()` → `Kernel`.
- `element_rho` is a plain mesh field (`mFields.NonDof`), seeded directly; `gTbulk` must be
  non-NaN or `h_ghost` demands a nodal `"T"` field (the dead-T branch found 2026-08-27
  executes in any run without a bulk temperature).

### The contract itself (probe numbers the test's thresholds are calibrated against)

- K is 12×12; `2m == n` (facet vs layer nedelec dofs) is silently assumed by the kernel's
  `Dm`/`Ds` assembly and now pinned by the test.
- "Block sparsity" as the scaffold imagined it is WRONG — no block is structurally zero.
  The truth is **interface locality**: far-from-interface columns (master bottom, slave top)
  carry only the `rho_harm`-scaled gradient coupling (~7e-8 at rho = 1e-8), near columns the
  k_reg-dominated penalty (~3e-3). Ten-folding rho scales far columns by exactly 10.0000 and
  moves near columns ≤ 0.14% — the test asserts the mechanism, not the numbers.
- The annihilation patch test holds at machine zero (max|K·q| ≈ 2e-19 vs negative control
  1.4e-3) with q built as edge line-integrals in the edge object's own node order. This
  resolves the 2026-08-27 round's open physics caveat by execution.
- Edge flips: the only LEGAL flip reverses local edge k in **all rows of both prisms** —
  ThinShellFactory extrudes every row from the same mid-surface edges. A partial flip that
  still satisfies the DEBUG bitset assert (`cl_FEM_Element.cpp:810-824`) but breaks
  co-location leaves an rho-scale annihilation residue: that assert has a blind spot
  (report-only observation). The test runs canonical and one legal flip variant.

### Latent production defect found (NOT fixed — tests-only session): uninitialized `Kernel::mController`

`Kernel::mController` is not in the constructor's initializer list (`cl_FEM_Kernel.cpp:46-52`);
only `Kernel::set_controller` (`:1129`) ever assigns it. `IWG_Maxwell::link_to_group`
unconditionally does `controller()->thermal_kernel()` (`cl_IWG_Maxwell.cpp:266`), and the
`BELFEM_ASSERT` in `controller()` cannot catch garbage-non-null. Reproduced tonight: the
fixture's `link_to_group` call passed solo and SEGV'd after 21 sibling tests (heap-layout
dependent). Any kernel that never receives a controller — unit contexts today, and plausibly
future library consumers — walks this path. Same family as the gantry `mFunMKF` null jump
(dl20260827_gantry_thermal_nulljump). Filed as **DR-129**; the fixture documents and skips
`link_to_group` (the contract test calls `h_ghost` directly and needs nothing it sets).

### The symmetry discrepancy is now measured on real edge functions (adjudication owed)

`h_ghost` assembles an **exactly symmetric** K (max|K−Kᵀ| = 0.0 on the fixture — same
accumulation order both triangles; 2026-08-27's random-matrix check gave 1.8e-15), while
`src/fem/maxwell/doc/ghost_penalty_stabilization.md:44-46` states the sign pattern is
nonsymmetric (nonsymmetric Nitsche, Rivière/Wheeler/Girault 2001). Doc and code cannot both
be right; whether the kernel was accidentally symmetrized or the doc is stale is a
**formulation question for Christian** — the contract test deliberately asserts nothing about
symmetry, and `IWG_Maxwell`'s `SymmetryMode::Unsymmetric` is a no-op for this kernel either
way. Filed as **DR-130** `[RULING]`.

## Changes Made

- `tests/fem/support/cl_TS_TestStack.hpp` — `TS_TestGhostStack` fixture (+ the maxwell/mesh
  includes it needs); saves/restores `gTbulk`; `set_rho()` helper; legal-flip parameter.
- `tests/fem/test_InterfaceOrientation.cpp` — scaffold comment rewritten (stale `:1805-1937`
  cites dropped, per Codex), `DISABLED_GhostElementContract` → `GhostElementContract` with the
  three assertion groups, serial-only guard, two flip variants.
- `tests/fem/CMakeLists.txt` — `fem/maxwell` + `fem/maxwell/matrices` include dirs inside the
  existing `if( USE_MAXWELL )` guard.
- `todo/debt_register.md` — DR-46 closed (gate run), DR-129 and DR-130 filed.

## Open Questions

- **DR-130 (Christian):** symmetric kernel vs nonsymmetric doc — which is intended? Tonight's
  fixture makes the follow-up cheap: one assertion pins whichever answer is ruled.
- **DR-129:** one-line hardening (`mController = nullptr` default + a guard or documented
  contract in `link_to_group`) — src change, needs its own approved round.
- Grok's wrapper died at max turns in BOTH tonight's rounds (second and third consecutive
  failures counting 2026-08-27); its partial material was recovered from the log each time and
  verified, but the wrapper's turn budget needs raising (`GROK_MAX_TURNS`).

## Files Updated

- tests/fem/support/cl_TS_TestStack.hpp
- tests/fem/test_InterfaceOrientation.cpp
- tests/fem/CMakeLists.txt
- todo/debt_register.md
- devlog/README.md (index entry)
