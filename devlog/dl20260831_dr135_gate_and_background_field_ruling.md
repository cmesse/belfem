# Devlog 2026-08-31 — DR-135 gated, and the background-field cluster closed

**Date:** 2026-08-31
**Topic:** Run the DR-135 gate on an existing mesh; close DR-135, DR-136 and DR-137
**Module:** fem/kernel (SideSet BC, Controller), fem/maxwell (factory), doc
**AIs involved:** Claude
**Verification:** executable — five sandbox runs on the helix mesh, two on a disk-deck copy.
No source changed; one documentation file and the register did.

## Summary

DR-135's gate was believed to need a mesh that did not exist. It did not: the trigger is
reachable on the helix mesh already in the tree, and the run both fires the reroute and
discriminates the defect. DR-135, DR-136 and DR-137 are struck and archived.

## What the gate actually needed, and three wrong turns before it

The row's gate reads "a symmetry sideset crossing the periodic target face". Every clause
in that sentence is load-bearing, and I got each one wrong in turn:

1. **Not a `dirichlet` BC.** A `dirichlet { }` block in `boundary conditions` never reaches
   `SideSet::impose_dirichlet` on this mesh. The trigger is a *topology* domain: `AirSymmetry`,
   `BufferSymmetry` or `FerroSymmetry`, which `create_magnetic_kernel` routes through
   `impose_dirichlet( 0.0 )` (`cl_MaxwellFactory.cpp:851-856`) — "phi symmetry corresponds to
   a constant phi on the symmetry plane".
2. **"Crossing" is not loose phrasing.** `PeriodicityFactory::tag_periodic_sidesets`
   (`cl_Mesh_PeriodicityFactory.cpp:1754-1775`) OVERWRITES the domain type of any sideset whose
   nodes ALL lie on a periodic plane with `DomainType::Periodic`, which reaches neither the
   selection switch nor the impose_dirichlet branch. A symmetry sideset ON the target face is
   therefore silently swallowed and can never trigger this row. It must cross: some nodes on
   the plane, most off it.
3. **The current drive cannot coexist with the trigger.** An air-symmetry surface on the outer
   boundary takes the cohomology to 0 generators, and a current condition needs one, so the run
   aborts on the generator budget. Christian's redirect solved it: drive the model with a
   `background` field instead. It is imposed weakly through its own stiffness matrix and spends
   no generator, so the two coexist.

## The gate run

Deck: the helix periodic model, `topology { air symmetry { sidesets : 133 } }` (one of the four
outer-cylinder bands, z = 0..5, so it crosses both periodic planes), driven by
`background : bgfield { sidesets : 137, 141, 145 ; direction : 0,0,1 ; ramp 0.1 T over 1 s }`.
No transport current, no geometry authored.

```
SideSet 133: 16 node(s) are condensed onto sources outside the free set; pinning those source dofs instead
  unknowns 98948 / prescribed 143 ( was 2 ) / condensed 26978
  17/17 timesteps, 0 rejections, exit 0, full 100 ms, dt grown to the cap
```

**The run discriminates rather than merely executes.** Cross-checked against the mesh's own
`/hanging/nodes/` table in the `.bfm`: sideset 133 carries 143 nodes, 17 of them condensed,
whose sources are **16 outside the sideset and 1 inside**. The one interior source was already
fixed when the loop reached it, which is why the log counts 16 first flips. The pre-fix path
would have called `fix()` on 16 eliminated dofs and lost those constraints silently.

The imposed field tracked the analytic ramp to six figures at every step — B/mu0 with the -1 ms
offset, 159.155 / 238.732 / 358.099 / ... / 8037.325 A/m — which incidentally exercises the
Tesla-family conversion end to end under adaptive stepping. (Those values print in the `I_max`
slot with an `A` suffix; that mislabelling is a separate defect already being fixed elsewhere.)

What the run does NOT show, kept explicit so the strike does not overclaim: no field-level
equivalence to a direct source-face pin, the two deliberate abort branches (non-unit weight,
post-freeze first flip) were never exercised, and it was single-rank so the row's MPI residual
stands. Binary `build/bin/belfem`, 2026-08-30 19:37:19, debug.

## Why the equivalence clause was waived rather than deferred

The row's third observable is "identical fields to a direct source-face pin". The 16 sources are
the z = 0 rim of the ADJACENT band — a line of nodes. A deck can pin sidesets or bearing points,
not an arbitrary node set, so the twin is not expressible on this mesh, and the clause is
unconstructible rather than merely unrun. That is the DR-90 shape ("gate found INAPPLICABLE
rather than blocked"), one rung stronger, because here the trigger fired.

## DR-136: a design ruling, and where it now lives

Christian's ruling: **the background field is imposed weakly by design; the strong variant is
deliberately dropped.** Costed, it buys a handful of prescribed dofs for two real costs.
Removing the automatic gauge pins is the first — with phi prescribed over the whole boundary the
potential is already determined, so an autopin fixing an interior node to 0.0 becomes an
inconsistent extra constraint rather than a redundant one. The second is that the source value
is spatial, so each condensed target needs the value at its source's position under the periodic
transform: a gradient derivation over the phi domain, not a copied constant.

The reasoning was written into `doc/input_file_reference.md` §13 in the same session, together
with the standing veto that the dead branch must **not** be repaired by routing it through
`fem::pin_dirichlet_dof` (that helper pins a condensed dof's source to one constant, which is
right for a constant Dirichlet sideset and wrong for a spatial value). Documents state facts
inline, so it could not be left pointing at an archived register row.

## DR-137: struck as accepted exposure, and a reproducer that failed

Re-confirmed statically: the ceiling is checked before `watchdog_magnetic` with an early
`return` between them, so on any deck where `max iterations` <= `watchdog window` the watchdog
cannot fire. **The row's scope was narrower than the code** — the same ordering is in
`iterate_magnetic()` (ceiling `cl_FEM_Controller.cpp:2242`, watchdog `:2271`), which every
magnetic-only deck uses, not just `iterate_coupled()` (`:1846` / `:1914`).

Struck as accepted exposure because both paths end the attempt through `reset_timestep()`: the
ceiling winning costs burnt iterations, not a wrong answer.

**The attempted reproducer failed, and the failure is worth recording so it is not retried the
same way.** Stiffening the power law to n = 40 and driving the tolerance to 1e-16 does not
create the watchdog's precondition. The residual falls monotonically to a floor near -146 dB and
improves by a hair on EVERY iterate, so `watchdog_magnetic` takes its `mEpsilon < mBestEpsilon`
branch each time and the no-new-best window never elapses. Run at `watchdog window` 30 and 10,
the two logs were byte-identical over every Magnetic/rejected/Watchdog line — 44 rejections
each, 0 watchdog fires, both ending at the floor-retry limit. The watchdog was never *eligible*,
so nothing was proven either way. A discriminating deck needs a genuine plateau (best residual
frozen for the window, relaxation pinned or shrinking so the 2026-08-21 anti-cascade spare at
`:855` does not excuse it) — the omega-sawtooth limit cycle, not a residual floor.

The net below it was observed to hold: the ceiling cut all 44 attempts and the run stopped
cleanly on the floor-retry limit with a diagnostic message.

## By-catch, not filed

Christian ruled that no new DR rows are being opened, so these are recorded here only:

- **A Maxwell `dirichlet` BC naming a sideset the DofManager does not own is accepted and
  silently discarded.** `SideSetData::sideset()` (`cl_FEM_DofMgr_SideSetData.hpp:151-162`)
  returns `mEmptySideset` with the requested id stamped on it instead of aborting, so both the
  factory pre-fix loop and the per-timestep `impose_bc` iterate zero nodes. A deck naming
  sideset 9999 runs clean to convergence. Removing the `sideset` key aborts and an unknown
  section type aborts, so the section is genuinely parsed — only its target is thrown away.
  No shipped deck uses this BC type (`grep -rn dirichlet examples/` is empty).
- **Deck id spaces differ silently.** The helix deck's `periodic { target : 294, 282, 284 }`
  are gmsh surface tags; that mesh's BELFEM sidesets are `1-5, 29..145, 210, 211, 212`. Both
  spellings are accepted and only one resolves.
- **A symmetry domain declared on a periodic plane is silently retagged**, as above. Same
  failure mode: the deck says something and the code discards it without a word.

## Files Updated

- `doc/input_file_reference.md` — §13, why the background field is not imposed strongly
- `todo/debt_register.md` — DR-135, DR-136, DR-137 struck; `[P]` 9, `[W]` 4
- `todo/debt_register_closed.md` — the three rows archived (136 closed, 13 live)
- `todo/README.md` — register counts

## Status

DR-135 struck on evidence with one clause waived as unconstructible; DR-136 on a design ruling;
DR-137 as accepted exposure, never executed. **Struck is not verified** applies to all three.
No source was changed in this session.
