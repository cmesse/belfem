# DR-135 Fixed: Sideset Dirichlet Pins Survive Periodic Condensation

**Date:** 2026-08-28
**Purpose:** Fix DR-135 ( `SideSet::impose_dirichlet` silently loses pins on
periodically-condensed dofs ) through a full three-round jury; file DR-136.
**Module:** fem/kernel ( SideSet, DofManagerBase/DofManager ), fem/maxwell ( factory )

## The defect and the fix

`SideSet::impose_dirichlet` called `fix( aValue )` unconditionally; on a dof
condensed onto a periodic source the flag and value were historically
ineffective for elimination — the T-matrix reconstructs the alias from its
source, so the constraint was silently lost ( the DR-126 bearing mechanism ).

Landed ( Christian: "Let's fix 135 right away" ):

- **`fem::pin_dirichlet_dof`** ( declared `cl_FEM_SideSet.hpp`, defined
  `cl_FEM_SideSet.cpp` ) — the single shared shape test: a hanging dof with
  exactly one source and unit weight ( `BELFEM_EPSILON` tolerance ) gets its
  SOURCE fixed to value/weight; non-unit single-source weights, broken
  sources, and a first free-to-fixed flip after `initialize()` are loud
  `BELFEM_ERROR`s. Everything else returns to the caller's historical path.
- **`SideSet::impose_dirichlet`** calls the helper for non-duplicate nodes,
  announces reroutes once per sideset ( first flips only, so the
  per-timestep path stays quiet after step one ).
- **The Maxwell factory pre-fix loop** ( `cl_MaxwellFactory.cpp` ) calls the
  SAME helper under its existing non-duplicate filter, so the source set
  pre-pinned before graph freeze is exactly the set the per-timestep path
  reroutes.
- **`DofManagerBase::is_initialized()`** ( virtual, base false ) exposes
  `DofManager::mInitializedFlag` for the post-freeze guard.
- Tolerances unified on Christian's mid-round direction: `BELFEM_EPSILON`
  everywhere, including redefining the bearing's `gBearingUnitWeightTol`
  from 1e-12 — safe because pair weights are stored literals.

Deliberate behaviour change: a 1-source non-unit condensation on a Dirichlet
sideset now aborts loudly at both sites instead of silently fixing the
eliminated target.

## The jury round ( plan+audit → code+audit → reconcile+audit )

- **Codex round 1 ( request changes )** found the two load-bearing gaps in the
  original plan: the Maxwell factory pre-fixed Dirichlet nodes through its own
  unconditional `fix(0.0)` loop, so the first reroute through the repaired
  SideSet would have come per-timestep, AFTER graph freeze; and single-source
  unit-weight does not prove periodic provenance ( a periodic target whose
  master is itself hanging is multi-source ).
- **Round 2 split:** Codex approve; **Grok request changes** — correctly
  refusing the landed factory fork: a second, weaker copy of the shape test
  with a different node filter meant the factory pre-pin set ≠ the SideSet
  reroute set, and a duplicate sideset node would have turned the new guard
  into a production abort. Its thin-shell point was physical: duplicate-pair
  condensations can carry an inhomogeneous relation ( cut: phi_dup =
  phi_orig + I ) that a plain source pin would violate.
- **Round 3 unanimous approve** on the reconciliation ( shared helper,
  duplicate gate at both sites ).

Vendor infrastructure cost one evening's friction, all root-caused: Grok's
sandbox died pre-API on `/run/podman` permissions ( fixed by Christian ), and
two later Grok runs crashed with phantom bash syntax errors because
`ask_grok.sh` was being edited by a parallel session WHILE my instances
executed it — bash reads scripts incrementally, so a moving file shifts under
the byte offset. Worked around by running a scratchpad snapshot of the
wrapper with the repo path pinned.

## Corrections to the record

- Claude's DR-136 row initially claimed `"background dirichlet"` is
  deck-reachable, citing post-creation code. **Wrong** ( Grok round 3,
  confirmed by direct read ): the creation switch hard-refuses the type
  ( documented since 2026-08-11 ). Retracted in the row. Since Christian
  plans a background-dirichlet deck "soon, but not now", DR-136 is now the
  checklist for finishing that feature — including Grok's veto on routing it
  through `pin_dirichlet_dof` ( its value is spatial; the source needs the
  transform-aware value, not a copy ).
- Both auditors' round-1/2 "dead path" phrasing and Claude's counter were
  each half-wrong; the register row now states the mechanism precisely.

## Status and residuals

**Reviewed, not verified.** `g++ -fsyntax-only` green on all touched TUs with
module flags; no executable gate ran. The row keeps Grok's residuals verbatim:
multi-source periodic targets stay silently lost ( indistinguishable from
refinement ), thin-shell 1-source duplicates keep the historically-ineffective
pin, per-timestep value transport of an off-sideset source is not MPI-proven,
and the helper cannot itself enforce the non-duplicate contract. The owed gate
is a probe deck with a symmetry sideset crossing a periodic target face — it
does not exist yet. Register: DR-135 retagged `[RUN][P]`; DR-136 filed
`[CODE][W]`. Uncommitted.
