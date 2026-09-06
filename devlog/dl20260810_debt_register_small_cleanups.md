# Devlog 2026-08-10 — Debt-Register Cleanup Batch: DR-49, DR-47, DR-24

**Date:** 2026-08-10
**Topic:** Three small register items cleared in one pass — spline tests compiled out by a
wrong ifdef, a test disabled by a blocker that no longer exists, and an unguarded
prefix-order assumption in `set_currents`
**AIs involved:** Claude
**Claude Confidence:** high on DR-49 and DR-24 (both statically verifiable); medium on DR-47
— the flip is justified, but the test has never executed, so its expectation is unproven
**Verification:** source-level only — **nothing here has been compiled or run.** Three source
files changed; the build and `make check-fast` are Christian's.
**Register rows:** DR-49 (fixed), DR-47 (follow-up done), DR-24 (residual guarded), DR-59
(gather half deleted after a jury round — §6; packing half still open in live code)

## Summary

Three items, chosen because none of them needs a run to *do*, only to confirm. One turned out
bigger than its row suggested, one smaller, and one had to be scoped back to what is honestly
checkable.

## 1. DR-49 — spline tests were gated on a dependency the code does not use

Every spline construction and evaluation test in `tests/math/test_Spline.cpp` sat behind
`#ifdef BELFEM_SUITESPARSE`. SuiteSparse is **off by default** (legacy, not BSD-3 clean), so
the whole §3–§7 block compiled out of the normal build: **spline coverage in the fast gate was
zero**, and DR-42's planned `aCol` tests would have been invisible too.

The gate was also wrong on its own terms. Production `Spline::update_data` solves with
**SuperLU**, not UMFPACK (`cl_Spline.cpp:427-429`, with a comment saying exactly why:
"SuiteSparse is legacy here and is not BSD-3 clean, so it is off by default and cannot be
relied on"). The tests were gated on a library the code under test had stopped using.

Re-gated on `BELFEM_SUPERLU` at all four sites (62/76, 262/712, 731/766), which is ON by
default (`CMakeLists.txt:69`). Section headers corrected from "require UMFPACK", and the
§3–§7 header now carries the reason the old gate was wrong so it does not come back.
`tests/sparse/test_Solver.cpp` keeps its SuiteSparse gates — those really are UMFPACK solver
tests. `tests/doc/tests_09_spline.md` updated to match.

**Expected effect: spline coverage goes from nothing to the full §3–§7 set.** That is also the
risk — these tests have not run in the default configuration, so the first `make check-fast`
after this change is the real test of the change.

## 2. DR-47 — `Hex8TbUnitCirculation` re-enabled

Disabled 2026-08-05 because instantiating the kernel chain over a HEX8TB block died with "no
lagrange function available". That blocker is gone: the `ElementType::HEX8TB` case landed with
hex8tb R4 on 2026-08-08 (`cl_IF_InterpolationFunctionFactory.cpp`, folded into the existing
HEX8/HEX8TS case).

Before flipping it I checked the dimensions rather than assuming them, because the test's
comment describes four wall dofs while the helper loops over `number_of_edges()`:

- `TS_TestWall::element()` returns the **HEX8TB** wall (the fixture also builds a HEX8TS layer
  next to it — easy to misread).
- `number_of_edges(HEX8TB)` = **4** (`meshtools.cpp`), against 8 for HEX8TS.
- `EF_HEX8TB` sizes `mE`/`mC` as 3×**4** (`cl_EF_HEX8TB.cpp:36-37`).

So the helper's loop bound, the edge-function column count and the test's 4×4 assertion all
agree, and the `HEX8TS` argument only selects the 8-node corner-coordinate table. My initial
worry — that the helper would index 8 or 12 columns into a 4-column matrix — was **refuted**.

What is *not* resolved is the expectation itself: the test still carries `EXPECTED: pending
Christian sign-off`. That caveat is shared with the already-enabled `Penta6TsUnitCirculation`
and `Hex8TsUnitCirculation`, so enabling it is consistent with how its siblings are treated.
**Its first green run is the sign-off; if it fails, that is information about the wall element,
which is the point of having it.** It is currently the only automated check the HEX8TB wall
has.

## 3. DR-24 — `set_currents` guarded for what is actually checkable

`IWG_Maxwell::set_currents` (`cl_IWG_Maxwell.cpp:92-104`) fixes `aI(k)` onto abstract node dof
`k`, walking `mAbstractNodeDofs` in container order. The caller
(`hphirun.cpp:72-73`, `:99-103`) builds `tI` from the current-BC list in *its* order. Nothing
tied the two orders together.

The register asked for "a one-time assert". Written honestly, that splits in two:

- **Checkable, now checked:** a `BELFEM_ERROR` rejects `aI.length() > mAbstractNodeDofs.size()`
  — the old loop guarded with `if (tCount < aI.length())` and silently discarded any surplus,
  which is a wrong-physics-without-a-message failure. Plus a `BELFEM_ASSERT` for a null dof in
  the fixed prefix (the container is `set_size(..., nullptr)`-initialised).
- **Not checkable here, and the comment now says so:** the *ordering* correspondence is
  established by the factory — current-BC order versus the order `collect_abstract_node_dofs`
  walks `mAbstractNodes` — and there is no key on the dof side to re-derive it from.
  Asserting it for real means threading the BC identity through to the dofs, which is a design
  change, not a one-time assert. Recorded rather than faked.

`BELFEM_ERROR` (not `ASSERT`) for the length check because it is O(1) and `set_currents` runs
once per timestep, not per element — the error tier follows the call rate.

## 4. Still open: DR-59, and it needs a decision not a patch

> **Superseded the same day — see §6.** Christian ruled, a jury round verified, and the bodies
> were deleted. This section is left as written because it is what was true when the batch was
> handed over, and because its "two landmines" framing turned out to be half right: only one of
> the two died with the code.

The retired `populate_rho_database_parallel` bodies are now **unreferenced** in both
`Material_Metal` (`:781`) and `Material_Alloy` (`:556`) — roughly 100 lines each — and they
carry two landmines: a receive loop whose read counter never advances
(`cl_Material_Metal.cpp:828-833`, identical in Alloy), and a dense owned-node packing the
`Database` projector reads back by `node->index()`, correct only while one rank owns
everything. Delete them, or keep them as the starting point for a genuinely distributed build?
If kept, both defects should be fixed now, because the next person to add a caller inherits
them silently. **Christian's call** (`todo/example_deck_and_material_db_repair.md` §8).

## 5. Handoff

Three source files changed, none compiled:

- `tests/math/test_Spline.cpp` — ifdef re-gate (+ `tests/doc/tests_09_spline.md`)
- `tests/fem/test_InterfaceOrientation.cpp` — `DISABLED_` prefix removed
- `src/fem/maxwell/cl_IWG_Maxwell.cpp` — two guards in `set_currents`

The meaningful gate is `make check-fast` with `USE_TEST=ON` (the shared tree has it OFF —
toggle, run, restore). Two of the three changes are *expected* to make previously-silent tests
execute for the first time, so a red result there is a finding, not necessarily a regression.

---

## 6. Addendum — DR-59 resolved by deletion, after a jury round (same day)

Christian authorised deletion conditionally: *"if we are sure that populate_rho_database_parallel
is never called, we can remove it. same with Hex8TbUnitCirculation. I am not even sure if this
even refers to the current HEX8TB with four edges, or an older version we had with eight."*
A `--jury` round was run rather than deleting on my own reading
(`tmp/ai_exchange/review_dead_code_deletion.md`; pre-registration frozen before dispatch,
including a declared bias — I had written both register rows and re-enabled the test that
morning).

**Outcome: delete A, keep B.** Codex and Grok agreed independently on both, and both caught the
same defect in my own brief.

**Candidate A — deleted, 222 lines.** No call site; both members `private` and non-virtual with
`Alloy` not derived from `Metal`, so no base-pointer path; `create_database_mesh` and the
two-arg `populate_rho_database` stay reachable from the serial builders. My self-flagged weakest
claim (orphaning) survived a deliberate attack from Grok, which re-ran it and reported the claim
holds.

**Candidate B — kept.** Christian's recollection was **correct and my pre-registration was wrong
to dismiss it.** I had written that "an 8-edge HEX8TB never existed"; the documentation records
an earlier experimental `HEX8TB`, then an `HEX8TS` wrap, both removed in June 2026 as unphysical
because they could not represent binormal-H at the shell/connector fold
(`src/fem/interpolation/doc/nedelec_thinshell.md:7`). Three generations, not one. What resolves
the actual question is the edge count of each, obtained by walking every historical version of
`src/mesh/cl_Element_HEX8TB.hpp`: `23184fc9` and `abf47b99` (April) are
`ElementTemplate< 8, 8, 4, 0, 0 >`, the file is **absent** at `55e66d29` (deleted 2026-04-17),
and `b42044da` / `d3231b7a` (July/August) are `ElementTemplate< 8, 8, 4, 6, 1 >`. **Every
HEX8TB ever written had four edges;** the eight-edge element was the `HEX8TS`-typed wrap
(`dl20260422_thinshell_facet_renumbering_break.md:68`), and `HEX8TS` still has eight edges today
as live thin-shell machinery — the two sit one case apart in the factory
(`cl_Element_Factory.cpp:189-193`). The test asserts a 4x4 block against a HEX8TB fixture and
postdates the current element by a week, so it cannot be addressing a retired variant.

**Both auditors, independently: my Alloy line numbers were wrong** — the stuck receive counter
is at `cl_Material_Alloy.cpp:598-606`, not `:568-573` (which is the `partition` call). The same
stale citation had propagated into `debt_register.md` DR-59; corrected there.

**DR-59 is only half closed, and the register now says so.** The gather defect died with the
code. The other half — dense `aRho( tCount++ )` packing versus the projector's
`F( tElement->node( i )->index() )` read (`cl_DatabaseProjector.cpp:253-257`) — is in the
**shared** two-arg helper the serial path still calls. It is correct today only because every
node carries owner 0 on the lockstep path. Raised by Grok, verified at source, and it is now a
latent defect in *live* code rather than in a dead branch, which is a worse place for it to sit.

Nothing here was compiled or run; the deletion is source-level.

---

## 7. Correction (2026-08-11): the DR-24 guard was wrong and broke every parallel run

Christian ran helix and hit the guard added in §3:

```
set_currents: 1 currents for 0 abstract node dofs — the surplus would be dropped silently
This error occured on proc 1
```

He suspected a false positive. It was, and worse than that — the check aborts **every parallel
run that has a current boundary condition**.

**Mechanism.** `MaxwellFactory` calls `set_abstract_nodes` inside a `mCommRank == 0` guard
(`cl_MaxwellFactory.cpp:587-592`), so `mAbstractNodes` — and therefore `mAbstractNodeDofs` — is
**empty on every worker rank**. The caller builds the current vector from the global BC count
(`hphirun.cpp`) and calls `set_currents` on all ranks. So on any worker: `aI.length() > 0`,
`mAbstractNodeDofs.size() == 0`, and a `BELFEM_ERROR` fires.

**What I got wrong.** I read `if ( tCount < aI.length() )` as sloppiness that silently discards
a surplus, and "hardened" it. In fact the loop is *deliberately rank-tolerant*: a rank with no
abstract dofs must fix none. The condition I turned into an error is the normal distributed
case, not a defect.

**Fix:** the `BELFEM_ERROR` is removed and the invariant written into the source instead, so the
next reader does not repeat the mistake. The null-dof `BELFEM_ASSERT` stays — it is cheap and
still meaningful.

**Lesson, and it is the same one twice in two days.** I added an always-active check to a path I
had only reasoned about, on the strength of reading a single call site (`hphirun.cpp`) that
happens to run on rank 0. `BELFEM_ERROR` on a path whose rank behaviour has not been traced is
not hardening; it is a new failure mode. The register row for DR-24 now records the reversal
rather than the original claim.
