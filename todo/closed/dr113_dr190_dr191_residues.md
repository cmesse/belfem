# DR-113, DR-190, DR-191: the residues of the small-fixes cluster

**Date:** 2026-09-30
**Purpose:** retire the three rows opened or narrowed on 2026-09-29 by paying their recorded
residues. DR-190: after a caught `create_field` refusal the equation still points at the freed
`DofManager`; an `IWG::detach_field()` in a `try`/`catch` around `set_equation` ends that window.
DR-191: the MUMPS wrapper hands a fresh Fortran-`allocate`d instance to `JOB = -1`, and the 5.9.1
driver reads `KEEP(40)` from it before initializing anything; the pool zeroes that one field at
allocation. The Kohler switch vectors of the five metals that still size-then-store move to the
brace assignment the other four already use, so an unset entry becomes impossible instead of
detectable; the transverse-branch test gains the two near-switch samples that pin the stored
value. DR-113: no code — its `make check` leg ran on 2026-09-29 and the row said "owed" until this session's text edit.
**Module:** `src/fem/kernel`, `src/fem/iwg` (DR-190); `src/sparse`, `src/physics/materials`,
`tests/physics`, `tests/fem` (DR-191); `todo/debt_register.md` (DR-113)
**AIs involved:** Claude (plan + implementation), Codex `gpt-6-astra` + Grok `grok-4.7` (juries;
the kernel puts the round on the astra row)
**Status:** ✅ COMPLETE 2026-09-30 — plan jury and code jury both Codex astra/high + Grok 4.7/high, reconciled (plan round: four P2 adopted, Silver widened; code round: no P0/P1/P2, four P3 applied). Landed: `IWG::detach_field` + the `try` through the push + three `BELFEM_ERROR` guards (DR-190); `KEEP(40)` zeroed on a fresh MUMPS pool, the near-switch test samples, brace-assigned switch vectors in five metals, Silver's low-field cubics repaired (DR-191); DR-113's verified `make check` leg recorded. Verified on a private non-MKL Armadillo debug tree (positive controls on HEAD, clean with the fix, both mutations red on the predicted sample, `test_physics` 159/159, `test_fem` 300/300). DR-190 and DR-191 struck (Christian's rulings O2–O6). `make check` in the shared tree, Christian, 2026-09-30: 17/17 (fast 8, mpi 4, 21.24 s real).

> **Scope guards:**
> - DR-190: the detach and the catch in `create_field` only. No change to what is refused, to
>   `set_equation`'s order, to the production (MPI-abort) reaction, or to the callers that never
>   catch. No null guards added to `set_omega`, `set_penalty` or `set_wetted_sidesets`: a caller
>   that reuses an equation after a refusal gets a null dereference in every build instead of a
>   read of freed memory, and no in-tree caller does either.
> - DR-191 (MUMPS): one field of the fresh pool slot, cited to the driver line that reads it. No
>   change to the `JOB = -1` argument set, to the reuse path, or to the rollback.
> - DR-191 (Kohler): the switch vectors only — same sizes, same values, same store sites. The
>   polynomial vectors, the Bezier objects, and the `kohler()` bodies are untouched. The tin
>   longitudinal cubic's sign (§1.5) is a physics ruling, not a step of this plan (O2).
> - DR-113: register text only. Gate ( a ) stays blocked on the Tier 2 FEM fixture
>   (`todo/test_gap_backlog.md` F1 / T-62) and is not touched here.

---

## 1. Current Behaviour and How It Fails

### 1.1 DR-190: the dangling `IWG::mField` after a caught refusal

`Kernel::create_field` (`cl_FEM_Kernel.cpp:590-613`) owns the `DofManager` in a `unique_ptr`
(`:596`) until the push succeeds (`:611-612`). `set_equation` (`cl_FEM_DofManager.cpp:85-130`)
calls `aIWG->set_field( this )` at `:95` before `create_dofs` and `create_blocks`, so every
refusal inside those leaves `IWG::mField` (`cl_IWG.hpp:138`, set at `cl_IWG.cpp:710-714`)
pointing at the object the `unique_ptr` then frees. Readers of that pointer:
`set_wetted_sidesets` (`cl_IWG.cpp:119`, `:131`, `:137`), `set_omega` (`:992`), `set_penalty`
(`:1002`), and `IWG_Timestep::set_timestepping_method`, which null-guards it (`cl_IWG_Timestep.cpp:226`)
and then walks the field's blocks through a `reinterpret_cast` (`:236-237`). `set_field` also
copies the mesh pointer (`:713`); the mesh outlives the field (the fixture deletes kernel, then
parameters, then mesh, `test_ThinShellSelection.cpp:352-358`), so `mMesh` stays valid and is not
part of the window. No derived IWG caches the field separately (`grep "DofManager *\* *m"` over
`src/fem/iwg`, `src/fem/maxwell`, `src/fem/thermal` headers: only the two factories'
`mMagneticField` / `mThermalField`, which are assigned from `create_field`'s return and therefore
never hold a refused object).

No in-tree caller catches a `create_field` refusal and reuses the equation (18 call sites:
`staticheat.cpp:71`, `poisson.cpp:48`, `hangingnodes.cpp:59`, `cl_ThermalFactory.cpp:173`,
`cl_MaxwellFactory.cpp:786`, `cl_CutFactory.cpp:2474`, and twelve test fixtures — the only catch
is `expect_refusal`). Under the production abort policy the window
does not exist; under the debug throw policy it is what `expect_refusal`
(`test_ThinShellSelection.cpp:471-485`) leaves behind in the two `*Refused` tests (`:626-655`).

### 1.2 DR-191 (1): the `dmumps_` valgrind context

`mumpstools_create_solver` (`src/sparse/mumpstools.f90`) allocates the pool
`gSolvers( gMaxNumSolvers )` of `DMUMPS_STRUC` (`:65-70`), which Fortran does not zero, sets
the four fields the user guide requires for the initialization call — `COMM`, `PAR`, `SYM`,
`JOB` (`:104-107`; guide 5.9.0: `JOB`, `PAR`, `SYM` "must be initialized by the user … before" the
call, nothing else) — and calls `DMUMPS` (`:113`).

The 5.9.1 driver (`MUMPS_5.9.1/src/dmumps_driver.F`, from the tarball in
`~/codes/belfem_transfer/tmp/`, the version `libdmumps.so` in `/opt/scls/mkl` reports) does this
first on `JOB = -1`:

```
      IF ( JOB .EQ. -1 ) THEN
        id%INFO(1)=0
        id%INFO(2)=0
        OLDJOB = id%KEEP( 40 ) + 456789                 ! :942
        IF ( OLDJOB .EQ. -1 .OR. OLDJOB .EQ. 1 .OR.
     &       OLDJOB .EQ. 2 .OR. OLDJOB .EQ. 3 ) THEN     ! :943-946
           id%INFO(1) = -3
```

`KEEP(40)` is the driver's own "last successful phase" tag, written at the end of every phase
(`:1066`, `:1074`, `:1303`, `:1564`, `:2125`, `:2210`) and never by the user. On a fresh
`allocate`d instance it is garbage: the comparison at `:943-946` is the conditional jump valgrind
reports, and if the garbage happens to satisfy it the first `JOB = -1` returns `INFO(1) = -3`
("instance not terminated") for no reason — four values in 2³², so never observed, but not zero.
Nothing else is read before the defaults are set: `MUMPS_PROPINFO` (`bcast_errors.F`) receives
`ICNTL(1)` and never uses it, and `DMUMPS_INI_DRIVER` (called at `:1044`) writes `NPROCS`,
`MYID`, and every `KEEP`/`ICNTL`/`CNTL` default through `DMUMPSID` before reading anything
beyond `COMM`, `PAR` and `SYM`. `KEEP` is a fixed `INTEGER :: KEEP(500)` component
(`/opt/scls/mkl/include/dmumps_struc.h:173`), so it is addressable before the call.

A reused slot is not affected: the free path runs `JOB = -2`, which stores `-2 - 456789`
(`:1074`). So the clear belongs at pool allocation, not at every create — zeroing on every
create would also hide the driver's own reuse check. **Recorded, out of scope (O6):** the
defaults routine `DMUMPSID` (`:61`, from `DMUMPS_INI_DRIVER`) zeroes `KEEP` (`:1151`) and then
writes `KEEP(40) = -1 - 456789` (`:1333`), so a `JOB = -1` that fails after reaching the defaults
leaves `OLDJOB = -1`, which the entry check refuses on the next `JOB = -1`; the wrapper's rollback
(`mumpstools.f90:133-143`) runs no `JOB = -2`, so such a slot refuses every later re-init. Pre-existing.

Confidence high on the field (driver source, this session); this replaces the "recollection, not
verified" of the previous plan.

### 1.3 DR-191 (2): the test window

`MetalFieldDependence.KohlerTransverseBranchesMeetAtTheStoredSwitch`
(`tests/physics/test_MetalFieldDependence.cpp:165-247`) checks the low-field cubic at `x/4` and
the log polynomial at `2x` (`:213-216`). A stored switch anywhere in `( x/4, 2x )` passes both.
The branch gap at `0.9x` and `1.1x`, computed from the same match (this session):

| metal | `0.9x` | `1.1x` |
|---|---|---|
| lead | `+7.2e-5` | `-5.7e-5` |
| tin | `+3.6e-5` | `-2.7e-5` |

relative to the log polynomial — five orders above the `1e-10` tolerance, so samples there pin
the switch to `( 0.9x, 1.1x )` and a `1e-6` self-check that the two oracles differ at each sample
keeps the test honest if the coefficients ever move.

### 1.4 DR-191 (4): the size-then-store Kohler switches

Nine metals carry Kohler switches. Four assign the whole vector at once and cannot leave an
entry unset: Chromium `mBSKohlerSwitch = { x, x1 }` (`cl_Material_Chromium.cpp:284`), Aluminum
(`:243`), Iron (`:336`), Indium (`cl_Material_Indium.cpp:195`). Five size first and store later:

| metal | size | stores |
|---|---|---|
| Copper | `cl_Material_Copper.cpp:279` | `:285-286` |
| Nickel | `cl_Material_Nickel.cpp:187` | `:195-196` |
| Silver | `cl_Material_Silver.cpp:257` | `:330`, `:336` |
| Lead | `cl_Material_Lead.cpp:205`, `:207` | `:247-248`, `:270` |
| WhiteTin | `cl_Material_WhiteTin.cpp:212`, `:214` | `:263-264`, `:286` |

All entries are written today (the previous round wrote the two that were not). The
size-then-store form is what made DR-191 possible; the brace form removes the class. Both
backends' `operator=( std::initializer_list )` produce a vector of the list's length with its
values, including the one-element case (`cl_AR_Vector.hpp:338-353`, `cl_BZ_Vector.hpp:394-411`:
one element is `set_size( 1, value )`).

### 1.5 Recorded, not in scope: the tin longitudinal cubic

`WhiteTin::create_kohler` matches the longitudinal low-field cubic `p0` to the log polynomial at
`BS = e^{2.5} ≈ 12.2` with value, slope and curvature (`cl_Material_WhiteTin.cpp:242-246`); the
resulting linear coefficient is negative, so `p0` is below zero for `BS` under about 4, where
`BS = B · S` and `S = rho_ref / rho_0T` (`cl_Material_WhiteTin.hpp:76`, `cl_Material_Metal.cpp:1013`;
the register's `-2.9e-4` is at `B = 1 T`, `S = 1`). The database build evaluates every field node
without a cutoff, so the corner is reached. A magnetoresistance ratio below zero is unphysical;
the remedy is a choice between clamping, a constrained re-match, and leaving a `≤ 4e-4` relative
defect. That is O2, Christian's call, not a step.

### 1.6 Found by the plan jury (Grok), confirmed on the source: Silver's low-field cubics

`Silver::create_kohler` (`cl_Material_Silver.cpp:255-345`) builds both low-field cubics at the
stationary point `w0` of `q1 − p1` — where the two log polynomials are farthest apart, not
where they cross (`:271-282`; the jury caught my "crossing"); `w0 = -5.3446`, so the cubic branch is
`BS < 4.77e-3`). Block 1 (`:281-302`) computes the longitudinal cubic from `p1` and stores it in
`q0`; block 2 (`:304-317`) computes the transversal cubic from `q1` but reuses block 1's `dgdx` and
`d2gdx2` (`:310`); block 3 (`:319-328`) recomputes the longitudinal cubic with `d2dgw2` where
`d2gdx2` belongs (`:323`), writes `p0( 0 )`, and then `q0( 1 )`, `q0( 2 )`, `q0( 3 )`. Net:
`p0( 1..3 )` are never written, `q0( 1..3 )` carry longitudinal coefficients, and `q0`'s curvature
was wrong before that. The unwritten entries are zero, not garbage: `Cell::set_size` is
`std::vector::resize` (`cl_Cell.hpp:199-202`), which value-initializes each `Vector`, and
`Vector() = default` (`cl_AR_Vector.hpp:51`) is not user-provided, so the whole object including
Armadillo's local buffer is zeroed first — the gate's valgrind pass reports nothing for silver
while the gtest goes red (run evidence, this session). So the shipped low-field longitudinal
curve was `p0( 0 ) · BS³` alone: deterministic, wrong, and invisible to valgrind. Silver is a tape-stack component in seven shipped decks. The mechanical
repair is block 1 → `p0`, block 2 recomputing its own derivatives, block 3 deleted, plus a Silver
case in the far-point test; every `silver_RRR*.hdf5` cache written before it keeps the wrong table
(`cl_Material_Metal.cpp:893` naming, same acceptance weakness as DR-191's). Whether it rides this
round or is filed is O4.

**Bottom line:** three small, independent mechanisms — a pointer that outlives its object, one
integer the library reads before anyone writes it, and a construction pattern that permits an
unset entry — each with a valgrind- or test-observable gate, plus one register cell that lags a
run that already happened.

## 2. Approach

**DR-190.** `IWG::detach_field( const DofManagerBase * aField )` clears `mField` only when it is
`aField`, and nothing else (the mesh pointer stays, §1.1; a throw before the link at
`cl_FEM_DofManager.cpp:95` — from `:87-91` — must not clear a pointer to some other live field).
`create_field` protects everything from the link to the completed hand-over, because
`Cell::push` is `std::vector::push_back` (`cl_Cell.hpp:265-267`) and can throw after
`set_equation` returned (Codex, plan round):

```cpp
try
{
    tField->set_equation( aEquation );
    mDofManagers.push( tField );
}
catch( ... )
{
    // set_equation links the equation to the field before it can refuse;
    // the field dies with tOwner, so the equation must not keep pointing at it
    aEquation->detach_field( tField );
    throw ;
}
tOwner.release() ;
```

`tOwner` is declared outside the `try`, so the field is still alive while the handler runs and is
freed when the exception leaves the function — the detach never touches it. The `check_mesh`
refusal (`:599-603`) precedes the link and needs no detach. `catch( ... )` with a rethrow is the
house form (`cl_MaxwellBoundaryConditionFactory.cpp:36`, `cl_ThermalBoundaryConditionFactory.cpp:41`,
`cl_Material_Alloy.cpp:245`). The three setup-time readers of `mField` (`set_wetted_sidesets`,
`set_omega`, `set_penalty`) get a `BELFEM_ERROR( mField != nullptr, … )` if O5 is ruled yes
(both auditors; `coding_philosophy.md:639` puts once-per-run preconditions on the always-active
tier). Rejected either way: an early return on null — a silent no-op.

**DR-191 (1).** In the allocation block of `mumpstools_create_solver` (`:65-70`), beside the
`gOccupied` fill:

```fortran
        ! the 5.9.1 driver reads KEEP(40), its last-phase tag, at the top of
        ! JOB = -1 (dmumps_driver.F, the OLDJOB check) before it initializes
        ! anything, and a fresh allocate leaves it undefined
        forall( k=1:gMaxNumSolvers ) gSolvers( k )%KEEP( 40 ) = 0
```

Zero decodes to `OLDJOB = 456789`, which none of the four refused values match. One field, not
the whole array: the driver sets the other 499 itself and a broader clear would suggest a
contract the guide does not state. Rejected: a clear on every create (hides the driver's reuse
check, §1.2) and a clear of `KEEP` wholesale.

**DR-191 (2).** Two more samples in the transverse test, each against its own branch oracle at
`1e-10` relative, plus a `1e-6` relative self-check that the oracles disagree there.

**DR-191 (4).** The five metals assign their switch vectors with braces at the site where the
last entry is known today; the `set_size` lines go. Silver's first store (`:330`) moves down to
join the second (`:336`) — `x` is not reassigned in between (`:331-335` compute `w1` from `p1`).
Lead and tin: `mBSKohlerLongSwitch = { std::exp( w1 ), std::exp( w2 ) }` and
`mBSKohlerTransSwitch = { x }` at the existing store sites. Rejected: the NaN fill plus finite
assert (O1 of the previous plan) — a detector for a class the brace form abolishes.

**DR-113.** The `make check` of 2026-09-29 (commit `31be8ee`, 17/17 in one of the two shared trees, the other not run) ran the
17 registered suites, which include `commmpi_np2` and `commmpi_np4` (`ctest -N` in
`cmake-build-debug`, this session); `DrainCheck.CleanFabricIsSilent` and
`DrainCheck.StrayScalarFiresOnEveryRank` landed in `52465d8` (2026-09-23) and were in that
binary. Gate ( b ) is verified at 2 and 4 ranks; the row's "make check owed" clause is rewritten
to say so. Text only, exempt from the jury.

## 3. Gap Table

| # | item | class | handled by |
|---|---|---|---|
| 1 | `IWG::mField` after a caught refusal | (c) | R2: detach in the catch |
| 2 | `IWG::mMesh` after the same refusal | (a) — the mesh outlives the field | no change (§1.1) |
| 3 | element/facet flags left by `collect_element_indices` after a refusal | out of scope, recorded on DR-190 | — |
| 4 | a throw inside `DofManager`'s own constructor | out of scope, recorded on DR-190 | — |
| 5 | `KEEP(40)` on a fresh pool slot | (c) | R3 |
| 6 | `KEEP(40)` on a reused slot | (a) — `JOB = -2` wrote it | no change (§1.2) |
| 7 | transverse switch pinned to `( 0.9x, 1.1x )` | (c) | R4 |
| 8 | size-then-store switch vectors, five metals | (c) | R5 |
| 9 | stale `Pb*Sn*_RRR*.hdf5` caches | operational, recorded | devlog + register |
| 10 | tin longitudinal cubic sign | (b) | O2 |
| 11 | DR-113 gate ( b ) run evidence | text | R1 |

## 4. Ordered Steps

- [x] **R1** plan jury (`review_dr113_dr190_dr191_residues`, Codex `gpt-6-astra`/`high`, Grok
  `grok-4.7`/`high`) — done 2026-09-30, reconciled; the DR-113 register cell was rewritten in the
  working tree during the round (text; it lands with the approval).
- [x] **R2** DR-190: `IWG::detach_field( const DofManagerBase * )` (`cl_IWG.hpp` beside
  `set_field`, `cl_IWG.cpp`), the `try`/`catch` through the push in `Kernel::create_field`, the
  three `BELFEM_ERROR` guards if O5 says so; the two `*Refused` tests call
  `set_timestepping_method( EulerMethod::BackwardDifference1 )` on the equation after the refusal —
  the field walk behind the null guard becomes a no-op after the fix and is an invalid read of the
  freed field before it (valgrind is the observer; the gtest line only proves it does not crash).
- [x] **R3** DR-191 (1): the `KEEP( 40 )` fill in `mumpstools_create_solver`.
- [x] **R4** DR-191 (2): the `0.9x` / `1.1x` samples and their self-check in the transverse test.
- [x] **R5** DR-191 (4): brace assignment of the switch vectors in Copper, Nickel, Silver, Lead,
  WhiteTin.
- [x] **R5b** (O4: this round, Christian 2026-09-30) Silver's low-field cubics repaired per §1.6, a Silver
  case in the far-point test (both branches, at `x/4` and `0.9x` against the test's own C2 match).
- [x] **R6** code jury on the own-files diff, same slug and depth (`git diff HEAD -- <files>` to
  `tmp/ai_exchange/<slug>.diff`).
- [x] **R7** gates (§7) on a private non-MKL debug tree (`/opt/scls/debug`, Armadillo, in the
  session scratchpad — none exists on disk today; all three trees are MKL); mutation checks
  recorded; Christian's `make check-fast` on the shared trees.
- [x] **R8** close-out (rulings 2026-09-30: O2 leave and record, O3 strike with the print retired, O4 repaired here, O5 guards added, O6 on the closed row), conditional on the rulings: DR-191 struck only if O2 is ruled (leave =
  a recorded by-design closure, not a residual), O4 is either repaired here or filed elsewhere,
  and O6 has a home; DR-190 struck only if O3 retires the print and the recorded items (flags
  after a refusal, a throw inside the constructor, the rank-0 barrier) are ruled recorded-not-owed
  — else both rows stay open, narrowed. DR-113 cell current; devlog, `todo/README.md`, R6 of
  `todo/dr190_dr191_small_fixes.md` ticked.

## 5. Open Design Questions

- **O1** — **RESOLVED 2026-09-30 → brace form.** brace form (R5) or NaN fill plus finite assert? Proposal: brace form; it abolishes the
  class rather than detecting it, costs nothing at run time, and is the form four metals already
  use. Jury's view welcome; Christian decides.
- **O2** — **RESOLVED 2026-09-30 (Christian) → leave, record as accepted.** the tin longitudinal cubic below 4 T (§1.5): clamp at zero in `WhiteTin::kohler`,
  re-match with a through-origin constraint, or leave and record. Physics; Christian's call.
  Proposal: leave for this round, keep the register line.
- **O3** — **RESOLVED 2026-09-30 (Christian) → strike, print retired, the three items recorded not owed.** DR-190 after R2: the only run residue is the per-rank print of the gate halves on
  `tapestack3d` at np 4, whose premise §1.1 of the previous plan refuted from the source (the
  gate is built from rank 0's shared list). Proposal: retire the print with the row
  (register preamble, the DR-42/49 exception) and strike, with the three recorded items
  (§3 rows 3–4 and the rank-0 barrier) ruled recorded, not owed. Christian's call.
- **O4** — **RESOLVED 2026-09-30 (Christian) → repaired in this round (R5b).** Silver's low-field cubics (§1.6): repair in this round as R5b (widening DR-191 by the
  extend-don't-branch policy; a value change in a shipped material with a stale-cache consequence),
  or file and leave DR-191 open with it. Proposal: repair here — it is the same mechanism as the
  lead/tin match, the fit coefficients are untouched, and the far-point test pins it.
- **O5** — **RESOLVED 2026-09-30 (Christian) → yes, added.** `BELFEM_ERROR( mField != nullptr, … )` in `set_wetted_sidesets`, `set_omega`,
  `set_penalty` (both auditors). Proposal: yes, three once-per-run checks.
- **O6** — **RESOLVED 2026-09-30 (Christian) → DR-191's closed row, recorded.** the MUMPS rollback residue (§1.2): once DR-191 is struck, does it live on DR-191's
  closed row as recorded, or on a new sparse-wrapper row? Proposal: closed row, recorded; a
  `JOB = -2` on rollback is a one-line wrapper change that can wait for a caller that hits it.

## 7. Definition-of-Done Checklist

- [x] `valgrind --leak-check=full --track-origins=yes test_fem --gtest_filter='ThinShellSelection.*Refused'`
  on the private tree: 0 invalid reads, 0 `DofManager` leak records (before R2, with the new
  test lines: invalid reads inside `set_timestepping_method` — the positive control).
- [x] `valgrind --track-origins=yes test_physics --gtest_filter='AlloyRRR.ARhoDatabaseForAnotherRRRIsRefused'`
  on the private tree: 0 `dmumps_` contexts (before R3: 1).
- [x] `test_physics --gtest_filter='MetalFieldDependence.*:MetalRRR.*:AlloyRRR.*'` green; red
  with either transverse store shifted to `0.8x` and, separately, to `1.2x` (mutation of both
  new samples).
- [x] if R5b: `valgrind --track-origins=yes test_physics --gtest_filter='MetalFieldDependence.*'`
  0 `Silver::kohler` contexts (before: the cubic-branch reads of `p0( 1..3 )`).
- [x] `test_physics` and `test_fem` full runs green on the private tree, both backends if the
  scratchpad budget allows (the switch vectors are read on both; the brace path is exercised by
  the four brace metals on both already).
- [x] Christian: `make check` in the shared tree, Christian, 2026-09-30: 17/17 (fast 8, mpi 4, 21.24 s real).
- [x] Register: DR-113 cell rewritten; DR-190 and DR-191 per O3; the stale-cache note kept on
  DR-191's closed row.

## 8. Audit Trail

`tmp/ai_exchange/review_dr113_dr190_dr191_residues.md`; the driver excerpt and the user-guide
text extraction in the session scratchpad (`mumps/`); the branch-gap numbers of §1.3 from a
Python evaluation of the test's own match formulas.
