# Conditioning Diagnostic: the PETSc-Native Estimate

> **CANCELLED 2026-09-03** (todo/ currentness sweep, round 3): low relevance — a PETSc-native conditioning metric for a diagnostic whose MUMPS path already works; plan v2 was never revised after its audit. Status lines and checkboxes below are as they stood at closure and are not maintained.

**Date:** 2026-08-29
**Purpose:** Let a PETSc-backed field report a conditioning number from the Krylov solve it is
already doing, instead of paying for a shift-invert eigen round. Mechanism:
`KSPSetComputeSingularValues()` before the solve, `KSPComputeExtremeSingularValues()` after.
**Module:** `src/sparse` (+ `src/fem/kernel` routing)
**Environment:** PETSc **3.25.4** ( `/opt/scls/mkl/include/petscversion.h:5-12`, verified
2026-08-29 ). Earlier notes in this repository say 3.25.2 — that is wrong for the installed tree,
and the difference matters because the audit's API findings were taken against 3.25.4.
**AIs involved:** Claude (draft), Codex (audit)
**Status:** **PLAN v2 — AUDITED 2026-08-29, verdict "revise before implementation", six material
corrections folded in below.** No source modified.

> **Scope guards:**
> - This is R7 of `todo/conditioning_shift_invert_fallback.md`, unblocked because shift-invert now
>   works and prints every timestep. It is an OPTIMISATION of a working diagnostic, not a fix.
> - The eigen path stays. It is the only route to a true `kappa_2` and the only route at all for a
>   backend with nothing native.
> - Out of scope, explicitly ( Christian, 2026-08-29 ): an option to force an ARPACK round even
>   when the solver offers a native number.

---

## 1. The one thing that must not be fudged: it is a DIFFERENT NUMBER

`KSPComputeExtremeSingularValues` returns `sigma_max / sigma_min` of the **PRECONDITIONED**
operator. It is not `kappa_2( A )`. Confirmed against PETSc's own documentation by the audit.

> **CITATION RETRACTED.** v1 wrote "`M^-1 A` for left preconditioning, which is what this deck runs
> ( `cl_SolverPETSC.cpp:450` )". That citation is **to the wrong file**: line 450 of
> `cl_SolverPETSC.cpp` is vector allocation. The string "left preconditioning" appears at line 450
> of `cmake-build-debug/tapestack3d/out_serial.txt` — a RUN LOG — and I carried the line number
> across without noticing. The deck does appear to use left preconditioning, but the evidence is
> PETSc's own `KSPView` output in that log, NOT the source.
>
> And BELFEM never sets the side at all: there is no `KSPSetPCSide` anywhere in `src/`, so the side
> is PETSc's per-KSP-type default and could change with the method. **The implementation must call
> `KSPGetPCSide()` and attach the side to the reported quantity**, rather than assume it. PETSc's
> left convention is `B ~ A^-1` giving `BA`, not the `M^-1 A` this plan wrote.

That is not a technicality, it is the whole point of a preconditioner: a good `M` makes that ratio
SMALL while `kappa_2( A )` is unchanged. On the thermal system, `kappa_2( A ) = 6.44e7` measured;
a well-behaved ASM-preconditioned GMRES might report `1e2`. **Both are correct. Printing them under
one label would be a lie**, and a worse one than the MUMPS-COND1-versus-kappa_2 confusion that
started this whole line of work, because the gap is larger and the number looks healthier.

So the footer must name it. The end-of-step block already names the quantity per field
( "MUMPS COND1, rhs-dependent" / "kappa_2, spectral" ), and this adds a third. **Label corrected on
the audit's advice** from "sigma ratio, preconditioned" to **"preconditioned kappa_2 estimate"** —
the standard name, and it carries the two facts the old label dropped: that it is a kappa_2, and
that both singular values are Krylov ESTIMATES rather than computed values.

One qualification the audit added and I had wrong by implication: "a good preconditioner drives the
ratio toward 1" is an idealised tendency, not a definition. For a nonsymmetric or non-normal
operator, fast GMRES convergence does not require a small Euclidean singular-value ratio, so a large
reported ratio does not by itself mean the preconditioner is doing badly.

## 2. What each PETSc configuration can actually supply

BELFEM resolves `AUTO` to GMRES for an iterative PC and PREONLY for `PCLU`
( `cl_SolverPETSC.cpp:625-634` ). The three cases are genuinely different:

> **CORRECTED 2026-08-29 ( audit §2 ). v1's "GMRES-family only" was FALSE**, and the CG row was
> wrong on top of it. PETSc 3.25.4 ships `KSPComputeExtremeSingularValues` implementations for
> **CG, GMRES, MINRES and FETI-DP** — CG installs
> `KSPComputeExtremeSingularValues_CG` when singular-value computation is requested, so it needs
> NO switch to `KSPComputeEigenvalues`. O3 rested on a false premise and is struck.

| deck | KSP | what is available |
|---|---|---|
| `preconditioner : asm` etc. ( the thermal deck ) | GMRES | supported — the target of this plan |
| `krylov method : cg` | CG | **supported directly**, same API |
| ( if ever used ) | MINRES | supported |
| `preconditioner : lu` | PREONLY | **no Krylov estimate.** But note it reports `its = 1`, NOT zero: PREONLY applies the preconditioner once and sets `KSP_CONVERGED_ITS`. **An iteration floor therefore cannot identify PREONLY** — capability must be gated on the KSP TYPE. If the factor package is MUMPS, `MatMumpsGetRinfog( F, 10 )` is the separate route |

**The failure behaviour is the dangerous part, and v1 did not specify it** ( audit §2, from PETSc's
`itfunc.c` ):

- **never armed** → returns `PETSC_ERR_ARG_WRONGSTATE`. It does not quietly return zeros.
- **armed but the active KSP has no implementation** → returns **SUCCESS** and writes
  `emax = emin = -1.0`.

So the naive ratio is `( -1 ) / ( -1 ) = 1.0` — **a plausible, healthy-looking, entirely false
conditioning number**, produced with a success code. That is precisely the silent-wrong class this
campaign exists to remove, and it would have shipped. NaN-on-unavailable is therefore sufficient
ONLY with all of:

- both outputs initialised before the call;
- the returned `PetscErrorCode` checked ( note `KSPGetIterationNumber` at
  `cl_SolverPETSC.cpp:223-224` currently ignores its status — do not copy that );
- non-finite rejected;
- `smin <= 0`, `smax <= 0` and `smax < smin` rejected — this is what catches the `-1` case;
- capability gated on the RESOLVED active KSP type, queried after `AUTO` resolution
  ( `cl_SolverPETSC.cpp:628-637` ), not on the deck enum;
- convergence required as `reason > 0` from `KSPGetConvergedReason` ( already read at
  `cl_SolverPETSC.cpp:219-242` ), not merely "not a known divergence".

## 3. The estimate's quality is bounded by the iteration count — and that is a real problem here

`KSPComputeExtremeSingularValues` estimates from the Krylov space GMRES actually built. Its accuracy
therefore rises with the ITERATION COUNT, and:

- a **well-preconditioned** solve converging in 5 iterations gives a Krylov space of dimension 5
  and a correspondingly poor estimate — the better the preconditioner, the worse the estimate;
- BELFEM sets no GMRES restart, so PETSc's default of 30 applies. **O1 is now CONFIRMED and v1
  UNDERSTATED it** ( audit §3 ). PETSc's own documentation says: *"Disable restarts ... otherwise
  this estimate will only be using those iterations after the last restart."* So the estimate uses
  the iterations **since the last restart** — at most 30, and possibly a handful if convergence
  happens just after one. A 90-iteration solve does not give a 90-dimensional estimate.

This is the opposite of the eigen path's failure mode, and worth stating plainly: shift-invert gets
*more* accurate the harder the problem, this gets *less*. A deck whose thermal solve converges in a
handful of iterations will get a number with very few significant figures, and nothing in the API
reports that uncertainty.

**Consequence for the design, corrected:** `KSPGetIterationNumber()` is already called
( `cl_SolverPETSC.cpp:224` ) but **the total iteration count is the WRONG quality proxy** — it can
be 90 while the estimate rests on 3 directions. What is needed is the LAST-CYCLE dimension:
`KSPGMRESGetRestart()` and `its % restart`, handling exact-cycle termination ( where the modulus is
0 but the cycle was full ). That effective dimension is what should be logged and what any floor
should test.

Simplest alternative worth weighing: **disable restarts for the armed solve**, as PETSc's own
documentation suggests. That buys a full-dimension estimate at the cost of GMRES memory growing
with the iteration count — acceptable for one solve per timestep on this deck, but it changes the
solve's own behaviour, which a diagnostic should be very reluctant to do. Recorded as O5, undecided.

## 4. Where it hooks — mirror the MUMPS pattern exactly

The MUMPS shape is already proven in this deck and should not be reinvented: arm at
`initialize_timestep`, capture and disarm on the FIRST solve of the step
( `Controller::arm_conditioning_*` / `capture_conditioning_*` ). PETSc gets the same lifecycle:

```
arm      : KSPSetComputeSingularValues( ksp, PETSC_TRUE )   before KSPSolve
capture  : KSPComputeExtremeSingularValues( ksp, &smax, &smin )   after a CONVERGED solve
disarm   : KSPSetComputeSingularValues( ksp, PETSC_FALSE )
```

> **BLOCKING CORRECTION ( audit §1 ): this cannot simply mirror MUMPS.** `KSPSetComputeSingularValues`
> is a SETUP-TIME flag, not a per-solve toggle like MUMPS's `ICNTL(11)`, so arming and disarming it
> around one iterate of each timestep does not have the clean lifecycle the MUMPS pattern relies on.
> It is also logically COLLECTIVE in MPI — every rank must take the same arm/disarm path, which the
> first-iterate pattern must therefore guarantee rather than assume.
>
> Two consequences the implementer has to resolve BEFORE writing R1:
> 1. what the correct point to set the flag actually is, given KSP setup lifetime; and
> 2. whether leaving collection ENABLED for the whole run is in fact cheaper and simpler than
>    toggling it — in which case the sample can be taken from the most informative CONVERGED solve
>    of the step ( largest effective Krylov dimension ) rather than mechanically from the first,
>    which the audit points out may have the smallest Krylov space of any solve in the step.
>
> That second option is strictly better information for the same cost, and it dissolves O2: the
> floor stops being an accuracy certificate and becomes a display policy.

## 5. The architectural piece: retire the type-switch

`Controller::compute_conditioning()` currently decides the route with
`solver()->type() != SolverType::MUMPS`, and `Solver::set_mumps_error_analysis()` reaches the
wrapper by `reinterpret_cast` to the concrete type ( `cl_Solver.cpp:284-295` ). With a SECOND native
backend that pattern stops scaling — every new backend would add another `||` to a condition in the
Controller and another cast in `Solver`.

`get_cond0()` is **already a `Wrapper` virtual** ( `cl_SolverWrapper.hpp:201` ), so most of the
mechanism exists. The proposal is to finish it:

```
virtual bool supports_native_conditioning() const ;   // default false
virtual void arm_native_conditioning( bool ) ;        // default no-op
virtual real get_cond0() const ;                      // exists ; PETSc overrides
```

The Controller then asks `supports_native_conditioning()` instead of comparing an enum, MUMPS keeps
its current behaviour behind the same interface, and the `reinterpret_cast` in `Solver` can go.
**This is a refactor of a working path**, so it is sequenced AFTER the PETSc estimate is proven,
never bundled with it — see the step order.

> **CORRECTED ( audit §5 ): do NOT overload `get_cond0()`.** I asked whether reusing one accessor
> for MUMPS COND1, a preconditioned kappa_2 estimate and a true spectral kappa_2 repeats the
> original sin, and the answer is yes. Three different quantities behind one `real` accessor is
> exactly how "conditioning" came to mean two incomparable things in the first place. The interface
> must return a TYPED result carrying the quantity's IDENTITY alongside its value — so the footer
> label is derived from the data rather than re-decided by the caller, and a new backend cannot be
> added without declaring what it measures.

## 6. Ordered steps

- [ ] **R1** — PETSc wrapper: `arm_native_conditioning()` storing a flag, applied via
      `KSPSetComputeSingularValues` at the point the KSP is configured; `get_cond0()` override
      calling `KSPComputeExtremeSingularValues` and returning `smax / smin`. Gated on GMRES-family
      AND a converged solve AND an iteration count above the floor; returns `BELFEM_QUIET_NAN`
      otherwise, which the Controller already treats as "no sample".
- [ ] **R2** — availability reporting: `supports_native_conditioning()` on the PETSc wrapper answers
      for the ACTIVE configuration, not for PETSc in general — PREONLY must answer false.
- [ ] **R3** — Controller: route a PETSc field to the native number when available, eigen otherwise;
      name the quantity in the footer as a third label. No behaviour change for MUMPS or for a
      backend with nothing native. (after: R1, R2)
- [ ] **R4** — the executable gate. **v1 called the two numbers "cross-checkable"; they are NOT**
      ( audit §6 ): `kappa_2( A )` and an estimate of `kappa_2( P( A ) )` SHOULD disagree, so
      comparing their magnitudes proves nothing and would be a misleading test. What the deck gate
      can validate is routing, labelling, finiteness and reproducibility. **Numerical correctness
      needs a separate unit test**: a small matrix where the preconditioned operator is formed
      explicitly and its kappa_2 computed directly, then compared against PETSc's estimate. Plus
      focused cases for each failure mode named in §2 — unsupported KSP ( the `-1, -1` trap ),
      unarmed, zero- and one-direction Krylov spaces, restarted GMRES, and CG.
- [ ] **R5** — retire the type-switch ( §5 ), as a separate change with its own round. (after: R4)
- [ ] **R6** — `doc/input_file_reference.md` + `doc/input_schema.yaml` if any deck key changes, and
      the support matrix in `sparse_usage_guide.md`: which solver yields which quantity.

## 7. Open questions

- **O1 — what does a GMRES restart do to the estimate?** BELFEM sets no restart, so PETSc's default
  30 applies. If the Hessenberg is discarded per cycle, a long solve reports from its last cycle
  only and the estimate is *worse* than the iteration count suggests. Verify in the PETSc source.
- **O2 — the iteration floor. REFRAMED by the audit: this is a DISPLAY POLICY, not an accuracy
  certificate.** No threshold certifies the estimate, so the honest framing is "below this
  effective dimension we do not show a number", and it must be tested on the LAST-CYCLE dimension,
  not the total iteration count ( §3 ). **Not decided** — and note the §4 option of sampling the
  step's most informative converged solve largely dissolves it.
- ~~**O3 — CG.**~~ **STRUCK**: rested on the false premise that CG needs `KSPComputeEigenvalues`.
  PETSc supports the singular-value API on CG directly, so CG is simply one of the supported types.
- **O5 ( NEW ) — disable restarts for the armed solve?** It buys a full-dimension estimate, at the
  cost of changing the behaviour of the solve being measured. A diagnostic altering the thing it
  observes needs an explicit decision. **Not decided.**
- **O4 — three labels in one footer cell.** The cell is six characters
  ( `conditioning_string`, `cl_FEM_Controller.cpp` ) and the quantity name sits beside it. With
  three possible quantities across two fields, is the end-of-step box still readable, or does the
  naming want a legend printed once per run?
