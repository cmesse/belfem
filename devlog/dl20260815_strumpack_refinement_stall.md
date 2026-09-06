# STRUMPACK Refinement Stall: a Relative Target on a Zero RHS

**Date:** 2026-08-15
**Purpose:** Diagnose why the quench run hung for an hour on a converged
step, and fix the tolerance plumbing that caused it (DR-76)
**Module:** sparse

## The symptom, and the false suspect

The tapestack3d quench run went silent for ~1 h at step 96 (BDF2,
Δt = 25 ms) — on an iterate where both fields were *converging beautifully*
(magnetic −83.63 dB, thermal −52.59 dB at iterate 3). First suspect was the
same-day BC-globals broadcasts; gdb backtraces on all four ranks refuted
that and put every rank inside STRUMPACK `IterativeRefinementMPI` /
`backward_multifrontal_solve`.

## The mechanism

Three facts compose into the stall:

1. `strumpacktools.cpp` applied `set_rel_tol( relative_tolerance() )`
   unconditionally, and nothing ever called `set_maxit` — so STRUMPACK's
   refinement cap stood at its library default of 5000.
2. The shared default was raised 1e-8 → 1e-10 on 2026-08-13 to cure the
   thermal drift (DR-70/DR-72) — a change aimed at the PETSc thermal
   block. The deck's `linear magnetic { library : strumpack ; }` states no
   tolerance, so it inherited the tightening unnoticed.
3. STRUMPACK is a direct solve; refinement is a safety net. On the iterate
   where Newton has essentially converged, the linear RHS is numerically
   zero — and no *relative* criterion is meetable against a zero
   denominator. The net becomes a trap: refinement grinds toward 5000
   sweeps of forward/backward substitution on a 2.5M-dof factorization.

Same disease as DR-72: a relative criterion applied to a quantity that has
become numerically zero. The failure fires precisely on the *best* iterates.

## The fix (approved: "make it so")

The root defect was one global number serving two masters — a Krylov
stopping criterion (PETSc, where the 1e-10 default is the drift cure and
stays) and a refinement target for an already-exact direct solver. Split
by explicitness:

- `SolverParameters` gained `mHaveRelativeTolerance` — set by the deck
  parser and by `set_relative_tolerance()`, copied in the copy ctor, and
  carried in `synchronize()` (payload widened 7 → 8 uints; STRUMPACK
  options are applied per rank, so a divergent flag would configure the
  distributed solver inconsistently).
- `strumpacktools` applies `set_rel_tol` **only when a tolerance was
  stated**. Absent, STRUMPACK keeps its own default and refinement
  converges in one or two sweeps behind the exact factorization.
- `set_maxit( 50 )` bounds the net either way. Deliberately the *second*
  line of defense: the wrapper escalates any non-SUCCESS return to
  `BELFEM_ERROR`, so capping alone would have converted the hour-long
  stall into a run abort on a perfectly solved step. With the conditional
  rel_tol in front, the cap only fires for a deck-stated target that is
  genuinely unmeetable — which now fails in minutes with STRUMPACK's named
  NO_CONVERGENCE message naming the knob to loosen.

Bridge for the currently running binary (which applies the default
unconditionally): the tapestack3d deck now states
`relative tolerance : 1e-8` in `linear magnetic` — the exact value
STRUMPACK received for months before 2026-08-13 — so the next restart is
cured without a rebuild.

Input contract updated in the same session: reference §4.1 row and schema
note now document the stated-only STRUMPACK semantics and the cap.

## By-catch

The A/B question that started the evening (2 vs 4 OMP threads) was
answered along the way: restarting from the same memdump, iterate 1 agreed
to all printed digits and the difference grew from iterate 2 — the
signature of reduction-order roundoff under a different thread count, not
of the day's code changes. The residual differences are ~3 parts in 1e4 on
a divergence, well under the formulation's own noise floor.

Status: **reviewed, not verified** — syntax gates green on all four sparse
TUs (with `BELFEM_STRUMPACK` active, so `set_maxit` compiled against the
real API). Executable gate: rebuild, then the previously stalling step 96
must complete in bounded time; A/B with and without the deck line.
Verify round: `tmp/ai_exchange/strumpack_refinement_fix.md`.

## Follow-up: the headroom gate, and two review catches

Christian, reading the deck after the fix, applied the headroom rule to the
new line — "the relative tolerance in linear can't be larger than the
nonlinear one" — which is a theorem for an iterative solver and a category
error for a direct one (the STRUMPACK number is a refinement exit
threshold, not the delivered accuracy). Rather than leave that distinction
to a comment, it is now enforced where it is true:

- `Solver` exposes `parameters()` (const);
  `Controller::check_iterative_solver_headroom()` hard-errors at setup
  when a field's linear solver is PETSc and its relative tolerance is
  looser than that field's nonlinear tolerance ("impossible tolerance
  pairing..."). Called from set_params for both fields and repeated in
  set_thermal_kernel for the late-link path. Direct solvers exempt, and
  the header comment says why. Deck sweep: no shipped deck uses PETSc;
  garber is MUMPS; the tapestack3d thermal pairing passes with three
  decades of headroom.
- The tapestack3d deck comment on the 1e-8 line now explains
  demand-vs-delivery instead of merely citing DR-76.

Codex round-1 verify of the DR-76 fix: five PASS, one FAIL worth having —
the `SolverParameters` copy ctor silently dropped `mUseMetisNodeNDP`, and
since `Solver` takes its parameters by value, a deck's
`metis nodendp : false` never reached STRUMPACK. Fixed (copy ctor +
carried in `synchronize()`, payload 8 → 9). Input contract updated for the
gate (reference §4.1, schema cross-key note).

Round 2 verify: all five questions PASS — collectivity of the gate,
solver-exists-before-set_params in every execution mode, the parse order
that keeps the late thermal check from comparing against the 1.0e-6
constructor default, the widened synchronize payload, and doc accuracy.
Its one wording note earned a code comment: "direct solvers are exempt" is
shorthand for the real predicate, `type() != PETSc`. That distinction has
teeth, because STRUMPACK under `compression scheme : blr` runs a genuine
GMRES over a lossy factorization — the one configuration where the
headroom rule binds and the gate stays silent. Left ungated on purpose
(compression is off by default and no deck enables it), but the header
comment, the reference row and the schema note now say so out loud rather
than implying full coverage.

Executable gates owed at the next rebuild: `make check-fast`; the
previously stalling step 96 completing in bounded time; and a deliberate
bad-pairing deck (PETSc linear 1e-6 against nonlinear 1e-8) aborting at
setup with the named message.
