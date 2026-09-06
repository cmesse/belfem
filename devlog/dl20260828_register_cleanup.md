# Register Cleanup: DR-95, DR-85, DR-126 Closed; DR-134 and DR-135 Filed

**Date:** 2026-08-28
**Purpose:** Close DR-95 ( comment demotion applied at all three sites ), DR-85 ( Christian's ruling; stale-doc residue split into DR-134 ), and DR-126 ( superseded by the autopin feature; SideSet sibling split into DR-135 ).
**Module:** fem/kernel ( Controller ), fem/maxwell ( postproc + factory, record only ), todo records

## What was done

DR-95's remaining half — the record-hygiene half was fixed 2026-08-26 — was the
comment demotion itself, comment text only, no behavioural change. Applied on
Christian's approval ( "make it so" ), wording taken from the register's own
suggested resolution, which was the 2026-08-19 jury round's only 3/3 finding
plus Grok's verified V6 broader-state correction, so the wording was
pre-audited; no fresh vendor round was run for a comment-only landing.

1. **`cl_FEM_Controller.hpp:53-62`** ( was `:53-59`; register cited the stale
   `:53-56` ): the `mDeltaTimeInitial` doc comment no longer asserts as fact
   that the solver "derives its permutation from the FIRST Jacobian and keeps
   it for life". It now states the measured invariant — the dt of the first
   post-restart step decides the restart's outcome and cannot be repaired
   afterwards — and hedges the mechanism as "Working hypothesis, not proof",
   using the `.cpp`'s wording: STRUMPACK solver-entry state ( ordering,
   matching, equilibration, tiny-pivot threshold ).
2. **`cl_FEM_Controller.cpp:3662`** ( register cited the stale `:3503-3504` ):
   the parse comment's "a warm restart re-enters here, never at the dumped
   step" was false as written — the cap is `if ( mDeltaTime >
   mDeltaTimeInitial )` at `:4390`, so a dump at or below the initial step IS
   adopted verbatim. Replaced with the schema's accurate contract wording:
   "keeps a dumped delta_time only up to this value".
3. **`todo/closed/dr92_restart_step_ramp.md:5-9`**: the framing note "The
   defect is which matrix STRUMPACK's permutation is derived from" demoted to
   the same hypothesis-plus-invariant form with an inline DR-95 date stamp;
   the historical body stays as written.

The already-hedged mechanism block at `cl_FEM_Controller.cpp:4373-4378`
( "Current working explanation, not proof" ) was the template and is untouched.

## Register

DR-95 struck in `todo/debt_register.md`; the `[P]` open count recounted
mechanically, 27 → 26. Status: **reviewed, not verified** in the trivial sense
that no executable gate applies to comment text; all three sites and the cap
conditional were read directly this session.

## Second closure, same session: DR-85 struck, DR-134 filed

DR-85 ( J/Jc lookup used the unprojected air-average b ) closed and struck on
Christian's ruling ( "I think we can close DR-85" ). The fix landed 2026-08-16
and was re-confirmed on disk at the strike: `bn = dot( bn, n ) * n` present at
`cl_MaxwellPostprocessor.cpp:550` and `:867`; it had also traced clean in the
2026-08-28 landing-claim sweep. The surviving A/B one-frame gate ( JJC with and
without the projection ) never ran and is **waived by the ruling** — recorded
as waived, not passed.

The row's stale-doc twin did not close with it: `h_ts_metal` / `h_ts_metal_t` /
`h_ts_hts` / `h_ts_hts_t` and `mt_maxwell_h.cpp` exist nowhere under `src/`
( grep-verified ), yet `dof_manager_usage_guide.md` §20 still documents the
family as current code at ~6 sites, including the §20.3 walkthrough that shows
the very unprojected average DR-85 fixed in the live path. Split into
**DR-134** ( P3, doc-only ) so the residue survives the strike. Register
bookkeeping: `[P]` count 26 → 27 ( DR-134 in, recounted mechanically ), the
untagged-rows note now names DR-100 alone.

## Standing rule recorded mid-session

Christian: **"comments are free to be changed, code changes have to be audited
by the jury."** Both closures this session are comment/doc/register text only,
so no jury round was owed; the rule is saved to persistent memory.

## Third closure, same session: DR-126 struck as superseded, DR-135 filed

DR-126 ( the exactly-singular magnetic system: unpinned φ constant per
thin-shell buffer patch, floating air level despite the bearing ) closed and
struck on Christian's ruling: "with today's autopin method, DR-126 is
obsolete." Superseded by the autopin feature in `04c33664` ( "add autopins" ):
`MaxwellFactory::find_autopins` ( `cl_MaxwellFactory.cpp:3154`, using the new
`fn_Graph_multibfs` ) labels every connected φ component and selects one
safely pinnable node per component — ERROR when a component offers none — and
the factory fixes each pin to 0.0 at setup, chasing a periodically-condensed
candidate to its NODE source first; the pin set persists in the `.bfm`, with a
loud master-rank warning on pre-autopin files. That retires the flagged
residue (b) ( per-patch anchors ) by construction and defuses residue (i)
( gauge pin after graph freeze ) in the common path, since the factory-time
pins that residue named as its requirement now exist on every component.

**Reviewed, not verified:** the row's own witness — condest baseline dropping
by orders plus a frame-stable air/buffer φ constant on the same deck — has not
yet run on an autopin build; it stays the natural check at the next tapestack3d
launch and is recorded in the struck row.

Residue (ii) is not covered by autopins and split into **DR-135** ( P3,
`[CODE]` `[W]` ): `SideSet::impose_dirichlet` ( `cl_FEM_SideSet.cpp:288` )
calls `fix( aValue )` unconditionally, so on a periodically-condensed dof the
pin is silently discarded by the T-matrix — the exact mechanism the 2026-08-28
bearing fix rerouted around. Mechanism confirmed by reading; reachability
unaudited ( no shipped deck known to put a symmetry sideset on a periodic
face ).

Register bookkeeping was merged with a parallel session's edits ( DR-88 and
DR-100 struck there ): after both sessions, `[P]` = 25 and `[W]` = 7, both
recounted mechanically.

## Fourth closure: DR-89 struck under the DR-42/DR-49 pattern

DR-89 ( magnetic solve exiting at STRUMPACK's library-default absolute
tolerance 1e-10, turning small-||b|| startup convergence into a GMRES
lottery ) struck on Christian's ruling. The state that justified it, verified
on disk this session, not from the devlog trail:

- the fix landed 2026-08-17 and is in the tree: `absolute tolerance` is a full
  deck key, `mAbsoluteTolerance = 1e-14` default
  ( `cl_SolverParameters.hpp:162` ), STRUMPACK always applies it, PETSc only
  when deck-stated; the key appears at three sites in
  `doc/input_schema.yaml`; suite-tested at landing ( ctest 14/14 )
- the original X2 gate ( t=50/100 ms memdump replay ) became unrunnable
  2026-08-24 when the dump was deleted; the analog run on the surviving
  t=7100 ms dump passed the DR-89-relevant parts ( restart mechanics clean,
  ~3 iterates/step, no startup lottery, no crawl )
- the analog's one anomaly, a -80 dB floor at deep quench, was booked to the
  DR-88 conditioning class — and DR-88 has since been struck as superseded by
  DR-126's gauge singularity, itself closed by the autopins ( `04c33664` ),
  so the residual anomaly now points at a closed, plausibly fixed cause. The
  stale "OPEN DR-88" clause in the row was updated in the same edit.

The `[F]` half was already settled: the user-visible contract ( deck key +
default ) was decided and landed a week ago. What remained was exactly the
DR-42/DR-49 shape — design finished, residue a single run — and Christian
took that option over keep-open-retag-`[P]` and strike-and-retire-gate.

**Surviving gate ( struck is not verified ):** the next campaign cold start,
where the small-||b|| abs_tol path is genuinely exercised — startup endgame
must reach < -115 dB deterministically with no iterate crawl and no delta-t
collapse.

The `[F]` list is now DR-65 alone.
