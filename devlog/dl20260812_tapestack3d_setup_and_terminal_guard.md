# tapestack3d setup: periodic mesh, terminal definition, and a descriptive guard for over-specified terminals

**Date:** 2026-08-12
**Purpose:** Record four setup defects found while bringing the 8-tape soldered stack case up, and the one source change made in response
**Module:** `src/homology`, plus the `cmake-build-debug/tapestack3d` case files

## Summary

Four failures in sequence, three of them case-configuration errors and one a source
robustness gap. Each aborted at a different stage of `MaxwellFactory::create_magnetic_kernel`.
Only the last one produced a source change.

| # | Symptom | Stage | Cause | Fix |
|---|---|---|---|---|
| 1 | `The Curve seems to branch or is not properly connected: 1: 0 2: 0 1+2: 1` | `CurveFactory::sort_end_nodes` | terminal curves for inner tapes intersected with the air end face | intersect with the solder end face instead |
| 2 | `Number of facets does not match` | `PeriodicityFactory::map_facets` | mesh not periodic; gmsh meshed the two end planes independently | `Periodic Surface` constraints in the `.geo`, remesh |
| 3 | `Row index 2 out of bounds, which must be smaller than 2` | `Cohomology::updatekGeneratorsFromHomology` | 8 current conditions requested, topology provides 2 generators | one condition over the whole soldered stack |
| 4 | same message again | same | `a:b` outside brackets makes each id its own condition | bracket the range |

## 1. Terminal curves must meet a face that carries the tape edge

In 3-D a `topology{curves}` entry `a @ b` is the node set shared by sidesets `a` and `b`
(`MaxwellFactory::create_curves`, `CurveFactory::intersect`). The case used the air end
face (49 at z=0, 50 at z=L) for every tape. That face is the air cross-section with the
stack cut out of it, so its boundary runs along the two *outer* tapes only. For the six
inner tapes the intersection collapsed to the single tip node where the tape meets the air
face, giving a curve with one node and no segments — exactly the `1: 0 2: 0 1+2: 1` counts.

The end face that carries an inner tape's edge is the end face of the solder volume below
it: 42…48 at z=0 and 28…34 at z=L. With those, every one of the 32 curves resolves to 20
nodes / 19 segments, and the two halves of each tape share exactly the centre node.

## 2. The periodic planes have to be meshed as copies

`PeriodicityFactory::map_facets` collects every sideset lying entirely in the source plane
and every sideset in the target plane, drops one temporary node per facet centroid, and
requires the counts to match before it can pair them with the kd-tree. The `.geo` carried no
`Periodic Surface` statement, so gmsh meshed the two planes independently: 5300 facets at
z=0 against 5110 at z=L, and the air end faces had **zero** coincident centroids. Even the
solder faces, whose counts coincided at 158 because their boundary discretisation is shared,
only agreed on about half their centroids.

Adding

```
For k In {1:numTapes-1}
	Periodic Surface { s2+k } = { s4+k } Translate { 0, 0, domainLength };
EndFor
Periodic Surface { s } = { s-1 } Translate { 0, 0, domainLength };
```

and remeshing gives 5300 = 5300 with a worst centroid mismatch of 3.9e-8 mm.

## 3. One galvanically connected conductor carries one condition

The tapes are thin shells; the volumes between them are solder. Every solder volume touches
the tapes on both sides, so the conductor is one connected slab and the air region is a
periodic cylinder minus one wire: `H^1` has rank 2, one generator encircling the conductor
and one free axial generator from the periodicity. That leaves exactly one condition, but
the case declared eight (`[1,2],[3,4],…`).

The terminal is now the solder cross-section given as **sidesets**, not curves. The key
chosen in `MaxwellBoundaryConditionFactory` selects the algorithm: `input curves` sets
`tIsThinShell` and `Homology::suggest_Homology` then builds the generator as the plain sum
of the listed curves' segments, whereas `input terminals` takes the boundary of the terminal
2-chain — the envelope of the whole cross-section, which encircles the tapes *and* the
solder. The curve form is wrong here twice over: it misses the solder current, and the six
inner tape curves lie inside the conductor, so they are not in the air complex at all.

The per-tape curves stay defined in `topology{curves}`; the thin shell claims them by
sideset in `MaxwellFactory`, independently of the boundary-condition list.

## 4. Brackets are load-bearing in terminal lists

`Section::get_id_groups` opens a group on `[` and flushes it on `]`. Outside a bracket every
id is flushed as its own group, **including each member of an `a:b` range**. So
`input terminals : 42:48 ;` is not one terminal of seven sidesets, it is seven conditions,
and it reproduced failure 3 with an identical message. The correct form is
`input terminals : [42:48] ;`.

## Source change: a descriptive guard instead of a backend assertion

`updatekGeneratorsFromHomology` sizes `tTransInv` as `n x n` with `n` the number of
cohomology generators, then writes one row per row of `tPseudoInv`, which has one row per
condition. With more conditions than generators the loop walked off the end of the matrix and
surfaced as a raw Blaze bounds assertion — a message that names neither the input key
responsible nor the quantity that is too large.

The guard that covers the related failure (`w.n_cols()+tCount==n`, "Inconsistency in the
definition of the terminals…") sits *after* that loop, so it never had the chance to fire.

Added a `BELFEM_ERROR( aNumConditions <= n, … )` immediately after `aNumConditions` is
computed, naming both counts and stating the three things a user needs: that each
galvanically connected conductor carries one condition, that a periodic model spends one
further generator on the free axial loop, and that ids outside a bracket each become their
own condition. `BELFEM_ERROR` rather than `BELFEM_ASSERT` because this is setup code that
runs once and the condition is an input error that must survive a release build.

### The audit found a second overrun site

A read-only Codex audit against a pre-registered claim set
(`tmp/ai_exchange/cohomology_terminal_guard.md`) refuted two of the five claims, and the
refutation was load-bearing rather than cosmetic: **`aNumConditions <= n` is not the complete
safety condition.**

`smithForm( tF )` returns `t`, the number of nonzero Smith pivots, i.e. the rank of the
terminal coupling matrix. Step 6 writes one row per condition and then one row per free cut,
`tQ_.n_cols() - t` of them. Since `tF` is `n x tNumUnique`, `smithForm` sizes `tQ`/`tQ_` from
`aMat.n_rows()`, so `tQ_.n_cols() == n` and the total row count is `aNumConditions + ( n - t )`.
That is `<= n` exactly when `aNumConditions <= t`. A deck with few enough conditions but
*dependent* ones therefore still overran, one loop further down, and the final
`w.n_cols()+tCount==n` check — the one that names linear dependence — was still unreachable
for it.

Re-derived from source before acting on it (`fn_Smith.cpp:144-148` for the `tQ_` sizing,
`:185-196` for `t` as the pivot count) rather than accepted on the auditor's report. Finding
confirmed; a second `BELFEM_ERROR( aNumConditions <= t, … )` now sits immediately after
`smithForm`, naming the rank and the three configurations that produce dependence.

The two guards are deliberately separate: `<= n` means "more conditions than the topology can
carry", `<= t` means "the conditions you gave are not independent". Different causes, different
remedies, so they get different messages.

The existing rank/independence check is untouched and is **not** made dead by either guard:
with `aNumConditions <= t` the assembled rows are in bounds but can still fail to be
independent, which is precisely what it catches.

## Evidence

**Verified by execution:**

- The corrected deck clears the cohomology: the run wrote `tapestack3d.bfm` (538 MB), which is
  only produced after the cuts and cohomology are computed. Setup failures 1–4 are gone. The
  transient itself was not run to completion.
- The periodic remesh: 5300 facets on each plane, worst centroid mismatch 3.9e-8 mm.
- **The first guard fires and renders.** A copy of the deck with the brackets removed
  (`input terminals : 42:48`) aborts with the intended message, wrapped correctly in the error
  box, reporting "7 current or voltage conditions, but the topology only provides 2 cohomology
  generators". That also **confirms n = 2 empirically**, which until then was a topological
  prediction plus a matrix bound read off a traceback.

**Reviewed, not verified:**

- The second guard (`aNumConditions <= t`). It compiles and rests on the same arithmetic as
  the first, but has not been observed to fire. Constructing a deck that reaches it is awkward
  in this model: with `n = 2` the capacity is one condition, so a two-condition deck trips the
  pre-existing kernel check rather than the rank guard. A model with more generators would be
  needed to exercise it.

## Confidence

- Failures 1, 2 and 4: **high**, each verified directly against the mesh file or the parser
  source, and 1 and 2 confirmed by the case getting past those stages afterwards.
- Failure 3 diagnosis (n = 2, one condition available): **high** — `n` is read straight off
  the traceback, and it is what the topology predicts.
- `input terminals : [42:48]` being the complete terminal fix: **medium**. The bulk branch
  also folds in every `topology{curves}` curve whose edges the terminal sideset flags, i.e.
  all 16 front tape curves on top of the envelope. Whether that is the intended handling for
  a thin shell piercing a bulk terminal or a double count depends on sign conventions in
  `Chain::getBoundary()` and `Segment::edge_direction()` that were not resolved statically. If
  it double counts, the symptom is a factor on the imposed current, not a crash — worth a
  sanity check on the total current in the first output step.

## Addendum (late evening): the DR-08 freeze reproduced, and the fix

The overnight throughput experiment (8×2 vs 4×4 ranks×threads, warm-started from
the same `memdump.hdf5`) reproduced the DR-08 coupled freeze with the cleanest
possible control — identical deck, mesh, binary and state, differing only in
decomposition. The thermal Picard-2 residual came out 0.000100 on 8×2 and
0.000102 on 4×4: ordinary partition-order roundoff, straddling
`tolerance switch : 1e-4`. The 8×2 run promoted to Newton and converged in six
iterates; the 4×4 run stayed Picard and froze bit-flat for four. **The parallel
decomposition selected the algorithm.** Deterministic on replay, which rules out
the entire race/leak/corruption family — a logic defect.

Two defects behind one symptom: the promotion gate deadlocks (Newton requires
`eps2 < epsSwitch2`, but when Picard cannot reduce the residual, Picard is
exactly what must get below the switch for Newton to be allowed to run), and
`KSPSolve`'s convergence outcome was never read (error code only — a diverged or
maxits-terminated solve returned 0 and was accepted silently; DR-45's defect
class). Leading hypothesis for the flatness itself: the KSP relative tolerance
is met by the nonzero initial guess, so the solve returns its own input after
zero iterations.

Fix applied same evening (Fable session, per the critical-refactor convention),
three parts, `todo/thermal_picard_freeze_dr08.md`: **R1** converged-reason +
iteration-count instrumentation with hard error on divergence; **R2** the
magnetic `try_escalate_to_newton` precedent ported as its thermal twin, fired
from the existing stagnation detector — whose own comment already recorded this
exact signature ("thermal frozen at −39.5 dB for 100+ iterations") and answered
it with only a warning — at two flat iterates, once per attempt; **R3** demotion
hysteresis, one decade. Deviation from the plan as pre-registered: the
escalation clears the `mJustPicard2` chatter latch rather than respecting it,
because the magnetic twin does the same and the once-per-attempt cap replaces
the latch's protection. The staggered `iterate_thermal` handoff is deliberately
untouched (plan O1). All TUs syntax-clean both configs; **reviewed, not
verified** — gate is the 4×4 reproducer, which must now promote by the third
flat iterate.

The 8×2-vs-4×4 timing numbers are confounded by the algorithm split (6 vs 9
iterations) and the experiment must be repeated; per iteration the difference
was ~8 %, not the ~39 % the step times suggested.

## Addendum 2 (night): tape 1 upside down — root cause, and signed sidesets

Christian read the first visualization frames and saw the bottom tape's layer
stack mirrored against the other seven. Traced to the end, the chain is:

1. `fix_facet_masters` normalizes cross-type facet masters by DOMAIN TYPE
   (`Air=1 < Conductor=5` ⇒ the conductor side becomes master — the enum
   header says explicitly that these values "control the order in which
   master and slave priority is given").
2. That rule is **mirror-symmetric across a tape stack**: it puts the master
   on the solder side of BOTH air-adjacent tapes, which are geometrically
   opposite sides. A translationally uniform stacking direction is
   unrepresentable — one outer tape always flips. With air below the stack,
   that is tape 1.
3. `Facet::flip()` rewrites the facet node winding to the master traversal
   (`cl_Facet.cpp:36-55`), and `ThinShellFactory::process_nodes_tri3` builds
   the layer normals from that winding — flipped master ⇒ flipped layer stack.
4. The `.geo` is innocent: all 16 tape surfaces carry identical gmsh windings
   (measured), and the as-loaded master (smaller element id,
   `cl_Mesh_ConnectivityCalculator.cpp:227`) is not uniform either — it puts
   tape 8's master on the air side. Neither state is usable as-is.

**Two designs considered and one rejected.** A geometric uniformization pass
(flip whole tapes to the majority mean normal) was drafted and killed on
Christian's objection: a corc wrap has no meaningful mean normal, so the
heuristic could silently flip a currently-correct corc tape. Which side the
layers face is user intent — no local geometry rule can derive it.

**Implemented instead: gmsh-style signed sidesets** (Christian's design).
`sidesets : -5, -6, 7:20 ;` — a negative id flips the orientation of every
facet of that sideset, applied after `fix_facet_masters` (which would undo it)
and before the cut/thin-shell pipeline reads the windings. Pieces:
`Protoshell::mFlippedSideSets` (+accessor), `MaxwellFactory::
read_signed_sidesets` (sign-aware parse; signs bind to single ids, signed
ranges fatal; `sidesets()` stays unsigned so all downstream consumers are
sign-agnostic), `MaxwellFactory::flip_thin_shell_sidesets` (the flip pass,
one Default-level log line per flipped sideset), call site between
`fix_facet_masters()` and `create_cuts_sub_master()`. Unsigned decks take the
exact code path they took before. Input contract updated in both artifacts,
including the sign-whole-sheets-together pitfall (a half-signed tape tears at
the shared centre line and the node-normal averaging cancels to zero there).

**Reviewed, not verified**: `-fsyntax-only` clean in both configurations;
gate is a remesh-free rerun of tapestack3d with `sidesets : -5, -6, 7:20 ;`
and the visualization showing eight uniform tapes. Left open: the
tape-8-vs-loaded-state observation above means the domain-type rule and the
element-id rule disagree about the air-adjacent tapes generally — any model
with air on both sides of a stack needs signs on one outer tape.

## Addendum 3 (night): solver defaults — measured, not assumed

Two default changes, both driven by numbers that only existed once `KSPSolve`'s
converged reason and iteration count were instrumented (R1 of the DR-08 plan).
Step time over the evening, same deck and mesh: **4:01 → 0:30, about 8x.**

### PETSc preconditioner: GAMG/JACOBI → ASM

The old default was `comm_size() > 1 ? GAMG : JACOBI`. Measured on the thermal
block (337k dofs, rtol 1e-8): GAMG needed **~74 Krylov iterations on the two
Picard iterates and ~700 on every Newton iterate of the same step** — a 9x
penalty appearing exactly where the algorithm flips. Cause: smoothed
aggregation assumes a symmetric elliptic operator, which describes the Picard
system (M/dt + K, SPD) but not the Newton tangent, which carries the
non-symmetric mixed term from `mt_thermal_h`. ASM costs slightly more on the
Picard operator (~107) and does not degrade on the Newton one.

**The wall-clock mechanism was NOT the one I first assumed.** The thermal solve
was 520-672 ms against STRUMPACK's 24-33 s for the magnetic factorisation — it
was never the bottleneck. The 5x came from the coupled iteration converging in
**2 iterates instead of 9**, because a better-conditioned thermal correction
means fewer sweeps, and each sweep costs a magnetic factorisation. Recorded
because "we changed the thermal preconditioner and the run got 5x faster" is
true but the obvious reading of it is wrong.

Second, independent reason to drop the conditional: the old default made the
preconditioner **a function of the rank count**, so serial and parallel runs of
one deck solved with different numerics — the same defect class as DR-08.

### STRUMPACK: METIS_NodeNDP enabled by default (`metis nodendp` opts out)

STRUMPACK's own diagnostic recommended it after detecting a 56-level tree and
warning about stack-overflow segfaults. Verified after the change:

| | NodeND | NodeNDP |
|---|---|---|
| levels | 56 | **14** |
| separators | 160,475 | 103,939 |
| supernodal tree | "built from etree" | "from METIS_NodeNDP" |
| factor nonzeros | 703,896,611 | 710,951,079 (+1%) |
| factor memory | 5631 MB | 5688 MB (+1%) |
| nd time | 7.4 s | 9.0 s |

**The fill is unchanged; only its SHAPE changed** — and that is the whole point,
so the register must not record this as a memory reduction. Per-rank RSS at
peak went from one rank at 23 GiB against siblings at 2.8 GiB, to 4.79 GiB on
rank 0 and **1.84-2.06 GiB across all seven workers** (<12% spread). Both of
the evening's OOM kills were peak-on-one-rank, not aggregate, so a balanced
tree was worth more here than a smaller factor would have been. The ordering
phase pays 1.6 s more per setup for it.

The mechanism is exactly what the API doc states: NodeNDP returns the separator
tree, NodeND does not, so STRUMPACK rebuilds one from the etree — and that
reconstruction was what produced the depth.

### Open

- `factor nonzeros` unchanged means the fill itself is still the memory floor;
  a genuinely smaller factor would need a different ordering (SCOTCH untested)
  or compression (BLR — off by design, stalls Newton).
- The ASM default rests on one deck. GAMG still wins on symmetric operators and
  a deck that knows this should ask for it explicitly.
- Keyframe trims cost a full step for near-zero simulated time when Δt lands
  just short of a save point (observed: 28 s for 0.0208 ms). Harmless once Δt
  reaches the cap and every step lands on a keyframe. A stretch-instead-of-
  shorten rule in the trim would remove it; the Δt stash it would need already
  exists.

## Addendum 4 (early morning 2026-08-13): the DR-07 |B| tangent, audit-gated

Christian's three-task order executed while the quench run grinds through
1.0–1.24 Ic: deck retune ( target iterations 50→3, max iterations 20/15, the
magnetic min-relaxation floor raised to 0.01, save every 50 ms ), the silent
rejection message restored ( printed AFTER the Δt update so it cannot lie —
the cut is a clamp, not a halving, at the floor ), and the DR-07 |B|
derivative chain implemented under the full cross-review cycle:
pre-registration → Codex + Grok design audits → implementation → Codex + Grok
result audits, unanimous pass.

The audits earned their place twice. First: both voices independently refuted
the plan's assumption that the existing Newton beta channel could be reused —
`add_rho_field_tangent` differentiates bj_angle ( field–current, metal
Kohler ) while HTS uses bn_angle ( field–tape-normal ), so binding
mFundRhodBeta would have applied the wrong ∂β/∂q rows. β stays return_zero
with the reason recorded at the dispatch site. Second: Codex caught that
`rho_piecewise` returns the RAW power law in its PL regime — no parallel
combination — so the piecewise B-derivative follows its own residual rather
than the powerlaw convention; the same look found the EXISTING
drho_piecewise_dJ carrying a parallel factor its residual does not ( mild,
recorded ), and flagged the compute_drhodT_hts reconstruction b = a·c/(c−a)
as possibly sign-wrong, which is why the T-leg is deferred rather than rushed.
The fact-check before pre-registration also caught my own 2026-08-12 plan
note stating a spurious ln10 on the B-leg chain factor ( it cancels; both
auditors confirmed ) — now corrected in the plan file.

Landed ( 576 insertions, 7 files, -fsyntax-only clean in both
configurations ): JcFunction::deval_dB/_dbeta/_dT with zero defaults that do
not abort ( zero is exact for constants, and the Newton consumer early-outs
on it — constant-jc decks are bit-identical by construction ),
clamp-consistent JcFunctionDatabase overrides, null-safe djc/dn routing
helpers, Material::drho_powerlaw_dB and drho_piecewise_dB with defect
overloads ( jc AND its derivative modulated by the same D ), 8 MaxwellData
wrappers mirroring their compute_rho_* preambles, and mFundRhodB bound in all
8 HTS dispatch paths.

**Reviewed, not verified.** Gate: the quench deck's 14–50-iterate Newton
steps at 1.05–1.2 Ic — 35 % of tonight's wall clock for 2.8 % of its
simulated time — must drop substantially on the recompiled binary, with the
overnight run as the before-measurement.

## Addendum 5 (2026-08-13 morning): the J/Jc asymmetry — three angle conventions, one ruling needed

Christian read the first quench visualization and asked why J/Jc is 1.4 on the
−x tape half but 1.02 on +x when B and J themselves look symmetric, suspecting
the postprocessor. Three-way audit (Claude probe + Codex + Grok, independent
inventories in agreement): the plot and the solve use DIFFERENT angle
conventions, and a third one hides in bulk HTS.

- **Solve:** `Calculator::bn_angle` folds the field-to-tape-normal angle to
  [0, 90°] via `abs(dot(n,b))` — used by every thin-shell residual, tangent,
  defect variant, and (via compute_rho) the thermal Joule heating.
- **Plot:** the thin-shell postprocessor uses the signed dot, [0, 180°].
- **Bulk HTS:** dummy π/2 everywhere — always the ab-plane column of a table.

The fold is Kim-era and documented (`todo/closed/beta_angle_normal_sign.md`),
valid for ModifiedKim (even in θ) and invalid for the measured sp-ap table,
which spans [0, 180°] with real asymmetric pinning: jc(θ)/jc(180−θ) up to
1.42 at 82–86 K / 0.3–0.5 T (BELFEM's own spline; Codex spot check 1.461).
Christian's observed 1.37 sits inside that band — the plot is the honest
lookup of the attached table; the solve symmetrizes measured physics and
samples only one lobe of the data.

Registered as **DR-69 (P1, ruling)**: switching the solve to signed is the
physically consistent move for table-based Jc, but it makes the SIGN of the
tape normal material input — the gmsh-style signed sidesets become
load-bearing, old unsigned decks can change physically, Kim decks must stay
folded (the material has to advertise whether its Jc(θ) is even), and a
future β-tangent needs a signed bn_angle chain of its own. Not implemented —
this is a physics-contract decision, not a bug fix.

**Same morning, the DR-07 gate passed on live evidence:** after Christian's
recompile, the warm-restarted run crossed 1.30 → 1.38 Ic in 21 consecutive
steps under 6 iterates with one rejection — the same regime that cost the old
binary 16–37 iterates per step and a rejection every few steps. Throughput
through the knee roughly 3×. The |B| tangent does what the audits said it
would.

---

## Addendum 6 (2026-08-13): seven debt-register rows closed against the campaign

Christian's ruling: the 8-proc coupled quench campaign in
`cmake-build-debug/tapestack3d` — the first production run to exercise the
whole stack at once — closes **DR-07, DR-08, DR-09, DR-33, DR-34, DR-36 and
DR-38**. Each row's status cell in `todo/debt_register.md` records exactly
what the run did and did not cover; the honest deltas against the original
gates:

- **DR-07** — |B| tangent verified by execution (the ~3× knee, addendum 5).
  β channel and T-leg remain open in `powerlaw_jc_n_field_derivatives.md`,
  not in the register.
- **DR-08** — R1 verified live (the KSP counts drove the ASM/GAMG finding);
  R2/R3 in the running binary, freeze never recurred. R1b/R4 retired unrun;
  R5/O1 stand as scope notes in the plan file.
- **DR-09** — superseded by the full 4×4 + 8×2 campaign.
- **DR-33** — BDF1→BDF5 ramp and BDF5 proper ran all night, coupled, with
  T-dependent ρcp and savepoint restores; the fixed-Δt Newton-count A/B is
  waived with the row. Unit coverage (nine cases) stays in `check-fast`.
- **DR-34** — hundreds of STRUMPACK factor+solve cycles at 4 and 8 ranks on
  1.6 M dofs; serial and 2-rank never separately run, accepted.
- **DR-36** — closed as nothing-actionable: MPICH is out of scope by standing
  decision. The run is Open MPI and adds no MPICH evidence.
- **DR-38** — repeated warm restarts from memdump mid-quench exercised the
  cold-start-on-load path; the bitwise A/B vs an uninterrupted run and the
  CORC circuit coupling are waived with the row.

Plan files updated in the same pass: `thermal_picard_freeze_dr08.md` and
`bdf_nonlinear_mass_verification.md` marked CLOSED with their retired boxes
struck; `powerlaw_jc_n_field_derivatives.md` status upgraded to verified for
the |B| channel. `todo/README.md` bullets refreshed. `check_doc_claims.py`
32/32 throughout.
