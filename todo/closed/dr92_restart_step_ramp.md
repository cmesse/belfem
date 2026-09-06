# DR-92: A Warm Restart Must Not Let the Solver Be Born at a Large Timestep

**Date:** 2026-08-19
**Purpose:** Theory of the DR-92 restart anomaly and the proposed re-entry cap.
[ Title and framing corrected 2026-08-19: this is NOT a growth limit. The
working hypothesis, not proof, is which matrix STRUMPACK's solver-entry state
( ordering, matching, equilibration, tiny-pivot threshold ) is derived from at
process birth; the measured invariant is that the first post-restart dt decides
the outcome. Wording demoted 2026-08-28 ( DR-95 ). ] Step 1 of Christian's five-step plan ( theory+audit -> design+audit
-> code+audit -> test -> document+prose sweep ).
**Module:** fem/kernel ( Controller ), fem/iwg ( IWG_Timestep )
**Status:** **COMPLETE 2026-08-19 — fix landed and VERIFIED BY EXECUTION ( all four
acceptance criteria ); make check 14/14. Five-step round finished: theory+audit,
design+audit, code+audit, test, document.** Historical detail below is kept as
written, including the refuted T3/T4.

**Original status line:** THEORY AUDITED by Codex and Grok 2026-08-19. **T3 and T4 REFUTED by a
follow-up measurement Grok's audit prompted; the mechanism is now O8 ( solver
initialisation ), see §3a.** Numeric errors in the original T3 corrected. Fix not
designed, no code written.
**Evidence:** `tmp/ai_exchange/bnorm_anomaly.md` ( full measurement record ),
register rows DR-92, DR-93.

---

## 1. The phenomenon

A warm restart that re-enters at a large timestep produces a pathological first
solve and usually loses the step, while a live run at the identical step, time
and Δt is healthy. Measured today in a sandbox replica ( Christian's directory
untouched ), all from dumps the campaign itself wrote:

| restore | first Δt | `GMRES it. 0` | Krylov its | outcome |
|---|---|---|---|---|
| memdump_2350 | 0.154 ms | 0.0297 | 16 | healthy |
| memdump_2300 | 7.09 ms | — | — | healthy |
| memdump_2700 **clamped** | **5 ms** | **0.903** | **16** | healthy, 2 iterates, −138.8 dB |
| memdump_2700 | 50 ms | **10851** | 50 = maxit | sick, grinds, rejects |
| memdump_2800 | 50 ms | **21245** | 50 = maxit | sick |
| campaign, **live** | 50 ms | **0.019–0.043** | ~16 | healthy |

`GMRES it. 0` is the LEFT-PRECONDITIONED residual ‖M⁻¹b_scaled‖
( source-verified in `GMResMPI.cpp:71-99`, `SparseSolverMPIDist.cpp:282-320` ),
so it measures solution/correction scale after preconditioning, not ‖b‖.

**The controlled comparison:** `memdump_2700` restored at 50 ms is sick and at
5 ms is healthy. Same dump, same binary, same rank layout, one deck line
changed. [ The original draft concluded "the restored state is therefore not
damaged"; Grok correctly refused that — the clamp shows the dump is not GROSSLY
corrupt, not that it is identical to a live continuation. ]

## 2. What is NOT the cause ( all excluded by measurement, not argument )

| candidate | how it died |
|---|---|
| MPI field-collect clobber ( DR-93 ) | a SERIAL restart reproduces the anomaly; `FieldData::collect` is `mCommSize > 1`-guarded ( `cl_FEM_DofMgr_FieldData.cpp:235` ) |
| partition / aura / ghost ownership | same serial result |
| reduction order, multi-rank front tree | same serial result |
| BDF integrator restore ( `mH` ONLY ) | probe shows `mH` is exactly the dump's `bdf_h` shifted with `bdf_last_dt` prepended. **NARROWED after Codex: the broader 'q-history as consumed by assembly' is NOT excluded** |
| φ/φ0 clobber at abstract nodes ( C1–C6 ) | present IDENTICALLY in healthy restores, including the clamped one |
| imposed current | probe shows `max |fixed value|` = I( t+Δt ), correct |
| step-size RATIO Δt/h₀ | a ratio-1.00 restore at 50 ms ( 21245 ) is sicker than a ratio-10.83 one ( 10851 ) |
| BDF order | BDF1 restore at 50 ms drops it.0 to 4.17 but still burns maxit and reaches only −70 dB — an amplifier, not the cause |

## 3. Theory

**T1. The restored state is not grossly corrupt.** [ CORRECTED after Codex
audit — the original wording "every dumped quantity was verified" and "a
damaged state could not do" were overclaims. ] Every quantity we MEASURED was
restored correctly: fields, `mH`, `bdf_step_count`, `bdf_last_dt`, imposed
current, dof values. That is not the same as verifying all dumped state, nor
identity with a live continuation. What the clamped restore shows is that the
same dump converges in 2 Picard iterates to −138.8 dB at 5 ms.

**T2. A restored process cannot absorb a LARGE first step.** [ CORRECTED —
"the deficiency is not in any dumped quantity" was an overclaim; `mFieldValues`
and the q-history AS CONSUMED BY ASSEMBLY are not excluded. ] Something a live
sequence carries makes a large step survivable. Codex supplied source evidence
for the leading candidate: BELFEM reuses initialized STRUMPACK solver state
across live steps ( `cl_SolverSTRUMPACK.cpp:247-260`, `:288-323` ) and STRUMPACK
itself reuses ordering / matching / equilibration after value-only updates
( `strumpack-local/src/src/sparse/SparseSolver.cpp:153-163`,
`SparseSolverMPIDist.cpp:188-198` ) — so a live sequence solves with a matching
computed from an EARLIER matrix, while every restart computes one fresh.
Controller adaptive state ( ω, ε₀, best-ε, coupling factor,
`mDeltaTimeTemporary` ) is also unpersisted. **Honest gap: we know THAT a large
first step fails, not WHY.**

**T3. ~~The deficiency is transient — one or a few small steps cure it.~~**
**REFUTED 2026-08-19.** Grok's audit noted the grow-back ladder was a SECOND
PROCESS, so it never tested state-healing. The test that does — reduce Δt
WITHIN one process — is already on disk: the sick 50 ms restore rejected step
222 and retried at 25 ms, where it.0 was **9888**, essentially unchanged from
10851. Halving Δt in-process does not cure it.
[ The original T3 also mis-cited its own ladder numbers ( 1.48 at 27.9 ms
instead of 1.63; "5.6×, +37%" instead of 6.5×, +56%; "two snaps" instead of
four ) — all corrected in the exchange record. ]

**T4. ~~The defect is a re-entry / process-state discontinuity amplified by a
large first Δt.~~ REFUTED 2026-08-19, superseded by §3a.**

**T5. This is one instance of a general pattern in the Controller: Δt values
are installed WITHOUT passing through the growth limiter.** Two sites do this:
- `load_memdump` adopts the dumped Δt directly, clamped only against the deck
  min/max ( `cl_FEM_Controller.cpp:4167-4179`, `:4278-4282` — line numbers
  corrected by Codex ).
- The save-point snap path computes a growth-limited value ( `:2517-2535` ) and
  then overwrites it from `mDeltaTimeTemporary` ( `:2538-2542`; snap setup and
  restore at `:2565-2570` ) — producing jumps up to 10.8× against a documented
  BDF5 limit of 1.2×.
Codex found no OTHER post-start Controller path that installs a larger outer Δt
while bypassing the limiter ( other direct assignments are initial setup,
rejection shrink, or IWG history restore ).
The doc already states the intended limits ( `src/fem/doc/timestepping_strategy.md:44-57`,
citing Grigorieff 1983 and Hairer & Wanner ): 1.5 for orders 1–2, 1.4 for 3,
**1.2 for ≥ 4**.

## 3a. The mechanism ( O8, confirmed by measurement 2026-08-19 )

**The permutation, matching and ordering are computed from a process's FIRST
Jacobian and reused for its entire lifetime.** Source ( Grok, high confidence ):
matching is computed once in `reorder_internal` ( `SparseSolverBase.cpp:321-329` ),
`reorder()` no-ops afterwards, and `update_matrix_values` reapplies the stored
`matching_` / `equil_` / permutation ( `SparseSolver.cpp:152-163`,
`SparseSolverMPIDist.cpp:188-198` ), clearing only `factored_`. BELFEM
initialises once and value-updates thereafter ( `cl_SolverSTRUMPACK.cpp:247-260` );
`Wrapper::free()` is never called in the timestep loop ( `cl_SolverWrapper.cpp:61, 83` ).

BDF5 assembles `J = α M + Δt K` ( `cl_IWG_Timestep.cpp:995-1013` ). A small
first Δt weights the mass term; a large one weights the stiffness term. The
measurements say a permutation derived from the small-Δt matrix serves the
whole run — including later 50 ms steps — while one derived from a 50 ms matrix
does not, and cannot be rescued by shrinking Δt afterwards.

[ Honest limit, per Grok: "well-conditioned" is our adjective, not a measured
quantity. STRUMPACK's own equilibration diagnostics are identical in both arms
( `r_cond = 1, c_cond = 1, type = N` ), so the log carries no direct signal. ]


## 3b. Why "just reconstruct the state exactly" cannot work ( Christian, 2026-08-19 )

Christian's question — if the restore were exact, the next step would be
identical and no cap would be needed — has a definitive answer, and it is not
the one the earlier drafts assumed.

**Neither saved scalar is the step the controller would naturally take next,
because every dump sits exactly on a save point.** `save_memdump` runs only at
save points, so the last completed step is ALWAYS the snap-truncated one, and
`mDeltaTimeTemporary` ( the pre-snap intent ) is not persisted at all.
Measured across every dump on disk:

| dump | on save point | `bdf_last_dt` | dumped Δt | mH ( truncation artifacts in bold ) |
|---|---|---|---|---|
| 1800 | yes | 0.05 | 0.05 | 0.05, 0.05, 0.05, 0.05 |
| 2300 | yes | 0.0071 | 0.0071 | 0.0071, **0.019**, 0.017, 0.010 |
| 2350 | yes | 0.000154 | 0.000154 | ×3, **0.0033** |
| 2700 | yes | **0.0046** | 0.05 | 0.045, **0.012**, 0.038, **0.018** |
| 2800 | yes | 0.05 | 0.05 | 0.05, 0.05, 0.05, **0.0046** |

`bdf_last_dt` is systematically a keyframe artifact; the dumped Δt is the
restored pre-snap intent. **This kills `bdf_last_dt` as a cap basis for a second
and deeper reason than O1** — not merely permissive in one case, but an artifact
in every case.

**And the decisive point: exact reconstruction would NOT fix DR-92.** The LIVE
run chooses 50 ms at that state and is healthy. So choosing the "right" Δt is
not the issue. A restart with perfect controller state ( PID terms, pre-snap
intent, iteration history ) would correctly choose 50 ms — and would still be
sick, because the sickness is a property of which matrix STRUMPACK was
initialised from, not of which Δt is correct.

**Therefore design A is not a workaround for un-persisted state; it is a
targeted fix for a solver-birth property, which incidentally sidesteps the
un-reconstructable controller state.** "Save more state" — the obvious next
idea — is recorded here as REFUTED so nobody spends a week on it.

## 4. The fix ( design AUDITED 2026-08-19 by Codex and Grok; both rank A first )

**A — re-entry cap.** `Controller::load_memdump` must not adopt the dumped
`delta_time` as the first step. Re-enter at

    mDeltaTime = min( dumped delta_time, deck "initial timestep" )

and let the existing order-clamped growth climb back. Rationale for the bound
( Codex ): `initial timestep` is the step the deck already declares it can start
cold from; it is required and validated ( `doc/input_schema.yaml:603`,
`cl_FEM_Controller.cpp:3595` ), so no new deck key and no input-contract
obligation. On tapestack3d it is 5 ms — exactly the value measured healthy
( it.0 = 0.903, 16 Krylov, Picard-2 −138.8 dB ).

**Rejected, with reasons ( both auditors ):**
- **B** ( free + re-init after the first accepted step ) — *unreachable*: on the
  sick path the first step is REJECTED, so "after the first accepted step" never
  runs on the case it is meant to cure ( Grok ). Also re-derives the same bad
  ordering from the next 50 ms matrix.
- **C** ( seed the ordering from a mass-weighted matrix ) — *a silent no-op as
  sketched*: `STRUMPACK::initialize` leaves `reordered_ = false`, so
  `update_matrix_values` skips matching and `reorder()` runs lazily at factor
  time on whatever matrix is current ( `SparseSolver.cpp:154-162`,
  `SparseSolverBase.cpp:559-565` ). Making it real needs new wrapper API, a
  dummy solve wasting a 2.6e9-nonzero numeric factor, and it contaminates the
  frozen tiny-pivot threshold ( `SparseSolverBase.cpp:346-349` ). Research
  prototype, not a step-2 fix ( Grok ).
- **D** ( document only ) — leaves a foot-gun that has already cost days.

**Cost, stated honestly ( Grok ):** recovery is `initial × 1.2^n` growth steps,
arithmetic not measurement — ~13 steps on tapestack3d, ~49 on
`examples/helix` ( initial 0.01 ms, maximum 0.1 s ). Each still numeric-factors.
Cheap against a rejected maxit grind; not free.

**This is a documented-policy reversal and step 5 must rewrite the prose**
( Grok ): `save_memdump` ( `cl_FEM_Controller.cpp:4054-4056` ) and
`src/fem/doc/timestepping_strategy.md:186-190` currently state that the dump's
Δt is persisted SO THAT a restart does not re-enter at `initial timestep` — the
DR-38 "restart cliff". We revert the Δt half and keep the BDF-history half,
which is what actually prevents the order cliff; they are separable. Shipping
the code without the prose leaves the docs asserting the bug.

**Acceptance test ( corrected — the original was wrong for A, Grok ):**
1. Restore `memdump_2700` with NO manual deck clamp. The restart banner must
   print `delta t = 5.0000 ms`, not 50.
2. First magnetic solve: `it.0` = O( 1 ) and ~16 Krylov its, i.e. the clamped
   arm's 0.903 / 16 — the 5 ms number, NOT the live-50 ms number ( 0.019–0.043 ),
   because `it.0 = ||M⁻¹ b_scaled||` scales with the solve.
3. Step accepts in ~2 iterates.
4. **Then continue the SAME process until Δt grows back to 50 ms and verify that
   solve is still O( 1 )** — this is what distinguishes a real first-ordering fix
   from "a small Δt happened to pass" ( Codex ).
Picard-1 relative residual is NOT a gate: it already matched live while the run
was sick.

## 5. Open questions ( logged, not decided )

- **O1. A pure growth-limiter is NOT sufficient, and this is the sharpest open
  point.** `memdump_2800` has `bdf_last_dt` = 50 ms and dumped Δt = 50 ms —
  ratio exactly 1.00, so any rule of the form "limit Δt to g·h₀" permits it —
  yet that restore is the SICKEST measured ( 21245 ). So the fix cannot be
  "apply the growth limiter to the restored Δt"; it must be a genuine ramp that
  ignores the dumped Δt when re-entering. What sets the ramp's starting value?
- **O2. Why does a correct state fail on a large first step?** T2's gap. Does
  the answer change the fix? ( The mitigation does not depend on it, but a
  wrong theory could hide a second defect. )
- ~~**O3. Does the ladder reach 50 ms healthy?**~~ **RESOLVED 2026-08-19 → YES.**
  The healed lineage climbed 5 → 50 ms under the controller's own 1.2× limit
  with ZERO rejections and it.0 rising only 1.08 → 2.20 ( the healthy live
  scale ). **A ramp started small is sufficient to reach full step.**
  [ **Record hygiene, 2026-08-26 ( debt register DR-95 ).** This bullet used to
  end "T3 and T4 hold; a ramp is sufficient" — wrong, and wrong in the way that
  matters: **T3 and T4 are struck REFUTED above** ( §T3/§T4 ), and nothing here
  revives them. What the ladder measured is a NEW process *born* at 5 ms
  climbing cleanly. T3 claimed **in-process** healing — that a sick process
  cures itself after a few small steps — and measurement killed it: the sick
  50 ms restore rejected step 222 and retried at 25 ms with it.0 = 9888,
  essentially unchanged from 10851. Those are different claims about different
  objects, and the old sentence invited T3 to be quoted as a survivor. The ramp
  conclusion stands on the ladder alone. ]
- **O4. Should the save-point snap-restore bypass ( T5, `:2520-2527` ) be fixed
  in the same change?** It is a real defect against BELFEM's own documented
  table, independent of DR-92.
- **O5. Deck key or hard-coded ramp?** A key costs the full input-contract
  obligation; a constant costs flexibility.
- **O6. Codex's strongest alternative theory ( logged, NOT refuted ):** the
  defect is not "large restored first step" but *"the first post-restart
  algebraic/solver entry state differs from a live continuation"* — unpersisted
  STRUMPACK ordering/matching/equilibration, `mFieldValues`, or q-history not
  proven identical. It fits every measurement, including why a ratio-1.00
  restart can still be sick. **If this is right, a ramp is a workaround that
  happens to work because small steps are easy for any preconditioner, and the
  real fix might be to persist or rebuild the solver entry state.** The
  decisive test remains the live-vs-restored assembled-system / preconditioned
  residual comparison via the existing dump hook
  ( `cl_FEM_Controller.cpp:1845-1864` ).
- ~~**O7. Cheap partial test of O6:** restart the sick case with
  `matching : off`.~~ **RUN 2026-08-19 → INCONCLUSIVE, and the method is
  unusable on this deck.** The run aborts instead of comparing: all ranks hit
  `PETSc KSPSolve did not converge ( DIVERGED_NANORINF after 0 iterations )`
  ( journalctl -t belfem ), i.e. the magnetic solve produces NaN without MC64
  and poisons the thermal solve. This re-confirms the standing ruling that
  matching is load-bearing for h-φ; it says nothing about O6. A different
  instrument is needed — see O8.
- **O8. The evidence that DOES bear on O6, from the grow-back ladder:** the
  ladder's STRUMPACK initialised at 5 ms and REUSED that matching / ordering
  all the way to 50 ms, staying healthy ( it.0 2.20 at 50 ms ), while restarts
  that initialise AT 50 ms are sick. So the discriminator is not "fresh vs
  inherited" but **which matrix the matching was computed from**: small-Δt
  ( mass-dominated, well-conditioned ) transfers upward; large-Δt
  ( stiffness-dominated ) does not. This is consistent with every measurement
  including `memdump_2800` ( ratio 1.00, sickest ). **Status: hypothesis, not
  proven** — the equilibration diagnostics printed by STRUMPACK are identical
  in both cases ( `r_cond = 1, c_cond = 1, type = N` ), so the log carries no
  direct signal. If true, the ramp fix works for a deeper reason than
  "small steps are gentle": it forces the preconditioner to be built from a
  well-conditioned matrix.


## 6. Literature

Pierre 2022 ( Calcolo 59:36, DOI 10.1007/s10092-022-00479-0 ) gives, for
gradient flows of semiconvex F, gradient stability under a timestep restriction ( `c_F·Δt < 2β_k`; the theory doc originally wrote `β_k` and asserted 'iff' — both corrected after Codex checked the source ), computes
β₄ and β₅, and proves BDF6 is not quadratically stable ( BELFEM offers bdf1–bdf5
only, `doc/input_schema.yaml:599-617` — that specific instability class is avoided, though bdf1–bdf5 still carry their own timestep and model assumptions ). **It cannot explain
DR-92**: an absolute Δt bound cannot separate a live 50 ms step from a restored
one. It DOES support capping Δt by order, and quantifies why β_k shrinks with k.
Applying it numerically would require our coupled h-φ+thermal system to be a
gradient flow of one semiconvex functional ( it is not obviously ) and an
estimate of `c_F` for an n≈25–30 power law ( intractable, and would give an
absurdly small bound contradicted by the healthy live 50 ms steps ).
