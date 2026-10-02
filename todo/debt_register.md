# Debt Register

**Purpose:** the live set of open items across all campaigns, with the blocking-1.0 column as
the release lens. Retired rows move to [debt_register_closed.md](debt_register_closed.md). The
canonical examples cited below now live there.

## Freeze lens

Every open row now carries a second tag beside its `[CODE]`/`[RULING]`/`[RUN]` kind:

- **`[F]`** — settle **before the design freeze**. Reserved for rows whose resolution changes a
  **user-visible contract**: a default, an input key's meaning, which examples ship, or a
  documented convention. After the freeze these become expensive or breaking to change.
- **`[P]`** — post-freeze. Internal defects, performance, test coverage, verification runs.
  A freeze locks design, not bug fixing, so these are safe to carry.
- **`[W]`** — watch. Latent hazards with no current trigger and no reachable caller; no action
  is owed unless something starts firing.

**The `[F]` list is EMPTY.**

Everything else open is `[P]` (3 rows) or `[W]` (3 rows), with no untagged rows remaining.
Nothing in those two groups needs a decision this week; they need time, which is the one thing
the freeze does not take away. The per-recount history is in git.

**Maintenance (protocol §11):** a session that resolves or creates an open item updates this
register in the same session. IDs are stable and are never reused. **Prefer extending an existing
row over filing a new one ( policy set 2026-08-29, Christian ):** when a session's own gate or audit
turns up a finding in the same contract or subsystem as the row being worked, widen that row's scope
cell and put the work in its plan, rather than branching a fresh ID. New rows are for genuinely
separate subsystems. The register had been inflating one ID per finding.

What follows is the operating manual for the register. It is not history — the per-pass narrative
this preamble used to carry is in git, which is what git is for.

## How to read a row

A ~~struck ID and description~~ means the row is **retired** — no design/code/ruling work remains.
The status column stays unstruck because closure evidence, surviving gates, and waivers live
there. The strike itself is not verification evidence: **struck is not verified.**

When a row's description and status cell disagree, **the status cell is the verdict** on whether
the item is open — but a stale description is still currentness work, not permission to read past
it. A dead clause inside a still-live row is struck in place. **The row supersedes examples in this
preamble:** where an example and the row conflict, the row wins.

Writing `closed` in the status cell and striking the ID and description are **two separate
actions**. A row marked closed in status but never struck is live work for the next sweep.

## When to strike

Strike the ID and description only when the work is finished with **no residual**. A row stays live
when it is merely downgraded, part-done, or still waiting on a decision, a design question, or an
unwritten test. A residual living in another session or artifact is still a residual.

**Exception, Christian's call — the DR-42/DR-49 pattern:** when the *design* work is finished and
the residue is a single run, Christian may strike the row anyway. Record the exception; the gate
survives in the status column, and **struck is not verified** applies with full force. Christian
may instead retire the gate *with* the row — write that in the status cell, never infer it.

A static sweep may recommend a strike, or a by-design / already-fixed / nothing-actionable closure;
it does not issue the ruling. **"Pending run" is not a row a static sweep can close.**

When a row's amendment history outgrows human readability, **compact it** and say so in the row.
Git carries the generations; the register carries the current state.

## How a sweep gathers evidence

**For what the source does today, the evidence is the tree — never the plan and never the devlog.**
Verify on disk, not by inference. Rulings and waivers are exceptions: the recorded decision is the
evidence.

If the source still matches an unstruck row's claim, leave it. Classifying what rows are *waiting
on* from their status columns is a different pass from a currentness sweep — useful, but say which
one you ran, and never report one as the other.

A stale row is neither reliably pessimistic nor reliably optimistic. Sweeps have found debt already
paid, debt **larger** than recorded, and written claims that were simply false. If a claim is
refuted, rewrite what survives — it is often narrower and different. Do not close a false row and
lose the real defect underneath it.

Uncommitted tree state is not closure. Commit state can change a row's *argument* — unreferenced
dead code rather than abort-shadowed — without changing its severity.

A cleanup row that is ready to execute pins its exact change set. Before closing a plan, lift its
leftovers into this register. Detailed residuals stay in the owning plan; this register is not a
second copy of it.

## Closure evidence (protocol §11)

Use this ladder, weakest first, for register evidence: **AI agreement → single-voice static source
trace/review → standalone probe → clean compile/link → named suite gate → named
production/end-to-end run that can observe the row.** A strike does not move a row up this ladder;
only the status evidence does. Say **reviewed** for static-only work and **verified** only when an
executable gate passed.

**Agreement is the weakest tier.** In this register, it has repeatedly done no work: promote a
finding on a re-read of the tree, not on a vote. Flag single-raiser rows as such, and keep
confidence in a mechanical file:line trace separate from confidence in any arithmetic riding on it.

**Pick a gate that can observe the defect.** For the canonical residual-invariant Jacobian/history
case, use Newton-iteration count, not a physics-output comparison.

**A status cell records what a run did *not* cover** as well as what it did. A waived gate is named
there.

Honest closures can be something other than "the test went green": already fixed in tree, stale or
moot decision, by-design, or out-of-scope. Record which one — never dress a decision as run
evidence.

A core-controller termination change is not applied on a static finding. It takes plan+audit →
code+audit with both vendors, plus an A/B.

## Standing traps

1. **A grep-shaped absence is not a functional gap.** Prove the thing is unavailable, not merely
   absent from the one function you read; trace the alternative routes first. Canonical: DR-19:
   the connector record is deliberately absent from the shell record because the coating block
   types and wall-to-facet links ship by two other explicit channels. That false positive survived
   a jury round, a row edit, and two independent reads before the code refuted it. Retract a false
   positive in the row; do not quietly delete it.
2. **Follow the ID you cite.** A plausible ID copied forward without being followed is the same
   defect class this register exists to catch. The always-active `set_currents` guard that broke
   every parallel run is **`61bbfa11`** — not DR-24, the periodic free-cut `mPrescribedCurrents`
   row. That wrong citation stayed in this register for two passes.
3. **When a row names N sites, count them in source.** Canonical: DR-45 — three `INFO` sites, not
   two, with the Matrix-RHS `solve` overload still gated on `rank() == 0` after the commit that
   supposedly fixed it.
4. **A commit can pass its own acceptance test and still leave the neighborhood wrong** — an
   adjacent input-contract key, a comment, or a document. The acceptance test clears the change,
   not what sits beside it.
5. **A source comment that asserts as fact what every other artifact in the same commit correctly
   hedges is itself debt.** The register row is allowed to be the better-calibrated statement,
   and has been.

## Triage

Every live row waits on a **run**, a **ruling**, or **code work** — or a named mix. Classify before
planning: a ruling costs a sentence from Christian, runs batch together, and only the code-work
residue needs a full audit round.

Each live row carries that classification as a tag in its ID cell:

| tag | meaning |
|---|---|
| `[RUN]` | the only residue is executing a named gate. Inputs exist, no code, no decision — autopilot-ready |
| `[RUN-BLOCKED]` | a run gate whose inputs are gone or whose upstream is unfixed. **Not** autopilot-ready; the row says what died |
| `[RULING]` | needs a sentence from Christian — a strike, a commit, a scope call. No machine time |
| `[CODE]` | code remains, so it takes the standing plan+audit → code+audit round |
| `[MIXED]` | a named combination; the row states the order |

The tags are a **classification pass from the status columns as written** — not a currentness
sweep, and no source was re-verified to place them. Re-derive a tag whenever a row's status
changes; a tag is the cheapest thing in the row to let go stale.

**`[RUN-BLOCKED]` is the tag that earns its keep.** Three gates have already died because the
dump or the deck they named was deleted between the fix and the run. When a row's gate depends on
an artifact, name the artifact in the status cell and say it must be kept — see the evidence rule
above, and expect the default to be deletion.

| ID                          | area | severity | status | reproducer | blocking-1.0? |
|-----------------------------|---|---|---|---|---|
| DR-29 **[CODE]** **[P]**    | homology: coreduce worklist port (#1/#2a) + stale-adjacency cleanup (#10) — sequenced after rectification guards/tests | P2 | deferred — **confirmed open 2026-08-13, and #10 is more specific than recorded.** ~~The live path `coreduce_complexCCR()` → `pGeneralizedCocombine(p)`~~ **live-path clause corrected 2026-08-30: production is `coreduce_complexPellikkaGeneralized()` (`cl_MaxwellFactory.hpp:48` default + `mSuggestHomologies` hardcoded true `cl_CutFactory.hpp:91`), which re-enters `pCoreduce(p+1)` after each cocombine phase — CCR is opt-in.** `pGeneralizedCocombine` (body `:1119-1228`) is a single forward pass with an inner backward scan — no worklist, no queue, no neighbor re-enqueue. The Pellikka path `pCoreduce` (`cl_SimplicialComplex.hpp:338-400`) is the older full-rescan form. **#10 pinned:** the stale-adjacency cleanup blocks (`:1175-1196`) sit inside `if ( tNeighborIDs.size() > 0 )` (`:1175`), yet `remove_kcochainFromMap` runs unconditionally (`:1223-1224`) and only deletes the cochain and erases the map entry (`hpp:320-332`) — so when there is no valid neighbor the removed IDs stay referenced in adjacent boundary maps. **Three-way analysis round 2026-08-30 (Claude + Codex sol/xhigh + Grok 4.6/xhigh, blind jury, reconciled):** #10 confirmed as an elementary coreduction pair whose bookkeeping is skipped; a production-reachable abort cascade was constructed independently by both auditors (3-D, zero-neighbor skip leaves dead face in a tet boundary → later merges shrink it to a dead raw singleton → `pCoreduce` raw `size()==1` picks the dead ID → `Map::operator()` abort `hpp:369`, all build types) — order-dependent, the "every change breaks it" mechanism: safety today is by phase ordering only, unchecked and undocumented. **Refuted and must not be "fixed":** the Smith matrices from `createMatrixFromCoboundaryMap` are CORRECT — δ∘δ=0 makes the live-restriction exactly the retracted operator (Codex derivation, independently re-verified); do not mirror merges into (p−1) coboundary rows and do not switch the matrix source to boundary maps. Fix plan (design only, nothing touched): (0) debug invariant — live↔live row/column agreement + mD·mD=0; (1) #10 cleanup unconditional AND find()-guarded incl. the `(p+2)` walk `:1190-1194`; (2) `pCoreduce` pair test live-and-unit (Mrozek §4 ±1 check missing); (3) worklist port last, FIFO+membership set, dequeue revalidation, component-aware omit seeding — `pCocombine` is NOT a template (LIFO `Cell::pop`, sorting `unique()`, unguarded pop). Report for Gregory: session 2026-08-30, `devlog/dl20260830_dr29_coreduce_analysis.md`. **HANDED OFF 2026-08-30 (Christian): the report went to Gregory Giard, who develops this module and decides what happens to it. The row stays DEFERRED pending his call. The fix plan above is a proposal for him, not a queued task — no BELFEM-side session should self-start it, and the analysis does not need redoing.** **AI EDIT BAN 2026-08-31 (Christian + Gregory):** the cohomology core (`cl_Cohomology`, `cl_Homology`, `cl_SimplicialComplex`, `cl_Chain`, `cl_Cochain`, `fn_Smith`) is closed to AI edits regardless of session edit approval — `doc/ai_collaboration_protocol.md` §7.1. This row is read-and-report only; steps (0)-(3) are executable by Gregory, or under his named authorization, and by nobody else. | Betti/cut-count comparison; step-0 invariant green on annulus×3 algorithms + both periodic 3-D fixtures, debug AND release — gates apply only if Gregory adopts the plan | no |
| DR-71 **[MIXED]** **[W]**   | homology **the `mSuggestHomologies == false` branch of `CutFactory::compute_cohomologies` is unreachable and broken in two independent ways.** Found 2026-08-13 while building the DR-23 layer fixture. The flag is hardcoded `true` (`cl_CutFactory.hpp:91`, comment: "this needs to be turned on for the incidence matrix to work!"), so the `reduce_complex*` calls inside `if (!mSuggestHomologies)` never execute in production. When forced: (1) running reduce AND coreduce together zeroes the SNF input — reduce guts the chain side (all top cells omitted), coreduce guts the cochain side, and `createMatrixFromCoboundaryMap` then sees only zero maps, so every surviving 1-cochain becomes a "generator"; (2) `pGeneralizedCombine` carries an unbounded iterator walk (`cl_SimplicialComplex.cpp:873-878`): the `while` advances `it2` with no `end()` bound and the `if (it2 == tSimplicesMap.end()) continue;` below it is unreachable-as-guard — when no boundary simplex qualifies, `++it2` walks past `end()` (UB; observed as a 99%-CPU infinite loop). The sibling cocombine loop (`:1139-1148`) IS bounded (breaks at `begin()`), so the defect is one-sided. NOT triggered from any production path today. Per-symptom reproducers: **(1) reproduces on any well-formed mesh** — force `mSuggestHomologies = false` so reduce and coreduce both run, and the SNF input zeroes by composition; **(2) was observed only through the malformed complex** (the DR-23 `aComputeConnectivities = false` fixture bug, which empties every edge coboundary so no boundary simplex ever qualifies and `++it2` walks past `end()`) — on a well-formed complex the loop may happen to terminate, but the missing bound and the dead guard are structural either way. Options when touched: bound the loop like its sibling, or retire the dead branch outright. **SPLIT BY THE AI EDIT BAN 2026-08-31** (`doc/ai_collaboration_protocol.md` §7.1): defect (2), the unbounded `++it2` at `cl_SimplicialComplex.cpp:873-878`, is INSIDE the banned core and is Gregory's alone; defect (1) and the "retire the dead branch" option live in `cl_CutFactory.cpp:278-281,295-298,329-332`, which the ban does NOT cover. Do not treat this row as wholly closed to AI work, and do not treat the banned half as actionable | P3 | open — latent; no production path reaches it | (1): any mesh with `mSuggestHomologies` forced false; (2): the `tests/homology/` annulus built with `aComputeConnectivities = false`, then `reduce_complexCCR` | no |
| DR-105 **[MIXED]** **[P]**  | examples/deck-envelope: **the tape_quench deck is tuned and its quench-front abort has a landed remedy; three things are left — a long rerun, an unexplained parallel wall, and one physics ruling.** Split out of DR-02 at its closure ( Christian's ruling, 2026-08-24 ). The "µs-scale Δt contraction radius" it was filed for turned out to be a COLD-START artifact: the same deck that rejects a cold 50 µs first step grows into 20 µs steps once BDF5 history is built, so ramp-in matters and step size per se does not. The two physical drivers are the 100× Gaussian jc-suppression defect under n≈37 and the measured current table's first jump at t≈0.77 ms | P2 | open — **recipe applied, remedy landed, rerun owed.** T5 recipe ( `initial timestep : 0.002 ms`, `maximum timestep : 0.02 ms`, strumpack magnetic / petsc thermal, Anderson on, Newton via the 1e-4 switch, chi off ) came out of a six-variant matrix at equal 12-15 min wall budgets: T5 reached 0.954 ms and rode the cap with near-zero rejections, straight through the 0.77 ms wall, while T6 found the natural ceiling at ≈16-30 µs with rejected probes at 31 µs — capping just under it beats letting the PID hunt ( 4.8 vs 3.2 ms/hour ). **The serial 9.8 h discriminator settled the np=8 question and then hit a new defect:** serial ran to t = 9.14 ms at ~1.0 ms/h — 24× past the np=8 wall — so that wall is a parallel effect, still unexplained; the run then aborted at the quench front on `rho_piecewise`'s `n > 1` assert ( `powerlaws.hpp:780`; the floor-retry abort above it is `cl_FEM_Controller.cpp:2700-2706` ) inside the 90 K < T ≤ 92.5 K window, where the plugin's own fits have fallen back to 1.0 but `T_crit` keeps the superconducting branch alive. Debug-only guard; a release build would divide by zero in the flux-flow branch. **That abort now has a landed remedy:** the 2026-08-29 usermat round single-sourced `hts_Tcrit`, made the fit cutoff inclusive ( `T <= hts_Tcrit` ) and added a deck-side `critical temperature` key read in `MaterialFactory` ( `cl_MaterialFactory.cpp:283` ) that overrides what a plugin's `_init` set; the shipped deck carries `examples/tape_quench_usermat/input.conf:76` = `critical temperature : 90.0 K ;`. **Residue: (1) the long rerun** — 0.3 s ≈ 60 h serial at measured throughput, MPI ranks would shorten it; **(2) the np=8 residual floor** ( re-measured 2026-08-31 from the surviving logs; **the earlier "floor-collapse at step 264" framing was wrong and is retracted** ) — nothing happens at step 264. np=8 was pinned at the 1 µs floor essentially throughout: 229 floor-retry warnings over 264 accepted steps ( 87% ) against serial's 47 over 4191 ( 1.1% ), and serial sits at the floor at its OWN step 264 too ( Δt = 1.1 µs ). Step 264 is merely where 20 consecutive retries ran out ( `cl_FEM_Controller.cpp:2267` ). The real discriminator is the achievable magnetic residual at comparable Δt: serial descends 1.5e-3 → 4e-7 in 9 Picard iterations, while np=8 sits FLAT at 1.6-2.2e-4 through every iteration of all twenty retries — ~200× above the deck's 1e-6 tolerance, insensitive to both Δt and relaxation, i.e. a ~500× residual-floor difference from rank count alone. ( Full traces and the thermal machine-zero second symptom: `devlog/dl20260831_riva_n_guard_and_dr105_wall.md` §3. ) **One confound is now ELIMINATED, not merely doubted:** both runs used the identical binary ( Build Date `Aug 27 2026 03:25:34`, commit `dcca23cc`, DEBUG on, `hphiTrun` ), verified from both log headers, so "debug binary" is common-mode. Deck identity still rests on the 08-27 devlog's record rather than on the artifacts, because `belfem.log` does not echo solver settings. Mechanism remains UNDIAGNOSED — ranked hypotheses, not results: rank-dependent residual assembly / MPI-reduction inconsistency over owned-vs-ghost dofs ( 2e-4 is ~1e9× machine epsilon, far too large for roundoff on a 2133-node problem ), STRUMPACK's distributed factorization path, thin-shell / cut partitioning at 267 nodes per rank, and Anderson depth 3 ( ranked LOW — the residual is already 2e-4 at Picard iteration 1, before Anderson has history ). Bearing CLEARED as a suspect ( `link_bearings_with_dofs`, `cl_FEM_DofManager.cpp:1472`: broadcast table, every rank links wherever the node exists ). Cheap discriminators, ~15 min each since np=8 walls fast: `compute conditioning : true` at np=8 ( the abort message's own hint — a modest κ rules conditioning out and points at assembly ), np=8 with `library : mumps`, and np=2 / np=4 to see whether the floor scales with partition count or is a cliff. Log at `cmake-build-claude/tape_quench_dr105/belfem_np8_wall0375ms.log`; **(3) whether to soften the defect or the current table at all**, which is physics and Christian's call — and T4 already REFUTED the obvious lever there: the 175 A ramp-BC variant is a STRONGER drive, not a table-jump control. Compacted 2026-08-30; residue (2) re-measured and its framing retracted 2026-08-31 | **Two different decks share this row's name — do not confuse them.** The T5 envelope lives ONLY in the reconstructed `cmake-build-claude/tape_quench_dr105` ( strumpack/petsc, `anderson depth : 3`, timestep 0.002 → 0.02 ms ). The shipped `examples/tape_quench_usermat` is NOT the T5 deck and, as of 2026-08-31, is no longer even the same MODEL: it was rebuilt onto a single-tape mesh ( `tape.geo`, replacing the two-tape `tapestack3d.geo` ), moved to `resistivity type : riva`, sidesets 1,2 and one in/out curve pair, mumps at `:16`, timestep 0.5 → 1 ms at `:53-54`, still no Anderson. Running it discharges no residue of this row. **Two consequences for the row's diagnosis:** its `riva` selection removes the `n > 1` abort CLASS independently of the `critical temperature` remedy — `rho_riva` has no such precondition anywhere and `riva_rho_pl` takes n ≤ 1 as the finite ohmic limit `ec/jc` ( `powerlaws.hpp:2481-2512` ) — so this row now has two independent remedies where it records one; and its defect rewrite ( depth 0.25 → 0.1, re-sized against an explicit ~0.95 mm MPZ argument, switched on at 20 ms with a C2 quintic inside the current table's 5-155 ms plateau ) is residue (3)'s physics ruling TAKEN IN CODE, addressing both named drivers — though on the new model, not by softening the T5 deck. **CURRENTNESS 2026-08-31: the reconstructed deck is REPAIRED and residue (1) is unblocked.** It previously aborted at setup on its `custom { }` material blocks ( the factory hard-errors on that spelling since the 08-29 rename — a deliberate diagnostic, NOT an alias, `cl_MaterialFactory.cpp:252-262` ); all six material blocks now read `usermat` and `ybco1` carries `critical temperature : 90.0 K` ( `input.conf:66-72` ). It remains `resistivity type : piecewise`, so its quench-front protection is the `critical temperature` key alone. **The 2026-08-31 riva n-source gate is INERT for this deck** — every arm is conditioned on `ResistivityLaw::Riva` ( `Material::check_riva_n_source`, `cl_Material.cpp:790` ), so a piecewise deck takes none of them and nothing in that change can alter this run's behaviour. Serial evidence `belfem_serial.log`, parallel evidence `belfem_np8_wall0375ms.log`, both in the deck dir | no |
| DR-113 **[CODE]** **[P]**   | comm: **the point-to-point fabric has no structural defense against the DR-100 class.** Every exchange between a rank pair shares two tags ( `comm_tag` base for sizes/scalars, `+1` for payloads ), so correctness rests entirely on both ranks executing their exchanges in the same order — one skipped or surplus message silently poisons the NEXT recv on that tag ( truncate, garbage, or hang ), with the crash surfacing in whatever innocent collect runs later, an arbitrary distance from the producer. DR-100 took a month from defect to first observed crash and two probe campaigns to walk back to the producer. Architectural hardening options, in rising cost: ( a ) a debug-build drain check ( `MPI_Iprobe` sweep like the DR-100 probe ) at exchange-boundary chokepoints that aborts loudly at the FIRST stray, naming source/tag/size; ( b ) per-exchange tag separation ( widen `comm_tag` with an exchange-sequence component mod the tag limit ) so a stray can never match a later exchange's recv. The contract itself is now documented on `comm_tag()` ( 2026-08-27 ) | P3 | **gate ( b ) VERIFIED 2026-09-29 — Christian's `make check` 17/17 in one shared tree ran `commmpi_np2` and `commmpi_np4` with the DrainCheck tests in the binary ( they landed in `52465d8`, 2026-09-23 ); the other shared tree not run. Gate ( a ) still owed, blocked on the Tier 2 FEM fixture ( `todo/test_gap_backlog.md` F1 / T-62 )** — `DrainCheck.CleanFabricIsSilent` and `DrainCheck.StrayScalarFiresOnEveryRank` in `tests/comm/test_CommMPI.cpp`, registered through the existing `commmpi` suite at 2 and 4 ranks ( `tests/comm/CMakeLists.txt` ): a clean ordered exchange leaves the check silent; an unmatched one-int `send` makes every rank throw ( debug policy, collective fold ), rank 0 drains it, a final check is silent; the stray test skips in release where the check compiles out. Plan jury ( Codex terra/high + Grok 4.6/high ) and code jury ( Codex terra/high + Grok 4.7/high ) both converged; the code round caught the tests landing before the FIRST `#endif // BELFEM_MPI` instead of the last, fixed. Verified 2026-09-29 as recorded in the status head; the second shared tree remains unrun. Residuals recorded: the unmatched send relies on Open MPI eager completion ( not an MPI guarantee; a rendezvous transport fails into the ctest TIMEOUT, never a silent pass ); the payload tag `+1` is exercised but not asserted, since no wrapper can post a payload-tag stray without its size message. Gate ( a ), a chokepoint at more than one rank, needs the Tier 2 FEM fixture ( `todo/test_gap_backlog.md` F1 / T-62 ). Previous: **option ( a ) landed 2026-08-31, gates owed** — `comm_drain_check()` in commtools ( barrier → `MPI_Iprobe` sweep of both pair tags → allreduced verdict → `BELFEM_ERROR` on EVERY rank, throw policy intact ) at five chokepoints: `DofManager::solve`, `solve_from_residual`, `postprocess`, `SolverData::collect_matrices` entry, `Controller::finalize` exit. No-op in release/serial; finalize is now collective in assert-active builds. Codex + Grok code audits accept at xhigh ( `tmp/ai_exchange/dr113_fabric_hardening.md` ). Gates: build + `make check` GREEN 2026-08-31, but that does **not** exercise this check — only `tests/comm` registers `TESTRANKS` ( 2 4 ) and that binary never builds a `DofManager`/`Controller`, while every other suite runs at one rank where the function returns at its `size() < 2` guard, so **no test reaches a chokepoint with more than one rank**. Still owed: a parallel gate deck ( RLC np=4 / 2D_Tapestack np=2, debug ) with zero drain aborts, and a deliberate-stray run proving the tripwire fires. Residuals: `comm_check` throw mid-sweep is not collective ( broken-MPI path ); a stray born inside `collect_matrices` before `collect_rhs_vector` is not covered. Option ( b ) tag separation stays deferred by design — touches every comm call site, `MPI_TAG_UB` care needed | any future stray now aborts at the nearest chokepoint naming source/tag/size; the check IS the archived DR-100 drain-probe pattern made permanent | no |
| DR-122 **[CODE]** **[W]**   | homology: **order-2 cut-pipeline residue, split from DR-103 at strike time — three known gaps, all unreachable behind the deliberate order-1 hanging-edge guard** ( `BELFEM_ERROR( max_element_order() == 1 )` in `MaxwellFactory::create_hanging_edges_and_facets`, `cl_MaxwellFactory.cpp:1413` ): ( 1 ) the diagonal thick-cut cases ±5/6/7 rely on chord-midside neighbor coverage that `src/homology/doc/thick_thin_cuts_and_conjugate_edges.md` documents as unenforced; ( 2 ) `CutSet::create_duplicates` sets the duplicate's node source to the original unconditionally — no cascade when the original is itself hanging ( unlike `InterfaceProcessor::duplicate_nodes` and the ThinShellFactory ), which at order 2 trips the `is_hanging` assert at `cl_FEM_DofMgr_DofData.cpp:3513` in debug and writes a wrong T-matrix row in release; ( 3 ) no TRI6/TET10 homology fixture exists. NOT in scope: the former `fix weights for second order` todo — ruled stale at DR-103's strike ( unit nodal weights are exact at any Lagrange order ) **( 4 ) WIDENED 2026-09-18 ( Codex, DR-133 round ): `CutFactory::close_terminal_loops()` builds the duplicate segment with `insert_node( tNode, n-1-k )` ( `cl_CutFactory.cpp:1598` ), which for a LINE3 `[e0, e1, mid]` yields `[mid, e1, e0]` — the midside node lands in an endpoint slot — while the same function's arclength code at `:1688` treats node 2 as the center. Unreachable today behind the same order-1 guard; wrong the day order 2 is admitted** | P3 | **open ( filed 2026-08-27 at DR-103's strike )** — watch only: no reachable trigger while the order-1 guard stands; flips to live work the day order-2 h-φ is scoped in | order-2 inductor mesh ( `scratchpad/dr31_order2/` ), once the order-1 guard is lifted | no — order-2 h-φ is not a 1.0 feature ( DR-103 ruling ) |
| DR-177 **[RULING]** **[W]** | homology ( closed core, report to Gregory ): **`Chain::getSimplicesMap` is declared as `OrderedMap< id_t, int > &` but returns the member `OrderedMap< index_t, int > mSimplicesMap`** ( `cl_Chain.hpp:44`, `:107-108`, `:157-160` ). The two typedefs coincide on the default build, so it compiles; under `BELFEM_INT64` ( `USE_MKL_64BIT_API`, `CMakeLists.txt:155`, `typedefs.hpp:47` ) `index_t` widens and the return is ill-formed, and `Cochain::operator()` ( `cl_Cochain.cpp:76` ) plus `tests/homology/test_Cohomology.cpp` call it. Not an AI edit: the six core units are closed | P3 | **open ( filed 2026-09-22 from the test-gap Tier 1 run )** — Gregory's call which type the map keys are; the fix is a one-word signature change in a closed file | configure with `-DUSE_MKL_64BIT_API=ON` and `make belfem_homology` | no — an opt-in build |
