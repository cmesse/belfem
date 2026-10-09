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

Everything else open is `[P]` (3 rows) or `[W]` (2 rows), with no untagged rows remaining.
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
| DR-199 **[CODE]** **[P]**  | fem/kernel, fem/maxwell, fem/thermal: **the `gauge` boundary condition fixes its dofs only at the first timestep, after the free/fixed split, and an untyped `gauge` block calls an uninitialized source-function pointer.** Filed 2026-10-09 on Christian's instruction ( via the Server session ) from the two code audits of the autopin fix. Established in c31ea43: a `gauge` is created like a bearing on POINT ids ( `cl_MaxwellBoundaryConditionFactory.cpp:92-100`, `cl_FEM_DofMgr_BearingData.cpp:67-89` ) but is pre-fixed by no factory — Maxwell pre-fixes `Bearing` and `Dirichlet` only ( `cl_MaxwellFactory.cpp:898-941` ), thermal likewise and then initializes ( `cl_ThermalFactory.cpp:220-236` ); its first imposition is `impose_bc` → `Bearing::impose_dirichlet` → `Dof::fix` at the first `compute_boundary_conditions` ( `cl_FEM_PhysicalBoundaryCondition.cpp:254-261`, `cl_FEM_Bearing.cpp:194-197` ), after `shift_fields()` has run the one-time `initialize()` → `split_dof_container` ( `cl_IWG_Timestep.cpp:365`, `cl_FEM_DofManager.cpp:300-302, 367-372` ); no call re-splits. Executed on a two-component strip deck with a point element: debug aborts in `assemble_jacobian` ( `cl_FEM_DofMgr_SolverData.cpp:1517`, the fixed column addressed with its free index, `cl_SpMatrix.hpp:843` ), release writes one past the CSR values there ( `cl_SpMatrix.hpp:881-913`, no check ) and then throws from Armadillo's bounds check in the fixed-value gather ( `:1952-1961` / `:2143-2153`, counts frozen by reference `cl_FEM_DofMgr_SolverData.hpp:85-90` ); two ranks: SIGSEGV on rank 1 ( workers link bearing dofs for local nodes, `cl_FEM_DofManager.cpp:1329-1391`, and flip locally after the only Dirichlet sync ). An untyped block ( `nodes` only ) crashes earlier: `mFun` has no initializer ( `cl_SourceFunction.hpp:89`, `cl_SourceFunction.cpp:72-106` ) and `impose_bc:231` calls it. Read from source only: the hanging-dof reroute ( the LIMITATION comment `cl_FEM_Bearing.cpp:144-150` ), the rank-1 crash site. Predates the autopin fix; 12a8eec did not change its reach ( a gauge was never a contact there, `cl_MaxwellFactory.cpp:3954-3961` ). The key is documented ( `doc/input_file_reference.md:912, 1027-1029, 1058, 1067-1074` ); no shipped deck uses it ( grep of `examples/*/input.conf`: none; the two `gauge` hits are a commented `coulomb gauge penalty` ) | P2 | open — diagnosis only, jury round 1 ( Codex astra/high, Grok 4.7/high, `tmp/ai_exchange/review_gauge_late_flip.md` ): both verdicts REAL for the late flip and for the untyped crash; smallest containment agreed by both: refuse a first free→fixed flip after `is_initialized()` in `Bearing::impose_dirichlet` ( as `pin_dirichlet_dof` does, `cl_FEM_SideSet.cpp:80-90` ) and refuse an absent `type` at deck read; the feature-preserving repair ( pre-fix gauge points at factory time, both factories ) needs Christian's rulings on `tHaveBearing`, on the autopin contact set and on the thermal order. No code change until Christian has seen the result | the strip deck `tmp`-only ( 18 nodes, point tag 11, `gauge { nodes : 11 ; type : constant ; amplitude : 0 ; }` ): aborts at the first assembly in debug and release, serial and on two ranks; the same deck with `bearing { nodes : 11 }` passes assembly. Not yet a test in the tree | no for the shipped decks ( abort, not a wrong answer; no deck uses the key ); Codex recommends blocking 0.9.3 on the containment, Grok calls it a product decision, Claude recommends the containment before 0.9.3 — Christian's call |
| DR-193 **[MIXED]** **[P]**  | examples/tape_quench: **the usermat plugin's property tables end below the temperatures the T5 quench deck reaches.** Active tables end at 199 K ( `k_Hast` ), 282 K ( `cp_HTS` ), 300 K ( `cp_Cu`, `k_HTS` ), 303 K ( `cp_Hast` ), 390 K ( `rho_Hast` ), endpoint held above, copper ρ clamped at 1358 K; the envelope is the minimum over every active property, 199 K, reached at ≈ 47 ms, so the 1770 K hot spot that both Picard and the hybrid reproduce is the model outside its tables. Re-filed from DR-105 residue ( 3 ) | P3 | open — ruled 2026-10-03 ( O3 ): report results against the 199 K envelope; table extension is a separate task ( NIST cryogenic and high-temperature compilations for Cu and Hastelloy C-276, DOI beside each number ) | the reconstructed T5 deck, `cmake-build-claude/tape_quench_dr105`, past ≈ 47 ms | no |
| DR-177 **[RULING]** **[W]** | homology ( closed core, report to Gregory ): **`Chain::getSimplicesMap` is declared as `OrderedMap< id_t, int > &` but returns the member `OrderedMap< index_t, int > mSimplicesMap`** ( `cl_Chain.hpp:44`, `:107-108`, `:157-160` ). The two typedefs coincide on the default build, so it compiles; under `BELFEM_INT64` ( `USE_MKL_64BIT_API`, `CMakeLists.txt:155`, `typedefs.hpp:47` ) `index_t` widens and the return is ill-formed, and `Cochain::operator()` ( `cl_Cochain.cpp:76` ) plus `tests/homology/test_Cohomology.cpp` call it. Not an AI edit: the six core units are closed | P3 | **open ( filed 2026-09-22 from the test-gap Tier 1 run )** — Gregory's call which type the map keys are; the fix is a one-word signature change in a closed file | configure with `-DUSE_MKL_64BIT_API=ON` and `make belfem_homology` | no — an opt-in build |
