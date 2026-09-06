# Lessons Learned {#doc_lessons_learned}

**Date:** 2026-08-20
**Purpose:** Trigger-indexed operating rules distilled from this project's own failure record, so that a session about to repeat a known mistake is stopped at the moment it matters
**Module:** cross-cutting

Distilled from 298 devlog entries and the debt register spanning 2026-03-18 to 2026-08-19,
which yielded **537** clustered incidents in 20 clusters. The catalog itself is 539 rows: the
clustering and the cards were locked as a pre-registration before INC-538/539 were added and are
deliberately not retrofitted (`lessons_learned_evidence.md:873-877`). Addenda through 2026-08-31
(INC-540 … INC-566) are appended in the evidence file and likewise not folded into the counts. Every rule cites the incidents
that paid for it. The `INC-NNN` identifiers index
[the incident evidence](lessons_learned_evidence.md), where every row cites the dated `devlog/`
entry it was mined from.

**Status:** proposed. The rules are written and cited; two now carry a mechanical check —
L-08 partially, through `scripts/check_doc_claims.py`, and L-21 through
`scripts/check_wrapper_policy.py`. The rest are not yet enforced, and each card names the
enforcement target it is waiting on.

---

# LAYER 1 — THE TRIPWIRE LAYER

*This is the part a session loads. Budget: 250 lines. It ends at the LAYER 2 marker.*

## How to use this

The error rate is not the target. **The error survival rate is.** Errors will recur; this
file exists to keep each one's survival time near zero. A lesson that cannot fire mid-session
is decoration.

Two ways in, and you will normally use one of them, not both:

- **Starting an activity?** Go to the *trigger index* and read the two or three lines under it.
- **Staring at a symptom?** Go to the *failure signatures* table and take the discriminator.

Cards in Layer 2 carry the evidence and the reasoning. You do not need them to act.
This file does not explain *why we work this way* — that is `doc/ai_workflow_best_practices.md`,
and it is not restated here.

## Vocabulary

The evidence ladder is `doc/ai_collaboration_protocol.md` §11, strongest first: end-to-end
reproducer, focused regression, compile/link, numeric probe, static source trace, literature
consistency, **AI reviewer agreement — weakest**.

- **verified** = an executable gate ran. Name the gate. Nothing else earns the word.
- **reviewed** = a static read is complete. Most work is this. Say so.
- Concurring audits do not lift a claim past a source trace. Three voices agreeing that code
  is live has been wrong here repeatedly; one smoke run settled it in seconds.
- **A build passing is not a verification of anything but the build.**

## Iron rules

1. The evidence is **the tree** — never the plan, never the devlog, never the register row.
2. "Verified" requires an executable gate that ran; agreement ranks below a source trace.
3. Every number, ID, hash and line anchor is re-opened at its source in the session that
   cites it. Recognition is not lookup.
4. Every universal or negative claim carries its enumeration or its search scope in the same
   sentence, or it is not made.
5. **Presence of code is not evidence of execution.** Prove the path runs before reasoning
   about it.
6. On a numerical regression, establish what the number *measures* before touching the physics.
7. Endangered artifacts are copied by the agent proposing the destructive step, in the same
   message that proposes it — never delegated.
8. **Struck is not verified.** Retiring a design does not retire its run gate.

## Trigger index — by activity

**Before designing anything**
- Two viable designs still standing, or the question is physics → stop, route to Christian. `L-13`
- Building a model from one sample, one deck, one metal, one mesh → widen before it sets a
  threshold or a severity. `L-14`
- Extrapolating an anchor from fitted data → compare the target tolerance to the spread of the
  inputs the extrapolation must cancel; if the inputs' spread dwarfs the target, MEASURE the
  anchor instead. Two designs died this way in one night (INC-553). `L-14`
- State the degenerate-point physics of the domain boundary before choosing a functional form —
  at |B| = 0 there is no field direction, so nothing may stay angle-dependent there; the
  constraint was free and it outranked the fit both times. `L-13` (INC-553)

**Before changing a default, a constant, or a dispatch table**
- A default change needs a **run gate, not a review gate** — on a deck that omits the key. `L-04`, `L-02`
- Diff **default member initializers and class-level constants** as their own explicit step
  when hunting a regression; call-path reading does not surface them. `L-02`
- Changing a shared convention (orientation, node order, master/slave, argument order)?
  Enumerate every consumer; all-`real` signatures let stale call sites compile silently. `L-06`
- Audit a default against the DOCUMENTATION'S taught examples, not only in-tree callers — a
  guide's `class X : public Base` snippet is a consumer, and in-tree safety by accident of a
  shadowing derived default is one derivation away from a silent wrong answer (INC-552). `L-06`

**Before calling a third-party library — MPI, LAPACK, PETSc, MUMPS, HDF5, anything vendor-prefixed**
- NEVER directly. Only through the dedicated wrapper layer (`commtools.hpp`, `fn_*` LAPACK,
  the `Solver*` classes, io). A missing operation extends the wrapper first. Hard policy,
  Christian 2026-08-30. The wrappers carry the 32/64-bit and lifecycle contracts a raw call
  re-derives wrongly. `L-21`

**Before editing the cohomology core — six units in `src/homology/`**
- `cl_Cohomology`, `cl_Homology`, `cl_SimplicialComplex`, `cl_Chain`, `cl_Cochain`, `fn_Smith`
  (`.cpp` and `.hpp`) are CLOSED to AI edits. **The rest of the directory is not** — CutFactory,
  CutProcessor(Manual), CutData, CutSet, InterfaceProcessor, BeltedTree, Topology,
  `en_CutAlgorithm`, `CMakeLists.txt` and the module's `doc/` follow ordinary rules. Session edit
  approval does not reach the six; lifting needs Gregory Giard's named authorization **relayed by
  Christian and recorded in a devlog with scope and date**. Hard policy, Christian + Gregory
  2026-08-31, `doc/ai_collaboration_protocol.md` §7.1.
- The reason is **novelty, not difficulty**: a *modified* Pellikka reduction, unpublished, so the
  nearest prior is textbook Pellikka and every deliberate modification reads as a defect. Read,
  report to Gregory, do not repair. `L-22`
- The attested case: `coreduceOmit` (`cl_SimplicialComplex.cpp:1236-1257`) removes ONE 0-cochain
  per stall and re-runs the whole `pCoreduce` cascade each time. Batching it — it removes them
  all anyway — was proposed and struck as unsafe: it strands 1-cells whose boundary is empty, and
  `pCoreduce` only ever pairs a cell whose boundary holds **exactly one** entry
  (`cl_SimplicialComplex.hpp:364`), so they survive into cocombine as spurious generators.
  **The loop looks trivially hoistable and is not.** `L-22`
- Generalize the test, not the file list: an algorithm the models have never seen inverts the
  usual prior. **Do not treat concurrence as coverage** — the errors are correlated, so a
  majority verdict is worth no more than one voice (2026-08-30: 2 of 3 read a deliberate
  invariant as a defect; the third refuted them). A hard but widely published algorithm is safer
  ground than a simple unpublished one. `L-22`, `L-03`

**Before adding an assert or a `BELFEM_ERROR`**
- List the legal states it must accept: empty containers, worker ranks, NaN, degenerate
  geometry, `n == 4`. A new guard is a new failure mode. `L-16`
- Retryable state → return status to the controller. Wrong-answer state → hard-fail in
  release. Never silent, never abort what the controller could have rescued. `L-18`
- A DIAGNOSTIC (conditioning estimate, eigen probe) may never abort the run it measures —
  degrade to "n/a" with a strike counter, as the ARPACK path does (INC-555). `L-18`

**Before deleting something as redundant or dead**
- Name what would break if it were load-bearing, then check that specific thing. `L-17`

**Before touching restart, serialization, or any save/load path**
- The first step after a restore is a **different program**. Seeding, first Δt, integrator
  history and solver birth state are each restored or explicitly re-derived. `L-12`
- A cache is validated against everything that determines its contents, not against its
  input file. `L-02`

**Before claiming a path works**
- Prove it executes: a work count, a probe, or a deliberate break that must fail. `L-09`
- "The run completed" is not "the path ran" — check the deck actually enables it. `L-09`, `L-01`
- Executing is still not working: a guarded fallback can run the new path and return the OLD
  answer, silently. The probe must assert the new EFFECT differs from the old behavior — the
  one probe written that way caught a fix shipping 100 % inert (INC-549). `L-01`, `L-16`
- Consuming foreign data (vendor columns, Python-written HDF5, embedded schemas)? Read the
  contract from the artifact itself — declared metadata, the dataset's own type — never assume
  the convention. Units (A/cm read as A/m) and string csets (UTF-8 read as ASCII) each burned a
  session the same evening (INC-549, INC-550). `L-02`, `L-06`

**Before merging**
- Serial correctness says nothing about ownership. Which rank owns each entity, and what does
  a non-owner do? `L-11`
- Both matrix backends, if the change touches `data()`, strides, or layout. `L-11`
- The fix worked — now look for the second defect it was hiding. `L-10`

**Before editing a shared checkout with other live sessions**
- `git status` before EVERY edit, not at session start — a peer's half-finished refactor can be
  sitting in the exact file your one-line brief names (INC-556). `L-05`
- Stage explicit paths; a broad `git add` swept three sessions' WIP into other people's commits
  in one day. State in the message whose work rode along (INC-556).
- Prose that duplicates a row-derived fact (a count, a tally) rots silently under multiple
  writers — machine-check it or do not write it. The register's tag counts are now enforced by
  `check_doc_claims.py`, whose first run caught a real drift within minutes (INC-556). `L-05`
- **Never `git checkout -- <file>` in a shared tree. Use `git stash`.** `checkout --` is the one
  git verb with no undo; `stash` is reversible and costs nothing when you turn out to be right
  that the file was a no-op. A session closing out ran it on `devlog/README.md` believing the
  change was empty and destroyed another session's index line (INC-558). `L-05`
- Before concluding "modified but no content change", read `git diff --stat`, not a grep of
  `git diff`. The grep returned nothing and `--stat` said `1 insertion(+)` — the contradicting
  evidence was on screen and the weaker signal was the one acted on (INC-558). `L-01`

**Before writing a `file:line`, a count, or "the only" / "none"**
- Re-open the **named file**. A line number can be genuine, freshly read, and from another
  document entirely — a run log and a source file both have a line 450. `L-08`
- Universal or negative claim → print the hits and count the printout. Never report the
  number alone; both such claims made this way were wrong. `L-08`
- Citing a source tree as authority → match its version to the linked binary first.
  `/opt/tpls_old` is one release behind what BELFEM links. `L-08`
- Generated or bulk-renamed code → grep it for `file:line` before trusting it. A rename
  turns checked citations into invented ones and nothing compiles-checks a comment. `L-08`
- Committing a comment that states WHY a change is right → verify the mechanism against the
  primary source, not the reporter. A correct change with a false reason survives every review
  and poisons the next design that consults it; one was caught only because its author
  volunteered the retraction (INC-551). `L-08`

**Before writing a jury brief**
- The "established facts — do not re-litigate" section may contain **only firsthand,
  re-read facts.** One secondhand number there wasted a whole round. `L-03`, `L-08`
- Pre-register predictions and the decision rule before any data exists. `L-03`
- State what result would falsify the hypothesis. If none would, it is not an experiment. `L-01`

**Before closing a DR or a todo item**
- Read the tree, not the plan. `L-05`
- Striking the design does not retire the run gate — the status column is never struck. `L-05`
- A stale row is **not** reliably pessimistic; it has been found worse than recorded. `L-05`

**Before releasing**
- Every `file:line` in a shipped document is advisory; anchors are greppable tokens. `L-08`
- Input-key changes land in `doc/input_file_reference.md` **and** `doc/input_schema.yaml`
  in the same session — including when only the *behavior* changed. `L-05`

## Failure signatures — symptom to first discriminator

| symptom | look here first | cheapest discriminator |
|---|---|---|
| **A numerical result regressed** | The yardstick, before the physics. Every regression chased in this corpus was a measurement or contract artifact. | Establish what the printed number *is*: linear exit residual vs nonlinear residual, display clamp vs true value, dB vs linear. `L-04` |
| **Wrong only in parallel** | Ownership, not physics. | Which rank owns it; what a non-owner does with the map lookup. `L-11` |
| **Wrong only after a restart** | The first post-restore step. | Compare against a *cold* run at the same state — the artifact that fires once per process is not about the restore. `L-12` |
| **Only a restart clears it** — in-run retries wall, Δt-insensitive, but a warm restart passes | The soft reset restores less state than the load path. | Diff the reject path against the load path member by member; first check whether the residual is evaluated against un-rewound state (dofs vs fields). `L-19` |
| **Rejections cluster while the residuals were improving** | The guard, not the physics. | Was the killed attempt at/near its best-ever and within reach of tolerance? Count guard-attributed cuts separately from convergence failures before tuning anything else. `L-16`, `L-20` |
| **A fix changed nothing** | The path may never have run. | Work count, or break it deliberately and confirm the failure. `L-09` |
| **Convergence stalls at a suspiciously round floor** | A tolerance or exit test, not the solver. | Which tolerance is active, against which norm, with whose default. `L-04` |
| **"Converged" in very few iterations** | The residual may be an algebraic identity. | Recompute the residual from a fresh assembly; check what the iterate was multiplied against. `L-04` |
| **A guard under test does not fire** | The binary, before the guard. A stale object makes a guard ABSENT, and absent is indistinguishable from broken. Probes fail OPEN. | Diff object mtimes against source mtimes; pair the guard test with its inverse, so a stale binary cannot show green in both directions. `L-08`, `L-09` |
| **The suite is green but the thing is broken** | Coverage, not correctness. | Is the test in the `check` dependency list; is it behind an `#ifdef`; does its assert fire if you break the code. `L-01`, `L-09` |
| **A step reads as an off-by-one, a redundant re-run, or a trivially hoistable loop** — in code implementing a published or unpublished algorithm | Your prior, not the code. For a *modified* algorithm the nearest famous one is the wrong reference, and the mismatch presents as a defect. | Find the primary source's side conditions before touching it. No paper, or the paper defers the step to a reference you cannot read → it is not yours to change. `L-22` |
| **An assert fires on a healthy case** | The guard, not the data. | Read the guard against the legal state set: empty, worker rank, NaN, degenerate. `L-16` |
| **Same error message, second time** | Two unrelated causes, one message. | Diff the two failing inputs; do not assume recurrence. `L-10` |
| **A quantity is zero where physics says it should not be** | Something upstream never ran or never wrote. | Trace back to the producer; confirm it has a consumer *and* a producer. `L-09` |
| **Results changed after only a doc/comment edit** | A stale build artifact or an orphaned duplicate file. | Confirm the edited file is the one the target compiles. `L-02` |
| **Editing an input changed nothing** | A cache keyed on the wrong thing. | Check what the cache checksum covers versus what the file actually contains. `L-02` |
| **Two AIs agree and the run disagrees** | The run is right. | Stop reviewing; the next step is the gate, not another opinion. `L-03` |
| **A grep found nothing, so the feature is missing** | The premise is fine; the conclusion is not. | Look for the other channel before promoting an absence to a defect. `L-07` |
| **An audit "verified" something that later broke** | The word, not the auditor. | Re-read what gate actually ran. `L-03` |
| **A performance cost is "obviously" in X** | Attribution without measurement. | Profile before attributing. (annex, N5) |
| **A pooled resource exhausts at a suspiciously exact count** | The allocator, not the consumer. High-water counters never recycle freed slots while one instance stays alive. | Read the free path: does it rescan, or only reset when ALL are free? Count = held + per-step burns (INC-555). `L-16` |
| **A fix landed, runs, and changed nothing measurable** | Its own guard or fallback may be swallowing it — executing is not working. | Assert the new effect differs from old behavior at one point where they must differ (INC-549). `L-01` |

## Escalation triggers — stop and route to Christian

More analysis is the wrong move when:

1. **Measurement contradicts a mechanism all voices endorsed.** The measurement wins, and the
   next decision is his. This has resolved correctly every time in the corpus.
2. **The question is physics, formulation, or a modeling convention.** Never settled by model
   vote — a wrap element passed six weeks of audits before a basis-completeness objection
   killed it.
3. **Two viable fix designs remain and both are defensible.** He has killed two at once with
   one data-layout insight.
4. **A premise came from him and the evidence now contradicts it** — say so directly; several
   of the most valuable corrections in the corpus began this way.
5. **A "portability gap" or "missing support" in the build or MPI layer.** Some of these are
   deliberate, validated choices that look like omissions.
6. **The next step would destroy or overwrite a reproducer.** Copy first, ask second.

---

# LAYER 2 — LESSON CARDS

*Evidence and reasoning. 22 cards. Read one when Layer 1 sends you here.*

## L-01 State the falsifier before you run the check
**Domain:** PROC
**Trigger:** Designing any experiment, probe, control arm, A/B, or acceptance test; or reading a green result as evidence.
**Rule:** Before running a check, write down which outcome would falsify the hypothesis. If no outcome would, it is not a discriminator — redesign it or drop it. When a check comes back clean, confirm it *could* have failed before believing it.
**Failure signature:** A check passes and the defect persists → the check was structurally unable to report it (wrong tree, wrong deck, path not enabled, assert unreachable, test not in the dependency list) → break the code deliberately and confirm the check goes red.
**Missed signal:** In nearly every case the entry already contained the refutation — a probe counter reading zero, a deck without the key, a suite that never built the binary. The information was present and unread.
**Evidence:** 83 incidents, 73 entries, 28 threads, 10 pre-distilled. INC-020 (a "successful periodic run" whose deck was a non-periodic variant — zero diagnostic calls), INC-344 (two probes vacuous: resumed time past the horizon, zero steps assembled), INC-290 (a history-consistency probe trivially zero by construction, comparing a field against itself right after the shift defined it), INC-336 (a polynomial-exactness test at h=3e-4 where a wrong-order scheme is indistinguishable from roundoff), INC-154 (a whole test block behind an off-by-default `#ifdef`), INC-141 (`tests/core` registered with ctest but absent from the `check` dependency list), INC-384 (an alarming result that was a `sed` in the fixture never matching), INC-176 (a review round where the auditor was handed the charter instead of the artifact). Positive instances, kept because the rule is easier to trust when it has paid out: INC-549 (a probe asserting the new path's EFFECT caught a fix shipping 100 % inert behind its own guard), INC-557 (a temporary forced-reset probe drove a branch no healthy deck reaches, confirming a jury-flagged hang fix by execution).
**Enforcement:** `doc/ai_collaboration_protocol.md` — add a required "falsifier" line to any probe or experiment proposal, and an auditor duty to check experiment design, not only conclusions. `Status: proposed`
**Applicability:** Universal.
**Exceptions:** None. A check with no failure mode is not cheaper than no check; it is more expensive, because it produces false confidence.

## L-02 Validate stale-able state against what determines it
**Domain:** ENG
**Trigger:** Reading a cache, a build artifact, a persisted file, a default value, or any state the current process did not just produce.
**Rule:** Validate cached or persisted state against everything that determines its contents, not against its nominal input. When validation fails, say so loudly — never fall back silently.
**Failure signature:** Editing an input changes nothing, or a doc-only edit changes results → a cache keyed on the wrong thing, or an orphaned duplicate file the target does not compile → compare what the checksum covers against what actually feeds the artifact.
**Missed signal:** The `.bfm` case is the clearest: the file caches the *enriched* mesh — cuts, thin-shell layers, coating walls, periodicity, renumbering — and was reused on a checksum of the *raw* geometry alone. Everything that made the file worth caching was outside its validity check.
**Evidence:** 66 incidents, 57 entries, 22 threads, 8 pre-distilled. INC-083 (a regression that was two default-initialiser flips; an exhaustive call-path audit of ~1500 lines across 57 files missed it, and diffing default member values found it immediately), INC-118 (a `!= BELFEM_QUIET_NAN` guard always true because NaN compares unequal to everything), INC-254 (`BELFEM_ERROR`'s check compiled into every build but its *reaction* did not), INC-257 (a test pinning a default tolerance the source had moved past), INC-137 (a Doxygen layout path silently unresolved because the target ran in a different working directory), INC-193 (Doxygen merging two backend headers into one entity and picking one brief).
**Enforcement:** ENG — extend the `.bfm` config-tag check to every cached artifact; add a startup assert that a loaded cache's tag covers the enrichment inputs. PROC — add "diff default initializers and class-level constants" as an explicit numbered step in the regression-hunting checklist. `Status: proposed`
**Applicability:** Universal; sharpest for `.bfm`, build trees, and vendor defaults.
**Exceptions:** None.

## L-03 Name the rung
**Domain:** PROC
**Trigger:** Writing "verified", "confirmed", "clean", or "high confidence"; assembling a jury brief; deciding whether to run another review round.
**Rule:** State which rung of the evidence ladder supports the claim. "Verified" requires a named executable gate that ran. Concurring audits never lift a claim past a static source trace. A jury brief's "established facts" section may contain only firsthand, re-read facts.
**Failure signature:** A reviewed-clean change fails on first run → the claim rested on agreement or on a build → stop reviewing and run the gate; the stop condition is already in the protocol.
**Missed signal:** The word itself. In INC-090 three-AI agreement on a sign structure was recorded as verification while five errors were caught in that same round; in INC-262 three read-only passes signed off on code sitting below an early `return`, and a smoke run caught it immediately.
**Evidence:** 62 incidents, 62 entries, 27 threads, 12 pre-distilled — the widest thread spread in the catalog, so this fails everywhere rather than in one corner. INC-218 (a jury round wasted because its brief's do-not-re-litigate section carried a secondhand number that was two kernels' counts added together; retracted the same day when Christian did the arithmetic), INC-002 ("`make hphirun` passed" used as the verification line while eight semantic defects sat in the merge), INC-198 (an incompatibility flagged in devlog, plan and commit message, with flagging treated as handling — the code did not compile on that backend), INC-360 (a mechanism "verified by execution", retracted, then re-confirmed — the retraction rested on a stale viewer frame), INC-255 (a pre-registered default refuted by an auditor citing 52 dependent sites).
**Enforcement:** `doc/ai_collaboration_protocol.md` §11 — add the firsthand-only rule for brief "established facts", and require every "verified" to name its gate inline. `Status: proposed`
**Applicability:** Universal.
**Exceptions:** None. The vocabulary rule is cheap and it did more work than any other single practice in this corpus.

## L-04 Audit the yardstick before the physics
**Domain:** ENG
**Trigger:** Any numerical regression, stall, unexpected convergence, or "the number looks wrong".
**Rule:** Before touching a formulation, establish exactly what the printed number measures — which residual, against which norm, with whose default, through which display path.
**Failure signature:** Convergence stalls at a round floor, or converges implausibly fast → an exit test, a display clamp, or a library default, not the solver → recompute the quantity independently from a fresh assembly and compare.
**Missed signal:** This is the corpus's single most transferable lesson, stated in the source feedback document: *every* regression chased turned out to be a measurement or contract artifact, not wrong physics.
**Evidence:** 51 incidents, 55 entries, 16 threads, 6 pre-distilled. INC-231 (the printed nonlinear residual **is** the linear solve's exit residual; at rtol 1e-8 GMRES stopped on a ~4e-9 plateau, Picard "stalled", the controller promoted to Newton on exit-test noise — Christian's A/B refuted Claude's pre-registered floor model on all three predictions), INC-230 (the floor attributed to factorization accuracy when occasional 1e-16 finals proved the factor was better than the floor), INC-275 (`9.000000` was a display clamp while the dB column carried ~5.8e4), INC-226 (a relative refinement criterion against a numerically zero RHS, grinding hardest on the *best* iterates), INC-321 (a 5× speedup whose obvious causal reading was wrong — the thermal solve was 520-672 ms against a 24-33 s magnetic factorization and was never the bottleneck), INC-467 (a vendor workspace query returning below the reference minimum), INC-287 (defaults flipped as "safely reversible" that silently moved keyless decks to a different scheme).
**Enforcement:** ENG — one-line provenance in the residual print path naming which residual it is. PROC — "what does this number measure" as the mandatory first question in any regression round; a run gate, not a review gate, for every default change. `Status: proposed`
**Applicability:** Universal; strongest for solver tolerances, residual reporting, and vendor defaults.
**Exceptions:** None.

## L-05 The evidence is the tree
**Domain:** PROC
**Trigger:** Reading a plan, a register row, a devlog, an index, or a status line; closing a DR; running a currency sweep.
**Rule:** Treat every planning artifact as a claim to check against the tree, never as a source. Closing a row requires reading the code. Striking a design does not retire its run gate.
**Failure signature:** Work is planned around a documented state that no longer exists → the artifact outlived the behavior → grep the tree for the named symbol before believing the row.
**Missed signal:** Staleness is **not reliably pessimistic**. The corpus recorded a four-for-four run of stale rows being conservative, then broke it twice: two rows were found *worse* than recorded.
**Evidence:** 47 incidents from only 30 entries across 20 threads, 9 pre-distilled — the count is inflated by sweep entries, which is why the entries column exists. INC-144 (five plans finished without being closed, nine status lines actively misdescribing the tree), INC-145 and INC-506 (rows written as open that were already fixed — one of them fixed in the same devlog's own addendum), INC-496 ("struck is not verified" — the run gate survives in the status column, which is why that column is never struck), INC-500 (a diagnostic recorded as landed whose functions exist with zero consumers tree-wide), INC-133 (the devlog index missing entries on disk), INC-504 (an input key whose *behavior* changed while both contract artifacts still described the reverted behavior).
**Enforcement:** Register protocol — closure requires two new fields, "why prior tests missed it" and "residual risks", plus a lessons check naming the violated card. `scripts/check_doc_claims.py` extended in the code→document direction. `Status: proposed`
**Applicability:** Universal.
**Exceptions:** None.

## L-06 One source of truth per convention, tested at every consumer
**Domain:** ENG
**Trigger:** Touching orientation, node ordering, master/slave roles, facet numbering, argument order, or an index space; adding a second consumer of an existing convention.
**Rule:** A convention has exactly one authority and a test at every consumer. When you change one, enumerate all of them — a stale consumer will still compile and still look plausible.
**Failure signature:** Results plausible but wrong in one configuration only → a consumer reading the convention differently → compare the two consumers side by side against the authority, not against each other.
**Missed signal:** The TET10 case is the sharpest: the wrong node map lived in one line of a generator while its *sibling cross-check script carried the correct map with a "careful, these are swapped" comment*. The convention was being consciously fought in one place and silently violated in another, with no gate between the generator and the hand-written C++.
**Evidence:** 46 incidents, 43 entries, 15 threads, 5 pre-distilled. INC-115 (the T-first `(T,B,beta)` reorder applied inconsistently; because every argument is `real`, stale sites compiled silently and evaluated resistivity at T = |B|), INC-245 (the TET10 generator trap), INC-243 (`EF_TET4::E()` computing `η∇ξ − ξ∇ζ` where the Whitney form is `η∇ξ − ξ∇η`, contradicted by the file's own comment and its separately-coded curl), INC-237 (a mesh-level edge-direction flip fixing HEX8 while HEX20/27/64 kept the old convention), INC-268 (HEX27 derivative rows swapped, breaking ΣdNᵢ = 0), INC-284 (two consumers indexing layer nodes with and without `original()` normalization).
**Enforcement:** ENG — extend the interface/orientation battery in `make check-fast` to every element family that shares a convention; a generator that emits tables emits a checksum the consumer asserts. `Status: proposed`
**Applicability:** Elements, meshes, material APIs, index spaces.
**Exceptions:** None.

## L-07 Carry the scope of every absence
**Domain:** AGENT
**Trigger:** Writing "no", "never", "every", "all", "nothing", "missing", "unwired", "not present"; concluding a defect from a grep.
**Rule:** A universal or negative claim carries its enumeration or its search scope in the same sentence. Before promoting an absence to a defect, look for the other channel.
**Failure signature:** A confident gap report that the person running the code contradicts within the hour → the grep premise was true and the conclusion did not follow → ask what else could carry the thing before concluding it is absent.
**Missed signal:** In INC-170/INC-502 the grep was correct — the connector record genuinely does not appear in the shell record — but it is *deliberately* absent because the data ships by two other explicit channels. The premise was sound; only the inference was wrong.
**Evidence:** 42 incidents, 46 entries, 17 threads, 5 pre-distilled. INC-197 ("nothing includes `fn_eigen.hpp`, therefore it has no build coverage" — written twice in one session, refuted later that session by a test that includes it and five existing tests), INC-199 ("zero `@brief` across all 44 files" — the count excluded the best-documented directory in the module), INC-001 (a branch audited as carrying machinery that lived only on an older unmerged branch, written at "medium-high" confidence), INC-105 (a case statement believed functional whose enum value the factory omitted at both selection sites).
**Enforcement:** `CLAUDE.md` standing instruction — a negative claim states its pattern and its searched paths inline; an unscoped absence is not a finding. `Status: proposed`
**Applicability:** Universal, and specifically to grep-derived conclusions.
**Exceptions:** None. The scope clause costs a few words.

## L-08 Re-open every literal
**Domain:** AGENT
**Trigger:** Writing any number, ID, commit hash, `file:line`, tracker row, attribution, or measured value — especially late in a long session.
**Rule:** Open the source and read the literal in the session that cites it. Recognition is not lookup. A plausible ID reads correct at a glance and that is exactly the danger.
**Rule, extended 2026-08-29 — three shapes that pass the rule as originally written:**
- **Re-open the NAMED FILE, not just a file.** A line number can be genuine, freshly read, and from the wrong file. `cl_SolverPETSC.cpp:450` was cited from a run log's line 450 read minutes earlier (INC-543). A lookup happened; it was of the wrong document.
- **A source tree is not evidence about a binary until its version is matched to the linked library.** Sources under a plausible path can be a different release from the one that runs (INC-546).
- **A mechanical rename converts verified citations into invented ones**, and nothing compiles-checks a comment. Grep generated text for `file:line` patterns before trusting it (INC-547).
**Failure signature:** A claim survives review rounds and then dissolves when someone follows the pointer → the pointer was confirmed by recognition → re-follow it before defending the claim.
**Second signature (2026-08-29):** a COUNT or an enumeration stated as a fact about the tree, produced by a search whose scope was narrower than the claim. Both instances were caught by others, and one was nearly reported back as a discrepancy in a peer's correct work (INC-544, INC-545). Discriminator: if the claim is universal ("the only", "none", "N rows"), print the hits and count the printout — never report the number alone.
**Missed signal:** The measurement that forced the existing anchoring rule: of 87 `file:line` citations in one document, 26 were stale and 24 had shifted by the same offset — and the document and the code had been edited *in the same commits*. Discipline was being followed and still failed, because a line number is already wrong by the time it is committed.
**Evidence:** 37 incidents, 25 entries, 15 threads, 6 pre-distilled. **INC-543 to INC-547 (2026-08-29) are five violations in ONE session, by an agent that had this card available** — the rule did not fail, it failed to FIRE; that is the datum, and it is why the enforcement line below now names a mechanical check rather than a standing instruction. Earlier: INC-166 and INC-501 (an incident attributed to the wrong row, wrong on both halves — the revert is `61bbfa11` dated 2026-08-11, and the cited row is unrelated; it survived a jury round plus two independent reads), INC-216 (three references confirmed by recognition, correct at one commit and stale at the next), INC-390 (a map asserted to exist from session memory and refuted — it had been removed), INC-531 (an entire exchange conducted on reported numbers, where the decisive bit-identical value was the *previous run's file still open in the viewer*), INC-405 (text from the shell/`rg` channel mangling words, so quoted identifiers had to be Read-confirmed).
**Enforcement:** `CLAUDE.md` standing instruction — any specific number or ID in a final document carries a same-session re-read; without one it is UNVERIFIABLE by definition. Machine contracts anchor on greppable tokens, never line numbers. **Partially mechanized 2026-08-29:** `scripts/check_doc_claims.py` now recounts the debt register's open `[P]`/`[W]` rows and diffs them against the header sentence — the one count in this corpus that had drifted three times in a day. It caught the author of the rule getting it wrong within the hour of writing it. The general case (arbitrary counts, cross-file line anchors) remains unmechanised. `Status: partially enforced`
**Applicability:** Universal; risk rises with session length.
**Exceptions:** None.

## L-09 Presence is not execution
**Domain:** ENG
**Trigger:** Reasoning about what a code path does; auditing new wiring; claiming a feature works; reading a plan that says a mechanism exists.
**Rule:** Before reasoning about a path, prove it executes — a work count, a probe, or a deliberate break that must fail. A path with no consumer, no test and no run is unknown code, not working code.
**Failure signature:** A fix changes nothing, or a quantity is zero where physics says it should not be → the path never ran → count the work it claims to do (sources set, edges tied, facets tagged) and check the count is non-zero.
**Missed signal:** Read-only audits are structurally bad at this, and the corpus says so repeatedly: three read-only passes signed off on code below an early `return` (INC-262), and the first smoke run caught it in seconds. Reading tells you what code *would* do; only running tells you whether it *does*.
**Evidence:** 24 incidents, 29 entries, 12 threads. This cluster **was not in the pre-validated seed list** and was proposed independently by three of six blind assigners working on disjoint slices. INC-005 (a flagging call placed before the thing it flags is created), INC-025 (loops gated on flags nothing sets until after the loop), INC-041 (weights computed and never handed to `set_sources`), INC-073 (a helper declared and defined, called by nothing), INC-064 (a call behind an always-aborting `default:`), INC-403 (a Newton path forced to Picard every timestep), INC-438 (a branch hardcoded off and broken when forced), INC-141 (`tests/core` never built by `make check`), INC-523 (a committed test calling an API that never existed), INC-458 (a `getri` wrapper whose final call re-factorized instead of inverting), INC-150 (a parallel database path that silently produced nothing).
**Enforcement:** ENG — new wiring lands with an assert on its work count, not only on its inputs; every test target appears in the `check` dependency list and CI-less discipline is replaced by a `check` completeness assertion. PROC — "has this executed?" as a standing audit question with a named answer. `Status: proposed`
**Applicability:** Universal; highest yield on newly wired paths and vendor wrappers.
**Exceptions:** None.

## L-10 A root cause is not the sole cause
**Domain:** PROC
**Trigger:** A fix works; a symptom recurs after a fix; the same error message appears twice.
**Rule:** After a confirmed fix, look for the second defect it was hiding. Treat a recurring message as two causes until the inputs are diffed.
**Failure signature:** The identical message returns after a successful fix → two unrelated causes sharing one message → diff the two failing inputs rather than assuming recurrence.
**Missed signal:** Masking is usually *symmetric* — two sign errors canceling at the surface (INC-367), a crossed table row that coincides with the correct value at evaluation point 0 (INC-246). Both were invisible to any run that only looked at the place they cancelled.
**Evidence:** 19 incidents, 22 entries, 15 threads. INC-317 (two unrelated causes behind one identical message, one run apart), INC-058 (a second fuse flag missed by the first switch-off, found by an end-to-end A/B), INC-367 (two independent sign errors, one masking the other at the outer surface), INC-035 (a symptom attributed to a known defect whose gates were already in tree and functional), INC-233 (one counter serving three incompatible roles, its latent bugs surfacing only once split).
**Enforcement:** Register protocol — closure requires a "why prior tests missed it" field, which forces the masking question. `Status: proposed`
**Applicability:** Universal.
**Exceptions:** None.

## L-11 Serial correctness says nothing about ownership
**Domain:** ENG
**Trigger:** Shipping a change validated in serial; touching maps, gathers, duplicates, postprocessing, or anything with an owner; changing matrix transfers or backend-dependent layout.
**Rule:** Before shipping, state which rank owns each entity and what a non-owner does. Run the second backend if the change touches `data()`, strides, or layout.
**Failure signature:** Correct in serial, wrong or hanging in parallel → ownership advertised but not used, or a hidden collective → check for a rank guard around the producer and what the consumer does without it.
**Missed signal:** "Serial works" was repeatedly read as "the file and the load path are correct". In INC-391 the serial load was *silently wrong too* — 16761 facets instead of 12285 — and the parallel crash was the first thing loud enough to notice.
**Evidence:** 17 incidents, 16 entries, 10 threads. INC-046 (a map populated on root and dereferenced on every rank), INC-070 (a bug parallel-only: all duplicates of an owned original inserted into the owned-node bitset without an ownership check), INC-143 (Matrix MPI transfers shipping a stale `capacity()` because Blaze never shrinks it), INC-299 (a duplicate whose own-side elements live on a non-owner rank claimed by nobody, leaving its row silently 0.0), INC-534 (a fix distributing all history levels while verifying only the read ones — a near-miss caught by audit), INC-517 (an error rejecting the legitimately empty worker-rank table).
**Enforcement:** ENG — a parallel smoke run on 2 and 4 ranks as a merge gate for anything touching ownership, transfers, or persistence; extend the existing battery to run under both backends. `Status: proposed`
**Applicability:** MPI paths, persistence, postprocessing, matrix transfers.
**Exceptions:** Genuinely serial-only utilities, stated explicitly.

## L-12 Restart is not continuation
**Domain:** ENG
**Trigger:** Any save/load, memdump, warm start, or resumed run; changing what a dump contains.
**Rule:** Treat the first post-restore step as a different program. Seeding, first Δt, integrator history and solver birth state are each explicitly restored or explicitly re-derived — and the choice is written down.
**Failure signature:** Healthy live run, sick resumed run at the same state → something the process earned over time was not restored → compare against a **cold** run at the same state; an artifact that fires once per process is not about the restore.
**Missed signal:** The contract was documented in four places as "BDF spacing deliberately not persisted" while the dump carried hundreds of megabytes of field history — the expensive half was saved and the three integrator members that made it usable were not.
**Evidence:** 11 incidents, 11 entries, 4 threads. INC-339 (a restart discarding its full history and re-anchoring at BDF1), INC-441 (`load_memdump` called before `set_circuit`, so the load was a silent no-op), INC-393 (a topology map built pre-enrichment and reused post-load), INC-012 (a restore path unusable on one construction route because Hesse forms were absent). The thread's final state caps the re-entry step at the deck's `initial timestep` and deliberately reverses the Δt half of the earlier restart-cliff policy — with the mechanism still only a leading hypothesis, which is itself worth knowing before someone "fixes" it again.
**Enforcement:** ENG — a restart regression: dump at step N, restore, and require the first step's iterate count to match the live run's within a stated tolerance. `Status: proposed`
**Applicability:** All persistence and restart paths.
**Exceptions:** None.

## L-13 Know which forks are not yours
**Domain:** AGENT
**Trigger:** Physics or formulation questions; two viable designs; measurement contradicting a unanimous mechanism; anything that looks like a portability gap in build or MPI.
**Rule:** Stop and route to Christian rather than analyzing further. More analysis is the wrong move when the missing input is a domain prior, not a fact about the code.
**Failure signature:** Rounds of audit converge on a mechanism and the run disagrees → the question was never a code question → present the measurement and the two candidate readings, and ask.
**Missed signal:** The corpus's highest-value corrections all arrived this way and all arrived *fast*: a tolerance hypothesis Claude had argued against yielded a large throughput gain when Christian tested it; a basis-completeness objection killed a design six weeks of source audits had not; "OpenMPI is the only validated MPI" turned a written-up portability gap into a deliberate choice.
**Evidence:** 9 incidents, 10 entries, 5 threads, 2 pre-distilled. INC-164 (two hardcoded OpenMPI sites written up as a gap; the premise corrected — PETSc has always crashed on MPICH), INC-530 (a rho-database index-space "trap" ruled to be design, not defect), INC-037 (a surplus cut dof believed spurious, physically legitimate and required to stay free), INC-039 (the wrap unphysical at the basis level), INC-269 (an attribution corrected by Christian the same day).
**Enforcement:** `CLAUDE.md` standing instruction listing the six escalation triggers from Layer 1. `Status: proposed`
**Applicability:** Universal.
**Exceptions:** None — and the cost of a wrong escalation is minutes, against days for a wrong analysis.

## L-14 Widen the sample before the model sets a threshold
**Domain:** PROC
**Trigger:** Building a model, a threshold, or a severity from one deck, one metal, one mesh, one hardware, one run.
**Rule:** A model built from one sample is a hypothesis about that sample. Widen the sample before it sets a threshold, a severity, or a decision.
**Failure signature:** A clean first port followed by a grossly wrong second → the first sample was benign, not representative → run the widest cheap sweep available before generalizing.
**Missed signal:** The benign-first-sample trap is the sharpest form. Copper's Bézier spans happen to be balanced, so the first audit pass came back clean; Silver and Indium then exposed a root-finder converging on an extrapolated second root, and a census found five of seventeen in-tree bases non-monotone on the wider bracket.
**Evidence:** 8 incidents, 7 entries, 7 threads, 2 pre-distilled. INC-308 (a mesh-format claim written into the register, two index files, a plan and a README, refuted by `head -3` on a working counterexample that had been in the tree all along), INC-379 (a crossover measured as non-existent on one machine, contradicted by an independent re-measurement on different hardware), INC-418 (a "six orders of magnitude apart" threshold measured on too small a sample), INC-132 (a flagged fit defect withdrawn once the comparison widened across metals), INC-532 (a linear drift-accumulation model that holds only without conduction or source).
**Enforcement:** PROC — a threshold or severity derived from measurement states its sample size and span inline. `Status: proposed`
**Applicability:** Fits, thresholds, performance crossovers, severity ratings.
**Exceptions:** None.

## L-15 Protect the evidence and the machine before the experiment
**Domain:** AGENT
**Trigger:** Any restart recipe, cleanup, overwrite, concurrent job, log watcher, or subprocess you are tempted to kill.
**Rule:** Copy endangered artifacts yourself, in the same message that proposes the destructive step — never delegated. Budget memory before launching concurrent jobs. Anchor a log watcher to a unique event marker. Low CPU is not a hung process.
**Failure signature:** A reproducer is gone, or a job died, or a test "failed" and did not → an operational shortcut, not a code defect → re-read the raw artifact by hand before acting on any watcher's verdict.
**Missed signal:** In INC-350 the recipe was correct and its first line *was* the backup; it was delegated, skipped under a mistaken belief about the state's age, and the only reproducer of a live anomaly is permanently gone. The standing rule that came out of it is Iron Rule 7.
**Evidence:** 8 incidents, 8 entries, 6 threads, 4 pre-distilled. INC-350 (the lost t=3100 dump), INC-347 (a false test failure from a watcher whose baseline reset when one poll returned 0), INC-377 (`mktemp` templates that on BSD created literal shared files, so concurrent runs collided), INC-436 (a validated 18-check harness that lived in a scratchpad and was swept, shipping the production code with no test at all), INC-345 (hand-launched runs missing the thread environment the scripted ones set), INC-407 (an auditor truncated by a turn budget, succeeding on re-run with a raised one).
**Enforcement:** `CLAUDE.md` standing instructions — the backup rule, the subprocess patience norm (never kill on CPU% or elapsed time below a stated floor; check output growth and RSS), and the unique-marker rule for watchers. `Status: proposed`
**Applicability:** Universal, and acute whenever a long-running job holds the machine.
**Exceptions:** None.

## L-16 A new guard is a new failure mode
**Domain:** ENG
**Trigger:** Adding or tightening a `BELFEM_ASSERT` or `BELFEM_ERROR`, especially one written from reading a single call site.
**Rule:** Before adding a guard, enumerate the legal states it must accept — empty containers, worker ranks, NaN, degenerate geometry, boundary parameter values. Then check the comparison actually has the sense you intend.
**Failure signature:** An abort on a case that used to work → the guard, not the data → read the guard against the legal state set, and check NaN and empty separately since both defeat naive comparisons.
**Missed signal:** These guards were added *as hardening*, which is why they were not scrutinized. INC-155's own verdict states the missing rule verbatim: it is "not hardening; it is a new failure mode".
**Evidence:** 4 incidents, 5 entries, 3 threads, 1 pre-distilled. INC-096 (an inverted assert firing on every healthy interface in debug while leaving the release divide-by-zero path completely unguarded — wrong in both directions at once), INC-122 (an assert on a discriminant that legitimately vanishes at n = 4), INC-155 (an error rejecting the deliberately rank-tolerant empty worker table), INC-126 (a sharing condition checking `Py(3) == Qy(3)` where the curves share the last P and the *first* Q).
**Enforcement:** ENG — a guard added to a hot or shared path ships with a test that exercises its legal-but-extreme states. `Status: proposed`
**Applicability:** All assertion and error sites; highest risk in shared framework code.
**Exceptions:** None.

## L-17 Name what breaks before deleting the redundant
**Domain:** ENG
**Trigger:** Removing anything as dead, redundant, double-counted, or pure overhead — a term, a loop, an include, a branch.
**Rule:** Before deleting, name the specific thing that would break if it were load-bearing, and check that thing. "Dead" and "inert" are different properties.
**Failure signature:** A deletion compiles clean and results move → the construct had a side effect its name did not advertise → restore it and bisect the observable it was holding up.
**Missed signal:** All three cases were *reasoned* deletions with plausible arguments, and the reasoning was the problem. In INC-426 the commit message even asserted a supporting number — reasoned from the code, never measured — and every cubic gas silently returned an absolute entropy one standard state too low until `make check` caught it at r² = −56.4.
**Evidence:** 3 incidents, 3 entries, 3 threads, 1 pre-distilled. INC-426 (`mSref` removed as a "double count" while it was compensating the departure spline's standard-state anchor), INC-396 (a loop deleted as "pure overhead" that was resetting a key count from the raw facet count to the post-`unique()` count, leaving connectivity ~2× oversized), INC-444 (includes that were dead but not inert, transitively supplying a `using namespace` and a constant).
**Enforcement:** ENG — a deletion-only change is a run gate, not a review gate, when it touches a numerical path. `Status: proposed`
**Applicability:** Universal; highest risk for terms in thermodynamic and constitutive expressions.
**Exceptions:** Deletions covered by a test that would fail if the construct mattered.

## L-18 Match the failure tier to the recovery above it
**Domain:** ENG
**Trigger:** Deciding how a failure reports: assert, error, status return, or silence; adding a scope guard to a vendor return code.
**Rule:** A retryable state returns a status to the controller above it. A wrong-answer state hard-fails in release. Neither is ever silent, and neither aborts a run the controller could have rescued.
**Failure signature:** A run dies on a condition the controller handles, or continues on one it cannot → the tier is wrong → ask what policy exists above this call and whether it can act on the signal.
**Missed signal:** Vendor return codes rarely map one-to-one onto tiers. MUMPS returns warnings as a *sum* of flags, so a `> 0` test conflates an expected `+8` on ill-conditioned HTS systems with a genuine fault — both auditors raised the scope guard independently.
**Evidence:** 2 incidents, 2 entries, 1 thread — the thinnest card here, promoted deliberately because **the rule already exists in `CLAUDE.md`** ("Third category — expected algorithmic failure, never an abort") and these are that rule being violated anyway. A documented rule that is still violated needs enforcement, not restatement. INC-215 (the MUMPS warning-sum scope guard), plus the PETSc soft-fail work that turned a `DIVERGED_ITS` abort into a controller timestep cut.
**Enforcement:** ENG — audit every vendor return-code check for flag-sum semantics; `doc/coding_philosophy.md` already states the tier rule, so the gap is a check, not text. `Status: proposed`
**Applicability:** Solver wrappers, vendor interfaces, controller-adjacent code.
**Exceptions:** None.

---


## L-19 A restart that fixes it is a diagnosis, not a workaround

**Trigger:** Any wall that in-run retries cannot cross but a warm restart can; any manual
restart ladder; any state-rollback code (reset, revert, un-shift) reviewed or written.

**Failure signature:** Rejection cascades non-responsive to Δt, cleared by restarting from a
dump → the reject path restores a SUBSET of what the load path restores → diff them member by
member; the residual being evaluated against un-rewound state is the first suspect (fields
restored, dof values not — `reset_timestep` vs `load_memdump`'s `seed_dof_values`). A
non-finite iterate that survives the rollback repeats on every retry. Corollary: every
rollback needs a completeness argument against its corresponding cold path, and "same dump,
different outcome" run pairs measure FP nondeterminism, not mechanism — count consecutive
in-run failures against restart successes instead.

**Evidence:** INC-540 (dl20260828_quenchfront_night §1: the seed hole — three runs consumed,
a manual restart ladder, 15-reject cascades at 2 µs; one-line fix, validated in production the
same hour). The jury's differential audit also refuted the dispatching session's own prime
suspect — the BDF un-shift slot — because the retry's re-shift heals it before assembly.

## L-20 A stack of half-working mitigations marks one undiagnosed defect

**Trigger:** About to add another tuning knob, guard, or scheme change for a recurring
instability that already carries several; any knob that "helps a little" for reasons nobody
can state.

**Failure signature:** Multiple mitigations each move a wall without removing it, and their
benefits resist attribution → suspect a common root defect beneath the stack; find it before
adding layer N+1 → after the root fix, re-evaluate every layer on clean ground — the simple
configuration may now win. The compensations themselves become load-bearing folklore ("BDF5
handles the front better") that inverts once the defect is gone.

**Evidence:** INC-541 (dl20260828_quenchfront_night §4: a watchdog window tuned tight enough
to execute attempts 3 dB from convergence — calibrated against the pre-fix pathology), INC-542
(§5: BDF5+Newton ratcheted to 0.1 ms where plain BDF1 held 0.7–1.9 ms through the same front
once the seed hole was fixed — the tree's own guidance and the paper's original scheme
vindicated over the accumulated apparatus).
## L-21 Third-party libraries are accessed only through their dedicated wrappers

**Domain:** ENG
**Trigger:** About to write a call whose name carries a vendor prefix — `MPI_`, LAPACK's
`dgesvd_`, `mumps_`, `KSP`/`Mat`/`Vec`, `SP_`, HDF5's `H5*` — anywhere outside the module that
wraps that library. Also: reviewing a diff that adds one.

**Rule:** **HARD POLICY (Christian, 2026-08-30): third-party libraries are NEVER called
directly. Every access goes through the dedicated wrapper layer** — the comm module for MPI
(`commtools.hpp` for the collectives; `cl_Communicator.hpp` declares `comm_abort`, which core
must reach without pulling linalg), the `fn_*` LAPACK interface, the `Solver*` classes for the
sparse suites, the io module for HDF5 — in the fashion the codebase already applies
meticulously elsewhere. A missing wrapper is extended, not bypassed: add the function to the
wrapper layer first, then call it.

**Why the rule earns its keep here, not just aesthetics:** the wrappers are where this tree's
hard-won contracts live. `comm_type<T>` resolves the 32/64-bit `int_t` suite split that a raw
`MPI_INT` silently gets wrong; the abort path's `MPI_Initialized`/`MPI_Finalized` guards and its
deliberate use of `MPI_COMM_WORLD` over `gComm.world()` (contract at
`src/comm/cl_Communicator.hpp:235-258`, definition at `src/comm/commtools.cpp:67`) encode a
bug fixed on 2026-08-30 that any fresh direct call can re-introduce; chunking, tag discipline,
and empty-`data()` semantics live in `commtools.hpp` and nowhere else. A direct vendor call
starts from zero of that knowledge, every time.

**Failure signature:** A vendor-prefixed call sits in a file whose module does not own that
library → the call is re-deriving (or missing) a contract the wrapper already carries → move it
behind the wrapper; if the wrapper lacks the operation, that is the work item.

**Evidence:** INC-565, INC-566. The 2026-08-30 sweep found `MPI_Allreduce` open-coded in
`cl_SolverPETSC.cpp`, `cl_SolverSTRUMPACK.cpp` and `test_commmpi_main.cpp`, and `MPI_Abort`/
`MPI_Initialized`/`MPI_Finalized` in `assert.cpp` — four sites, three of them written that same
month, one of them (the test main) written the same day the policy was stated, by a session that
knew the tree well. The violations were not ignorance of the codebase; they were the default
behavior of writing "the normal way" under time pressure. That is what makes this a tripwire
rather than a style note. All four sites were closed the same day, through `allreduce` /
`comm_abort` in the comm module — after a jury round on the wrappers themselves found the first
drafts uncompilable and, in the abort case, regressing the very contract being wrapped.

**Enforcement:** ENG — mechanically checked by `scripts/check_wrapper_policy.py`, which matches
**calls** (a vendor-prefixed identifier followed by an open parenthesis, after comments and
string literals are stripped) for MPI, HDF5, PETSc, MUMPS, STRUMPACK, SuperLU, LAPACK/BLAS and
ARPACK, and prints every exclusion it applies by name. A plain prefix grep cannot do this job:
the recipe this card shipped until 2026-08-31 returned fourteen `MPI_` hits outside `src/comm/`,
**none** of them a call — comments, a serial-build `typedef int MPI_Comm`, parameters and
members, a Blaze macro, a PETSc error enumerator. Vendor *constants* crossing a module boundary
are out of the sweep's reach (not syntactically distinguishable from a mention) and remain review
work. `Status: hard policy, mechanically checked, green on the tree as of 2026-08-31`

**Applicability:** The C++ tree. The exception is the wrapper module itself, whose entire job is
the direct call.

**Exceptions:** Three, all ruled by Christian on 2026-08-31 and all encoded in the sweep as
**named, printed** entries rather than silent allowlist lines — an exception nobody can see is
indistinguishable from a hole.

- **The Fortran drivers** (`src/sparse/*.f90`). For MUMPS, PARDISO and ARPACK those drivers
  **are** the wrapper the C++ tree uses; their collectives stay inside the vendor package's own
  communication, and their contracts are the vendor's (`MPI_INTEGER`, not `comm_type<T>`). The
  wrapper requirement is specific to C++. This ruling makes a testable prediction, and it holds:
  MUMPS and ARPACK have **zero** direct C++ calls anywhere in the tree.
- **PETSc reaches into `src/comm/`.** PETSc's requirements are unusual enough that its scope
  extends past the sparse module: `PetscInitialize` is paired with `MPI_Init` in
  `Communicator::init` and cannot sit behind a `Solver*` class. `comm/` is a PETSc owner by
  ruling, not by module.
- **The HDF5 database probes in `src/physics/materials/`.** `cl_JcFunction_Database.hpp` and
  `fn_rho_database_is_current.hpp` open a file to test for a group that is a *sibling* of the one
  the io wrapper opened, and the wrapper exposes no root handle to probe with. Both say so in
  place. Scoped to HDF5 in that directory and nothing wider: an `MPI_` call in
  `physics/materials/` is still a finding, and an `H5*` call in `physics/gasmodels/` is still a
  finding.

Nothing else. "It is just one call" is the failure mode, not an exception, and a wrapper that
lacks the operation you need is the work item.

## L-22 An algorithm with no published implementation inverts your prior
**Domain:** PROC
**Trigger:** Editing, reviewing, or reasoning about code that implements a *modified*, in-house, or unpublished algorithm — anything whose nearest neighbor in the literature is not what the code does. In this tree: `src/homology/`'s cohomology core, and by the same test any future kernel whose paper is a draft.
**Rule:** Familiarity with the nearest published algorithm is a liability, not an asset, once the code deliberately departs from it. Consult the primary source and its side conditions, or derive from first principles, before concluding anything is a defect. Where neither is possible, report and do not repair. **Do not treat concurrence as coverage:** the missing prior is shared across vendors, so agreement between reviewers carries no more weight than one voice.
**Failure signature:** A step in unfamiliar code reads as an off-by-one, a redundant re-run, or a trivially hoistable loop → the prior is completing toward the textbook version → open the primary source and establish what the deviation *is* before calling it either deliberate or defective. Both answers occur (INC-560 vs INC-561); guessing which is the failure.
**Missed signal:** The deviation is *unknown until checked*, not presumed deliberate — INC-561 is the case where the in-tree deviation carried no correctness weight at all, and only the paper's own theorem established that. What is constant is that the reviewer could not tell which it was from the artifact. The strongest instance is the one where the vote went the wrong way: two of three voices agreed a correct operator was defective, and only an independent derivation broke it. The corpus contains no case where re-reading the code harder resolved one of these; it contains three where the primary source or a derivation did.
**Evidence:** 6 incidents. INC-559 (Claude + Grok concurred that a correct Smith-stage matrix was missing a merge projection; Codex refuted with a Schur derivation — a majority vote would have corrupted a working kernel), INC-560 (batch-omitting the 0-cochain removals in `coreduceOmit`, struck before it landed — and note that the "Mrozek §5 licenses one omission per component" gloss carried in the 2026-07-03 devlog and copied into this card's first draft is itself a near-miss error: the paper adds an ∅ cell of dimension −1, one per connected component, and licenses no 0-cell omission at all), INC-561 (the in-tree full rescan flagged as a deviation from the textbook worklist; the paper reversed the direction of the finding and Thm 6.2 dissolved its severity), INC-562 (a prior session's "we already have a worklist" reused as a template and refuted as a category error 8 weeks later), INC-563 (the policy document written to describe this failure mode exhibited four instances of it, all found by a cross-vendor jury). Apparent counter-instance, corrected: INC-564 (two AI edits *inside* the same code did land and hold — but one preserved a latent defect it did not understand, and the other was heavily human-guided surgery in which the human supplied the constraint that BOTH the SPFA and greedy algorithms are needed; only the dispatch wrapper was cut, and an AI framing it alone would have deleted the live rectifier, which runs at `cl_Cohomology.cpp:548` inside `clean_spfa`).
**Enforcement:** For the cohomology core, ENFORCED BY POLICY rather than by tooling — `doc/ai_collaboration_protocol.md` §7.1 closes **six units** (not the directory) to AI edits outright, liftable only by the module owner's named authorization relayed by Christian and recorded in a devlog. Elsewhere the rule is unenforced; a pre-edit tripwire on files whose module documentation cites a preprint or "in preparation" would generalize it. `Status: proposed (policy-enforced for the six-unit cohomology core only)`
**Applicability:** Any code whose algorithm is unpublished, in draft, or deliberately modified from a published one. **Not** a general license to avoid unfamiliar code — it is specifically about the *near miss*: the hazard peaks where the code closely resembles something famous and differs in a small deliberate way. Genuinely alien code produces visible uncertainty, which is safe because it gets flagged.
**Exceptions:** Consulting the primary source and deriving independently are the paths that have worked. Removing "dead" code is **not** one of them: the one clean-looking instance was human-guided, and the thing that made it safe — knowing that the greedy rectifier is load-bearing inside the SPFA path — is exactly what an AI reading the same code does not have (INC-564). The rule bites on reasoning alone, and "this looks redundant" is reasoning.

# ANNEX — candidates not promoted to cards

One line each. These met the ≥ 2 incident threshold but lost on spread or budget. They are
recorded so a second pass can find them, not because they are wrong.

- **N5 cost attribution unmeasured** (2 incidents, 2 entries, 2 threads) — *profile before
  attributing a cost.* INC-369 (an O(N²) container scan attributed to "graph construction is
  inherently expensive"), INC-394 (a 157 MB HDF5 read blamed for what was under 2% of the
  run). Enforceable as PROC; lost on count alone and is a reasonable first promotion.
- **N6 backend asymmetry** (2 incidents, **1 entry**, 1 thread) — *a green result on one
  matrix backend does not generalize to the other.* INC-249: reference tests indexing
  `Matrix::data()` linearly, correct under Armadillo and wrong under Blaze's padded columns.
  Both incidents come from a single entry, which is exactly the inflation the entries column
  exists to expose. The underlying rule is already in `CLAUDE.md` under the column-major
  backend caveat, and it is folded into L-11's enforcement rather than standing alone.
- **From the source feedback document, with no seed counterpart** — *cheap reproducibility as
  a deliberate investment.* Once the bug was reproducible in a sandbox from a state file, each
  hypothesis cost about four minutes to kill and most wrong ideas died the same hour. This is
  a positive practice rather than a failure mode, so the failure-shaped seed list has no home
  for it; it belongs in `doc/ai_workflow_best_practices.md`, not here.
- **From the source feedback document, suggestion (b)** — *treat your own hypotheses as
  adversaries by default.* Pre-registration and refutation-seeking came from Christian's
  process, not from the model's defaults. Partly covered by L-01 and L-03; recorded because
  the general habit is broader than either card.
- **36 uncertain assignments** were flagged by the assigners themselves during clustering, each
  recorded with its rival cluster and the question that would settle it. None of them changes
  whether a card exists.
