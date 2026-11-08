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

Everything else open is `[P]` (4 rows) or `[W]` (1 rows), with no untagged rows remaining.
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
| DR-200 **[CODE]** **[P]**  | linalg: **the `Matrix`/`Vector` wrapper pays for two defects of its own** — six `virtual` members nothing derives from (2.64× Blaze / 1.28× Armadillo on a 3×3 kernel through `const Matrix &`, parity once removed) and overload-set gaps that force temporaries (Armadillo `Matrix` `+=`/`-=`/`*=` and ctor from expressions, Blaze `Vector` `+=`/`-=` from any vector expression, Blaze `Matrix` compound ops only under `BELFEM_GCC`), which 110 `matrix_data()`/`vector_data()` hatch lines in 13 consumer files work around by hand. Plus the release check policy: Armadillo checks bounds and shapes in release, Blaze does not (decided 2026-11-07, Christian: on in debug, off in release; Blaze's shape exceptions stay). Plan `todo/linalg_wrapper_hardening_plan.md` v3 (two plan-jury rounds, Codex astra + Grok), rounds A de-virtualize + `final` → B check policy (with R2.0: `fn_AR_eigen.hpp:33` discards `eig_gen`'s status, an out-of-bounds read once checks are off) → C overload sets (counter + numerical gates, `strans` for complex, O12 runtime shape branch) → D hatch retirement (BDF2–5 tests first). Investigation: `todo/linalg_wrapper_design_investigation.md` v6 | P2 | **open — rounds A, B and C landed and verified 2026-11-08** (Christian's `make check`, both backends, every suite green): de-virtualized `final` wrappers; checks follow the build mode; complete expression overload sets on both backends with the O12 orientation contract, pinned by the executed allocation-parity counter, the overload value tests and the BDF1–5 assembly baseline; dead matrix `dot` retired (O14). Round D: inventory done (110 lines, 41 are `Cell`/HDF5/`Spline` accessors and stay); R4.1 (`cl_IWG_Timestep.cpp`, 70 live lines) and R4.2 (kernel/postproc/eigen, 12 lines) replaced and code-juried 2026-11-08, no P0/P1 — the eigen residual's Blaze branch was `sqrt` of the 1-norm and now matches Armadillo's 2-norm (a Blaze stop-test change, Grok); Blaze `make check` green, Armadillo owed, uncommitted; R4.3 with belfem-5e | `tests/linalg` tripwires (`is_polymorphic`, `sizeof`) + `make check` both backends, both modes; round C: `test_LinalgAllocations` parity + positive control | no |
| DR-193 **[MIXED]** **[P]**  | examples/tape_quench: **the usermat plugin's property tables end below the temperatures the T5 quench deck reaches.** Active tables end at 199 K ( `k_Hast` ), 282 K ( `cp_HTS` ), 300 K ( `cp_Cu`, `k_HTS` ), 303 K ( `cp_Hast` ), 390 K ( `rho_Hast` ), endpoint held above, copper ρ clamped at 1358 K; the envelope is the minimum over every active property, 199 K, reached at ≈ 47 ms, so the 1770 K hot spot that both Picard and the hybrid reproduce is the model outside its tables. Re-filed from DR-105 residue ( 3 ) | P3 | open — ruled 2026-10-03 ( O3 ): report results against the 199 K envelope; table extension is a separate task ( NIST cryogenic and high-temperature compilations for Cu and Hastelloy C-276, DOI beside each number ) | the reconstructed T5 deck, `cmake-build-claude/tape_quench_dr105`, past ≈ 47 ms | no |
| DR-201 **[CODE]** **[P]**  | build/tests: **with `USE_NONFREE=ON`, `make check` runs the nonfree tests without building them.** ctest runs every registered test, but the `check` target's dependency loop (`CMakeLists.txt`, the `foreach( TESTDIR core containers … circuit )` block) lists only the open-source suites, so `test_nonfree-suite-a`, `test_nonfree-suite-b` and the `the-nonfree-python-extension` extension are never built by it. Observed 2026-11-07 19:11 in `cmake-build-debug` (`USE_NONFREE=ON`, `USE_PYTHON=ON`): `nonfree-suite-a` and `nonfree-suite-b` "Unable to find executable", `the-nonfree-python-extension` "No module named 'the-nonfree-python-extension._the-nonfree-python-extension'" (`the-nonfree-python-extension_REQUIRE_EXTENSION=1` turns a missing extension into a failure, `nonfree/python/CMakeLists.txt`). With `USE_NONFREE=OFF` none of the three is registered (the `nonfree/` subdirectory is never added, `CMakeLists.txt` `if ( USE_NONFREE )`), which is correct. A plain `make` builds `test_nonfree-suite-a`/`test_nonfree-suite-b` (they are in ALL) but not the extension | P2 | open — filed 2026-11-08 by belfem-5e, diagnosis read-only, nothing edited. Proposed fix, not yet approved: `Add_Test.cmake` appends `test_${TESTNAME}` to a global property (e.g. `BELFEM_NONFREE_CHECK_TARGETS`) when the test source lies under `${BELFEM_NONFREE_SOURCE_DIR}`, `nonfree/python/CMakeLists.txt` appends `the-nonfree-python-extension`, and the root `check` block adds `add_dependencies( check … )` over that list — no nonfree target name in the open-source `CMakeLists.txt`, empty list when nonfree is off. Alternatives offered: `if( TARGET the-nonfree-python-extension )` in the root loop, or keep `python`-labelled tests out of `make check` (`ctest -LE python`) | `make check` in a tree with `USE_NONFREE=ON` (and `USE_PYTHON=ON` for the the-nonfree-python-extension case), without a prior `make` | no |
| ~~DR-203~~ **[CODE]** **[P]**  | ~~mesh: **`TensorMeshBSpline.{Quad16ReproducesCubics, Quad16TranslatedAnisotropic, Hex64ReproducesCubics}` fail: every node reports `number_of_sources() == 0`**~~ (`tests/mesh/test_TensorMeshBSpline.cpp:82`). The tests were added 2026-10-22 (`9abf660`) with the devlog recording "test not run" (`dl20261022_guide_sweep_code_fixes.md:47`); they build the mesh in memory (`Mesh( 3, tNumNodes, tStep, tOrigin )`), read no file, so this is independent of DR-202. Source population is `cl_TensorMeshFactory.cpp:1265` (`tNode->set_sources( tSources, tWeigts )`). Seen in the 2026-11-07 `make check` on this Mac; not caused by the DR-200 round-A diff (no semantic change; clean rebuild) | P2 | open — filed 2026-11-07; **root cause found 2026-11-07 (read-only):** the tensor-mesh constructor `Mesh::Mesh( aOrder, aNumNodes, aStep, aOrigin )` (`cl_Mesh.cpp:97-120`) never sets `Connectivity::Compute`, unlike its two siblings (`:92`, `:132`), so the `finalize()` inside `populate_tensor_mesh` (`cl_TensorMeshFactory.cpp:85`) skips `connect_nodes_to_elements()` (`cl_Mesh.cpp:753-756`); `create_bsplines` → `connect_nodes_to_bsplines` then loops `tNode->number_of_elements()` = 0 (`:1178`), collects no control points and calls `set_sources` with empty containers. `tests/mesh/test_Mesh.cpp:250-264` documents the missing flag and works around it with `unfinalize(); set_connectivity( Compute ); finalize()`; the B-spline test does not, and no cubic tensor mesh built through the constructor gets its node weights. Fix candidate: set `Connectivity::Compute` in the tensor constructor before `populate_tensor_mesh` (one line; node-to-element, node-to-node and facet links as the siblings compute — edges and faces stay lazy, `finalize_edges` acts only once edges exist) and retire the `test_Mesh.cpp` workaround; owes a jury round (`src/mesh`). **Jury 2026-11-07 (Codex astra/high + Grok 4.7/high): no P0/P1; two P2s to fix with the diff — the control-point connectivity passes inside `finalize()` (`cl_Mesh.cpp:784-788`) run before `create_bsplines` creates the control points and set their flags on empty containers (rebuild after the B-splines); `aComputeConnectivities = false` with `aOrder > 1` silently restores the defect (refuse at construction). Pre-existing, logged: `FacetToElement` set after an empty facet loop and never reset by `unfinalize()`; `mEdgesAreFinalized`/`mFacesAreFinalized` true while empty.** **Owner: Claude (linalg session), assigned by Christian 2026-11-07; fix in the tree (`cl_Mesh.hpp`, `cl_Mesh.cpp`, `test_Mesh.cpp`, +12 −28), syntax-gated with the build's flags in both modes, jury running (`tmp/ai_exchange/review_dr203_tensor_mesh_connectivity.md`); gate `test_mesh` is Christian's.** **RETIRED 2026-11-07: fix committed; Christian's `make check` passed `test_mesh` on Blaze (VTK on) and Armadillo — the three B-spline cases, the connectivity cases without their workaround, and the new refusal case. Pre-existing bookkeeping gaps Grok recorded (`FacetToElement` sticky, `*AreFinalized` on empty containers) stay as a [W] note here, not a new row** | `test_mesh --gtest_filter=TensorMeshBSpline.*` green on both backends | no |
