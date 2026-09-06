# Devlog 2026-08-31 — Cohomology Core Closed to AI Edits

**Date:** 2026-08-31
**Purpose:** Record the ruling that bans AI edits to the cohomology core, its exact scope, and where it is written down
**Module:** src/homology (policy: cross-cutting)
**AIs involved:** Claude (author), Codex `gpt-5.6-terra`/high, Grok `grok-4.6`/high (jury round)
**Claude Confidence:** high on the ruling, scope and wiring; medium on the training-data premise (not checkable); the first draft's supporting examples were wrong and are corrected below
**Codex Audit Confidence:** high (REVISE — 3 high, 1 medium; all confirmed)
**Grok Audit Confidence:** high on tree-checkable claims; low on the human-history claims, which it correctly declined to verify
**Literature References:** Pellikka et al. 2013; Mrozek & Batko 2009 §5 (the empty-cell trick — see the round-2 correction below, it is not what an earlier draft said it was); Giard et al., *Generalized Pellikka algorithm for cohomology computation*, in preparation (drafted 2025, resumed 2026; not in `./literature/`)
**Verification:** static source trace only — `check_doc_claims.py` 37/37 on branch `main`, which does not cover prose and does not corroborate any claim in this entry. **Reviewed, not verified.**
**Exchange thread:** `tmp/ai_exchange/review_cohomology_edit_ban.md`

## Summary

Christian and Gregory Giard, who owns `src/homology`, ruled that no AI — Claude, Codex, Grok,
or a successor — may edit the cohomology core. The ruling was prompted by DR-29: the module's
own developer read the analysis and could not extract its main message, and Christian's own
experience is that he "has broken it and had to restore it often enough in the last few
months." The policy is stricter than the ordinary edit-safety default, which the user can lift
per session; this one cannot be lifted that way.

## Scope, as ruled

Banned, in `src/homology/`: `cl_Cohomology`, `cl_Homology`, `cl_SimplicialComplex`,
`cl_Chain`, `cl_Cochain`, `fn_Smith` — both `.cpp` and `.hpp`, six units, twelve files.

Not banned, ordinary edit-safety rules apply: the rest of the directory (`cl_CutFactory`,
`cl_CutProcessor`, `cl_CutProcessorManual`, `cl_CutData`, `cl_CutSet`, `cl_InterfaceProcessor`,
`cl_BeltedTree`, `cl_Topology`, `en_CutAlgorithm`, `CMakeLists.txt`), the module's `doc/`, and
every caller outside the directory. `cl_Topology` was checked against the boundary before being
left out: it classifies mesh regions by domain type and carries none of the chain-complex
algebra.

Analysis stays open and stays wanted — read, trace, run, and report into `./devlog/`, a
register row, or `./todo/`, addressed to Gregory. What is banned is touching the source and
self-starting a fix. A fix plan in a register row is a proposal for the module owner, not a
queued task.

The override is Gregory's named authorization relayed by Christian, plus a devlog line naming
who authorized what. Christian's ordinary edit approval does not reach these files.

## Key finding — the reason is novelty relative to training data, not difficulty

Christian's framing, and it is the load-bearing part of the policy. What runs here is a
*modified* Pellikka reduction — Gregory's own work, and **unpublished**. The first draft of this
entry called it "published work"; that is false against every citation in the tree (both
auditors), and the correction *strengthens* the argument — a draft in preparation is better
ground for a training-data claim than a publication would be, because there is no published
implementation to have learned from.

**Provenance settled by Christian, 2026-08-31, after the audit rounds flagged it as
irreconcilable.** The work was drafted in 2025, shelved, and resumed in 2026 — which is why
`cohomology_algorithms.md` said "2025 preprint" and `dl20260830_dr29_coreduce_analysis.md:9` said
"2026 draft": one work, two stages, not two citations. **And the author is Giard**, as the source
headers say everywhere (`src/homology/*.hpp`, "Developers: Christian Messe, Gregory Giard"). The
module documentation had spelled it **Giarda** in nine places across four files; all are corrected
in this session, and `doc/literature_references.md` — which carried no entry for the work at all —
now has one, with a note that a reappearing "Giarda" is a typo.

That inverts the usual situation. Given BDF5, RCM or a CG loop, a model has read thousands of
implementations and carries a strong, largely correct prior about what the code should look
like; a real deviation stands out as one and the prior is an asset. Here the nearest object in
the weights is *textbook* Pellikka, so every deliberate modification reads as a defect, and
"fixing" it regresses the algorithm toward the published version it deliberately is not. The
model cannot separate "this differs from the algorithm I know" from "this is wrong", because the
algorithm it knows is the wrong reference.

**The canonical case — corrected against the source, and the correction matters.** The first
draft wrote it as Gregory phrased it in conversation: "for cohomology, and not for homology, a
node is removed from the simplicial complex so the process can start." Read against the live
function, that compresses two things. `coreduceOmit()` (`cl_SimplicialComplex.cpp:1236-1257`)
first runs `pCoreduce(0..2)`, then **loops**: while any 0-cosimplex remains, remove *one*
0-cochain, zero its coboundary's references, and re-run the entire `pCoreduce(0..2)` cascade.
One node per stall, each removal unblocking a further coreduction. It answers the problem
Mrozek & Batko 2009 §5 poses — a simplicial complex admits no elementary coreduction pair at all,
since a boundary is never of size one — but it is **not** Mrozek's answer, and the first draft of
this entry said it was. The paper adds ∅ as a simplex of dimension −1 carrying every vertex in
its coboundary, one per connected component (`mrozek2009.txt:444-464`); it licenses no 0-cell
omission, and `coreduceOmit` adds no −1 cells. *(Correction found in round 2, and it is the
most on-theme error of the session — see below.)* `reduceOmit()`
(`:939-961`) is the strict dual on the chain side and removes **top-dimensional** chains, which
is what "and not for homology" means. Grok read the original sentence as claiming homology has
no omit step at all and refuted it; that refutation is itself wrong — but **it was produced at
high effort with the source open, which is the finding.** A sentence a careful auditor inverts is
a failed sentence whether or not it is true. Also worth knowing: `pellikka2013.txt:164-165`
defers the reduction details to its own reference [28], so citing Pellikka 2013 for the omit rule
was not supported either.

**What an agent actually does with that loop is on the record.** On 2026-07-03 the proposal was
to batch the omission — it removes every 0-cochain in the end, so why re-run the cascade each
time? Struck as unsafe: batch omission leaves empty-boundary edges that `pCoreduce`'s
`size() == 1` test can never remove, and the junk survives into cocombine as spurious generators
(`devlog/dl20260703_coreduce_risk_reeval.md:25-27`). The loop looks trivially hoistable. It is
not. That is the failure mode in one line — code that compiles, runs, and is wrong — and it is
attested rather than imagined, which is why it replaced the original example.

Two consequences worth stating separately. **A review round does not cover this**, because the
errors are *correlated*: models missing the same prior fail the same way, so a majority verdict
carries no more weight than a single voice. The first draft put this too strongly — "three
models will concur, confidently, on the same wrong reading" — and the record refutes it: on
2026-08-30 Claude and Grok concurred on a false defect and **Codex broke the tie** with a Schur
derivation (`dl20260830_dr29_coreduce_analysis.md:41-46`), and in this session's own jury both
auditors found real errors in this document. The defensible claim is Grok's wording: *do not
treat concurrence as coverage* — not *concurrence is inevitable*. And **the generalizable test is
novelty, not difficulty**: a hard but widely published algorithm is safe ground, a simple
unpublished one is not.

That is also why this is a scope boundary and not a directory ban. The cut plumbing carries no
such property — it is ordinary bookkeeping over the core's output, the kind of code the priors
do cover. **The first draft of this entry cited DR-111 and "the rectification work" as the
evidence for that, and both were wrong** (Grok, confirmed): DR-111 is `ThinShellFactory`
clockwise winding in `fem/kernel` (`todo/dr111_quad4ts_winding.md:1-7`), not homology at all,
and rectification (`clean_spfa`) lives in `cl_Cohomology.cpp` — *inside* the ban. The IDs that
actually fit are DR-107 (CutFactory hygiene), DR-108 (`CutSet::create_duplicates`) and DR-133
(`CutFactory::orient_terminal_curves_2D`; note that the landed sign fix in that row was
`Homology::reorient_generators`, which is banned).

This is not a quarantine of a broken module. Gregory's position on DR-29 is that the module has
not been an issue for a while, and that the non-unit-coefficient behaviour is a
thin-cut/mathematical obstruction rather than a defect — a coarse mesh where a conductor loops
under itself admits no unit representative. That is separate from how the implementation handles
it, which `src/homology/doc/thin_cut_nonunit_rectification.md:22-30` documents as a real
limitation ("the current code cannot handle this … the unit-only restriction is structural"),
and `cl_Cohomology.cpp:880-886` hard-errors on any non-unit coefficient that survives pocket
removal. Obstruction and handling are two claims; the first draft ran them together (Codex).
The ban protects working code whose correctness conditions live outside it.

## Changes Made

- `doc/ai_collaboration_protocol.md` — new §7.1 "Protected Module — the Cohomology Core"
  (the authoritative text: file table, exclusions, rationale, what stays allowed, override)
- `AGENTS.md` — Edit Safety Rule gains the ban and points at §7.1
- `CLAUDE.md` — homology bullet under Specialized Modules; new line in the Minimum Session
  Compliance Checklist beside the read-only rule
- `doc/lessons_learned.md` — Layer 1 trigger index gains "Before editing anything in
  `src/homology/`": one heading + four bullets, plus a failure-signature row, taking Layer 1 to
  **224 lines** after the round-2 corrections
  against its 250 budget (heading `:19` to the LAYER 2 marker). The first draft of this entry
  said "2 entries, 6 lines" and "210 lines" — both wrong, both caught by both auditors. Cause,
  recorded because it is the generalizable part: the counts were measured after the first edit
  and never re-measured after the rationale rewrite added a bullet. An L-08 miss inside the file
  that exists to prevent L-08 misses
- `CLAUDE.md:53` — the compliance checklist described Layer 1 as "~140 lines"; now "~220"
- `src/homology/doc/README.md` — banner above the Overview, so module-first navigation hits it
- `todo/debt_register.md` — DR-29's status column records the ban and marks the row
  read-and-report only; DR-71's row records that the ban **splits** it (below) and its stale
  `cl_CutFactory.hpp:86` citation is corrected to `:91`

`check_doc_claims.py`: 37/37 — which corroborates nothing here. The checker probes build flags,
make targets and named symbols in `CLAUDE.md` and `doc/coding_philosophy.md`
(`scripts/check_doc_claims.py:5-16`); this policy is prose and is invisible to it.

## Lessons Extracted — the precedent sweep behind the ban

Christian asked for the devlogs to be checked for precedents once the jury returned. All 429
entries were searched for the module's core files; six incidents were extracted and filed as
**INC-559 … INC-564** in `doc/lessons_learned_evidence.md`, and distilled into a new Layer 2
card, **L-22 "An algorithm with no published implementation inverts your prior"**. Layer 1 gains
a failure-signature row pointing at it.

**The pattern is consistent and it is not "AI is bad at hard code".** Every failure in the
record is the same shape: a deliberate deviation from a published algorithm read as a defect,
and the "repair" moving the code toward the textbook version.

| # | what happened | what it teaches |
|---|---|---|
| INC-559 | DR-29: Claude **and** Grok concluded a correct Smith-stage matrix was missing a merge projection; Codex refuted it with a Schur derivation (`BA=0` ⇒ `αd+λD=0`, so the live restriction *is* the retracted operator). Grok's alternative fix was discarded | **A majority vote would have corrupted a working kernel.** The minority voice won on a derivation, not on a reading. This is the empirical core of the ban |
| INC-560 | Batch-omitting the 0-cells in `coreduceOmit` — struck as unsafe against Mrozek §5, which licenses one omission *per connected component* | The prior carries the published trick and **not its side condition**. This became the policy's canonical example |
| INC-561 | The in-tree full-rescan `pCoreduce` flagged as a deviation from Mrozek Algorithm 6.1's worklist. The paper reversed the framing — the worklist is the baseline — and Thm 6.2 then dissolved the concern entirely | Reading the primary source corrected the *direction* of the finding and cost it its severity. Deviation ≠ defect |
| INC-562 | The 2026-07-03 conclusion "we already have a worklist in-module" (`pCocombine`) reused as a template, refuted 8 weeks later as a category error: LIFO `pop`, sorting `unique()`, unguarded dereference, no ±1 check | Shape matched, semantics did not. **Conclusions about unfamiliar algorithms decay into false premises faster**, because nothing in the reader's prior contradicts them |
| INC-563 | This document's own first draft: four load-bearing errors, all found by the jury | The policy describing the failure mode exhibited it |
| INC-564 | The two AI edits *inside* the banned core that did land — `759a0801`'s `pGeneralizedCocombine` restructuring and DR-107's `clean_greedy()` deletion — **corrected 2026-08-31 (Christian): neither is the counter-instance it looks like.** The first preserved DR-29 #10's zero-neighbor skip, a latent defect it did not understand. The second was **heavily human-guided surgery**: both algorithms are needed, the greedy rectifier `cohomology::rectify_greedy_sweeps` is live inside `clean_spfa`, and only the standalone member plus the `mFunClean` switch were cut | The safety came from a human supplying the constraint the code does not state. An AI framing it alone — "two cleaners, one is dead" — deletes the live rectifier. **"This looks redundant" is reasoning, and reasoning is what fails here** |

**The discriminator L-22 ends up on is therefore not "did an AI touch it" but "was the
constraint supplied from outside the model's prior".** The corpus contains no case where
re-reading the code harder resolved one of these, and three where the primary source or an
independent derivation did.

**The DR-107 correction is what fixes the boundary, and it narrows the ban's cost rather than
widening it.** On first pass this entry filed the `clean_greedy()` deletion as a clean AI success
— dead code, removed under a double jury. Christian's correction: it was heavily human-guided
surgery, and the guidance was load-bearing. Both algorithms are needed — production is the SPFA
certifier *plus* the greedy rectifier, `cohomology::rectify_greedy_sweeps` called from inside
`clean_spfa` at `cl_Cohomology.cpp:548`, with `Cohomology::clean()` forwarding to `clean_spfa` at
`cl_Cohomology.cpp:189-192` — and
what the round actually deleted was the standalone member and the `mFunClean` function-pointer
switch. The greedy algorithm never moved. An AI holding only the code sees two cleaners and one
dead dispatch; the move that framing suggests is deleting the greedy algorithm and taking the
rectifier with it. So the residual "what still works" list is shorter than this entry first
claimed: **read the primary source, or derive it independently.** Deleting apparent redundancy
is not on it, because "this looks redundant" is reasoning, and reasoning is the thing that fails
here. That makes the ban cheaper than it appeared — it forbids less that was working.

**One refinement the record forces onto the "novelty" framing.** The hazard peaks at *near*
similarity, not at absence of a prior. Genuinely alien code produces visible uncertainty, which
is safe because it gets flagged; code that is 95 % isomorphic to something famous with a small
deliberate delta produces *confidence*. Every incident above sits in the second category. L-22's
Applicability section says so explicitly, because "novel ⇒ dangerous" would misroute a session
to the wrong files.

## Round 2 — the power sweep, and what it found in the lesson itself

Christian ordered a second round on the distilled lesson at raised depth: Codex `gpt-5.6-sol`
xhigh + Grok `grok-4.6` xhigh, jury, full file uncut (`DIFF_CAP=700`; L-22 sits past the default
400-line cap and would have been invisible). Both completed. Full record in the exchange thread.

**The finding of the round, and it is the failure the card describes, committed while writing the
card.** The claim "Mrozek & Batko 2009 §5's empty-cell trick licenses omitting one 0-cell per
connected component" is false. The primary text (`mrozek2009.txt:444-464`) says a simplicial
complex admits **no** elementary coreduction pair — a boundary is of size n+1 or 0, never 1 — and
the remedy is to treat ∅ as an additional simplex of dimension **−1** whose coboundary holds every
vertex, so each vertex becomes a free coface; for several components, one such −1 generator **per
component**. It licenses no 0-cell omission. BELFEM's `coreduceOmit` adds no −1 cells; it removes
0-cochains in a loop. I had mapped the code onto the nearest published trick — precisely what
L-22 warns against — while writing L-22.

**And I did not invent it.** The gloss is verbatim from
`devlog/dl20260703_coreduce_risk_reeval.md:25-27`, is preserved in the catalogue at
`lessons_learned_evidence.md:568` (INC-476), and I carried it into §7.1, into INC-560 and into the
card without opening the paper. A false literature attribution propagated unchallenged for two
months and then into policy — INC-562's shape applied to a citation. **What survives untouched is
the safety conclusion:** the batch-omission argument (stranded empty-boundary 1-cells that
`pCoreduce`'s exactly-one test at `cl_SimplicialComplex.hpp:364` can never pair) is BELFEM's own,
derived from this code in that devlog. Only the citation was borrowed. *(Neither auditor nor I
re-derived the cascade itself; recorded as an open gap, not as verified.)*

**Second self-inflicted citation error, same class.** L-22 and INC-564 anchored the greedy
rectifier at `cl_Cohomology.hpp:135`. That line is the *declaration* `void clean_spfa();`. The
forward is `Cohomology::clean()` at `cl_Cohomology.cpp:189-192` and the rectifier call is at
`:548`. I produced it by grepping the header, seeing `135: clean_spfa();`, and reading a
declaration as a call — then reporting it as verified. Both errors are INC-543's shape, inside the
material distilled to prevent INC-543.

**P0, and the reason this round could not wait: my Layer 1 heading over-banned.** "Before editing
anything in `src/homology/`" is broader than the policy it summarises — §7.1 bans six units and
leaves CutFactory, CutProcessor(Manual), CutData, CutSet, InterfaceProcessor, BeltedTree,
Topology, `en_CutAlgorithm`, `CMakeLists.txt` and the module's `doc/` under ordinary rules. The
bullets underneath were right; the heading is what fires, and it would have blocked legitimate
cut-plumbing work while contradicting the document that outranks Layer 1. L-22's Status line
repeated it. Both corrected to the six-unit boundary, and the override gate now states the full
condition (Gregory's authorization **relayed by Christian and recorded in a devlog**).

**Also fixed:** L-22 had been appended *after* `# ANNEX — candidates not promoted to cards`, so a
sequential reader would take the cohomology card for rejected material — it now sits directly
after L-21. "In every case the deviation was the point" was contradicted by this card's own
INC-561, where the deviation carried no correctness weight and only the paper's Thm 6.2
established that; the rule now reads *treat the deviation as unknown until checked against the
primary source*, since both outcomes occur and guessing which is the failure. The corpus header
said "539 incidents grouped into 20 clusters" where clustering was locked on 537
(`lessons_learned_evidence.md:873-877`). And a compression of mine inverted: "empty-boundary
edges `size()==1` can never remove" reads as though empty-boundary edges have `size()==1`.

**Nine pre-existing defects in other cards were found and deliberately NOT touched** — they are
other people's rules and outside this task. Listed in the exchange thread and summarised for
Christian: L-21's `assert.cpp:294-320` citation points at a contract that moved to
`cl_Communicator.hpp:235-258` (INC-543's shape inside the card meant to stop it, and
`coding_philosophy.md:631` is stale the same way); L-21's proposed grep excludes `.f90` and so
misses the densest remaining direct MPI in `mumpstools.f90` and `parpacktools.f90`, whose status
as an intended exception is nowhere written down; the header's "none is mechanically enforced yet"
is a false universal contradicted by L-08 in the same file, and "each card names its enforcement
target" is false for L-19/L-20; Layer 1 pins INC-555 to "the ARPACK path" when the row is the
MUMPS solver-ID pool; the `git stash` tripwire is unsafe as written in a shared tree; "the run is
right" is stated unconditionally against this file's own stale-binary entries; L-18 is broader
than `coding_philosophy.md:600-627`; L-09 still says "CI-less"; and "`checkout --` is the one git
verb with no undo" is false.

## Open Questions

- The policy is written but not mechanically enforced. A pre-edit hook matching the six file
  stems would make it a tripwire instead of a rule; not built here.
- DR-71 is **split**, not enclosed: defect (2), the unbounded `++it2` at
  `cl_SimplicialComplex.cpp:873-878`, is inside the ban and is Gregory's; defect (1) and the
  "retire the dead branch outright" option are in `cl_CutFactory.cpp:278-281,295-298,329-332`,
  which the ban does not cover. Recorded in the row.
- **Open for Christian and Gregory — the file list is not derived from the stated criterion.**
  `fn_Smith` is textbook Smith normal form, which the policy itself names as the example of a
  *safe* prior, yet it is banned (on coupling grounds: the SNF stage consumes the reduced complex
  directly). `cl_BeltedTree` holds `Cell<Cochain*>` cohomology generators and a
  `SimplicialComplex*` (`cl_BeltedTree.hpp:26-46`), has no literature citation, and is *outside*
  the ban with no boundary check recorded. By the novelty test the list has one member it cannot
  justify and lacks one it can. Not a production concern — the default path is
  `PellikkaGeneralized` — but the criterion and the list disagree.
- **Open for Christian — comment-only edits.** §7.1 says "edits" and "touching the source" and
  does not rule on comments. A peer file (`todo/doc_currentness_fixes.md:31-35`) has already
  inferred "including comments". Under that reading, adding a one-line note at `coreduceOmit`
  explaining the omit — which would *shrink* the hazard this ban exists for — is forbidden.
  Recorded as "treat as banned until ruled" in §7.1 and the module banner.
- `tests/homology/` is outside `src/homology/` and therefore outside the ban. An agent can still
  pin a wrong reading into `test_Cohomology.cpp`'s expectations. Accepted residual (Grok).
- **The policy's central empirical claim is not reconstructible from this repository.** Both
  git-visible reverts in the module — `6df7d7c9` "revert to old Cut Processor" (2025-06-18,
  Christian) and `34c44150` "reverting to the previous CutProcessor" (2025-08-27, Gregory) —
  predate this AI collaboration by about nine months and are both `cl_CutProcessor`, outside the
  ban. Two AI-authored changes *inside* the ban held up: `759a0801`'s restructuring of
  `pGeneralizedCocombine` (the pre-image was read: the zero-neighbor skip DR-29 #10 describes is
  pre-existing, and the claim of preserved semantics is true) and DR-107's deletion of
  `Cohomology::clean_greedy()` under a double jury. The restore episodes Christian cites are
  consistent with working-tree breakage, which git does not record. If the policy is to survive a
  later challenge, one named episode in §7.1 would carry more than the argument does.
- `doc/README.md:12` still describes Layer 1 as "141 lines" and Layer 2 as "18 operating rules".
  **Deliberately not fixed here** — that file is under concurrent modification by the doc
  currentness sweep, whose own `todo/doc_currentness_fixes.md` already queues the correction.
