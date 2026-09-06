# AI-Assisted Development: Practices That Earned Their Place {#doc_ai_workflow_best_practices}

**Date:** 2026-08-11
**Purpose:** The working practices of the BELFEM multi-AI development method, restricted to
those that were adopted, used in anger, and observably changed outcomes — with what
happened after each one, and what each one costs.

This is a companion to two other documents and duplicates neither:
`METHODOLOGY.md` maps our terms onto published vocabulary, and
`doc/ai_collaboration_protocol.md` is the normative rule set for this repository. This file
is the experience report: which rules paid off, why, and what a team adopting them should
expect. It is written to be readable outside this project.

**Context it came from.** A single-maintainer HPC finite-element framework (C++17, MPI,
~500 kLOC) developed with one primary AI (broad exploration and implementation) and two
auditor AIs from different vendors, over roughly five months and ~250 recorded sessions.
The domain matters: results are numerical, most defects are *silently wrong* rather than
loudly broken, and a full validation run costs minutes to hours. Practices that assume
cheap green/red feedback do not transfer unchanged.

---

## 1. Reviews and evidence

### 1.1 Independent audits from *different vendors*, run blind

The primary AI's findings are audited by two models from two other vendors, in parallel,
each unable to see the other's response ("jury"); a sequential variant where each auditor
sees the previous findings ("relay") is used when the goal is to deepen one line of
argument rather than to sample independent ones.

**What happened after adopting it.** The value showed up almost entirely in
*disagreement*, not agreement:

- A proposal that all three models endorsed was corrected by a single dissent that
  identified 52 existing code sites depending on the old behavior — the design shipped
  with the opposite default.
- An auditor found a second, unfixed overload of a function whose sibling had just been
  repaired; no one else had looked for a second overload.
- An auditor derived a tangent by hand and found it off by a factor, in eight overloads.
- An audit of *newly written hardening code* caught a zero-size allocation that would
  abort in a reachable production path.

Meanwhile several claims that all three models agreed on were later refuted by reading the
source: a suspected library-state bug, an ID-collision diagnosis, and a file-format
diagnosis that the maintainer contradicted from experience. Concurring audits are a
sampling of correlated priors, not a proof.

**Practical rules.** Use different vendors — the point is decorrelated failure modes, and
two instances of the same model fail the same way. Read the *deltas* first. Never record
unanimity as evidence.

### 1.2 Pre-registration: findings are locked before any auditor is invoked

The primary AI writes its complete findings — mechanism, file:line evidence, confidence —
into the shared thread *before* any auditor runs, and never edits that entry afterwards.
Resolutions are appended, not merged back.

**Why it earned its place.** Without it, agreement is unfalsifiable: a model that has seen
the auditor's answer will produce a consistent story every time. With it, three things
become measurable: whether the auditors found anything the primary missed, whether the
primary's confidence was calibrated, and whether a round was worth running at all. It also
makes the record honest — several sessions are on file where the pre-registered analysis
was simply wrong, which is information the next session needs.

### 1.3 Auditor output is testimony, not verdict

Every auditor claim carrying a file:line citation is re-checked against the actual source
before it is accepted, and labeled CONFIRMED / REFUTED / UNVERIFIABLE with the evidence
quoted.

**What happened after adopting it.** Roughly one auditor finding in five does not survive
verification — including confidently stated, mechanically plausible ones. Refuting them is
cheap (a source read) and skipping the step is expensive: an unverified finding becomes a
"fix" to code that was correct, which is the most costly failure mode of AI-assisted
review.

### 1.4 The evidence ladder: "reviewed" is not "verified"

A finding is *reviewed* when a static audit is complete. It is *verified* only when an
executable gate has passed. Written records must not say "verified" for read-only work —
they name the highest rung actually reached:

1. End-to-end reproducer (a full run reproduces, or stops reproducing, the behavior)
2. Focused regression (a permanent test locks the claim)
3. Compile / link
4. Numeric probe (a one-off computation confirms a value or identity)
5. Static source trace
6. Literature consistency
7. AI reviewer agreement

Lower rungs support, never replace, higher ones — and reviewer agreement is the *bottom*
rung by construction.

**What happened after adopting it.** Sessions that would previously have concluded at
"three audits agree" went on to execute something, and execution repeatedly changed the
answer:

- A new test suite passed for the wrong reason: at realistic step sizes the wrong scheme
  and the right one differed by less than round-off. Rescaled, it separated by fifteen
  orders of magnitude. A static review could not have found this.
- A feature built specifically to detect silently stale caches shipped with a silent miss
  of its own, found by running it against three real input decks.
- Hand-derived test expectations were wrong because an unqualified function call bound to
  a C library symbol rather than the project's own.

The ladder's second effect is negative and just as valuable: it ends rounds. Which leads
to —

### 1.5 A stop condition for review rounds

A further audit round on the same subject requires at least one of: a material code change
since the last round, new experimental evidence, reviewer disagreement on a load-bearing
point, or a safety boundary (ownership, parallelism, persistence, ABI, formulation).
Otherwise the next step is the executable gate, not another review.

**Why.** Review rounds are cheap to start and produce output every time; without a stop
condition they crowd out the runs that would actually settle the question. This rule
converted "let's get one more opinion" from a default into a decision.

### 1.6 Calibrated confidence on every non-trivial claim

Claims carry **high** / **medium (~N %)** / **low**. High means verified in code or
literature; medium means likely but unverified; low means hypothesis.

**What it changed.** Auditors spend their effort on the weak claims instead of re-scanning
everything, and the record stays useful months later. The subtler benefit: it makes it
legitimate to *split* confidence instead of averaging it — e.g. "that the wrapper handles
this return code is high (read at the cited lines); that the library actually emits it is
medium (source not installed)". Averaged into one number, that claim would have been
either overstated or useless.

---

## 2. Division of labor

### 2.1 Investigation is read-only by default; source edits require explicit approval

Writing to the scratch channel, the task plans, and the session log is always allowed.
Everything else waits for the maintainer to say go.

**What happened after adopting it.** Diagnosis sessions produce a *mechanism and a
proposed fix* rather than a speculative patch, which means the fix can be rejected without
anything having to be reverted — and it frequently was rejected in favor of something
better. It also makes the primary AI's incentive honest: nothing is gained by rushing to
an edit.

### 2.2 Physics, formulation and design are adjudicated by the human, never by model vote

Models settle questions of code. They do not settle questions of physics or of what the
project is for.

**What happened after adopting it.** The single cheapest source of refutations in the
project is the maintainer's one-line correction. Repeated pattern, several times a month:
an elaborate multi-model analysis is dissolved by a sentence containing information that
exists nowhere in the repository — that a code path had never worked and never needed to,
that a configuration had always crashed on a third-party library, that two directories
contained different models entirely. In one case the maintainer replaced a jury-designed
hashing scheme with "just check the numbers", which was simpler, sufficient, and had
survived none of the review because nobody had proposed it.

**The transferable rule:** route design questions to the human *early*, before the models
have built a structure worth defending.

### 2.3 Role separation between the AIs, stated explicitly

Primary: exploration, implementation, documentation, literature routing. Auditors:
standard compliance, memory safety, parallel correctness, edge-case logic, independent
re-derivation — explicitly *not* architecture decisions. A third voice is used to break
ties and to supply a different causal reading, not to add a vote.

**Why it matters.** Without stated scope, auditors drift into redesign, and their output
becomes a competing plan rather than a check. With it, an audit is a bounded artifact that
can be verified item by item.

### 2.4 The human owns the executable gates

In this project, runs and builds are the maintainer's: a shared build tree and hour-long
validation runs make unattended builds actively harmful. So sessions are written to end at
a *named gate* — the exact command or deck that will settle the question — rather than at
a conclusion.

**Consequence to plan for.** The bottleneck moves to the human. At the time of writing,
ten of fourteen release-blocking items need a *run*, not code. That is a healthy state
compared to the alternative, but it means batching runs deliberately, and it means an
AI-side "done" is never the end of an item.

---

## 3. Memory: what gets written down, and where

### 3.1 Two tiers, separated by audience, with distillation as the price of deletion

- **AI-only, ephemeral:** one file per topic for the AI-to-AI exchange. Machine-parseable,
  never committed, swept on a schedule.
- **AI+human, durable:** a dated session log ("what changed and why", backward-looking) and
  task plans ("what needs to happen", forward-looking). Both in-repo, curated for a human
  skimming for decisions.

The tiers are joined by one rule: before scratch is swept, its conclusions are *lifted*
into the durable record. **If it is not in the session log or a task plan, it is lost.**

**What happened after adopting it.** The exchange can be garbage-collected aggressively
without anxiety, which keeps it usable; and the durable corpus stays navigable because it
contains conclusions rather than transcripts. An early failure mode is worth naming: for a
while every audit pass produced its own session log, and a chain of eleven of them
documented one investigation. Chains like that should be collapsed into their closing
entry.

### 3.2 The session log records what was *rejected*, not only what was built

Superseded designs, refuted hypotheses, and approaches ruled out — with the reason — are
first-class content.

**Why this is the highest-value habit in the whole method.** AI sessions have no memory.
Without a written record of rejected approaches, the same idea is re-proposed indefinitely,
and each re-proposal costs a full round to re-refute. Entries in this project explicitly
carry lines of the form "rejected on file size — do not re-propose", and they work.

### 3.3 A compression layer over an append-only log

Session logs are append-only and grow without bound, so current state lives in two derived
artifacts: one ≤1-page **campaign page** per active workstream (accepted design, branch,
last passing reproducer, open P0/P1, superseded approaches) and one **debt register** table
of open items with a "blocks the release?" column. A session that opens or closes an item
updates both in the same session.

**What happened after adopting it.** Onboarding a session into a months-old workstream
went from reading a dozen dated entries to reading one page. But note the failure mode
below — derived pages rot, and they rot silently.

### 3.4 Periodic currentness sweeps of the forward-looking plans

Every few weeks, all active plans are re-read against the actual tree.

**What happened after adopting it.** Every sweep found real drift: plans that were finished
but never closed, status lines that actively misdescribed the code, a campaign page that
said "uncommitted" for five days after the work had landed, and a tracked item whose
reproducer had been deleted five days earlier. One generalization earned its keep:
**sequential plans survive; leapfrogged plans rot.** Every stale plan had been overtaken
*sideways* — superseded by a different plan, or by a decision taken elsewhere — rather than
abandoned. Nobody closes a plan that was never rejected.

### 3.5 Anchor references by searchable token, never by line number

Machine-readable contracts anchor on a greppable string. Line numbers may appear in prose
for human navigation, but they are advisory, not truth.

**The measurement that forced the rule.** In one reference document, 87 `file:line`
citations were audited: 26 were stale and 24 had shifted by exactly the same offset because
one file had grown. Crucially, the document and the code had been edited *in the same
commits* — the discipline was being followed and still failed, because a line number read
from a working copy is already wrong by the time it is committed.

**The generalization.** The dominant defect class in a maturing AI-assisted corpus is not
faulty reasoning; it is **pointers that were correct when written and were never
re-followed** — stale line numbers, a tracker ID naming the wrong item, an enumeration that
listed two of three call sites. None of these are visible to a reader who does not re-check
them, and all of them survive review rounds.

### 3.6 Mechanical checking in the code → document direction

A script verifies the checkable claims in the convention documents — build flags, make
targets, executable names, module build status, existence of named symbols — against the
tree, and is run after editing those documents or the build configuration.

**Why this direction specifically.** Human discipline reliably keeps *prose* current; it
has. What it cannot catch is a new key or flag added by someone who does not know the
document exists. That gap is structural, so it needs a mechanical pass, and the pass must
run from the code toward the document, not the reverse.

**Its limit, learned the expensive way.** In one investigation every mechanically checkable
fact was verified and correct, and the conclusion was still wrong, because the missing
premise — that one third-party configuration had never been validated and crashed a
dependency — was written down nowhere. *A mechanical check catches stale facts; it cannot
supply a posture nobody ever recorded.* The fix was to write the posture down and add a
configure-time guard that enforces it.

### 3.7 Derived documents rot faster than authoritative ones — so point, don't restate

Two instances, both costly. A test-writing guide inferred a rule from the authoritative
design document, inferred it wrongly, and the error survived for months precisely because
the authoritative document was right — nothing contradicted it where anyone looked. The
result was 72 tests that silently never ran in release builds. Separately, a bootstrap file
had grown a copy of another repository's directory tree, which then drifted; deleting the
copy was the fix, not repairing it.

**Rule:** a document should state facts it owns and link to the rest. Restating another
document's content creates a second copy that will diverge without warning.

---

## 4. Domain practices that generalize

### 4.1 Literature-first — but conditionally

Consulting primary sources is *mandatory* for algorithms, formulations, numerical methods,
correctness validation, and any explanation of *why* something is implemented as it is. It
is explicitly *not* required for style fixes, build changes, or mechanical patches.

**What happened after adopting it.** A verification pass comparing implementations against
their source papers found four production defects of a class that no test and no build
would ever surface: two swapped coefficient arrays, a missing crossover term, a sign typo
in one of forty fitted terms, and a double-counted reference state. Each was quantified
against the paper's own constraints before being reported.

**The conditional half is load-bearing.** Making literature routing universal turns it into
ceremony, and ceremony gets skipped wholesale — including where it mattered.

### 4.2 Prefer making a silent failure loud

The recurring defect shape in this codebase is silence: a matrix silently transposed by a
lazy expression template, a stale cache silently reused, plugin registrations silently
becoming no-ops after an ABI change, a test binary aborting mid-suite so that 72 skipped
tests reported as one failure. A large share of the highest-value changes were not new
behavior at all — they converted a silent wrong answer into a loud one.

**As a review question:** for any defect found, ask what *would* have made it announce
itself, and whether that guard is cheaper than the next occurrence.

### 4.3 Changes to defaults need a run gate, not a review gate

Several default changes passed review, landed, and were rolled back within a day once a
real deck exercised them — a scheme default, an acceleration default, a newly added guard
that aborted every parallel run. Reviews assess whether a change is *correct*; only a run
assesses whether it is *safe for the decks people actually have*.

### 4.4 Probes: get executable evidence without the full gate

Where a build is expensive or owned by someone else, a standalone probe — a small program
linked against the already-built libraries, or a one-off numeric computation — buys a rung
on the evidence ladder cheaply. Probes found the near-vacuous test, the parser bug in a
new detection feature, and several wrong hand-derived expectations. Probes are deleted once
their question is answered; a probe left in the tree becomes an unowned diagnostic.

---

## 5. What this costs, and what it does not fix

- **Token and wall-clock cost is real.** A full jury round on a non-trivial subject is
  substantial. The stop condition (§1.5) exists because the method's natural failure mode is
  reviewing instead of running.
- **Agreement is seductive.** Three concurring audits *feel* like proof and are not. This
  has to be written into the rules, because it will not be resisted in the moment.
- **The record is only as good as its sweeps.** Every derived artifact — campaign page,
  register, index — will drift. Budget for the sweeps or do not build the layer.
- **Nothing here substitutes for tests and runs.** The method is a way of producing
  well-founded hypotheses fast and refuting bad ones cheaply. The last rung is still an
  executable gate, and in a numerical code that gate is slow and often human-owned.
- **It does not scale down.** For a typo or a build tweak, all of this is overhead; the
  protocol says so explicitly, and that exemption is what keeps the rest credible.

---

## 6. Minimum viable adoption

If a team wants one afternoon's worth of this, take these five in order:

1. **Write findings before asking anyone** — pre-registration costs nothing and makes every
   later signal interpretable.
2. **Use a second model from a different vendor, and read the disagreements first.**
3. **Say "reviewed" until something has actually run**, and name the gate that would make
   it "verified".
4. **Keep a dated session log that records rejected approaches**, not just changes.
5. **Anchor every reference to something greppable.**

Everything else in this document is refinement on top of those five.
