# Methodology

BELFEM is developed with a multi-AI workflow: Claude acts as the primary development AI
(exploration, implementation, documentation), with Codex and Grok — models from two other
vendors — as independent auditors of its findings and changes. Coordination runs through
ephemeral, per-task exchange files under `tmp/ai_exchange/` (git-ignored); their
conclusions are distilled into the permanent session record in `devlog/` before the
scratch files are swept. The operational rules live in
`doc/ai_collaboration_protocol.md`; the tooling is `.claude/commands/cross-review.md` and
`scripts/cross_review.sh`. The method's distilled output is `doc/lessons_learned.md`: a
trigger-indexed rule layer derived from the session record itself, with every rule citing the
incidents that paid for it.

## Terminology mapping

| Our term | Published term |
|---|---|
| Independent Codex/Grok audits of Claude's findings | **Cross-model adversarial review**; deep lineage: **N-version programming / design diversity** (Avizienis) — reviewers chosen from different vendors so failure modes decorrelate |
| Parallel independent audit round | **Jury** |
| Sequential audit round, each auditor sees prior findings | **Relay** |
| Claude's findings locked in the exchange file before any auditor output is read | **Pre-registration** (blinded analysis) |
| `tmp/ai_exchange/` shared findings files — spoken of as **"the blackboard"** (synonym adopted 2026-08-20; the directory name stays) | **Blackboard architecture**. One deliberate deviation: during blind jury rounds the blackboard property is suspended — auditors must not read each other's postings until reconciliation |
| **Read-only diagnosis** sessions, probe before repair | **Hypothesis-driven debugging** |
| `devlog/` session logs | Engineering **daybook**; decision entries serve as **ADRs** (architecture decision records) |
| The evidence ladder's ranking of AI reviewer agreement as the **weakest** rung | **Correlated failure in multiversion programming** (Knight and Leveson 1986) — an experimental test of the failure-independence assumption behind N-version programming rejected it. Our ladder applies that result to LLM ensembles: vendor diversity decorrelates some failure modes, so agreement bounds error rather than establishing truth |
| A jury brief's **"established facts — do not re-litigate"** section carrying an unverified claim | No published analog in software engineering; the nearest is **unblinding / protocol contamination** in clinical trials. A fact placed where auditors are told not to challenge it propagates into every voice and is invisible in the reconciliation, so a round can return unanimous and wrong |
| **Blind slice assignment** — independent agents cluster disjoint slices of one corpus, and convergence on a category is read as evidence about the taxonomy | **Independent coding / inter-rater reliability** in qualitative research methodology. Our use is inferential rather than metric: convergence is treated as a signal that a category is real, not as a reliability coefficient |
| The **blindness trust boundary includes the tooling** | No published analog found. Blinding protocols assume the channel between parties is what needs controlling; here the environment breached the blind with no party acting |
| **Recurrence-triggered promotion** — a rule that fires repeatedly while still unenforced is escalated from written rule to mechanical check | No published analog found. Related in spirit to defect-density thresholds, but the trigger is the rule's failure to hold, not the defect count |

Note for readers from the computational-science literature: our **"verified"** means an
executable gate ran. It is narrower than *verification* in the V&V sense (solving the equations
right, as opposed to *validation*, solving the right equations) and narrower again than
verification in formal methods (proof against a specification). Where we mean either of those,
we say so explicitly.

## Deviations from common practice

These are deliberate, and they are the point:

1. **Pre-registered findings.** The primary AI writes its complete findings before any
   auditor is invoked, and never edits them afterward. Agreement is only meaningful when
   it cannot be produced by anchoring.
2. **Source-grounded verification.** Every auditor claim carrying a file:line citation is
   re-checked against the actual source before acceptance and labeled CONFIRMED /
   REFUTED / UNVERIFIABLE with the evidence quoted. Reviewer output is testimony, not
   verdict.
3. **Physics and literature outrank reviewer agreement.** Questions of physics,
   formulation, or design are never settled by vote among models — they are routed to the
   human maintainer and adjudicated against the published literature
   (`doc/literature_references.md`).
4. **The devlog is durable, not scaffolding.** Exchange files are disposable by design,
   but the distilled record of what was confirmed, refuted, or left open — with citations
   and confidence levels — is kept permanently in `devlog/`.
5. **Agreement is bounded evidence, and the protocol says so.** Reviewer diversity across
   vendors decorrelates some failure modes, not all. Knight and Leveson (1986) tested the
   independence assumption behind N-version programming directly and rejected it:
   independently developed versions fail together far more often than independence
   predicts. We treat that result as applying to LLM ensembles, which is why inter-AI
   agreement sits at the bottom of the evidence ladder, below a static source trace. In
   this project's record, unanimous reviewer agreement has been wrong where a single
   executed probe was right. Diversity bounds error; it does not establish truth.
6. **The blind includes the tooling.** Blinding is a property of the whole environment, not
   of participant discipline. This was learned directly: a notification hook injected a
   completing auditor's output into the primary AI's context before its findings were
   locked. The breach was disclosed in place, the affected agreement marked
   non-independent, and the remaining auditor dispatched after the lock so one voice stayed
   provably clean. Standing rule since: harness hooks and automation are audited as part of
   the review protocol, and a breach is recorded and worked around, never averaged away.
7. **Mechanical checks carry their own falsifiers.** A check that cannot fail produces false
   confidence at the cost of a real check. The project's own distillation run demonstrated
   both sides in one night: sampled citation audits held at 49 of 50, while a
   duplicate-detection check keyed on shared defect IDs reported a reassuring rate and
   missed every duplicate that cited none. Before any check's result is used, the question
   "what would this have reported if the answer were bad?" must have a concrete answer.
8. **Written rules graduate into mechanisms.** A documented rule that is violated anyway
   needs enforcement, not restatement. Each rule in `doc/lessons_learned.md` names an
   enforcement target at birth and carries a status (proposed, adopted, automated,
   execution-verified); a rule that fires three times without graduating forces a promotion
   review toward a test, an assertion, or a script. The firing signal is designed to come
   from closure fields in the debt register's protocol; those fields are proposed alongside
   the rules and adopt with them.

## The method, measured

The record this methodology produced was itself distilled in August 2026: 298 devlog entries
and the debt register yielded a 539-incident catalog and 18 operating rules
(`doc/lessons_learned.md`), extracted under the same protocol it documents: pre-registered,
blind-audited across vendors, with sampled citation verification.

Two calibration results are worth stating. First, when cross-referenced against the
scientific-software-engineering literature, four of the extracted rules independently reproduce
published results — among them the oracle problem, single-source-of-truth, and the
verification-before-validation ordering — which is evidence that the extraction machinery finds
real signal. Six have no analog in the sources consulted, though adjacent literature exists
for two of those and was not read. The remainder concern epistemic hygiene in long-running,
multi-model collaborative development, and that is the territory this document exists to
describe.

Second, the audits found the record's weaknesses to be systematic rather than random: sampled
citations survive verification at high rates, while counts, denominators and deduplication do
not survive without checks designed to fail. The reader adopting this protocol should expect
the same distribution.

The governing target throughout is not the error rate but the error survival rate: under this
protocol errors remained frequent, and their survival time approached zero. The distinction,
and the scaffolding that achieves it, is the substance of `doc/lessons_learned.md`.
