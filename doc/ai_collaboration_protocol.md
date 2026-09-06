# BELFEM AI Collaboration Protocol {#doc_ai_collaboration_protocol}

**Date:** 2026-03-13
**Purpose:** Concrete roles, communication standards, quality gates, and file conventions for Claude/Codex cooperation in the BELFEM codebase.
**Status:** Official extension of `AGENTS.md` and `CLAUDE.md`. This file is the authoritative protocol.

---

## 0. Instruction Precedence

If guidance conflicts across collaboration docs, apply this order:

1. `doc/ai_collaboration_protocol.md`
2. `AGENTS.md`
3. `CLAUDE.md`

---

## 1. Roles

**Claude Code — Primary AI (broad exploration)**

- Navigates the full codebase, identifies architectural issues, proposes refactors
- Writes implementation plans, documentation, and code changes
- Handles literature lookups and cross-referencing
- Scope limit: flag MPI edge cases, thread-safety / hybrid MPI+OpenMP hazards, and C++ standard subtleties for Codex audit rather than assuming they are correct

**Codex — Secondary AI (precision audit)**

- Audits Claude's findings for C++17 standard compliance, MPI/OpenMP correctness, and edge-case logic
- Double-checks critical claims: memory safety, undefined behavior, standard library guarantees
- Validates low-level correctness that broad exploration may miss
- **Must prefer independent verification** over simple agreement and ground audits in `doc/coding_philosophy.md`
- Scope limit: does not restructure architecture or make broad design decisions unilaterally

---

## 2. Audience Tiers & the Communication Channel

### Audience Tiers

Collaboration artifacts are organized by their **reader**, and lifetime follows audience:

- **AI-only artifacts** — read only by AIs (the exchange). Ephemeral, machine-parseable, disposable. No human ever navigates them, so they can be garbage-collected aggressively.
- **AI+human artifacts** — read by both AIs and humans (`./todo/` progress files, `./devlog/` session logs). Durable, tracked in-repo, curated for a human skimming for decisions.

An artifact is ephemeral **iff no human ever needs to return to it**. The boundary between the two tiers is a **distillation** step, not a copy: at thread/session close the signal is *lifted* from the AI-only scratch into the durable AI+human record (the devlog, and `./todo/` progress files as applicable) before the scratch is swept. That distillation is precisely what makes aggressive GC of the AI-only tier safe.

**Format follows audience.** AI-only files stay machine-parseable — keep the `# AI`, `## Audit`, and `Confidence:` markers so the receiving model can orient — but need no human-prose polish. AI+human files are curated and distilled for a human reader.

### Communication Channel — `./tmp/ai_exchange/<slug>.md` (AI-only, ephemeral)

The AI-to-AI exchange is a set of **sharded, ephemeral, per-task files** under `./tmp/ai_exchange/`: one file per topic, `./tmp/ai_exchange/<slug>.md`, where `<slug>` is a short `lowercase_underscore` topic tag (e.g. `periodic_thin_cut`). **One topic per file.** Either AI may create, read, or append to the file for the current task at any time.

**"The blackboard" is a spoken synonym for this channel** (adopted 2026-08-20, after the published term it maps to — see `METHODOLOGY.md`). An instruction like "see blackboard" or "put it on the blackboard" means the current task's `./tmp/ai_exchange/<slug>.md` file. The directory name stays `ai_exchange`; only the vocabulary gained the alias. One deliberate deviation from the published pattern: during **blind jury rounds the blackboard property is suspended** — auditors must not read each other's postings until reconciliation.

These files are AI-only scratch: machine-parseable, disposable, and **never committed** (`./tmp/ai_exchange/` is git-ignored). The durable record of any conclusion lives in the devlog (§6), not here. Do not archive them — they are swept, not retained (see §10).

### Nonfree Module Exception

For work scoped to `./nonfree/` or to porting code from `./tmp/manta/` into
`./nonfree/`, do **not** create, read, or append the AI-only exchange at all —
suppress the AI-only tier entirely. The shared exchange is routed through
external AI vendors and is therefore a vendor-bound leak surface; nonfree work is
Christian's proprietary project and must not cross it. Keep coordination in chat
and record session outcomes only in the local human-readable record
(`./nonfree/devlog/`).

### Entry Format

Copy this format exactly for every entry:

```markdown
# CLAUDE 2026-03-13 14:20:00 PDT
## Query / Finding
Confidence: medium (~65%)

Body of the message. Reference specific files and line numbers.
Invite the partner AI to challenge the claim.

---

# CODEX 2026-03-13 14:25:30 PDT  (model=gpt-5.6-terra, effort=high)
## Audit & Verdict
Confidence: high

Codex's analysis. Checklist items addressed (with N/A noted).
Counter-points or confirmation with independent evidence.

---

# CLAUDE 2026-03-13 14:27:15 PDT
## Resolution / Next Action

Summary of outcome. What to do next.

---
```

**Rules:**
- Always start with `# AI_NAME YYYY-MM-DD HH:MM:SS TZ`
- Auditor entries additionally carry the depth that produced them, `(model=…, effort=…)`, appended
  by the wrapper — never written by hand. A knob the caller did not choose is marked `[defaulted]`.
  Claude's own entries carry no stamp; the asymmetry is deliberate, since the point is to attribute
  an auditor's *findings* to the depth that produced them (see §9.1)
- One topic per file (`./tmp/ai_exchange/<slug>.md`); start a **new file** for a new topic
- No archiving and no size threshold — these files are ephemeral. Distil the conclusion into the devlog (§6) before the file becomes GC-eligible (see §10).

---

## 3. Communication Style — Calibrated Uncertainty

AIs are typically trained to sound confident. In this project, we do the opposite: **communicate uncertainty honestly** so the partner AI and the user can assess claims independently.

### Three-Tier Confidence

| Tier | Meaning | Action for receiving AI |
|------|---------|------------------------|
| **high** | Strong evidence, verified in code and/or literature | Spot-check is sufficient |
| **medium** (~N%) | Likely correct but not fully verified; add approximate % when helpful | Verify before acting on it |
| **low** | Educated guess, needs investigation | Treat as hypothesis, investigate independently |

### Examples

```
The race condition is in cl_Comm_Table.cpp:245 where the buffer is reused
before MPI_Wait completes. [confidence: high — confirmed by reading MPI-3.1 §3.7.3]

I think the segfault is triggered by the Jacobian assembly, but it could also
be the material lookup. [confidence: low — stack trace is ambiguous]

The thin-shell N=1 case should converge for single tapes but will miss
top/bottom losses for stacks. [confidence: medium (~75%) — matches Alves et al.
2022b Table 1, but I haven't verified the current code path]
```

**Why this matters:** When Claude says "I'm ~70% sure the bug is here," Codex can focus its audit energy on verifying that specific claim rather than re-scanning everything. When Codex finds a C++ standard violation, it should say whether it's guaranteed UB or merely implementation-defined.

---

## 4. Codex Audit Checklist

Before signing off on a Claude finding, Codex addresses **only the relevant items** from this checklist. Items that don't apply to the change under review should be marked "N/A" — there is no requirement to comment on all 8 for a localized fix.

When a deviation from `doc/coding_philosophy.md` is found, Codex flags it with a confidence level and proposes the corrected pattern.

| # | Item | What to check |
|---|------|---------------|
| 1 | **Naming conventions** | `a`/`t`/`m`/`g` prefixes in framework code; mathematical-kernel exception for standard notation |
| 2 | **Container selection** | `Cell<T>` for arrays, `Vector<T>`/`Matrix<T>` for linear algebra only, `DynamicBitset`, `ShiftRegister` |
| 3 | **Memory management** | Manual `malloc`/`free` on hot paths, ownership via `Cell<T*>` (non-owning) vs `Cell<T>` (owning), alignment |
| 4 | **Error handling** | `BELFEM_ASSERT` for debug-only checks, `BELFEM_ERROR` for always-active runtime errors, release-build behavior |
| 5 | **MPI patterns** | 64 KB chunking, `data()` safety on empty `Cell`, alignment, distributed operation correctness |
| 6 | **Thread safety / OpenMP** | BELFEM is MPI-first and deliberately not thread-safe internally. Flag any accidental internal thread assumptions or hybrid-model hazards. OpenMP is possible but requires external synchronization (`#pragma omp critical`) |
| 7 | **Performance-first rules** | No hidden allocations in loops, zero-abstraction penalty, preallocate buffers |
| 8 | **Literature compliance** | When applicable: did Claude follow the routing tables? Do citations match the implementation? (Use `CLAUDE.md` citation format) |

---

## 5. Literature-First Rule (Conditional)

Literature routing is **mandatory** when the task involves:
- Algorithms or formulations (h-φ, thin-shell, cohomology, MPFA, etc.)
- Numerical methods or solver strategies
- Correctness validation or debugging of physical behavior
- Documentation of *why* something is implemented a certain way

Literature routing is **not required** for:
- Localized code style fixes (typos, formatting, naming)
- Build system changes
- Null checks, memory leak fixes, or other mechanical patches
- Documentation edits that don't touch algorithmic content

**When literature is used:**
1. Claude consults the routing tables (`literature/papers/fem/index.md`, `literature/books/index.md`, or `doc/literature_references.md` as fallback)
2. Codex verifies that Claude's citations are correct and that the implementation matches the cited equations/sections
3. Both AIs use the citation format defined in `CLAUDE.md`: "Author et al. YYYY, Section X" or "Bathe, §X.Y". The old `paperN` aliases are retired — do not write new ones; `doc/literature_references.md` decodes the ones in older records.

---

## 6. Session Documentation — Devlogs in `./devlog/`

At the end of every meaningful session, create a devlog summary in `./devlog/` and update `./devlog/README.md` with a one-line entry linking the new file.

**The devlog is the durable AI+human distillation of the ephemeral exchange (§2).** The per-task `./tmp/ai_exchange/<slug>.md` files are AI-only scratch and are swept. At thread/session close, **lift the conclusion** — what was confirmed, refuted, or left open, with file:line evidence and confidence — into the devlog (and into the relevant `./todo/` progress file when there is one) *before* the exchange file becomes GC-eligible. If it is not in the devlog or a `./todo/` file, treat it as lost once the scratch is swept.

**Distinction from `./todo/`:** Devlogs record *what was changed and why* (backward-looking). Task files in `./todo/` document *what needs to be done* (forward-looking). Both are AI+human, durable, and tracked in-repo — in contrast to the AI-only exchange.

**Todo-file prose (standing rule, revised 2026-08-20):** `./todo/` files get **no** Codex prose sweep. The earlier rule required one on every draft or substantial edit; it is withdrawn as not worth its cost. `./todo/` files are working artifacts for the AI-plus-human team, not documentation a user reads, and the sweep spends vendor budget to smooth prose whose density is more useful than its polish. This changes nothing about **technical** audits of a todo file's content, which are unaffected and still expected where the plan warrants one. The prose sweep is now reserved for user-facing documentation — see `CLAUDE.md` §"Prose Gets a Language Sweep" for the scope and the exclusions.

**Todo-file checkboxes (standing rule):** whenever a `./todo/` file lists steps, **always give each step a GitHub checkbox** (`- [ ]`). Keep the boxes live:
- Tick a box (`- [x]`) the moment its step is completed — the file is the running progress record, not just the initial plan.
- When a step becomes obsolete, **strike it through** (`~~step text~~`) rather than deleting it, so the history of what was abandoned and why stays visible.
See `todo/closed/periodic_thin_cut_continuity_fix.md` for the reference style (checked, unchecked, and struck-through steps), and `todo/plan_template.md` for the current skeleton.

### Nonfree Module Exception

For work scoped to `./nonfree/` or to porting code from `./tmp/manta/` into
`./nonfree/`, write session devlogs to `./nonfree/devlog/` instead of
`./devlog/`. Do not update `./devlog/README.md` for nonfree-only sessions.
If `./nonfree/devlog/` has its own README or index, update that local index;
otherwise the dated devlog file is sufficient.

### Filename Format

`dlYYYYMMDD_topic.md`

Examples:
- `dl20260313_solver_refactoring.md`
- `dl20260315_mpi_deadlock_investigation.md`
- `dl20260320_thin_shell_convergence_fix.md`

### Required Header

Fill only the applicable fields:

```markdown
# Devlog 2026-03-13 — Solver Refactoring

**Date:** 2026-03-13
**Topic:** Brief description
**AIs involved:** Claude, Codex
**Claude Confidence:** medium (~70%)
**Codex Audit Confidence:** high
**Literature References:** Messe et al. 2023 §2.7, Bathe 2016 Ch. 8
**Verification:** focused regression — `make check-fast` green on sideconnectors @ <commit> (see §11; "reviewed" ≠ "verified", name the highest evidence level reached; include branch/commit + command for probe level or above)

## Summary

Brief description of findings and outcomes.

## Key Findings

- Finding 1 (with file:line references)
- Finding 2

## Changes Made / Proposed

- List of files modified (if any), with brief rationale

## Open Questions

- Anything unresolved, for future sessions

## Files Updated

- path/to/file1.cpp
- path/to/file2.hpp
```

### When to Write a Devlog

- End of any debugging or investigation session
- After completing a refactoring task
- When resolving a significant exchange thread from `./tmp/ai_exchange/<slug>.md` (distil the conclusion into the devlog before the scratch is swept)

---

## 7. Edit Safety Rule

- **Investigation and review = read-only by default.** The working name for this mode is a
  **read-only diagnosis** (published analog: hypothesis-driven debugging — see
  `METHODOLOGY.md`): probe before repair, and no source edit until the diagnosis is accepted.
- Writing to `./tmp/ai_exchange/` (AI-only exchange), `./todo/`, and `./devlog/` is always allowed (AI exchange, task plans, devlogs).
- For nonfree-only sessions, use `./nonfree/devlog/` for devlogs and suppress the
  AI-only exchange (`./tmp/ai_exchange/`) entirely unless explicitly requested.
- **Source-code edits only after the user explicitly says editing is approved.**

### 7.1 Protected Module — the Cohomology Core

**Hard policy, 2026-08-31 (Christian and Gregory Giard, who owns `src/homology`). No AI —
Claude, Codex, Grok, or any successor — edits the cohomology core. This is not the ordinary
read-only default: it is not lifted by the user approving edits for the session.**

Under the ban, in `src/homology/`:

| unit | files |
|---|---|
| cohomology | `cl_Cohomology.cpp` / `.hpp` |
| homology | `cl_Homology.cpp` / `.hpp` |
| complex reduction | `cl_SimplicialComplex.cpp` / `.hpp` |
| chains | `cl_Chain.cpp` / `.hpp` |
| cochains | `cl_Cochain.cpp` / `.hpp` |
| Smith normal form | `fn_Smith.cpp` / `.hpp` |

**Not under the ban** — ordinary edit-safety rules apply: the rest of `src/homology/`
(`cl_CutFactory`, `cl_CutProcessor`, `cl_CutProcessorManual`, `cl_CutData`, `cl_CutSet`,
`cl_InterfaceProcessor`, `cl_BeltedTree`, `cl_Topology`, `en_CutAlgorithm`, `CMakeLists.txt`),
the module's `doc/`, and every caller outside the directory.

**Why — the algorithm is novel, and that is the whole of it.** This is not a vague claim that
the code is hard. What runs here is a *modified* Pellikka reduction — Gregory's own work, and
as of this writing still unpublished. It was drafted in 2025, shelved, and picked up again in
2026, so citations in this tree carry both years for one work — see
`doc/literature_references.md`, "Giard, G., et al., *Generalized Pellikka algorithm for cohomology
computation*, in preparation". A draft in preparation is *stronger* ground for this policy than a
publication would be: there is no published implementation for a model to have learned from.

That inverts the usual situation. Meeting BDF5, RCM, or a CG loop, a model has read thousands of
implementations and carries a strong and largely correct prior about what the code should look
like — a genuine deviation stands out as one, and the prior is an asset. Here the nearest object
in the weights is *textbook* Pellikka. Every deliberate modification therefore reads as a
defect, and "fixing" it regresses the algorithm toward the published version it deliberately is
not. The model cannot tell "this differs from the algorithm I know" from "this is wrong",
because the algorithm it knows is the wrong reference.

**The canonical case, and it is attested rather than imagined.** `coreduceOmit()`
(`cl_SimplicialComplex.cpp:1236-1257`) runs `pCoreduce(0..2)`, then loops: while any 0-cosimplex
remains, it removes *one* 0-cochain, zeroes its coboundary's references, and re-runs the whole
`pCoreduce(0..2)` cascade. One node per stall, each removal unblocking a further coreduction.
It is BELFEM's own device for the problem Mrozek & Batko 2009 §5 identifies — a simplicial
complex admits no elementary coreduction pair at all, because a boundary is never of size one —
though *not* Mrozek's solution to it: the paper adds ∅ as a simplex of dimension −1 with every
vertex in its coboundary (one per connected component), and licenses no 0-cell omission
(`mrozek2009.txt:444-464`). Do not conflate the two; an earlier draft of this section did. `reduceOmit()` (`:939-961`) is the strict dual on the chain side and removes
**top-dimensional** chains instead, which is what Gregory means by "for cohomology, and not for
homology, we remove a node". Nothing at either site says any of this.

What an agent does with that loop is on the record. On 2026-07-03 the proposal was to batch the
omission — it removes every 0-cochain in the end, so why re-run the cascade each time? It was
struck as unsafe: batch omission strands 1-cells whose boundary is empty, and `pCoreduce` only
ever pairs a cell whose boundary holds exactly one entry (`cl_SimplicialComplex.hpp:364`), so the
junk survives into cocombine as spurious generators
(`devlog/dl20260703_coreduce_risk_reeval.md:25-27` — that entry's own Mrozek attribution is
wrong, see above; its mechanical argument is not). The loop looks trivially hoistable. It is
not. That is the failure mode in one line: code that compiles, runs, and is wrong.

The test that generalizes is **novelty relative to training data, not difficulty.** A hard but
widely published algorithm is safe ground; a simple but unpublished one is not.

**Do not treat concurrence as coverage.** Reviewer agreement is worth least exactly here,
because the errors are *correlated*: models missing the same prior fail the same way. They will
not always agree — on 2026-08-30 Claude and Grok both read a deliberate invariant as a missing
merge projection and Codex refuted it with a Schur derivation
(`devlog/dl20260830_dr29_coreduce_analysis.md:41-46`) — but a majority verdict on this module
carries no more weight than a single voice, and two-thirds of one jury was wrong.

Note what this policy is **not**. It is not a quarantine of a broken module — it has not been a
source of defects for a long time. Gregory's position on the non-unit-coefficient behavior is
that the obstruction is mathematical (a coarse mesh admits no unit representative), not a bug;
that is separate from the implementation's handling of it, which
`src/homology/doc/thin_cut_nonunit_rectification.md:22-30` documents as a real limitation
("the current code cannot handle this"). The ban protects working code whose correctness
conditions live outside it.

**Scope caveats, both open and neither settled by this policy.** `fn_Smith` is textbook Smith
normal form — the novelty test alone would not ban it; it is in the list on coupling grounds,
because the SNF stage consumes the reduced complex directly. `cl_BeltedTree` is a cohomology
algorithm (`cl_BeltedTree.hpp:26-46` holds `Cell<Cochain*>` generators and a
`SimplicialComplex*`) with no literature citation, and it is *outside* the list. By the stated
criterion the boundary has one member it does not justify and lacks one it does. Neither is a
production concern today (the default path is `PellikkaGeneralized`), and both are Gregory's and
Christian's to rule on. **Comment-only edits are likewise unruled** — treat them as banned until
someone decides, and propose the comment in a devlog instead.

**Still allowed, and still wanted:** reading, tracing, running, and reporting. Findings go into
`./devlog/`, a debt-register row, or a `./todo/` file, addressed to Gregory. What is banned is
touching the source and self-starting a fix — a fix plan in a register row is a proposal for
the module owner, never a queued task.

**Lifting it for one session** requires Gregory's named authorization, relayed explicitly by
Christian (e.g. "Gregory approved *this specific change* on *date*"), and a devlog line
recording who authorized what. Christian's ordinary edit approval does not reach these files.


---

## 8. Quick Start for Every Session

1. Read `AGENTS.md` → this file (`doc/ai_collaboration_protocol.md`).
2. Read the relevant module README and applicable `CLAUDE.md` sections.
3. Identify the question type → decide literature-first or code-first (see Section 5).
4. Open or append to the AI-only exchange for the current task,
   `./tmp/ai_exchange/<slug>.md`, except for nonfree-only sessions where this
   channel is suppressed unless explicitly requested.
5. State your confidence and wait for the partner AI's audit when the exchange
   channel is in use.

---

## 9. Claude-Initiated Codex Audit

Claude can invoke Codex directly — without waiting for the user to relay a prompt — when an independent check is warranted. This is the primary mechanism for covering blind spots that arise from training on different data and different company incentives.

### When Claude should call Codex

Call Codex proactively when:

- A claim about **C++17 standard compliance, MPI correctness, or UB** could go either way and the consequences are non-trivial
- A **memory-safety or ownership** invariant is asserted but not proven by reading the code
- The **cohomology, thin-shell, or periodic-BC logic** is intricate enough that a second read-through adds real value
- An implementation was **derived from literature** and independent verification of the translation is useful
- Claude's own confidence is **medium or below** on a safety-critical path

Do **not** call Codex for trivial fixes, naming corrections, build-system changes, or anything where Claude is already high-confidence and the risk of error is low.

### How to invoke

Choose the depth first — see §9.1. These recipes carry it; an invocation without it still runs,
but is recorded as `[defaulted]`.

```bash
CODEX_MODEL=gpt-5.6-terra CODEX_EFFORT=medium .claude/scripts/ask_codex.sh "Your audit prompt here"
```

Or pipe a longer prompt:

```bash
echo "Audit cl_CutProcessor.cpp:420-486 for sign-blind admission correctness..." \
  | CODEX_MODEL=gpt-5.6-terra CODEX_EFFORT=medium .claude/scripts/ask_codex.sh -
```

The script:
1. Resolves the per-task AI-only exchange path `./tmp/ai_exchange/<slug>.md`
   (slug from `$AI_EXCHANGE_SLUG`, else a sanitized session tag `sess_<id>` from
   `$CLAUDE_CODE_SESSION_ID`, else `scratch`), creating the directory/file if absent,
   and instructs the auditor to read the existing thread in that file first
2. Prepends the role preamble (Codex reads `AGENTS.md` and this protocol)
3. Resolves and validates the depth (`$CODEX_MODEL`, `$CODEX_EFFORT`; see §9.1), rejecting an
   unknown value before any prompt is read, so a typo costs no billed call
4. Runs `codex exec --sandbox read-only -m <model> -c model_reasoning_effort=<effort>` with
   working root at the BELFEM repo
5. Captures Codex's final response via `--output-last-message`
6. Appends a `# CODEX <timestamp>  (model=…, effort=…)` entry to the resolved
   `./tmp/ai_exchange/<slug>.md`
7. Echoes the response to stdout so Claude sees it inline

To pin a topic file explicitly, set `AI_EXCHANGE_SLUG=<topic>` (e.g.
`AI_EXCHANGE_SLUG=periodic_thin_cut`) before invoking the wrapper. The same
resolution applies to `ask_grok.sh`.

### 9.1 Depth Selection

Both wrappers take a model and a reasoning effort, and stamp both into the exchange entry. Choose
them before dispatch; do not let them default. The defaults exist so a bare invocation still works,
not as a tier — they were chosen to reproduce what the vendor configs gave on the day the wrappers
pinned them, and an unchosen knob is marked `[defaulted]` in the record precisely so it is not
mistaken for a choice. That correspondence is an observation of a local config, not something this
repository can check; if the vendor config drifts, the wrapper defaults stay put.

| Subject | Codex | Grok |
|---|---|---|
| Prose sweep of an ordinary guide or README; citation and doc-claim mechanics | `gpt-5.6-luna`, `medium` | not used |
| Prose sweep of a dense technical document (the input reference, the coding philosophy) | `gpt-5.6-terra`, `medium` | not used |
| Single narrow claim, round 1 | `gpt-5.6-terra`, `medium` | `high` |
| Plan audit or code-diff audit, round 1 | `gpt-5.6-terra`, `high` | `high` |
| Round ≥ 2, a round-1 split verdict, or a safety-boundary subject (ownership, lifetime, MPI collectives, ABI) | `gpt-5.6-terra`, `xhigh` | `xhigh` |
| Unattended post-commit `--quick` | `gpt-5.6-luna`, `medium` | `grok-4.6`, `medium` |

```bash
CODEX_MODEL=gpt-5.6-terra CODEX_EFFORT=high AI_EXCHANGE_SLUG=<slug> .claude/scripts/ask_codex.sh - < prompt.md
GROK_MODEL=grok-4.6 GROK_EFFORT=high        AI_EXCHANGE_SLUG=<slug> .claude/scripts/ask_grok.sh  - < prompt.md
```

Four things the table encodes, each of which cost something to learn:

- **The key is the subject's scope and the round number, not the confidence of the claim under
  audit.** The claim's confidence is the thing being tested; ordering a cheap round because you feel
  sure is exactly backwards when you are confidently wrong.
- **The second rung raises effort, not model.** Moving both knobs at once makes it impossible to
  learn which one helped. `gpt-5.6-sol` is allowlisted as an escape hatch but is named by no row,
  pending a measured comparison on a real BELFEM diff.
- **Effort is capped at `xhigh`.** This is a cost-and-conservatism policy, not a proven
  constraint: `max` has shown no measured benefit here, and `ultra` is documented by the vendor as
  delegating to subagents — a different execution shape than the single-agent read-only audit the
  wrappers assume, and one no in-repo evidence confirms. Raise the cap only with a measurement.
- **Effort is not a substitute for a sharp prompt.** The narration-stub failures documented in
  `doc/lessons_learned_evidence.md` were a prompt-shape defect; no effort setting would have
  touched them.

`scripts/cross_review.sh --jury` and `--relay` refuse to run unless all four variables are set, and
name this table when they do. `--quick` pins the cheap row itself and forces it onto the auditors,
so the unattended post-commit hook cannot be retuned by whatever shell happened to commit.

### Workflow

1. **Write a Claude query entry** to the per-task `./tmp/ai_exchange/<slug>.md` first (standard `# CLAUDE` header, Section 2 format), stating your claim and confidence.
2. **Call the script** with a focused audit prompt referencing the specific files/lines, and with
   the depth from §9.1 set explicitly.
3. **Read the response** from stdout (or from the delta hook on the next prompt).
4. **Write a resolution entry** (`# CLAUDE … ## Resolution`) summarizing what was confirmed, refuted, or left open.

### Prompt quality

A good audit prompt is specific:

> "Audit `src/homology/cl_CutProcessor.cpp:420-426`: does positive-only `tCase` admission correctly implement single-imposition for all interior shared faces, or can two adjacent elements both get positive cases for the same face? Reference the self-cancel block at lines 476-486."

A poor prompt is vague:

> "Check if the cut logic is correct."

---

## 10. File Summary

| Path | Audience | Purpose | Lifetime |
|------|----------|---------|----------|
| `./tmp/ai_exchange/<slug>.md` | AI-only | Live, per-task AI-to-AI exchange (sharded by topic) | Ephemeral; safe to delete anytime — recommended sweep of files older than 14 days |
| `./todo/ai_exchange.md` (legacy) | AI+human | Pre-2026-06-18 single-file exchange — historical record only, not written anymore | Persistent (frozen) |
| `./todo/ai_exchange_archive_*.md` | AI+human | Archived legacy exchange threads | Persistent |
| `./todo/*.md` (other) | AI+human | Task planning and implementation plans | Persistent |
| `./devlog/dlYYYYMMDD_topic.md` | AI+human | Session summaries — the durable distillation of the exchange | Persistent |
| `./devlog/campaigns/<name>.md` | AI+human | Current-state page per ACTIVE campaign (§11) | Persistent, rewritten in place |
| `./todo/debt_register.md` | AI+human | **Live** open-item register with the blocking-1.0 lens (§11); its preamble is the operating manual for both files | Persistent, one table |
| `./todo/debt_register_closed.md` | AI+human | Archive of retired (struck) register rows, same columns | Persistent, append-on-strike |
| `./nonfree/devlog/dlYYYYMMDD_topic.md` | AI+human | Session summaries for nonfree-only work | Persistent |
| `./doc/ai_collaboration_protocol.md` | AI+human | This file — the authoritative protocol | Persistent |

---

## 11. Evidence Hierarchy

### "Reviewed" is not "verified"

A finding or change is *reviewed* when a static audit is complete. It is *verified*
only when an executable gate has passed. Devlogs, todo trackers, and exchange entries
must not write "verified" for syntax-only or read-only work — write "reviewed", or
name the actual gate that ran.

### Evidence ladder

Strongest first:

1. **End-to-end reproducer** — a full model run reproduces (or stops reproducing) the behavior
2. **Focused regression** — a permanent test locks the claim (e.g. the interface/orientation battery via `make check-fast`)
3. **Compile / link** — the affected configuration builds
4. **Numeric probe** — a one-off computation confirms a value or identity
5. **Static source trace** — file:line reading of the code
6. **Literature consistency** — the claim matches the cited derivation
7. **AI reviewer agreement** — concurring independent audits

Lower levels support, never replace, higher ones. Reviewer agreement is the weakest
tier: three concurring audits do not lift a claim past level 5. Physics and design
questions are adjudicated by Christian against the literature regardless of level
(see the cross-review rules).

### Review stop condition

A new audit round on the same subject requires at least one of:

- a **material code change** since the last round,
- **new experimental evidence** (a failing or passing run, a probe result),
- **reviewer disagreement** on a load-bearing point,
- a **safety boundary**: ownership, parallelism, persistence, ABI, formulation.

Otherwise the next step is the executable gate, not another review.

### Campaign pages and the debt register (same-session rule)

`devlog/campaigns/<name>.md` (one page per ACTIVE campaign: current accepted design,
branch, last passing reproducer, open P0/P1, superseded approaches, links) and
`todo/debt_register.md` (the live rows, one table: `ID | area | severity | status |
reproducer | blocking-1.0?`, with `todo/debt_register_closed.md` carrying the retired
rows in the same columns) are the compression layer over the append-only devlog. A
session that resolves or creates an open item updates the campaign page and the register
**in the same session**; the dated devlog entry links to them and does not restate their
content. The register's own preamble governs how a row is read, when it may be struck,
and what counts as closure evidence — read it before editing a row.
