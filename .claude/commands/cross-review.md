---
description: Three-AI review round — Claude pre-registers, Codex + Grok audit headless, then verification + reconciliation
argument-hint: "[--jury|--relay] [path]"
---

Run one round of the frozen three-AI cross review on: $ARGUMENTS

Target: no path → the working-tree diff (`git diff HEAD`; clean tree → the last commit). A path argument → review that file (plan, spec, header) instead. Mode: `--jury` (default) = parallel, blind, independent audits; `--relay` = sequential, each auditor sees the thread so far. Never mix the two in one round.

Follow these steps IN ORDER — the ordering is the protocol (`doc/ai_collaboration_protocol.md` vocabulary applies):

1. **Resolve the target** and pick a topic slug `review_<topic>` (lowercase_underscore). The exchange file is `tmp/ai_exchange/<slug>.md`.

2. **Pre-register.** Review the target yourself, completely, BEFORE any auditor is invoked. Write your full findings to the exchange file as a `# CLAUDE <timestamp>` entry in the house format: every finding with file:line citations and per-finding confidence (high / medium / low). This entry is FROZEN once written — never edit it afterwards; all later material is appended below it.

3. **Choose the depth, then dispatch.** `cross_review.sh --jury/--relay` refuses to run without all four depth variables — a round whose depth is not recorded cannot be compared to any other round. Pick the row that matches the subject from the depth-selection table in `doc/ai_collaboration_protocol.md`:

   | Subject | Codex | Grok |
   |---|---|---|
   | single narrow claim, round 1 | `gpt-5.6-terra` / `medium` | `grok-4.6` / `high` |
   | plan or code-diff audit, round 1 | `gpt-5.6-terra` / `high` | `grok-4.6` / `high` |
   | round ≥ 2, a split verdict, or a safety-boundary subject (ownership, lifetime, MPI collectives, ABI) | `gpt-5.6-terra` / `xhigh` | `grok-4.6` / `xhigh` |

   Only these three rows apply here. A jury round always runs **both** auditors, so the
   Codex-only prose-sweep rows of §9.1 have no meaning for this command — reach for
   `ask_codex.sh` directly when you want a one-voice prose pass.

   Run as its own single background Bash command (never nested in a compound/piped command) and wait for it to finish:
   `CODEX_MODEL=gpt-5.6-terra CODEX_EFFORT=high GROK_MODEL=grok-4.6 GROK_EFFORT=high AI_EXCHANGE_SLUG=<slug> scripts/cross_review.sh --jury [path]` (or `--relay`).
   The auditors run headless and read-only via the existing `ask_codex.sh` / `ask_grok.sh` wrappers and append their entries to the exchange file themselves, each stamped with the model and effort that produced it.

4. **Verification pass.** Read the auditor entries from the exchange file. Re-check EVERY file:line citation against the actual source. Label each finding CONFIRMED / REFUTED / UNVERIFIABLE, quoting the evidence line under each label. Append this as a new `# CLAUDE <timestamp>` entry with a `## Verification` section.

5. **Reconciliation table.** Append at the end of the exchange file:
   `| finding | raised by | verdict | severity (P0–P2) | evidence | agreement |`
   The **evidence** column names the highest evidence-ladder level backing the verdict (`doc/ai_collaboration_protocol.md` §11: reproducer > regression > compile/link > probe > source trace > literature > AI agreement). Rules: two-plus independent raisers = high confidence, but agreement is the WEAKEST evidence tier and never upgrades a verdict past source-trace level. Single-raiser findings are flagged "needs human adjudication" — never auto-resolved. Physics and design questions are NEVER settled by vote — route them to Christian explicitly.

6. **Report to the user:** a short summary plus the P0/P1 list only; the full record stays in the exchange file. NEVER apply fixes as part of this command — fixes are a separate, explicitly approved step.
   **Stop condition (protocol §11):** do not propose another review round unless there is a material code change, new experimental evidence, reviewer disagreement on a load-bearing point, or a safety boundary (ownership, parallelism, persistence, ABI, formulation). Otherwise the next step is the executable gate (`make check-fast`, a probe, or a reproducer), not another review.
