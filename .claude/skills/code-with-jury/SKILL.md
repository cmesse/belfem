---
name: code-with-jury
description: BELFEM's standing loop for any change to executable code beyond a trivial fix — plan, plan jury (Codex + Grok), implement, code jury, close out. Load it before writing a plan or touching a source file in src/, tests/ or nonfree/, including a one-line fix to a factory, kernel, solver or IWG, and whenever the user says "plan+audit->code+audit", "jury", "audit this", or asks for a fix, refactor or new feature in C++/Fortran/CMake. Also names what is exempt (comments, docs, registers, diagnostic probes) so the round is not spent where it buys nothing.
---

# plan+audit → code+audit

Christian's standing rule: **every** code change beyond a trivial fix runs the full loop, and
every audit launches **both** Codex and Grok. Three vendors exist to cover each other's blind
spots — a second blind Claude shares Claude's training and is not a substitute.

The round mechanics are frozen elsewhere and are not restated here:
`.claude/commands/cross-review.md` (one round, step by step) and
`doc/ai_collaboration_protocol.md` §9–§11 (depth table, evidence ladder, stop condition).
This skill is the **outer loop** those rounds sit inside.

## Is the round owed?

| Change | Round owed? |
|---|---|
| Anything that compiles into the binary — C++, Fortran, CMake logic | **yes**, even one line |
| Comment text, Doxygen, `doc/*.md`, `todo/`, `devlog/`, debt register | no — lands on approval. User-facing docs still take the Codex **language sweep** (CLAUDE.md) |
| Diagnostic probe: temporary print, counter, dump | no — land it, strip it when its question is answered. Piggyback it on a fix's round when both are in flight |
| Input deck, `.geo`, example data | no, unless it changes a key's contract — then also `doc/input_file_reference.md` **and** `doc/input_schema.yaml` |

A task that mixes text and code **splits**: land the text, route the code through the jury.

## The loop

Two stops. Everything else runs unattended.

### 1 — Plan

Write the plan as a `todo/<topic>.md` file per `todo/plan_template.md` (Dn/Rn/On step IDs, gap
table, live checkboxes, Status line), or inline for a small fix. It must name: the defect or
goal, the files and functions touched, the mechanism of the fix, what could break, and the
**gate** that will decide it (a reproducer, `make check-fast`, a probe — name it now, not later).

### 2 — Plan jury

Pick the slug `review_<topic>`, pre-register your own complete findings on the plan to
`tmp/ai_exchange/<slug>.md`, then dispatch on the plan file's path. Follow
`.claude/commands/cross-review.md` from its step 2 — pre-registration before dispatch is the
protocol, not a formality.

### 3 — ⛔ STOP: plan approved, edits approved

Report the reconciled plan findings and the P0/P1 list. Ask for approval of the plan **and** for
permission to edit source — the session is read-only by default. Do not implement before both.

### 4 — Implement

Apply the plan as approved. Keep the diff to the plan's scope; a discovery that widens it comes
back here as a plan amendment, not as extra hunks. Syntax-gate what you write
(`-fsyntax-only` with `-std=gnu++17` and the real warning flags — below the compile rung, say so),
and never call it verified on that basis.

### 5 — Code jury

**Shared checkout:** `git diff HEAD -- <your files> > tmp/ai_exchange/<slug>.diff` and dispatch
on that path. The default `git diff HEAD` would hand the auditors every peer's half-finished
edit as your own. Pre-register that the diff is own-files-only.

Same slug as the plan round, so the plan and code entries sit in one thread. Then the mandatory
verification pass: re-check **every** auditor `file:line` against source, label
CONFIRMED / REFUTED / UNVERIFIABLE with the evidence quoted, and append the reconciliation table.

`git status` after any Grok round — Grok has edited source during a read-only audit before.

### 6 — ⛔ STOP: fixes are a separate approved step

Report the reconciliation. Never apply a finding inside the round that raised it. Then hand
Christian the gate to run — **he runs builds and binaries**, not you.

### 7 — Close out (same session)

Devlog `devlog/dlYYYYMMDD_topic.md` + `devlog/README.md`. Tick the plan's checkboxes and refresh
its Status line. An item opened or closed also updates `devlog/campaigns/<name>.md` **and**
`todo/debt_register.md` in this session — register rows are exactly 6 cells, rebuilt cell by cell.
Work scoped to `nonfree/` logs to `nonfree/devlog/` and suppresses the exchange tier entirely:
no exchange file, no jury, coordination in chat.

## Depth

`doc/ai_collaboration_protocol.md` §9.1 owns the table and `.claude/commands/cross-review.md`
carries the operational excerpt — read one of those rather than a third copy. Set all four
variables explicitly; `cross_review.sh --jury` refuses without them, so a round's depth is always
in the record.

The row this loop reaches for most often: **critical framework code** — `src/mesh`,
`src/fem/kernel`, dof manager, sidesets, blocks, calculator, basis, distributor, `src/fem/maxwell`
factories — puts `gpt-6-astra` on the Codex seat for the plan audit and the code audit alike.
Name the model in the exchange entry either way.

## When a voice is missing

If either vendor's quota or budget runs out, **say so immediately** and let Christian decide — he
can top up credits. Do not silently proceed on two voices, and do not substitute a blind Claude.
Note the missing voice in the exchange file with the reduced-independence caveat.

## Hard limits

- **Cohomology core is closed to AI edits** — `cl_Cohomology`, `cl_Homology`,
  `cl_SimplicialComplex`, `cl_Chain`, `cl_Cochain`, `fn_Smith`. Session edit approval does not
  lift it. Read, report to Gregory.
- **Physics and design questions are never settled by vote.** Route them to Christian and the
  literature. Two auditors agreeing is the weakest rung of the evidence ladder and never lifts a
  verdict past source-trace.
- **Say "reviewed", not "verified"**, until an executable gate has actually run.
- **Stop condition:** no further round on the same subject without a material code change, new
  experimental evidence, reviewer disagreement on a load-bearing point, or a safety boundary.
  Otherwise the next step is the gate.
