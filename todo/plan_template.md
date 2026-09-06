# Todo Plan Template

**Date:** 2026-07-01
**Purpose:** Canonical structure for substantial `./todo/` implementation plans, extracted from
`closed/meshfile_refactor_plan.md` (the `.bfm` save/load refactor — the reference example of a plan
that stayed navigable through a multi-week, tri-AI, 25-defect campaign). Use this template for
any plan expected to span multiple sessions, audits, or AIs. Small single-session bug notes do NOT
need every section, but they DO need the header, live checkboxes, and file:line citations.
**Module:** meta (documentation convention)

> **How to use:** copy the skeleton in §T below into `todo/<lowercase_with_underscores>.md`, delete
> sections that do not apply (keep the remaining sections in order), and register the new file in
> `todo/README.md`. When the plan completes, move it to `todo/closed/` and update the README entry
> with a DONE summary.

---

## Conventions (apply to every plan)

- **Live checkboxes** for every actionable item: `[ ]` open, `[x]` done, `[◐]` partially done /
  in progress. Tick the box **in the same turn** the fix lands in code. Never delete an obsolete item —
  strike it through (`~~text~~`) and say why.
- **Stable item IDs**, so audit threads and devlogs can reference them across sessions:
  - `Dn` — defects found during implementation/audit (D1, D2, …), each with severity
    (CRITICAL / HIGH / LOW), a root cause, and — once resolved — a **"Fixed YYYY-MM-DD"** line
    stating *what* changed and *who verified it*. Retracted findings stay in place, marked
    **FALSE POSITIVE (retracted YYYY-MM-DD)** with the reason — negative results prevent re-audits.
  - `Rn` — ordered plan steps (R1, R2, …), with dependencies noted as `(after: …)`.
  - `On` — open design questions (O1, O2, …), **logged, not silently decided**; when resolved,
    annotate **RESOLVED YYYY-MM-DD → decision** in place.
- **Every claim carries a code citation** (`cl_File.cpp:123-145`) or is explicitly marked as an
  assumption. Attach confidence (high / medium / low) to non-trivial claims per the collaboration
  protocol.
- **Historical decay is marked, not deleted.** When a section goes stale (code deleted, design
  superseded), prepend a `> **HISTORICAL (reconciled YYYY-MM-DD).**` blockquote saying what still
  holds and where the current truth lives. Dead citations stay but are flagged dead.
- **The Status line in the header is the single source of truth** for where the plan stands —
  update it whenever the state changes materially (audited / blocked / n of m steps done / DONE).
- **Attribution:** name who found, fixed, decided, and verified each item (Claude / Codex / Grok /
  Christian). Decisions by the user are recorded as such (`decided YYYY-MM-DD, Christian`).
- Audit threads live in `tmp/ai_exchange/<slug>.md` (ephemeral); distill them into the plan and
  devlog before they are swept, and list them in §8 so the trail can be reconstructed.

---

## §T — Skeleton

```markdown
# <Title: What Is Being Built/Fixed>

**Date:** YYYY-MM-DD
**Purpose:** One paragraph: the user-visible goal, and the mechanism in one sentence.
**Module:** `src/<module>` (+ secondary modules)
**AIs involved:** Claude (exploration + plan), Codex (audit), Grok (third voice) — as applicable
**Status:** OPEN | PLAN — audited, pending approval | IN PROGRESS (Rn of Rm) | ✅ COMPLETE (date)
    On completion, rewrite Status as a self-contained summary: what landed, what was verified
    (which tests / procs), and the residual follow-ups with pointers to their new todo files.

> **Scope guards (from the task brief):**
> - What is explicitly OUT of scope (e.g. "parallel distribution is out of scope").
> - Compatibility promises kept or explicitly dropped.
> - The entity/state list that IS in scope.

---

## 1. Current Behaviour and How It Fails

What the code does today, with citations. Classify the failure modes in a table:

| Failure | Mechanism | Evidence |
|---|---|---|
| Silent data loss — X | why it happens | `file.cpp:lines` |
| Hard error — Y | … | … |

End with a **"Bottom line:"** paragraph — the one-sentence diagnosis the rest of the plan rests on.

## 2. Architecture: Why <Approach> Is the Right Spine

The chosen mechanism and why (prefer reusing a battle-tested engine over inventing one).
Name the alternative(s) rejected and the deciding reason.

## 3. Gap Table

One row per piece of state/behaviour in scope. Classify each:
- **(a)** deterministically rebuildable (state *from what*),
- **(b)** ambiguous / open question (link the On),
- **(c)** must be handled explicitly.

| # | State | Needed for | Handled today? | Class | Citation / rationale |
|---|---|---|---|---|---|

Follow with **§3.1 Cross-cutting findings** — risks that span rows (e.g. "index vs ID keying"),
called out as the most important correctness changes.

## 4. Ordered Steps

Steps `R1…Rn`, each independently testable, with dependencies noted as `(after: …)`. Include an explicit
round-trip / end-to-end **test step** as the final gate.

### 4.0 Implementation Progress (updated YYYY-MM-DD)

Add this subsection once implementation starts; it becomes the living heart of the plan. It contains:
- the **"Implemented and audited"** checklist (what landed, with the key design facts a future
  reader needs — not just "done"),
- the **Dn defect tracker** (see Conventions), grouped by audit round with the round's date and
  which AIs participated,
- a **"Still unimplemented"** / "out of scope by design" paragraph so absence is distinguishable
  from omission.

## 5. Open Design Questions (not silently decided)

`O1…On`. Each states the question, the options, and — when resolved — the decision + date in place.

## 6. <Schema / Format / Interface Design>

The concrete artifact being designed (file format, API contract, data layout), precise enough to
audit against the code. When the artifact ships, move the authoritative reference into the module
doc (`src/<module>/doc/…`) and mark this section historical.

## 7. Definition-of-Done Checklist

- [ ] Every gap-table row mapped to a step or an open question.
- [ ] Each claimed gap backed by a citation, not assumption.
- [ ] Ordered steps with dependencies.
- [ ] Open questions logged, not decided.
- [ ] End-to-end test passes (name the reproducer and configurations).

## 8. Audit Trail

- Exchange thread(s): `tmp/ai_exchange/<slug>.md` (swept per protocol §10 after distillation).
- Per-AI findings summary: what each auditor caught, and the note that findings were
  independently re-verified against the cited code before inclusion.

## Appendix <X> — Decision: <topic>

One appendix per non-obvious decision worth defending later: the question, the answer, why the
tempting alternative breaks (enumerated reasons), precedent in the codebase, and a decision-summary
table. Mark **SUPERSEDED IN PART** with a pointer when later work revises it.
```

---

## Lifecycle

1. **Draft** (Claude) → **audit** (Codex, optionally Grok third voice) → status `PLAN — audited,
   pending user approval; no source modified`.
2. **Approval** (Christian) → implementation; §4.0 appears and the Dn tracker starts.
3. Each session: tick boxes, log new Dn/On in place, update the Status line, distill the exchange
   thread. Have Codex polish the prose after substantial edits (read-only pass).
4. **Done:** rewrite Status as the completion summary, split residual follow-ups into their own
   todo files (link them), move the file to `todo/closed/`, update `todo/README.md` (move the entry
   to the Closed section with the DONE summary), and record the session in `devlog/`.
