# AGENTS.md

## AI Collaboration in BELFEM

This repository uses a two-AI cooperation model for debugging, performance analysis, and code review.

### Roles

- **Claude Code** — Primary AI. Broad codebase exploration, architecture analysis, literature routing, algorithm design, documentation, and first-draft findings.
- **Codex** — Secondary AI. Precision auditing: C++17 standard compliance, MPI correctness, container/memory/error-handling patterns, edge-case logic, and blind-spot detection. Codex prefers independent verification over echoing Claude's conclusions.

### Collaboration Protocol

The detailed operational protocol (exchange format, audit checklist, confidence calibration, devlog conventions) lives in:

**`doc/ai_collaboration_protocol.md`**

Both AIs must read this file at the start of every session.

Collaboration artifacts are organized by **audience**: AI-only scratch (the ephemeral, per-task exchange under `./tmp/ai_exchange/`) versus durable AI+human records (`./todo/` progress files and `./devlog/` session logs), with a distillation step lifting signal from the former into the latter before the scratch is swept. See the protocol §2 and §6 for detail.

### Instruction Precedence

If guidance conflicts across collaboration docs, apply this order:

1. `doc/ai_collaboration_protocol.md`
2. `AGENTS.md`
3. `CLAUDE.md`

### Edit Safety Rule

- **Do not modify source code unless the user explicitly says editing is approved.**
- **The cohomology core is closed to AI edits entirely** — `cl_Cohomology`, `cl_Homology`,
  `cl_SimplicialComplex`, `cl_Chain`, `cl_Cochain`, `fn_Smith` in `src/homology/`. Session-level
  edit approval does not reach them; only Gregory Giard's named authorization does. Reading and
  reporting stay open. See the protocol §7.1.
- Investigation and review tasks are read-only by default.
- Writing to `./todo/` and `./devlog/` is always allowed (task plans, AI exchange, session devlogs).

### Key Principle: Calibrated Uncertainty

Both AIs communicate uncertainty honestly rather than defaulting to confident agreement. Claims carry a confidence level (high / medium / low) so the partner AI and the user can assess and challenge them independently. See the collaboration protocol for details.

### Project Context

- For codebase architecture, coding standards, and literature routing, see `CLAUDE.md`.
- For documentation organization conventions, see `doc/documentation_guidelines.md`.
- For the HPC design rationale, see `doc/coding_philosophy.md`.
