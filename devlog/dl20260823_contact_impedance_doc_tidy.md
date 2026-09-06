# Devlog 2026-08-23 — contact_impedance_theory.md reworked into a consistent theory note

**Date:** 2026-08-23
**Topic:** `src/fem/maxwell/doc/contact_impedance_theory.md` restructured
at Christian's request: one consistent present-tense theory document
instead of a revision history of corrected mistakes.
**AIs involved:** Claude (rework), Codex (language sweep, 30 edits + one
flagged reference overstatement, all applied after review)
**Claude Confidence:** high (no technical content changed)

## Changes

- Header `**Correction:**` line removed; `**Related:**` no longer cites
  the todo design note (docs-don't-cite-todos rule).
- §4 rewritten from "Correction: rho/h, not rho*h" into a direct
  statement of the two reciprocal scalings (R_ct = ρ·h measurable,
  ρ/h in the matrix; TSA does not change conditioning).
- §9's "Codex review correction (…devlog…)" preamble replaced by a plain
  scope statement (also removed a forbidden devlog citation); §10/§11
  freed of "earlier discussion / originally assumed / we introduced"
  framing; "honest cap" unified to "resistivity cap"; Milestone-C
  bookkeeping replaced by inline statements of what is outstanding.
- Duplicate `## 11` fixed (References is §12); §12's Juntunen & Stenberg
  parenthetical narrowed to the scalar-elliptic scope §3.4 actually
  claims (Codex flag).
- Codex sweep applied (CAPS emphasis, correction-work connectives,
  em-dashes removed); every edit checked against the text before
  keeping.

No derivation, number, table value, or uncertainty marker changed. The
worktree copy under `.claude/worktrees/debt-register-sweep/` was left
untouched.
