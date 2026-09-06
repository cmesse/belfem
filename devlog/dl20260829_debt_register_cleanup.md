# Devlog 2026-08-29 — Debt Register Cleanup

**Date:** 2026-08-29
**Topic:** Archive struck debt-register rows and audit apparent solved rows
**AIs involved:** Codex
**Codex Audit Confidence:** high
**Verification:** mechanical register audit

## Summary

Reconciled the live register after the day's parallel closure rounds. Every struck row is now in
`todo/debt_register_closed.md`; the live register contains no struck rows.

The status cells of the fixed but unstruck rows were audited against the register's closure rules.
None is silently complete: DR-97, DR-119, DR-120, DR-135, DR-138, and DR-139 retain named run
gates; DR-144 retains G10/G11 code work; and DR-146 explicitly retains its `expk` design half.
They therefore remain live unless Christian applies the DR-42/DR-49 single-run exception to a
run-only row.

## Changes Made

- Reconciled the archive after DR-86, DR-90, DR-117, DR-118, DR-121, DR-123, DR-124, DR-128,
  DR-140, DR-141, DR-142, DR-143, DR-145, and DR-148 were struck during the day's completed
  rounds.
- Preserved all closure evidence and surviving-gate text in the archived rows.
- Confirmed the live-register preamble counts against the actual table.

## Verification

- `rg '^\| ~~DR-' todo/debt_register.md` returns no rows.
- No ID is duplicated within the live table or the closed table, and no ID occurs in both files.
- The live table contains 27 rows: 2 `[F]`, 18 `[P]`, and 7 `[W]`.
