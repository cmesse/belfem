# Debt Register: Archive Sweep of the Day's Strikes

**Date:** 2026-08-27
**Purpose:** move struck rows out of the live debt register into the archive
**Module:** cross-cutting (todo/)

## What was done

Six rows carried a struck ID in `todo/debt_register.md` but had never been moved to
`todo/debt_register_closed.md`. They were moved verbatim, in register order:

| row | why it was struck |
|---|---|
| DR-40 | costheta run + both clauses verified 2026-08-23 |
| DR-79 | MUMPS BLR compression closed 2026-08-26, gate executed |
| DR-82 | HDF5 memory-datatype corruption fixed, probe gate passed 2026-08-27; hardening spun off as DR-116 |
| DR-83 | `BELFEM_DEBUG` → `DEBUG` rename closed 2026-08-26, all four sites executed |
| DR-107 | homology hygiene bundle, `make check` gate discharged three times |
| DR-114 | default-config build break fixed by Christian, `make` clean |

The live register is now **32 rows**, the archive **84**.

## How it was done

Rows moved with their cell padding collapsed to match the archive's formatting, nothing else
touched — no strike markup was added or removed, no status cell edited. A whitespace-normalized
set diff of all table rows before and after the move is empty, so the split is lossless.

Three of the six (DR-40, DR-79, DR-83) have a struck ID but an unstruck *description*. They were
moved as written rather than tidied: the archive already contains rows in that mixed form, and
the register's own rule makes the status cell the verdict on whether an item is open. All three
status cells read closed.

## Trackers refreshed in the same pass

- `todo/debt_register.md` preamble: the freeze-lens counts now read 19 `[P]` (was 20, DR-114
  archived), 6 `[W]`, 3 `[F]`, plus the four untagged rows (DR-85, DR-94, DR-100, DR-102) that
  the sentence had never named. DR-100 is flagged as fixed but still live here — it is not struck.
- `todo/run_gate_batches.md`: DR-40, DR-83 and DR-107 removed from the "strike owed" lists
  (DR-75 and DR-92 still owe one); DR-79's `[RUN-BLOCKED]` row in §6 struck in place, with its
  garber A/B residue noted as riding the closed row.
- `todo/README.md`: the archive bullet now states the current 84/32 split and names the six rows.

## Not done

No source was read and no gate was run — this was a bookkeeping pass over the register files
only. **Struck is not verified** still applies to every row moved; the surviving gates (DR-79's
garber A/B in particular) live in the archived status cells and are waivable only by Christian.
