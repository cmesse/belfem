# DR-78 Strike and DR-112 Loader Hardening

**Date:** 2026-08-27
**Scope:** `todo/debt_register.md`, `todo/debt_register_closed.md`, `todo/run_gate_batches.md`
**Status:** DR-78 closed; DR-112 opened

## Summary

Christian ruled to strike DR-78 and split the remaining loader-policy concern into a new debt row.
The DR-78 crash fix and its np>=2 gate are done. The separate policy question is that, in a 0.9
pre-release tree, obsolete warm-restart memdumps do not need a compatibility fallback.

## DR-78 Closure

DR-78 covered the BDF warm-restart crash when a dump carried integrator state but the field side did
not guarantee the numbered history fields. The practical high-risk path was parallel: workers could
enter BDF5 with empty `qold` storage and either crash under NDEBUG or assemble a wrong RHS.

The closure evidence is the 2026-08-24 np>=2 gate: a coupled corc scratch deck wrote a 235 MB
history-bearing dump at full BDF5, then a fresh two-rank warm restart resumed at full order with
the banner and advanced cleanly through steps 14-16 with zero errors.

The original serial `qold(1)` bounds throw remains mechanism-unpinned, but it is not kept as live
debt. The restored path now synchronizes and verifies the numbered history fields before accepting
BDF state, so the dangerous condition becomes a named loader error rather than a raw crash or
silent wrong RHS.

## DR-112 Split

DR-112 now tracks the stricter loader contract: a warm restart under BDF2-5 should reject an
obsolete or incomplete no-`bdf_*` memdump instead of cold-starting the BDF order ramp at BDF1.
Partial triples already hard-error; the new row is about removing the missing-triple fallback.

No source code was touched.
