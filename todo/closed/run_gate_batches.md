# Run-Gate Batches

> **CLOSED 2026-08-30 — OBSOLETE, moved out of the active set by the register/todo currentness sweep.**
> This file was derived from a classification pass over **44** live register rows; the live table now holds
> **24**, and every batch it schedules is discharged or struck (Batch A / DR-107, Batch B, Batch C / DR-30,
> Batch D / DR-87; DR-90 and DR-79 are archived). Its tag census and its "waiting on Christian" list no
> longer describe the register. Kept for the record, not for execution.
>
> **One thing here has NOT been institutionalized and should outlive this file:** §6's evidence-decay rule
> — *say which output files a gate will need, and until when, before the campaign runs; a gate whose recipe
> is recorded survives the deletion of its artifacts.* It is not in `doc/lessons_learned.md`. Promoting it
> there is a `doc/` edit and was therefore left for Christian rather than taken in a todo-only session.

**Date:** 2026-08-26
**Purpose:** the executable residue of `debt_register.md`, batched for an autopilot session.
**Module:** cross-cutting

**Status:** ready to run. Derived from a classification pass over the 44 live register rows —
**from their status columns as written; no source was re-verified.** Every row carries its tag in
the ID cell (`[RUN]`, `[RUN-BLOCKED]`, `[RULING]`, `[CODE]`, `[MIXED]`).

## The headline, before the batches

Of 40 live rows, **one is genuinely run-only** — DR-107, and it is scheduled. The register's executable residue is much thinner
than the row count suggests, and the reason is worth stating up front:

| tag | rows | what it means for tomorrow |
|---|---:|---|
| `[RUN]` | 1 | autopilot discharges these |
| `[RUN-BLOCKED]` | 3 | a gate exists but its input is gone or its upstream is unfixed |
| `[RULING]` | 14 | **the largest class** — a sentence from Christian, zero machine time |
| `[CODE]` | 13 | needs the standing plan+audit → code+audit round |
| `[MIXED]` | 9 | the row names the order; usually a ruling before the run |

*Updated 2026-08-27, evening.* **The 02:00 `at` job ABORTED itself** (a solver was running at the
time — the guard worked). The homology bundle has since been **committed anyway** (it is in
`d8e5fa20`/`dcca23cc`), so DR-107's `make check` gate is now owed as a *post-commit regression
check*, not a pre-commit one. Tonight's autopilot (`~/belfem_autopilot/night_20260827.sh`,
launched by Claude on Christian's go: "do these runs on auto pilot", 8 MPI procs approved) runs:
Batch A (build at HEAD + `make check`), Batch B (RLC_Circuit — now also the post-merge
re-confirmation of the committed DR-100 fix `7a211d59` — and circuit), then a **reconstructed
DR-105 run**: the deck was rebuilt at `cmake-build-claude/tape_quench_dr105` from the surviving
`tmp/examples/Tape_Quench/CustomMat` copy + the three dl20260823 plugin repairs + the T5 recipe
recorded in the register row, so DR-105 has left §6 — the §6 "recipe survives deletion" lesson
did its job. Verdicts land in `~/belfem_autopilot/night_SUMMARY.txt`.

*Updated 2026-08-26, 19:30.* **DR-107's gate was scheduled: `at` job 19, Thu 2026-08-27 02:00**, running
`~/belfem_autopilot/dr107_gate.sh` (build the uncommitted bundle, then `make check`; verdict appended to
`~/belfem_autopilot/dr107_SUMMARY.txt`). **DR-87 and DR-31 were struck by ruling the same evening** — see
§4 and the archive. Earlier the same evening:* A parallel session fixed and struck **DR-108** (residual spun out as
**DR-110**, a ruling) and then **ran and struck DR-30** — its gate was discharged on a
bulk-superconductor substitute deck at 1 and 2 ranks, since the 2-rank thermal-coupled corc proved
unreachable. DR-30 briefly appeared here as `[RUN]`; it is retired and **Batch C below is
already done.**

**Fifteen rows are waiting on Christian, not on a machine.** Several are "fixed, gates green,
uncommitted" and need only a review and a strike. That is the cheapest yield available and it does
not need a run at all — see §4.

## 1. Batch A — the test suite (discharges DR-107) — **DONE, DR-107 STRUCK 2026-08-27**

Gate discharged post-commit by three green runs (overnight ctest 14/14 with homology passed;
morning `make check` after the DR-114 fix restored the dep chain; Christian's own `make check`).
Struck on Christian's confirmation. Kept below for the record.

Tree: `cmake-build-debug` (configured, `USE_TEST=ON`, `USE_DEBUG=ON`, `USE_PARDISO=ON`; `hphirun`
and `hphiTrun` present). **`build/` no longer exists** — anything that named it is dead.

```
cd cmake-build-debug && make -j<cores> && make check
```

- **DR-107** `[RUN]` — owes `make check` with **homology 9/9** before commit + strike. This is the
  whole gate. The homology cleanup (−1201/+135) is landed but uncommitted, so this run is also the
  regression check for it.
- Free by-catch: `make check-fast` covers the order-1 regression that DR-102 and DR-103 already
  gated green; a re-run confirms the uncommitted bundle still passes as a whole.

## 2. Batch B — example decks on the rebuilt binary

One rebuild, then in order of cost. None of these is a formal gate for a `[RUN]` row; they are the
regression battery for the **uncommitted bundle** (DR-94, DR-102, DR-103, DR-107 all sit
in the working tree unreviewed).

```
cd examples/RLC_Circuit && ./Allrun          # np=4 by default
cd examples/circuit     && ./Allrun
```

- `RLC_Circuit` at np=4 is **expected to crash** (DR-100, `[CODE]`, MPI_ERR_TRUNCATE in the SIZES
  collect). Reproducing it is not progress — the mechanism is already pinned. Run it only to
  confirm the fix round has not been overtaken; do not spend the session on it.

## 3. Batch C — DR-30 — **DISCHARGED, do not queue**

Struck 2026-08-26 by a parallel session. Kept here only so nobody re-queues it: the row's named
gate (2-rank thermal-coupled corc) was **unreachable** — no such deck exists and the scratch one
was gone — so the gate was executed on a **bulk-superconductor substitute deck, serial and 2
ranks**. Result: the aura path is live but `element_exists()` never fails, because the thermal
group contained every element the maxwell group did across all 144 716 aura reads. Corc terminal
voltages agree between 1 and 2 ranks to 1e-8…1e-6 on all six tapes. Fix judged correct, harmless
and defensive-in-practice.

Residuals stay open in the archived row (T absent from `mPostprocessorSourceFields`; `gTbulk` bias
if the groups ever diverge; the `norm(N*q)` vs `dot(Nvec,q)`+clamp mismatch).

The question this raised — whether the silent-keep in `link_element_maxwell_thermal` is unreachable
and could be an assert — is **ruled and closed: it stays.** It is the sanctioned fallback for a
maxwell element with no thermal counterpart, kept because coupled models added later may need it,
and the consumer discharges the obligation by comparing element ids. Latent is not dead. The intent
now lives in the comment at the link site, so the branch cannot be mistaken for an oversight.

**The lesson worth carrying:** an unreachable gate was discharged by substituting a deck that
exercises the same code path, rather than by waiting for the deck the row named. That is the
counter-move to the evidence decay in §6.

## 4. Batch D — DR-87 — **STRUCK BY RULING, do not queue**

Struck 2026-08-26 (Christian: "strike 87 on the existing evidence"), a recorded DR-42/49-style
exception. Its gate had already PASSED and been verified by execution; what remained was a **watch
with no defined pass threshold**. Three findings while preparing the run settled it:

1. **The deck is no longer in the configuration DR-87 was measured in.** The measurement was
   STRUMPACK magnetic + PETSc thermal; `input.conf:15-18` now reads
   `linear magnetic { library : mumps }`, switched during the 2026-08-24 four-run comparison that
   closed DR-02. A run would have measured a different allocator story and walked into DR-106.
2. **Memory:** 32 GiB available against a 47 GiB peak — the exact shape that killed DR-77 at
   51.6 GiB with swap exhaustion.
3. **Disk at 97%.**

**The watch survives as a standing caution, not as debt:** if in-step PEAK RSS creeps past ~53 GiB
the ratchet is back; next rungs are `mallopt(M_MMAP_THRESHOLD)` or jemalloc. PEAK is the metric —
the trough is sampling-limited at 60 s. It rides the next long STRUMPACK campaign. And the row's
standing warning stays live: **do NOT "fix" a symptom here with `krylov method : preonly`** — the
factorization is degraded to ~4 digits and the outer GMRES is what reaches 1e-11.

**DR-31 was struck in the same ruling** and is the cleanest exception yet: its gate is not merely
unrun but **unrunnable by design**. Order-2 h-φ terminates at an always-active
`BELFEM_ERROR( max_element_order() == 1, "Not implemented for higher order" )` in
`create_hanging_edges_and_facets` (`cl_MaxwellFactory.cpp:1399`), and every in-tree 3D h-φ deck
computes cohomologies. The fix stays dormant-but-committed; **struck ≠ verified applies with full
force — that code has never executed.**

## 5. Not a run at all — the cheapest yield (11 rows)

These need a sentence, not a machine. Grouped by what the sentence is:

- **Fixed, gated, uncommitted — review and strike:** DR-94, DR-102, DR-103 (~~DR-107~~ — struck
  2026-08-27 once Batch A went green, and archived to `debt_register_closed.md`).
- **Gate discharged, strike owed:** DR-75, DR-92 (~~DR-40~~ and ~~DR-83~~ struck and archived to
  `debt_register_closed.md` 2026-08-27).
- **Scope calls:** ~~DR-53~~ (**struck 2026-08-26** — fuse flags are experimental and not deck-settable), DR-65, DR-89, DR-110 (add a capacity member to `Basis` or document the invariant as caller discipline — DR-108 measured the bytes as free) (accept the analog or let the gate ride a future cold
  start), ~~DR-31~~ — **struck 2026-08-26**, see §4.

## 6. Blocked, and why — do not queue these

| row | gate as written | why it cannot run |
|---|---|---|
| ~~**DR-90** `[RUN-BLOCKED]`~~ ( **STRUCK 2026-08-29**, see the amendment below ) | rel-tol A/B, ~2 h replay from the **t = 2300 ms** memdump | that dump is gone. The only surviving dump is `cmake-build-debug/tapestack3d/memdump.hdf5`, and it reads `timestamp = 7.1` s, `running_timestep = 3349` — the t = 7.1 s campaign dump, not the one the gate names |

> **DR-90 amended 2026-08-29 — the blockage was the smaller half of the problem.** The gate is not
> merely unrunnable, it is **inapplicable**: the magnetic field now runs `library : mumps`
> (`cmake-build-debug/tapestack3d/input.conf:18`), a direct solver, and the deck's
> `relative tolerance : 1e-10` is annotated `// for PETSC` and never reaches it. There is no
> iterative magnetic tolerance left to A/B. Meanwhile the conditioning half of the question is
> already answered on disk — `out.txt` from the 2026-08-29 run reports a magnetic conditioning
> number every timestep from two independent estimators (spectral ratio 6.03e9 to 7.4e10, MUMPS
> COND1 4.25e6 to 3.43e8), so the deleted dump is no longer needed for it. And the stall class the
> row describes does not reproduce: both rejections in that run are thermal-limited. Full reasoning
> and the proposed re-specified gate are in the DR-90 status cell of `debt_register.md`. **This is
> the fourth row lost to evidence decay and the first one recovered from it** — recovered only
> because the run logs happened to survive, which is exactly the habit §6 asks for below.
| ~~**DR-105** `[RUN-BLOCKED]`~~ | `build/tape_quench` at the tuned recipe (`initial timestep : 0.002 ms`, `maximum timestep : 0.02 ms`); ~60 h serial | **UNBLOCKED 2026-08-27** — deck reconstructed at `cmake-build-claude/tape_quench_dr105` from the recorded recipe; overnight np=8 run launched (see the update at the top) |
| ~~**DR-79** `[RUN-BLOCKED]`~~ ( struck and archived 2026-08-27; the residue below rides the closed row ) | garber pre/post A/B on MUMPS compression | needs a **pre-fix binary** to A/B against, and is only meaningful if compression was binding on that deck. Gate (a) already PASSED; this residue may be worth waiving rather than rebuilding history |

A **fourth** instance surfaced with DR-30 and is the recoverable kind, worth contrasting: its gate
needs a thermally-coupled corc deck, and the scratch one is gone — but the **recipe** was written
into DR-06's status cell, so the deck can be rebuilt in minutes. A gate whose recipe is recorded
survives the deletion of its artifacts. That is the cheap habit the rule below is asking for.

**DR-89 is the precedent, and it is the same failure three times over:** its X2 gate died when the
50/100 ms dumps were deleted. Three of this register's run gates have now been lost to evidence
decay rather than to anything technical.

**Therefore, before any campaign in §1–3: say which output files the gates will need and until
when.** The default on this machine is deletion.

## 7. Ordering for the session

1. Batch A. If `make check` is not green, stop — everything below assumes the uncommitted bundle is sound.
2. Batch B, cheapest first; expect the `RLC_Circuit` np=4 crash and move on.
3. Skip Batch C — discharged.
4. Take §5 to Christian — it is 15 rows for a handful of sentences and no machine time.
5. Nothing else needs the machine: Batches C and D are discharged or struck.
