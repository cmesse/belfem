# The rows were live and the sentences were wrong: a currentness sweep of the register and `todo/`

**Date:** 2026-08-30
**Purpose:** record a devlog-driven currentness pass over the live debt register and the `todo/`
directory — what this week's work actually discharged, and what the register still says about it
**Module:** `todo/` (no source touched)
**AIs involved:** Claude (pre-registration + tree checks + edits), Codex (`gpt-5.6-terra`, high),
Grok (`grok-4.6`, high) — blind, same brief, `tmp/ai_exchange/register_currentness_20260830.md`
**Verification:** static. Every claim below that concerns source state was re-read in the tree;
every claim that concerns an *event* (a gate ran) is cited to the devlog that recorded it.
**Nothing here is run-verified, and no row was struck.**

## What the session was asked to settle

For every open register row except DR-106, DR-144 and the homology rows: has this week's recorded
work made it solved or obsolete? Then the same question for `todo/*.md`. Seventeen rows in scope.

## The answer is not the one the question expects

**No row is strikeable, and that is not the interesting result.** The interesting result is that
**seven of seventeen rows are fixed in the tree while their description still states the defect in
the present tense**, and **eight `file:line` citations now name a different function than the one
the row is about**.

A previous sweep the same morning (`dl20260830_debt_register_strike_sweep.md`) asked "does any row
deserve a strike?", answered no, and stopped. That was the right answer to the wrong question. The
register's own rule says it explicitly:

> When a row's description and status cell disagree, **the status cell is the verdict** on whether
> the item is open — but a stale description is still currentness work, not permission to read past
> it.

A strike sweep and a currentness sweep are different passes. Only the second one finds this.

### The seven pre-fix descriptions

DR-97, DR-105, DR-110, DR-119, DR-120, DR-135, DR-139. Each now carries a dated **CURRENTNESS**
clause appended to its description naming what survives:

| row | the description says | the tree says |
|---|---|---|
| DR-97 | thermal budget "is silently inert in fully-coupled mode" | enforced — latched `tThermalBudgetSpent` at `cl_FEM_Controller.cpp:1860-1866` |
| DR-119 | thermal watchdog "structurally unable to fire" | coupled-mode clamp at `:4526-4559` |
| DR-120 | both keys "stored into `uint` without validation" | validated before the cast at `:3997-4018` / `:4198-4219` |
| DR-135 | `impose_dirichlet` "is the unrepaired sibling" | repaired 2026-08-28; the row is already `[RUN][P]` |
| DR-139 | `~ElectricalCircuit` "never deletes `mMNA`" | guarded `delete mMNA` at `cl_ElectricalCircuit.cpp:71-74` |
| DR-110 | "either add a capacity member **or** leave it documented" | the documented half is already taken, in the source at `cl_Mesh_Basis.cpp:315-321` |
| DR-105 | the quench-front `n > 1` abort is an open diagnosis | it has a landed remedy — see below |

The descriptions were **not** deleted. What each one records is the diagnosis, which is history
worth keeping; what was wrong is that a reader had no way to know the sentence was past tense. The
clause is appended rather than substituted for exactly that reason.

DR-139 also carried a note that its fix "sits inside the uncommitted circuit-cluster WIP that
session left behind at its exit". That WIP is committed (`b99ba890`, 2026-08-29) and
`git diff HEAD -- src/circuit/` is empty. The clause is struck in place.

### DR-105 is narrower than it reads

Its 9.8-hour serial run aborted at the quench front on `rho_piecewise`'s `n > 1` assert, in the
90 K < T ≤ 92.5 K window where the plugin's own fits have fallen back to 1.0 but `T_crit` keeps the
superconducting branch alive. The reproducer cell concludes "the deck-consistent value is
T_crit = 90.0".

The 2026-08-29 `usermat` round built exactly that: `hts_Tcrit` single-sourced, the fit cutoff made
inclusive (`T <= hts_Tcrit`), and a deck-side `critical temperature` key read in the shared tail of
`MaterialFactory` (`cl_MaterialFactory.cpp:283`) that overrides what a plugin's `_init` set. The
shipped deck carries it: `examples/tape_quench_usermat/input.conf:76` reads
`critical temperature : 90.0 K ;`. So DR-105's residue is a **rerun**, not an open diagnosis —
plus the still-unexplained np=8 floor-collapse, which is a separate thing the row also names.

### The citations that had rotted into other functions

Both auditors, independently, found load-bearing anchors pointing at unrelated code. Re-anchored:

- DR-97 `cl_FEM_Controller.cpp:3434` → a BLR-compression warning. Parse is at `:4198-4204`.
- DR-120 `:3368` → blank; `:3542` → a thermal print format.
- DR-135 `cl_FEM_SideSet.cpp:288` → inside the Element construction loop. `impose_dirichlet` is `:352`.
- DR-138 `tests/circuit/test_ElectricalCircuit.cpp:411` → the test is `:535`, its 1e-3 bound `:590`.
- DR-147 `cl_MaterialFactory.cpp:379` → a duplicate-label error. `data_file()` is `:558/:567/:580`.
- DR-111 cited `src/fem/maxwell/cl_ThinShellFactory.cpp:1517-1520`. **That path does not exist** —
  the factory moved to `src/fem/kernel/`; the winding is `:1543-1546`.
- DR-110's three `cl_MaxwellFactory.cpp` callers had drifted to `:1581,1830,1901` (+ `:1996,2102`).
  **See the correction below — my first reading of this one was wrong.**
- DR-39's four anchors, all shifted; `cl_SourceFunction.hpp` also moved to `src/numerics/sources/`.

This is the failure mode `doc/input_schema.yaml` already anchors against by token instead of line
number. The register does not, and cannot cheaply — its citations are prose. The mitigation that
does work is what happened here: a sweep that re-reads them.

DR-39's `.subckt` clause is struck in place — item (c) was discharged 2026-08-27 when the ngspice
parser shipped. Items (a) and (b) are confirmed still absent in the tree and are that plan's
Phases 4 and 6.

## DR-138 is Christian's call, not the sweep's

Its unit gate **ran and passed**: `dl20260829_circuit_cluster.md` records the
`RLCRingAcrossRejectedStep` rerun green in the full 15/15 suite, under 1e-3, exactly as
pre-registered. But the row itself still names "a coupled reject-retry A/B" as owed. The register
says a static sweep cannot close a pending run, and that a residue living in another artifact is
still a residue. So this is precisely the DR-42/DR-49 single-run exception — which is a ruling.
The sweep records the candidate and does not take it.

## `todo/`: four files moved, none deleted

Active set 48 → 44.

**Solved.** `ai_wrapper_tier_selection.md` (thirteen steps, R10 gate green, zero open boxes) and
`tapestack3d_jjc_noise_2125ms.md` (its own Status has read CLOSED and run-verified since 08-29;
residue is DR-128's own file).

**Obsolete.** `run_gate_batches.md` was derived from a classification of **44** live register rows;
the table now holds **24**, and Batches A–D are all discharged or struck. One thing in it is worth
more than the file: §6's rule that *a gate whose recipe is recorded survives the deletion of its
artifacts*, and that a campaign should say which output files its gates need before it runs. That
is **not** in `doc/lessons_learned.md`. Promoting it is a `doc/` edit, so it is flagged in the
file's closure banner and left for Christian rather than taken in a todo-only session.

`restart_timestep_consistency.md` is obsolete and **this is the find of the session, and it is
Grok's**. Both halves are gone: (a) the missing-BDF-triple hole was fixed and DR-112 struck earlier
today; (b) the Δt inconsistency was retracted as a false positive on 08-29, and the archived DR-112
row says in terms that this file "is superseded by that retraction". The plan's Status still read
"OPEN — R1 (localization) not started". Neither Claude nor Codex caught it; both of us went looking
in the *files* and Grok went looking in the *closed register*, which is where the supersession was
recorded.

**Not moved, flagged.** `make_doc_repair.md` split the auditors — its one open box is O5, a policy
question routed to Christian ("should the gate be able to fail a build?"). Grok read that as
SOLVED, Codex as still live. A sweep does not issue the ruling, so it stays active.
`current_sign_2d_fix.md` says DONE but owes O5 (DR-133) and the O7 voltage-driven run.
`shared_library_and_install_plan.md` and `gauge_newton_tangent_plan.md` both carry unticked boxes
their own prose says landed.

`todo/dr128_spap_low_field_axis.md` existed since 08-29 and was **never registered** in
`todo/README.md`. Now indexed. `todo/` has zero broken links and zero unindexed files.

## What the two auditors were each worth

They were not redundant, and the split is the useful record:

- **Codex** was the one that reframed the question. My pre-registration asked "is this row still
  live?" and answered correctly for every row; Codex asked "is this row's *text* still true?" and
  found the seven pre-fix descriptions I had walked past.
- **Grok** was the one that read outside the files it was pointed at — the closed register, not
  just the live one — and that is where `restart_timestep_consistency.md`'s death certificate was.
  It also caught the largest set of rotted citations.
- Both were wrong about something. Codex attributed a quote to `weekly_report_20260828.md`
  ("nothing owed from this file") that **is not in that file**; its §8 is titled "What is owed".
  I claimed Grok's replacement line numbers for DR-110 (`cl_MaxwellFactory.cpp:1581`, …) were
  invented. **That claim was false and it was my error — see the correction below.**

One near-miss worth recording: a naive cell count flagged DR-39's row as malformed markdown
(8 cells, not 6). It is not — the row contains `netlist\|subckt\|spice`, correctly escaped, and the
splitter was wrong. Verified against the pre-edit backup before touching anything.

## Files updated

- `todo/debt_register.md` — seven CURRENTNESS clauses, one struck WIP note, one struck `.subckt`
  clause, eight citation blocks re-anchored. Row count unchanged at 24; `[P]` 17 / `[W]` 7.
- `todo/README.md` — sweep section, four moved-file links, five stale Status lines struck and
  corrected, archive/live counts recounted (121 / 24), `dr128_spap_low_field_axis.md` registered,
  a currentness note on the legacy "Active Tasks" digest.
- `todo/closed/` — four files moved in with closure banners.
- `devlog/README.md`, this entry.

`scripts/check_doc_claims.py`: 37/37.

## Correction, same day: I called a true finding invented

The DR-110 caller list above is wrong in the first version of this entry, and so is what it says
about Grok.

Grok gave `cl_MaxwellFactory.cpp:1581, :1830, :1901, :1996, :2102` as the live `add_source` sites.
I ran `grep -rn 'add_source(' src/fem/kernel/cl_MaxwellFactory.cpp`, got nothing, ran a tree-wide
fallback piped through `head -6`, saw only `cl_Mesh_Basis`/`cl_ProtoMesh`/`cl_Mesh_Periodicity`
hits, and concluded the file had no callers at all. I then wrote that into the register as a
correction and into this devlog as a vendor error.

**Both of my greps were broken, in different ways.**

- `src/fem/kernel/cl_MaxwellFactory.cpp` **does not exist**. The file is
  `src/fem/maxwell/cl_MaxwellFactory.cpp`. `grep -rn` on a nonexistent path prints a warning to
  stderr and exits — which looks exactly like "no matches" when you are reading stdout.
- The fallback was truncated by `head -6` before it reached the `fem/maxwell/` entries, because
  directory order put `mesh/` first.

Grok's five line numbers are exact. Codex found this in the compaction round by checking the same
claim independently.

Two rules out of it, both about the shape of the mistake rather than the grep:

- **An empty grep result is not evidence until the path is confirmed to exist.** A wrong path and
  a real absence are indistinguishable on stdout. The same trap as `git diff | grep` swallowing a
  real change in INC-558, four hours earlier in this same tree.
- **Never conclude "zero occurrences" from a command with `head` in it.** `head` cannot tell you
  what it truncated. If the question is *does this exist anywhere*, the answer must come from
  `grep -c` or an untruncated list.

And a third, about the round rather than the tool: **I was quicker to believe a vendor had
fabricated a citation than that my own search was broken** — after a session whose entire subject
was stale citations. The finding was the kind Grok had been right about all day.

---

# Phase 2 — the compaction, and what the deletion audit caught

Christian's follow-up: *"Would it make sense to bring them up to date? Remove what is obsolete, so
that we can focus on the work that needs to be done. We have git and the devlogs, nothing will be
lost."* The register's own rule already licences it — *"When a row's amendment history outgrows
human readability, compact it and say so in the row."*

Eleven rows rewritten to state what REMAINS plus the minimum evidence needed to do that work and
judge its severity. Table text **79,656 → 68,888 characters**; DR-97 alone went 10,774 → 4,719.

## The audit found more in the deletions than the deletions removed

Codex (`gpt-5.6-terra`, xhigh) and Grok (`grok-4.6`, xhigh) got the before/after pairs and one
question: what did this drop that a future session would need? They came back with **eight false
statements in my compacted text**, seven of which I had inherited from the old rows and carried
forward without re-reading — the exact failure compaction was supposed to end.

| row | I wrote | the tree says | found by |
|---|---|---|---|
| DR-97 | "`run_coupled` requires BOTH residuals under target" | **there is no `run_coupled` in `src/`.** The accept is `( ! tSolveM ) && tThermalDone` (`:1747-1752`), and `tThermalDone` includes `mThermalStalled` — the thermal need NOT be under target, and the code's own comment calls that acceptance UNCERTIFIED | Grok |
| DR-105 | rerun the reconstructed deck | **it aborts at setup.** Its material block still says `custom { }` (`input.conf:66-70`) and the factory hard-errors on that spelling since the 08-29 rename (`cl_MaterialFactory.cpp:252-262`) | Grok |
| DR-105 | one deck | **two decks share the row's name.** T5 lives only in `cmake-build-claude/tape_quench_dr105`; the shipped `examples/tape_quench_usermat` has the `critical temperature` fix but is mumps at 1 ms with no Anderson — running it discharges nothing | Grok |
| DR-110 | "every caller allocates `number_of_nodes()`, i.e. 2-3" | false — `cl_Mesh_BfmFile.cpp:1524+` allocates a persisted `n`. The real invariant is capacity-matched allocate-then-append | both |
| DR-39 | Pulse/PWL have "zero tree-wide hits" | the netlist parser names and refuses both (`cl_NgspiceCircuitFactory.hpp:62`, `.cpp:432`) | both |
| DR-39 | the ε₀ placement is "one real behavioural drift" | structural only — an unrefuted equivalence argument in the parser plan shows the ω update is identical | Grok |
| DR-139 | "the class is implicitly copyable" | all four copy/move ops deleted at `cl_ElectricalCircuit.hpp:119-122`; DR-141 closed | Codex |
| DR-147 | call sites `:472/:328/:312` | those are `key_exists("units")` checks; the `read_user_defined` calls are `:484/:340/:324` | Codex |

**The DR-105 one is the session's most valuable finding in wall-clock terms.** The row's residue is
a ~60-hour serial rerun. It would have died in the first second on a rename diagnostic, and nothing
in the row said so.

## Tags were wrong in three places

`[RUN]` means "the only residue is executing a named gate — inputs exist, no code, no decision,
autopilot-ready". Measured against that: DR-105 still holds a ruling and an unexplained parallel
wall (`[MIXED]`), DR-111 has no fix at all (`[CODE]`), DR-135's probe deck does not exist
(`[RUN-BLOCKED]`). DR-119 and DR-120 earned their `[CODE]`→`[RUN]` promotion; Grok called DR-120
"the best `[RUN]` in the set". The kind-tag census is now CODE 11, RUN 7, MIXED 3, RUN-BLOCKED 2,
RULING 1.

## Two rows got LONGER, and that is correct

DR-39 (+68) and DR-110 (+234) grew, because the audit showed their old text was false and the true
statement is longer. Compaction is not a byte target; it is a currency operation that happens to
shrink most rows.

## Restored on purpose

Deletion is the risky direction, so the brief asked specifically for constraints, refuted variants,
measurements, artifact paths and named residuals. What came back and went in:

- **DR-97:** Grok's refutation of the unconditional budget form (self-caught before dispatch) — the
  reason the landed predicate is conservative, and the thing that gets re-proposed without it. Plus
  the A/B deck path, the 44 min 13 s cost, and both measured step shapes as gate criteria. The deck
  carries a trap: `examples/3D_tapestack` already sets `watchdog window : 8`, so the window must be
  defaulted to reproduce the envelope at all.
- **DR-105:** the np=8 numbers and log name, and T4 as the *refuted* softening lever.
- **DR-111:** `4a42d982` as a "not that bug" marker, and the 3D 10-step clean run behind "2D-specific".
- **DR-119:** the notice is gated on rank 0 *and* thermal-section presence, because a magnetic-only
  deck defaults `mIsFullyCoupled = true` (self-caught).
- **DR-138:** the coupled `reset_timestep()` → `initialize_*` `shift()` sequence — that *is* the
  recipe for the owed A/B, and the row had lost it.

Both auditors independently praised the same judgement: **delete votes, keep constraints.** Every
"do not re-propose this" and "do not overclaim that" sentence survived — DR-135's five residuals,
DR-138's snapshot-at-shift ban, DR-119's miss-is-a-spare-clause reading, DR-147's layering blocker.
What was cut was roster-and-vote chronology, which the register's own evidence ladder ranks last.

## Files updated

`todo/debt_register.md` — eleven rows rewritten, three retagged, eight false statements corrected.
Row count unchanged at 24; `[P]` 17 / `[W]` 7; `check_doc_claims` 37/37.
Exchange: `tmp/ai_exchange/register_compaction_20260830.md`.
