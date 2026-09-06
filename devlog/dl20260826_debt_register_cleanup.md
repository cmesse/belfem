# Debt Register Cleanup: Archive Split, Preamble → Operating Manual, and Waiting-On Tags

**Date:** 2026-08-26
**Purpose:** Record the split of `todo/debt_register.md` into a live register plus an archive,
and the replacement of its 233-line pass-by-pass preamble with a rules document.
**Module:** cross-cutting (todo/, doc/)

## What Christian asked for

> "Please move all closed DR entries into a new file `debt_register_closed.md`. And use the help
> of Codex and Grok to tidy up the lines 11 to 243. Codex+Grok to distill the essentials that we
> keep and Codex to do a prose sweep. **We are not interested in a changelog of the debt
> register, that's why we have git. What we need is a clear set of rules that we can work with.**"

That last sentence is the whole brief. The register's preamble had grown into a chronological
pile — "Correction pass 2026-08-09", "Maintenance pass 2026-08-10", "Strike-through pass
2026-08-11", "Second strike-through pass 2026-08-11, later the same day", "Currentness sweep
2026-08-13", "Cross-review pass 2026-08-19", "Triage pass 2026-08-23" — each recounting which
rows were touched, by whom, and how many. Real rules were in there, but only as asides inside
narrative paragraphs.

## The split

| | before | after |
|---|---|---|
| `todo/debt_register.md` | 353 lines, 108 rows, 233-line preamble | ~160 lines, **43 live rows**, ~115-line rules preamble |
| `todo/debt_register_closed.md` | — | 86 lines, **65 struck rows** |

Selection rule: a row moves iff its ID is struck (`~~DR-nn~~`), which is the register's own
definition of retired. Verified that no *unstruck* row is closed in its status cell — three
looked like candidates on a naive grep (DR-78, DR-93, DR-108) and all three were false hits, the
word "closed" appearing in the description or referring to a closed *audit round*.

**Round-trip check, run before and after assembly:** `sorted(orig) == sorted(live + closed)` →
`True`, 108 = 65 + 43. No row lost, reordered within its file, or mangled.

Row order is preserved as-is in both files, which means it stays roughly chronological rather
than numeric (the live table still runs `… DR-76, DR-79, DR-78, DR-77, DR-75 …`). Flagged to
Christian, deliberately not resorted.

## The three-way distillation

Codex and Grok each got the verbatim extract (`tmp/ai_exchange/debt_register_header_extract.md`)
and the same brief: inventory every rule buried in the narrative, verdict KEEP/DROP/MERGE on
each, propose a replacement. Claude drafted independently first, so the reconciliation was a
merge and not a rubber stamp.

**What each voice actually added** — the useful part of the record, because the split was
lopsided:

- **Claude's draft** had the four-section shape and the load-bearing rules (evidence is the tree;
  "pending run" is not closable by a static sweep; grep-shaped absence; agreement is weakest).
- **Codex (36 items)** added the closure-evidence *rung ladder*, "gate on the behaviour the fix
  actually touches" (a Jacobian fix is gated by iteration count, not physics output), "a status
  cell says what a run did **not** cover", and the observation that the *"fixed-but-still-gating
  rows stay live"* rule had **stopped being absolute** once DR-42/DR-49/DR-45 were struck as
  exceptions. Its restatement — live iff design/code/ruling residue remains — is what landed.
  It also flagged that "evidence is the tree, never the plan and never the devlog" is too broad
  as written: rulings and waivers *are* historical decisions, and the record is where they live.
- **Grok (57 items)** did what neither of the other two did: it checked the old preamble's claims
  **against the live table** and found two of them false.

## The two false claims Grok caught (both verified before acting)

1. **`[seeded — confirm]`.** Header lines 4-6 assert that *"every row is `[seeded — confirm]`
   until Christian's correction pass"*. `grep` finds the tag at exactly two places in the file —
   lines 5 and 13, **both inside the preamble prose**. **Zero table rows carry it.** Claude's own
   draft had faithfully carried this dead rule forward as a live one. Dropped.
2. **The 2026-08-23 triage census** (`14 run-only / 5 ruling-only / 14 mixed / 14 genuine code`,
   plus a blocking-row ID list and a cheapest-batch plan) is a dated snapshot of status columns
   *as written that day, explicitly without re-verification*. Several rows it names are since
   struck — DR-02, DR-06, DR-16, DR-17, DR-18, DR-19, DR-22, DR-23. The **taxonomy** survives as
   the triage lens; the census is dropped and deliberately not replaced with a fresh count.

Grok also caught that **lines 1-10 were themselves rot** and would have kept lying above the new
rules. Claude had explicitly planned to leave them untouched as out of the stated 11-243 scope.
They were rewritten.

**Agreement did no work in this round, again.** Every one of the highest-value findings came from
a *single* voice — the ladder and the not-absolute-anymore rule from Codex alone, both false
claims from Grok alone. Protocol §11's ranking of AI agreement as the weakest evidence tier held
literally, which is now itself one of the rules in the file.

## What the preamble says now

Seven sections, all imperative: **how to read a row** (struck vs live, status cell is the verdict,
the row supersedes the preamble, closing status and striking are two actions) · **when to strike**
(no residual; the DR-42/DR-49 recorded-exception pattern; compact a row whose history outgrows
readability) · **how a sweep gathers evidence** (the tree, not the plan, not the devlog; a stale
row is neither reliably pessimistic nor reliably optimistic; rewrite a refuted claim to the debt
that survives) · **closure evidence** (the rung ladder, reviewed vs verified, single-raiser flags,
pick a gate that can observe the defect) · **five standing traps** · **triage** (run / ruling /
code).

Four IDs are kept as canonical examples of durable lessons, and only those: **DR-19**
(grep-shaped absence — the false positive that survived a jury round, a row edit and two
independent reads), **DR-42/DR-49** (the strike-exception pattern the register names after them),
**DR-45** (count the sites in source: three `INFO` sites, not two), and the commit **`61bbfa11`**
with **DR-24** as the explicit anti-citation. All four now live in the archive file, which the
preamble points at.

## The prose sweep (Codex, read-only)

The finished preamble went back to Codex for the readability pass, told what it may **not**
change: technical claims, identifiers, commit shas, DR IDs, the markdown structure, the em-dash
house style — and *density beats readability where they conflict*, because this file is read by
us, not by users of the code.

It came back with four flags rather than pure polish, and three were accepted as written:

1. **"The row supersedes this preamble" was too broad.** The rule it came from
   (*"see the row itself, which supersedes this table"*) was about the pass-narrative table, not
   about operating rules. Letting a row override the rules would be new legislation. Narrowed to
   *"the row supersedes **examples** in this preamble"*.
2. **The rung ladder's top rung was over-citable.** A bare "production run" invites any campaign
   log to be waved at any row; tightened to *"named production/end-to-end run **that can observe
   the row**"*, which matches the "pick a gate that can observe the defect" rule two paragraphs
   below it. Codex also added *"a strike does not move a row up this ladder; only the status
   evidence does"* — a third restatement of struck-is-not-verified, which is the one confusion
   this file exists to prevent.
3. **The Newton-iteration-count sentence is an example, not a universal rule**, and now says so.

One flag is **carried to Christian rather than resolved**, because Codex is right on the facts:

4. **The compaction rule is new legislation.** *"When a row's amendment history outgrows human
   readability, compact it and say so"* is **not** in lines 11-243. It was lifted from the DR-02
   table row, which records exactly that happening — DR-02's four generations of amendment grew
   past readability and Christian had to run Codex over it to compact them. So it is real recorded
   practice, and it is the only rule in the new preamble that came from outside the range Christian
   scoped. Kept, and flagged here rather than smuggled.

One edit was **rejected**: Codex softened the rule to *"On its own, 'pending run' is not a row a
static sweep can close."* The brief forbade added hedges, and that rule earns its bluntness. The
hedge was removed before applying.

## Collateral

- `doc/ai_collaboration_protocol.md:352` and `:399-400` both described the register as "one
  table", which this change makes false. Fixed on Christian's approval: the file table gains an
  archive row, and §11's same-session rule now names both files and points at the register's
  preamble as the authority on reading and striking rows.
- `todo/README.md` — archive registered under *Added 2026-08-26*.

## Waiting-on tags and the run batches (same session)

Christian: *"tag all DR entries that are just pending a run. Then we can run them on auto pilot
tomorrow."*

All 44 live rows now carry a tag in the **ID cell** — not a seventh column, because eight rows have
unescaped `|` in their text and a new column would render inconsistently for exactly those. Tags:
`[RUN]`, `[RUN-BLOCKED]`, `[RULING]`, `[CODE]`, `[MIXED]`, with the legend in the preamble's Triage
section. Placed by reading every live row's status and reproducer cells — **a classification pass
from the status columns as written, no source re-verified**, which the register's own rules require
be said out loud.

| tag | rows |
|---|---:|
| `[RUN]` | **2** — DR-107, DR-87 |
| `[RUN-BLOCKED]` | 3 — DR-79, DR-90, DR-105 |
| `[RULING]` | **15** |
| `[CODE]` | 15 |
| `[MIXED]` | 9 |

**The headline is that the autopilot list is two rows long**, and the largest class is waiting on a
sentence from Christian rather than on a machine — including five rows that are fixed, gated green
and merely uncommitted. That is the cheapest yield in the register and it needs no run at all.

**Three gates are blocked by evidence decay, not by anything technical.** Verified on disk rather
than taken from the rows:

- **DR-90** names a replay from the **t = 2300 ms** memdump. The only surviving dump is
  `cmake-build-debug/tapestack3d/memdump.hdf5`, and reading its HDF5 metadata gives
  `timestamp = 7.1` s / `running_timestep = 3349` — the t = 7.1 s campaign dump. The one the gate
  names is gone.
- **DR-105** names `build/tape_quench` with a tuned recipe measured over a six-variant matrix.
  **`build/` no longer exists**, so the deck and its logs are gone.
- **DR-79** needs a pre-fix binary to A/B against; its gate (a) already passed.

DR-89's X2 gate died the same way when the 50/100 ms dumps were deleted. That is **four** run gates
lost to deletion, which is why `[RUN-BLOCKED]` is a tag rather than a footnote, and why the legend
carries the instruction to name a gate's artifacts in the status cell and expect deletion to be the
default.

`todo/run_gate_batches.md` holds the session plan: Batch A (`make check`, discharges DR-107 and
regression-guards the uncommitted bundle), Batch B (three example decks, with the `RLC_Circuit`
np=4 crash expected and not worth the session), Batch C (DR-87's peak-RSS watch, which rides any
long run), the 15-row ruling list, and the blocked table.

Tree note: `cmake-build-debug` is configured `USE_TEST=ON`/`USE_DEBUG=ON`/`USE_PARDISO=ON` with both
binaries present, so Batch A and DR-101's PARDISO deck are live there. Nothing was built or run.

## Owed

- **`doc/lessons_learned_evidence.md:809`** cites `todo/debt_register.md:218-221` for INC-512
  (DR-01/03/04 struck with "premise dissolved" wording). Those rows are now in the archive, so
  the citation is stale — a fresh instance of the exact `file:line`-rot failure `CLAUDE.md`
  documents under the input contract. The durable fix is a searchable anchor, not a corrected
  line range. Not edited: outside the approved scope, and `lessons_learned*` is sweep-excluded.
- `doc/lessons_learned_evidence.md:731` has `### B15 — debt register` reading `(section
  missing)`. Pre-existing, unrelated to this change, noted while grepping.
- **Eight rows render as broken markdown** — unescaped `|` inside cell text, so the row splits
  into extra columns for every reader. Live: DR-84, DR-88, DR-89 (23 pipes against a 7-pipe
  header), DR-92, DR-102, DR-104, DR-108. Archived: DR-98. **All pre-existing** — verified by
  re-running the pipe count against `ac73de0d`, where the same eight are already broken (the
  other thirteen hits there were the preamble's own sub-tables, which this change deleted). This
  is INC-173's failure class recurring: the table is assumed to render because nobody reads the
  raw markdown. Not fixed — it means editing the description text of rows outside this session's
  scope. The fix is `\|` in each, and it is mechanical.
- Nothing was compiled or run this session. Every change is markdown.
