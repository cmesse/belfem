# DR-39 Netlist Question → ngspice Parser Plan Promoted, Phase-0 Decisions Taken

**Date:** 2026-08-27
**Purpose:** Record the promotion of `todo/ngspice_parser_plan.md` out of `deferred/` and the
Phase-0 design sign-off. No source touched; markdown only.
**Module:** `src/circuit` (planning)

## What happened

Christian asked whether DR-39's "(c) there is no netlist parser at all" means BELFEM has no
SPICE parser (it does mean that — verified 2026-08-13: circuits are built exclusively from the
`circuit{}` deck section via `ElectricalCircuitFactory::read_circuit()`), and asked to start
planning one, with the goal of the circuit example reading from a spice-like file. The "md file
with ideas" is `todo/deferred/ngspice_parser_plan.md`, the full 2026-06-28 proposal parked in
the 2026-08-26 sweep for lack of consumer pressure. This request is that consumer pressure, so
the plan was promoted back to the active set (`git mv` → `todo/ngspice_parser_plan.md`).

## Phase-0 decisions (Christian, in-session)

1. **FEM linkage = 6b** — terminal pairs stay in `input.conf` beside `topology{curves{}}`; the
   `.cir` holds only lumped elements. Decision history worth keeping: Christian initially picked
   6a because he hoped a `.cir` would let ngspice visualize the circuit. Two clarifications
   changed the picture: (a) the ngspice-loadability argument does not separate 6a from 6b at
   all, since the Option-A extension lines are comments and both variants keep the file valid
   SPICE; (b) ngspice is a simulator, not a schematic viewer — a netlist buys cross-simulation
   of the lumped subset (with the terminal pairs vanishing from ngspice's view), not a circuit
   drawing (confidence: high on ngspice, medium on the wider netlist→schematic tool landscape).
   With the compatibility argument moot, 6b's separation of concerns decided it.
2. **BELFEM-only lumped elements = Option A** comment directives (`* belfem: superconductor …`,
   timed switch). Under 6b the `terminal_pair` directive form is unused.
3. **No `.subckt` in v1** — hard-error on `X` cards; re-scopes the `.subckt` leg of DR-39 the
   way the register itself suggested.
4. **v1 scope = Phases 1–3 + deck wiring** — `spice_number_to_si`, `NetlistParser`,
   `NgspiceCircuitFactory`, `circuit { file : ... }` in `hphirun`, and converting
   `examples/circuit`'s lumped part to a `.cir`. `circuitrun` (Phase 4) and Pulse/PWL (Phase 6)
   come later. Ground stays a remap layer (plan recommendation, unchallenged).

## Plan refresh done at promotion

- Header rewritten: Status OPEN, decisions recorded, `RESOLVED` annotations placed in §5, §6,
  §10, §11 per the template convention.
- **DR-39 2026-08-13 verification folded in:** two Newton-loop copies, not three —
  `Controller::solve_circuit()` (`cl_FEM_Controller.cpp:1981-2039`) and
  `electricalCircuit.cpp:131-159` — with the one real drift (controller updates `tEpsilon0` at
  the top of the loop, executable at the bottom) flagged for the Phase-4 `CircuitSolver`
  extraction to resolve deliberately. Older line citations in the plan marked advisory.
- **DR-115 wired in as a sequencing gate:** `examples/circuit` currently aborts at np=4 on the
  first timestep rejection (`shift_back()` `Vector` overrun); it must be fixed or consciously
  waived before that example becomes the netlist showcase.
- Noted that `tmp/ngspice-manual.pdf` was in ephemeral `tmp/` — re-fetch at implementation start.
- `todo/README.md`: promotion entry added under 2026-08-27, registration entry re-statused, the
  two live links to the old `deferred/` path repaired (historical sweep prose left as written).

## Audit round 1 (same session, Christian's go)

Codex + Grok ran blind in parallel on the refreshed plan (pre-registration first, exchange
`tmp/ai_exchange/ngspice_parser_plan_audit.md`). **Verdict: not implementation-ready as
promoted; now repaired in place.** Both independently found the same three P0s, all adopted
into the plan's new §12 as O-items: **O5** — deck-side node references (`node +`, and the
`output{voltages}` sibling) are raw ints today (`cl_ElectricalCircuitFactory.cpp:76-86,632-636`)
and become ambiguous against the netlist's name→index remap; adopted contract: under
`circuit{file:}` every deck node reference is a netlist node *name*, and component labels =
SPICE instance names with one case-fold policy. **O8** — SPICE cards cannot carry the BDF
`order : 2` the showcase sets on both L and C; open for Christian (directive vs deck block vs
BDF1 waiver; directive recommended). **O7** — SPICE `.tran tstep` is the *print* increment,
not the initial Δt (the plan's mapping was wrong SPICE); coupled v1 ignores `.tran` with a log
line, `solver{timestep{}}` owns time. Plus P1s: hybrid factory ownership/sequencing (O6,
including `hphiTrun` needing the same branch), deterministic node-packing + unknown-current
order frozen for restart (O9 — `load_state` validates dof count only), a parser
error-reporting contract, the units seam (`100F` = 100 *femto*farad — the showcase capacitor
must be written `100`, not `100F`), and v1 hard-errors on all value-changing SPICE features
(resolves §11's warn-vs-error question).

Stale plan text struck in place: §3 `S`/`W`→`create_switch` and `X`→flatten rows, §9.2's 6a
"conf *or* cir" language, §7's "three Newton copies", the header's "pulls `cl_MaxwellFactory`"
claim (the executable is include-clean FEM-free since the Cut-A cleanup), and the `GROUND`
ground alias (disputed between auditors; v1 hard-errors pending the re-fetched manual). The
2026-08-13 "tEpsilon0 behavioural drift" was downgraded: Grok's unrefuted analysis shows the
top-vs-bottom assignment sites are behaviourally equivalent for the ω update.

**By-catch, all Claude-tree-verified before filing:**
- **DR-115 root cause diagnosed statically** (Grok, unprompted): `mX` is sized only in
  rank-0-only `compute_MNA_matrix()` (`cl_ElectricalCircuit.cpp:544`,
  `cl_FEM_Controller.cpp:244-249`) while `reset_timestep` runs `shift_back()` on every rank
  (`:2368-2370`) into `update_components()`'s `mX(j)` reads (`:252-255`) — empty `mX` on
  non-root ranks, matching "three of four ranks, serial clean" exactly. Register row updated;
  the "history depth" hunch was the wrong neighborhood. Fix still owes its own round; the plan
  now records DR-115 as a FIX gate (not waivable — `Allrun` defaults to multi-rank).
- **DR-117 filed:** `set_periodic` phase offset `(aPhase/2*constant::pi)*aPeriod`
  (`cl_SourceFunction.cpp:175`) is off by π²; latent (sine ignores `time_offset`), live for
  square/triangle/sawtooth with nonzero phase.
- **DR-118 filed:** `electricalCircuit.cpp` reject path advances time on failed steps (no
  `continue`, `:160-180`), aborts on a first-step rejection (`ShiftRegister::revert`), and
  re-stamps order-2 L/C with BDF1 companions on Δt changes (`cl_Inductor.cpp:81-85` vs
  `:55-56`). The plan's old "FAITHFUL" pseudocode had silently repaired the missing `continue`
  while claiming fidelity — now a correction banner. Standalone executable only; the
  controller path stamps correctly.

Grok stayed read-only (git status checked). Reviewed, not verified — no gate ran; the
polarity question (BELFEM current-source stamp vs SPICE `n+`→`n-` convention) is explicitly
unproven and assigned to the Phase-3 cross-simulation test.

## Sign-off (same session)

Christian resolved all five O-items in one pass: **O8 = the `* belfem: order <instance> <n>`
directive** (per-instance, next to its card, `.cir` stays vendor-valid), and **O5/O6/O7/O9
approved as adopted** (netlist node names in deck references, one factory owner with
lumped-in-topology as hard error, `.tran` ignored by `hphirun`, first-appearance node packing
frozen for restart). All RESOLVED annotations are in the plan's §12; the design is settled.

## Phase 1 implemented + code-audited (same session, Christian's "let's begin")

The manual was re-fetched first (master git clone → `tmp/ngspice-manuals/manual.lyx`) and every
round-1 disputed semantic settled at the source: ground is `0` + auto-converted `gnd` ONLY (the
`ground` alias refuted — Codex was right), `mil` = 25.4e-6, atto `a` in the main table,
`.tran tstep` = print increment verbatim, `$`/`;`/`//` comments confirmed.

**Landed:** `src/circuit/fn_spice_number.{hpp,cpp}` (ngspice number grammar → SI: validated
prefix scan + locale-independent `std::from_chars` per the `fn_GT_parse.hpp` precedent,
longest-match case-folded suffixes, trailing unit letters, hard errors on everything else) and
a new `tests/circuit/` gtest suite (89 expectations / 21 tests, `fast`), wired into `check`
AND `check-fast` foreach lists with `set( LIBLIST circuit )`.

**Code-audit round (Codex + Grok, blind):** all findings verified against the tree and fixed
same session — (1) post-scale overflow hole (`1e300T` → `+inf` without error; `isfinite` moved
onto the product) found by both; (2) **missing `LIBLIST` — the suite would not have LINKED;
Grok caught it, Codex had passed CMake as clean, and the hand-linked probe structurally could
not see it**; (3) `strtod` locale hazard closed via `from_chars` (+ guarded `de_DE` regression
test, ran un-skipped); (4) magnitude policy pinned toolchain-independently after an empirical
surprise — libstdc++ 11's `from_chars` reports subnormals as out-of-range while GCC 12+ would
accept them, so the contract (zero or normal double range, before and after scaling; subnormals
reject loudly) is now enforced at our level; (5) A-is-atto-not-Ampere twin footgun documented +
tested (`1A` = 1e-18, not 1 Ampere); (6) QUIRK comments demoted to "BELFEM policy forced by the
manual's rules, not traced against an ngspice binary" — the Phase-3 cross-simulation gate owns
that residual. Exchange: `tmp/ai_exchange/ngspice_phase1_code_audit.md`.

**Verified by execution:** scratchpad probe mirroring the suite 1:1 against the compiled
function — round 2 **88/88 pass**, negative control demonstrated (probe binaries in
`tmp/probe_spice/`). Both test TUs syntax-clean against real gtest headers.

## Phase 1 gate + Phase 2 (same session, after Christian's `make check`)

Christian's in-tree `make check` came back **15/15 green** with the new `circuit` suite in the
fast set — Phase 1 ticked, verified by execution on the named gate.

**Phase 2 landed and was code-audited the same way:** `src/circuit/cl_NetlistParser.{hpp,cpp}`
(netlist text → IR: title line, all comment forms, `+` continuations, paren/comma/`=`
normalization, case-folding, first-appearance node list, `* belfem:` directives incl. O8
`order`, §12.2 hard errors with `source:line: card` context). Gate innovation this round: the
REAL gtest TUs were compiled standalone against the real sources and executed — no mirror
probe. Audit round (Codex + Grok, blind): the substantive defect was **Grok's F1 — nodes
introduced only by a directive's `n+`/`n-` (superconductor/switch) never reached
`node_names()`, silently breaking the O9 restart packing contract**; fixed with a
written-order harvest + a regression test that fails on the old code. The one genuine
disagreement — Codex demanded continuations NOT survive interleaved comment/blank lines,
Grok called the superset acceptable — was resolved by re-reading the manual firsthand: it
is self-contradictory (":2690 must immediately follow" vs ":2722 empty lines ignored",
":5908 comments anywhere"), so the superset stays, now documented as a BELFEM extension
with the tension named. Also fixed: `.end nonsense` no longer silently terminates,
`.endc`/`.ends`/`.endif` locked by test (Grok flagged the hole as P0-blast-radius if the
matcher is ever "simplified"), `* belfem :` colon-space recognized, `r==5`/duplicate kwargs
hard-error via a shared `store_kwarg`, error-message contract asserted, on-disk `Ascii`
ctor smoke test. Final: **47/47 pass** (26 parser + 21 Phase-1) on the standalone build of
the real TUs. Phase-3 obligations recorded in the exchange (unknown directive kinds, missing
values, `ic=`/`uic` refusal, deck-side label folding). Exchange:
`tmp/ai_exchange/ngspice_phase2_code_audit.md`.

## Phase 3 (same session, after the Phase-2 in-tree gate went green 47/47)

`cl_NgspiceCircuitFactory` landed: O9 node map ("0"/"gnd" one node, last index; others
first-appearance, logged), components created in netlist line order (elements + directives
merged), labels = folded instance names, `node_index()` as the O5 lookup,
produce-and-hand-over via `unique_ptr`. Gate: **two SOLVED analytic circuits** — a DC divider
(v(mid) = 5.0, pinning V-source first-node-is-plus) and a current-source polarity circuit
(+1000 V, confirming the SPICE n+→through→n− convention and closing the Phase-2 polarity risk
by execution). Audit round was complementary: **Codex caught the parser letting directive
nodes jump the O9 queue** past a still-pending card (Grok had verified that ordering as clean
— vendor disagreement resolved by the tree; directives are now statements that flush the card
and end continuation chains, with a failing-on-old-code regression); **Grok caught a
pre-existing P0 in `ElectricalCircuit::compute_MNA_matrix`** — the V-source unknown-current
column is counted only past other V-sources, so a switch created before a V shifts the stamp
onto the switch's dof (Claude verified: one SWITCH case in the file, jacobian only; MNA
`default: break`). **Filed as DR-121** (reachable from the deck path too; latent — no shipped
deck uses a switch); the netlist factory refuses switch-before-V layouts with an error naming
the row until the MNA fix gets its own round. Also fixed: exception-safe non-copyable
ownership, diode/superconductor positivity guards, single-digit order token (uint wraparound),
case-folded `node_index()`. Suite 66/66 standalone. Exchange:
`tmp/ai_exchange/ngspice_phase3_code_audit.md`.

## Final gates (same session)

The first in-tree `make check` after Phase 3 caught what the standalone probe structurally
cannot: missing include paths and LIBLIST entries in `tests/circuit/CMakeLists.txt` (the probe
passes `-I` flags and archives by hand — code verified, wiring not; same shape as the Phase-1
LIBLIST miss). Fixed (ode/sources/fem-kernel includes; `LIBLIST circuit ode sources`), then
**Christian's in-tree `test_circuit` ran 66/66 green — Phases 2 and 3 both closed, verified by
execution,** with the SUPERLU solves of the two analytic circuits visible in the run log.

**Post-gate flake, root-caused and fixed:** a rerun flipped `FileConstructorSmoke` red. The
mechanism is the pre-existing `Ascii` quirk Grok had flagged in the Phase-2 audit: relative
paths resolve against `getenv("PWD")` (`cl_Ascii.cpp:36`), a SHELL variable that goes stale
when ctest/make change the working directory — the test's `ofstream` wrote to the real cwd
while `Ascii` looked where PWD pointed. The test now uses an absolute path
(`::testing::TempDir()`), verified by running the suite under a deliberately poisoned
`PWD=/nonexistent`: 66/66. The `Ascii` PWD dependency itself is pre-existing, tree-wide, and
left untouched — worth a small DR if it ever bites a production reader. **Final gate after
the fix: Christian's `make check`, 15/15 suites green (fast label 10 tests, 8.7 s) — the
campaign's Phases 1-3 are closed on the canonical gate.**

## DR-115 fix round (continued session, Christian's "let's continue")

The Phase-3b gate itself: `Controller::reset_timestep`'s `shift_back()` was the ONLY
unguarded `mCircuit->` site in the controller (all 20 sites enumerated firsthand); the fix is
the one-line rank guard mirroring the documented rank-0 circuit contract, with a
do-not-widen rationale (hardening `ElectricalCircuit::shift_back` would silence the debug
tripwire — both vendors independently agreed). Pre-registered with falsifier, TU
syntax-gated with real build flags, combined light plan+code round: **Codex CONFIRM + Grok
CONFIRM, zero blocking findings** — no collective below `shift_back`, no pre-MNA path on
rank 0 (all 17 `reset_timestep` callers tabulated), `hphiTrun` shares the controller,
`mReset` keeps ranks exiting loops together. Sharp Grok detail: the factory stores
`mCommRank` and never reads it, so every rank holds a live circuit — which is why workers
crashed instead of seeing null. **By-catch filed: DR-123** (switch latch not reverted on a
rejected step) and **DR-124** (first post-restart `compute_MNA_matrix` zeros the freshly
loaded `/circuit` state — the restore may be a silent no-op; static-only, cross-referenced
into `todo/restart_circuit_verification.md` as the expected R6 failure). The hunk is
independent of the concurrent thermal-watchdog work in the same file (both vendors checked).
**Owed: the run gate — np=4 `./Allrun` in `examples/circuit` through a timestep rejection.**

## Phase 3b: hybrid deck wiring (same continued session, alongside the DR-115 gate)

`circuit { file : <netlist> }` landed INSIDE `ElectricalCircuitFactory`
(`read_circuit_from_netlist`) — so `hphirun` and `hphiTrun` both get netlist support with
zero executable changes (O6 by construction). The classic TERMINALPAIR case and output block
were extracted into shared helpers with alias-preserved bodies; hybrid node references
resolve as netlist node names (O5), `.tran` is logged-and-ignored on rank 0 (O7), and the
guards (`number of nodes` forbidden, lumped-in-topology refused, unknown names refused) all
hard-error. Input Contract updated in BOTH artifacts same session (schema `file` key +
reference §11.1 with the worked example). Audit round (Codex + Grok, converged): the real
catch was a **one-sided case fold** — currents folded but deck pair labels did not, so the
documented `Z1` example would have failed at save time; hybrid mode is now case-insensitive
throughout, with cross-source duplicate labels refused. Also fixed: `unique_ptr` factory
members (ctor-throw leak), rank-guarded node-map logging, mode-conditional schema machine
fields, and a rank-0 voltage-column decoder line. **74/74 on the standalone real-TU build**
(the mixed-case `save_timestep` round trip is the regression that would have caught the P1).
Showcase conversion drafted and HELD until the DR-115 np=4 gate frees `examples/circuit`.

## DR-115 struck; showcase converted (the day's finish line)

Christian's np=4 `./Allrun` on the classic deck **reached and survived a timestep
rejection** — the exact former crash point — discharging DR-115's gate; the row is struck
and archived with the fix evidence (by-catch DR-123/DR-124 stay open). `make check` then
went green with the 74-test suite and the fem-kernel LIBLIST chain. With both gates in
hand, `examples/circuit` was converted: the lumped block became `tapestack.cir` (SIN
source, L/C with `* belfem: order` directives) and `input.conf` keeps only the two
terminal pairs (node NAMES `live`/`0`) and the output list. Pre-flight probe against the
installed production deck: **identical node-index mapping to the old deck** (live→0,
ground→1), same five components, two BCs in deck order — FEM coupling and CSV columns
continuous; only the labels fold to lower case. The `./Allrun` smoke of the converted
showcase is the campaign's final v1 gate.

## v1 SHIPPED (the day's last line)

After a stale-binary false alarm (`make check` rebuilds tests, not solver binaries — the
first showcase attempt aborted in the OLD constructor; timestamps proved it, one `make`
fixed it), **Christian's `./Allrun` on the converted showcase runs: production `hphirun`
builds its circuit from `tapestack.cir`.** DR-39's "there is no netlist parser at all" was
asked as a question in the morning and discharged as shipped code by night: five audit
rounds, 74 tests, one pre-existing parallel crash fixed and struck (DR-115), five register
rows filed as by-catch, and a demo deck stock ngspice can read back. Item (c) of DR-39 is
annotated discharged; items (a)/(b) remain as the plan's Phases 4 and 6.

## What was owed at v1 close

Remaining post-v1: **Phase 4** (`circuitrun` + shared `CircuitSolver`, discharging DR-118),
**Phase 6** (PULSE/PWL), the classic terminal-pair deck test, cross-simulation against a
real ngspice binary, DR-121/123/124's own rounds. Superseded text below kept as written: (the boxes above tick on green (the Phase-2 box ticks on green — `test_circuit` is new,
no new test directory this time, but two new source files need the reconfigure); then
Phase 3 (`NgspiceCircuitFactory` + O5/O6 wiring) under the same discipline. The DR-115 fix round before Phase 3b. Phase 3b additionally
owes the Input Contract update (`doc/input_file_reference.md` + `doc/input_schema.yaml` for
the new `circuit { file }` key and the O5 grammar change) in the same session that lands it.
