# Test-gap sweep, DR-113 drain-check tests, grok-4.7

**Date:** 2026-09-21
**Purpose:** Spend the pre-reset budget on tests rather than a second documentation sweep:
close DR-113's deliberate-stray gate with a permanent test, confirm DR-157's state, and turn the
incident record into a ranked backlog of missing tests.
**AIs:** Claude Fable 5.1 (primary; five blind sweep agents), Claude Fable via web API (two
pre-sweeps, Christian's dispatch), Codex `gpt-5.6-terra`/high and Grok `grok-4.6`/high (DR-113
plan jury), Codex `gpt-5.6-terra`/high and Grok `grok-4.7`/high (DR-113 code jury, in flight at
the time of writing)

## Summary

- **Docs vs tests.** Christian asked whether a second documentation sweep or more tests was the
  better use of the remaining budget. The Roles sections of `AGENTS.md` and protocol §1 do omit
  Grok, a one-hour fix, but a general sweep has a low find rate after this month's rounds. Martin
  2025 settled the rest: an LLM gets roughly 80 % of a module right and the check is tests it
  did not grade itself (Ch. 13 "It Will Be Wrong", p. 287); intent must be stated twice, once as
  tests (Ch. 17 "Infancy", p. 324); coverage measures execution, not testing (Ch. 32, p. 529).
  The protocol's evidence ladder (§11) already ranks a focused regression second and reviewer
  agreement last, while the practice has been the reverse. Decision (Christian): finish DR-113
  and DR-157, then sweep the record for test gaps.
- **DR-157** was already struck on 2026-09-18: the `-DUSE_GASMODELS=ON` token is on the
  nightly's configure line. Only the observation is open, and the run is tonight. Two leftovers
  wait on it: `CLAUDE.md` still says the gas tests "have never run there (DR-157)", and the
  comment in `tests/physics/CMakeLists.txt` claims the gas tests skip themselves when the tables
  cannot be found, which no test does (they call the data-path lookup and open the file; the
  harness pins `BELFEM_DATA` to the source tree's `share/`, so CI is fine).
- **DR-113.** Two tests appended to `tests/comm/test_CommMPI.cpp`, registered through the
  existing `commmpi` suite at 2 and 4 ranks: `DrainCheck.CleanFabricIsSilent` (ordered pairwise
  exchange over both fabric tags, then `comm_drain_check` must not throw) and
  `DrainCheck.StrayScalarFiresOnEveryRank` (rank 1 posts an unmatched one-int `send`; every rank
  calls the check in a bounded collective retry loop, catches the debug-policy throw, checks its
  role-specific message, rank 0 drains the stray and a final check is silent; skipped in release
  where the check compiles out). Pre-registered in `tmp/ai_exchange/dr113_drain_check_test.md`
  with four predictions before any audit. Plan jury: Codex and Grok both REVISE, converging on
  every point and both ruling out a `src/` change. Applied: T1 walks a global `i < j` pair
  schedule (a rendezvous cycle was otherwise possible), post-catch checks are `EXPECT_*` so a
  failing rank still reaches the drain and barrier, and the reliance on Open MPI's eager
  completion of a one-int standard send (not an MPI guarantee; sm 4096 B, tcp 65536 B on the
  local 5.0.10) is stated in a comment with its failure mode, a ctest timeout rather than a
  silent pass. The stale "no throw expectation" comment in `test_commmpi_main.cpp` was
  refreshed in the same diff. Syntax-checked with the real `-Werror` flag set under the
  debug/Blaze and release/Armadillo trees. Code jury (Codex terra/high, Grok 4.7/high): both
  caught that my insertion had landed before the FIRST `#endif // BELFEM_MPI`, ahead of 30
  tests whose one-int collect on the same tag would have consumed the stray, exactly the
  failure the tripwire exists to stop. Fixed by rebuilding the file from the pristine copy;
  Grok's four smaller corrections (the "DR-100 shape" comment, a fatal length assert before a
  collective, the deadlock-freedom claim, the release comment) applied; Codex's doc point (the
  coding philosophy said every test main forces `throw_on_error`) closed with one sentence
  naming the Tier 2 exception. Open risk accepted and recorded: the payload tag is exercised,
  not asserted. **Reviewed, not verified:** `make check` in both trees is Christian's. Gate (a), a chokepoint at more than one rank,
  stays owed and needs the Tier 2 FEM fixture (backlog F1).
- **grok-4.7.** `grok models` lists it as the CLI default. Allowlisted in `ask_grok.sh` and
  `cross_review.sh`, named in the depth table (protocol §9.1, the cross-review command, the
  `--quick` pin); `grok-4.6`/`4.5` stay allowlisted for reproducing older rounds;
  `grok-4.7-build-fast` deliberately not listed. Two probes through the wrapper at `medium` and
  `xhigh` returned correctly stamped with the tree untouched. The wrapper default stays `grok-4.6`
  because it mirrors `~/.grok/config.toml`, which still says 4.6. No quality comparison against
  4.6 on a real brief was run; the table moved on the probe alone.
- **The sweep.** Corpus: the 601-row incident catalog (already the distillation of the devlogs
  to 2026-09-01, one sixth their size), the 49 September devlogs, the seven open register rows,
  and two web pre-sweeps Christian ran with the same brief. The brief
  (`tmp/ai_exchange/test_gap_sweep_prompt.md`) paraphrases the Martin and Thomas & Hunt
  criteria with chapter cites so the licensed texts never leave the machine, fixes an output
  schema, and caps each slice at 30 proposals. Five blind agents with read-only tree access
  covered the catalog in thread groups, the devlogs and the register. Reconciliation: every
  cited id resolves (352 incident, 73 register, 50 devlog cites, 0 missing); the "red today"
  claims at the top were re-read in source. Result: `todo/test_gap_backlog.md`, 77 tests in
  four tiers, 7 enabling fixtures, 9 rulings, and §2's 19 defects still in the tree; the same
  list flattened into `tmp/devlog_sweep_claude/fable_consolidated.md` as the Claude column of the
  cross-check against the Grok, Astra and Copilot sweeps (stable `T-nn` ids, a crosswalk from
  the pre-sweep rows, and every slice's covered/doubted sections verbatim).

## Key Findings

- **Two "latent" catalog rows are live defects.** `get_bottom_nodes`/`get_top_nodes` HEX20
  sizes 9 and reads `node( 21 )` on a 20-node element (INC-084, copied from the HEX27 case);
  `BlockData::collect_thin_shell_facet_ids` counts shells with a selected block and fills from
  every shell, an out-of-bounds write when a shell has no selected block (INC-072). Both were
  recorded "confirmed real but latent" and nobody re-grepped them. Both re-read in source today.
- Also confirmed alive by reading: `hdf5::group_exists`/`dataset_exists` both bare `H5Lexists`;
  `Map::operator()` release path aborts without the key; no QUARTIC branch in
  `interpolation_order`; the `.bfm` loader reads the writer's git hash and never compares it
  (DR-160); no plugin ABI guard; `Ascii` resolves relative paths against `$PWD`;
  `CutData::collect_coefficients` tests `== ±1` only. Reviewer-read, not re-read: the silent
  drop of `key : value // comment;` deck lines, `set_relative_tolerance` unvalidated, the
  `CurveFactory` two-chain sentinel, duplicate-id mesh globals, the HEX8TB twin-slot FIRST/LAST
  contradiction, the still-disabled absolute entropy test, and `Vector::operator=( expression )`
  resizing silently.
- **Fixtures decide the ranking.** No test constructs a `Controller`, a `Postprocessor`, a
  `BfmFile`, a `CutFactory` or a `MaxwellFactory`, and no Tier 2 suite builds a mesh. About a
  third of the P0 proposals sit behind one of seven fixtures; whether each can be built at all
  is the first question of its session. `test_DofSeeding.cpp` already builds a serial
  Mesh + Kernel + DofManager and `test_GhostSwitch.cpp` has a deck helper, which lowers several
  costs the web pre-sweep had priced as L.
- **The catalog's own "misaligned rows" were a false alarm of mine.** Rows in two other tables
  of the evidence file also start with an incident id (a verdict sample and the addendum); the
  catalog table is clean. The four catalog agents got a one-line correction mid-run.
- **The web pre-sweep's first pass could not see the tree** and marked seven reachability
  questions unknown; its second pass, with tree access, resolved all of them and withdrew two
  proposals (one closed by guard, one already covered). That is the behavior the schema was
  designed to force, and the reason the in-tree agents were briefed to grep before proposing.

## Changes Made

- `tests/comm/test_CommMPI.cpp`: `<stdexcept>` include, `DrainCheck.CleanFabricIsSilent`,
  `DrainCheck.StrayScalarFiresOnEveryRank` (last in the TU by design).
- `tests/comm/test_commmpi_main.cpp`: throw-policy comment refreshed.
- `.claude/scripts/ask_grok.sh`, `scripts/cross_review.sh`: `grok-4.7` allowlisted; `--quick`
  pinned to `grok-4.7`/`medium`.
- `doc/ai_collaboration_protocol.md` §9.1, `.claude/commands/cross-review.md`: depth table and
  examples name `grok-4.7`.
- `todo/test_gap_backlog.md` (new), `todo/README.md` (registered, active count 5).
- `tmp/ai_exchange/`: `test_gap_sweep_prompt.md`, five `test_gap_sweep_claude_slice_*.md`,
  `dr113_drain_check_test.md` (plan, two audits, reconciliation), `dr113_drain_check_test.diff`,
  `grok47_probe.md`.

## Open Questions

- Code jury verdict on the DR-113 diff, then `make check` in `cmake-build-debug` (Blaze, debug)
  and `build` (Armadillo, release); prediction P2 (every rank catches within a handful of
  sweeps) is the one that can fail for a reason worth knowing.
- Tonight's nightly: the ctest log must name `gastables` and `gasmodels` in all four matrix
  jobs; then the two DR-157 leftovers above.
- Rulings O1–O9 in the backlog, and whether L1–L6 get register rows this week.
- Whether the depth table should wait for a measured 4.6-vs-4.7 comparison before naming 4.7.

## Files Updated

`tests/comm/test_CommMPI.cpp`, `tests/comm/test_commmpi_main.cpp`, `.claude/scripts/ask_grok.sh`,
`scripts/cross_review.sh`, `doc/ai_collaboration_protocol.md`, `.claude/commands/cross-review.md`,
`todo/test_gap_backlog.md`, `todo/README.md`, `devlog/README.md`, this file.
