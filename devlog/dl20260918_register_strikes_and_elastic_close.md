# DR-131, DR-155 and DR-156 struck; the elastic-moduli plan moves to closed

**Date:** 2026-09-18
**Purpose:** Record three register strikes on Christian's ruling and the relocation of a completed
todo, with what each strike waives.
**Module:** `todo/` (register and plans); no source touched
**AIs involved:** Claude (register read-back and edits); no audit round — register and todo text
only, per the "comments free, code takes the jury" rule.

## Rulings

Christian, reading the register: DR-131 can be closed; DR-155 also looks closed; DR-156's
benchmark is not the register's scope. Separately: `todo/elastic_moduli_quasiharmonic.md` can move
to `closed/`.

## What each strike rests on

- **DR-131** (plugin templates and example decks on Darwin). Every gate ran on 2026-09-01 and the
  residual was discharged the same day: templates 4/4 green, all four decks fixed and rebuilt,
  `disk_pulse` run to completion on macOS with its plugin confirmed mapped into the live process.
  Nothing was owed. Struck on evidence.
- **DR-155** (`matvec_csc` OpenMP stack overflow, fixed by the default-OFF `BELFEM_OMP` token).
  Executed: the token-OFF debug build with `make check` 15/15, the fix confirmed present by symbol
  with a falsification control, and the Darwin decks without `OMP_STACKSIZE` under DR-131.
  **Waived with the strike, never run:** the `USE_BELFEM_OPENMP=ON` build, and a Linux
  `-DUSE_MKL=ON` eigenvalue run of the forwarded MKL path. The nightly job configures without MKL,
  so CI has not covered that path either. Switching the token ON re-arms the defect; that is
  documented in `doc/parallel_execution.md` and is not carried as debt.
- **DR-156** (ParMETIS / PT-Scotch never benchmarked). Struck by scope ruling: the wiring is landed,
  green at np 2/4, D5 fixed and gated. Whether the two orderings pay for their dependency is a
  measurement question and stays with W-R5 of `todo/closed/parmetis_ptscotch_wiring.md`.

## Edits

- `todo/debt_register.md`: three rows removed; `[P]` header count 8 → 5, recounted 2026-09-18
  with the reasons above. Live register now 9 rows (5 `[P]`, 3 `[W]`, one `[cdP]`, `[F]` empty).
- `todo/debt_register_closed.md`: the three rows appended cell for cell, ID and area struck, the
  status cell prefixed with the strike note and otherwise kept. Archive now 144 rows.
- `todo/closed/parmetis_ptscotch_wiring.md` and `todo/closed/matvec_csc_omp_stack_overflow.md`:
  Status lines no longer point at a live DR row; the matvec plan's "not yet compiled or tested"
  was stale against the 2026-09-01 gate and now says so.
- `todo/elastic_moduli_quasiharmonic.md` → `todo/closed/` (`git mv`). It landed as `2ffebcc`.
  Its four unticked definition-of-done boxes were ticked: the first three are plan-quality items
  the completed plan satisfies, and the gate box is covered by `make check` 17/17 with
  `CopperMatchesLedbetter1981Isothermal` and `PoissonRatioIsMonotoneAndBounded` in
  `tests/physics/test_MetalElastic.cpp`. `todo/README.md` entry re-pointed and marked CLOSED.
- `scripts/check_doc_claims.py`: 38/38 after the edits.

## Same session: DR-157's one token landed

Christian: running the gas models on the nightly is an easy patch, and the default stays off.
`.gitlab-ci.yml` now appends `-DUSE_GASMODELS=ON` to the nightly's configure line, with a comment
stating that `USE_GASMODELS` stays OFF in `CMakeLists.txt` for users. `USE_COMBUSTION` is untouched
and still OFF; it is the only gas-side option with a nonfree dependency. DR-157 retagged
`[MIXED]` → `[RUN]`; the strike waits on the next scheduled pipeline, whose ctest logs must name
`gastables` and `gasmodels` in all four matrix jobs. Risk noted in the row: the job's 1 h timeout at
`-j2` has never been measured against the two gas suites.

## Same session: DR-158 fixed (option B), jury-audited, configure gate run

Christian chose option B: take the Tier 2 launcher from the installation the compiler wrapper
belongs to, rather than refuse on a `--version` banner (which cannot tell two Open MPIs apart).
`config/scripts/Add_Test.cmake` now searches `MPI_HOME/bin`, then the wrapper's directory and its
symlink target's, then CMake's default search, each alone before the next, and prints the choice
once. The jury (Codex gpt-6-astra/high, Grok 4.6/high, blind, on a restricted diff of my three
files because another session was editing the kernel) converged on three findings, all confirmed
and fixed: the wrapper lookup used CMake's default search, which visits cmake's own install prefix
before `PATH` and so does not resolve the bare `mpicxx` the way `make` does (my own configure gate
had caught the same thing minutes earlier and it was recorded before the round returned, not fixed
mid-round); `MPI_HOME` proves nothing about what was linked, so its launcher is warned about when it
is not beside the wrapper; and my two documentation hunks described the order I had first written,
not the one I committed to. Grok's remaining points (warning wording, `NO_CACHE`, empty quotes in
the fatal message, a REALPATH trap under ccache, the Related list) were all taken. The selection
gate ran as four scratch configures of the full tree with a fake foreign `mpirun` first on `PATH`;
the record with the reconciliation table is `tmp/ai_exchange/review_dr158_launcher.md`. Row
retagged `[CODE][W]` to `[RUN][P]`; owed is `make check` in a real tree reconfigured with this file.
One lesson worth a line: a `pkill -f` whose pattern matches the relaunching shell kills the relaunch.

## Same session: DR-119's notice half ran in the sandbox

Christian asked for a proposal on DR-119 and then for the runnable half to be run here. The clamp
notice fires at setup and needs a coupled deck with a magnetic ceiling of at most 33 and a defaulted
thermal window; no shipped deck has that. `examples/tapestack3d` (magnetic ceiling 30, thermal min
iterations 2) was copied to the scratchpad, remeshed with `gmsh -3`, and its thermal `watchdog
window : 8` line removed. The release binary printed the notice, clamped to 27 against the ceiling
of 30, at the default info level and went on into the first BDF1 step before the run was stopped.
The log is `tmp/dr119_notice_run.log` (lines 79 and 80). The trip half cannot be scheduled: a stall
long enough to fire the watchdog is a property of a deck at a step, and forcing one would test the
probe. Row retagged `[RUN]` to `[RULING]`: the proposal is to strike on this run plus the existing
code approval, with the trip half waived as unproducible. A first attempt at the run produced an
empty log because the binary block-buffers stdout into a file; `stdbuf -oL` fixed that.

## Same session, later: DR-157 and DR-158 struck

Christian's ruling after the two fixes landed. DR-157's token is in; the pipeline observation is
the nightly's own business. DR-158's `make check` wiring run after a reconfigure is waived with the
strike, since the chosen launcher on this host is unchanged. Archive 146 rows, live register 7.

## Same session: DR-133 audited (astra + Grok) and probed — the diagnosis was wrong

Christian asked what can be done about DR-133 and whether it risks the cohomology behavior. I
pre-registered a counter-analysis: the facet's node order equals the master element's CCW edge
(`Facet::set_master` relinks, `flip()` re-links), the original nodes stay with the master, so
"antiparallel to the facet direction on the master side" is the counter-clockwise test the row said
was missing. Both auditors confirmed that reading and refuted the row's conclusion, and both
independently designed the same discriminator: negate both halves of one tape in
`tapestack2d_layered`, fresh cuts, judge by integrated `Jz`. That run had already been executed in the
sandbox on the Blaze release binary: tape 6's integrated current kept its sign to 0.06 % and the far
field was identical to 0.08 %. What the jury found instead: four never-written `Vector` components
(UB under Blaze, benign in the probe), a debug-only assert guarding a dereference on a setup path,
a first-segment lookup that depends on the first listed sideset id and the curve's start node, and
that the row's "BFS-deterministic" mechanism is false in 2-D (propagation is 3-D only; the maxwell
orientation doc is stale on that). Codex also found a latent LINE3 closure defect in
`close_terminal_loops`, unreachable behind the order-1 guard, now widened into DR-122. Risk to the
cohomology behavior: none from the proposed hygiene diff, which stays inside the dimension-guarded
2-D function; the one action that would break it is rewriting the 2-D test's sense, because the
global −1 in `Homology::reorient_generators` (closed core) is calibrated to the present sense. Row
retagged `[MIXED][P]`; the diff needs Christian's approval and a code jury before the strike. Record:
`tmp/ai_exchange/review_dr133_orientation.md`.

## Same session: the DR-133 hygiene diff, astra + Grok code round

Christian approved the diff. Landed in `cl_CutFactory.cpp`, 2-D function only: fill-constructed
3-vectors, `BELFEM_ERROR` for the null facet and for two zero lengths, a `BELFEM_ERROR` that the
matched facet's master element carries the segment's nodes, and the sense read off the master's
edge nodes; the contract in one Doxygen block at the declaration. Both auditors passed the C++ and
failed the doc hunk in `thin_shell_facet_orientation.md`: "per-block numbering gives a uniform
master" is an input assumption (lower element id, contiguous ids per block), a pre-existing
paragraph promised a BFS consistency that seedless components never get in either dimension, and my
insertion had shifted the doc's Stage 2 line range onto the wrong function. All fixed, plus Grok's
comment trims, a guard on the second normalization, and Christian's correction that the lower-id
rule applies only to same-type facets. Syntax-checked with the tree's flags; not built. Row retagged
`[RUN][P]`: owed are `make check` on the rebuilt tree and the tape-6 probe rerun on the patched binary.

## Same session: DR-133 struck

`make check` green on both backends (Christian). The tape-6 probe reran on the patched debug Blaze
binary with assertions on: no new check fired, per-tape currents identical to the pre-patch table,
tape 6 sign held, far field unchanged. Row struck and archived; the 3-D sibling's two hygiene twins
are recorded in the archived row, not fixed. Archive 147 rows, live register 6 (4 `[P]`, 2 `[W]`).

## Not done

One C++ change, the DR-133 hygiene diff in `cl_CutFactory.cpp`, built and `make check`-green on both backends by Christian, probe-verified here. The CMake change was configure-gated in scratch trees, not built; `make check` after a reconfigure is owed and Christian runs builds. The two DR-155 gates named above stay unrun by ruling; if
`USE_BELFEM_OPENMP` is ever defaulted ON, they come back as the precondition.
