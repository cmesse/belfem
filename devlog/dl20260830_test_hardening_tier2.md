# Test Hardening: the Tier 2 Round, and What the Devlogs Said About Our Own Tests

**Date:** 2026-08-30
**Purpose:** Record the jury round that picked BELFEM's three costliest test gaps from its own
539-incident history, the plan that came out of it, and the first three steps landing in code
**Module:** `config/scripts`, `tests/comm`, cross-cutting

## The task

Christian asked for a detective round: sweep the devlogs, find what breaks most often, and use it to
name gaps in the test suite. Three AIs, then plan+audit → code+audit. Two books on scientific
software engineering in `tmp/` as the source of selection criteria.

The framing worth keeping: **"broken most often" is not the same question as "test gap."** A
subsystem can fail repeatedly and still be untestable at a price we would pay. The round only worked
because that second filter was applied explicitly rather than assumed.

## The evidence already existed

`doc/lessons_learned.md`'s evidence file already held **539 incidents mined from 298 devlogs**, one
row per incident, each citing its source devlog, with a `fix class` and a `how discovered` column.
Re-sweeping the devlogs would have measured narrative attention — a subsystem that failed loudly
three times and earned one long devlog looks colder than one that produced five short notes. The
catalog is already normalized.

Two numbers framed everything after:

- **280 of 539 incidents were code defects. 38 were closed by adding a test — 7%.**
- Discovery: AI/static audit 356, production solver runs 158, hand probes 90, the suite ~38.

Both auditors correctly noted the second number measures recorded *remediation class*, not
detection. Accepted. It does not survive being halved either.

Coverage correlated: `src/mesh` 52,318 source lines against 915 test lines, `src/fem` 90,730 against
5,719 — the two thinnest-covered modules carry the largest concentrations of code defects, while the
best-covered (`containers`, `math`, `circuit`) are near-silent in the incident record.

## What the round produced

Three gaps, both auditors converging on the first and third:

- **G1** — the periodic/thin-shell pipeline is never driven through its production factory. The
  suite does not merely miss it, it *documents* missing it (`tests/mesh/test_MeshEntities.cpp:11`),
  and the one periodic test bypasses `PeriodicityFactory` **on purpose**
  (`tests/homology/test_CohomologyPeriodic.cpp:26-34`). False comfort, worse than a zero.
- **G2** — the suite is **serial by construction**. `Add_Test.cmake:29` registered every test as a
  bare command; `mpirun`/`MPIEXEC` appeared in no CMake file in the tree. Meanwhile
  `tests/comm/test_CommMPI.cpp` — 19 KB, written three days earlier, its own `main()`, 29 rank
  guards — was referenced by no `CMakeLists.txt` at all. A Tier 1/Tier 2 vocabulary already existed
  in comments across three files and had never been built.
- **G3** — the reject/warm-restart state contract. Layer 1 carries two failure signatures for it
  (L-12, L-19) and the suite implements neither; `load_memdump`, `reset_timestep` and `seed_dof`
  have no call sites under `tests/`.

Christian's scoping ruling: build the parallel capability, but run multi-rank only where the
*subject* is communication. Verified against the tree — the `fem` matches are `comm_size() != 1`
serial skips, and `sparse`'s are comments and enum names.

## What the process caught, which is the actual story

**Sixteen defects were recorded, four of them mine, and none reached a build.**

- **D1** — the sealed pre-registration's second gap rested on INC-533, whose DR-73 had been
  **retracted the same day it was filed**: the 8-rank dof figure was magnetic + thermal combined.
  Both auditors caught it independently. The retraction record names the exact mechanism — a
  secondhand number placed in a brief's "established facts" section — and this round reproduced it.
- **D2** — 94% of physics test code sits behind `USE_GASMODELS`, default OFF. Default-build physics
  coverage is 0.7%, not the 10.5% first reported.
- **D5/D6/D12** — three planned assertions would have failed against *correct* code or passed with
  their cited defect live. R5 asserted a post-restore Δt equal to a cold run, but restart
  deliberately caps that step and the cap *is* the INC-351 fix; the assertion tested for the bug.
- **D9** — the Tier 2 design **failed open.** A `GTEST_SKIP` returns 0 from `RUN_ALL_TESTS()` and
  ctest reads the exit code, so a launcher silently degrading to one rank would have produced a
  green line: 29 skips and one trivial pass. The campaign's own headline failure mode, reproduced
  inside its fix.
- **D11** — six assertions across R3/R4/R5 would have passed with their defect live, forcing G1
  apart into three fixtures.
- **D15** — the second revision introduced four consistency defects of its own, caught by a third
  read. Root cause worth keeping: *a rewrite reconciles the text it replaces, not the text that
  referenced what it replaced.*
- **D16** — the code audit found the landed `set_throw_on_error( true )` carried a comment false of
  its own file, and would convert a release `MPI_Abort` into a caught throw.

## The disagreement worth recording

On D16 the two auditors prescribed **different remedies**. Codex required forcing
`set_throw_on_error( false )`, reasoning that a rank-local throw strands peers in a collective.

That is precisely the argument in `feedback_debug_must_throw_in_parallel` — implemented, audited
twice by both vendors, and **reverted by Christian on 2026-08-30, the same day as this work**. It
was refused. Grok, auditing independently, reached the remedy that had already been applied: delete
the call and let `gThrowOnError = BELFEM_ASSERTIONS_ACTIVE` stand, so debug throws and release
aborts.

**An auditor finding a real defect does not make its prescribed fix right.** The defect was real and
both found it; one of the two remedies would have re-litigated a ruling hours old.

## What landed

R1, R1b, and most of R2. `Add_Test.cmake` gained an opt-in `TESTRANKS` path that registers
`${TESTNAME}_np<N>` under the launcher and **suppresses the serial twin**, with a non-skippable rank
sentinel, `MPI_Allreduce( MAX )` over the test result, `TIMEOUT`, `PROCESSORS`, `--oversubscribe`,
per-rank label/environment propagation, and per-suite variable reset. `test_CommMPI.cpp` was
rewired into the tree's first Tier 2 binary and gained three `share`/`receive` tests — that function
had no test call site anywhere despite being the mandated path for large payloads.

**VERIFIED — two executable gates ran on Christian's build.**

`make check`: **17/17 passed**, with `commmpi_np2` and `commmpi_np4` registered and running, and a
label summary of `mpi = 2 tests` / `fast = 10 tests` — the D13 leak closed, no multi-rank work in
`check-fast`, and the 13 pre-existing suites unchanged. A green `commmpi_np4` is itself the
evidence, not merely the absence of failure: the sentinel cannot be skipped, so it means
`comm_size()` really was 4 and CTest's `ENVIRONMENT` forwarded `BELFEM_TESTRANKS` through `mpirun`.

Then the failure-propagation probe, which produced exactly its predicted **two-sided** result:
**`commmpi_np4` failed while `commmpi_np2` passed.** Rank 3 was the only failing rank, so rank 0's
own `RUN_ALL_TESTS()` returned 0 and the red verdict could only have arrived through the
`MPI_Allreduce`. That is the difference between the probe proving the reduction works and merely
proving that some failure somewhere fails the run — and it is why the probe failed one rank rather
than all of them. O4 is now settled by measurement; the ~60%-confidence assumption about Open MPI
folding non-root exit statuses is no longer load-bearing anywhere. Probe removed the same session.

The mechanism is at the top rung of the evidence ladder. That says nothing about the fixtures still
to come.

## Residue

- `todo/test_hardening_campaign.md` — R3a/R3b/R3c, R3d, R4, R4b, R5, R6, R7 open
- **O5 open and gating R4:** "force a z-slice partition" has no known mechanism —
  `mesh::Partitioner` is METIS-only, with no injection point for a prescribed partition
- `CLAUDE.md` and `doc/coding_philosophy.md` corrected: the nightly GitLab CI exists, but the
  pipeline definition is not in the tree, there are still no sanitizers, and a green nightly
  measures the suite rather than the code
- The two assert reactions were split rather than treated alike, on Christian's ruling: the syslog
  hook is for *a calculation that crashed after hours*, whose terminal output is gone, and a test
  never has that problem. `set_syslog_on_error( false )` restored — not by analogy to the other test
  mains, but because the nightly would become the highest-frequency writer to the `belfem` syslog
  identity and would dilute `journalctl -t belfem` for the production crash it exists to diagnose.
  `throw_on_error` stays unset so debug throws and release aborts
