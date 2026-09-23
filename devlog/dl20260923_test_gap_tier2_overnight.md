# Test-gap Tier 2: the overnight run

**Date:** 2026-09-23
**Purpose:** Session record of the unattended Tier 2 run of the test-gap campaign
(`todo/test_gap_tier2_execution.md`), under the standing approvals of the Tier 1 run
**Module:** `tests/` (core, containers, linalg, math, sparse, physics, fem, circuit, io, mesh,
homology); `src/containers`, `src/fem/kernel`, `src/sparse`, `src/io`, `src/comm`, `src/mesh` for
jury-confirmed local fixes
**AIs involved:** Claude Fable 5.1 (B1 gate and slice-1 reconciliation, before the reload), Claude
Opus 5.5 (the loop from the reload on, and the implementation agents), Codex (`gpt-6-astra`/high on
kernel and mesh rows, `gpt-5.6-terra`/high otherwise, `gpt-5.6-luna`/medium for prose sweeps), Grok
`grok-4.7`/high

## Summary

Christian started the run at ~22:00 ("make this an overnight run"). It worked off the Tier 2 batches
in the order of the tracker's handoff, one batch in flight at a time, with Opus agents writing the
test-only rows and the loop gating, mutating and juring every line. **Landed on `devel`:** B1's five
test rows, slice 1 (DR-185 and the map message), B6 (REBCO table, λ pins, magnetization), B3 (sparse
wrapper contracts plus the relative-tolerance and max-iterations setters), B5a and B5b (fem, circuit,
physics contract rows), and B2's T-06 and T-07. **Held:** T-04 (HDF5 existence helpers) with the
evidence below, B4 and B7 (peer work in progress in `src/mesh` and `tests/comm`), and the
ruling-gated rows. Every commit passed `make check` 17/17 in both trees before it landed.

## What changed in `src/`

| Commit | Row | Change |
|---|---|---|
| `25dbc53` | DR-185 | `Kernel` deletes its copy and move members; four traits pinned by `static_assert` for Kernel, Solver, ElementMapper, ElectricalCircuit |
| `a32ad4e` | T-05 / N-20a | `Map` / `OrderedMap::operator()` name the missing key in release builds too; "unknown" for a key type without a formatter |
| `0d6b726` | T-09 / T-20 | `set_relative_tolerance` refuses non-positive / non-finite values, and the deck key goes through it; new `set_max_iterations( int )`; input reference and schema updated |
| `fe52394` | T-18 by-catch | `MUMPS::warning_message` reported an unknown bit only when no known bit was set; now always |
| `1b92533` | T-07 | `Ascii` resolves a relative path against the process cwd, not `$PWD` (ruling O1); the banner's run directory too |
| `9cb58e3` | T-06 | a `.bfm` stamped by another commit loads with a warning naming both hashes (ruling O9) |
| `90f8901` | T-15 | the `gRhoMin < gRhoMax` precondition stated at the declaration |

## Findings the measurements produced

- **T-19:** the degree-8 polyfit bound the implementing agent chose (1e-5) failed on Armadillo at
  1.2e-5. Measured worst coefficient: Armadillo/MKL 1.2e-5, Blaze/MKL 5.9e-6, reference LAPACK
  6.8e-6, OpenBLAS 3.6e-6; the bound is 1e-4. A normal-equations solve misses by 7.8e+2 (the comment
  had said 7e+1).
- **T-29 → DR-189:** the REBCO low-field bridge engages at every angle only from 76 K up on
  `sp-ap.hdf5`; between 56 and 72 K the monotone limiter leaves slope kinks of 0.2–1, and from 52 K
  down the guard falls back to the flat clamp. The bridge comment's "the measurements satisfy
  J0 ≥ jc( Bmin, θ )" is false there. Two `DISABLED_` pins; the enabled cases assert the C1 join and
  the angle-free J0 at 76, 80, 84 K. Monotonicity in B holds on a sampled grid to roundoff
  (1.6e-15 jc, 8.2e-15 n), refuting a jury worry of a 1e-10 rise.
- **T-30:** the model sits at −1.7 % (Ni, vs 90.7) and −1.8 % (Fe, vs 83.2) at 300 K; the pins stay
  disabled because neither reference was read from the cited table.
- **T-04 held:** a mechanical classification of every literal existence probe against how its
  label is written found `group_exists` used on datasets (`abstract`, `pinned`, `orphaned` in
  `cl_Mesh_BfmFile.cpp`, `extra` in `cl_Spline.cpp`) and 16 BfmFile probes on names written as both.
  A type-checked `group_exists` would silently drop that data on every `.bfm` reload; the plan
  inventory had missed it in two rounds.
- **By-catch, recorded not fixed:** `BfmFile::save` throws on a mesh without facets; a Maxwell
  `dirichlet` section refuses any dimensioned amplitude; the refusal paths of three factories and of
  `JcFunctionDatabase` leak what they built before the throw (DR-175 class).

## How the loop ran

- **Mutations on scratch copies only.** A harness (session scratchpad, `mut/run.sh`) copies one
  source file, mutates it, and recompiles the real test sources against the built library with the
  mutant first on the include path (or compiled in for a `.cpp`). The shared tree was never mutated,
  so juries could read while mutations ran. Two lessons: a header that another header includes
  first is masked by its include guard (force it with `-include`), and a define carrying quotes
  needs `eval` on the compile line.
- **Staging beside peer WIP.** Besides `stage_line.py` (added lines), a new
  `tmp/ai_exchange/stage_replace.py` applies an exact replacement to the index copy of a file, so a
  one-line edit in a file carrying a peer's uncommitted hunks is committed without them.
- **Agents.** Five Opus implementation agents (B1 remainder was gated by the loop; B3, B5a, B5b, B2
  tests); every line they wrote was read before its jury. Their reports were verified against the
  tree; two plan errors they found (a protected decoder, an equivalent mutant) were real.
- **Jury rounds:** slice 1 plan + code, B6, B3 plan + code, B5a (Grok rerun alone after it hit
  100 turns), B5b, B2 plan (two rounds on T-04) + T-06/T-07 code, B2 tests.

## Open for Christian

All in `todo/test_gap_tier2_execution.md` §3g, "For the morning quiz": DR-189, the T-30 reference
values, the metis setter, the duplicate-column row, DR-175 widening, the ElementMapper Newton
policy, DR-187's narrowed ruling (the thermal default arm carries by-design exclusions), T-04's
options, and the two by-catch register rows.
