# Devlog 2026-08-11 — SPFA Solver Test Port (DR-23, T8 part one)

**Date:** 2026-08-11
**Topic:** The difference-constraint solver behind the thin-cut rectifier was production code
with zero checked-in coverage; its test harness existed only in a swept scratchpad. Ported.
**Module:** `tests/math`, covering `src/math/graph/fn_Graph_spfa`
**AIs involved:** Claude
**Claude Confidence:** high that the tests are correct as written (they verify answers
independently, not return values); **the tests have not been executed** — `USE_TEST` is OFF in
the shared tree
**Verification:** `mpicxx -fsyntax-only` under the project's own `-Wall -Werror
-pedantic-errors`. **Not run.**
**Plan:** `todo/thin_cut_nonunit_rectification_implementation.md` T8 · **Register:** DR-23

## Why this was worth doing first

`spfa_difference_constraints` decides whether a set of thin-cut coefficients can be rectified
to unit, and when it cannot, produces the negative cycle that names the throat. It is called
from `Cohomology::clean_spfa()` on every non-unit generator. It had **no** checked-in test.
The harness that once exercised it — 18 checks, per `dl20260715_spfa_clean_implementation.md`
— lived in a session scratchpad and was swept.

It is also the easy half to port: a free function over plain `Cell`s, no mesh, no MPI.

## What landed

`tests/math/test_GraphSpfa.cpp`, 16 tests, wired into the existing `math` suite (`fast`
label). `graph` is already in `BELFEM_LIBLIST_BASE`, so no link changes were needed.

**The tests verify the answer, not the verdict.** This is the point of the file:

- a *feasible* verdict is accepted only after re-testing **every** constraint against the
  returned theta, so a solver that returns `true` with a wrong solution cannot pass;
- an *infeasible* verdict is accepted only after verifying the returned certificate is a
  genuine closed walk whose weights sum **strictly** negative.

Coverage: no arcs; isolated vertices staying at zero; a simple chain; a tight unit chain
(the production shape, both directions pinned); disjoint components; self loops both signs;
zero-weight and balanced cycles; a negative triangle; the greedy-hang 3-cycle at c = 2; a
negative cycle buried in a larger feasible graph; 200 randomized feasible-by-construction
instances; 100 with an injected negative cycle; a 2000-vertex chain; a 500-vertex chain closed
by a negative return arc.

**Two of these are regression guards, not coverage padding**, and both correspond to defects
the original scratch harness caught:

- **A zero-weight cycle must read FEASIBLE.** It only forces the potentials around the cycle
  to be equal. An early solver accepted a zero-weight predecessor cycle as proof of
  infeasibility; the audit caught it.
- **Large feasible instances must survive the spurious `cnt > n` trigger.** The classic SPFA
  update-count test can fire *before* the predecessor graph holds a cycle. An early version
  treated that as an internal error and aborted on a feasible system — the scratch test caught
  it before it reached the repo, and the randomized feasible set is what would catch it again.

The randomized instances are feasible **by construction** — potentials are drawn first, then
each arc weight is set no smaller than the difference it must admit — so an infeasible verdict
there is a false negative in the solver, not a fixture artefact. The PRNG is a fixed-seed
xorshift rather than `std::rand`, whose sequence is implementation-defined, so a failure is
reproducible on any machine.

## What is still owed on DR-23

**The Cohomology layer has no coverage at all.** `clean_spfa()`,
`rectify_greedy_sweeps()`, `fire_node_coboundary()` and `remove_cut_pockets()` take mesh
entities and need a fixture — exactly what `dl20260715` predicted would be the hard part. The
solver underneath them is now guarded; the code that calls it is not.

**T9** — the single-layer corc periodic regression — is a run, not a unit test.

So DR-23 stays open and blocking-1.0, with its scope now sharpened rather than closed.

## Handoff

`make check-fast` with `USE_TEST=ON` is the gate, and it is a real one: these tests have never
executed. Alongside them, the spline suite (DR-49) and `Hex8TbUnitCirculation` (DR-47) also run
for the first time in that configuration. A red result in any of the three is information about
the code under test, not necessarily a regression from this work.
