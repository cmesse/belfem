# DR-23 Stage 1: Cohomology-Layer Test Suite on a Programmatic Annulus

**Date:** 2026-08-13
**Purpose:** Land the layer half of DR-23/T8 — checked-in coverage for the
`Cohomology` entry points (`clean_spfa`, rectification, pocket census) through the
production call sequence.
**Module:** `tests/homology` (new), `todo/`

## What landed

- **`tests/homology/`** (new test module, label `fast`, wired into
  `tests/CMakeLists.txt`): `test_Cohomology.cpp` plus main and CMake.
- The fixture builds a **programmatic annulus** — one ring of N=8 sectors between
  two radii, as QUAD4 or as 2N TRI3, built node-by-node like the `tests/mesh`
  tensor meshes and `TS_TestStack`; **no mesh files needed**. Christian's
  question — where checked-in test meshes should live — was ruled the same
  session: **`tests/meshes/`**, colocated with the tests rather than under
  `share/` (which carries user-facing runtime data). The directory now exists
  with a README carrying the conventions (`.geo` source of truth checked in
  next to the generated `.msh`; prefer programmatic meshes when possible).
- Six cases: {TRI3, QUAD4} × {Pellikka, CCR, PellikkaGeneralized} coreduction.
  Each drives the production `CutFactory::compute_cohomologies` sequence with
  `mSuggestHomologies` semantics: flag the region, `SimplicialComplex`,
  **coreduce only**, chain-side `Homology` SNF, `Cohomology` constructor (which
  runs `generatorsOfCohomology` + `clean()`, i.e. `clean_spfa`),
  `updatekGeneratorsFromHomology`, `clean()` again.
- Assertions: betti numbers of the annulus (H_0/H_1/H_2 = 1/1/0 on the chain
  side, one cohomology generator on the cochain side), the rectification
  contract (every surviving coefficient ±1), and the **winding pairing**: the
  generator paired with the inner-ring cycle must give ±1 — rectification may
  change the representative, never that pairing.

**Evidence:** compiled `-Wall -Werror -pedantic-errors` with the tree's flags and
run as a standalone probe against the prebuilt `.a`s — **6/6 green in 45 ms**.
The rectifier and pocket census *execute* (their log lines appear) but find
nothing to repair on this mesh — the generators come out unit-coefficient
already, so the assertions pin the contract, not a repair; actual repair
behaviour stays covered by the 16 solver-tier tests. The gtest suite itself has
not run (`USE_TEST=OFF` shared tree); the gate remains `make check-fast` with
`USE_TEST=ON`, now covering both DR-23 halves.

## The debugging story, kept because both lessons are load-bearing

The first fixture draft produced H^1 = 24 generators (one per edge), then a
99%-CPU hang inside `reduce_complexCCR`. Both symptoms briefly looked like
framework defects. The degeneracy was caused outright by **one fixture bug**,
and that same bug is what *exposed* the latent hang — the hang itself lives in
real (dead-branch) framework code, see DR-71 below. The fixture bug: the mesh
was built with `aComputeConnectivities = false` (copied from `TS_TestStack`,
which hand-builds its edges and needs it). With connectivity off, `finalize_edges`
skips `ConnectivityCalculator::connect_edges_to_elements`
(`cl_Mesh.cpp:975-982`), every edge cochain gets an empty coboundary, the
coreduction cannot collapse a single copair, and the cohomology SNF sees zero
maps — one "generator" per edge. The chain side stays intact (elements know
their edges), which is why `Homology` returned correct betti numbers throughout
and made the discrepancy diagnosable. Lesson one, now in the test header:
**fixture meshes for topology work must keep `aComputeConnectivities = true`.**

Lesson two: the fixture initially mirrored the `!mSuggestHomologies` branch of
`CutFactory::compute_cohomologies` (reduce **and** coreduce), because that reads
like the straightforward pipeline. That branch is dead in production — the flag
is hardcoded `true` (`cl_CutFactory.hpp:86`) — and it is broken when forced,
in two independent ways now recorded as **DR-71**: reduce+coreduce together
zero the SNF input, and `pGeneralizedCombine` carries an unbounded iterator
walk (`cl_SimplicialComplex.cpp:873-878`, `++it2` with no `end()` bound and an
unreachable guard below — the observed hang; its cocombine sibling at
`:1139-1148` is properly bounded). The production sequence is coreduce-only
with the homology-suggested generator update, and that is what the tests pin.

## Stage 2 (same day, Christian's go-ahead): the periodic quotient

Reference per Christian: the sidecoating / tapestack3d infinite-tape setup —
a periodic mesh whose end planes are gmsh `Periodic Surface` copies, with the
cohomology structure named in `dl20260812_tapestack3d_setup_and_terminal_guard.md`:
one generator encircling the conductor **plus one free axial generator from the
periodicity**, H^1 rank 2. Two 3D miniatures of that in the new
`test_CohomologyPeriodic.cpp`:

- **PeriodicBar** — 1×1×N HEX8 bar, z-fold: solid torus, betti (1,1,0); the
  single generator is the free axial one, paired ±1 with the axial cycle
  through the fold.
- **PeriodicBarHole** — 3×3×N with the center column removed, z-fold:
  (annulus × S¹), the "periodic cylinder minus one wire" skeleton — betti_1 = 2,
  betti_2 = 1; generators are basis-dependent, so the invariant assertion is
  the 2×2 pairing matrix against the axial and encircling reference cycles
  being **unimodular**.

The periodicity is **hand-wired** — distinct slave-plane nodes, symmetric
`set_periodic()` links, id-matched edge/face pair lists — bypassing the
geometric `PeriodicityFactory` on purpose: the subject is the fold inside
`SimplicialComplex`/`clean_spfa`, not the plane matcher. One production
convention had to be replicated for that to work: slave-plane boundary faces
must be **re-based from master to slave role**
(`PeriodicityFactory::fix_face_slaves` is private; `create_complex:414` reads
`face->slave()` for the fold orientation and null-derefs on the FaceFactory
default — the first probe run's segfault). And one fixture rule follows from
the quotient itself: a cycle through the fold closes only in the quotient, so
the pairing helper walks explicit edge steps with no implicit wrap.

**One API contract learned and recorded in the test:**
`updatekGeneratorsFromHomology`'s 3D branch expects the suggested homology as
terminal **in/out pairs** (`[in0, out0, …]`, one condition per pair — the shape
`suggest_Homology` produces from a maxwell deck). Raw SNF homology generators
do not satisfy that: fed anyway, the update builds one condition generator
(g_ax − g_hole) plus one free cut (g_ax + g_hole) — pairing determinant 2,
exactly as its own Step-4 comment predicts. The periodic tests therefore
assert on the **constructor** generators (the SNF + `clean()` path, where the
periodic `clean_spfa` handling runs); the update path stays covered by the 2D
annulus suite, whose 2D contract is one-generator-per-condition.

**Evidence:** full suite 8/8 green in 7.2 s (probe against the prebuilt
`.a`s, `-Wall -Werror -pedantic-errors`); the hole case dominates at ~5 s
(SNF on 184 edges), inside the fast-label criterion.

## What DR-23 still owes

- ~~**Stage 2** — the periodic-quotient case~~ — **landed the same day, see the
  Stage 2 section above** (3D periodic bar fixtures rather than the 2D strip
  first proposed, matching the production infinite-tape topology).
- A layer-tier case that **forces a repair**: all fixture generators come out
  unit-coefficient already, so the SPFA path runs but never rectifies here;
  repair behaviour stays covered by the 16 solver-tier tests.
- **The suite gate**: `make check-fast` with `USE_TEST=ON` has still never run —
  everything above is probe-tier evidence.
- **T9** — the corc periodic regression run (Christian).

Register: DR-23 row updated (layer half landed, gate now both halves); DR-71
filed. Plan file `thin_cut_nonunit_rectification_implementation.md` T8 → [◐]
with the same record.
