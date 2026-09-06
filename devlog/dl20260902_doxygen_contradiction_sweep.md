# Devlog 2026-09-02 — Doxygen Contradiction Sweep (comments and pages vs the code)

**Date:** 2026-09-02 (evening session, closed 2026-09-03 00:xx PDT)
**Topic:** Pre-release sweep of the documentation Doxygen renders — header/source doc comments and
the markdown pages — for statements the code contradicts; documentation and comments edited only
**AIs involved:** Claude (coordination, verification, core/containers/comm partition), 16 Claude
subagents (read-only finders) + 12 Claude subagents (verify-and-apply), Codex `gpt-5.6-terra`/medium
(independent finder, comm+containers+core+io), Grok `grok-4.6`/high (independent finder, the
`input.conf` contract)
**Claude Confidence:** high on every applied edit — each cited comment and code line was re-opened in
the working tree before the edit landed; medium on the completeness of coverage (three partitions
were audited but not applied, see Open Questions)
**Literature References:** N/A (documentation currency; one auditor attribution corrected to
NASA RP-1311 for the gas mixing rules, per the code's own citation)
**Verification:** **reviewed, not verified.** No build, no `make check`. Executable gates that ran:
`scripts/check_doc_claims.py` 37/37; `doc/input_schema.yaml` parses; a Doxygen 1.9.1 probe into the
scratchpad before and after (see below). No source token changed: a `git diff -U0` filter over every
`.hpp/.cpp/.f90` hunk shows comment lines only, plus one code line that is Christian's own
uncommitted work (`cl_MaxwellFactory.cpp` `std::abs(...)`, mtime 21:27, matching his 21:03 commit).

## Summary

Doxygen renders every `src/` doc comment plus `doc/`, `src/*/doc` and `examples/`. The tree was
partitioned into 16 read-only audits with a shared brief (contradictions only; quoted code evidence
required; missing docs, style and code bugs out of scope). Every finding was then re-opened at its
source by a second agent (or by Claude) before a comment was changed. **331 contradictions were
corrected** across 167 files; **about 245 audited findings remain unapplied** because the API session
limit killed the last applier wave — the reports are preserved (see Files) and the residue is
tracked in `todo/doxygen_contradiction_sweep.md`.

The dominant defect classes, in order of frequency: copy-paste drift between near-identical
overloads and element classes (a `Matrix` overload documented as "vector", `@param` names from the
sibling class, a HEX node count on a PENTA); comments that describe a design the code left behind
(a `Cell<bool>` memory comparison that `std::vector<bool>` refutes, "Parallel now aborts" for an
error reaction that is build-type-selected, `set_bh_curve()` taught where `load_bh_curve()` is the
call that reroutes μ); wrong tier vocabulary ("asserts" for an always-on `BELFEM_ERROR` and the
reverse); and defaults stated in prose that the code's initialisers contradict (spline solve on
every rank, gas tables born in SPLINE mode, the timestep controller's PID default).

## Key Findings

**The input contract had six defects Grok found and the tree confirmed.** The `custom { }`
material subsection is refused by name since 2026-08-29 but appeared in neither contract file (now a
`refused_subsections` entry in the schema and a sentence in the prose); the `RRR` table cell said
"discarded" where `check_unused_input` refuses it on non-`PureMetal` builtins; `critical temperature`
on a `curve` section was documented as ignored where it is fatal; the schema's
`required_unless_section` for `linear` OR-combined the two field sections while each factory falls
back to `section("linear")` per field; a schema note asserted something false about the prose; and
`KrylovMethod::AUTO`'s enumerator comment described neither backend correctly.

**Codex's 17 on comm/containers/core/io were all real.** Seven "vector" comments on matrix and raw
array overloads in `commtools.hpp`; the `DynamicBitset` "8x saving over `Cell<bool>`" claim
(`Cell<T>` wraps `std::vector<T>`, so `Cell<bool>` is already bit-packed — also corrected in
`CLAUDE.md` and `coding_philosophy.md`, `check_doc_claims.py` still 37/37); `random_seed()`'s serial
branch draws a seed and never calls `std::srand`; `gTbulk` is NaN until the executable or deck sets
it; `InputFile::section(string)` looks up by type, not label; `XML::next_sibling_of_same_name`
advances rather than counts; `create_group`/`select_group` briefs said "dataset"; the HDF5
existence probes test link presence only.

**The 2D thin-shell, mesh and kernel partitions carried the most rot** (mesh 72, kernel 45, fem
physics 39 corrections). Representative: `Element::compute_edge_directions` cited from the mesh
namespace where the method lives on `fem::Element`; the periodicity guide had the k-d tree on the
wrong side and the wrong entity; `cl_FEM_Tmatrix.hpp` documented `B = T·A` for a `Tᵀ·A` product;
`nonlinear_controller_theory.md` taught the legacy `sqrt(target/used)` controller while the code
default is the log-space PID that `bdf_timestepping_theory.md` already describes;
`Kernel::material()` documented as creating a missing material where it hard-errors; the kernel
README listed an enumerator from a different enum.

**Two refutations went the other way.** The linalg auditor flagged a `parpacktools.f90` comment
saying `dsaupd` returns `INFO = 3`; the ARPACK source (`dsaupd.f:254`) documents exactly that, so the
three sites claiming the symmetric driver "has no info = 3" were the stale ones and were corrected
— and `check_saupd` has no case for it (code-side, below). The mesh-entity auditor claimed
`compute_edge_directions` does not exist; it does, on `fem::Element`.

**Doxygen mechanics.** A 1.9.1 probe (`GENERATE_HTML=NO`, output into the scratchpad) found five
warning sites outside the quarantined `lessons_learned_evidence.md` that
`cmake/report_doxygen_warnings.cmake` states do not exist: an `<aLabel>` placeholder read as an HTML
tag, an `@param T` on a function without `T` in `examples/block3d/src/matlib.cpp`, `::Boundary`
link requests in the homology handoff, an out-of-tree markdown link in the materials README, and
four undocumented `ThinShellFactory` constructor parameters. All five fixed; the re-probe shows
**0 warnings outside the evidence file and the layout-version noise** (the 32 `DoxygenLayout.xml`
lines are the raw template against 1.9.1 — `make doc` normalises the layout, the probe did not).

**Shared-checkout incident.** Every file this session had edited before 20:19:07 carried that exact
mtime afterwards and none of the edits: a bulk `checkout`/`reset` in the shared tree (Christian
committed `cfaaa01b` at 20:24) reverted roughly an hour of work — the five Doxygen fixes, the Codex
batch, the bitset claims. Detected only because the closing Doxygen re-probe reported the *same*
five sites; the io applier had also reported my io edits "not present" and I had misread that as
equivalent rewording. All lost edits were re-applied from the recorded strings. Lesson: in a shared
checkout, `git diff --stat` on your own files before the closing summary, not just `git status`.

## Changes Made / Proposed

Applied (documentation and comments only), by partition, count of findings corrected:

| partition | corrected | source |
|---|---|---|
| comm / containers / core (incl. two usage guides) | 20 + 17 | Claude subagent + Codex |
| io / circuit / visualizer / executables | 33 | subagent |
| linalg / sparse (+ 3 ARPACK `INFO = 3` sites) | 36 + 3 | subagent + Claude |
| mesh core/IO, entities, algorithms (+ 5 mesh doc pages) | 40 + 16 + 16 | subagents |
| fem/kernel: controller group, dof/calculator group | 23 + 22 | subagents |
| fem interpolation / iwg / maxwell / postproc / thermal | 39 | subagent |
| physics/gasmodels | 23 | subagent |
| numerics / database / gastables | 20 of 23 | subagent (interrupted) |
| physics/materials | 6 of 38 | subagent (interrupted) |
| top-level pages (`README.md`, `doc/`, `examples/`) | 6 of 22 | subagent (interrupted) |
| `input.conf` contract (both files) + `en_SolverEnums.hpp` | 6 | Grok |
| Doxygen warning sites; `todo/debt_register.md` `[P]` count 8→7 | 5 + 1 | Claude |

Not applied — audited reports preserved under `tmp/ai_exchange/doxygen_sweep_reports/`:
`homology_math.md` (42; five in the closed cohomology core, report-only for Gregory),
`module_docs_nonfem.md` (116), `module_docs_fem.md` (36), plus the remainders of materials (32),
top-level pages (16) and numerics (3). Each carries quoted evidence and proposed text; the apply
brief (`apply_brief.md`, same directory) is the procedure.

**Code-side observations relayed, none acted on** (comment-only mandate; a code decision each):

- `src/core/constants.hpp` `y2k = 2462502.5` is the JD of 2030-01-01; the comment said 2000. Comment
  now states the value; no consumer in `src/`. Which side was meant is Christian's.
- `src/core/random.hpp` serial `random_seed()` never seeds `std::rand()`.
- `src/sparse/arpacktools.cpp` `check_saupd` has no `case 3`; `dsaupd` documents it.
- `src/fem/interpolation/nedelec/cl_EF_PENTA6TS.cpp:36` `mSumW = 8.0` (HEX value) for a wedge whose
  weights sum to 1; `tests/fem/test_InterfaceOrientation.cpp:201` recovers an area from
  `abs_det_J() * sum_w()`.
- `Controller::check_thermal_diagnostics()` consumes its latch and emits neither warning the magnetic
  twin emits; docs now say "magnetic field only".
- Error-message string literals that name the wrong thing: "OldMaxwellFactory" (`cl_IwgFactory.cpp:46`,
  `cl_FEM_Kernel.cpp:561`), `compute_volume` inside `compute_surface`
  (`fn_Mesh_compute_surface.cpp:30`), "tape roller" (`cl_FEM_Element.cpp:1424`), the `gas` and
  `gastable` help texts (a `-v` flag with no parser branch; program name and default verbosity).
- `cl_GM_Helmholtz.hpp:417` declares `check_parent()` that is never defined.
- Mesh: `faces_are_finalized()` returns `mIsFinalized`; `VtkWriter::write_node_fields()` increments
  `mNumberOfNodes` twice; `GmshReader` `delete_unused_nodes_and_elements()` commented out and
  `read_elements_v41()` never sets physical tags; `element_type_from_gmsh` is a plain `static_cast`
  (gmsh 32 is `MSH_TET_22`, not QUAD16 — `GmshDefines.h:114,118`); const/non-const `max_element_id()`
  scan different sets; `unfinalize()` resets connectivities twice; `cl_Mesh_BfmFile.hpp` forward-declares
  `belfem::ProtoMesh` in the wrong namespace; `Mesh::memory()` undercounts owned vertices.
- Periodicity k-d prune compares a squared plane offset with an unsquared best distance;
  `OrderConverter::check_input_mesh()` reads `mMesh` before the constructor assigns it.
- Dof manager: `reorder_dofs()` gates on METIS/SCOTCH but only calls `symrcm`; `SolverData::jacobian()`
  returns `mSystemMatrix`; `run_shift_invert()` error text names `run_parpack()`; `link(Facet*)`
  hard-codes a TET slave index; `set_field_index` `uint` vs `index_t` mismatch;
  `dof_manager_usage_guide.md:2133` still names `SolverReordering::RCM`, which no enum has.
- `fn_compute_permutation.hpp` header guard still says `RCM` while the routine runs METIS.
- `stringtools.cpp` `muV` scale `1e-3`, `MN` scale `1e3`, unreachable second `cd` branch;
  `to_real("")` returns 0 where the header promises NaN; `Timer` truncates at ~49.7 d;
  `DynamicBitset::to_int_fail` returns `BELFEM_UINT_MAX`, not `gNoIndex` under `BELFEM_INT64`;
  `assert.cpp` reaches `comm_abort` with no declaration in a non-MPI build.

## Second round — auditing the sweep itself (2026-09-03)

Christian asked for a deep Codex + Grok round on what the first pass might have missed. Depth per
§9.1 for a code-diff audit, round 1: Codex `gpt-5.6-terra`/high, Grok `grok-4.6`/high. Predictions
and a falsifier were pre-registered in the exchange before dispatch. The tree was byte-identical to
the patch snapshot before and after both runs (no auditor write, and no edit of mine could produce a
phantom finding).

**Codex, Part A — the sweep's own regressions: zero.** All 450 hunks across 177 paths, checked for a
new false statement, a deleted live caveat, broken Doxygen markup, a cross-file contradiction, or a
non-comment source token. The only non-comment line is the owner's `std::abs(...)` at
`cl_MaxwellFactory.cpp:2999`. That is the pre-registered prediction confirmed, and it is the result
that matters most: 331 edits landed without introducing a defect.

**Codex, Part B — four findings in areas no partition covered**, all confirmed. The largest is a
numerical one and is reported, not fixed: the tetrahedron quadrature tables label "order"
inconsistently across two families. A numeric probe over every table (monomial exactness against the
reference tetrahedron) settles it:

| table | comment | exact through degree | reachable |
|---|---|---|---|
| `gauss_tet10` | "4th order" | none — see below | no call site |
| `gauss_tet20` | "6th order" | 5 | no call site |
| `gauss_tet35` | "7th order" | 6 | requested order 7 |
| `gauss_tet46` | "8th order" | 8 | requested order 8 |
| `gauss_tet56` | "9th order" | 8 | requested order 9 |
| `gauss_tet220` / `tet236` | "13th"/"14th" | 13 / 14 | yes |

The Witherden and Vincent tables label themselves by exactness degree; the Shunn and Ham tables label
themselves one higher. `fn_intpoints.cpp` dispatches on the requested order and treats the two as
interchangeable, so a request for order 7 or 9 receives a rule exact to degree 6 or 8. Comments now
state the measured degree; the dispatch is untouched and is Christian's call.

`gauss_tet10` is worse and dead: its `aWeights` array repeats the barycentric coordinate values, so
the weights sum to 2.5 instead of the reference tetrahedron's 1/6 and the rule integrates nothing
exactly. It has no call site. The comment now says so rather than claiming an order. Wiring it in
would silently corrupt every integral — worth a look before release even though it is unreachable.

A parser artifact worth recording, because it nearly became a false report: several live tables
(`tet35`, `tet165`) close their weight sum with `aWeights(k) = 1./6. - sum(aWeights)`. A literal-only
parse misses that line and makes a correct table look badly normalised. The first probe flagged two
live tables; reading the files refuted it.

**Grok — 16 findings on the cross-partition axis, all confirmed, all applied but three.** The
partitioned first round could not see any of these, which is the structural point:

- `CLAUDE.md` presented ε < 10⁻¹¹ as a BELFEM implementation choice. It is Messe et al. 2023's
  recommendation; the class default is `1e-6` and the shipped decks use it. Corrected in both
  `CLAUDE.md` and `doc/literature_references.md`.
- `doc/coding_philosophy.md` cited `CMakeLists.txt:98` for `USE_DEBUG` and `:101` for `USE_TEST`.
  Those lines are `USE_SUITESPARSE` and `USE_PARDISO`, both default OFF — a reader trusting the
  citation over the sentence concludes tests are off. Replaced with greppable tokens.
- The build-flag statement was GCC-only in three documents: GCC Fortran debug is `-O0 -fcheck=bounds`,
  not `-Og`, the native flag is Linux-x86 only, and ICC is `-O1 -xHost` / `-O0 -g` throughout.
- `CLAUDE.md` said the example decks' plugins are SHARED. All five are MODULE.
- `core_usage_guide.md` told the reader that release builds "may abort" on a `BELFEM_ASSERT`
  violation, contradicting its own later paragraph and `assert.hpp`.
- First-hour path: `examples/README.md` ran the solver from `cmake-build-debug/` after the README had
  the reader build in `build/`; its deck table listed five of seventeen decks as "the decks"; four
  decks load a plugin that the two-step recipe never builds (with a maintainer home directory in the
  `disk_pulse` build comment); `--slurm` is Lawrencium-specific without saying so; `getting_started`
  gave a `belfem-conf` path that only works from the repository root and never named PyYAML;
  `python/README.md` counted fifteen decks and cited the deleted `sidecoating`.
- `src/fem/thermal/doc/README.md` still sent readers to `hphiTrun` in its See Also, ten lines after
  calling it retired — the interrupted applier's fingerprint.
- The homology usage guide and README taught a six-argument `CutFactory` (the real one takes five,
  with no Poisson flag) and `Topology::analyze()` (the methods are `run()` and
  `run_on_enriched_mesh()`). Fixed in both, including twelve call sites in snippets.

**Three Grok findings deliberately not acted on, all Christian's call:**

- `doc/mainpage.md` publishes `ai_collaboration_protocol`, `ai_workflow_best_practices`,
  `lessons_learned` and `lessons_learned_evidence` as `@subpage`s, so the public Doxygen site carries
  the AI process documents and the incident catalogue. This was raised as an open question on
  2026-08-14 and never ruled on; it is a publishing decision, not a defect.
- `src/homology/doc/handoff_for_gregory_20260831.md` is a named-person working note living in a module
  doc directory, so Doxygen renders it as module documentation. Moving it is Gregory's and Christian's
  call.
- `cohomology_theory_and_implementation.md` walks an algorithm whose functions are not in the build.
  The page already carries a banner saying so; restructuring a careful document to survive a skim is
  an editorial decision.

**New code-side observations from this round** (added to the tracker): the `tet10` table above; the
order-convention mismatch feeding `fn_intpoints.cpp`'s dispatch; `scripts/check_doc_claims.py`
inventories only `STATIC|SHARED`, so it cannot see the `MODULE` the templates actually build, and its
`rglob` walks `tmp/`, where the only literal `SHARED` in the tree now lives — the probe passes on the
strength of a gitignored scratch directory; `examples/disk_pulse/src/CMakeLists.txt` defaults
`BELFEM_DIR` to a maintainer home path in a `CACHE` entry (code, not comment).

**Gates after the round:** `check_doc_claims.py` 37/37, `input_schema.yaml` parses, Doxygen 1.9.1
re-probe shows **zero warnings** outside `doc/lessons_learned_evidence.md` and the raw-layout noise.
188 files changed; the only non-comment source line in the whole diff is still the owner's.

## Open Questions

- The three unapplied partitions and three remainders (above) — ~245 audited findings with proposed
  text, waiting on a verify-and-apply pass. Tracked in `todo/doxygen_contradiction_sweep.md`.
- The auditors' code-side list above: each is a code decision, several look like real defects
  (`mSumW`, `check_saupd`, `faces_are_finalized`, the gmsh type cast).
- No Codex pass ran over the final diff (comment edits need no audit round by standing rule; the
  session limit made the optional pass moot). A build over the tree is owed as usual.
- The one uncommitted code line in `cl_MaxwellFactory.cpp` is Christian's; nothing in this session
  touched it and it must not be swept into a documentation commit.

## Files Updated

- 167 files, comment and markdown lines only (`git diff --stat`); by area: `src/mesh` 28,
  `src/fem/kernel` 24, `src/sparse` 18, `src/fem/*` 30, `src/physics/*` 16, `src/core|comm|containers|io|circuit` 30,
  `doc/` 5 (`input_file_reference.md`, `input_schema.yaml`, `coding_philosophy.md`, …), `CLAUDE.md`,
  `examples/block3d/src/matlib.cpp`, `examples/README.md`, `examples/scripts/README.md`,
  `todo/debt_register.md`
- tmp/ai_exchange/doxygen_contradiction_sweep.md (pre-registration, Codex, Grok, resolution)
- tmp/ai_exchange/doxygen_sweep_reports/ (16 audit reports, 12 apply logs, both briefs)
- todo/doxygen_contradiction_sweep.md (residue tracker), todo/README.md
- devlog/dl20260902_doxygen_contradiction_sweep.md (this file), devlog/README.md
