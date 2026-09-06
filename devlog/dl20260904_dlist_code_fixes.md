# D-list code fixes from the Doxygen contradiction sweep

**Date:** 2026-09-04
**Purpose:** Apply the twelve code-side observations (D1–D10, D13, D14) the 2026-09-02/03
documentation sweep routed to Christian, after his ruling, with a plan jury and a code-diff jury;
settle the publishing items P2 and P3
**Module:** cross-cutting (core, containers, sparse, mesh, fem/kernel, fem/interpolation,
physics/gasmodels, physics/gastables, scripts, examples)

Plan, both jury rounds, verification passes and reconciliation tables:
`tmp/ai_exchange/dlist_code_fixes.md`. Rulings and boxes: `todo/doxygen_contradiction_sweep.md`.

## What changed

| id | file | fix |
|---|---|---|
| D1 | `cl_EF_PENTA6TS.cpp` | `mSumW` 8 → 1 (wedge weights sum to 1; the 8 was HEX8TS's). New test `Penta6TsVolumeIdentity` |
| D2 | `arpacktools.cpp/.hpp/.f90` | `check_saupd` decodes `info = 3` (no shifts applied, raise ncv) as a warning like the nonsymmetric twin; it aborted the run after the Fortran driver had already extracted |
| D3 | `constants.hpp` | `y2k` = 2451544.5 (2000-01-01 00:00 UT); it held 2030's JD. No consumer |
| D4 | `random.hpp` | serial `random_seed()` now calls `std::srand`; the `/dev/urandom` read goes through `memcpy` with a short-read check instead of a `reinterpret_cast` of a `char[]` |
| D5 | `cl_FEM_Controller.cpp/.hpp`, `input_file_reference.md` | `check_thermal_diagnostics()` emits the two no-MUMPS warnings its magnetic twin emits; the doc hedge "not yet emitted" removed |
| D6 | five literals | `OldMaxwellFactory` → `MaxwellFactory` (2 sites + 2 guide quotes); `compute_surface` assert names itself, "parralell" fixed in both postproc files; "tape roller" → ThinShellFactory; `gastable` help heading; phantom `--verbose` lines dropped from both gas tools |
| D7a | `cl_Mesh.hpp` | `faces_are_finalized()` returns `mFacesAreFinalized` |
| D7b | `cl_Mesh_VtkWriter.cpp` | stray `mNodeMap`/`mNumberOfNodes++` in the owner loop removed (latent; output byte-identical) |
| D7c1 | `cl_Mesh_GmshReader.cpp` | `delete_unused_nodes_and_elements()` re-enabled, null-guarded, and clearing every reader container afterwards. Withdrawn on Christian's first reading, then reopened by him as a regression question ("there might be a reason why I commented that call out"): `git log -S` puts the comment-out in `173f3b31` (2026-03-02), the same hunk that replaced the vertex filter by an unconditional move of every vertex into the mesh, which the cleanup would then have freed; the filter returned on 2026-07-02, the cleanup did not. Codex SAFE and Grok SAFE-WITH-CONDITION (call site as proposed) on a dedicated regression brief; both walked every adoption path. Frees only isolated point elements and the nodes nothing flagged |
| D7c2 | `cl_Mesh_GmshReader.cpp/.hpp` | MSH 4.1 `$Entities` parsed into a per-dimension entity→physical-tag map; elements get their physical tag as the 2.2 path always did (first tag, magnitude). Withdrawn on Christian's first reading, restored on his second ("we can make the changes 2 and 3"); the code jury's bounds check on a declared-but-missing tag added |
| D7d | `meshtools.cpp`, `Mesh_Enums.hpp` | `element_type_from_gmsh` is an explicit switch (gmsh 36 → QUAD16, unknown → `BELFEM_ERROR` naming the number); inverse likewise. `QUAD16 = 32` untouched — it is serialised |
| D7e | `cl_Mesh_BfmFile.hpp` | `ProtoMesh` forward-declared in `belfem::mesh` |
| D7f | `cl_Mesh.cpp` | `Mesh::memory()` adds the vertices' footprint |
| D8a | `cl_Mesh_PeriodicityFactory.cpp` | k-d prune compares `|Δ|` with the Euclidean best, not `Δ²` |
| D8b | `cl_Mesh_OrderConverter.cpp/.hpp` | `check_input_mesh()` reads `mInMesh`; `mMesh` initialised to `nullptr` |
| D9a | `cl_FEM_DofMgr_DofData.cpp` | METIS/SCOTCH gate on `reorder_dofs()` deleted; the body only ever ran `symrcm` |
| D9b/c | `cl_FEM_DofMgr_SolverData.hpp`, `cl_FEM_DofManager.hpp`, `cl_FEM_DofMgr_EigenValues.cpp` | `jacobian()` → `system_matrix()` (it returned the system matrix); `run_shift_invert()` error text names itself |
| D9d | `cl_FEM_Calculator.cpp/.hpp` | facet link bounds the slave coordinate loop by the slave's node count and dispatches the slave integration index through `mFunSlaveIntegrationIndex`; the four index helpers take a `mesh::Facet *` |
| D9e | `cl_FEM_Dof.hpp` | `set_field_index` declared with `index_t` (broke under `BELFEM_INT64`) |
| D10a | `stringtools.cpp` | `MN` 1e3 → 1e6, `muV` 1e-3 → 1e-6. Tests for the force and volt families |
| D10b | `stringtools.cpp` | `to_real( "" )` returns NaN as the header promises; test flipped |
| D10c–e | `cl_DynamicBitset.cpp`, `fn_compute_permutation.hpp`, `cl_GM_Helmholtz.hpp` | `gNoIndex` sentinel; guard renamed; dead `check_parent()` declaration removed |
| D13 | `scripts/check_doc_claims.py` | library-kind probe sees `MODULE`, skips `tmp/`, `archive/`, `nonfree/`, `.claude/` and build trees; new row "plugin library kind". 38/38 |
| D14 | `examples/disk_pulse/src/` | last maintainer home path (in `bgpulse.cpp`'s build comment) replaced; the override comment its siblings carry added |
| P2 | `todo/handoff_for_gregory_20260831.md` | moved out of `src/homology/doc/`; nav and README follow |
| P3 | `cohomology_theory_and_implementation.md` | status banner under the title |
| P1 | `doc/mainpage.md`, `Doxyfile.in`, `doc/README.md` | the two lessons-learned pages no longer render on the site; the AI protocol and workflow pages still do (Christian's ruling) |

Behaviour changes a user can observe, all named in the plan and confirmed by both auditors:
serial ARPACK start vectors are no longer bit-reproducible (D4, MPI builds never were);
MSH 4.1 meshes now carry physical tags into the Exodus "physical tag" variable and the
`ProtoMesh` transfer (D7c2, parity with 2.2); decks using `MN` or `muV` were 1000× off and are
now right (D10a); an unsupported gmsh element type aborts with its number instead of a
misleading node-count message (D7d).

## Audit record

**Plan jury** (Codex terra/high, Grok 4.6/high): Codex asked for four corrections, all adopted —
`memcpy` for the seed read, a `return` after every `BELFEM_ERROR` in the new switches, clearing
the reader's alias containers too, and recorded its objection to deleting the reorder gate.
Grok raised a blocker: that `delete_unused_nodes_and_elements()` would free live sideset
geometry because sideset elements are never flagged. **Refuted on the source**:
`GmshReader::create_mesh()` flags every sideset facet's element and its nodes and aborts unless
every facet was reached; Grok had read `create_sidesets()` and stopped one function short.
Grok's other corrections were adopted: the second `check_saupd` site, `92 → HEX64` stated
explicitly, `2451544.5` under the name `y2k`, the four-helper `Facet *` shape for D9d, two
comment sites for the rename.

**Code jury** (Codex terra/xhigh; Grok 4.6/xhigh failed on turns, retried with a narrowed
brief and confirmed all six named hunks, refuting one comment of mine — gmsh does not sign
physical tags): Codex found the six guide snippets the rename orphaned, the stale precondition comment
on `reorder_dofs()`, the unconstrained seed template, three core-guide sentences contradicting
the new seeding, and a release-build out-of-bounds read on a malformed `$Entities` line — all
confirmed and fixed.

**Cleanup regression audit** (dedicated brief to both, terra/xhigh and 4.6/xhigh): Codex SAFE,
Grok SAFE-WITH-CONDITION where the conditions are the call site as proposed. Grok added the
load-bearing detail that `Mesh::finalize()` unflags block elements, so the call must never move
below it, and that the delete loops were not null-safe on a short read — guarded.

## Gates

- `python3 scripts/check_doc_claims.py` → 38/38 (was 37/37).
- Every touched translation unit compiled with the library flags of `cmake-build-debug`
  (`-Og -g -Wall -Werror -pedantic-errors …`, `-c -o /dev/null`); `arpacktools.f90` with
  `gfortran -fsyntax-only`. That is the compile rung, not the build, and not `make check`.
- **`make check` passes** (Christian, 2026-09-04, full suite on the rebuilt tree) — the
  regression rung: the three new tests and the flipped one run green, and the reader with the
  cleanup re-enabled reads every fixture the suite loads.
- **Still owed:** one MSH 4.1 deck read end to end with its Exodus "physical tag" variable
  inspected (D7c2); Valgrind on the non-MKL tree over one MSH 2.2 and one MSH 4.1 deck with
  sidesets for D7c1 — a green suite shows no crash, not the absence of a double free.

## Open, for Christian

- **D9a split verdict.** With the gate gone, a programmatic
  `KernelParameters::set_reordering_method( NATURAL )` would run RCM instead of aborting. No deck
  key reaches that knob and no FEM caller sets it. Deleted as ruled; a one-line
  `BELFEM_ERROR` excluding NATURAL/AUTOMATIC is the alternative.
- **P1** ruled by Christian the same day ("agreed"): the protocol and workflow pages stay on the
  site; `lessons_learned.md` and `lessons_learned_evidence.md` left the `@subpage` tree in
  `doc/mainpage.md` and are excluded in `Doxyfile.in`; `doc/README.md` marks them repository-only.

## Noticed, not touched (out of scope, relayed)

1. ~~`Element::memory()` returns `sizeof( Element )` for every derived element~~ — fixed the same
   day, see the follow-up section below.
2. `mBoundaryEdges` are called "references" in `Mesh::memory()`; no block owns 1-D boundary edges
   and the element delete loop in `~Mesh` is commented out — likely leaked.
3. `stringtools.cpp` second `cd` branch unreachable; ampere family has no `muA`.
4. `gastables::Arguments` parses no verbosity flag; wiring one is a feature.
5. `cl_FEM_Element.cpp` `compute_edge_directions_thinshell()` tests `physical_tag() == 1`, which
   after `MaxwellFactory::set_physical_tags_for_elements()` means "first material", not a
   direction flag. Single caller. Looks accidental.
6. `more/gmsh/tet56.msh` (gmsh type 31) still has no BELFEM enumerator; the new error names it.

## Follow-up: `Element::memory()` (same session, Christian: "let's fix the Element::memory() function")

Exchange: `tmp/ai_exchange/element_memory.md`. `Element::memory()` was non-virtual and started
from `sizeof( Element )`, so every element was charged the base object size; it also charged
the edge and face arrays unconditionally, though the template allocates them on demand, and the
facet array, which no element ever allocates. Both auditors confirmed and added the same three
corrections: gate the facet array like the neighbour array, count the source weights (sized by
the same count as the sources), and say plainly that the counts are fill counts, since no
capacity is stored.

Now: `Element::memory()` is virtual and returns `sizeof( Element )` plus a protected
`array_memory()` for the arrays the element and its bases own; `ElementTemplate::memory()`
overrides with `sizeof( *this )`, the base arrays, the node array, and the edge and face arrays
under their flags. Diagnostic only; nothing consumes the numbers and no test pins them. Compiled
in both trees through `cl_Mesh.cpp`, `cl_Facet.cpp`, `cl_Element_Factory.cpp`. Code jury (both):
no P0; the accounting comment was rewritten per array (fill counters vs stored capacity vs
facet-count sizing), and the dormant `allocate_facet_container()` now asserts one slot per facet.

**Second round, same session (Christian: "same treatment for Node::memory, Edge::memory and
Face::memory; reset_control_point_container should nullify the pointer"):** `mesh::Vertex` got a
protected `array_memory()` that charges its five owned arrays by their **stored capacity** —
never by `number_of_*()`, which `Facet`, `Face` and `Segment` override to describe the wrapped
element or the master/slave pair (Grok's load-bearing point: the fill-count route would have
double-counted every sideset facet) — plus the graph vertices, dofs, sources and weights as fill
counts. `Node`, `Edge`, `Face` and `Facet` (`Mesh::memory()` sums facets too, and their
facet-neighbour and thin-shell graph arrays are live) now build on it. Both auditors caught that
`number_of_duplicates()` hides the one slot a duplicate node holds for its original (counter
`-1`); the count is `std::abs( mNumberOfDuplicates )`. Hardening in the same files: the
duplicate reset and the control-point reset null their pointers (each was a double free on a
reset followed by destruction); `allocate_element_container()` got the `== nullptr` precondition
its siblings have; the five `Vertex` capacity allocators and the two narrow `Element` counters
(control points, elements) refuse a count their `uint8_t`/`uint16_t` field cannot hold — checked
**before** anything is allocated, after the code jury caught the first version checking after the
`malloc` and the truncating store, which on the throwing error hook would have leaked. The dormant
`Element::allocate_facet_container()` size check is always-active. Compiled in both trees across
mesh and homology consumers. No test pins a byte count: reviewed, not verified. One behaviour
change: an explicit request for more than 255 elements on a vertex aborts instead of truncating.

**Still not honest (relayed):** `ControlPoint::memory()` omits its `Basis` dofs, sources and
weights; a `ControlPoint` is a `Basis`, not a `Vertex`, so neither helper reaches it. A `Basis`-level
helper would close that and let `Element` and `Vertex` share the three base terms.
`Curve::memory()` charges only its segment-pointer array, not the segments it deletes nor the
elements they wrap. The duplicate array of a `Node` stores no capacity, so it is billed by fill count
like sources and dofs. `increment_*_counter()` on a `Vertex` are debug asserts, so a counter can
wrap in release before the new width check sees it.

Deleted on Christian's ruling: `src/fem/postproc/cl_Surface.{cpp,hpp}` was in no CMake list,
nothing in `src/`, `tests/`, `examples/`, `nonfree/` or `archive/` included it or named
`fem::Surface` (the `Surface` classes in `nonfree/manta` and `archive/satellite` are their own),
and it named a header that does not exist (`cl_Element_Penta6TS.hpp`). Last touched 2025-03-19.
Relayed, not touched: `Element::reset_control_point_container()` frees
the array without nulling the pointer (Grok); the destructor then frees it again if the count
was reset. Neither is part of this change.
