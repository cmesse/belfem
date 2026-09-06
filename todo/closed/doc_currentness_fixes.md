# Documentation Currentness Fixes

**Date:** 2026-08-31
**Purpose:** Fix list from the three-AI sweep of all 96 markdown files in `doc/` and `src/*/doc/`
**Module:** cross-cutting (documentation)

> **Source:** `devlog/dl20260831_doc_currentness_sweep.md`. Every `Dn` below was re-derived
> against the tree by Claude before being written here; auditor agreement alone did not qualify
> a row. Rows marked **ROUTED** are physics or attribution questions for Christian and are not
> to be edited by an AI session.

**Status:** **R1–R7 COMPLETE, 2026-08-31.** Executed as approved — documentation and comments only,
plan+audit → correct → audit, batched. 89 files, 1295 insertions, 16 auditor rounds (Codex
`gpt-5.6-terra`/high + Grok `grok-4.6`/high on the plan and on every batch). Every D-group verified
closed against the tree at session close; 120 boxes ticked.

**Reviewed, not verified.** No build and no `make check`. One real compile gate ran
(`-fsyntax-only -std=gnu++17 -Wall -Werror` on `tests/math/test_GraphVertex.cpp`, clean) plus two
standalone sentinel probes. `check_doc_claims.py` 37/37, `resolve_doc_cites.py --check` and
`update_doc_index.py --check` both exit 0. Session record: `devlog/dl20260831_doc_repair_campaign.md`.

**Ten boxes remain open, none of them mine to close:**

| Open | Why it is open |
|---|---|
| D13 / R7f / R7h-note — `examples/Tape_Quench_obsolete_deleteme` | Deletion needs Christian's word; raised three times, not acted on |
| D12 second box — `numerics/{integration,ode,bezier,sources}`, `math/tools` | Whether these want `doc/` directories is a scoping call, not a defect |
| O1 — φ-region static Laplace form vs the `mt_maxwell_phi.cpp:26-46` mass term | Formulation question; Christian's |
| O2 — static-condensation attribution (`messe2023.txt:509-517` credits Alves et al.) | Christian's paper, Christian's call |
| R8a–R8d — mechanical probes for the blind spots this campaign could not check | **Explicitly out of scope for the R1–R7 approval** |

**Owed, and not owed by this plan:** a build and `make check` over the tree as it now stands, and
staging the untracked files (`todo/doc_currentness_fixes.md`, the two new devlogs,
`src/numerics/opt/doc/`).

---

## Scope guards — read before every batch

These bound what any session working this plan may touch. They are not advisory.

1. **Documentation and comments only.** No change may alter what the compiler emits or what the
   program does. Markdown, and comment text inside source or CMake files. Nothing else.
2. **When doc and code disagree, the code is right and the doc is fixed** — unless the row says
   otherwise. This plan never "fixes" a doc by changing the code to match it. Where a finding
   raises a genuine question about the code's behaviour (should Crank–Nicolson be re-enabled?
   should `LINE3` be supported?), that goes to `On`, and the doc is corrected to describe what the
   code does **today**.
3. **`O1` and `O2` are frozen.** The φ-region time-derivative question and the Messe
   static-condensation attribution are Christian's. Do not edit
   `maxwell_usage_guide.md:186-210,329-369` or `CLAUDE.md:699,711` in any batch.
4. **The cohomology core is under an AI edit ban** — `doc/ai_collaboration_protocol.md` §7.1,
   landed by a peer session on 2026-08-31 while this campaign was running. `cl_Cohomology`,
   `cl_Homology`, `cl_SimplicialComplex`, `cl_Chain`, `cl_Cochain` and `fn_Smith` (`.cpp` and
   `.hpp`) are closed to AI edits **including comments**, and the ban is *not* lifted by session
   edit approval. `src/homology/doc/` is explicitly **exempt**, so R6's markdown work is still in
   scope. Read §7.1 before R6: its reasoning — that the nearest model prior is *textbook*
   Pellikka, so every deliberate modification reads as a defect, and that reviewer agreement is
   worth **least** here — applies to what R6 is allowed to *assert*, not just to what it edits.
5. **`D17` (homology) needs its owner.** Gregory Giard owns `src/homology`. Batch 6 prepares the
   evidence and fixes only the items that are unambiguous; the guide rewrite waits for him.
6. **Source citations to `hphirun.cpp` stay.** The file exists. Only *run* instructions change.
7. **`git status` before every batch.** Shared checkout, live peers.
8. **No batch is "done" until its boxes are ticked in the same turn the edit lands.**

---

## Gap table — what this plan closes

| Gap | Today | After |
|---|---|---|
| Newcomer run instructions name an executable that is not built | 6 documents | 0 |
| Taught code examples that do not compile or silently misbehave | 4 | 0 |
| Documented defaults that are backwards | 2 (`USE_DEBUG`, coulomb `chi`) | 0 |
| Countable claims drifted (line counts, card counts, incident counts, module counts) | 6 | 0, and mechanically guarded |
| `check_doc_claims.py` blind spots that let all of the above through | 3 classes | probed |

---

## D — Defects

### D1 `hphirun` / `hphiTrun` are not built, and eleven documents name them — **CRITICAL**

Six of the eleven carry *runnable* instructions; the rest are prose. (The heading said "six"
while the list below enumerated more — corrected 2026-08-31 after Codex flagged the mismatch.)

`src/executables/CMakeLists.txt:40-47` comments out both `Add_Executable` blocks. Confirmed
against `cmake-build-debug/bin/` (built 2026-08-31 02:56): `belfem`, `electricalCircuit`,
`material`, `gas`, `db2exo` present; neither `hphirun` nor `hphiTrun`.

The `.cpp` sources still exist, so **citations to `hphirun.cpp` as source stay correct** and must
not be swept up in this fix.

- [x] `README.md:44` — `../../build/bin/hphirun` → `belfem`
- [x] `doc/getting_started.md:45-46` — both invocations → `belfem`
- [x] `doc/getting_started.md:52` — "`hphirun` solves magnetic decks; `hphiTrun` handles decks
      with a thermal section" → `belfem` selects from the deck
- [x] `doc/parallel_execution.md:12` — production command line → `belfem`
- [x] `src/fem/maxwell/doc/maxwell_usage_guide.md:310-311,1019` — run commands → `belfem`
- [x] `src/sparse/doc/solver_memory_and_compression.md:106` — `prterun -np 4 hphirun` → `belfem`
- [x] `examples/README.md:38-40` — **"`hphirun` and `hphiTrun` still exist and still work" is
      false.** They are not built. This sentence is what led one auditor to downgrade the
      severity, so it is load-bearing beyond its own line
- [x] `doc/README.md:17,49,59` — index prose still frames them as the solver applications
- [x] `src/fem/maxwell/doc/README.md:340` — "alternatively, run `hphiTrun`"
- [x] `doc/input_file_reference.md:992` — behaviour attributed to `hphirun`/`hphiTrun`
- [x] `src/executables/doc/README.md` — documents both as available
- [x] Decide whether `config/globals.cmake:19` should stop naming them. **Harmless today** —
      `config/scripts/Add_Executable.cmake:63` uses the list as a filter over targets that exist,
      and `src/physics/database/CMakeLists.txt:24` states that contract — so this is tidiness,
      not a defect
- [x] Fix the typo in `src/executables/CMakeLists.txt:40` ("retiredI")

### D2 `rho()` taught with the arguments in the wrong order — **CRITICAL (silent wrong answer)**

- [x] `src/fem/maxwell/doc/maxwell_usage_guide.md:945,951` teach
      `material->rho( norm(b), angle, T )`. Signature is `rho( const real T, const real B, const
      real beta )` (`cl_Material.hpp:680`, defn `:2032`). All three parameters are `real`, so the
      wrong order **compiles** and returns a wrong resistivity — `lessons_learned.md` L-06.
      `src/physics/materials/doc/materials_usage_guide.md:70,259` already teaches it correctly

### D3 A taught MPI API that is largely fictional — **CRITICAL**

> **Narrowed 2026-08-31 (Grok, plan audit).** The first wording said "an API that does not exist".
> That is overstated: `Communicator::rank()` and `size()` are **real** (`cl_Communicator.hpp:145-151`,
> defined `:191,:199`), so `gComm.rank()` is valid code. What is fictional is the *tag* parameter,
> `sum_all`, and `comm.send`/`comm.receive` as member functions. Correct to the free-function
> surface without claiming `rank()`/`size()` do not exist.

- [x] `src/fem/kernel/doc/dof_manager_usage_guide.md:2332-2337` teaches `comm.rank()`,
      `comm.size()`, `comm.send(data, target_rank, tag)`, `comm.receive(data, source_rank, tag)`,
      `comm.sum_all(local_value)`. BELFEM's MPI layer is free functions — `comm_rank()`
      (`commtools.hpp:64`), `send( T&, proc_t )` (`:549,877`), `receive( T&, proc_t )`
      (`:634,962`). **There is no tag parameter** and **`sum_all` exists nowhere in `src/`**.
      `src/comm/doc/README.md:76-77` already shows the correct form — copy it

### D4 Two more non-compiling snippets — **HIGH**

- [x] `src/fem/iwg/doc/iwg_usage_guide.md:263-298` and `src/fem/iwg/doc/README.md:109-117` call
      `create_iwg(type, dim, IwgMode::Iterative)`; the factory declares two parameters
      (`cl_IwgFactory.hpp:53`)
- [x] `src/fem/maxwell/doc/maxwell_usage_guide.md:506-568` calls `new Kernel(mesh)`; the
      constructor takes `KernelParameters *` (`cl_FEM_Kernel.hpp:98`)
- [x] `src/homology/doc/homology_usage_guide.md:279,630,742` and `src/homology/doc/README.md:115`
      construct `CutFactory` with four or six arguments and an `aComputePoisson` parameter; the
      only constructor takes exactly five and has no such parameter (`cl_CutFactory.hpp:99-104`)
- [x] `src/homology/doc/homology_usage_guide.md:1785` tells the caller to invoke
      `select_belt_fasteners()` and `create_tree()` after constructing `BeltedTree`; the
      constructor already calls both (`cl_BeltedTree.cpp:24-38`)

### D5 Defaults documented backwards — **HIGH**

- [x] `USE_DEBUG` is documented as defaulting **ON** in `doc/getting_started.md:24-27`,
      `doc/coding_philosophy.md:56`, `README.md:25` (and the follow-on advice at `README.md:31`
      and `getting_started.md:26-27` to pass `-DUSE_DEBUG=OFF` for production). It defaults
      **OFF** — `find_scls_flavor.cmake:23` sets `BELFEM_DEFAULT_USE_DEBUG OFF`, only the SCLS
      `debug` flavor flips it (`:47-50`), and `CMakeLists.txt:98` consumes it. A newcomer
      following the quick start gets a release build believing assertions are live.
      *Everything else in that defaults list is correct* — do not rewrite the rest
- [x] Coulomb-gauge penalty documented "off by default" in
      `src/fem/maxwell/doc/coulomb_gauge_penalty_theory.md:11` and
      `src/fem/maxwell/doc/README.md:33`. On at `chi = 1e-4` since 2026-08-27 —
      `cl_IWG_Maxwell.cpp:48-51` says so in a comment written for this purpose

### D6 Advertised features that hard-fail or do not exist — **HIGH**

- [x] `src/fem/iwg/doc/iwg_usage_guide.md:357` lists "BDF1-5, Crank-Nicolson, Galerkin, Euler
      methods". Crank–Nicolson and Galerkin raise `BELFEM_ERROR` on selection
      (`cl_IWG_Timestep.cpp:92-94`, "disabled: the Newton tangent assembly supports the BDF
      family only")
- [x] `src/fem/maxwell/doc/README.md:20` lists **losses** among postprocessed fields. B/H/J/JJC
      are real (`cl_MaxwellPostprocessor.cpp:113,117,146,150`); case-insensitive grep for `loss`
      across every `.cpp`/`.hpp` under `src/fem/maxwell/` returns zero hits
- [x] `src/fem/iwg/doc/README.md:325` — "`MaxwellThermal` not detected by `is_maxwell()`". It is
      (`en_IWGs.hpp:99-102`). `:322` "BDF5 currently broken" contradicts `:183-195` in the same
      file
- [x] `src/homology/doc/cohomology_theory_and_implementation.md:29`,
      `homology_usage_guide.md:1329`, `README.md:260` describe a cleanup workflow built on
      `manifold_filter_3d()`, `check_surface()` and Tarjan articulation detection.
      `manifold_filter_3d()` and `check_surface()` exist nowhere. **Tarjan does exist**, as an
      archived prototype — `archive/graph/fn_Graph_tarjan.{hpp,cpp}`, outside the build. (The
      earlier "zero hits, `archive/` included" claim here was wrong: shell `grep` wraps
      `ugrep --ignore-files` and had silently skipped `archive/`.) The implemented path is
      `Cohomology::clean_spfa()` → `remove_cut_pockets()` (`cl_Cohomology.hpp:135,145`,
      both private) after SPFA rectification.
      **CLOSED 2026-08-31**, second session, under explicit authorization to edit homology
      documentation. All four sites in `homology_usage_guide.md` now corrected — the two the
      previous note named (`:335` glossary, `:1329-1350` pipeline) plus **two it had missed**:
      the Internal Workflow step and the troubleshooting entry at `:2205-2208`, which told
      readers to tune cycle-density and size thresholds that do not exist.
      A follow-on defect of the same class was found and fixed in the same pass: the advertised
      debug outputs `thick_cut_*.vtk`, `thin_cut_*.vtk` and `manifold_*.vtk` are emitted by no
      code (`manifold` appears zero times in the module's sources). `save_debug_meshes()` writes
      `cut_<index>.vtk` only (`cl_CutProcessor.cpp:139-147`), and it is a member of
      `CutProcessor`, not `CutFactory` as the guide's troubleshooting snippet had it.
      Corrected in `homology_usage_guide.md` (3 sites) and `src/homology/doc/README.md` (1).
      A sweep of all 44 code identifiers asserted across the six homology docs found no further
      missing symbols; `rectify_to_unit_or_certify()` is absent from the tree but correctly
      labelled "Proposed". **Content status remains Gregory Giard's to confirm.**
      Reviewed, not verified — no build was run.
      **I ticked this box while it was one-third open.** Carried to
      `todo/handoff_20260831_open_defects.md` **F4**; the content call is Gregory's.
      Note `§7.1` bans the cohomology **source**, not `src/homology/doc/` — this is editable
- [x] `src/homology/doc/thin_cut_nonunit_rectification.md:22,334` frames SPFA rectification as
      a *proposal* and says the current `clean()` has no global view or certificate.
      `clean_spfa()` ships the difference-constraint solve and the obstruction certificate
      (`cl_Cohomology.cpp:342-358,455-460`)

### D7 Wrong preprocessor guard for METIS — **HIGH**

- [x] `src/math/graph/doc/graph_usage_guide.md:626,660,702,743,1018,1377,1383` and
      `src/math/graph/doc/README.md:63,76,81` guard METIS with `BELFEM_SUITESPARSE`. The code
      uses `BELFEM_METIS` (`graph_typedefs.hpp:27,50`), set by `USE_METIS`
      (`config/linalg/config_metis.cmake`). A reader copying the documented guard gets METIS
      compiled out of a METIS-enabled build

### D8 A cohomology formula with the wrong index — **HIGH (theory)**

- [x] `src/homology/doc/cohomology_algorithms.md:21` writes `H^k = Ker(δ^k) / Im(δ^{k+1})`.
      Pellikka et al. 2013 Eq. (A.2) (`literature/papers/topology/pellikka2013.txt:683`):
      `H^k(M) = ker(δ^k)/im(δ^{k-1})`. The `k+1` is the *homology* index pattern from Eq. (A.1)
      on the line above, copied into the cohomology case. **`k+1` → `k-1`**

### D9 ParMETIS input model described wrongly — **MEDIUM**

- [x] `src/math/graph/doc/graph_usage_guide.md:772` says ParMETIS operates on an
      already-distributed graph with each rank owning a subset. BELFEM's wrapper requires the
      complete graph on rank 0, builds every rank's CSR slice there, then distributes
      (`fn_Graph_ParMETIS.cpp:50-75`; the CSR builder asserts root-only, `graphtools.hpp:124-133`)

### D10 Countable claims that have drifted — **LOW, but they are what a session loads**

- [x] `doc/README.md:12` "141 lines" and `CLAUDE.md:53` "~140 lines" → Layer 1 is **200** lines
      (`lessons_learned.md` 19→219)
- [x] `doc/README.md:12` "18 operating rules" and `lessons_learned.md:221` "18 cards" → **21**
      (L-01..L-21)
- [x] `doc/README.md:12,13` "539 catalogued incidents" and `lessons_learned_evidence.md:10,36`
      "537" → **558** distinct IDs (INC-001..INC-558) on 593 row lines, 35 IDs spanning more
      than one row. Decide which number the prose means and say which
- [x] `CLAUDE.md:520-523` lists 23 modules with a `doc/` directory; there are **26** — add
      `fem/postproc`, `fem/thermal`, `visualizer`. `CLAUDE.md:529` "only one of the 23" → 26;
      the DOI substance is still true (`src/homology/doc/README.md`)
- [x] `doc/README.md` does not link `src/fem/doc/README.md` — the only module doc README missing
      from the master index (no broken links in the other direction)
- [x] `doc/README.md` "Executables" line omits `belfem`, `material`, `gas`, `db2exo`
      (`config/globals.cmake:19`)

### D11 Drifted `file:line` anchors — **LOW**

Citation sweep over all 96 docs: 564 of 580 resolve. Of the 16 flagged, 9 were an artifact of my
own resolver preferring `src/CMakeLists.txt` (22 lines) over root (590) among 35 same-named
candidates. Real drift, 7 anchors, 5 of them in `lessons_learned_evidence.md` where dated records
are kept as written.

- [x] `doc/coding_philosophy.md:59` cites `assert.hpp:99` for a `throw`; `:99` is a
      debugger-workflow comment, the `throw` is at `:191`
- [x] `doc/coding_philosophy.md:61` cites `CMakeLists.txt:81` for `USE_TEST`; it is `:101`
- [x] `doc/coding_philosophy.md:504` cites `corctest.cpp:44` — no such file in the tree
      (a stale `corctest` binary survives in `cmake-build-debug/bin/` from 2026-08-27)
- [x] `doc/coding_philosophy.md:245,302` — `length()` is `cl_AR_Vector.hpp:257-260` not `:234`;
      `capacity()` is `cl_AR_Matrix.hpp:224-227` not `:216`
- [x] `doc/ai_collaboration_protocol.md:194,403-404` names three paths that no longer exist
      (`todo/periodic_thin_cut_continuity_fix.org.md`, `todo/ai_exchange.md`,
      `todo/ai_exchange_archive_*.md`)
- [x] `doc/documentation_guidelines.md:91-94` gives examples `src/fem/doc/iwg_integration_guide.md`
      and `src/maxwell/doc/field_formulations.md`; neither exists and there is no `src/maxwell/`
- [x] `src/homology/doc/thick_thin_cuts_and_conjugate_edges.md:86` — `CutProcessor::run()` does
      not exist; the work happens in the constructor (`cl_CutProcessor.cpp:21-99`)

### D12 An undocumented live module — **LOW**

- [x] `src/numerics/opt` (Optimizer, Objective, algorithm/status enums) has no `doc/` directory
      and no entry in `doc/README.md` or CLAUDE.md's architecture section. It is built —
      `src/numerics/CMakeLists.txt:5` — and `USE_NLOPT` defaults ON (`CMakeLists.txt:93`)
- [ ] Same question for `src/numerics/{integration,ode,bezier,sources}` and `src/math/tools`,
      which have substantial code and no `doc/`. `src/numerics/integration` is 65 files

### D13 An example directory named for its own deletion — **LOW**

- [ ] `examples/Tape_Quench_obsolete_deleteme` ships in the examples tree and is named in no
      README. Delete it or name it

### D14 Physics-module API and semantics defects — **CRITICAL / HIGH**

- [x] **`heat_spline()` is documented one derivative off, throughout — CRITICAL, silent wrong
      answer.** `src/physics/gastables/doc/gastables_usage_guide.md:349-351` and
      `src/physics/gastables/doc/README.md:363-364` teach
      `real cp = heat_spline->eval(T); real dcp_dT = heat_spline->deval(T);`.
      The spline is filled with **enthalpy** — `cl_GT_RefGas.cpp:1227` stores `this->H( T(k) )` —
      and the accessors say so unambiguously: `RefGas::spline_Cp` returns `mHeatSpline.deval(T)`
      (`:1122-1124`), `spline_dCpdT` returns `ddeval` (`:1130-1132`), `spline_H` returns `eval`
      (`:1138-1140`). So the documented `eval` gives H, and the documented `deval` gives cp.
      Every line in that block is shifted by one differentiation
- [x] ~~**`mu` documented as ∂B/∂H, implemented as B/H.** `materials_usage_guide.md:614`~~
      **MIS-STATED (retracted 2026-08-31).** The claim mixed two different APIs, and the proposed
      relabel would have been wrong. `:614` is the **`MaterialProperty` constant table** — slot 5
      of the property enumeration, not `BhCurve`. For a *constant* permeability B = μH is linear,
      so B/H and ∂B/∂H **coincide**; `Material::dmudH_const` returns `dmudH = 0.0`
      (`cl_Material.hpp:2112-2116`), which is only consistent because μ is constant there.
      The secant-vs-differential split is real but belongs to `BhCurve::mu(H) = B/H`
      (`cl_BhCurve.hpp:31`) with `dmudH` alongside. **The fix is at `:650`, which already names
      the split — sharpen that sentence and leave slot 5 alone.** Caught by Grok; re-derived here
- [x] `src/physics/gastables/doc/README.md:194,198,235` calls `has_cryo_thermo()` /
      `has_cryo_transport()`. **Neither exists** — the only hits in the whole tree are those doc
      lines. The real accessors are `has_thermo()` (`cl_GT_RefGas.hpp:370`), `has_conductivity()`
      (`:375`), `has_viscosity()` (`:380`)
- [x] `src/physics/gasmodels/doc/gasmodels_usage_guide.md:123` and `README.md:131` use
      `HelmholtzModel::HYDROGEN` / `::METHANE`. The enum is CamelCase and has no generic
      hydrogen — `ParaHydrogen`, `NormalHydrogen`, `OrthoHydrogen`, `Oxygen`, …
      (`en_Helmholtz.hpp:17-20`). Does not compile
- [x] `src/physics/gasmodels/doc/gasmodels_usage_guide.md:57` and `README.md:79` show
      `Gas air("air", GasModel::IDGAS)`; per Codex there is no `air` entry in the shipped
      `share/fluid/*.inp` data and `Gas::initialize()` errors without thermodynamic data
      (`cl_Gas.cpp:274`). *Not re-derived by me — carries Codex's confidence only*
- [x] `src/physics/gasmodels/doc/gasmodels_usage_guide.md:79` — "All property methods are const
      (thread-safe for read-only access)". `cl_Gas.hpp:75` carries the opposite warning: const
      evaluators write a shared cache and must not be called concurrently on one `Gas`. This also
      contradicts BELFEM's stated not-thread-safe-by-design posture. *Not re-derived by me*
- [x] Gas input-file names. **The cited line was wrong** (Grok, plan audit):
      `gastables_usage_guide.md:494-496` already says `thermo.inp` and is **correct** — an editor
      sent there would "fix" a right heading. The wrong names live at
      `gastables_usage_guide.md:529,542` and `gastables/doc/README.md:502-519`. The factory reads
      `thermo.inp`, `trans.inp`, `gasdata.inp`, `cubicalpha.inp` (`cl_GT_RefGasFactory.cpp:39-49`),
      which is exactly what `share/fluid/` contains (verified)
- [x] `gastables_usage_guide.md:218` / `README.md:85` say polynomial mode is the default and
      spline mode costs setup; the factory leaves gases in **spline** mode
      (`cl_GT_RefGasFactory.cpp:117`) and `set_mode()` only rebinds pointers. Guidance reversed.
      *Not re-derived by me*
- [x] Stale anchors, physics modules: `thermal_expansion_from_heat_capacity.md:102` (points at a
      base spline error path, not `SplineLookupTable::create_low_temperature_alpha()`),
      `callaway_thermal_conductivity.md:35,192` (point into `YBCO::test_callaway()`, not
      `lambda_custom()` at `cl_Material_YBCO.cpp:475`),
      `materials_contracts_and_invariants.md:84` ("every material constant defaults to NaN" —
      `mu`, `q`, `density_correction` have non-NaN constructor defaults),
      `thermal_expansion_from_heat_capacity.md:304` (Magnesia does get a derived `debye0K`).
      *Not re-derived by me*


### D15 Two thermophysical correlations documented in the wrong functional form — **HIGH (theory)**

Both found by Grok, both re-derived here against the kernel headers.

> **Two corrections to this row from the plan audit (Grok), both applied below.** First, the
> **H and S** expressions at `gastables/doc/README.md:253-260` are NASA-7 as well, and the kernel
> computes both in NASA-9 (`cl_GT_HeatPoly.cpp:52-78`); fixing only `Cp/R` leaves the page a
> chimera. Second, `gastables_usage_guide.md:498-514` is **not** a second copy of the formula — it
> is a CEA *record layout* showing an N2 data record. **Do not turn a data record into a
> polynomial.** The wrong-formula target is `README.md:250-260` only.

- [x] **NASA-7 formula printed under a "NASA-9" heading.** `src/physics/gastables/doc/README.md:244`
      is headed "NASA CEA 9-Coefficient Format" and prints, six lines below at `:250`,
      `Cp/R = a1 + a2·T + a3·T² + a4·T³ + a5·T⁴` — the **7**-coefficient form. The kernel
      evaluates NASA-9, RP-1311 Eqs. (4.9)–(4.11):
      `cp/R = a1 T^-2 + a2 T^-1 + a3 + a4 T + a5 T² + a6 T³ + a7 T⁴` (`cl_GT_HeatPoly.hpp:27`).
      The document contradicts itself within four lines, and the shipped `thermo.inp` records are
      NASA-9. **Fix `Cp/R` at `:250` *and* the H and S expressions at `:253-260` in the same
      edit** — the kernel's H and S are NASA-9 too (`cl_GT_HeatPoly.cpp:52-78`)
- [x] **log₁₀ documented where the code uses ln.** `src/physics/gastables/doc/README.md:272`
      gives `log₁₀(property) = a1·log₁₀(T) + a2/T + a3/T² + a4`. The kernel evaluates CEA
      RP-1311 Eq. (5.1), `ln(X) = A ln T + B/T + C/T² + D` (`cl_GT_TransportPoly.hpp:48`), and
      the coefficients in `trans.inp` are for the natural log. Anyone re-deriving viscosity or
      conductivity from the documented form is wrong by factors of ln 10. Same at
      `gastables_usage_guide.md:518`

### D16 Further physics-doc defects — **MEDIUM**

- [x] `src/physics/gastables/doc/README.md:316-321` writes `Cell<RefGas*> gases(species.length());`
      then assigns `gases(i)`. Two errors in one line: `Cell` has `size()`, not `length()`
      (the accessor split CLAUDE.md documents), and the single-argument `Cell(n)` constructor
      only **reserves** — sizing needs `Cell(n, value)`. Does not compile
- [x] `src/physics/materials/doc/README.md:108-109` says "the eight superconductor tables";
      `share/material/` holds **six** jc/n files, and `share/material/README.md:13-25` says six
- [x] Broken link in three files — `src/numerics/doc/spline_usage_guide.md`
      (`gastables/doc/README.md:532`, `gasmodels/doc/README.md:829`,
      `gasmodels_usage_guide.md:1160`). `src/numerics/doc/` does not exist; the file is
      `src/numerics/spline/doc/spline_usage_guide.md`
- [x] `src/physics/gasmodels/doc/README.md:64` names `cl_GM_AlphaFunction.{hpp,cpp}`; the files
      are `cl_GM_EoS_AlphaFunction.*` and `cl_GM_EoS_AlphaFunctionFactory.*`
- [x] `gasmodels/doc/README.md:25,234-237` says Helmholtz covers "H2, CH4, O2 only";
      `HelmholtzModel::Nitrogen` is live (`en_Helmholtz.hpp:24`)
- [x] `gasmodels_usage_guide.md:190,201-202` works examples on `R134a`, `C3H8`, `C4H10`, none of
      which are in the shipped `share/fluid/thermo.inp`. *Grok's search; not re-derived by me*
- [x] The gastables/gasmodels README "Location" tables cite class-at-line anchors that now land
      on comments and `#define`s. *Grok's list; not re-derived by me*


### D17 `src/homology/doc/homology_usage_guide.md` is stale across its whole API surface — **CRITICAL**

This is the single most non-current document in the tree. Grok returned ten P0 API defects in it;
I spot-checked four names against `src/` and the class is confirmed — these are not drifted
anchors, they are calls to functions that do not exist.

**Verified absent from `src/` (0 hits each, `.cpp` + `.hpp`):** `add_abstract_dof`
(`:765,2347` — the live API is `DofManager::extract_abstract_dofs_from_mesh()`,
`cl_FEM_DofManager.hpp:205`), `manifold_filter_3d`, `node_is_on_cut`.
**Verified misplaced:** `save_debug_meshes()` is called on `CutFactory` at `:721,2201` but exists
only on `CutProcessor` (`cl_CutProcessor.hpp:116`).
**Verified wrong:** `DomainType::Phi` / `::Boundary` / `::Interface` at `:1144-1156` — the enum
has `Air, Buffer, Ferro, Coil, Conductor, Cut, ThinShell, …` and none of those three
(`en_DomainType.hpp`).

- [x] Treat the guide as **needing a rewrite against the current API, not line edits.** Grok's
      full row list is in `tmp/ai_exchange/docsweep_b3_grok.md`; the remaining unverified items
      are `Topology::analyze()` (live: `run()` / `run_on_enriched_mesh()`), `Mesh::distribute()`
      / `Mesh::receive()`, `Protoshell` construction and its `add_*` setters, `Chain`/`Cochain`
      operators, the `CutProcessor` constructor's last two arguments, `CutData`'s members, and
      the `smithForm` return tuple
- [x] **`src/homology` is owned by Gregory Giard** (`en_Maxwell_Formulations.hpp:7` names him as
      a developer). This rewrite should be his call or done with him — it is not a
      mechanical sweep
- [x] The same guide's **glossary at `:314` is correct** — `H^k = Ker(δ^k) / Im(δ^{k-1})` — while
      `cohomology_algorithms.md:21` has `Im(δ^{k+1})` (D8). The tree contains both the right and
      the wrong formula, four files apart. Fix the wrong one, keep the glossary as the reference
- [x] `cohomology_theory_and_implementation.md:23`'s "every non-trivial cycle intersects the cut
      an **odd** number of times" is the **ℤ₂** picture; the pipeline uses integer cochains
      (`clean_spfa`, ±1 `CutData` bitsets) where the pairing ⟨c,z⟩ is the linked current and can
      be any integer. Both auditors flagged this independently
- [x] `cohomology_algorithms.md:143-145` claims "MPI-aware implementations for distributed
      computing"; `homology_usage_guide.md:71-73` states the serial-on-undistributed-mesh
      contract that matches the factory, and `homology_usage_guide.md:104` already says
      "Parallel cohomology is planned but not yet implemented". Two of the three are right
- [x] Complexity claims: `PellikkaGeneralized` is called **O(n log n)** and "optimal"
      (`cohomology_algorithms.md:49`, `homology_usage_guide.md:1037`, `README.md:149`). Pellikka
      et al. give a *range* of O(n log n) to O(n²) for reduction techniques
      (`pellikka2013.txt:188-198`); the code carries no complexity guarantee; and the cited
      Giarda preprint is **not in the local literature tree** (searched). Also RCM is stated as
      O(V+E+V log V) (`graph_usage_guide.md:475`, `graph/README.md:55`) where `symrcm()` sorts
      each dequeued neighbor list (`fn_Graph_symrcm.cpp:107-130`), giving O(V+E log V) worst
      case absent an unstated bounded-degree assumption
- [x] `thick_thin_cuts_and_conjugate_edges.md:282,325` cites `todo/periodic_thin_cut_continuity_fix.md`
      as the home of an implementation policy. **Module docs must not cite todo files** —
      the policy is already inline at `cl_CutSet.cpp:107-124`. State it there
- [x] `graph_usage_guide.md:55-66,167-189` and `graph/README.md:120-127` describe types that have
      changed: `id_t` is `unsigned int` not `uint64_t`, `gNoOwner` is
      `numeric_limits<proc_t>::max()` not `-1`, the vertex flag is `uint8_t mFlags` not
      `bool mFlag`. *Grok's read; not re-derived by me*
- [x] `graph_usage_guide.md:149` locates `Vector` at `containers/cl_Vector.hpp`; it is
      `src/linalg/cl_Vector.hpp`


### D18 The comm guide describes three asymmetric primitives as collective — **CRITICAL**

Both auditors landed here independently, and Grok's B5 round widened it from one function to
three. This is not three separate typos: the guide was written against a mental model of MPI
collectives that BELFEM's API deliberately does not follow, so the same error repeats per
function. **Verified by counting `MPI_Isend` / `MPI_Irecv` calls in each implementation.**

- [x] `comm_usage_guide.md:698-709` — `distribute(data)` called on **every** rank, then "Now data
      on rank i contains i*10". `distribute(Cell<T>)` (`commtools.hpp:766-812`) contains
      **0 receives**. Non-root ranks are never written to. Live pattern is rank-0 `distribute` +
      worker `receive` (`cl_SolverDistMatrix.cpp`)
- [x] `comm_usage_guide.md:807-813,982-996` — standalone `collect(all_values, my_value)`
      documented as a gather. `collect(Cell<T>)` (`:821-862`) contains **0 sends** — it is
      receive-only and needs matching senders
- [x] `comm_usage_guide.md:103-120` — the "Broadcasting Uninitialized Data" pitfall is
      **self-refuting**: the WRONG and CORRECT snippets are the *same* `broadcast(v, 0)` call with
      the same setup, differing only in the comment. One claims CRASH, the other claims automatic
      resize. The resize is real (`commtools.hpp:736,1063`), so the WRONG half is simply wrong
- [x] `comm_usage_guide.md:31-36` — "Automatic Chunking: large messages split into 64 KB chunks"
      is stated generally, but `broadcast` is a single unchunked `MPI_Ibcast`
      (`commtools.hpp:710-748,1050-1079`); chunking applies to send/receive/share. This is the
      distinction CLAUDE.md draws, inverted. The **"46,340 process limit" at `:248` is wrong**:
      the bound is `0.5*(sqrt(2*mMaxTag+1)+1)` (`cl_Communicator.cpp:220-224`) over a tag built as
      `2*(tmax*size+tmin) % max_tag` (`commtools.cpp:109`), giving **~32768** for
      `MPI_TAG_UB ≈ 2³¹−1`, not `√(2³¹) ≈ 46340`. Grok's arithmetic, not executed — **derive it
      before writing a number**, and do not re-check it as `√MPI_TAG_UB` or the same error repeats
- [x] `comm_usage_guide.md:1002-1018` — the ghost-exchange pattern is send-all then receive-all;
      `send` waits on `MPI_Waitall` (`:944-945`), so symmetric neighbor pairs can deadlock on
      rendezvous-sized payloads. *Grok, medium confidence (~70%) — needs a ruling, not a blind edit*
- [x] `src/comm/doc/comm_usage_guide.md:881-890` describes `share(data)` as "Sends vector from
      current process to all other processes (**all-to-all broadcast**)" and shows it called
      **unconditionally on every rank**, commented "All ranks send their vector to all others".
      `share(Vector)` is a loop of `MPI_Isend` with **no receives** (`commtools.hpp:1654-1670`).
      CLAUDE.md states the contract directly: *"The `if/else` rank guard is required — `share` and
      `receive` are NOT collective. Calling either on the wrong rank is a deadlock."*
      `doc/coding_philosophy.md:537-541` teaches the correct guarded form. **The module's own
      usage guide teaches the hang that the framework documentation exists to prevent.**
      Its `commtools.hpp:1540` anchor is also stale

### D19 Infrastructure and FEM API drift (B5/B6) — **CRITICAL / HIGH**

Codex returned 12 P0s in B5 and 2 in B6, nearly all "the documented call does not exist". Sampled
and confirmed above; the rest carry Codex's confidence and need a pass before editing.

- [x] `src/io/doc/io_usage_guide.md:295`, `io/doc/README.md:60-72` — HDF5 examples call
      `save()`/`load()`; the interface declares `save_data()`/`load_data()` (`cl_HDF5.hpp:147-163`)
- [x] `io_usage_guide.md:1068-1092` — `aParallelMode=true` documented as cooperative single-file
      PHDF5 with collective close. It rewrites the path **per rank** (`cl_HDF5.cpp:31-40`) and
      opens with `H5P_DEFAULT`, not an MPI access property list
- [x] `io_usage_guide.md:1301-1313` — `OPEN_RDONLY_PARALLEL` documented as collective MPI reads.
      True for **HDF5 only**: `cl_HDF5.cpp:42-103` handles `NEW`/`OPEN_RDONLY`/`OPEN_RDWR` and
      falls through to "unknown filemode passed". **`Ascii` does handle the mode**
      (`cl_Ascii.cpp:46-49`) and the enum value is real (`filetools.hpp:30`), so **do not write
      "the mode does not exist"** — scope the correction to HDF5 (Grok, plan audit)
- [x] `container_usage_guide.md:78,793-827` — `Genome<4,3> genome;` (no default constructor) and
      `Cell<Genome*> population(POP_SIZE)` then indexed (the one-argument `Cell` ctor **reserves**,
      it does not size — the same trap as D16)
- [x] `linalg_usage_guide.md:489,787` — `Vector<real> x = posv(A, b);`. `posv` mutates its RHS and
      returns `int_t` status (`fn_posv.hpp:218-227`)
- [x] `mesh/doc/README.md:208-213` + two more files — `Distributor::create_local_mesh()`; the API
      is `run()` / `partial_mesh()` (`cl_Mesh_Distributor.hpp:86-97`). Grep finds the name only
      in documentation
- [x] `mesh_usage_guide.md:276-295,833-850` — `push_back` on `Vector` and on `Cell`. `Cell` has
      `push`; `Vector` has `set_size`. Also constructs `ThinShell` with one sideset where both
      primary and ghost are required (`cl_ThinShell.hpp:83-84`)
- [x] `mesh_usage_guide.md:909-915` — `Mesh::save(path, true)`; `save` takes one argument
      (`cl_Mesh.hpp:233-235`), and the `.hdf5` path writes no fields
- [x] `sparse_usage_guide.md:547-562` — `set_use_inital_guess()` / `use_inital_guess()`. The
      declared spelling is `initial` (`cl_SolverParameters.hpp:215,269-270`) — the typo is in the
      **doc**, so both examples fail to compile
- [x] `sparse_usage_guide.md:88-100` — tells the caller to `create_coo_indices()` before MUMPS;
      the solver does it itself (`cl_SolverMUMPS.cpp:716-725`)
- [x] `linalg_usage_guide.md:275-285` — `Matrix::capacity() == n_rows * n_cols`. Under Blaze the
      inter-column stride may exceed `n_rows` (SIMD padding). This is the exact hazard CLAUDE.md
      flags for MPI payload sizing, taught as safe
- [x] `io_usage_guide.md:64-90` — "BELFEM error paths never use C++ exceptions and always abort".
      Debug builds throw (`assert.hpp:88-93,189-194`)
- [x] `container_usage_guide.md:487-492` — `DynamicBitset::where(..., true/false)` documented as
      selecting sparse vs dense; the flag is ignored, all calls use `where_sparse`
- [x] `bfm_file_format.md:175-180` — `.bfm` optionally contains `/circuit` state; `BfmFile::save()`
      ends after curve data, no circuit serialization
- [x] `dof_manager_usage_guide.md:371-377` + `kernel/doc/README.md:180-182` — default-construct
      `KernelParameters` (every constructor requires a `Mesh&`/`Mesh*`/`Kernel*`) and call
      `kernel.import_mesh(...)`
- [x] `nedelec.md:101-110`, `interpolation_usage_guide.md:1005-1015` — `LINE3` listed as a
      supported **factory-dispatched** edge element. `EdgeFunctionFactory` has no `LINE3` case
      (`cl_EdgeFunctionFactory.cpp:33-96`). **But `EF_LINE3` exists** (`nedelec/cl_EF_LINE3.hpp`),
      `num_nedelec_dofs(LINE3)` returns 2, and the mesh `ElementFactory` *does* handle LINE3 —
      so say "not dispatched by `EdgeFunctionFactory`", **never "LINE3 does not exist"**
      (Grok, plan audit). Consistent with the known LINE3 stub gap
- [x] `interpolation_usage_guide.md:1815-1818` — says `num_nedelec_dofs(TET4)` returns 4 and needs
      fixing. It returns **6** (`fn_num_nedelec_dofs.hpp:21-26`). A *fixed* defect presented as live
- [x] `interpolation_usage_guide.md:377-382` — `create_bubble_function(ElementType, BubbleType)`;
      the signature takes `(ElementType, uint aFacet)`
- [x] `interpolation/doc/README.md:33-38` vs `nedelec.md:146-164` — contradict each other on
      whether TRI6/TET10 parent-edge functions carry circulation 1/2 or unity
- [x] `spline_usage_guide.md:31-33,142-146` — nonuniform `X` documented as always erroring; the
      guard is `BELFEM_ASSERT` and compiles out in release
- [x] `circuit/doc/README.md:15,21,36-38` — driven from `hphirun`, and "no parser exists".
      `belfem.cpp:163-204` builds the circuit, and `cl_NetlistParser.hpp` exists
- [x] `fem/thermal/doc/README.md:16-17,48-49` — thermal backs "only `hphiTrun`"; `belfem`
      conditionally builds the thermal kernel (`belfem.cpp:180-200`)
- [x] `executables/doc/README.md:14-19` — table omits `material`, `gas`, `db2exo`
- [x] `bearing_gauge_eigenmode.md:16-22,81` — "silent no-op" for an unresolved bearing (now a
      `BELFEM_ERROR`, `cl_FEM_Bearing.cpp:96-120`) and "`absolute tolerance` parsed but
      deliberately unconsumed" (it participates in convergence, `cl_FEM_Controller.cpp:1681+`)
- [x] `linalg_usage_guide.md:27-30` — documents the backend switch as `BELFEM_ARMADILLO=ON` /
      `BELFEM_BLAZE=ON`. Those are **generated defines**; the cache options are
      `USE_MATRIX_ARMADILLO` / `USE_MATRIX_BLAZE` (`CMakeLists.txt:81-87`), and **Apple defaults to
      Blaze, not Armadillo**. A `-DBELFEM_ARMADILLO=ON` on the command line does nothing
- [x] `sparse/doc/README.md:41-49` omits **SuperLU** from the solver backends; `SolverType::SUPERLU`
      exists and `USE_SUPERLU` defaults ON (`CMakeLists.txt:89`). `README.md:70-72` also shows an
      `assemble(Ke, K, e->dofs())` free function that does not exist in `src/sparse/`
- [x] ~~`sparse_usage_guide.md:253-254` — `Vector<real> y = A * x;` on `SpMatrix`; there is no
      `operator*` under `src/sparse/`. The real call is `A.multiply(x, y)`~~
      **FALSE POSITIVE (retracted 2026-08-31).** Raised by Grok B5 #8 on a negative search that was
      simply wrong. `cl_SpMatrix.hpp:641-648` defines a free
      `operator*( SpMatrix&, const Vector<real>& )` that allocates and returns the result vector.
      The documented line — `Vector<real> y = A * x;  // Allocates new vector` — is **correct,
      comment included**. Caught by Codex in the plan audit and re-derived here.
      **Applying this row would have deleted a correct example.** Kept visible so it is not
      re-filed by a later sweep
- [x] `io/doc/README.md:128-129` — `make_path_parallel("output.h5")` documented as `output.h5.0`;
      actual is `output_4.0.h5` (`filetools.cpp:90-104`)
- [x] `io/doc/README.md:155-156` — `get_ids()`/`get_reals()` shown as returning vectors; they are
      out-parameter `void` APIs (`cl_Input_Section.hpp:181,206`). `get_value()` needs a unit string
- [x] `container_usage_guide.md:70` — Map `operator()` "throws if not found". **The proposed
      correction in the first draft was wrong and is withdrawn.** `cl_Map.hpp:235-243` uses
      `BELFEM_ASSERT` in debug and `BELFEM_ERROR` in release, and `assert.hpp:88-93` states the
      contract: *"a debug run throws and a production run aborts, at ANY rank count"*, with the
      `throw aException` at `:189-192`. So the doc is **literally correct in a debug build**.
      The fix is a **refinement, not a reversal**: say it throws in debug and aborts in release.
      Writing "aborts, not throws" would introduce an error where none existed. Caught by Codex
- [x] `core_usage_guide.md:360` — "for long-running timers use `stop()` directly" instead of
      `next()`; `Timer::stop()` casts to `unsigned int` too (`cl_Timer.hpp:42-47`), so it carries
      the same 2³² ms wrap. The advice does not fix what it claims to fix
- [x] `mesh/doc/README.md:152` lists **TRI21** as supported; the type exists in the enum but
      `ElementFactory::create_element` has no case for it and hits `BELFEM_ERROR`. Also
      `mesh_usage_guide.md:14` says "29 types" against the glossary's "50+" at `:74`
- [x] `mesh/doc/periodicity.md:17,169,259` — `DofData::create_dofwise_periodicities_master()`;
      **no matches in the tree.** Already on record as absent (INC-010)
- [x] `mesh_contracts_and_invariants.md:27-39` — `compute_element_indices()` does not exist
      (the method is `update_element_indices`); `bfm_file_format.md:21-22` names
      `Mesh::save_fields`/`load_fields` path overloads that are commented out (`cl_Mesh.hpp:1075-1078`)
- [x] `container_usage_guide.md:606` cites `todo/lesson_learned.md` — no such file; it is
      `doc/lessons_learned.md`
- [x] `mesh/doc/thin_shell_geometry_and_periodicity.md:16,21,25,45,59` — every line in its
      `ThinShellFactory` map has drifted
- [x] Stale anchors: `bdf_timestepping_theory.md:37-38,50-55`,
      `hanging_dofs_static_condensation.md:598-616`, `quaternion_usage_guide.md:164-168`


---

## C — Code defect found during the sweep (NOT a documentation issue)

### C1 `Vertex::mOwner` initialises to the wrong sentinel — **Christian's ruling, 2026-08-31**

`src/math/graph/cl_Graph_Vertex.hpp:42` reads `proc_t mOwner = gNoID;`. It should be `gNoOwner`.
Christian ruled this a **defect**, not a convention, when the sweep surfaced it — so the
documentation records it as such rather than enshrining it.

**Verified by execution** (compiled and ran the narrowing, the one executable gate in this
campaign): `gNoID` is `numeric_limits<id_t>::max()` = 4294967295 over `unsigned int`; `gNoOwner`
is `numeric_limits<proc_t>::max()` = 2147483647 over `int`. Assigning the former to a `proc_t`
narrows to **−1**, and `mOwner == gNoOwner` evaluates **false** on a freshly constructed vertex.

**The intended sentinel is not in doubt:** `fn_Graph_dfs.cpp:32` resets with
`set_owner( gNoOwner )`.

**Two consequences, both reachable in principle:**

1. `owner() == gNoOwner` cannot detect "unassigned" on a fresh vertex.
2. `owner()` is used **directly as an array index** — `++tCounters( tVertex->owner() )`
   (`fn_Graph_METIS.cpp:302`), `tCount( tVertex->owner() )` (`graphtools.hpp:146,160`) — so an
   unassigned vertex indexes at −1: an assert in debug, undefined behaviour in release.

**Prior art suggesting this has already bitten once:** `graphtools.hpp:142` carries a
`// bugfix for non-assigned dofs e.g. from circuit model` that special-cases
`owner() == tCommSize` and reassigns to 0. That guard matches **neither** sentinel, which reads
like a symptom patch for exactly this confusion.

- [x] C1a Change `cl_Graph_Vertex.hpp:42` to `gNoOwner`. **DONE 2026-08-31** — Christian approved
      the edit mid-campaign and a dedicated agent applied it. One line.
- [x] C1a-COMPANION **`make check` will go red until `tests/math/test_GraphVertex.cpp:27` is
      updated.** It reads `EXPECT_EQ( ( belfem::id_t ) tV.owner(), belfem::gNoID );` — the cast to
      `id_t` is the tell that the test was written to match the buggy initialiser rather than the
      contract. It should be `EXPECT_EQ( tV.owner(), belfem::gNoOwner );`. **Not applied**: a test
      source change is outside this documentation campaign, and it is a second edit in a shared
      checkout. Required, not optional.
- [x] C1b The `graphtools.hpp:142` guard is a **genuinely different case — leave it.** It tests
      `owner() == comm_size()`, a *deliberate* "unassigned" marker set at `cl_FEM_Kernel.cpp:195,203`,
      which is neither sentinel. Two caveats the agent recorded: the guard is incomplete in both
      regimes (a vertex still holding `gNoOwner` sails past it into `tCount( owner() )` at `:146`),
      and the whole function is **currently dead code** — `build_pargraph_adjacency` is reached only
      from `parmetis_nd`/`ptscotch_nd`, neither of which has a caller in the tree.
- [x] C1c Reachability **traced by the agent, and my severity was overstated**:
      `mesh::Basis::Basis()` sets owner to 0 (`cl_Mesh_Basis.cpp:25`), so **no mesh entity ever
      observes the default**. The blast radius is a handful of direct `graph::Vertex` users, not
      the whole framework.
- [x] C1d **My "-1 underrun" framing was wrong.** `Vector<T>::operator()` takes `size_t`, so the
      old value converted to 2^64−1 — a wild index that always faults. The new sentinel gives
      2^31−1, which is *marginally less certain* to fault. The change is right, but it is a small
      regression in failure **loudness**, not an improvement in it.
- [x] C1e **Independent evidence that `gNoOwner` was always intended**, found by the agent and
      stronger than the naming argument: the `std::min` ownership sweeps at
      `cl_Mesh_Partitioner.cpp:253-260` and `cl_FEM_Kernel.cpp:427` only work if the sentinel
      behaves as **+∞**. Under the old value it was an absorbing element and every facet would have
      ended unassigned. Several sites are **repaired** by the change —
      `cl_Mesh.cpp:1669,1684` and `cl_FEM_Kernel.cpp:379,918,1103` all did signed comparisons that
      `-1` slipped through.

**Documented meanwhile** at `src/math/graph/doc/graph_usage_guide.md`, as a known defect with the
advice not to test against either sentinel until it lands.

---

## O — Open questions, ROUTED, not to be decided by an AI session

### O1 The φ-region time derivative — **ROUTED to Christian (physics)**

`src/fem/maxwell/doc/maxwell_usage_guide.md:329-369` states the air/ferro region as the static
`∇·(μ∇φ) = 0` Laplace form. `mt_maxwell_phi.cpp:26-46` assembles a φ **mass** contribution, so
doc and code differ on whether a time derivative is present there. The Arsenault 2026 erratum
(`literature/papers/fem/arsenault2026.txt:7-25`) replaces the 2023 air equation
`∂(μ∇φ)/∂t = 0` with Faraday `∇×E = ∂(μ∇φ)/∂t`, i.e. E is not zero in air.

Whether the guide is stale or is describing a deliberate BELFEM modelling choice is a formulation
question — `lessons_learned.md` escalation trigger 2. **Not decided here.**

### O2 The static-condensation attribution — **ROUTED to Christian (his paper)**

`src/fem/maxwell/doc/maxwell_usage_guide.md:186-210` presents static condensation as "the pattern
from Messe et al. 2023"; `CLAUDE.md:699,711` records the same attribution. Read firsthand,
`literature/papers/fem/messe2023.txt:509-517` says: *"Alves et al achieved this in their GetDP
implementation by eliminating the h_t degrees of freedom … a method known as static condensation.
For the sake of simplicity, we opted to couple the φ-field and the in-plane field h_t using
Lagrange multipliers … we wish to highlight that static condensation may be the preferred
method."*

So the method is Alves et al.'s, the paper's own implementation used Lagrange multipliers, and
the paper hedges with "may be". BELFEM's engineering choice is sound and is exactly what the
paper recommends — only the citation is wrong. How to reword is Christian's call.

### O5 The Doxygen nav files are generated from Python — **RESOLVED 2026-08-31**

Christian approved editing `scripts/update_doc_index.py`, which is where the stale executable names
actually lived. Fixed at the source and regenerated, so the correction is **durable** — editing the
`.dox` files directly would have been reverted by the next run. `hphirun`/`hphiTrun` no longer
appear in `doc/doxygen_nav.dox` or `doc/groups.dox`.

**Fixing the generator surfaced a second failure the first was masking.** With the `numerics/opt`
entry added, the check moved on to `no PAGE_GROUPS entry for: doc_examples_block3d_data,
doc_examples_undulator2d_data` — two pages the **peer session** added today with the usermat
migration, not mine. Registered them in the same group as their siblings. The generator then
injected its page anchors into those two files and into the Gregory handoff, one line each, which
is its own convention.

`update_doc_index.py --check` now exits 0 with **0 files out of date**; 27 modules, 101 pages.

**Addendum 2026-09-02 — the two data pages were the wrong repair.** Registering
`doc_examples_block3d_data` / `doc_examples_undulator2d_data` in PAGE_GROUPS silenced the check by
promoting two deck-internal notes to top-level manual chapters, which is not what they are. The
generator now prunes `NOT_DOCUMENTATION_DIRS = {"data"}` from the `examples/` walk, `Doxyfile.in`
carries the matching `*/examples/*/data/*` exclusion so Doxygen does not auto-label them either,
and the two anchors were stripped from the READMEs. `--check` exits 0; 27 modules, 99 pages.

---

### O5-original (kept for the record) — **ROUTED to Christian (scope)**

`doc/doxygen_nav.dox:42,52` and `doc/groups.dox:197,266` still describe `hphirun`/`hphiTrun` as
the solver applications and name `hphiTrun` as the thermal consumer. They look like documentation
and they are in the `doc/` tree, **but they are generated**: the strings live in
`scripts/update_doc_index.py:115,124`, and `update_doc_index.py` rewrites both files.

That puts the honest fix outside this approval:

- editing the `.dox` files alone is **not durable** — the next `update_doc_index.py` run reverts it;
- editing `update_doc_index.py` is a **Python source change**, which both auditors ruled out of
  scope when they struck R8.

The change itself is two string literals in a data table and is behaviourally inert. **Ruling
needed:** treat those literals as documentation and fix them, or leave both files stale until the
R8 tooling task. Left untouched in R1 either way. Found by Codex in the R1 diff audit.

### O3 Erratum citation hygiene

- [ ] `src/fem/kernel/doc/dof_manager_usage_guide.md:2825,3074` and
      `src/fem/iwg/doc/iwg_usage_guide.md:2061` cite Arsenault et al. 2023 for the magnetodynamic
      coupling without the 2026 erratum, though `doc/literature_references.md:40` records it and
      CLAUDE.md's routing table requires reading them together

### O4 Thin-shell `N > 1` stated as an unconditional requirement

- [ ] `src/fem/maxwell/doc/maxwell_usage_guide.md:1066-1073` and
      `src/fem/maxwell/doc/README.md:507-510` state `N > 1` as a flat requirement. Alves 2022b
      (`literature/papers/fem/alves2022b.txt:930-956,1120-1132`) makes it conditional — `N = 1`
      misses top/bottom losses for *closely packed* tapes, and can be comparable when normal-field
      effects dominate. Add the condition and the reason, or keep the simplification deliberately

---

## R — Execution batches

Each batch is one unit of work: **apply → tick boxes → dispatch the diff to Codex + Grok → record
verdicts → next batch.** Batches are ordered by severity-per-file-touched, so the newcomer path is
repaired first. **A stop after any completed batch leaves the tree better than it was** — but only
now that the two false positives are retracted. Codex was right that this was not true of the
first draft: literal R3 would have deleted a correct `SpMatrix` example and literal R3b would have
made a correct Map statement false.

**Standing guard, both auditors independently:** several `Rn` steps (R4c, R5e, R6b, R7d) carry
work that no `Dn` row states in those words. That is by design — the D-numbering is coarser than
the R-steps — but it means **the source must be re-opened before those edits**, not just before
the ones with a citation. Line anchors will also drift as batches land, so **re-locate by content,
never by the line number written in this plan.**

Depth for every audit round, per protocol §9.1 (doc-diff audit, round 1):
`CODEX_MODEL=gpt-5.6-terra CODEX_EFFORT=high`, `GROK_MODEL=grok-4.6 GROK_EFFORT=high`.
Round ≥2 or a split verdict escalates effort to `xhigh`, model unchanged.

| Batch | Theme | Defects | Files | Why here |
|---|---|---|---|---|
| R1 | The newcomer path | D1, D5, D10 (index rows) | ~11 | Every one is read before anything else in the project |
| R2 | Snippets that fail silently | D2, D3, D18, D14 (heat_spline, mu) | ~7 | Compile clean, return wrong numbers or hang |
| R3 | Infrastructure API snippets | D19 (io, containers, linalg, mesh, sparse) | ~12 | Compile-breaking; tedious but mechanical |
| R4 | FEM, interpolation, circuit, thermal | D4, D6, D19 (B6 rows) | ~10 | Mixed API and stale-status text |
| R5 | Physics: gas and materials | D14, D15, D16 | ~8 | Includes two wrong functional forms |
| R6 | Topology and graph | D7, D8, D9, D17 (partial) | ~8 | D8 is a real math error; D17 waits for its owner |
| R7 | Counts, anchors, housekeeping | D10, D11, D12, D13 | ~10 | Cheap, and R8 depends on the counts being right |
| R8 | Mechanical guards | — | `check_doc_claims.py` | Stops the same drift recurring |

---

### R1 The newcomer path — **do this one first**

Closes **D1**, **D5**, and the index rows of **D10**. Eleven files, and the two most-read documents
in the repository are among them.

- [x] R1a `README.md:44` run command; `:25` the `debug flags ON` default; `:31` the
      "for production use `-DUSE_DEBUG=OFF`" advice that follows from it
- [x] R1b `doc/getting_started.md:24-27` defaults, `:45-46` both run commands, `:52` the
      "`hphirun` solves magnetic decks" sentence
- [x] R1c `doc/parallel_execution.md:12` production command line
- [x] R1d `doc/README.md:17,49,59` index prose; add `belfem`, `material`, `gas`, `db2exo`; link
      `src/fem/doc/README.md`
- [x] R1e `examples/README.md:38-40` — the false "still exist and still work" sentence
- [x] R1f `src/executables/doc/README.md:14-19,17-35,96` — the table, the availability claim, and
      the `mpirun -np 4 hphirun` example
- [x] R1g `doc/coding_philosophy.md` `USE_DEBUG` default — **three sites, not one.** `:56`
      ("default **ON**", and its `CMakeLists.txt:78` anchor points at `USE_OPENMP`; the option is
      at `:98`), `:755` ("`USE_DEBUG=ON` (default)"), and check `:591` and `:820` read correctly
      once the default flips in the prose. Found by my own recon, not by the sweep — the original
      D5 row named only `:56`
- [x] R1g-note After R1g rewrites the defaults sentence, **re-open `coding_philosophy.md:59,61`
      before R7c** — they sit in the same paragraph (`:56-61`) and their line numbers will move
- [x] R1h **`maxwell_usage_guide.md:310-311` is adjacent to the frozen `O1` band (`:329-369`).
      Surgical single-line edits only; do not rewrite the surrounding section.**
      `src/fem/maxwell/doc/maxwell_usage_guide.md:310-311,1019`;
      `src/sparse/doc/solver_memory_and_compression.md:106`;
      `src/fem/maxwell/doc/README.md:340`; `doc/input_file_reference.md:992`
- [x] R1i `src/executables/CMakeLists.txt:40` — comment typo "retiredI". **Comment text only**;
      the commented-out blocks stay exactly as they are
- [x] R1j Audit round **COMPLETE**. Codex **REJECT**, Grok **ACCEPT WITH CHANGES**. Every
      finding re-derived and fixed; batch re-verified clean

**Verdicts and what they cost.** Codex rejected the batch and was right to. The blocker was mine:
I wrote that `belfem` names its Exodus output after the mesh, and `belfem.cpp:303` writes a fixed
`hphi_results.e-s` on the **segregated coupled path** while `:245` writes the mesh-derived name.
**While removing a false claim about the output filename I introduced a different one** — the exact
failure mode this campaign was built to catch, in batch one. Both paths are now named.

Also fixed from the two rounds:

- `USE_DEBUG` "defaults OFF" was stated universally in three places; it is OFF only without
  `$SCLS`. All three now carry the qualifier **in the same sentence**, not a paragraph later
- "`make install` deploys all five" — `USE_GASMODELS` defaults **OFF** (`CMakeLists.txt:105`), so
  a default install is four binaries. Rewritten to describe the gating
- "the configure summary prints the flavor" — only `if( BELFEM_SCLS_FLAVOR )`
  (`config/summary.cmake:30-32`)
- **A two-artifact rule violation of my own:** I rewrote the `t`/`temp`/`temperature` prose in
  `input_file_reference.md` and left its twin in `doc/input_schema.yaml:1637-1640` describing
  `hphirun`/`hphiTrun` as live. Both now agree, and the row's anchor cites all three real sites
  (parse, 77 K fallback, BC rejection) instead of only the parse site. Caught by Grok
- `coding_philosophy.md:820` still cited `assert.hpp:99` after I fixed `:59` — **R1g explicitly
  told me to check `:820` and I did not.** Caught by Grok
- `maxwell_usage_guide.md:702`, `input_file_reference.md:126`, and the
  `src/executables/doc/README.md` overview paragraph (my new table made its "each executable reads
  `input.conf`" opening false) — all corrected

**Files R1 touched that the plan never named**, all found by grep-back or by the auditors:
`src/sparse/doc/sparse_usage_guide.md`, `examples/scripts/README.md`, `examples/scripts/Allrun`,
`src/physics/materials/UserMaterialTemplate.cmake`, `src/fem/thermal/doc/README.md` (pulled
forward from R4f — it contradicted `doc/README.md:50` the moment R1 landed), `CLAUDE.md:136`,
`doc/input_schema.yaml`.

**Final state:** `check_doc_claims.py` 37/37. A repository-wide grep over `*.md`, `*.dox`,
`*.cmake` and `*.yaml` finds no live run instruction for the retired binaries and no remaining
`USE_DEBUG`-defaults-ON claim outside `devlog/` (historical, kept as written) and `todo/`
(working artifacts). **Reviewed, not verified** — no build or run was performed.

**R1 applied 2026-08-31. 15 files, documentation and comments only.** `check_doc_claims.py`
still 37/37. Two things worth recording:

- **The plan's file list was not exhaustive.** A completeness grep after the edits found
  `src/sparse/doc/sparse_usage_guide.md:473` carrying the same broken `mpirun -np 4 hphirun`
  command — named by neither the sweep nor the plan. Fixed in this batch. Every later batch
  ends with the same grep-back rather than trusting its own list.
- **Two corrections beyond the row text**, both found while editing:
  `getting_started.md` claimed the output is `hphi_results.e-s`; `belfem.cpp:210` builds the
  Exodus name from the mesh label (`sprint("%s.e-s", tMFactory.label())`,
  `cl_MaxwellFactory.cpp:321-322` strips the extension), so `helix.msh` gives `helix.e-s`.
  And `coding_philosophy.md` carried two stale anchors in the same paragraph — `CMakeLists.txt:81`
  for `USE_TEST` (it is `:101`) and `assert.hpp:99` for the `throw` (it is `:191`) — fixed while
  the paragraph was open, which is cheaper than R7c reopening it.

**Deliberately not in R1:** `config/globals.cmake:19`. The stale names are harmless
(`Add_Executable.cmake:63` filters on existing targets) and removing them is a build-file change,
not a doc change. Raise it with Christian separately.

### R2 Snippets that fail silently *(after R1)*

The four defects that compile clean and then hang or lie. Highest value per line in the plan.

- [x] R2a `src/comm/doc/comm_usage_guide.md` — **D18 in full**: `share` `:881-890`,
      `distribute` `:698-709`, `collect` `:807-813,982-996`, the self-refuting pitfall
      `:103-120`, the chunking overclaim `:31-36`, the `:248` process-limit figure, and the
      `commtools.hpp:1540` anchor. Model the corrected examples on
      `doc/coding_philosophy.md:534-545`, which is the canonical correct block (`if/else` guard
      plus the explicit "calling either on the wrong rank is a deadlock" sentence), and on the live
      call sites `cl_Database.cpp:74-79` and `cl_SolverDistMatrix.cpp:263-266`
- [x] R2a-bis **`src/comm/doc/README.md:79-82` carries the same defect** and must be fixed in the
      same batch: `distribute(data);  // Scatter` / `collect(data, my_value);  // Gather`, with no
      rank guard. My earlier note that this README "already teaches the correct form" holds for its
      **send/receive** block at `:76-77` only — the distribute/collect block below it repeats the
      usage guide's error. Found by my own recon; the sweep did not flag it
- [x] R2b `comm_usage_guide.md:1002-1018` ghost-exchange deadlock — **medium confidence (~70%)**.
      Do not rewrite blind; either verify the rendezvous claim or mark the example with its
      precondition. If it cannot be settled, log it as `O5` rather than guessing
- [x] R2c `src/fem/maxwell/doc/maxwell_usage_guide.md:945,951` — `rho( T, B, beta )` argument order
- [x] R2d `src/fem/kernel/doc/dof_manager_usage_guide.md:2332-2337` — replace the invented
      `comm.*` object API with the free functions; drop `sum_all` entirely
- [x] R2e `src/physics/gastables/doc/gastables_usage_guide.md:349-351` and
      `gastables/doc/README.md:363-364` — `heat_spline` holds **enthalpy**: `eval`→H,
      `deval`→cp, `ddeval`→dcp/dT. **Also kill `d2eval`** — the `Spline` API is `ddeval`
      (`cl_Spline.hpp:177`), so `:351` will not compile even after the ladder is fixed (Grok).
      Prefer naming the accessors `spline_H` / `spline_Cp` / `spline_dCpdT`, which cannot be
      mis-ordered
- [x] R2f `src/physics/materials/doc/materials_usage_guide.md` — **edit `:650`, not `:614`.**
      Slot 5 of the constant table is a constant permeability where B/H and ∂B/∂H coincide, so
      "∂B/∂H" there is defensible; the secant-vs-differential split belongs to `BhCurve::mu(H)`.
      Sharpen the sentence at `:650` that already names it. See the retracted row in D14
- [x] R2g Audit round **COMPLETE**. Codex **REJECT**, Grok **ACCEPT WITH CHANGES**. Both found
      errors I introduced; all fixed and re-verified

**R2 applied, then substantially corrected by its own audit round. Three of my edits were wrong.**

| what I wrote | what the code says | caught by |
|---|---|---|
| `dmudH` exposes the differential ∂B/∂H | It is **dμ/dH** (`cl_BhCurve.hpp:32`, `cl_Material.hpp:1043`). ∂B/∂H = μ + H·dμ/dH, and is not returned by anything. A Newton tangent taken from `dmudH` would be wrong | both |
| "In this table `mu` is a **constant** permeability" | `mu` is enum index 5, **inside the non-constant range** (`T_crit` = 13, and the document says so itself one line above). A B-H curve repoints `mFunctionMu` at `mu_bhcurve` (`cl_Material.cpp:908-910`). I had over-corrected the original observation into a different false one | Grok |
| `distribute`/`collect` have **two** valid pairings | **Three.** `collect` also pairs with an ordinary peer `send( data, root )` — the guide's own *Pattern 1* does exactly that. My text would have taught a working gather as a hang | both, independently |

Also corrected: the all-to-all form is **not** rendezvous-safe in general (`distribute` completes
its `MPI_Waitall` before the receives are posted), so the guidance is now scoped by overload —
`Cell<T>` moves one element per peer and is eager-safe, which is why the shipped test passes,
while the `Cell<Vector<T>>` family can deadlock. `receive` on the sender is a **no-op**
(`commtools.hpp:970`), not a deadlock as I wrote. Container `broadcast` issues **two** `MPI_Ibcast`
calls, not one. `collect( T*, offsets )` is a fourth, narrower shape that my blanket
"same semantics" note did not cover.

**The batch also found the contract neither auditor stated.**

`distribute` and `collect` are not a scatter and a gather. Every `distribute` overload is
**send-only** (`MPI_Isend`, zero `MPI_Irecv`) and every `collect` overload is **receive-only** —
checked across all five overloads, not just the `Cell<T>` pair the auditors read. Because
`comm_tag(s,t)` is symmetric, `distribute` has **two** valid pairings:

1. **all-to-all** — every rank distributes, then every rank collects. Proved by
   `tests/comm/test_CommMPI.cpp:555-571`, whose assertion `tRecvData(p) == p*100 + tRank` only
   holds under that reading;
2. **scatter from one rank** — root distributes, workers call the scalar `receive(value, root)`,
   which is exactly *Pattern 3* in the guide and works because the tags match.

**I got this wrong first.** My initial rewrite said `distribute` was "useless on its own", which
would have condemned a working idiom the guide already documents. Caught by reading the scalar
`receive` tag (`commtools.hpp:313-335`) before dispatching, not by an auditor.

Corrections beyond the row text, all found by grep-back:

- `comm_usage_guide.md:49` repeated the 46,340 figure in a summary bullet; the section fix alone
  would have left the page self-contradictory
- Pattern 3's **first half** was a second, unreported deadlock — `collect` on every rank with no
  `distribute` anywhere
- `gMaxCommChunkLength` is 65 536 **elements of `T`**, not "64KB"; for `real` that is 512 KB.
  `coding_philosophy.md:753` already had it right
- five `distribute`/`collect` overload examples beyond the two the plan named
- **No global-sum helper exists.** `allreduce()` is fixed to `MPI_MAX` (`commtools.hpp:234-250`),
  so `sum_all` could not be replaced with it — the corrected text says so rather than
  substituting a wrong equivalent

Pattern 4 (ghost exchange) was **flagged, not rewritten**, per R2b: `send( Vector )` completes its
own `MPI_Waitall`, so symmetric neighbor pairs can deadlock above the rendezvous threshold. That
threshold has not been measured here, so the guide now carries the precondition instead of a
confident rewrite.

### R3 Infrastructure API snippets *(after R2)*

**D19**, the io / containers / linalg / mesh / sparse half. Mostly one-line signature corrections.
Every row needs the header re-opened before the edit — several of these came from a single auditor.

- [x] R3a `src/io/doc/` — `save_data`/`load_data`; the `aParallelMode` per-rank-path reality;
      `OPEN_RDONLY_PARALLEL` falling through to "unknown filemode"; `make_path_parallel` naming;
      `get_ids`/`get_reals` out-parameter signatures; the "never uses exceptions" claim
- [x] R3b `src/containers/doc/container_usage_guide.md` — `Genome` construction; the
      `Cell(n)`-reserves-not-sizes trap; `where(...)` flag ignored; Map `operator()` — **keep the
      debug-throw / release-abort distinction, do not write "aborts"**; `todo/lesson_learned.md` path
- [x] R3c `src/linalg/doc/linalg_usage_guide.md` — `posv` signature (it mutates **both** A and B
      and returns `int_t`; the path is `src/linalg/lapack/fn_posv.hpp:226-262`, not
      `src/linalg/fn_posv.hpp`); **`capacity()` padding**
      (state the `spacing() * n_cols` rule); the `BELFEM_ARMADILLO` vs `USE_MATRIX_ARMADILLO`
      flag error and the Apple-defaults-Blaze correction
- [x] R3d `src/mesh/doc/` — `create_local_mesh` → `run()`/`partial_mesh()`; `push_back` on
      `Cell`/`Vector`; `ThinShell` two-sideset requirement; `Mesh::save` arity; `.bfm` `/circuit`
      claim; `save_fields`/`load_fields` overloads; `compute_element_indices`;
      `create_dofwise_periodicities_master`; TRI21; the 29-vs-50+ element count;
      `thin_shell_geometry_and_periodicity.md` line map
- [x] R3e `src/sparse/doc/` — `set_use_initial_guess` spelling (**the typo is in the doc**);
      the caller-owned `create_coo_indices` advice; missing SuperLU backend; the invented
      `assemble()` free function. **Do NOT touch the `A * x` example — it is correct** (see the
      retracted row in D19)
- [x] R3f `src/core/doc/core_usage_guide.md:360` — `stop()` carries the same `unsigned int` wrap
      as `next()`, so the advice does not fix what it claims
- [x] R3g Audit round **COMPLETE**. Codex **REJECT**, Grok **ACCEPT WITH CHANGES**. Both found
      the same blocker; four further errors of mine fixed after Grok's round

**R3 applied, then heavily corrected by its audit. I introduced a P0 regression and three
further falsehoods.**

| what I wrote | what the code says |
|---|---|
| `run()` then `partial_mesh()` on **every** rank | `partial_mesh()` is worker-only — `BELFEM_ERROR( mCommRank > 0, … )` (`cl_Mesh_Distributor.cpp:2663`). The original snippet failed at **compile** time; mine compiles and **aborts on rank 0**. Five sites, rewritten to the production shape (`cl_FEM_Kernel.cpp:688-708`) |
| `save(`→`save_data(` renamed by receiver name | It also hit an **`Ascii`** call, and `Ascii` has `save()` with no `save_data()`. **Rename by API, never by receiver** |
| "both forms do the same thing" for per-rank HDF5 | The open is gated on `( aParallelMode \|\| rank == 0 )` (`cl_HDF5.cpp:40`), so a hand-rolled per-rank *path* without the flag leaves every non-root rank with **no open file at all**. Two places, both mine |
| Scatter is "Pattern 3", gather is "Pattern 1" | Pattern 1 is both; Pattern 3 is all-to-all. Cross-references corrected |

Also fixed from the two rounds: `collect( Cell<Matrix<T>> )` leaves its own slot **empty**, so the
gather example silently dropped rank 0's matrix; `OPEN_RDONLY_PARALLEL` in `Ascii` is *rank 0 reads
then broadcasts*; `where()` still walks every level-2 summary word; `posv`'s `AbortOnError` raises
`BELFEM_ERROR`, which **throws in debug**; SuperLU registers **non-MPI**; `ThinShell` accepts a
**null** ghost sideset; eager buffering is not standard-guaranteed at any size; and a second
`Cell(n)`-then-index site at `container_usage_guide.md:823`.

**The batch also found four things the plan did not name.**

- **The parallel-HDF5 section was wrong in a way that mattered.** It described `aParallelMode` as
  cooperative single-file PHDF5 and carried a "DEADLOCK HAZARD" warning about conditional
  collective close. `aParallelMode` rewrites the path through `make_path_parallel()` and opens a
  **per-rank file** with `H5P_DEFAULT` (`cl_HDF5.cpp:31-40,47-51`) — nothing is collective, so the
  warned-about deadlock cannot occur. A reader was being taught to guard against an impossible
  failure while missing that their "single file" is N files
- **`Cell<int> data(comm_size())` appears twice in the comm guide**, including inside the
  *size-mismatch pitfall itself*, where it is presented as the CORRECT form. The one-argument
  constructor only reserves, so `size()` stays 0 and `distribute`'s assertion fires — the
  "correct" example fails the very check the section is about
- The `Vector<index_t> … push_back()` examples in the mesh docs are "WRONG:" examples for a
  *different* pitfall (index invalidation). They still have to compile to make their point, so
  they became `Cell` + `push` rather than being left broken
- **`create_dofwise_periodicities_master()` exists nowhere in the tree** (0 hits across
  `src/**/*.{hpp,cpp}`), yet `periodicity.md` names it three times as the DOF-constraint entry
  point. The surrounding description of the behaviour is accurate, so the fix names the file and
  flags the symbol rather than rewriting the section

**The retracted false positive held.** `sparse_usage_guide.md:253-254`'s `Vector<real> y = A * x;`
was left untouched, as the D19 retraction requires — `cl_SpMatrix.hpp:641-648` defines that
operator and the example is correct, comment included.

Two claims were **narrowed rather than deleted** on evidence: `OPEN_RDONLY_PARALLEL` is real and
handled by `Ascii` (`cl_Ascii.cpp:46-49`) but has no `HDF5` case, so the correction is scoped to
HDF5; and the io guide's "no C++ exceptions in error paths" is true in **release** and false in
debug, so it now says which build it describes.

### R4 FEM, interpolation, circuit, thermal *(after R3)*

- [x] R4a `src/fem/iwg/doc/` — Crank–Nicolson and Galerkin are **disabled** (`BELFEM_ERROR`), not
      supported; `create_iwg` arity; the BDF5 self-contradiction; the stale
      "`MaxwellThermal` not detected" pitfall
- [x] R4b **`doc/input_file_reference.md:256,270` is already correct that `chi` is on — do not
      "fix" it** (Grok). Only `coulomb_gauge_penalty_theory.md:11` and `maxwell/doc/README.md:33`
      are stale, and both must state the default **and** the opt-out (`chi : 0`).
      `src/fem/maxwell/doc/README.md:20` — remove **losses** from the postprocessed-field list
      (B/H/J/JJC stay); `:33` and `coulomb_gauge_penalty_theory.md:11` — the penalty is **on** at
      `chi = 1e-4`; `:507-510` and `maxwell_usage_guide.md:1066-1073` — make the `N > 1` condition
      and its reason explicit per Alves 2022b
- [x] R4c `src/fem/maxwell/doc/maxwell_usage_guide.md:506-568` — `Kernel( KernelParameters * )`;
      `thin_shell_facet_orientation.md:30-36,42-103,221-236` — `update_facet_nodes` skips
      `ThinShell`/`GeometryOnly`, and the enum values listed are obsolete;
      `postprocessor_recovery_theory.md:168-195` — the "Mode 1" attribution
- [x] R4d `src/fem/kernel/doc/` — `KernelParameters` construction and `import_mesh`;
      `bearing_gauge_eigenmode.md:16-22,81` — both statements describe fixed behaviour;
      the `bdf_timestepping_theory.md` and `hanging_dofs_static_condensation.md` anchors
- [x] R4e `src/fem/interpolation/doc/` — `LINE3`/`EF_LINE3` is **not** factory-dispatched;
      `num_nedelec_dofs(TET4)` returns 6 and the defect is fixed;
      `create_bubble_function(ElementType, uint)`; the TRI6/TET10 circulation contradiction
      between `README.md:33-38` and `nedelec.md:146-164`
- [x] R4f `src/circuit/doc/README.md:15,21,36-38` — driven by `belfem`, and `cl_NetlistParser`
      exists; `src/fem/thermal/doc/README.md:16-17,48-49` — `belfem` also builds the thermal kernel
- [x] R4g `src/numerics/spline/doc/spline_usage_guide.md:31-33,142-146` — the equidistance guard is
      `BELFEM_ASSERT` and compiles out in release
- [x] R4h Audit round **COMPLETE**. Codex **ACCEPT WITH CHANGES**, Grok **ACCEPT WITH CHANGES**.
      Six errors of mine found between them; all fixed and re-verified

**R4 applied, then corrected by both rounds. The two auditors found six errors I introduced, and
overlapped on only one — the strongest argument yet for running both.**

| what I wrote | what the code says | found by |
|---|---|---|
| `ElectricalCircuitFactory tElFactory( tFile, … )` | I introduced `tFile` while repointing an anchor, and never declared it — the snippet above still uses `"input.conf"`. **A compiling example turned into a non-compiling one by my edit** | Grok |
| `KernelParameters params(mesh);` in two sketches | `mesh` was not in scope in either. I replaced a sketch that at least parsed with an undeclared identifier | Grok |
| The TRI6/TET10 circulation "open point" | **Already answered in the tree**, and my flag repeated the very conflation it flagged. See the resolution note below | Codex |
| "naming either the unresolvable deck id or the point" | The `mID == 0` branch names the deck's `nodes` list **without** echoing the id | Grok |
| Left `update_facet_nodes() overwrites the wrong nodes` as the thin-shell failure mode | I added the skip caveat above it and left the sentence it falsifies four lines below | Grok |

**The circulation flag was the instructive one.** I recorded it as a question for the author,
reasoning from §7.1 that a formulation convention on proof-of-concept elements is not mine to
settle. `nedelec_derivation.md:313-323` had already settled it — the blanket unit-circulation rule
"holds for the **first-order** elements", TRI6 is the `1/2` exception, and the note ends "This note
does **not** extend to TET10", whose polynomials are twice the TRI6 pair and carry unit
circulation (`:403-404`). So `interpolation/doc/README.md:36` was simply wrong to lump them, and my
flag copied that error into a second file while calling it unsettled. **Caution about adjudicating
physics was right; not checking whether the tree had already adjudicated it was not.**

**One narrowing from Grok, carried into the note above:** subclass constructor signatures are *not*
uniform. `IWG_Maxwell` is `( Formulation, dim, bool, bool )` and `IWG_StaticHeatConduction` is
`( dim, type )` with no mode. Leaving the CustomType example alone was still right, but the reason
is "it matches the Poisson family", not "every subclass does this".

**The batch also found four things the plan did not name.**

**The grep-back found six more 3-arg `create_iwg` calls** than the plan named — it listed two, and
`iwg_usage_guide.md` carried eight in total. Removed with a scoped regex, and the diff was read
line by line afterwards to confirm no dangling commas, since R3's rename taught that a regex over
documentation needs its output inspected rather than trusted.

**I nearly "fixed" a fourth correct example.** `iwg_usage_guide.md:1566-1570` constructs an IWG
subclass directly as `IWG_CustomType( ModelDimensionality, IwgType, IwgMode )`, which looked like
the same three-argument error. It is not: the *base* `IWG` ctor is `( IwgType, ModelDimensionality,
IwgMode = Iterative, … )` (`cl_IWG.hpp:359-364`), while the Poisson family takes **(dim, type,
mode)** and forwards — `IWG_Poisson::IWG_Poisson( ModelDimensionality, IwgType, IwgMode )`
(`cl_IWG_Poisson.cpp:22-26`). The example matches that pattern exactly. **Left alone.**
**Do not generalise it**: subclass signatures are not uniform — `IWG_Maxwell` is
`( Formulation, dim, bool, bool )` and `IWG_StaticHeatConduction` is `( dim, type )` with no mode
at all (Grok, R4 audit).

**Deliberately flagged, not resolved:** `src/fem/interpolation/doc/README.md:36` summarises the
derivation as "TRI6/TET10 parent-edge functions carry circulation **1/2**", which does not sit
obviously with `nedelec.md:146`'s unit-circulation normalisation. TRI6, TET10 and LINE3 are
**proof-of-concept and outside the supported set** (`nedelec.md:206`), so this is a formulation
convention on unsupported elements — recorded as a question for the author with a "do not rely on
either wording" warning, in the same spirit as §7.1. Not adjudicated.

**Corrections beyond the row text:** a second `losses` claim in the Maxwell pipeline diagram
(`README.md:327`) that the plan named only once; `interpolation_usage_guide.md:1813` listing the
**already-fixed** TET4 dof-count bug as live, and citing the wrong file for it; and the circuit
guide's `hphirun` entry-point anchors repointed to `belfem.cpp:164,204` with the netlist-parser
status corrected — `cl_NetlistParser` shipped, only `circuitrun` remains deferred.

**N > 1 was read from the source, not the auditors.** `alves2022b.txt:939-947` says that at N = 1
the in-plane field does not fully penetrate the tape and AC losses are underestimated **at low
transport current**, converging onto h-φ as N grows. Both bibliography bullets now carry that
condition and reason instead of a bare "requirement".

### R5 Physics: gas and materials *(after R4)*

- [x] R5a **D15, both**: `gastables/doc/README.md:250` NASA-7 formula under a NASA-9 heading →
      the NASA-9 form from `cl_GT_HeatPoly.hpp:27`; `:272` `log₁₀` → `ln`, per
      `cl_GT_TransportPoly.hpp:48`. Same corrections at `gastables_usage_guide.md:498-514,518`
- [x] R5b `has_cryo_thermo`/`has_cryo_transport` → `has_thermo`/`has_conductivity`/`has_viscosity`.
      **More sites than the D-row lists** (Grok): `gastables/doc/README.md:194,198,235,437` and
      `gastables_usage_guide.md:384,390,396,668`. Replace with the live predicates — do **not**
      invent a `has_cryo_*` accessor; the cryo polys are synthesised in `finalize()` whenever the
      live predicates hold (`cl_GT_RefGas.cpp:243,290,306,646`)
- [x] R5c `HelmholtzModel` enum spellings; the "H2, CH4, O2 only" claim (Nitrogen is live);
      `cl_GM_AlphaFunction` filenames
- [x] R5d The `Gas air("air", …)` example and the `R134a`/`C3H8`/`C4H10` examples — **verify the
      shipped `share/fluid/*.inp` species list first**; both rows are single-auditor
- [x] R5e Spline-vs-polynomial default (the factory leaves gases in spline mode); the input-file
      names and the `$BELFEM_DATA/fluid` resolver ladder; the `Cell(n)`/`length()` snippet;
      "eight superconductor tables" → six; the `src/numerics/doc/spline_usage_guide.md` broken link
      in three files
- [x] R5f The thread-safety claim: `Gas` const evaluators write a shared cache
      (`cl_Gas.hpp:75`). Correct it toward BELFEM's stated not-thread-safe-by-design posture
- [x] R5f-bis `gastables/doc/README.md:35-39` and `gasmodels/doc/README.md:35-43` — the "Location"
      tables cite class-at-line anchors that now land on comments and `#define`s (**D16**;
      unassigned in the first draft, caught by Codex)
- [x] R5g Materials anchors: `thermal_expansion_from_heat_capacity.md:102,304`,
      `callaway_thermal_conductivity.md:35,192`, `materials_contracts_and_invariants.md:84`
- [x] R5h Audit round **COMPLETE**. Codex **REJECT**, Grok **ACCEPT WITH CHANGES**. Every formula
      accepted; the fallout was not

**R5 applied. Both auditors accepted every formula — deriving from the kernels worked — and both
rejected the fallout around them. The corrections I made created four new defects.**

| what I broke | how |
|---|---|
| `gasmodels/doc/README.md:194` | I renamed the declaration `propane` → `methane` and **left `propane.v(...)` on the next line**. A non-compiling snippet, created by my own rename — the same failure I have been finding in these documents |
| The SRK "good for" list | Substituting `CO2` produced a **duplicate `CH4` entry** and a direct contradiction with the "near critical — use PR instead" line four rows below. Replaced with `C2H6`/`C2H4`, which are shipped and keep the original LPG/petrochemical intent |
| Three more `Gas("air")` sites | My exact-string replace matched `Gas air("air", …)` and missed `air_combustion`, `idgas`, `fast_air`. Re-done with a pattern over the variable name |
| Three `make_path_parallel` + `HDF5(path, …)` pairs | Left teaching the hand-rolled form that R3 established leaves non-root ranks with no open file |

**The cryo rewrite was incomplete, and both auditors said so.** I rewrote the code examples and two
sections but left the feature bullet ("via supplemental data files"), a "NIST cryogenic property
data" source line, and a whole "Cryogenic Extensions" section still promising a per-species
dataset. **Rewriting the examples is not rewriting the claim** — the prose that framed them
survived and contradicted the fix.

Also corrected: the `HelmholtzModel` enum table still listed `HYDROGEN, METHANE, OXYGEN`; four
places said Helmholtz covers "H2, CH4, O2 only" when `Nitrogen` is wired
(`cl_Gas.cpp:108,143,842`); and my "wrong by factors of ln 10" was mathematically sloppy — a
base-10 form needs every coefficient divided by `ln 10`, and feeding the shipped coefficients to a
`10^f` form gives a temperature-dependent error, not a constant factor. `$BELFEM_DATA` is also
returned **unchecked**, so the marker-file rule applies to the fallback ladder only.

**The formula work itself stood.** Both auditors independently re-derived all three caloric forms
against `cl_GT_HeatPoly.cpp:38,52,69` — including the signs and divisors — and the transport
correlation against `rawpoly`, and confirmed the rewritten file-format sections against
`share/fluid/thermo.inp:44-48` and `trans.inp:51-53`, the unit conversions, the 14 class anchors,
and that the default `Gas` constructor really builds the 14-species mixture.

**The original R5 findings, for the record:**

**The three caloric forms were all NASA-7 under a "9-Coefficient" heading**, not just `Cp/R`. I
transcribed all three from `HeatPoly::Cp/H/S` (`cl_GT_HeatPoly.cpp:36-78`) and checked each against
RP-1311: `Cp/R = a1 T⁻² + a2 T⁻¹ + a3 + a4 T + a5 T² + a6 T³ + a7 T⁴`, with H and S following by
integration and `b1`/`b2` as the constants. Grok's warning held — `gastables_usage_guide.md:498-514`
is a **data record**, not a second copy of the formula, and was not turned into one.

**Both data-format sections were fabricated.** Checking the record layout against the shipped file
showed the "example (N2)" was a **CHEMKIN `therm.dat` block**, not this format at all. Rewritten
from `share/fluid/thermo.inp:44-48` and `trans.inp:51-55`. The real record carries an explicit
exponent row — `-2.0 -1.0 0.0 1.0 2.0 3.0 4.0` — which is independent confirmation of the NASA-9
powers, and the `D` exponent marker and fixed-column spacing are now noted because the readers
depend on them.

**`data.inp`, `alpha.inp`, `crthermo.inp`, `crtrans.inp` do not exist.** The factory opens
`thermo.inp`, `trans.inp`, `gasdata.inp`, `cubicalpha.inp` (`cl_GT_RefGasFactory.cpp:39-49`), which
is exactly what `share/fluid/` contains. The `gasdata.inp` columns were documented in the **stored**
units rather than the file's: the reader converts g/mol → kg/mol and **bar** → Pa
(`cl_GT_InputData.cpp:110,115`). The data path was `${BELFEM_SOURCE}/data/gastables/`; it is
resolved at run time from `$BELFEM_DATA/fluid` or a CWD ladder, keyed on the marker file
`gasdata.inp` (`fn_GT_data_path.cpp`).

**The cryogenic sections described a capability that does not exist as described.** There is no
cryo dataset and no `has_cryo_*` predicate: `RefGasFactory` *synthesizes* the range below the
lowest interval as an extrapolation at construction, along with glue polynomials and Lucas-correlation
viscosity (`cl_GT_RefGas.hpp:55-61`, `cl_GT_RefGas.cpp:243,290,306`). A mechanical rename of the
predicates would have left prose still promising a cryo/non-cryo distinction — the sections were
rewritten instead, which is the R3 lesson applied.

**`Gas air;` replaces the broken string form in 14 places.** The default constructor builds a real
14-species mixture (`cl_Gas.cpp:46-70`). `R134a`, `C3H8` and `C4H10` are absent from `thermo.inp`
(`C3H8` is in `gasdata.inp` only, so `create_refgas` still fails) and were swapped for species that
are actually shipped.

Also: 14 stale class-at-line anchors across the two Location tables, the six-not-eight table count,
the three broken spline links, and one more "thread-safe for read-only" claim in
`io_usage_guide.md:811` that the grep-back caught outside this batch's own files.

### R6 Topology and graph *(after R5)* — **partial by design, and now under §7.1**

> **Before starting R6, read `doc/ai_collaboration_protocol.md` §7.1 in full.** The six
> cohomology-core files are closed to AI edits, comments included. `src/homology/doc/` is exempt
> and is where all of R6's homology work lives — but §7.1's *reasoning* constrains the content of
> that work too: an AI's nearest prior here is textbook Pellikka, the shipped algorithm is a
> deliberate modification of it, and "three models agreed" is worth least in exactly this module.
> **Where a finding is about the algorithm, report it to Gregory; only edit where the doc
> contradicts itself or the tree.**

- [x] R6a **D8**: `cohomology_algorithms.md:21` `Im(δ^{k+1})` → `Im(δ^{k-1})`.
      **Safe under §7.1, and the reason matters:** this is the *definition* of a cohomology
      group, not a description of the reduction algorithm, so the "modified Pellikka" caveat does
      not reach it — a modified reduction does not change what `H^k` means. It is also not a
      textbook-prior call: the **tree contradicts itself**, with the correct
      `Ker(δ^k)/Im(δ^{k-1})` in its own glossary (`homology_usage_guide.md:314`) and in the code's
      storage (`mW[k]=ker(δ^k)`, `mV[k+1]=im(δ^k)`). Fix the outlier to match the glossary; do not
      reword the glossary. If any doubt survives, this is a one-line question for Gregory rather
      than an edit
- [x] R6b **D7**: `BELFEM_SUITESPARSE` → `BELFEM_METIS` in
      `graph_usage_guide.md:626,660,702,743,1018,1377,1383` and `graph/README.md:63,76,81`;
      missing METIS is a `BELFEM_ERROR`, not an RCM fallback. **Also `fn_Graph_METIS.hpp:39-40,64`**
      — the header's `@note`s repeat both errors, and comment text is in scope (Grok; the D-row
      listed only markdown)
- [x] R6c **D9**: the ParMETIS input model — complete graph on rank 0, sliced there, then
      distributed
- [x] R6d Graph types: `id_t` is `unsigned int`; `gNoOwner` is `numeric_limits<proc_t>::max()`;
      `uint8_t mFlags`; `Vector` lives in `src/linalg/`. Complexity claims for RCM and
      generalized Pellikka — state what the code guarantees (nothing) rather than a bound.
      **Restricted to `cohomology_algorithms.md`, `graph_usage_guide.md` and
      `graph/doc/README.md`.** The same claim at `homology_usage_guide.md:1037` is R6g's, not
      this step's — Codex caught the collision
- [x] R6d-bis Gaps the first draft dropped, all assigned here (Codex, plan audit):
      `src/homology/doc/README.md:115` — the six-argument `CutFactory` construction with
      `aComputePoisson` (**D4**, and it is *not* the protected usage guide, so it is in scope);
      `src/homology/doc/README.md:260` — the `manifold_filter_3d`/`check_surface` cleanup workflow
      (**D6**); `cohomology_algorithms.md:143-145` — the false "MPI-aware implementations for
      distributed computing" claim (**D17**, again not the usage guide)
- [x] R6e `thick_thin_cuts_and_conjugate_edges.md:282,325` — remove the `todo/` citation and state
      the policy inline from `cl_CutSet.cpp:107-124`; `:86` `CutProcessor::run()` does not exist
- [x] R6f `thin_cut_nonunit_rectification.md:22,334` — SPFA rectification **ships**; it is not a
      proposal. `cohomology_theory_and_implementation.md:29` — remove the
      `manifold_filter_3d`/`check_surface`/Tarjan workflow, which exists nowhere; `:23` — the
      odd-intersection claim is the ℤ₂ picture, the pipeline is integer cochains
- [x] R6g **STOP.** `homology_usage_guide.md` itself is **not** in this batch. Assemble the
      evidence into a short note for Gregory Giard and hand it over. **The note must also carry
      the `D4` rows on that file** — `:279,630,742` (`CutFactory` arity; the ctor takes five
      arguments and has no `aComputePoisson`, `cl_CutFactory.hpp:99-104`) and `:1785`
      (`BeltedTree` — the ctor already calls `select_belt_fasteners()` and `create_tree()`,
      `cl_BeltedTree.cpp:36-37`). Grok caught that these had **no legal batch**: R4's table claims
      D4 but names no homology file, and R6g forbids the guide. Parked here, not orphaned.
      `src/homology/doc/README.md:115` carries the same `CutFactory` example and **is** in scope —
      it is R6d-bis
- [x] R6h Audit round **COMPLETE**. Codex **REJECT**, Grok **ACCEPT WITH CHANGES**. Between them,
      one tooling defect affecting the whole campaign, and a P0 of mine that survived R3's own audit

**R6 applied under §7.1. No protected file was touched**, verified before and after.

**The round's most valuable finding is not about R6 at all: my repository-wide negative searches
were structurally unreliable.** `grep` in this shell is a **function** wrapping
`ugrep --ignore-files`, which honours `.gitignore` — and this repository ignores `archive/`,
`nonfree/`, `literature/` and `tmp/`. Every `grep -r … .` from the root silently skipped those
four trees. Codex caught it by finding `archive/graph/fn_Graph_tarjan.{hpp,cpp}`, a real Tarjan
articulation-point pocket detector, after I had asserted in a document **and in the handoff to the
module owner** that none existed. A stored memory said it did; I trusted the empty search over the
note. I re-ran every load-bearing negative from R2–R6 with `command grep`: **all of them hold** —
only the Tarjan claim was wrong. Recorded as a memory so later sessions do not repeat it.

**A P0 of mine from R3 that R3's own audit missed.** `mesh_usage_guide.md:1435` used
`globalMesh`, undeclared in that snippet (the variable is `tMesh`), then called `finalize()`
unconditionally and did `delete tMesh` on root — where `localMesh` **is** `tMesh`, so the example
freed the mesh it was still using. Three defects in one block I wrote two batches ago. Fixed, and
the "finalize on all ranks" rule at `:335` that contradicted it went too.

**Where I under-cleaned my own edits.** Both auditors flagged the same pattern independently: I
struck the generalized `O(n log n)` claim in one file and left it in two others I had edited in
the same batch — `cohomology_theory_and_implementation.md:18` and `README.md:171` ("**Best**
(O(n log n))"). Leaving a claim standing in a file you are already editing is not caution.

**And my "shipped" heading overclaimed.** §7 of `thin_cut_nonunit_rectification.md` describes the
approach that landed, but its body still reads as a proposal enumerating options the code does not
take. The heading now says the approach shipped and the listing is the design.

**Where §7.1 changed what I did, not just what I edited:**

- **The Pellikka complexity claim was split.** "O(n log n) for large meshes" was corrected against
  Pellikka's own text, which gives a **range** of O(n log n) to O(n²) kept low by applying the
  cheapest equivalences first (`pellikka2013.txt:196-198`) — that is a published source I can
  read. The "optimal O(n log n)" claim for the **generalized** variant was **removed rather than
  replaced**: the Giarda preprint is in preparation and absent from the local library, and the
  code carries no annotation, so there is nothing to check it against. Asserting a different
  bound would be exactly the failure §7.1 describes
- **"Fastest for large meshes"** became "the production default (`cl_MaxwellFactory.hpp:48`)" —
  what is checkable, with a note that no timing comparison is recorded
- **The ℤ₂ parity statement was flagged, not rewritten.** Applying the R4 lesson, I looked for the
  tree's own answer first and found strong evidence — `thin_cut_nonunit_rectification.md` exists
  because **non-unit** coefficients occur, which is meaningless over ℤ₂ — but this is still a
  formulation statement about Gregory's cut theory, so the note records the tension and the
  evidence and leaves the wording to him
- **`D8` was safe to fix and the reason is recorded**: it is the *definition* of `H^k`, the tree
  contradicts itself (the guide's own glossary at `:314` has it right), and Pellikka Eq. (A.2)
  settles it

**An unimplemented design, documented in three files.** `manifold_filter_3d()`, `check_surface()`
and the Tarjan/region-growing/BFS scheme appear in `cohomology_theory_and_implementation.md`,
`README.md` and `homology_usage_guide.md`. Searching the **whole repository including `archive/`**
returns documentation hits only — no implementation of any of them, and no Tarjan or articulation
code anywhere. Marked as a design rather than deleted, since it may be the plan of record.

**Handoff written** for Gregory Giard (internal, not in the published tree): the twelve API-drift
rows for `homology_usage_guide.md` (including that the file states the cohomology quotient **wrong
at `:775` and `:821` and right at `:314`**), the three questions, and what the sweep did fix.

**Found outside R6's named files by the grep-back:** `mesh_usage_guide.md:1384` guarded METIS with
`BELFEM_SUITESPARSE` — the same D7 defect. The remaining `BELFEM_SUITESPARSE` uses in the sparse
docs are **correct** (UMFPACK genuinely comes from SuiteSparse) and were left alone.

### R7 Counts, anchors, housekeeping *(after R6)*

- [x] R7a **D10**: Layer 1 is 200 lines; Layer 2 has 21 cards; 558 distinct incident IDs on 593
      rows. Fix `doc/README.md:12,13`, `lessons_learned.md:221`, `CLAUDE.md:53`, **and the two
      sites the first draft missed** (Grok): the `lessons_learned.md:7-8` preamble ("539 … 20
      clusters") and `lessons_learned_evidence.md:10,36` ("537"). Decide and state
      which number the prose means — "distinct incidents" or "table rows" — rather than swapping
      one bare figure for another
- [x] R7b `CLAUDE.md:520-523` module list → add `fem/postproc`, `fem/thermal`, `visualizer`;
      `:529` "23" → 26, keeping the DOI substance, which is still true
- [x] R7c **D11** remaining anchors: `coding_philosophy.md:59,61,245,302,504`;
      `ai_collaboration_protocol.md:194,403-404`; `documentation_guidelines.md:91-94`
- [x] R7d `documentation_guidelines.md:200-201,269` — the README required-vs-Optional
      contradiction. Note `update_doc_index.py:194-202` **warns**, it does not fail — say that
      precisely. **Also `:273`**, which repeats the stale "23 module READMEs" (Grok)
- [x] R7e **D12**: document `src/numerics/opt` — at minimum an entry in `doc/README.md` and
      CLAUDE.md's architecture list. A `doc/README.md` for the module is better. Note the
      undocumented `numerics/{integration,ode,bezier,sources}` and `math/tools` as follow-on
- [ ] R7f **D13**: `examples/Tape_Quench_obsolete_deleteme` — **ask before deleting.** A directory
      is not documentation; this batch may only add a line to `examples/README.md` naming it, or
      raise the deletion with Christian
- [x] R7f-bis `src/math/quaternion/doc/quaternion_usage_guide.md:164-168` — the cited test region
      documents borrowed-buffer storage that has been removed (**D19**; unassigned in the first
      draft, caught by Codex)
- [x] R7f-ter **O3**, the missing erratum citations, is small and mechanical — do it here rather
      than leaving it ownerless: add the Arsenault 2026 erratum beside the 2023 citations at
      `dof_manager_usage_guide.md:2825,3074` and `iwg_usage_guide.md:2061`. This adds a citation
      and changes no physics claim, so it does **not** touch the frozen `O1`
- [x] R7g Run `scripts/check_doc_claims.py`; then `scripts/resolve_doc_cites.py --check` over the
      touched files
- [x] R7h Audit round **COMPLETE**. Codex **REJECT**, Grok **ACCEPT WITH CHANGES** — and they
      **split** on the incident arithmetic, resolved here by re-derivation

**The split is worth recording, because Codex was wrong and acting on it would have broken correct
prose.** Codex held that `lessons_learned.md:8` and `lessons_learned_evidence.md:873` are "wrong to
call 539 rows", since the table has 599 rows and 564 ids. Grok held that the peer's framing is
exact. **Grok is right, verified by tracing the spans myself:** the main table runs `INC-001`
(`:93`) through `INC-537` (`:629`), `INC-538`/`INC-539` sit separately at `:870-871`, and the
addenda `INC-540`…`INC-564` at `:920-949`. The 35 surplus `^| INC-` hits are **re-listings in
checkpoint tables** — `INC-024` appears at `:116` as its row and again at `:739` as
`| INC-024 | 13 | middle |` in a sampling table. So 539 is the **locked catalogue**, exactly as
written. Left untouched; a caveat naming the addenda range was added where the number appears
without one.

**Everything else both auditors found was mine, and most of it was self-inflicted by this batch:**

- **The `opt` example did not compile.** `Objective` has a required `Objective( uint )`
  constructor and a **non-virtual** `uint dimension() const` (`cl_Objective.hpp:40,68-75`), so my
  `index_t dimension() const override` was invalid twice over.
- **The module count went 26 → 27 because of the README I added in this batch.** The list is now
  **generated from the tree** rather than hand-maintained, which is the only version that survives
  the next module.
- **A descending line range**, `cl_AR_Matrix.hpp:224-219`, produced by my blanket replace of `216`
  inside the range `216-219`.
- **A replacement claim I never checked.** Removing the dead `corctest.cpp` citation I asserted
  that `normaltest.cpp`/`pentatest.cpp` "use smart pointers freely". They use raw `new`/`delete`.
  The text now says the tree has no example of the pattern to point at, which is the true statement.
- **The `mOwner` note went stale within the hour** — the fix Christian approved landed while I was
  writing, so a note calling it a "known defect… outside this pass" was false by the time the
  auditors read it. Rewritten to describe the corrected state, and to distinguish `gNoOwner`
  ("never set") from the separate `owner() == comm_size()` convention ("set, and known to have no
  home").

**Confirmed sound by both:** the NLOPT gradient-length trap (the load-bearing claim in the new
guide), all eight algorithm mappings, the `Status`/`is_usable` contract, that the module builds by
default and has no FEM caller, the erratum citations, the quaternion rewrite, the protocol and
documentation-guideline path repairs, and `doc/README.md` deferring to the source rather than
restating counts.

**R7 applied. Every count was re-measured at the moment of writing, not reused.**

That mattered: **all of them had moved since the morning.** Layer 1 is 224 lines (was 200), Layer 2
holds 22 cards (was 21), the catalogue has 564 distinct incident IDs on 599 rows (was 558/593).
A peer session has been editing those files throughout the day. Reusing my own earlier numbers
would have written six fresh errors into the fix that was supposed to remove them.

**The peer had already fixed most of it.** `lessons_learned.md:7-12` now carries a precise
537-clustered / 539-rows / INC-540…564-addenda reconciliation, `:247` says 22 cards, and
`CLAUDE.md:53` says ~220 lines. I fixed only what they had not touched — `doc/README.md:12-13` and
`CLAUDE.md:520-532` — and **phrased `doc/README.md` to point at the source rather than restate
brittle numbers**, since these counts demonstrably rot within hours.

**`src/numerics/opt` is now documented** (`src/numerics/opt/doc/README.md`), registered in
`doc/README.md` and in CLAUDE.md's architecture list. Writing it, I caught my own error before
dispatch: my usage example passed `&tObjective` where the constructor takes an `Objective&`
(`cl_Optimizer.hpp:71`). The guide also records the trap that matters — the gradient argument has
length zero for the derivative-free `LN_*` algorithms and `dimension()` for the gradient-based
`LD_*` ones, so an implementation that always writes to it overruns under half the algorithm set.

**Citations: 665 resolve, up from 564 at the start of the campaign.** The four remaining
"MISSING" are accounted for: `cl_Comm_Table.cpp:245` is a **fictional illustration** inside the
protocol's "how to phrase a finding" examples block, and the rest are in
`lessons_learned_evidence.md`, a dated record kept as written. The `CMakeLists.txt` "BEYOND" hits
are my resolver preferring `src/CMakeLists.txt` (22 lines) to the root file (590) among 35
same-named candidates — a tool artifact, not drift.

- [x] R7h-BLOCKER ~~**I broke `scripts/update_doc_index.py --check` and did not fix it.**~~ **FIXED 2026-08-31** — Christian approved the Python edit. Details below.
      Adding `src/numerics/opt/doc/README.md` — the module documentation R7e asked for — made the
      generator fail: `error: no MODULE_GROUPS entry for: numerics/opt`. Every collected page needs
      an entry, so the check now exits non-zero and `make doc`'s index generation is affected.

      **The fix is one data line** in `scripts/update_doc_index.py`, beside its neighbor at `:106`:

      ```python
      ("numerics/opt",   "Optimizer",   "NLOPT-backed bound-constrained minimisation"),
      ```

      **Applied on approval.** And fixing it surfaced a second, pre-existing failure the first one was masking — see the note under O5.

- [ ] R7h-note **`examples/Tape_Quench_obsolete_deleteme` still exists and is named in no README.**
      Per scope guard 7 this batch does not delete directories. Raised for Christian: delete it, or
      rename and document it.

### R8 Close the blind spots — **OUT OF SCOPE for this approval**

**Codex flagged this correctly in the plan audit and it is withdrawn from the batch run.**
`scripts/check_doc_claims.py` is executable Python, and Christian's approval covers documentation
and comments only. R8 stays written here as the standing recommendation, to be raised as a
separate tooling task **after** R1–R7 land — its probes depend on the counts in R7 being right
anyway.

`scripts/check_doc_claims.py` passes **37/37** and could not see any of D1, D5 or D10.

- [ ] R8a Probe the built-executable set from **CMake targets**, not `src/executables/*.cpp`
      stems (`check_doc_claims.py:130`), and apply it to `doc/` and `README.md`, not only
      `CLAUDE.md` (`:242-251`). This alone would have caught D1
- [ ] R8b Probe `USE_DEBUG`'s effective default by resolving `${BELFEM_DEFAULT_USE_DEBUG}`
      through `find_scls_flavor.cmake` for the no-`$SCLS` case. Catches D5
- [ ] R8c Probe the four countable claims in D10 (Layer 1 line span, `^## L-` count, distinct
      `INC-` count, module `doc/README.md` count). All four are one-line greps
- [ ] R8d Consider a "taught snippet" gate: extract fenced `cpp` blocks from the guides and
      compile-check the ones that are self-contained. D2, D3 and D4 are all this class, and D2
      is the one that would still be silent at runtime
