# DR-18x Cluster: Five Ruled Register Rows Landed as Five Small Rounds

**Date:** 2026-09-24
**Purpose:** Land the five register rows Christian ruled on 2026-09-24 — DR-184 (QUAD16 facet
table), DR-182 (periodic tag over a typed sideset), DR-186 (thin-shell recovery guards),
DR-187 (thermal group admission), DR-166 (graph counter limit) — each as its own diff with its own
gate, in one session, the way the mesh one-liner bundle went.
**Module:** `src/mesh`, `src/fem/maxwell`, `src/fem/thermal`, `src/math/graph`, `tests/mesh`,
`tests/fem`
**AIs involved:** Claude Fable 5.1 (read, probe, plan, implementation, gates), Codex
`gpt-6-astra`/high + Grok `grok-4.7`/high (plan jury and code jury; astra because four of the five
rows are `src/mesh` or a `src/fem/maxwell` factory)
**Status:** ✅ CLOSED 2026-09-24 — all five rows landed in the working tree: R1–R5 behind plan jury ( astra/high + grok-4.7/high, 12 findings ), code jury ( 8 findings, 1 refuted, P1 on DR-186 ) and a narrow round 2 on DR-186 at astra/xhigh + grok-4.7/xhigh ( approved, one message refinement ). Executed: scratch facet-table 8/8 and periodicity 10/10 red-then-green, thermal admission 2/2, four DR-186 deck legs. Owed: `make check` in both trees ( Christian ), a 2-rank deck run for the DR-186 lookup, the commit. Register rows DR-184 / 182 / 166 / 187 / 186 carry the outcome; the DR-187 thermal-sideset gap and the DR-186 residues are on their rows.

> **Scope guards:**
> - Each row is its own diff and its own commit; the jury reads them as one thread.
> - No kernel formulation is touched; DR-168 runs as a separate physics round.
> - The cohomology core is not touched.
> - Every claim below about what the tree does today was read from source on 2026-09-24; the
>   facet-table claim was additionally executed (§1.1).

---

## 1. Current Behaviour and How It Fails

### 1.1 DR-184 — `Element_QUAD16::get_nodes_of_facet( 1 )` returns node 8 for node 7

`cl_Element_QUAD16.hpp:61-67` lists `{ 1, 2, 6, 8 }` for facet 1, the ξ = +1 edge. The
interpolation puts node 7 at ( 1, +c ) and node 8 at ( c, 1 ) (`lagrange/cl_IF_QUAD16.hpp:75-80`),
so node 7 is on no facet and node 8 on two. The 2-D `get_nodes_of_edge` forwards to the facet
table (`cl_ElementTemplate.hpp:909-913`), so edge 1 carries the same defect.

**Executed evidence (probe, 2026-09-24, release library `build/lib/libbelfem.a`):** a scratchpad
oracle placed every Lagrange element type's nodes at their parameter coordinates and compared
`get_nodes_of_facet( f )` against the set of nodes inside the convex hull of
`get_corner_nodes_of_facet( f )`, and `get_nodes_of_edge( e )` against the nodes on the segment of
its first two. Over 32 types (LINE2…LINE5, TRI3…TRI15, QUAD4…QUAD16, TET4…TET35, PENTA6/15/18,
PYRA5/13/14, HEX8/20/27/64, and the six TS/TB types): **exactly one defect, QUAD16 facet 1 and edge
1** (`table { 1 2 6 8 }`, `hull { 1 2 6 7 }`). HEX64, the other cubic table, is clean. So the fix is
one token and the scope does not widen. Confidence high (executed).

By-catch, not a defect: `to_string( ElementType )` returns `"unknown"` for QUAD9TS and PENTA18TS
(`Mesh_Enums.hpp:234`), which only affects messages.

### 1.2 DR-182 — `tag_periodic_sidesets` stamps `Periodic` over a deck-typed sideset

`cl_Mesh_PeriodicityFactory.cpp:1744-1782`: every sideset whose nodes all lie on the master or the
slave plane gets `set_domain_type( DomainType::Periodic )`, whatever it was. It runs from
`MaxwellFactory::create_periodic` (`cl_MaxwellFactory.cpp:645`) at constructor line 134, **after**
`read_domain_types` (line 112, `:561-590`) has applied the deck's `topology` sections, so a sideset
the deck declares `air symmetry`, `conductor antisymmetry`, `cut`, `tape`, `curve` or
`background field` and that happens to lie on a periodic plane loses that type. Consumers that
then miss it: the magnetic sideset dispatch has no `Periodic` arm
(`set_block_types_in_magnetic_equation`, `:2216-2245`, default arm "pass"), and the homology puts a
`Periodic` sideset into its periodic list (`cl_Topology.cpp:493-498`). Nothing reads the old type
again. Confidence high (source trace; not reproduced on a deck — no shipped deck has the shape).

The only types a sideset can carry at tag time are the deck's own (`en_DomainType.cpp:domain_type`,
which refuses unknown strings, `cut` among them on purpose at `:127`) and `Periodic` itself. Of
those, a `tape` / `thinshell` sideset is restamped `ThinShell` by `read_thin_shell_data` (`:3294`,
from `create_cuts` `:1121`) before any consumer runs, and a `background field` sideset named by a
background condition is rewritten at `:179`; the types with no later writer, and so the ones the
stamp really loses, are `air/buffer/conductor/ferro symmetry` and `antisymmetry`, `curve`, and a
topology `background field` without a background condition (jury, confirmed). `Cut` (`:1202`) and
the `Inactive` stamps (`:1672-1686`, gated on interface and thin-shell types) never touch a
`Periodic` sideset. Residual, no code: a `.bfm` reload restores sideset types
(`cl_ProtoMesh.cpp:1356`) and the constructor still calls `create_periodic` when the mesh carries
no periodicity (`:130-134`), so a deck turned periodic after its cache was written would now be
refused if a restored typed sideset lies on the plane; recorded on the row.

### 1.3 DR-186 — both recovery admission guards are commented out

`cl_MaxwellFactory.cpp:2891-2917`, inside `assign_materials` (rank 0, per thin shell, per facet,
per volume side): a `/* … */` block holds two `BELFEM_ERROR`s — the conductor neighbour must be
linear (the recovery gathers one edge dof per edge into a vector sized for all Nedelec dofs) and
nonmagnetic (`is_constant( mu )` and equal to μ₀). The comment says they "throw on Lawrencium
without cause" and names neither deck nor predicate; it arrived with the release commit `f81dd84`.

What is live around it: the order-1 guard `BELFEM_ERROR( max_element_order() == 1 )` at `:1576`
runs at `:721` only under `mCommRank == 0 && mComputeCohomologies`, and `assign_materials` runs
later at `:1003`; so on a deck without cohomologies a higher-order conductor reaches the loop
unrefused. On the material side every builtin metal, alloy and HTS is constant μ₀
(`cl_Material.cpp:165`, `cl_Material_Alloy.cpp:204-217`), and so is a builtin `Ferromagnetic`
that never loaded a curve: `set_have( mu, false )` (`cl_Material_Ferromagnetic.cpp:23`) clears
only the have-bit, while `is_constant` tests the constant slot (`cl_Material.hpp:1669-1671`). The
slot goes NaN in `reset_dependencies` (`cl_Material.cpp:237-241`), reached by `load_bh_curve`
(`:951-958`, the factory's `curve` branch `cl_MaterialFactory.cpp:59-66`) and by `set_custom`
(`:402-404`), which the T-only user-material overload calls (`cl_Material_UserDefined.cpp:60-84`).
So `is_constant( mu )` is false for a B-H-curve material and for a plugin that registers `mu`
through a function, even a constant μ₀ one. No shipped plugin registers `mu`. Which predicate
fired on Lawrencium is unknown; the μ one is the one a correct deck can trip. Confidence medium on
the Lawrencium mechanism, high on the code (jury, confirmed).

### 1.4 DR-187 — the thermal dispatch: one arm unreachable, one loop dead

**Premise correction for the 2026-09-24 ruling.** The ruling assumed both `default:` arms carry
by-design exclusions (Coil, Inactive, GeometryOnly, Ghost, Cut blocks; interface, symmetry,
periodic, terminal sidesets). The tree says otherwise:

| arm | what reaches it today | evidence |
|---|---|---|
| blocks, `cl_ThermalFactory.cpp:194-208` | only the blocks `create_thermal_kernel` selected, and it selects exactly `Conductor`, `ThinShell`, `Ferro`, `Buffer` (`:75-135`, both the count and the fill loop); the dof manager builds blocks from that selection (`cl_FEM_DofMgr_BlockData.cpp:76`). The `Air` arm at `:198` is dead (never selected) and the `default:` arm is **unreachable** | source trace, high |
| sidesets, `:211-227` | **nothing.** `IWG_MaxwellThermal` never selects sidesets (its constructor `cl_IWG_MaxwellThermal.cpp:21-37` and no caller of `select_sidesets`/`set_sidesets` in `src/fem/thermal`), and `DofManager` builds sideset groups only when `selected_sidesets().length() > 0` (`cl_FEM_DofManager.cpp:103-105`). The loop iterates an empty container | source trace, high |

Two consequences the row must carry:

1. The `Dirichlet` / `Neumann` arms of the sideset switch could never match anyway: they test the
   **mesh** sideset's type, and no code path types a mesh sideset `Dirichlet` or `Neumann` — the
   deck string map has no such word (`en_DomainType.cpp`), and `impose_dirichlet` / `impose_neumann`
   / `impose_alpha` set the type on the **fem** `SideSet` (`cl_FEM_SideSet.cpp:365`, `:380`, `:407`).
   The enriched integration order `tOrder` (`:210`) is set nowhere else, so it has never been
   applied to a thermal boundary sideset.
2. A deck with a thermal Dirichlet condition on a sideset crashes today on a null dereference:
   `SideSetData::sideset` returns `mEmptySideset` on a miss
   (`cl_FEM_DofMgr_SideSetData.hpp:152-161`), a dummy allocated only inside `create_sidesets`
   (`.cpp:163-171`), so with no selection it is still `nullptr` (`:57`) when
   `cl_FEM_PhysicalBoundaryCondition.cpp:266` calls through it. A thermal Neumann condition is
   refused at parse (`cl_ThermalBoundaryConditionFactory.cpp:89`), any other sideset type by the
   default arm (`:91-101`). No shipped deck has a `boundary conditions / thermal` section (the
   three coupled decks — `corc_twolayer`'s thermal lines are comments — are adiabatic), which is
   why nothing has shown it. **Thermal
   sideset conditions are unsupported in the tree**, not silently inactive — a bigger gap than the
   row records; it is recorded on DR-187 (extend, don't branch) and its support is not in this
   round.

So "refuse an unexpected type by name" reduces to a tripwire on the block arm, and the sideset loop
is dead code whose design (switching on mesh types) could not work if it were live.

### 1.5 DR-166 — hard-coded `uint32_t` limits on the graph counter guards

`cl_Graph_Vertex.hpp:305` and `:314` compare `mVertexCounter` against
`std::numeric_limits< uint32_t >::max()`; the member's type is stated once, in the class, and the
guards would not follow a width change. Both are per-adjacency (hot path); ruled: stay
`BELFEM_ASSERT`, use `decltype`, document the release wrap.

**Bottom line:** five independent diffs (each a source file, and for three of them a test and a CMake line); DR-184's scope is proven by execution to be one
token, DR-187's ruling needs the premise correction above before its code is written.

## 2. Architecture

No new mechanism anywhere. Each fix is the smallest expression of the ruling in the file that owns
the contract, with a test that observes the contract where one can be built without a kernel.

## 3. Gap Table

| # | row | change | gate | class |
|---|---|---|---|---|
| 1 | DR-184 | `mNodes[ 8 ]` → `mNodes[ 7 ]` on facet 1 | new `tests/mesh/test_FacetTables.cpp`: the §1.1 oracle as a gtest over all 32 Lagrange types, facets and edges; RED on HEAD at QUAD16, GREEN after | (c) |
| 2 | DR-182 | refuse a sideset on a periodic plane whose type is neither `Default` nor `Periodic`, naming sideset, plane and type | three cases in `tests/mesh/test_Periodicity.cpp` on `PeriodicBrick`: caps tagged, a typed cap refused, a `Periodic` cap accepted, a side sideset untouched | (c) |
| 3 | DR-186 | both predicates evaluated once per ( shell, block ) pair (a setup-scope `Set< id_t >` per shell), the null-material guard kept, one always-active warning naming block, material, shell, element type and every predicate that failed, then `BELFEM_ASSERT`s; the commented block replaced by a short contract comment (linear gather, μ = μ₀), no ids or dates | syntax gate; sandbox A/B on the release tree: two shipped thin-shell decks (no warning) and one derived deck whose `Conductor` neighbour is `curve : RoxieIron` (warning, once per shell/block); the same negative deck on the debug tree (assert aborts naming the block). The order predicate is unreachable on a fresh run (`:1576` fires first; cohomologies are skipped only on a `.bfm` reload) and is not gated | (c) |
| 4 | DR-187 | one **public static** predicate `ThermalFactory::block_is_supported( DomainType )` = the SELECTION set (`Conductor`, `ThinShell`, `Ferro`, `Buffer`) used by both selection loops and the dispatch; the dispatch's `Air` arm goes with it; `default:` → `BELFEM_ERROR` naming block and type; the dead sideset loop and `tOrder` deleted (O1) with the enriched orders 5 / 13 recorded on the register row | a test in `tests/fem` calling the predicate on every enumerator 0..42 (the tripwire itself cannot fire on today's tree); the sideset gap recorded on the register row | (c), O1 |
| 5 | DR-166 | all three limits — the two header asserts and `init_vertex_container` (`cl_Graph_Vertex.cpp:49-52`) — compare against `std::numeric_limits< decltype( mVertexCounter ) >::max()`, one comment on the release wrap | existing `tests/math/test_GraphVertex.cpp` (70 000 count) unchanged; syntax gate | (c) |

### 3.1 Cross-cutting

- Every refusal added here is setup code (once per run) → `BELFEM_ERROR`; the two DR-186 predicates
  are the exception by ruling (assert plus warning), because the trigger is unknown and a production
  abort on Lawrencium is what the comment was avoiding.
- No `input.conf` key changes; the input contract is untouched.
- Shared checkout: the diff handed to the code jury is own-files only.

## 4. Ordered Steps

- [x] **R1 (DR-184)** Fix the token; write `tests/mesh/test_FacetTables.cpp` (the oracle, one
  `TEST` per element family so a failure names the type; hull tests as in the probe); register it in
  `tests/mesh/CMakeLists.txt`. Gate: the probe binary re-linked against the fixed header shows 0
  issues; the gtest RED/GREEN pair on the scratch build.
- [x] **R2 (DR-182)** The refusal in `tag_periodic_sidesets`; the `PeriodicBrick` fixture gains a
  lateral sideset (id 3, one side face of a layer-0 brick); three tests (caps tagged and the lateral
  sideset untouched; a typed cap refused naming the sideset and its type; a `Periodic` cap accepted),
  each after `set_translation_planes` and `finalize`. Gate: the tests on a scratch build;
  `test_mesh --gtest_filter='Periodicity.*'` green.
- [x] **R3 (DR-166)** The three limit expressions (two in the header, one in `init_vertex_container`) and the comment. Gate: syntax; `test_math` unchanged.
- [x] **R4 (DR-187)** The predicate, the tripwire, the deletion of the dead loop; the enumerator
  test. Gate: syntax; the test on a scratch build; `check-fast` fem suite.
- [x] **R5 (DR-186)** — landed, re-ruled after the code jury to always-active refusals with the material taken by label ( the dof manager's rank-0 blocks are incomplete in parallel ); round 2 at xhigh in flight. The restored guards. Gate: syntax; the three sandbox deck runs, warning
  grep on the logs.
- [x] **R6** Code jury — round 1 done ( 8 findings, 7 confirmed, 1 refuted, P1 on DR-186 ), fixes applied, round 2 narrow in flight; on the five diffs (one thread, slug `review_dr18x_cluster`), verification
  pass, reconciliation.
- [x] **R7** Close-out: register rows (DR-184, DR-182, DR-186, DR-187, DR-166 status cells), devlog,
  index line, Christian runs `make check` in both trees.

## 5. Open Design Questions

- **O1 (DR-187)** **RESOLVED 2026-09-24 → (a), Christian.** Given §1.4, the sideset half of the ruling has nothing to refuse. Options: (a)
  delete the dead loop and record the sideset gap on the row [recommended — it is dead code with a
  design that cannot work]; (b) keep it as a tripwire that refuses any sideset the dof manager ever
  holds; (c) leave it. Christian's call at the plan stop.
- **O2 (DR-186)** Warning deduplication: the loop is per facet, so an unguarded warning prints once
  per facet. ~~Plan: mark the block with its entity flag bit~~ **RESOLVED 2026-09-24 (jury):
  `mesh::Block` carries no entity flag bitset (it is a container, not a `graph::Vertex`); a
  setup-scope `Set< id_t >` of checked block ids per shell, both predicates evaluated and reported
  before the id is inserted.**
- **O3 (by-catch)** **RESOLVED 2026-09-24 → this round, `quad9` / `penta18`, Christian.** `to_string( ElementType )` lacks QUAD9TS and PENTA18TS arms. Two text arms in
  `Mesh_Enums.hpp`; in this round or its own? [recommended: this round, inside R1, it is where the
  test names types]. **Jury 2026-09-24:** the existing TS arms return the parent name (`quad4`,
  `penta6`, `hex8`), so the new arms return `quad9` and `penta18`; the test prints the enumerator's
  integer beside the name.
- **O4 (DR-182)** Allow-list is `Default` and `Periodic` only. The `AirPeriodic` … `FerroPeriodic`
  family is not reachable from a deck and is not on the list; the homology treats it as periodic
  anyway.

## 6. Definition of Done

- [x] five diffs, each own-files, each with its gate run and named in the devlog
- [x] `tests/mesh/test_FacetTables.cpp` red on HEAD's QUAD16, green after
- [x] the DR-187 register row carries the premise correction and the thermal-sideset gap
- [ ] `make check` in both trees — Christian
