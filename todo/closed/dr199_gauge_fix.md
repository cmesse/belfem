# DR-199: The `gauge` Condition Is Fixed at Setup, Guarded Against a Late Flip, and Refused Without a `type`

**Date:** 2026-10-09
**Purpose:** Make the `gauge` boundary condition work as the input reference describes it — a point condition with a time-dependent value — by fixing its dofs at factory time like a bearing, by making its point a contact of the automatic-pin rule, by refusing a first free→fixed flip after initialization inside `Bearing::impose_dirichlet` (the containment both auditors of DR-199 asked for), and by refusing any value-carrying condition that has no `type` at deck read. Both factories, Maxwell and thermal.
**Module:** `src/fem/kernel` (Bearing), `src/fem/maxwell` (factory, BC factory), `src/fem/thermal` (factory, BC factory), `tests/fem`, `doc/`
**AIs involved:** Claude Fable 5.1 (plan, implementation, gates), Codex gpt-6-astra + Grok grok-4.7 (plan and code juries)
**Status:** ✅ COMPLETE (2026-10-09). Landed on devel: gauge points pre-fixed in both factories, gauge point as autopin contact, late-flip refusal in `Bearing::impose_dirichlet`, `type` refusal in both BC factories, gauge arm imposes `mValue`. 19 tests in two suites (red/green), `check-fast` 8/8, strip reproducer serial, costheta with a point serial and two ranks (1.65e-9), pre-fix binary aborts on both. Code jury astra/xhigh + grok xhigh with one fix round (comment on worker bearings corrected: master-only). Residuals: thermal functionally untested; `mValue` regression needs a plugin; Dirichlet arm unscaled (pre-existing); periodic alias contact — see `devlog/dl20261009_dr199_gauge_fix.md`

> **Scope guards:**
> - Rulings (Christian, 2026-10-09, question form): (1) a pre-fixed gauge point is a CONTACT of
>   `impose_autopins` — only its own φ component loses the automatic pin; it does NOT set
>   `tHaveBearing`; (2) a value-carrying condition without `type` is refused at deck read in both
>   BC factories; (3) the thermal factory pre-fixes gauge points too.
> - The bearing's semantics are unchanged (a deck bearing still disables every automatic pin).
> - The double evaluation of the source function per gauge step (`impose_bc:231` and `:258`,
>   Grok's open risk) is NOT touched: secondary, no consumer depends on it.
> - No change to `find_autopins`, to the `.bfm`, to the input keys' names or units. The Maxwell
>   gauge amplitude stays dimensionless, the thermal one K.
> - devel only; never main.

---

## 1. Current Behaviour and How It Fails

DR-199 (`todo/debt_register.md`), verified in c31ea43 and executed on a strip reproducer
(`devlog/dl20261009_gauge_late_flip_investigation.md`):

| Failure | Mechanism | Evidence |
|---|---|---|
| Late flip | a gauge is pre-fixed by no factory (`cl_MaxwellFactory.cpp:898-941` handles Bearing and Dirichlet; `cl_ThermalFactory.cpp:218-236` likewise, then `initialize()`); its first `impose_bc` → `Bearing::impose_dirichlet` → `Dof::fix` (`cl_FEM_PhysicalBoundaryCondition.cpp:254-261`, `cl_FEM_Bearing.cpp:151-197`) lands after the one-time split (`cl_FEM_DofManager.cpp:300-302, 367-372` via `cl_IWG_Timestep.cpp:365`) | debug abort `cl_FEM_DofMgr_SolverData.cpp:1517` (`cl_SpMatrix.hpp:843`); release: raw one-past write there, then Armadillo throw in the gather (`:1952-1961`, `:2143-2153`); rank 1 SIGSEGV on two ranks |
| Untyped crash | `mFun` has no initializer (`cl_SourceFunction.hpp:89`, `cl_SourceFunction.cpp:72-106`); the factories set it only inside `if ( key_exists( "type" ) )` (`cl_MaxwellBoundaryConditionFactory.cpp:429`, `cl_ThermalBoundaryConditionFactory.cpp:215`); `impose_bc:231` evaluates it for every type but bearing and the circuit types | SIGSEGV at the first timestep |
| Over-constraint (never reached) | with the gauge dof pre-fixed, a gauge point in a component with no plane would carry the autopin AND the gauge: the DR-198/12a8eec over-constraint between two points | reasoning, both auditors |

**Bottom line:** fix the gauge dof before the split (then every per-step imposition is a value
update on a fixed dof), make the gauge point a contact so the component gets no second pin,
and put two always-on refusals where the two crashes started.

## 2. Architecture

Four small, independent pieces, each in the layer that owns the invariant:

- **Containment in the kernel** (`Bearing::impose_dirichlet`): the dof manager knows when the
  split is frozen (`DofManagerBase::is_initialized()`, `cl_FEM_DofManagerBase.hpp:224`); a first
  free→fixed flip after it is a wiring error at every call site, so the refusal belongs in the
  bearing, mirroring `pin_dirichlet_dof` (`cl_FEM_SideSet.cpp:80-90`). It fires in release.
- **Pre-fix in the factories**: the bearing loops already exist in both factories; a gauge is
  a bearing with a value, so it joins that loop with `impose_dirichlet( 0.0 )` — the value is
  overwritten by `impose_bc` before the first solve, as it is for a bearing whose value is 0.
- **Contact in `impose_autopins`**: the branch removed from 12a8eec on the code jury's P1s
  becomes valid once the dof is pre-fixed: the gauge point's node (through the mesh vertices,
  `mMesh->vertices()`, `tVertex->node( 0 )->original()`) is flagged on plane 1 like a
  non-duplicate sideset node. It does not set `tHaveBearing` (ruling 1).
- **Refusal at deck read**: both BC factories already own the `if ( key_exists( "type" ) )`
  branch; the `else` becomes a `BELFEM_ERROR` naming the section and the eight accepted
  type words. All 21 shipped decks carry a `type` in every value-carrying block (scanned
  2026-10-09, script in the devlog), so no deck changes behavior.

Rejected: initializing `mFun` to a constant 0 (a silent contract change; ruling 2); making the
gauge set `tHaveBearing` (reintroduces the floating-component hazard; ruling 1); retiring the
key (a documented thermal feature).

## 3. Gap Table

| # | State / case | Needed for | Handled today? | Class | Citation / rationale |
|---|---|---|---|---|---|
| 1 | Maxwell gauge point, plain dof | fixed before the split | no | (c) | R3 |
| 2 | Maxwell gauge point on a periodic alias (single source, unit weight) | source fixed before the split | no; the LIMITATION comment `cl_FEM_Bearing.cpp:144-150` calls for a factory bearing | (c) | R3 reroutes through the same `impose_dirichlet`; the comment is rewritten |
| 3 | gauge point in a component with no plane | no second pin | no | (c) | R4 contact |
| 4 | gauge point in a component with a plane / dirichlet | already no pin; gauge is a second fixed point in that component (two values on one potential, the user's choice as with two planes) | — | (a) | documented, not refused |
| 5 | thermal gauge point | fixed before `mThermalField->initialize()` | no | (c) | R5 |
| 6 | any later `impose_dirichlet` on a free dof after initialization (gauge today, any future caller) | loud refusal in release | no (`fix()` unconditional) | (c) | R1 |
| 7 | worker rank: bearing with linked dofs for a local node (`cl_FEM_DofManager.cpp:1329-1391`) | the pre-fixed flag arrives through `synchronize_dirichlet_bcs()` before the split; the per-step call is a value update | yes once R3 fixes on the master before the sync | (a) | R1 must not fire on a worker whose dof was fixed by the sync — it is fixed, so it does not |
| 8 | worker rank: empty bearing | returns before the refusal | yes | (a) | `cl_FEM_Bearing.cpp:104-121` |
| 9 | untyped `current` / `voltage` / `background` / `dirichlet` / `neumann` / `gauge` block | refusal at read | no | (c) | R2, both factories; `bearing` and the circuit types never evaluate the function |
| 10 | gauge point id that is no point of the mesh | refusal | partly: the bearing refuses on the master at impose time (`:110-113`) | (a)/(c) | R3's factory call hits the same refusal at setup; R4's vertex scan refuses too |
| 11 | warm restart (`cl_FEM_Controller.cpp:5436-5442`, the split frozen even earlier) | gauge fixed before it | yes with R3 (factory precedes restart) | (a) | the restart reads the fixed set from the factory state |
| 12 | `gauge` with `type` and a unit | Maxwell: refused (dimensionless contract); thermal: K | unchanged | (a) | `cl_MaxwellBoundaryConditionFactory.cpp:546`, thermal `:165-172` |

### 3.1 Cross-cutting findings

- The order inside `create_magnetic_kernel()` is load-bearing: bearing/gauge/dirichlet
  pre-fixes (`:898-941`) → `impose_autopins()` (`:944`) → … → the kernel's first `initialize()`
  happens in `shift_fields()` at the first timestep. R3 and R4 both sit inside that block.
- R1's refusal message must distinguish the empty-bearing case (already refused) from the
  late-flip case and name the point and node ids, so a user who hits it through a future
  caller sees which deck line is involved.

## 4. Ordered Steps

- [x] **R1 — containment** (`cl_FEM_Bearing.cpp`): in both branches of `impose_dirichlet` —
  the rerouted source and the direct dof — `BELFEM_ERROR( tTarget->is_fixed() ||
  ! mParent->is_initialized(), "Bearing %lu ( node %lu ): the dof was free when the system was
  initialized and cannot be fixed now — fix it at setup ( a gauge is pre-fixed by the factory )" )`
  before the `fix()`. The LIMITATION comment is replaced by the contract: "a first fix must
  precede initialize(); later calls update the value". `mParent` is a `DofManagerBase *`
  (`cl_FEM_Bearing.hpp:32`) and `is_initialized()` is on the base (`cl_FEM_DofManagerBase.hpp:224`).
- [x] **R2 — `type` refusal** in `cl_MaxwellBoundaryConditionFactory.cpp` (the case group at
  `:410-414`, right after the `expk` refusal) and `cl_ThermalBoundaryConditionFactory.cpp`
  (`:198-215`): `BELFEM_ERROR( tSection->key_exists( "type" ), "boundary condition '%s': 'type'
  is required ( ramp, sigmoid, sine, square, triangle, sawtooth, constant, userdefined )" )`.
- [x] **R3 — Maxwell pre-fix**: `cl_MaxwellFactory.cpp:900-909` gains `else if ( tBC->type() ==
  BoundaryConditionType::Gauge ) { for ( id_t tID : tBC->domains() ) bearing( tID )->impose_dirichlet( 0.0 ); }`
  without touching `tHaveBearing`.
- [x] **R4 — contact**: in `impose_autopins()`, after the dirichlet sidesets, for every Gauge BC
  and every point id: find the vertex in `mMesh->vertices()` (`BELFEM_ERROR` if absent: "gauge
  target invalid: point %lu does not exist on the mesh"), flag `tVertex->node( 0 )->original()`
  on plane 1; the skip message says `( gauge point %lu )` when the contact node is on no contact
  sideset. (after: R3)
- [x] **R5 — thermal pre-fix**: `cl_ThermalFactory.cpp:218-236` gains the same `else if` before
  `mThermalField->initialize()`.
- [x] **R6 — tests** (`tests/fem`):
  - `test_AutoPinGauge.cpp`: (a) deck `gauge { nodes : 11 ; type : constant ; amplitude : 0 ; }`
    → after `create_magnetic_kernel()` node 11's dof is fixed, the right pin is free with the
    `( gauge point 11 )` line, the left pin fixed; (b) an untyped `gauge` deck and an untyped
    `dirichlet` deck are refused at factory construction with a message containing `'type' is
    required` (the `expect_refusal` pattern of `test_DirichletPin.cpp:52-69`); (c) the typed
    gauge deck with the `.bfm` reused (cache case) — same outcome as (a).
  - `test_DirichletPin.cpp` or a new `test_BearingLateFlip.cpp`: the kernel-level fixture
    (`GluedStrip`, `:88-225`) extended with one VERTEX element on node 7
    (`mesh::ElementFactory::create_element( ElementType::VERTEX, 20 )`, node inserted, pushed to
    `mMesh->vertices()` before the Kernel is built — to verify against `cl_Element_Factory.hpp:39`
    and the gmsh reader's `create_vertices`, `cl_Mesh_GmshReader.cpp:1014-1034`): after
    `mField->initialize()`, `mField->bearing( 20 )->impose_dirichlet( 1.0 )` on the free dof is
    refused with "was free when the system was initialized"; the same call on a dof fixed before
    `initialize()` updates the value. Red on HEAD: (a), (b), the refusal.
- [x] **R7 — documentation**: `doc/input_file_reference.md` (`:1058` `type` row: "absent =
  fatal" instead of "treat as required"; the bearing/autopin paragraph `:960-975`: a gauge point
  is pre-fixed and is a contact; the thermal paragraph `:1067-1074`), `doc/input_schema.yaml`
  (`:1689-1698` `effectively_required` → `required: true`, `absent: fatal`; `bearing_note`),
  `CHANGELOG.md` [Unreleased] `### Fixed` (one entry: gauge pre-fixed, contact, two refusals;
  no shipped deck changes), DR-199 row: status cell updated, struck only when the gates ran.
- [x] **R8 — gates**: `make check-fast` in the private release tree (new cases red on HEAD,
  green after); the strip reproducer `typed0` passes assembly on one and two ranks (the zero-load
  deck then hits the step controller as its controls do — "no abort in assembly" is the gate,
  stated as such); the untyped deck refuses at read; costheta 1 s on the fixed binary against
  12a8eec: identical (no gauge, no untyped block → nothing changes). No thermal deck fixture
  exists: R5 is covered by R1's refusal and by reading, stated as residual.

### 4.0 Implementation Progress (updated 2026-10-09)

Round-1 plan jury, every item re-read against the tree before folding:

- **D1 (HIGH, Codex + Grok)** — the R1/R2 sketches lacked their printf arguments. Implemented
  with both ids cast to `long unsigned int` and the section label or type `.c_str()`.
- **D2 (HIGH, Codex)** — the test fixture `gDirichlet` had no `type`; R2 would refuse it. Migrated
  to `type : constant ; amplitude : 0`; a separate untyped fixture feeds the refusal test.
- **D3 (HIGH, Codex + Grok)** — no gate observed a nonzero, time-dependent gauge value after
  initialization through the real dispatch. Added: a constant-3 gauge whose dof reads 3 after
  `initialize()` + `compute_boundary_conditions( 0.5 )`, and a ramp whose values at 0.25 and
  0.75 differ; plus the alias case in the bearing test (source fixed, value updated, late first
  flip refused with "source dof was free").
- **D4 (HIGH, Grok)** — the gauge arm imposed `compute( aTime )` unscaled while `mValue`
  (`compute * mScale`) is what the global publishes; a `userdefined` Maxwell gauge in a declared
  unit was imposed raw. The arm now imposes `mValue`. Scope widened by one line; stated here.
- **D5 (HIGH, Grok)** — the skip message needs a side record: gauge contact nodes and point ids
  are kept in two parallel `Cell`s and the message says `( gauge point N )`; the node is taken
  from `dofmgr()->bearing( id )->node()`, keyed by vertex id, not from a node-id scan.
- **D6 (HIGH, Grok)** — the stale "worker bearings never hold dofs" sentence is replaced with
  the contract (bearings hold dofs on every rank whose local node carries them; the first fix
  precedes `initialize()`; later calls update the value).
- **D7 (HIGH, Grok)** — R1 guards both branches (source and direct); the `pin_dirichlet_dof`
  mirror covers the hanging branch only. Implemented on both.
- **D8 (MEDIUM, Codex + Grok)** — thermal has no deck-driven test fixture
  (`test_ThermalBlockAdmission.cpp` tests a static predicate); R5 is **statically reviewed,
  functionally untested** (protocol §11 wording). The containment R1 is its safety net.
- **D9 (LOW, Grok)** — R7's line numbers corrected: `type` row `:1068`, thermal paragraph
  `:1077-1085`, bearing prose `:973-998`, schema `bearing_note :1598-1620`; `required: true`
  is the schema's convention. The 21-deck `type` scan was executed (section scan, no hit).
- Grok's open risks taken: vertices pushed before `finalize()`; the untyped refusal test wraps
  the constructor; the thermal refusal sits inside the Gauge/Neumann/Dirichlet case, never on
  the bearing case.

Code round (Codex astra/xhigh, Grok 4.7/xhigh), one fix round:

- **D10 (HIGH, Codex + Grok)** — my replacement comment said bearings hold dofs on every rank;
  `connect_dofs_to_mesh` runs on the master only, so worker bearings are empty. Corrected; the
  plan round's adjudication is reversed and recorded.
- **D11 (LOW, Grok)** — refusal text overclaimed; reworded. Ramp test asserts 0.5 / 1.5;
  `<cmath>` included; dead point check → assert; schema contact sentence names the gauge point.
- **D12 (MEDIUM, Codex)** — no regression distinguishes `mValue` (needs a plugin). Residual.
- **D13 (MEDIUM, Grok)** — Dirichlet arm still unscaled (pre-existing); periodic alias contact
  untested. Residuals.

## 5. Open Design Questions

- O1 RESOLVED 2026-10-09 → Christian: gauge point is a contact, not a `tHaveBearing`.
- O2 RESOLVED 2026-10-09 → Christian: absent `type` is refused at deck read.
- O3 RESOLVED 2026-10-09 → Christian: thermal factory included.
- O4 — a gauge point whose node is a cut duplicate: `Bearing::impose_dirichlet` refuses a
  multi-source hanging dof (`:151-156`); a gauge there is refused at setup by R3, which is the
  bearing's existing behavior. Left as is; noted in the docs.

## 7. Definition-of-Done Checklist

- [x] R1–R7 landed on a private worktree branch from the devel tip, rebased onto cfaf7e8, code-juried
- [x] R8 gates run and recorded in the devlog with numbers
- [x] DR-199 row updated (status; struck only on Christian's word)
- [x] Server session told the SHA

## 8. Audit Trail

- `tmp/ai_exchange/review_gauge_late_flip.md` — the diagnosis round (DR-199)
- `tmp/ai_exchange/review_dr199_gauge_fix.md` — plan and code rounds on this file
