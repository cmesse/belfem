# No Automatic Gauge Pin in a φ Component That a Symmetry Plane Already Fixes

**Date:** 2026-10-09
**Purpose:** The automatic gauge pin (one node per connected φ component, since 2026-08-28) is
imposed even when a symmetry plane or a `dirichlet` sideset already fixes φ in that component.
That second condition over-determines the Neumann problem: it forces the magnetomotive force
between the plane and the pin to zero, and the solver answers with a nodal magnetic charge at the
pin. The fix decides, at imposition time and on the full mesh of rank 0, whether the pin's
component already holds a fixed φ, and applies the pin only when it does not.
**Module:** `src/fem/maxwell` (`cl_MaxwellFactory`), `tests/fem`, `doc/`
**AIs involved:** Claude Fable 5.1 (plan, implementation, gates), Codex + Grok (juries)
**Status:** ✅ COMPLETE (2026-10-09). Landed on devel: `impose_autopins()` with the mesh-level contact walk; five-case test (three red on HEAD); `make check-fast` 8/8 release; shipped-deck gate on gantry, costheta, undulator2d (control) and racetrack_usermat plus a two-rank run of gantry and costheta — see `devlog/dl20261009_autopin_symmetry_gauge.md`. Residuals: `gauge` late-flip hazard (pre-existing), no periodic fixture, two-plane MMF note for the guide

> **Scope guards (from the brief relayed by the Server session, Christian's words 2026-10-09):**
> - devel only, never main (pinned at v0.9.2); `src/homology` untouched
> - `find_autopins` selection stays as is: one candidate per component, farthest from any sideset
> - antisymmetry sidesets fix nothing and keep their pin; enclosed buffer patches keep theirs
>   unless they touch a fixed plane; a deck-written `bearing` still disables every automatic pin
> - Christian's ruling (via Server): "the cache caveat is acceptable" — a selection-time fix was
>   allowed. The imposition-time design below is chosen on merit (§2); it has no cache caveat.
> - register: extend a row rather than open one; if a row is needed at all, ask Christian

---

## 1. Current Behaviour and How It Fails

`MaxwellFactory::create_magnetic_kernel()` (`cl_MaxwellFactory.cpp:669` ff.):

| step | lines | what it does |
|---|---|---|
| antisymmetry | `:850-860` | the four AntiSymmetry sideset types are set Inactive; nothing is fixed |
| symmetry | `:874-880` | AirSymmetry / BufferSymmetry / FerroSymmetry: `impose_dirichlet( 0.0 )` on every node of the plane |
| deck Dirichlet | `:915-940` | `dirichlet` / `background dirichlet` BCs: every non-duplicate node of the sideset is fixed to 0 (source-rerouted through `pin_dirichlet_dof` when hanging) |
| autopin | `:944-975` | no `bearing` → every node of `mMesh->autopins()` has dof 0 fixed to 0 (hanging: its node source), unless that very dof is already fixed |
| selection | `:3709-3942` `find_autopins` | one node per connected φ component (air, buffer, ferro blocks; linked through node duplicates `original()` and periodic partners), the node farthest from any sideset (`graph::multibfs`); refuses a component whose nodes all lie on sidesets (`:3917`); persisted as dataset `pinned` in the `.bfm` (`cl_Mesh_BfmFile.cpp:371-379`, reloaded `:423-432`) |

| Failure | Mechanism | Evidence |
|---|---|---|
| Over-constraint of a gauged component | the pin's component already has φ = 0 on a plane; the pin at an interior node x_p adds ∫ h·dl = 0 between plane and pin, a condition on the field | `:874-880` + `:966-972`; mathematics in `tmp/ai_exchange/autopin_symmetry_gauge_claim.md`, literature Brenner & Scott 2008 §5.2, Monk 2003 Thm. 3.37, Bathe 2016 §4.2.2 Eqs. 4.42-4.45, Arsenault et al. 2021 Appendix, Dular et al. 1999 §VI Fig. 3 |
| Observed | Gregory Giard's undulator benchmark with a symmetry plane: "something weird"; a hand-written bearing ON the plane made it go away — because any deck bearing sets `tHaveBearing` (`:898-908`) and the whole autopin loop is skipped (`:944`), so the interior pin is absent. (The plan's first draft credited `:966-972`; corrected after the juries.) | Gregory's report 2026-10-09 |
| Shipped decks exposed | `gantry` (air + ferro symmetry), `costheta` (air symmetry), `racetrack_usermat` (air symmetry `sym1`); none sets a bearing | `examples/*/input.conf` |

**Bottom line:** the automatic pin is a gauge device, and a component whose gauge is already set by
a Dirichlet condition on φ must not receive one; the decision needs the full mesh and the final
set of φ-fixing sidesets, which exist together only at imposition time on rank 0.

## 2. Architecture: Why an Imposition-Time, Mesh-Level Rule

Two places could host the rule.

**Selection time (`find_autopins`)**, skipping gauged components: the symmetry types are known
there (`read_domain_types` precedes it), and the deck Dirichlet BCs are readable from
`mBoundaryConditionFactory`. Rejected for two reasons: (i) the pins are cached in the `.bfm`; the cache tag
(`fn_mesh_config_tag.hpp:330-336`) hashes the whole `topology` section, so a changed symmetry
plane does rebuild the cache, but it deliberately excludes the boundary conditions, so adding or
removing a `dirichlet` BC under an unchanged mesh would silently keep stale pins — not a one-time
"delete the .bfm" caveat but a standing trap (Grok, round 1); (ii) the selection would then depend
on BC semantics that live in the kernel-creation code, two thousand lines away.

**Imposition time (`:944`), mesh level** — chosen: at the autopin loop, rank 0 holds the full
mesh (`Kernel::mesh()`, `cl_FEM_Kernel.hpp:429-432`) and `mMesh->autopins()` is non-empty only on
rank 0 (the Distributor ships no pins; workers never enter the loop with entries). The rule walks
the pin's φ component on the mesh and asks whether it contains a node of a sideset that fixes φ.
The set of such sidesets is collected right where they are imposed (the symmetry switch and the
Dirichlet BC loop), so the rule cannot drift from the impositions.

**Why mesh level and not "is any dof in the component fixed":** the fem sideset groups hold only
the facets owned by the rank (`cl_FEM_DofMgr_SideSetData.cpp:86,115,132`), so at the autopin
site on rank 0 the dof-level fixed state is partial: `DofData::synchronize_dirichlet_bcs()`
(`cl_FEM_DofMgr_DofData.cpp:2619`) runs later, inside `DofManager::init_matrices()`
(`cl_FEM_DofManager.cpp:367`, before `split_dof_container()` at `:370` — D5 corrected the first
draft, which had the order reversed). A dof-level test at `:944` would miss a symmetry plane
owned entirely by other ranks. **The same holds for the sideset ids:** the symmetry switch at
`:862-890` iterates `dofmgr()->sidesets()`, which exist only for sidesets with a locally owned
facet, so the φ-fixing sidesets must be read from the rank-0 mesh by domain type
(`mMesh->sidesets()`, types set by `read_domain_types`, `:563-590`), not collected in that loop
(D1). Mesh sidesets on rank 0 are complete (`mesh::SideSet::nodes()`, collected in
`Mesh::finalize`, `cl_Mesh.cpp:779-782`).

Confidence high on both paragraphs (re-read today).

## 3. Gap Table

| # | State / case | Needed for | Handled today? | Class | Citation / rationale |
|---|---|---|---|---|---|
| 1 | component touching an Air/Buffer/Ferro Symmetry sideset | no pin | no — pinned | (c) | `:874-880` |
| 2 | component touching a `dirichlet` sideset (`background dirichlet` is refused at creation, `cl_MaxwellBoundaryConditionFactory.cpp:227-239`) or a `gauge` node (`cl_FEM_PhysicalBoundaryCondition.cpp:254-259`, fixes φ at impose time, does not set `tHaveBearing`) | no pin | no — pinned unless the pin node itself is on the sideset | (c) | `:915-940`; `dirichlet` reachable (`cl_MaxwellBoundaryConditionFactory.cpp:218-225`); the dirichlet value is time-dependent, the factory pre-pins 0 |
| 3 | component with antisymmetry only | keep the pin | yes | (a) | `:850-860` sets Inactive; nothing fixed |
| 4 | enclosed buffer patch not touching any fixing sideset | keep the pin | yes | (a) | one component per patch, `find_autopins` |
| 5 | pin whose dof hangs (periodic) | existing source reroute | yes | (a) | `:948-965`; the component test runs before the branch |
| 6 | periodic deck: component glued through partners | walk crosses partners AND the partner's duplicates | must mirror `find_autopins` | (c) | `:3805-3848` (partner elements `:3805-3825`, partner's duplicates `:3827-3848`) |
| 7 | interface duplicates (cond\|air, ferro\|air) | walk on originals (`:3732`), expand through duplicates (`:3783-3803`) | must mirror `find_autopins` | (c) | the four neighbor sources: own elements, duplicates' elements, partner's elements, partner's duplicates' elements |
| 13 | a φ-fixing sideset node that is a duplicate | its alias fix is inert (`cl_FEM_SideSet.cpp:337-352`); the Dirichlet loop skips duplicates (`:917`) | seed only from NON-duplicate sideset nodes, mirroring the imposition predicate (D2) | (c) | a plane whose contact with a component is duplicate-only must not suppress that component's pin |
| 8 | reloaded `.bfm` with stored pins | the rule applies unchanged | by design of §2 | (a) | `cl_Mesh_BfmFile.cpp:423-432` |
| 9 | MPI: pins only on rank 0, fix propagated by the Dirichlet sync | unchanged | yes | (a) | `cl_FEM_DofMgr_DofData.cpp:2619-2740` |
| 10 | ConductorSymmetry | not a φ fix | not in the set | (a) | conductors carry no node φ |
| 11 | pin on a node of a fixing sideset itself | already covered | yes | (a) | `:966-972` keeps a fixed dof; the walk also reports the component as gauged, same outcome |
| 12 | every pin skipped | no "unpinned" warning must fire | the warning keys on `autopins().size() == 0`, untouched | (a) | `:982-990` |

### 3.1 Cross-cutting findings

- The walk's adjacency must be the one `find_autopins` uses (elements of Air/Buffer/Ferro blocks,
  node `original()`, duplicates, periodic partners); a narrower walk would under-report gauging on
  periodic or interface-rich meshes and leave the over-constraint in place there. Mirror, do not
  reinvent: same element flag seeding, same three neighbor sources.
- Flag discipline: element flags and node flag planes 0 and 1 are scratch; cleared before and
  after (house rule). Plane 1 marks the nodes of the fixing sidesets (on their originals), plane 0
  marks visited nodes.

## 4. Ordered Steps

- [x] **R1 — the autopin block becomes `MaxwellFactory::impose_autopins()`** (private; the
  existing hanging/regular logic moves in unchanged; the pre-feature warning stays at the call
  site). It returns at once when `mMesh->autopins()` is empty (every worker), asserts rank 0
  otherwise, and prepares the scratch planes ONCE for the whole pin loop: element flag 0 on the
  elements of Air/Buffer/Ferro blocks (as `find_autopins:3717-3738`, but node flags are NOT
  seeded from it — D3); node flag 1 on every non-duplicate node of each φ-fixing contact —
  the mesh sidesets of type Air/Buffer/FerroSymmetry (D1: read from `mMesh->sidesets()` by
  domain type), the sidesets named by `dirichlet` BCs (`tBC->domains()`, the deck list on every
  rank), and the nodes named by `gauge` BCs (`mMesh->node( id )`). Clears all three planes before
  and after the loop (house rule).
- [x] **R2 — `MaxwellFactory::find_fixed_phi_contact( mesh::Node * aPin, Cell< mesh::Node * > &
  aVisited )`** (private): flag-plane-0 BFS from `aPin->original()` with `aVisited` as the queue
  (a `Cell`, reused across pins); neighbors from the four sources of gap 6/7 through flagged
  elements, `original()` on every pushed node; returns the first visited node carrying flag 1,
  else `nullptr`. The CALLER unflags plane 0 on `aVisited` after each search (D3: the per-pin
  search never touches the shared seed on plane 1). The sideset id for the message is found by
  scanning the collected contact sidesets for the returned node (one scan per skipped pin).
- [x] **R3 — apply in the loop**: before the hanging/regular branches, skip a pin whose
  component has a contact, with one `message( InfoLevel::Default )` line per pin on rank 0:
  `automatic gauge pin: node N` or `no automatic gauge pin for the phi component of node N: phi
  is already fixed at node M ( sideset S )`. Skipping is the recoverable outcome; no
  `BELFEM_ERROR` on it (Grok).
- [x] **R4 — test `tests/fem/test_AutoPinGauge.cpp`** (registered under `USE_MAXWELL` in
  `tests/fem/CMakeLists.txt`): a 2-D TRI3 strip, air cells x∈[0,2) and x∈[3,5), copper x∈[2,3)
  across the full height (two φ components), sidesets 24 (x = 0) and 25 (x = 5), written to a
  scratch directory as `.msh` + deck (the `ScratchDir` / `WorkingDir` pattern of
  `test_ThinShellAlignment.cpp:60-139`); `air symmetry` under `topology`, `dirichlet` and
  `bearing` under `boundary conditions`; the edges y = 0 and y = 2 untagged. Each case calls
  `create_magnetic_kernel()`, keeps the returned kernel, and reads the fixed state through the dof
  manager by dof id (`dofmgr()->dof( dofmgr()->calculate_dof_id( node, 0 ) )->is_fixed()`): the
  factory disconnects the dofs from the mesh at `:1009` (D4), so `node->dof( 0 )` is gone. The
  component of a pin is read as x < 2.5 (left) / x > 2.5 (right) and each case asserts exactly
  one pin per side (D4: the pin may sit on the copper interface line x = 2 or x = 3, which are
  air nodes too):
  - no symmetry: both pins fixed;
  - `air symmetry { sidesets : 24 }`: left pin free, right pin fixed, the skip line in the log;
  - `dirichlet { sidesets : 24 }`: as the symmetry case;
  - cache independence (Codex): the no-symmetry deck run first in a directory, then the
    `dirichlet` deck in the SAME directory reusing the `.bfm` (log says reusing): left pin free;
  - `bearing { nodes : 11 }` (node (4,1) of the right strip): no autopin fixed, node 11 fixed.
  Pair with the inverse: the test is built against the unfixed factory first and must go red on
  the symmetry, Dirichlet and cache cases. Not covered by a test, stated as residual: a periodic
  contact reached through partner adjacency (no periodic strip fixture exists; the walk is
  mirrored line by line from `find_autopins` instead).
- [x] **R5 — documentation:** `doc/input_file_reference.md` bearing paragraph (`:951-966`),
  `doc/input_schema.yaml` `bearing_note`, `examples/README.md:154` paragraph, and a `### Fixed`
  entry in `CHANGELOG.md` [Unreleased] naming gantry, costheta and racetrack_usermat as the decks
  whose results change. No `.bfm` deletion is needed (imposition-time rule).
- [x] **R6 — gates:** (a) `make check-fast` in the private release tree, the new test red before
  and green after; (b) shipped decks on the reference binary (devel 53b67ec) and the fixed one,
  serial, same mesh, no `.bfm`: gantry (1 s), costheta (1 s), racetrack_usermat (100 s), and
  undulator2d (20 s, no symmetry — must be unchanged: same pin, identical results). Compared: the
  pinned node ids (`.bfm` dataset `pinned`, identical by construction — the signal is the new
  skip line and the field), the terminal history, and φ / h at the old pin node, in its
  neighborhood and mesh-wide from the Exodus output; (c) the MPI case (Codex P1 / Grok 1): gantry
  and costheta on the fixed binary with `mpirun -np 2` — the skip line must appear and the
  terminal history must match the serial fixed run; the symmetry planes' ownership is whatever
  the partitioner gives, which is the production path.

### 4.0 Implementation Progress (updated 2026-10-09)

Round-1 plan jury (Codex gpt-6-astra/high, Grok grok-4.7/high), all CONFIRMED by re-reading:

- **D1 (HIGH, Codex + Grok)** — R1 collected the symmetry ids from the rank-partial
  `dofmgr()->sidesets()` loop; a plane owned entirely by workers would be missed. Fixed in the
  plan: read `mMesh->sidesets()` by domain type on rank 0.
- **D2 (MEDIUM, Codex; Grok open risk)** — seeding from every sideset node over-reports: a
  duplicate's alias fix is inert. Fixed in the plan: seed non-duplicate nodes only.
- **D3 (HIGH, Codex + Grok)** — the per-pin helper cleared the shared seed plane and, copied from
  `:3717-3738`, would have pre-flagged every φ node as visited. Fixed in the plan: caller-owned
  preparation and cleanup, per-search visited cleanup through the queue, no node seeding from the
  element pass.
- **D4 (HIGH, Grok)** — the test read `node->dof( 0 )` after the factory disconnected the dofs,
  and its x-filters could pass with zero pins. Fixed in the plan: dof manager lookup by id, one
  pin per side asserted, x < 2.5 split.
- **D5 (LOW, Codex + Grok)** — §2 had the sync/split order reversed; `background dirichlet` is
  not deck-reachable; the cache tag covers symmetry (topology) but not Dirichlet BCs. Text
  corrected.
- **Gates added (Codex)** — a two-rank deck run and a `.bfm`-reuse test case. A periodic
  fixture is not added (residual, stated in R4).

## 5. Open Design Questions

- **O1** — Should the rule also skip a pin when the component's φ is fixed by a `bearing` of
  another component? No: a deck bearing disables all autopins already (`:944`). Closed by the
  existing branch.
- **O2** — Message level for the per-pin line: Default (one line per component, at setup) —
  consistent with the DR-135 reroute announcement. Decide in the jury if Verbose is preferred.

## 7. Definition-of-Done Checklist

- [x] R1-R5 landed on branch `autopin_symmetry`, rebased on devel (f6618b9), fast-forwarded into devel
- [x] code jury on the own-files diff; findings reconciled (gauge branch removed, light pass clean)
- [x] R6 gates run and recorded in the devlog with numbers
- [x] Server session told the SHA, the file list and the deck comparison

## 8. Audit Trail

- `tmp/ai_exchange/review_autopin_symmetry_gauge.md` — the mathematics (narrow claim, round 1)
- `tmp/ai_exchange/review_autopin_gauge_fix.md` — plan and code rounds on this file
