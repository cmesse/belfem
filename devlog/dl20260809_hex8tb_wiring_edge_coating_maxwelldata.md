# Devlog 2026-08-09 — HEX8TB FEM Wiring, Edge-Coating Input, MaxwellData Migration

**Date:** 2026-08-08 / 2026-08-09 (one continued session)
**Topic:** Side-connector wall elements wired end-to-end: fem::Element generation
and MPI distribution, opt-in `edge coating` input keys, kernel migration onto
MaxwellData dispatch, legacy purge, first smoke runs
**AIs involved:** Claude (implementation), Codex + Grok (three jury rounds),
Christian (design rulings, builds, runs)
**Claude Confidence:** high on compile/link and startup wiring; medium on
end-to-end behavior until the smoke runs finish, and on physics until the
R7 battery runs
**Codex/Grok Audit Confidence:** high (all three rounds converged; every
cited reference rechecked)
**Literature:** Messe et al. 2023 (paper1) §4 cited for the checkerboard
tolerance context; no new formulation content this session
**Verification:** compile/link green (Christian's builds); serial + parallel
smoke runs launched — parallel running at session close after the width fix.
"Reviewed" applies to the physics; "verified" only to compile and startup.

## Summary

The HEX8TB edge-coating wall element went from "mesh-level construction only"
to a fully assembled path in one continued session: fem::Elements are
generated, linked to their recovery facets, activated with the conductor dof
table, dispatched to a MaxwellData-backed kernel pair, and fed by a per-tape
opt-in input key, with first MPI distribution fixes in place (P9a/P9d remain
open). Three frozen-protocol
jury rounds ran (edge-coating policy, wiring batch, MaxwellData draft); the
third caught four P0 defects in the draft *before* implementation. A ~1,470
line legacy purge (`SideLayerOld`) and two smoke-run defects (const
accessors, off-root block thickness) closed out the session with a parallel
run in flight.

## Key changes

### R4 geometry interpolation settled; D8 fusing gates demoted (2026-08-08)

HEX8 Lagrange shape functions are reused for the wall via case labels
(`cl_IF_InterpolationFunctionFactory.cpp` and the `nedelec_data_linear_h`
group in `cl_FEM_Calculator.cpp`); spatial derivatives stay on the exact-cuboid
metric. The deliberate variational crime (interpolation on true node coords,
analytic metric for derivatives) is recorded with its rationale and two
cautions in `side_connector_wall_element.md` §4.5. D8 closed; its fusing
gates were demoted to `mFuseEdges` re-enable preconditions in
`todo/side_edge_fusing_cut_aware_plan.md`.

### P1–P3 — facet wiring, activation, dispatch (2026-08-08)

`DofMgr_BlockData` links every wall element to its recovery facet (id = wall
id + 1, lives in the hidden sideset's facet map, NOT in `elements()`) and to
its master layer-block element, in both the serial and the parallel variant
(rank-0 `share()` of connector element ids; guard ladder on consumers).
Activation rides `set_block_types_in_magnetic_equation` + FieldList
(connectors share the Conductor edge-h dof table, Christian's edit);
`link_to_group` dispatches `Left/RightCoating`. The jury's "connector
blocks never selected into the kernel" blocking claim was REFUTED twice by
code trace: the Maxwell path selects ALL mesh blocks via
`IWG_Maxwell::initialize`; the metis-feeding `tBlocks` switch is partitioning
only. Facet ownership inheritance (wall ← facet ← master) was confirmed by
source trace as already implemented by the D6 pass in
`Kernel::partition_mesh`.

### Edge-coating policy + P4/P5 — input keys (jury round 1, wired 2026-08-09)

Christian's policy, refined by the jury: `edge coating : on` is OPT-IN per
tape; the guard triad (≥ 3 layers, identical top/bottom material — label
equality covers RRR — identical thickness with relative tolerance) lives as
loud `BELFEM_ERROR`s in the live `create_side_connectors`; wall material and
width DERIVE from the outer stabilizer layer (no separate material key;
`mConnectorWidth` 0.0 = derive sentinel, `edge coating width` overrides).
Parsing sits in `read_thin_shell_data` (per-tape on the `Protoshell`, fatal
on 2-D, width-without-flag, non-positive width);
`ThinShellFactory::create()` adopts the flags per tape.
`doc/input_file_reference.md` §8 synced same turn; the Codex prose pass
corrected the bool-token list (only `true`/`on`/`yes`/`1` are true — anything
else parses silently as false) and flagged the then-standing WIP stop.
`assign_materials` gates the coating on `MaterialType::PureMetal` (YBCO
subclasses `Metal` but is typed `HTS`) and warns when it is not copper.

### Frequency rule codified (2026-08-08)

Christian's policy from the wiring-batch fixes: the ASSERT-vs-ERROR tier
follows the call rate — hot paths keep `BELFEM_ASSERT`, once-per-run
setup/init code is generous with `BELFEM_ERROR`, no per-case deliberation.
Codified in `doc/coding_philosophy.md`; applied to the jury's F1–F3
hardening fixes (owned-wall-missing-facet, null facet master, zero
thickness).

### SideLayerOld purge (~1,470 lines, 2026-08-09)

The dead 13-argument `create_side_connectors` overload, its nine helpers,
the uncalled `compute_side_indices`, the `SideLayerOld`/`SideFacet` structs,
eight dead-only members, the commented old-signature call blocks, and five
orphaned includes were removed from `ThinShellFactory`. Kept deliberately:
the ctor's Air/Ferro element flag-1 loop (global mesh state that
CutFactory/EdgeCutter may read) and `mConnectorsForAllLayers` (guarded in
the live builder). **Incident, disclosed:** the first deletion script
assumed `//----` separators exist between all functions — this file has only
16 at column 0 — and one span swallowed ~1,500 lines including live
functions. Because the swallowed region existed in HEAD and appeared as pure
`-` lines in `git diff`, and none of the session's uncommitted edits were
inside it, the region was restored losslessly from the diff and the
removal redone with asserted span boundaries (def counts + brace balance per
span). Lesson recorded: separator-anchored span scripts are unreliable here;
anchor on definition lines and assert content.

### h_side_connector → MaxwellData dispatch (jury round 3, implemented)

Christian's diagnosis: the kernel could not reuse `mx->compute*` (T, B, β
computed differently than PENTA6TS) and the ρ tangents were unwired; his
proposal — connector-specific function pointers — was drafted, jury-audited,
then implemented. Core insight (C1, confirmed 3/3): the metal ρ family
consumes only cached `compute_b`/`compute_j`/`compute_T` + `bj_angle`, so
only the field-assembly layer needed connector variants; the tangent layer
is reused verbatim and the Newton part became the existing
`add_rho_field_tangent`. The jury caught four P0s pre-implementation, all
fixed: `mFundMudH` left null (crash in the Newton path); thermal runs
adopting the EmptyBlock calculator as peer (ThermalFactory never selects
connector blocks; `allocate()` skipped MaxwellData while the `link_maxwell`
rebuild created it with a bogus peer — fixed by a `block_exists` check in
the ctor, a connector exception in `allocate()`, and non-thermal link
dispatch for connectors); the master-`hn` copy left implicit (wrong physics
if omitted); `EF_PENTA6TS` friendship breaking on the helper move (the free
function moved to `cl_FEM_Calculator.cpp` under its own name instead).
Christian's ruling: constant-μ gate (`BELFEM_ERROR`) instead of porting the
dMdx blocks — plated walls are copper. Workspaces relabeled
`normal`/`tangent`/`binomial` (+`Tseam`), removing both documented aliasing
hacks; kernels are now a Picard/Newton pair (`h_side_connector`,
`h_side_connector_newton`) mirroring the conductor dispatch, with
`save_resistivity` and the |B|/β tangent channels (dρ/dj ≡ 0 for metals —
the old broken `drho` block, which also did not compile, is gone). Per-tape
temperature: seam-node interpolation with `gTbulk` fallback, under the
`compute_T_fem` clamp contract.

### Smoke-run fixes (2026-08-09)

1. `Cell<T>` gained const `first()`/`last()` overloads — the material
   derivation and the guard triad call them on const refs (build error).
2. First parallel start aborted in `EF_HEX8TB::link` (width NaN): the
   distributor's `send/receive_block_data` shipped no thickness, so off-root
   blocks lost it (P9b). Proof the rank-0 mesh was healthy: the `.bfm` saved
   by the same run carries `widths = 2.5e-6` read from the block objects.
   Fix: `proto::GroupData::mThickness`, broadcast appended to the block-data
   exchange, applied at ProtoMesh block materialization — walls AND
   thin-shell layers now carry thickness off-root through one channel.
3. The `.bfm` WIP stop (P6) removed once the path compiled; the sidecoatings
   deck gained `edge coating : on` (triad satisfied: copper 2.5 µm top and
   bottom; derived wall width 2.5 µm).

## Open items

- **R6 tail:** thin-shell type lists in `cl_FEM_Postprocessor.cpp`
  (~:246-250, :297-301) do not know HEX8TB — output-side only.
- **R7 battery:** numeric handedness (O1), wall current vs analytic r′,
  R0.3 probe, twisted-helix + cut-station regressions; r′ calibration dial
  (effective ρ vs interface-resistance term) decided then.
- **P7b:** the seam T-copy writes the mesh "T" field from assembly (MPI
  hazard) — move to a postproc/setup pass or gate on serial.
- **P8:** `side_connector_wall_element.md` §5/§6 still need sync to the
  implemented master-block/normal-calculator recovery route.
- **P9a:** `send/receive_thinshell_data` still ships no connector record —
  off-root `ThinShell::side_connector_blocks()` is empty (the dof-side facet
  pass does not need it; anything else iterating it off-root does).
- **P9d:** aura shipping of the recovery-facet master to neighbor ranks.
- Parallel smoke run in flight at session close.

## Files updated

- `src/fem/kernel/cl_ThinShellFactory.{cpp,hpp}` — P5 adoption, guard triad
  (prior day), SideLayerOld purge, WIP stop removal
- `src/fem/kernel/cl_FEM_DofMgr_BlockData.cpp` — facet wiring (both variants)
- `src/fem/kernel/cl_FEM_Calculator.{cpp,hpp}` — MaxwellData connector
  branch, field variants, frame prep, moved friend helper, thermal-peer +
  dispatcher + allocate fixes, case labels
- `src/fem/maxwell/cl_MaxwellFactory.cpp` — input keys, material derivation
  (create + reload), `assign_materials` PureMetal gate
- `src/fem/maxwell/cl_IWG_Maxwell.cpp` — dispatch pair, workspace relabel
- `src/fem/maxwell/matrices/mt_maxwell_h.{cpp,hpp}` — kernel pair, helper
  removal
- `src/fem/interpolation/cl_IF_InterpolationFunctionFactory.cpp` — HEX8TB case
- `src/mesh/cl_Protoshell.hpp` — edge-coating flag + width
- `src/mesh/st_ProtoMesh.hpp`, `src/mesh/cl_ProtoMesh.cpp`,
  `src/mesh/cl_Mesh_Distributor.cpp` — block thickness off-root
- `src/containers/cl_Cell.hpp` — const `first()`/`last()`
- `doc/input_file_reference.md`, `doc/coding_philosophy.md` (frequency rule),
  `src/fem/maxwell/doc/side_connector_wall_element.md` §4.5
- `todo/hex8tb_phase2_fem_wiring.md` (master tracker),
  `todo/side_edge_fusing_cut_aware_plan.md`, `todo/README.md`,
  `todo/closed/` (two supersessions)
- `cmake-build-debug/sidecoatings/input.conf` (deck key)

## Exchange threads distilled (GC-eligible)

`tmp/ai_exchange/review_edge_coating_policy.md`,
`review_hex8tb_wiring_batch.md`, `review_sideconnector_maxwelldata.md`,
`edge_coating_policy_spec.md`, `sideconnector_maxwelldata_draft.md`,
`hex8tb_wiring_batch.patch`.
