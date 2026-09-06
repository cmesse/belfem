# Devlog 2026-04-07 — Layered patch test (copper / silver) trace

**Date:** 2026-04-07
**Topic:** Static read-only trace of the DOF / interface wiring used by a copper-silver thin-shell tape, prompted by a ~20 dB error in the assembled solution that the user suspects is a *patch-test failure*.
**AIs involved:** Claude (this trace)
**Confidence:** see per-finding tags
**Status:** Trace in progress; this devlog is the snapshot left for the user when they return tomorrow.

## Test setup

`cmake-build-debug/input.conf`, with the material stack reduced to:

```
layers : tape
{
    copper    : 2.5 mum ;
    silver    : 1.5 mum ;
}
```

Two `thinshell : tape` blocks (sidesets 5-10, 11-16). 150 A current ramp BC through curves 1-2 → 5-6 and 3-4 → 7-8. Cohomology cuts (`generalized pellikka`) and current/voltage BCs are still active. **Important:** this is *not* a pure patch test with uniform `H_∞`; it's the full simulation reduced to two normal metals at 77 K so that all of the HTS / Hastelloy contrast story is removed.

What this rules out compared to previous runs:
- No power-law nonlinearity (Picard / Newton iteration is trivially well-conditioned).
- No harmonic-mean penalty collapse — both ρ are similar (Cu RRR=30, Ag RRR=100, both ~10⁻⁸ Ω·m at 77 K), so `eta * (km + ks)` gives a healthy α ≈ 0.07-0.1 with current code.
- Current Nitsche penalty (`mt_maxwell_h.cpp:2086-2088`) is `eta * (rho_m/hm + rho_s/hs)`, with the SIPG corrections (verified yesterday) on lines 2149-2152.

So the residual error is structural, not a tuning problem. Quoting the user: "the error is still in the 20 dB region… I think this means we fail the patch test."

## What's been verified to be correct (high confidence)

### 1. The hang_thinshell_edges_*_bottom / *_top functions are structurally consistent

Read `cl_MaxwellFactory.cpp:1364-1680`. There are four functions:

| Function | Layer | Vol side | Reorder? |
|---|---|---|---|
| `_on_nodes_bottom` (1364) | first layer (PENTA face 0) | master | no — master direct |
| `_on_nodes_top` (1424)    | last layer (PENTA face 1)  | slave  | yes via `to_master_orientation` |
| `_on_edges_bottom` (1490) | first layer (PENTA face 0) | master | no — master direct |
| `_on_edges_top` (1582)    | last layer (PENTA face 1)  | slave  | yes via `to_master_orientation` for both nodes and edges |

The asymmetry between bottom and top (master direct vs slave + reorder) is correct because by `update_facet_nodes` (`cl_Mesh.cpp:683-702`) the *facet element's* nodes are populated from `master->get_nodes_of_facet(index_on_master())`, i.e., always master local order. So the bottom variant talks to the master in master order natively, while the top variant has to reorder slave nodes/edges into master order first.

For copper/silver (both ends air), only `_on_nodes_bottom` and `_on_nodes_top` are called (the `_on_edges_*` paths are exercised when the master/slave is a conductor). The node-source path is what is actively running for this test.

### 2. `to_master_orientation` for TRI3 orientation 1 is correct

Read `fn_to_master_orientation.cpp:45-80`. For TRI3 orientation 1 (the "no rotation, just look-at-from-the-other-side" case), the mapping is:

```
master[0] = slave[0]
master[1] = slave[2]
master[2] = slave[1]
```

This swaps local indices 1 and 2, which is exactly the "look at the triangle from the other side" reflection. Combined with the master TET4's outward CCW vs the slave TET4's outward CCW, the geometry and the mapping are consistent.

### 3. ThinShellFactory layer-interface edge wiring is consistent

Read `cl_ThinShellFactory.cpp:1463-1623`. For the order==1 path in `link_elements_with_edges`:

```cpp
for ( Block * tBlock : aBlocks )
{
    tBottom = aLayers(l)->hasDuplicates ? EdgeDuplicates : Edges ;
    tTop    = aLayers(++l)->Edges ;
    ...
}
```

For `[copper, silver]` with `aMaterials(0) != aMaterials(1)`, `aLayers(1)->hasDuplicates = true`. Therefore:

| Block | Bottom edges | Top edges |
|---|---|---|
| 0 (copper) | `aLayers(0)->Edges` | `aLayers(1)->Edges` |
| 1 (silver) | `aLayers(1)->EdgeDuplicates` | `aLayers(2)->Edges` |

So copper's *top* face references `aLayers(1)->Edges` and silver's *bottom* face references `aLayers(1)->EdgeDuplicates`. **Different `Edge*` objects, different DOFs.** The two free edge DOF sets at the layer interface are properly decoupled. ✓

### 4. The ghost facet creation `set_master/set_slave` is consistent

In `create_ghost_facets` (`cl_ThinShellFactory.cpp:1656-`):

```cpp
tGhost->set_master( tMasterBlock->element( f ), 1 );  // copper, face 1 (top)
tGhost->set_slave ( tSlaveBlock->element( f ), 0, 1 ); // silver, face 0 (bottom), orient 1
```

PENTA6TS face 1 is `[mNodes[3], mNodes[4], mNodes[5]]` (top triangle) and face 0 is `[mNodes[0], mNodes[1], mNodes[2]]` (bottom triangle). Master copper top + slave silver bottom is the correct geometric pairing for the layer interface.

The hard-coded orientation `1` matches the integration table convention checked against `slave_integration_penta` (the `- 1` makes it 1-indexed → table index 0 → identity mapping). ✓

### 5. Edge canonical orientation is consistent across `Edges` and `EdgeDuplicates`

`create_edges_on_layers` (`cl_ThinShellFactory.cpp:1463-1504`) constructs both `Edges` and `EdgeDuplicates` from the same `tOrg` template, with the same `tDup->insert_node(...)` order. So copper's `aLayers(1)->Edges[i]` and silver's `aLayers(1)->EdgeDuplicates[i]` have:

- the same canonical direction (lower-original-id-first, inherited from temp edge construction);
- different `Node*` containers because the underlying layer-1 nodes are reused (yesterday's "shared interface node" issue), but…
- …same node order, so the edge directions agree.

Therefore `mS_copper[3] == mS_silver[0]` (etc.) under `Mesh::compute_edge_directions`, and the EF_PENTA6TS basis at `mE(:, 3)` for copper (zeta=+1) gives the *same vector value* as `mE(:, 0)` for silver (zeta=-1) at the same `(L1, L2)` integration point. The Nedelec basis at the layer interface is geometrically consistent across master and slave. ✓

### 6. SIPG `h_ghost` is symmetric and the Dm/Ds construction does the right thing for uniform H

The current `mt_maxwell_h.cpp:2147-2152` reads:

```cpp
Kmm += trans(Em)*Em*awdS - trans(Em)*Dm*rwdS - trans(Dm)*Em*rwdS ;
Kms -= trans(Em)*Es*awdS + trans(Em)*Ds*rwdS - trans(Dm)*Es*rwdS ;
Ksm -= trans(Es)*Em*awdS - trans(Es)*Dm*rwdS + trans(Ds)*Em*rwdS ;
Kss += trans(Es)*Es*awdS + trans(Es)*Ds*rwdS + trans(Ds)*Es*rwdS ;
```

Symbolically `trans(Kms) == Ksm`. Symmetric. ✓

For a uniform H represented as `q = [q0, q0, q0, q0, q0, q0]`, hand-evaluating `Dm·q` from the construction at lines 2128-2142:
- columns 0-2 of `Dm` are `−Em(:, 3..5) / hm`
- columns 3-5 of `Dm` are `+Em(:, 3..5) / hm`

`Dm·q = q0 · ( −Em(:,3) − Em(:,4) − Em(:,5) + Em(:,3) + Em(:,4) + Em(:,5) ) / hm = 0`.

So the consistency term vanishes on uniform H. ✓ Same for `Ds`. The penalty term also vanishes when `q_master = q_slave`. ✓

### 7. The PENTA6TS bulk K matrix in `h_ts_metal` is the standard `∫ ρ C^T C dV` form

Read `mt_maxwell_h.cpp:90-176`. This is the same code path that has been used for many previous test cases. Nothing has been changed here in the recent commits. The mass and stiffness matrices use `aCalc->E(k)`, `aCalc->C(k)`, and `aCalc->dV(k) = mThickness * mSurface` from `EF_PENTA6TS::link`. For uniform H with C·q ≈ 0, the contribution to the residual is 0. ✓ on paper.

### 8. T-matrix construction for `edge ↔ node-pair` (LINE2 case in `DofData::create_dofwise_t_matrices_master`)

Read `cl_FEM_DofMgr_DofData.cpp:3563-3760`. For a shell edge with two source NODES (the air phi DOFs), the coefficients are `[1, -1]` and `tWeights = tNodeWeights * tCoefficients`. The result is:

```
edge_DOF = phi(node 0 of shell edge) - phi(node 1 of shell edge)
```

which is `∫_(0→1) (-∇φ) · dl = ∫_(0→1) H · dl` for the H-formulation in air. Sign and magnitude both correct. The shell edge's `node(0)` and `node(1)` correspond to air volume nodes at the same `xy` (verified by walking the temp-index pipeline). ✓

## What I haven't verified yet (open items)

### A. The chain in `DofData::create_dofwise_t_matrices_master` for the *full* shell-edge population

Beyond the simple "edge has 2 node sources" case, I haven't traced the case where a shell edge has source nodes that are themselves *hanging* (e.g., when the air-side node was already duplicated by `CutFactory`). Lines 3656-3717 handle this recursively, and there's potential for an off-by-one or wrong-side index.

### B. How `h_ts_metal` builds `bn` from the air-side phi DOFs

Lines 117-125 of `mt_maxwell_h.cpp`:

```cpp
Calculator * tCalc = aCalc->get_normal_calculator( phi_m, phi_s );
bn = tCalc->Bm(0) * phi_m + tCalc->Bs(0) * phi_s ;
bn *= -0.5 * constant::mu0 ;
```

`tCalc->Bm(0)` is `mInvJm * dNdXi(0)` (verified at `cl_FEM_Calculator.hpp:1244`), i.e., the **full 3D gradient** of the master volume basis at integration point 0. So `bn` is the full averaged `B` vector, *not* the normal projection only. Then `b = bt + bn` adds the in-plane field from the shell edges to the *full* averaged volume B, which seems to *double-count* the in-plane part of `bn`.

For copper/silver this is irrelevant because `mat->rho()` for a normal metal ignores `b`, so the K matrix is unaffected. **But it's a smell** I want to come back to.

### C. The interaction of cohomology cuts with the thin shell duplicate node setup

`CutFactory::link_node_duplicates_and_originals` (`cl_CutFactory.cpp:2620-2717`) builds chains via `tOrg->original()->add_duplicate(tDup)` and `tDup->set_original(tOrg->original())`. If `tOrg` already had a `set_original` from `CutFactory::duplicate_nodes_on_face_sidesets`, the chain has length 2 and the recursion in `create_dofwise_t_matrices_master` may not see the right "original".

### D. The volume air block boundary condition application

The current 150 A BC is applied via the cut surfaces. The cut DOF substitution and the air block phi DOFs interact at the interface to the tape. I haven't verified that the BC is consistent with the shell edge's `[1, -1]` T-matrix coefficient convention.

### E. `slave_integration_penta`'s `sNodes[3][6]` table for orientation 1

Hard-coded as identity (`{0,1,2,3,4,5}` first column for slave facet 0 orientation 0). The `-1` in `slave_integration_penta` makes orientation index 1 → table column 0 → identity. **For our specific test all ghost facets use this identity branch**, so the orientation table can't be wrong here. But for any case where `Facet::compute_orientation` returns 2 or 3, the `_on_nodes_top`-style reorder must agree with `slave_integration_penta`. Worth a separate check.

## My current best hypotheses for the 20 dB error

I have not pinpointed the bug. These are *ranked candidates*, not findings:

1. **`bn` double-count in `h_ts_metal`** (item B above). For copper/silver this only affects the post-processed B, not the K matrix. So this *cannot explain* the 20 dB error in the assembled solution. **Demoted.**

2. **A wrong sign in the T-matrix `[1, -1]` coefficient when the shell edge's canonical orientation disagrees with the air node ordering it's hung on.** The convention is "shell edge node 0 → temp index 0 → master volume node 0", but the shell edge's *canonical* direction is set by `original()->id()`, while the air volume face's local order is set by `master->get_nodes_of_facet`. These can disagree if the master volume's first face node is *not* the lower-id one. **High suspicion**, would explain a sign error per-edge that compounds globally. Worth a focused trace tomorrow.

3. **Missing `mS` correction at the T-matrix level.** The basis function evaluation `EF_PENTA6TS::E()` applies `mS[i]` to flip the basis sign when the local element edge order disagrees with the canonical edge direction. But `DofData::create_dofwise_t_matrices_master` writes the T-matrix in terms of the *canonical* edge DOF without re-applying `mS`. If the local element's edge 0 has `mS = -1`, the K matrix entries get flipped, but the T-matrix doesn't know. **Plausible factor-of-2 → factor-of-many** error. **High suspicion.**

4. **Wrong scaling between layer thickness and mesh units.** Mesh is in `mm`, layer thickness in `mum`. SI conversion in `MaxwellFactory` looks correct. Probably not it.

5. **Cohomology cut sideset interaction with thin-shell nodes.** Item C above. Possible but harder to trace without a test run.

6. **`get_normal_calculator` returns wrong volume neighbor when both master and slave are air with duplicate nodes from `CutFactory::duplicate_nodes_on_face_sidesets`.** The `tPhi(tMaster->node(k)->index())` indexing is on the global mesh, not on the duplicated set, so the slave's "phi" might be read from the master's nodes if the master/slave node arrays still alias. Need to verify after `CutFactory` runs.

## Concrete next-day actions (in order of cost / value)

### Step 1 — Cheap, no rebuild

Look at `DofData::create_dofwise_t_matrices_master` `EntityType::EDGE` branch (`cl_FEM_DofMgr_DofData.cpp:3646-3760`) and verify by hand that **for a single shell edge with two non-hanging source nodes**, the resulting `tWeights` is correct *under the actual canonical direction of the shell edge* — i.e., walk through what happens when the lower-id air node is at facet position 1 instead of facet position 0. If the temp-index lookup `tNodesOnVolume(tEdge->node(k)->index())` ends up handing the wrong air node to the wrong temp slot, that's hypothesis #2 confirmed.

### Step 2 — Cheap, no rebuild

Walk `EF_PENTA6TS::E()` on a concrete example:
- A flat triangle with vertices at `(0,0,0), (1e-3,0,0), (0,1e-3,0)` (1 mm scale).
- Layer thickness 2.5 μm.
- Apply uniform `H = (1, 0, 0) A/m` and write down what the 6 edge DOF values should be (in canonical direction).
- Plug those 6 DOF values into `E(k) * q` for `k=0..ng` and check what field comes out. If it's not `(1, 0, 0)`, the basis evaluation is the bug.

### Step 3 — Low cost, single build needed

Add a debug-only assertion in `DofManager::compute_jacobian_and_rhs` (or wherever the assembled K is finalized) that checks **`‖K · q_uniform‖_∞ < 1e-10 · ‖K‖_∞`** for a programmatically-constructed `q_uniform` corresponding to a uniform `H_∞` patch. If it fails, we know exactly which row is broken (just print the worst row and its DOF id).

### Step 4 — Bigger build, but very informative

Strip the `homology` block from `input.conf`, drop the `current` BCs, and replace with a uniform tangential `H_∞` boundary BC at the top and bottom of the air region. **No cuts, no current, no time-stepping** — just a static uniform-H drive. If this still gives 20 dB error, the cohomology / current path is innocent. If it gives 0 dB error, the bug is in the cut interaction.

### Step 5 — When all else fails

Compile with `-DDEBUG_DUMP_K`, run the smallest possible 1-element-tall stack `[air-tet | copper-tri | silver-tri | air-tet]` and dump K to a CSV. Compare master rows vs slave rows for the layer-interface DOFs by hand. If the dump shows a non-symmetric or non-zero entry for uniform H, the bug is identifiable from the dump alone.

## Files touched (read-only)

- `src/fem/maxwell/cl_MaxwellFactory.cpp` — `hang_thinshell_edges_on_*_bottom/top` and the call site at lines 1220-1303
- `src/mesh/cl_ThinShellFactory.cpp` — `create_temporary_edges`, `create_edges_on_layers`, `link_elements_with_edges`, `create_ghost_facets`, `create_nodes_on_layers`
- `src/mesh/cl_Element_PENTA6TS.hpp` — face/edge node tables
- `src/mesh/fn_to_master_orientation.cpp` — TRI3 orientation table
- `src/mesh/cl_Mesh.cpp` — `compute_edge_directions`, `update_facet_nodes`
- `src/mesh/cl_EdgeFactory.cpp` — `grab_nodes`, `edge_key` (canonical direction by *index*, not id)
- `src/fem/interpolation/nedelec/cl_EF_PENTA6TS.cpp` — basis function evaluation
- `src/fem/interpolation/fn_IF_initialize_integration_points_on_facet.cpp` — `intpoints_penta_ts` master and slave overloads
- `src/fem/maxwell/matrices/mt_maxwell_h.cpp` — `h_ghost`, `h_ts_metal`
- `src/fem/kernel/cl_FEM_Calculator.{hpp,cpp}` — `Bm`, `Bs`, `Bscalar_master/slave`, `get_normal_calculator`
- `src/fem/kernel/cl_FEM_DofMgr_DofData.cpp` — `create_dofwise_t_matrices_master`

## Source edits made

**None.** This trace is read-only.

## Open question for the user

**What metric is "20 dB error"?**

- If it's `20·log10(‖H_computed − H_analytic‖ / ‖H_analytic‖)` then 20 dB means the error is **10× the signal**, which is catastrophic.
- If it's `20·log10(signal/noise)` or similar SNR, then 20 dB means a **10% error**, which is bad but not catastrophic.
- If it's a residual norm comparison, the interpretation is again different.

This affects how I prioritize hypotheses tomorrow. Knowing the metric tells me whether to chase a sign flip (gives ~200% error) or a unit-conversion factor (gives 1000× error) or a basis function magnitude bug (gives factor-of-N error).
