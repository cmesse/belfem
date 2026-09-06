# Coil Interface Decoupling + Side-Local Field Recovery (Interfaces Crisp, Cuts Continuous)

**Date:** 2026-07-03
**Purpose:** Implement the coil decouple path (R3) and the postprocessor visualization split
(R5) from `todo/interface_node_duplication_coil_ferro.md`; fix the parallel ownership hole the
split exposed. Serial-validated by Christian on the tape/coil model (three ParaView checks).
**Modules:** `src/homology`, `src/fem/kernel`, `src/mesh`
**AIs involved:** Claude Opus (analysis + edits), Fable (continuation, recovery-consistency
fixes), Codex (parallel ownership audit)
**Status:** Serial-validated; parallel (R7/R8) pending.

---

## 1. Decisions (Christian)

- **Air-coil = full decouple (option A):** coil nodes duplicated, coil elements relinked like
  the conductor path, but duplicates carry **no sources** and are **not registered** in the
  original's `->duplicate()` container — invisible to solve and postprocessing. Coil domains
  and φ-coil interfaces stay excluded from the computation.
- **Ferro-air = weight-1 tie, split in viz only:** perpendicular B and in-plane H are
  physically continuous at the iron surface; in-plane B and perpendicular H jump (μ jump).
  Nodes tied in the Maxwell solve, decoupled in postprocessing.
- Master-side correction to Christian's lower-ID intuition confirmed in code:
  `fix_facet_masters()` (`cl_MaxwellFactory.cpp:1069-1106`) makes the **higher** `DomainType`
  master → ferro/coil side is duplicated; facet normal points out of the iron.

## 2. Changes

| File | Change |
|---|---|
| `cl_InterfaceProcessor.hpp` | `enum class InterfaceTreatment { TieWeight1, Decouple }`; `mTreatment` member + accessor on `InterfaceSet` |
| `cl_InterfaceProcessor.cpp` | ctor: any `DomainType::Coil` block on the interface ⇒ `Decouple`; `add_duplicates()`: decoupled dups skipped in `tOriginalMap` (⇒ never sourced); sourcing loop: single `find()`, skip on miss; `unify_duplicates()`: final `add_duplicate()` registration guarded by `is_duplicate()` — fixes a latent null-container write for sourceless duplicates (`original()` returns self, `add_duplicate` on unallocated `mDuplicates`) |
| `cl_FEM_Postprocessor.cpp` | Node claiming (both branches) and `recover_fields()` accumulation made **side-local** (raw `tElement->node(k)`, no `original()` canonicalization); then extended by the **same-PP sibling merge** (§3); `compute_node_matrices()` walk restored to original+duplicates (guard-filtered) |
| `cl_Mesh.cpp` | `update_ownerships()`: per-side ownership — each node (original AND duplicate) owned by min rank over **its own** elements; duplicate falls back to original's owner when its element list is empty. Previously the pair shared one union-min owner |

## 3. The recovery rule (final design)

> **A node's recovery patch = its own elements ∪ the elements of every registered sibling
> selected in the same postprocessor.**

Both halves of the projection use the same membership predicate — `is_flagged(0)` at
Vandermonde build (`compute_element_coeffs()` internal guard, `cl_FEM_Postprocessor.cpp:738`)
and `mNodeMatrices.key_exists()` at accumulation — provably the same set, so `V` and `B` stay
consistent by construction.

| Pair class | Sibling in same postprocessor? | Result |
|---|---|---|
| Cohomology cut (air-air) | yes | full-disc patch, identical fit both copies → continuous B/H |
| Ferro-air interface | no (other PP) | side-local → crisp physical jump |
| Coil duplicate | unregistered | untouched, fields stay at init (φ = NaN) |
| Plain node | fast path | unchanged |

Physics rationale: same-domain duplicates are **gauge cuts** (φ jump = transport current,
H continuous — the correct viz is the two-sided patch, also better conditioned); cross-domain
duplicates are **material interfaces** (fields genuinely jump — side-local is correct). The
distinction needs no new flags: postprocessor membership encodes it.

This resolves **O3** definitively: the "existing conductor-side split mechanism" hunted by the
plan never existed — the old `recover_fields()` wrote every element contribution to the
original + all registered duplicates, so both copies got cross-side-blended values and
last-writer-wins masked everything (including a long-standing interface blur).

## 4. Debugging arc (three iterations, each serial-tested by Christian)

1. **Side-local claiming + accumulation** → ferro-air crisp, coils dead ✓, but cuts showed a
   "broken projection": `compute_node_matrices()` still assembled `V` over original+duplicates
   while `B` was one-sided — for cut pairs (both sides air, same PP, all elements flagged)
   `inv(V_both)·B_one` is inconsistent. Interfaces were accidentally consistent (the element
   flag gate excluded the other domain from `V`).
2. **Vandermonde made side-local** → projection consistent, but B showed a jagged **seam along
   the cohomology cuts** (`cohomology.png`). Static analysis verified the whole input chain
   (cut dup sources = {λ abstract dofs, original} all weight 1, `cl_CutSet.cpp:47-76`;
   dof-level type-matched tie `cl_FEM_DofMgr_DofData.cpp:3515-3546`; field row =
   `node->index()` `cl_FEM_Dof.cpp:25`; raw per-element φ gather `cl_IWG.cpp:790`). Evidence
   (`FieldPhi.png`): φ jumps cleanly ⇒ inputs correct ⇒ the seam was the **one-sided fit
   itself** — half-disc patches of 1-2 huge far-field triangles are rank-marginal for the
   recovery polynomial (systematic offset, not noise).
3. **Same-PP sibling merge** (§3) → all three behaviors correct simultaneously. Confirmed by
   Christian ("Whohoo!").

## 5. Codex audit: parallel ownership hole (verdict: real, fixed)

Exchange: `tmp/ai_exchange/postproc_sidelocal_parallel_ownership.md` (full findings + Claude
verification). Key results (all citations re-verified against the tree):

- The distributor ghosts duplicate siblings and their elements (`cl_Mesh_Distributor.cpp:363-417`),
  but those arrive as **aura** elements, and FEM blocks split owned/aura
  (`cl_FEM_DofMgr_BlockData.cpp:111-133`); `Group::elements()` returns owned only
  (`cl_FEM_Group.hpp:505-517`) and aura FEM elements carry no dofs — aura cannot rescue
  side-local recovery.
- Under shared pair-ownership, a duplicate whose own-side elements live on a non-owner rank is
  claimed by **nobody** → its field row silently stays 0.0. Collectives stay protocol-safe
  (that's what makes it silent).
- Fix (Codex preference (a), applied): per-side ownership in `update_ownerships()`. Solver
  coupling checked: `compute_hanging_dofs()` is master-only over global data
  (`cl_FEM_DofMgr_SolverData.cpp:2526`), and thin-shell duplicates (never `set_original()`ed)
  already live under per-side ownership — not a new regime. Soft consequence accepted by
  Christian: some edge/facet owners shift → different but consistent dof partition.

## 6. Evidence (serial, Christian's tape/coil model)

- `cohomology.png` — the intermediate cut seam (iteration 2), since resolved.
- `FieldPhi.png` — φ gauge regions with clean jumps across cuts (inputs correct).
- `FieldB.png` — B magnitude continuous across cuts, crisp at interfaces, coils dead.

## 7. Open items

- **R4 re-scoped:** the ferro strip (`cl_MaxwellFactory.cpp:1882-1900`) is still active in the
  working tree, yet the serial run is physical — likely because dof-level T-matrices are built
  from mesh-level sources during kernel init, *before* the strip runs, so it only clears
  mesh-level bookkeeping. Verify the ordering before deciding to remove or keep it.
- **R7/R8:** 2/4-proc validation — specifically the per-side ownership fix (failure signature
  it prevents: silent zero stripe on interface/cut duplicate rows) and solve-identical
  residual histories.
- **R6:** triple-junction guard (node already duplicated by `CutProcessor` hit by
  `InterfaceProcessor`) still unimplemented.
- **O4:** coil nodal fields other than φ default to 0.0, not NaN (`cl_Mesh_Field.cpp:40`);
  full-NaN "super clean" option still optional/waived.
- **O5:** Buffer≡Air membership in the admission predicate undecided.
