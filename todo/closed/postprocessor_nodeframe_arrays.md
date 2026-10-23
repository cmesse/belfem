# Postprocessor `NodeFrame`: fixed-size arrays instead of `Vector< real >`

**Date:** 2026-10-22
**Purpose:** The recovery postprocessor keeps one `NodeFrame` per selected node for the whole run,
and the struct held its origin, normal and two tangents as four `Vector< real >`. Those are 2- or
3-component geometric quantities; the only operations on them are element access, `+=`, `/=`, two
dot products and a norm. That is the container misuse CLAUDE.md names ("`Vector`/`Matrix` are for
linear algebra ONLY"). Under Armadillo each `Vector` is an `arma::Mat` with a 16-double local
buffer, so the struct was 784 B (probe, 2026-10-22; 112 B under Blaze, which instead heap-allocates
every sized member) and was itself heap-allocated once per node per postprocessor. The fix stores
the four quantities as `real[3]`, writes the two dot products and the norm as plain loops, and
stores the frames by value.
**Module:** `src/fem/kernel` (`cl_FEM_Postprocessor.hpp/.cpp`), `tests/fem`
**AIs involved:** Claude Fable 5.1 (exploration, plan, implementation, gates), Codex gpt-6-astra/high and Grok grok-4.7/high (plan jury and code jury)
**Status:** ✅ COMPLETE 2026-10-22 — R1–R5 landed; `NodeFrame` is 14 doubles, trivially copyable, held by value in `Map< mesh::Node *, NodeFrame >`; four 2D recovery cases added (R8b). Verified in Claude's private sandbox trees (Debug, GNU 16.1.0, SCLS): positive control of the new cases on the pre-change code 15/15; change 15/15 on Blaze and on Armadillo; whole `test_fem` 337/337 on Armadillo; mutation `T1 := Normal` in 2D turns three new cases red; restore 15/15. Residual: O5 (oblique-normal 3D fixture, tilted-T1 2D case) and R8c (deck-level VmHWM) not done; D1 (pre-existing re-`initialize()` crash after a throw) recorded, not fixed.

> **Scope guards:**
> - In scope: the `NodeFrame` struct, its storage map, and every access to its fields:
>   `compute_node_frame`, `compute_local_poly`, `accumulate_node_rhs`, plus the allocate/free
>   sites in `compute_node_matrices` / `delete_node_matrices` and the `add_node_gram` call.
> - Out of scope: the per-node inverse `Map< mesh::Node *, Matrix< real > * > mNodeMatrices`
>   (same storage pattern, pre-dates 0.9.3 — O4); the scratch members `mOffset`,
>   `mShellDirection`, `mLocalPoint` (sized once in `select_recovery_mode`, accepted
>   member-scratch pattern, and `compute_poly` takes a `Vector`); the recovery mathematics
>   (the 2D `T1` and 3D `T2` are deliberately not normalised — do not "fix" them); the
>   pre-existing re-`initialize()` lifecycle defect (D1); the `mesh::Vertex` counter width
>   (`tmp/ai_exchange/review_memory_092_093.md`).
> - Results move only at round-off (see §1 "what moves" for the honest statement).

---

## 1. Current Behaviour and How It Failed

`cl_FEM_Postprocessor.hpp:80-90` before the change, in the class's `private:` section (`:31`):

```cpp
struct NodeFrame
{
    Vector< real > Origin ;
    Vector< real > Normal ;
    Vector< real > T1 ;
    Vector< real > T2 ;
    real H   = 1.0 ;
    real Tau = 1.0 ;
};
Map< mesh::Node * , NodeFrame * > mNodeFrames ;
```

| Failure | Mechanism | Evidence |
|---|---|---|
| Container misuse | four `Vector< real >` for fixed-size geometry. The hot `dot` on them was already Armadillo's inline generic loop for `n ≤ 32` (`op_dot_meat.hpp:156`); only the once-per-node `norm` reached BLAS (`nrm2`, `op_norm_meat.hpp:559-584` via `fn_norm.hpp:35-41`). Nothing here *needs* the linear-algebra layer | `cl_FEM_Postprocessor.cpp:924-1027`, `:1075-1084` (pre-change lines) |
| 784 B per node under Armadillo | `Vector< real >` wraps `arma::Mat< double >` 1:1 (`cl_AR_Vector.hpp:38,45`), 192 B with its `mem_local[ 16 ]` (`Mat_bones.hpp:51`); 4 × 192 + 16 = 784. Blaze: `DynamicVector` is 24 B (`DynamicVector.h:531-541`), 4 × 24 + 16 = 112, but every `set_size` heap-allocates (`cl_BZ_Vector.hpp:271` → `DynamicVector.h:1751`) | sizeof probe 2026-10-22, standalone against `/opt/scls/include` |
| Heap chunks per node | Armadillo: `new NodeFrame()` + the `unordered_map` node = 2. Blaze: those 2 + one per sized member — non-shell 3, 2D shell 5, 3D shell 6 | pre-change `:693-694`, `:924`, `:932`, `:998`, `:1024` |
| Persistent for the run | `accumulate_node_rhs` reads the frame at every recovery; freed only in `delete_node_matrices()` from the destructor or before a refill | `:1122-1125`, `:719-732`, `:59`, `:686` |
| Multiplied by postprocessor instances | Air always; Ferro under h-φ with ferro blocks; Conductor / SuperConductor; ThinShellConductor / ThinShellSuperConductor; SideConnector — one instance per element type, each with its own map of *its* selected nodes | `cl_MaxwellFactory.cpp:2751-2767`, `:2816-2820`, `:2836-2845`, `:2899-2915` |
| Three of four vectors unused outside shell mode | `Normal`, `T1`, `T2` were `set_size`d only under `mShellLocal`; an unsized `arma::Mat` still costs its 192 B | `:932`, `:998`, `:1024` |

Magnitude (estimate, medium confidence): ~0.5–1 GB resident on a 3.5 M-dof 3D run under
Armadillo (~0.6 M nodes, the Air instance covering most of them, ~0.7 KB saved per frame). The
largest code-side memory delta found between v0.9.2 and v0.9.3 (introduced by `5ab4d41`,
`todo/closed/thinshell_nodal_recovery.md:11`) — single-raiser finding of a one-auditor round,
`tmp/ai_exchange/review_memory_092_093.md`. Not the dominant term of such a run (the direct
factorization is), but a guideline violation with a measurable cost and a local fix.

**What moves.** Three things, all at round-off for well-formed patches:
1. The hot dot products: Armadillo's `n ≤ 32` path is already a sequential loop in index order,
   so a hand loop is the same arithmetic; Blaze may pair terms differently.
2. `norm( Normal )` and `norm( T1 )`: BLAS `dnrm2` scales by the largest component before
   summing squares; a plain `sqrt( Σ x² )` does not. For a sum of a few dozen unit directions
   overflow/underflow is impossible, so this is an **algorithm change with round-off
   consequences**, not merely summation order. The patch-radius loop already used the plain form.
3. A frame perturbed at round-off can move a singular value across the hard cut in
   `recovery_block_inverse` (`fn_FEM_recovery_inverse.hpp:108`, `:119`) on a near-degenerate
   patch, so final-field differences are **not universally bounded to the last bit**. The
   recovery tests compare to manufactured fields with `EXPECT_NEAR` tolerances of 1e-13…1e-9.

**Bottom line:** four `Vector< real >` where `real[3]` belongs cost ~800 B and 2–6 heap chunks per
node per postprocessor; the fix is one struct and the three functions that touch it.

## 2. Architecture

`real[3]` is the house pattern for a fixed-size geometric quantity: `mesh::Node` stores its
coordinates as `real mCoords[3]` (`cl_Node.hpp:32`). The struct is 14 doubles (112 B), trivially
copyable, no heap, backend-independent. The dot products and the norm are written out in place —
the cross products in `compute_node_frame` were *already* written out component-wise and the patch
radius was already a plain sum of squares.

The justification is persistent footprint and allocation count. Runtime: the hot path replaces one
3-term loop by another; it introduces no allocation and no call, and is no slower in any way the
source shows — **unmeasured** (`doc/coding_philosophy.md:765`, "Benchmarks We Owe Ourselves").

Rejected: keeping `Vector< real >` and shrinking by `set_size( 0 )` — the Armadillo header is 192 B
regardless; `std::array< real, 3 >` — equivalent (`typedefs.hpp:71` has a `std::array` precedent),
but `Node::mCoords` is the closer one and the choice should not fork.

## 3. Gap Table

| # | State / behaviour | Needed for | Was | Class | Now |
|---|---|---|---|---|---|
| 1 | `Origin` | every fit | `Vector`, sized `mNumDimensions` | (c) | `real[3]` |
| 2 | `Normal` | shell mode | `Vector`, zero-sized in non-shell | (c) | `real[3]`, zero |
| 3 | `T1`, `T2` | shell mode | `Vector`; `T2` only in 3D | (c) | `real[3]`, zero |
| 4 | `H`, `Tau` | every fit / shell | `real` | (a) | unchanged |
| 5 | `Normal += mShellDirection` | shell mode | `Vector` op | (c) | loop to `mNumDimensions`; `mShellDirection` stays `Vector` |
| 6 | `norm( Normal )`, `Normal /= l`, `T1 /= norm( T1 )` | shell mode, once per node | `arma::norm` → `nrm2` | (c) | plain loops; 3D `T1` normalised from its two written components |
| 7 | `dot( mOffset, T1 )`, `dot( mOffset, T2 )` | shell mode, **hot** | Armadillo generic loop | (c) | one loop to `mNumDimensions` accumulating `tS1`, `tS2`; `mOffset` stays `Vector` (`compute_poly` consumes it) |
| 8 | storage | lifetime of the postprocessor | `Map< Node *, NodeFrame * >`, `new`/`delete` | (c) | `Map< Node *, NodeFrame >` by value (O1) |
| 9 | lookup in `accumulate_node_rhs` | recovery | `*mNodeFrames( aTarget )` | (c) | `const NodeFrame & tFrame = mNodeFrames( aTarget )`; `Map::operator()` returns a reference (`cl_Map.hpp:226-242`) |
| 10 | frame cleanup before a refill | no leak | `delete_node_matrices()` deleted frames | (a) | `mNodeFrames.clear()`; **the retry path itself is broken independently** (D1) |
| 11 | 2D | `Normal[2]`, `T1[2]`, `T2[*]` unused | n/a | (c) | contract: a loop that pairs a component with mesh data runs to `mNumDimensions` (the `Vector` operands are sized to it); the zeroing loop and the 3D axis search run to 3; the arrays' zero third component is **load-bearing for T1** (the 3D branch writes `tA`, `tB` only) |
| 12 | documentation | `postprocessor_recovery_theory.md` describes the frame mathematically, not its storage | (a) | no change |

## 4. Ordered Steps

- [x] **R1** — `NodeFrame` with `real Origin[ 3 ] = { 0.0, 0.0, 0.0 }` etc.; `H`, `Tau`;
  `Map< mesh::Node *, NodeFrame > mNodeFrames`; `#include <type_traits>`;
  `static_assert( std::is_trivially_copyable< NodeFrame >::value, … )` after the brace; the
  contract comment.
- [x] **R2** — `compute_node_matrices`: `NodeFrame & tFrame = mNodeFrames[ tNode ]`, no `new`,
  `add_node_gram( tNode, tFrame )`. `delete_node_matrices`: `mNodeFrames.clear()`; the
  `mNodeMatrices` delete loop stays.
- [x] **R3** — `compute_node_frame`: every field access with `[ ]`; explicit zero of
  Normal/T1/T2 in shell mode at entry; `Normal +=` loop; plain norm and divide loops; 2D
  rotation; 3D T1 from two components, two-term norm; T2 cross product. `BELFEM_ERROR` on the
  zero normal kept.
- [x] **R4** — `compute_local_poly`: `Origin[ d ]`; `tS1`/`tS2` loops.
- [x] **R5** — `accumulate_node_rhs`: bind by reference. Assert kept.
- [x] **R6** — inline loops, no named helpers (O2).
- [x] **R7** — `-fsyntax-only -std=gnu++17 -Wall -Werror` with the debug tree's flags on both
  TUs: clean.
- [x] **R8 — executable gate (3D)** — `PostprocessorRecovery.*` 15/15 on Blaze/Debug and
  Armadillo/Debug (private sandbox trees); whole `test_fem` 337/337 on Armadillo.
- [x] **R8b — executable gate (2D)** — `Shell2D( aCurved )` (two QUAD4TS layers, LINE2 tape
  facets, `mesh::ThinShell`; arc or flat tape) and `FarTriPatch` (8 TRI3, Air, 1e3 from the
  origin) added; cases `ThicknessJumpSurvivesCurvature2D`, `ConstantReproducedOnCurvedShell2D`,
  `TangentLinearFieldReproducedOnFlatShell2D`, `LinearFieldFarFromOrigin2D`. Positive control: all
  four pass on the pre-change code. Mutation `T1 := Normal` in the 2D branch: the first three go
  red (6, 6 and 15 truncated patches, plus value failures on the flat tape), the fourth and the
  3D cases stay green. Christian's ruling 2026-10-22: the fixture helper takes **no**
  dimensionality parameter; it derives `ModelDimensionality` from `mMesh->number_of_dimensions()`.
- [ ] **R8c (optional, not done)** — deck-level VmHWM after `Postprocessor::initialize()`.
- [x] **R9** — devlog `devlog/dl20261022_nodeframe_arrays_and_memory_hunt.md`, README lines,
  this file retired to `closed/`.

### 4.0 Implementation record (2026-10-22)

**Implemented and audited.** The struct, the map and the three consumers as in §4. Plan jury
(Codex gpt-6-astra/high, Grok grok-4.7/high): two compile-breaking omissions in R2/R4 (`*tFrame` at
the `add_node_gram` call, `Origin( )` at two sites), the 2D gate gap, the `nrm2` algorithm change,
the Blaze/O4 accounting, six line-citation corrections, D1 — all folded in before implementation.
Code jury (same seats): no blocking defect; four overclaiming comments (mine) corrected; the
`tS2`-in-2D nit left as a single loop; O5 logged.

- **D1 — HIGH (pre-existing, out of scope, not fixed).** `initialize()` deletes and nulls
  `mNodeBitset` (`:105-106`) and `mElementBitset` (`:128-129`) *before* `compute_node_matrices()`
  (`:132`), and sets `mIsInitialized` only at `:136`. If `compute_node_matrices` throws (a
  `BELFEM_ERROR` on a degenerate patch), a second `initialize()` re-enters
  `select_elements_and_owned_nodes()` and dereferences the null `mElementBitset` (`:228`). The
  comment at `:685` describes a path that crashes earlier. Found by Codex, confirmed by Claude.
  The by-value frames are leak-free under unwinding (the map clears from the destructor) but
  establish no retry safety.

**Still open by design.** R8c; O4; O5.

## 5. Open Design Questions

- **O1 — RESOLVED 2026-10-22 → by value in the `Map`** (both auditors concur). A `Cell` indexed by
  a per-postprocessor node slot would remove the map node too, but there is no such slot: the node
  index is borrowed during `compute_node_matrices` and restored by `update_node_indices()`.
- **O2 — RESOLVED 2026-10-22 → inline loops**, no named helpers.
- **O3 — RESOLVED 2026-10-22 → `static_assert( is_trivially_copyable )`** after the brace.
- **O4 — `mNodeMatrices` (open, logged).** One `Matrix< real >` per node. Shell mode is 4×4 = 16
  doubles, which *fits* Armadillo's `mem_local[ 16 ]` (`arma_config.hpp:29`, `Mat_meat.hpp:401`) —
  the cost is the ~192 B header plus the `new` and the map node. Volume order 2 is 10×10 and does
  heap-allocate. A flat `Cell< real >` plus a slot map would remove the header; ratio unmeasured.
  Separate plan if the memory matters after R8c.
- **O5 — recovery-frame coverage (open, logged by the code jury).** Both shell fixtures are
  axis-aligned extrusions (3D radial about z; 2D arc/flat), so the least-aligned-axis branches of
  the 3D tangent construction are exercised for one axis only, and a 2D `T1` *tilted* by a modest
  angle passes every case (the curved 2D cases catch a `T1` parallel to the normal through the
  truncation counter, the flat case catches a zero or normal-parallel `T1`; none pins the
  rotation). A rotated 3D fixture with manufactured fields along both tangents, and a 2D case that
  reads the tangent's direction, would close it. Raised by Codex and Grok from different angles.

## 6. Interface

No interface change. `NodeFrame` is a **private** nested type (`hpp:31`); `grep -rn NodeFrame src`
finds only `cl_FEM_Postprocessor.{hpp,cpp}`. `MaxwellPostprocessor::initialize()` only calls the
base (`cl_MaxwellPostprocessor.cpp:298-300`); `Gradient` never calls `initialize()`.

## 7. Definition-of-Done

- [x] Every gap-table row mapped to a step (1–11 → R1–R5; 12 → none needed).
- [x] `grep -n "Vector< real >" cl_FEM_Postprocessor.hpp` shows no `Vector` inside `NodeFrame`.
- [x] `PostprocessorRecovery.*` green including the new 2D cases, both backends (R8, R8b).
- [x] Layout: 14 doubles = 112 B by construction (`static_assert` on triviality in the header).
- [x] D1 recorded here; debt-register row at Christian's call.

## 8. Audit Trail

- Memory hunt that found it: `tmp/ai_exchange/review_memory_092_093.md` (one-auditor round, Codex
  quota; Grok grok-4.7/high).
- Plan and code rounds: `tmp/ai_exchange/review_nodeframe_arrays.md` — Codex gpt-6-astra/high,
  Grok grok-4.7/high on both. Every citation re-verified by Claude against the tree before
  inclusion; the only refutation in either round was of Claude's own fixture comment (the mesh
  checker never reorders a QUAD4TS: `Pipette::measure_quad4ts` returns `std::abs`).
- Distilled into `devlog/dl20261022_nodeframe_arrays_and_memory_hunt.md`.
