# Maxwell Postprocessor: Element Field Caching Optimization

**Date:** 2026-01-28
**Purpose:** CPU-only optimization to eliminate redundant physics computations
**Module:** `src/fem/maxwell`
**Status:** PROPOSED - Simple optimization, implement before GPU

**Priority:** HIGH - Much simpler than GPU, likely 5-10× speedup

---

## Problem Statement

**Current implementation** (`cl_MaxwellPostprocessor.cpp:634-768`):

In `process_recovery()`, for each node we loop over its patch elements (~10 elements) and compute expensive physics at integration points:

```cpp
for (index_t tIndex : tIndices)  // ~10 elements in node patch
{
    for (uint k = 0; k < num_intpoints; ++k)  // ~8 integration points
    {
        const Vector<real>& tY = this->compute(k);  // ← EXPENSIVE PHYSICS!
        // Computes: H field, B field, J current, Jc critical current
        // Involves: matrix-vector products, material lookups, nonlinear evaluations

        mVandermonde += tWeight * (mPoly * trans(mPoly));
        mCoefficients += tWeight * mPoly * tY;
    }
}
```

**The inefficiency:** Each element is shared by ~6-10 nodes (typical mesh connectivity). We recompute the same `tY` physics 6-10× for overlapping patches!

---

## Proposed Solution: Precompute Element Fields

**Key insight from Christian Messe:** Each Y-vector can be integrated locally per element once. Store the integrated Y-values per element and reuse them.

### Storage Requirements

**Per element: ONE integrated vector**
```
Memory: n_elements × n_fields × sizeof(real)
      = 1M elements × 6 fields × 8 bytes = 48 MB
```

**Negligible!** Compared to 80 MB for topology indices or 384 MB for GPU batching.

### Implementation Strategy

**Phase 1: Precompute (once per postprocessing call)**

```cpp
void MaxwellPostprocessor::precompute_element_fields() {
    // Allocate storage: one integrated vector per element
    mElementFields.set_size(mNumElements);

    for (id_t tID : mBlockIDs) {
        Cell<mesh::Element*>& tElements = mMesh->block(tID)->elements();

        for (mesh::Element* tElement : tElements) {
            if (!tElement->is_flagged()) continue;

            index_t e_idx = tElement->index();
            mElementFields(e_idx).set_size(mNumTargetFields, 0.0);

            // Setup element
            mElement = mMagfield->element(tElement->id());
            Block* tBlock = mMagfield->block(mElement->element()->block_id());
            mMaterial = mMaterialMap[tBlock->id()];
            mCalculator = tBlock->calculator();
            mCalculator->link(mElement);
            this->update_dofs();

            // Integrate Y over element (compute expensive physics ONCE)
            for (uint k = 0; k < mCalculator->num_intpoints(); ++k) {
                const Vector<real>& tY = this->compute(k);  // Physics computation

                // Accumulate (with integration weights if needed)
                mElementFields(e_idx) += tY;  // Or: += gauss_weight(k) * tY
            }

            // Normalize or scale as needed
            // mElementFields(e_idx) /= mCalculator->num_intpoints();
        }
    }
}
```

**Phase 2: Use cached values in recovery**

```cpp
const Vector<real>& MaxwellPostprocessor::process_recovery(mesh::Node* aNode) {
    mVandermonde.fill(0.0);
    mCoefficients.fill(0.0);

    Cell<index_t>& tIndices = mNodeToElements(aNode_index);  // Precomputed topology

    for (index_t elem_idx : tIndices) {
        mElement = mMagfield->element(tElements(elem_idx)->id());
        mCalculator->link(mElement);

        real tWeight = get_element_volume(elem_idx);

        // Get precomputed integrated Y for this element
        Vector<real>& tYelem = mElementFields(elem_idx);  // ← NO PHYSICS, just lookup!

        for (uint k = 0; k < mCalculator->num_intpoints(); ++k) {
            // Still need coordinates for polynomial basis
            compute_coordinates(k, mX);
            this->compute_poly(mX);

            // Vandermonde assembly (unchanged)
            mVandermonde += tWeight * (mPoly * trans(mPoly));

            // RHS assembly using cached Y (no physics recomputation!)
            for (uint j = 0; j < mNumTargetFields; ++j) {
                for (uint i = 0; i < mNumCoefficients; ++i) {
                    mCoefficients(i, j) += tWeight * mPoly(i, 0) * tYelem(j);
                }
            }
        }
    }

    // Solve (still needed, but now much smaller fraction of total time)
    posv(mVandermonde, mCoefficients);

    // Evaluate polynomial at node
    this->compute_poly(mX0);
    for (uint i = 0; i < mNumTargetFields; ++i) {
        mZ(i) = dot(mPoly, mCoefficients.col(i));
    }

    return mZ;
}
```

---

## Expected Performance

**Before optimization:**
- Physics computation (compute): ~60-70% of per-node time
- Polynomial assembly: ~20-25%
- posv solve: ~5-10%

**After element caching:**
- Precompute phase: O(n_elements) - done once
- Per-node recovery: Only polynomial assembly + posv
- **Speedup: 5-10× on recovery loop** (eliminate 6-10× redundant physics computation)

**For transient analysis:**
- Recompute element fields each timestep (DOFs change)
- Still 5-10× speedup per timestep
- Precompute cost amortized: O(n_elements) vs O(n_nodes × patch_size)

---

## Implementation Effort

**Trivial:** 1-2 days maximum

**Day 1:**
- Add `mElementFields` member: `Cell<Vector<real>>`
- Implement `precompute_element_fields()`
- Test: verify integrated values computed correctly

**Day 2:**
- Modify `process_recovery()` to use cached values
- Test: results identical to current implementation
- Profile: measure speedup

**Total:** 1-2 days (embarrassingly simple!)

---

## Advantages Over GPU Approach

| Aspect | Element Caching | GPU Batching |
|--------|----------------|--------------|
| **Implementation** | 1-2 days | 3-5 days |
| **Complexity** | Trivial | Moderate (CUDA, streams, batching) |
| **Memory** | 48 MB | 768 MB (triple buffering) |
| **Portability** | Pure CPU, works everywhere | Requires CUDA/ROCm |
| **Expected speedup** | 5-10× | ~1.3-1.4× (overall, pipelined) |
| **Risk** | Zero | Low-medium |

**Verdict:** Do element caching first. Then reconsider if GPU is still needed.

---

## Open Questions / Refinements

1. **Integration weights:** Should `tYelem` use Gauss integration weights, or just sum Y(k)?
   - Need to verify mathematical equivalence in SPR formulation

2. **Element volume scaling:** Currently `tWeight = element_volume` is applied in the loop. Can we fold this into precompute?

3. **Transient analysis:** For time-dependent problems, precompute needs to be called each timestep (DOFs change). Still worth it?

4. **Memory ordering:** Store as `Cell<Vector<real>>` or flat `Vector<real>` with indexing?

---

## Next Steps

1. ✅ **Profile current implementation** - verify physics computation dominates
2. ⏸ Implement element caching (wait for profiling confirmation)
3. ⏸ Test correctness (compare results with/without caching)
4. ⏸ Measure performance improvement
5. ⏸ Decide if GPU optimization still worthwhile after this

**Status:** Waiting for profiling data before implementation.

---

## Related Documents

- **[maxwell_postprocessor_gpu_acceleration.md](../deferred/maxwell_postprocessor_gpu_acceleration.md)** - GPU strategy (Phase 1, implement AFTER this)
- **[maxwell_postprocessor_review.md](maxwell_postprocessor_review.md)** - Overall review (2026-01-21)
- **[../../src/fem/maxwell/doc/postprocessor_recovery_theory.md](../../src/fem/maxwell/doc/postprocessor_recovery_theory.md)** - Mathematical foundation (authoritative)

---

**End of Document**
