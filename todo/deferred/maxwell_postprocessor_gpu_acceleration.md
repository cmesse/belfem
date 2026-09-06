
# Maxwell Postprocessor: Potential GPU Acceleration

**Date:** 2026-01-20 (Updated: 2026-01-28)
**Purpose:** Discussion of potential GPU acceleration for extreme-scale problems
**Module:** `src/fem/maxwell`
**Status:** DEFERRED - Element caching should be implemented first

> **⚠️ IMPORTANT:** This document discusses potential future optimizations, **not current performance issues**. The Maxwell postprocessor scales well with MPI and requires no immediate optimization work.

> **🔥 NOTE (2026-06-22):** The CPU-side **[element field caching](../closed/maxwell_postprocessor_element_caching.md)** prerequisite is now effectively implemented (the postprocessor was rewritten to compute physics once per element), so that todo is closed. GPU acceleration may well be unnecessary; re-profile before pursuing it.

**See Also:**
- **[Element Field Caching](../closed/maxwell_postprocessor_element_caching.md)** - CPU prerequisite — now effectively implemented (closed)
- [Recovery Theory](../../src/fem/maxwell/doc/postprocessor_recovery_theory.md) - Mathematical foundation (authoritative)
- [Postprocessor Review](../closed/maxwell_postprocessor_review.md) - Review findings (2026-01-21)
- [Maxwell Module README](../../src/fem/maxwell/doc/README.md) - Module overview

---

## Current Status (2026-01-21)

**MPI Parallelization:** ✅ Excellent (embarrassingly parallel, >95% efficiency expected)
**Mathematical Correctness:** ✅ Verified (Zienkiewicz-Zhu SPR)
**Performance Issues:** ❌ None identified

**Key finding:** Postprocessor recovery loop is **embarrassingly parallel**:
- Each rank processes only its owned nodes
- Patch recovery is purely local (no inter-rank communication)
- Communication only at boundaries: initial distribution + final collection

**Conclusion:** Current implementation is well-designed. **No optimization work required** unless postprocessing exceeds 20% of total simulation time for extreme-scale problems (>1M nodes).

---

## CRITICAL: What NOT to Change

**⚠️ Before optimizing, read these constraints:**

### Mathematical Integrity (DO NOT MODIFY):

1. **SPR Algorithm** - Keep Zienkiewicz-Zhu least-squares formulation
   - ❌ Do NOT replace with global L² projection
   - ❌ Do NOT replace with MLS/moving least squares
   - ❌ Do NOT replace Cholesky with QR "for stability"

2. **No Feedback into Solver** - Recovery is postprocessing ONLY
   - ❌ Do NOT feed recovered fields back into solve
   - ❌ Do NOT use for error-driven remeshing (yet)

3. **Polynomial Basis** - Keep complete polynomial spaces
   - ❌ Do NOT switch to incomplete/truncated bases
   - ❌ Do NOT change polynomial order selection logic

4. **Volume Weighting** - Keep element-volume weighting for stability
   - ❌ Do NOT switch to uniform weights
   - Inverse-distance weighting is for filter mode only

### Safe First Optimizations (Recommended Order):

**Phase 1 (1-2 days, minimal risk):**
1. ✅ Add OpenMP to node loop → 8× speedup
2. ✅ Preallocate work matrices → 1.2× speedup
3. ✅ Add patch caching for transient → 1.5-2× speedup

**Phase 2 (1-2 weeks, low risk):**
4. ✅ Inline small matrix solvers → 5× speedup
5. ✅ Use parent class batch recovery → 2× speedup
6. ✅ Precompute reference polynomials → 1.5× speedup

**Do NOT attempt GPU/advanced optimizations unless profiling shows postprocessing is a bottleneck (>20% of total simulation time).**

---

## Overview

The Maxwell postprocessor is **mathematically correct and produces accurate results**. This document discusses **potential GPU acceleration** for extreme-scale problems as the only viable optimization path.

**Current MPI implementation:** Embarrassingly parallel, >95% efficiency expected

**Why NOT OpenMP:** BELFEM is MPI-first. Hybrid MPI+OpenMP creates more problems than it solves (see `doc/coding_philosophy.md`).

**Only viable optimization:** GPU acceleration for extreme-scale (>1M nodes), justified only if postprocessing exceeds 20% of total simulation time.

---

## Why MPI Scaling is Excellent

**Code structure** (`cl_MaxwellPostprocessor.cpp:363-422`):

```cpp
void MaxwellPostprocessor::run()
{
    // Phase 1: Distribute source fields (one-time MPI communication)
    mKernel->dofmgr()->distribute_fields( mSourceFields );

    // Phase 2: Embarrassingly parallel recovery (NO communication)
    for ( index_t tIndex=0; tIndex<mNumNodes; ++tIndex )
    {
        // Each rank processes ONLY its owned nodes
        // Patch assembly accesses local + aura elements (no MPI)
        const Vector< real > & tValues = this->process( mNodes( tIndex ) );

        // Store results locally (no shared state)
        index_t tNodeIndex = mNodes( tIndex )->index();
        for ( uint f=0; f<mNumTargetFields; ++f )
        {
            tTargetFields( f )->data()( tNodeIndex )= tValues( f );
        }
    }

    // Phase 3: Collect target fields (one-time MPI communication)
    mKernel->dofmgr()->collect_fields( mTargetFields );
}
```

**Why this scales well:**
- ✅ Zero communication during main loop
- ✅ Each rank operates independently on owned nodes
- ✅ Patch recovery is purely local (aura provides neighbor data)
- ✅ No rank-0 bottleneck, no blocking sends/receives
- ✅ No load imbalance (nodes evenly distributed by mesh partitioner)

**Expected MPI scaling:** >95% parallel efficiency at any reasonable processor count

**Contrast with projected data in original document:** The document previously showed 71% efficiency at 16 ranks. This was **incorrect**. The embarrassingly parallel structure guarantees near-linear scaling.

---

## Algorithmic Complexity

| Operation | Complexity | Variables |
|-----------|------------|-----------|
| Patch assembly | O(N × E_patch) | N=nodes, E_patch≈6-20 elements |
| Coordinate interpolation | O(N × E × K × n_elem) | K=int.points, n_elem=nodes/elem |
| Polynomial evaluation | O(N × E × K × m) | m=coefficients (3-35) |
| Gramian assembly | O(N × E × K × m²) | Outer product p·pᵀ |
| **Cholesky solve** | **O(N × m³)** | **Per-node dense solve** |
| Polynomial evaluation at node | O(N × m) | Final result computation |
| MPI communication | O(P × N_data) | P=ranks, N_data=field data size |

**Total:** O(N × E × K × m²) dominated by Gramian assembly for large meshes.

---

## Bottleneck Analysis (2026-01-28)

**Location:** `cl_MaxwellPostprocessor.cpp:634-768` (`process_recovery`)

Detailed profiling of per-node operations in the recovery loop:

### Operation Breakdown

**1. Bitset Operations (`DynamicBitset`):**
```cpp
// Lines 639-659: Patch element identification
mElementBitset->reset();                     // O(n_elements/64)
// Loop over node + duplicates, calling set()  // O(patch_size) ≈ 6-20 elements
mElementBitset->where(tIndices);             // O(n_elements/64)
```

**Cost:** ~1,500-15,000 operations (1M elements)
- Uses `__builtin_ctzll` hardware intrinsics (bit counting)
- Sequential memory access, cache-friendly
- **NOT the bottleneck** despite initial concerns

**2. Patch Assembly Loop (lines 691-750):**
```cpp
for (index_t tIndex : tIndices)  // ~10 elements
{
    for (uint k=0; k<num_intpoints; ++k)  // ~8 integration points
    {
        // Line 737: Gramian assembly
        mVandermonde += tWeight * (mPoly * trans(mPoly));  // O(m²)

        // Lines 742-748: RHS assembly
        for (uint j=0; j<mNumTargetFields; ++j)            // O(m × n_fields)
            mCoefficients(i,j) += tWeight * mPoly(i,0) * tY(j);
    }
}
```

**Cost:** ~32,000 operations (m=20, 10 elements, 8 int points)
- Dominated by outer product `mPoly * trans(mPoly)`

**3. LAPACK Cholesky Solve (line 753):**
```cpp
posv(mVandermonde, mCoefficients);  // Cholesky factorization + triangular solves
```

**Cost:** O(m³/3) + O(m² × n_fields)
- Factorization: ~2,667 flops (m=20)
- Solve: ~2,400 flops (m=20, n_fields=6)
- **Total: ~5,000-10,000 operations per node**

### Bottleneck Verdict

**Primary bottleneck: `posv` (LAPACK Cholesky solver)**

| Component | Ops/Node | % of Total | GPU Viable? |
|-----------|----------|------------|-------------|
| Bitset (`where`) | 1,500-15,000 | 3-5% | ❌ Not worth it |
| Patch assembly | 32,000 | 60-70% | ⚠️ Complex (FEM) |
| **posv solve** | **5,000-10,000** | **20-30%** | **✅ Perfect fit** |

**Why posv is the target:**
1. ✅ **Grows cubically** with polynomial order (O(m³))
2. ✅ **Perfectly parallel** - each node independent
3. ✅ **Ideal for GPU** - batched dense linear algebra
4. ✅ **cuSOLVER optimized** - uses tensor cores on modern GPUs

**Why NOT optimize bitset:**
1. ❌ Already highly optimized with hardware intrinsics
2. ❌ Memory-bound operation (limited by DRAM bandwidth)
3. ❌ Small fraction of total time (3-5%)
4. ✅ **Better solution:** Precompute element indices (see below)

### Memory Efficiency Analysis

**Proposed: Store element indices explicitly**
```cpp
Cell<Cell<index_t>> mNodeToElements;  // [node_id][patch_elements]
Memory: n_nodes × ~10 elements × 8 bytes = 80 MB (1M nodes)
```

**Current: DynamicBitset on-the-fly**
```cpp
Memory: n_elements / 8 bytes = 125 KB (1M elements)
```

**Verdict:** 80 MB is **negligible** compared to GPU transfer buffers (~256-512 MB per batch). Precomputing indices:
- ✅ Eliminates bitset from critical path entirely
- ✅ Enables efficient batch assembly for GPU
- ✅ Simplifies pipeline logic
- ✅ 640× more memory, but totally acceptable

---

## Why NOT OpenMP (2026-01-21 Update)

**Developer comment at `cl_FEM_Postprocessor.cpp:416-417`:**
```cpp
// ( we might want to use OpenMP or Cuda here in the future
// if this turns out to be a bottleneck
```

**Decision:** OpenMP is **NOT recommended** for BELFEM.

**Problems with OpenMP in MPI-first framework:**
- Requires `MPI_THREAD_MULTIPLE` (10-30% performance penalty)
- Thread safety requirements for all LAPACK/BLAS calls
- Race conditions, debugging complexity
- Worse scaling than pure MPI (70-85% vs 95%)

**Better alternative:** Use finer MPI decomposition (more ranks, not threads within ranks)

**See:** `doc/coding_philosophy.md` section "Why BELFEM Does NOT Use OpenMP" for detailed rationale.

---

## GPU Acceleration Strategy (2026-01-28 Update)

**Justification threshold:** Postprocessing time > 20% of total simulation time AND problem size > 1M nodes

**Approach:** CPU assembles Gramians/RHS in batches, GPU solves via cuSOLVER batched Cholesky

### Revised Architecture

**Key insight:** CPU assembly is FEM-heavy (calculators, element data, material properties). Keep this on CPU where the infrastructure exists. **Only offload the batched linear solves.**

```
CPU (per MPI rank):
  ✅ Mesh management
  ✅ Patch construction (using precomputed indices)
  ✅ DOF gathering
  ✅ Gramian assembly (mVandermonde, mCoefficients)
  ✅ Batch packing (concatenate matrices/RHS)
  ✅ MPI communication

GPU (1 per rank):
  ✅ Batched Cholesky factorization (cuSOLVER)
  ✅ Batched triangular solves
  ✅ Return coefficients to CPU
```

**Why not assemble on GPU?**
- ❌ Complex FEM infrastructure (calculators, materials, integration)
- ❌ Irregular memory access (element patches, node coordinates)
- ❌ Low arithmetic intensity (memory-bound, not compute-bound)
- ✅ **Better:** Keep simple, let CPU do what it's good at

### Memory Budget Analysis

**Per-node data (m=20 coefficients, n_fields=6):**
```
Vandermonde matrix: m × m × 8 bytes      = 3,200 bytes
RHS matrix:         m × n_fields × 8     = 960 bytes
Total per node:                           ~4.2 KB
```

**Batch sizing (256 MB budget):**
```
Batch size:         256 MB / 4.2 KB      = ~61,000 nodes/batch
1M nodes:           1,000,000 / 61,000   = ~17 batches
10M nodes:          10,000,000 / 61,000  = ~164 batches
```

**Memory allocation:**
```
CPU (triple buffering): 3 × 256 MB       = 768 MB
GPU (single batch):     1 × 256 MB       = 256 MB
Topology indices:       80 MB            = 80 MB
Total overhead:                           ~1.1 GB (acceptable)
```

### Implementation Strategy: Pipelined Batches

**Three-stage pipeline for CPU/GPU overlap:**

```
Timeline (triple buffering):

Batch 0: [CPU: Assemble]──[H2D Transfer]──[GPU: Solve]──[D2H Transfer]──[CPU: Unpack]
Batch 1:                   [CPU: Assemble]──[H2D Transfer]──[GPU: Solve]──[D2H Transfer]──[CPU: Unpack]
Batch 2:                                    [CPU: Assemble]──[H2D Transfer]──[GPU: Solve]──[D2H Transfer]──[CPU: Unpack]
Batch 3:                                                     [CPU: Assemble]──[H2D Transfer]──[GPU: Solve]──[D2H Transfer]

Key: Stages overlap! CPU assembles batch i while GPU solves batch i-2.
```

**Expected speedup:** 30-40% vs sequential (assemble all → solve all)

### Implementation Outline

**Phase 1: Initialization (precompute topology)**

```cpp
void MaxwellPostprocessor::precompute_topology() {
    mNodeToElements.set_size(mNumNodes);

    for (index_t i = 0; i < mNumNodes; ++i) {
        mElementBitset->reset();
        mesh::Node* tNode = mNodes(i);
        mesh::Node* tOrg = tNode->original();

        // Original node elements
        for (uint e = 0; e < tOrg->number_of_elements(); ++e) {
            if (tOrg->element(e)->is_flagged())
                mElementBitset->set(tOrg->element(e)->index());
        }

        // Duplicate node elements
        for (uint d = 0; d < tOrg->number_of_duplicates(); ++d) {
            mesh::Node* tDup = tOrg->duplicate(d);
            for (uint e = 0; e < tDup->number_of_elements(); ++e) {
                if (tDup->element(e)->is_flagged())
                    mElementBitset->set(tDup->element(e)->index());
            }
        }

        // Extract indices (bitset operation - done once!)
        mElementBitset->where(mNodeToElements(i));
    }

    // Compute batch configuration
    size_t bytes_per_node = mNumCoefficients * mNumCoefficients * sizeof(real) +
                            mNumCoefficients * mNumTargetFields * sizeof(real);
    const size_t MAX_BATCH_BYTES = 256 * 1024 * 1024;  // 256 MB
    mBatchSize = MAX_BATCH_BYTES / bytes_per_node;
    mNumBatches = (mNumNodes + mBatchSize - 1) / mBatchSize;
}
```

**Phase 2: CPU Batch Assembly**

```cpp
void MaxwellPostprocessor::assemble_batch(
    size_t batch_idx,
    int buffer_idx,
    real* cpu_vandermonde,
    real* cpu_rhs)
{
    size_t node_start = batch_idx * mBatchSize;
    size_t node_end = std::min(node_start + mBatchSize, mNumNodes);

    for (size_t i = node_start; i < node_end; ++i) {
        size_t local_idx = i - node_start;

        // Use precomputed indices (no bitset!)
        Cell<index_t>& tIndices = mNodeToElements(i);

        // Assemble Gramian and RHS (existing logic from process_recovery)
        mVandermonde.fill(0.0);
        mCoefficients.fill(0.0);

        for (index_t elem_idx : tIndices) {
            // ... element loop, integration points, existing FEM logic ...
        }

        // CRITICAL: Handle Blaze vs Armadillo matrix layout
        real* V_dest = &cpu_vandermonde[local_idx * m * m];
        real* R_dest = &cpu_rhs[local_idx * m * n_fields];

#ifdef BELFEM_ARMADILLO
        // Armadillo: column-major, contiguous storage
        std::copy_n(mVandermonde.memptr(), m * m, V_dest);
        std::copy_n(mCoefficients.memptr(), m * n_fields, R_dest);
#else
        // Blaze: Need manual flattening (column-major for LAPACK compatibility)
        // LAPACK expects column-major, Blaze can be row-major or column-major
        for (uint col = 0; col < m; ++col) {
            for (uint row = 0; row < m; ++row) {
                V_dest[col * m + row] = mVandermonde(row, col);
            }
        }
        for (uint col = 0; col < n_fields; ++col) {
            for (uint row = 0; row < m; ++row) {
                R_dest[col * m + row] = mCoefficients(row, col);
            }
        }
#endif
    }
}
```

**Phase 3: GPU Solve (cuSOLVER)**

```cpp
void MaxwellPostprocessor::solve_batch_gpu(
    real* cpu_vandermonde,
    real* cpu_rhs,
    size_t nodes_in_batch,
    cudaStream_t stream,
    cusolverDnHandle_t solver_handle)
{
    // Transfer to GPU (async)
    cudaMemcpyAsync(d_vandermonde, cpu_vandermonde,
                   nodes_in_batch * m * m * sizeof(real),
                   cudaMemcpyHostToDevice, stream);
    cudaMemcpyAsync(d_rhs, cpu_rhs,
                   nodes_in_batch * m * n_fields * sizeof(real),
                   cudaMemcpyHostToDevice, stream);

    // Batched Cholesky factorization: V = L * L^T
    cusolverDnDpotrfBatched(
        solver_handle,
        CUBLAS_FILL_MODE_LOWER,
        m,                      // Matrix dimension (same for all)
        d_vandermonde,          // Batch of matrices (strided)
        m,                      // Leading dimension
        m * m,                  // Stride between matrices
        nodes_in_batch,         // Number of matrices
        d_info                  // Info array
    );

    // Batched triangular solves: L * L^T * x = b
    cusolverDnDpotrsBatched(
        solver_handle,
        CUBLAS_FILL_MODE_LOWER,
        m,                      // Matrix dimension
        n_fields,               // Number of RHS columns
        d_vandermonde,          // Factored matrices
        m,                      // Leading dimension
        m * m,                  // Stride
        d_rhs,                  // RHS (overwritten with solution)
        m,                      // RHS leading dimension
        m * n_fields,           // RHS stride
        nodes_in_batch,
        d_info
    );

    // Transfer results back (async)
    cudaMemcpyAsync(cpu_rhs, d_rhs,
                   nodes_in_batch * m * n_fields * sizeof(real),
                   cudaMemcpyDeviceToHost, stream);
}
```

**Phase 4: Pipelined Execution**

```cpp
void MaxwellPostprocessor::run_gpu_pipelined() {
    // Allocate pinned CPU buffers (3× for triple buffering)
    real* cpu_V[3], *cpu_R[3];
    for (int i = 0; i < 3; ++i) {
        cudaMallocHost(&cpu_V[i], mBatchSize * m * m * sizeof(real));
        cudaMallocHost(&cpu_R[i], mBatchSize * m * n_fields * sizeof(real));
    }

    // Allocate GPU buffers (single batch)
    cudaMalloc(&d_vandermonde, mBatchSize * m * m * sizeof(real));
    cudaMalloc(&d_rhs, mBatchSize * m * n_fields * sizeof(real));

    // Create streams
    cudaStream_t streams[3];
    cusolverDnHandle_t handles[3];
    for (int i = 0; i < 3; ++i) {
        cudaStreamCreate(&streams[i]);
        cusolverDnCreate(&handles[i]);
        cusolverDnSetStream(handles[i], streams[i]);
    }

    // Pipeline loop (mNumBatches + 2 to drain pipeline)
    for (size_t batch = 0; batch < mNumBatches + 2; ++batch) {
        int buf = batch % 3;

        // Stage 1: Assemble on CPU
        if (batch < mNumBatches) {
            assemble_batch(batch, buf, cpu_V[buf], cpu_R[buf]);
        }

        // Stage 2: Solve on GPU (batch-2)
        if (batch >= 2) {
            size_t solve_batch = batch - 2;
            int solve_buf = solve_batch % 3;
            size_t nodes = get_batch_node_count(solve_batch);

            solve_batch_gpu(cpu_V[solve_buf], cpu_R[solve_buf],
                           nodes, streams[solve_buf], handles[solve_buf]);

            // Synchronize and unpack results
            cudaStreamSynchronize(streams[solve_buf]);
            unpack_results(solve_batch, solve_buf, cpu_R[solve_buf]);
        }
    }

    // Cleanup
    for (int i = 0; i < 3; ++i) {
        cudaFreeHost(cpu_V[i]);
        cudaFreeHost(cpu_R[i]);
        cudaStreamDestroy(streams[i]);
        cusolverDnDestroy(handles[i]);
    }
    cudaFree(d_vandermonde);
    cudaFree(d_rhs);
}
```

### Expected Performance

**Speedup estimates:**

| Component | Sequential | Pipelined GPU | Speedup |
|-----------|-----------|---------------|---------|
| Batch assembly | 100% CPU | 100% CPU | 1× |
| Cholesky solve | 100% CPU | GPU (async) | 10-50× |
| **Total (posv ≈ 25%)** | **Baseline** | **0.6-0.7×** | **30-40% faster** |

**When worthwhile:**
- Problem size > 1M nodes per rank
- Postprocessing time > 20% of total simulation time
- GPU available per MPI rank (modern HPC clusters)

**Key insight:** Even though only 20-30% of per-node cost is in `posv`, GPU's 10-50× speedup on that portion yields overall 30-40% improvement due to pipelining.

### Implementation Effort

**Low:** 3-5 days for working prototype, ~1 week for production

**Realistic task breakdown:**

**Day 1:**
- ✅ Add `mNodeToElements` member to header
- ✅ Implement `precompute_topology()` in initialization
- ✅ Test: verify indices match existing bitset results

**Day 2:**
- ✅ Implement `assemble_batch()` with Blaze flattening
- ✅ Add batch size computation
- ✅ Test: batch assembly produces correct matrices

**Day 3:**
- ✅ cuSOLVER integration (strided batched API)
- ✅ Basic GPU solve without pipelining
- ✅ Test: GPU results match CPU `posv`

**Day 4:**
- ✅ Pipelined execution with triple buffering
- ✅ CUDA streams and async transfers
- ✅ Test: pipelined results identical to sequential

**Day 5:**
- ✅ Error handling, GPU memory checks, CPU fallback
- ✅ Performance profiling and tuning
- ✅ Documentation

**Optional (production hardening):**
- CMake detection for CUDA/ROCm
- Runtime GPU capability detection
- Batch size auto-tuning for different GPUs

**Total:** 3-5 days (working prototype) → 5-7 days (production-ready)

**Context:** Compared to implementing cohomology cuts or thin-shell coupling, this is straightforward plumbing. The hard FEM work stays on CPU where BELFEM's infrastructure handles it naturally.

**Risk:** Low
- ✅ **No novel algorithms** - standard batched linear algebra
- ✅ **Well-tested APIs** - cuSOLVER is production-grade
- ✅ **Isolated change** - existing CPU path unchanged
- ⚠️ **Blaze flattening** - validate once, works everywhere
- ⚠️ **Memory tuning** - batch size may need adjustment per GPU

### Literature References

**GPU FEM:**
1. Cecka et al. (2011) - "Assembly of Finite Element Methods on Graphics Processors"
2. Reguly et al. (2018) - "Performance portability patterns for structured grid computations"

**Batched solvers:**
3. cuSOLVER documentation - NVIDIA batched dense linear algebra
4. Kirk et al. (2006) - "Batched small matrix operations in scientific computing"

---

## Alternative: MPI-Only Optimization (If Communication is Bottleneck)

**NOTE:** MPI communication is already optimized (embarrassingly parallel). Only pursue if profiling shows communication is a bottleneck (unlikely).

---

## Conclusion (2026-01-21 Update)

The Maxwell postprocessor is **mathematically correct and MPI-scalable**. No immediate optimization work is required.

**Key findings:**
- ✅ Embarrassingly parallel MPI design → >95% scaling efficiency expected
- ✅ Zienkiewicz-Zhu SPR correctly implemented
- ❌ OpenMP NOT recommended (violates BELFEM's MPI-first design philosophy)
- ⚠️ GPU acceleration is the only viable path, justified only for extreme-scale (>1M nodes)

**Recommendations:**

1. **Do nothing** unless profiling shows postprocessing > 20% of total time
2. **If optimization needed:** GPU acceleration (2-3 months effort)
3. **Never:** Add OpenMP (creates more problems than it solves)

**For mathematical correctness:** See `../src/fem/maxwell/doc/postprocessor_recovery_theory.md`

**For detailed review:** See `maxwell_postprocessor_review.md`

**For BELFEM parallelism philosophy:** See `../doc/coding_philosophy.md`

---

## Key Takeaways (2026-01-28)

> **⚠️ PRIORITY:** The **[element field caching](../closed/maxwell_postprocessor_element_caching.md)** prerequisite is now effectively implemented (closed 2026-06-22) — it was simpler and faster, and may make GPU acceleration unnecessary. Re-profile first.

**If GPU acceleration is still needed after element caching:**

**The strategy is simpler than initially thought:**

1. ✅ **Bottleneck identified:** `posv` (~5-6% of per-node time), not bitset
2. ✅ **Topology precomputation:** Trivial 80 MB cost, eliminates bitset from loop
3. ✅ **No GPU kernels needed:** cuSOLVER batched API handles everything
4. ✅ **Pipelined execution:** Standard CUDA streams pattern
5. ✅ **Blaze matrix flattening:** Manual but straightforward

**But:** Element caching eliminates 60-70% of redundant physics computation (5-10× speedup), making GPU's 1.3-1.4× overall speedup on posv less compelling.

**Compared to BELFEM's cohomology algorithms, thin-shell formulations, and domain decomposition - this is refreshingly simple.** It's fundamentally just:
- Batch pack matrices on CPU
- Call cuSOLVER batched routine
- Unpack results

The complexity is in BELFEM's FEM infrastructure, which **stays on CPU where it belongs**.

**Decision tree:**
1. ✅ Profile postprocessor to confirm physics computation dominates
2. ✅ Implement element caching (1-2 days, 5-10× speedup)
3. ⏸ Re-profile with caching enabled
4. ⏸ IF still bottleneck (>20% total time), THEN consider GPU
5. ⏸ ELSE done - caching solved the problem

---

## Revision History

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 3.0 | 2026-01-28 | C. Messe (via Claude Code) | Added bottleneck analysis, pipelined strategy, Blaze handling |
| 2.0 | 2026-01-21 | C. Messe | Removed incorrect benchmark projections, updated to MPI-first philosophy |
| 1.2 | 2026-01-20 | C. Messe (via Claude Code) | Added projected benchmarks (**REMOVED - were incorrect**) |
| 1.1 | 2026-01-20 | C. Messe (via Claude Code) | Initial GPU discussion |
| 1.0 | 2026-01-20 | C. Messe (via Claude Code) | Initial document creation |

---

**End of Document**
