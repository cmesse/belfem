# Assembly Profiling — Sparse Scatter and Field-Lookup Costs

**Date:** 2026-07-17
**Topic:** Profile of hphirun system-matrix assembly (follow-up to
`dl20260717_bfm_load_profile.md`); A/J double-search fix
**AIs involved:** Claude (analysis), Christian (fix)
**Literature References:** Messe et al. 2023 (paper1), §2.7 (hybrid Picard/Newton)

## Profile facts (1,328 samples, debug build `-Og`, no `-DNDEBUG`)

| Block | Share |
|---|---|
| `IWG_Maxwell::compute_mkf` (h_ts_metal 13%, h_ghost 10%, h_ts_hts 4%, Calculator::link 6%) | ~34% |
| Sparse scatter (`assemble_jacobian` + `assemble_newton` + SpMatrix binary searches) | ~23% |
| `Calculator::qold` (string-keyed field lookups per dof per element) | ~8.6% |
| `Tmatrix::project` (small dgemms) | ~7.5% |
| `Matrix::operator` accessor flat cost (debug-inflated) | ~10% |

`std::lower_bound` alone was 7.5% flat — every CSR write ran
`SpMatrix::operator()` → `index_csr_zero_based` → binary search over the row.

## Key insight: why A (system) and J (jacobian) both exist

They share the sparsity graph but NOT the values — this is the hybrid
Picard/Newton strategy (Messe et al. 2023, paper1, §2.7):
`assemble_jacobian` writes the base operator K into both; `assemble_newton`
adds the consistent-linearization terms into J only
(`cl_FEM_DofMgr_SolverData.cpp:1694`). Newton branch: residual `r = A·x − b`
(line 2146) but solve with J (line 2152). Picard branch: solve with A, J unused.

## Change made

`SolverData::assemble_jacobian` free-free branch
(`cl_FEM_DofMgr_SolverData.cpp:1629-1643`, by Christian): since A and J share
the pattern, compute the CSR position once via the existing `SpMatrix::index()`
and write both raw `data()` arrays — halves the dominant scatter searches.
Claude added the debug bounds assert (`tPos < nnz`) that `operator()` used to
provide — without it, an out-of-pattern (i,j) would silently write
`data()[nnz]` out of bounds. Syntax-checked with production flags.

A considered alternative (assemble A only + one bulk `J += A` axpy after the
element loop, which would also skip J entirely during pure Picard) was set
aside in favor of the minimal in-place fix.

## qold acceleration — Christian's pointer-table design, audited

Christian implemented the `qold` hoist himself: `mQold` table of
`Vector<real>*` keyed `s * mMaxDofFieldIndex + fieldIndex`, built in
`init_qold_table()` (called from `Calculator::link(Group*)`), plus
`IWG::timestepping_order()` (base returns 0, `IWG_Timestep` returns mOrder).
Audited by Claude + Codex (thread `tmp/ai_exchange/qold_field_cache.md`;
Grok unavailable — 3 narration-stub runs, standing CLI issue).

**Verdict: design sound, two blockers before it can run.**

- [x] **F0 (Codex, high):** the hook never fires for BLOCKS — block assembly
      goes through `IWG::link_to_group()` (sets `mCalc` only); only the SideSet
      ctor calls `Calculator::link(Group*)`. Block calculators get an empty
      table → loud Map miss on first `qold()`. Move the build to a lifecycle
      point covering all calculators after the timestep method is set, and
      rebuild if `set_timestepping_method()` reconfigures.
- [x] **F1 (Claude, Codex-confirmed, high):** stride collision —
      `mMaxDofFieldIndex` is the INCLUSIVE max, so key(s, max) == key(s+1, 0);
      silent wrong-history at order ≥ 2. Stride must be max+1; better, replace
      the `Map` with a pre-sized `Cell<Vector<real>*>` (also removes the
      per-dof integer hash).

**Fixed 2026-07-17 (Claude, Christian-approved):** F1 → `mQold` is a
`Cell<Vector<real>*>` with stride `mMaxDofFieldIndex + 1` (aliasing comment in
code); F0 → `init_qold_table()` moved from `link(Group*)` to `allocate()`
(after the `mIsAllocated` guard — reached by block AND sideset calculators via
`DofManager::init_work()`), with early-outs for kernel-less groups and
`timestepping_order() == 0`. `qold()` now direct-indexes the Cell with a
nullptr debug assert. Syntax-checked. **Codex fix-verification: CONFIRMED
(high confidence)** — allocate() coverage proven for concrete blocks and
sidesets (thin shells are DomainType::ThinShell BLOCKS, same path; dummy/empty
groups never assemble), no in-tree qold() precedes init_work(), stride
arithmetic verified. Two caveats accepted: (1) a set_timestepping_method()
call AFTER calculator allocation would not rebuild the cache — all in-tree
callers configure before initialization; (2) missing-entry diagnostics are
debug-only (deliberate: internal invariant on a per-dof hot loop =
BELFEM_ASSERT tier).

**Cleared:** F2 ramp-up (standard Controller paths call `shift_fields()` before
any read — `cl_FEM_Controller.cpp:134-145,248-258,332-349`; the old lazy-copy
branch was never the normal warm-start), F3 (`mOrderActive ≤ mOrder` holds),
F5 (standard old-field creation already set write-flag false; memdump/restart
serializes all fields regardless), F6 pointer stability (heap `Field` objects,
`Cell<Field*>` — cached pointers cannot dangle).

Remaining open:
- [ ] Release-profile recheck: `Matrix::operator` accessor cost and OpenBLAS
      small-gemm buffer overhead (`blas_memory_alloc`, ~2.5%) need a `-O3
      -DNDEBUG` profile before acting.
- [ ] `compute_mkf` kernels (34%) = real physics; any tuning there is a
      separate literature-guided task, not a mechanical pass.
