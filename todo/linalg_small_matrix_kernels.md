# Small Matrices in the Calculator (O3): Fixed Sizes Are Not the Win, Armadillo's BLAS Dispatch Is

**Date:** 2026-11-07
**Purpose:** Settle O3 of the linalg wrapper investigation — would fixed-size small matrices
(`arma::Mat<T>::fixed<3,3>`, `blaze::StaticMatrix<T,3,3>`) pay in BELFEM's calculator? Measured on
the per-integration-point Jacobian chain (`J = dNdXi·X`, `inv3`, `dNdX = invJ·dNdXi`, hex8) on both
backends: on Blaze a fixed `J`/`invJ` buys 12 % of the chain and a fully static chain 30 %; on
Armadillo fixed sizes buy **nothing**, because the non-square 3×8 and 8×3 products are dispatched to
OpenBLAS `dgemm` per integration point, which makes the Armadillo chain **3.1× slower than Blaze**
(thread count pinned; 5.3× with OpenBLAS's pool unpinned). The actionable item is therefore not a
second matrix type but an **Armadillo-side** small-product kernel behind the existing backend split
(a hand loop loses to Blaze's own kernel at hex27), gated by a measurement of the chain's share in a
real assembly.
**Module:** `src/linalg` (new `fn_*` kernels, Armadillo implementation + Blaze forward),
`src/fem/kernel/cl_FEM_Calculator.hpp` / `.cpp` (nine product sites, two scratch sizes)
**AIs involved:** Claude (probe, draft); Codex `gpt-6-astra` + Grok `grok-4.7` (plan jury)
**Status:** PLAN v2 — three probes done 2026-11-07; **plan jury round 1 done** (Codex `gpt-6-astra`/high: 5 findings, 1 P1 —
`mdN` scratch sizing; Grok `grok-4.7`/high: 12 challenged claims + open risks, all CONFIRMED except one REFUTED). **Awaiting
Christian's decision**; the next step is **R1**, his profile. No source modified. Successor of O3 in
`todo/linalg_wrapper_design_investigation.md`; sequenced after round A of `todo/linalg_wrapper_hardening_plan.md`, independent
of rounds B–D (K wins before round B; round B adds to it).

> **Scope guards:**
> - Both backends; the two defaults (Linux Armadillo, Apple Blaze) must agree in results **to a stated
>   tolerance** (relative 1e-13 on `J`, `invJ`, `B`): a hand loop, `dgemm` and Blaze's kernel round
>   differently, and today's two backends already differ at that level.
> - The MORIS idea under test was `Matrix< F33RMat >` (Eigen fixed 3×3; MORIS's **Armadillo** branch
>   typedefs `F33RMat` to a *dynamic* `arma::Mat<real>`, `linalg_typedefs.hpp:69` — MORIS never had a
>   fixed-size Armadillo matrix either). 5 MORIS files use `F33RMat` (15 counting `F31RVec`/`F31RMat`),
>   mostly geometry rotations.
> - Out of scope: changing the default backend (O1 is decided: both stay), `ARMA_DONT_USE_BLAS`
>   (large products need BLAS), a second wrapper type for the whole tree.

---

## 1. What was measured

Probe `o3/bench.cpp` (source in `tmp/ai_exchange/review_linalg_small_matrices.md`): the hex8 chain per
integration point, 2·10⁶ repetitions, best of 5, GCC 15 `-O2 -DNDEBUG`, the tree's own defines and
headers (round A applied, i.e. the wrapper is already non-virtual). `X(0,0)` is perturbed every
iteration. **This first table was run with OpenBLAS's thread pool unpinned and observes only `detJ`
and `dNdX(0,0)`** — §1.2 repeats it pinned and fully observed; the ratios survive, the Armadillo
absolute times do not. Variants: **W** `belfem::Matrix<real>` (today's calculator), **D** raw
backend dynamic, **F** backend fixed `J`/`invJ` with dynamic `dNdXi`/`X`/`dNdX` (what a `SmallMatrix`
member in `calculator::MatrixData` would give), **A** everything fixed (needs the element type at
compile time — an upper bound, not a design).

| backend / configuration | W wrapper | D dynamic | F fixed J, invJ | A all fixed | W/D | F/D | A/D |
|---|---|---|---|---|---|---|---|
| Blaze (`BELFEM_NETLIB` → `BLAZE_BLAS_MODE 0`) | 161 ms | 157 ms | 138 ms | 97 ms | 1.03 | **0.88** | **0.62** |
| Armadillo, checks on (today's release) | 826 ms | 842 ms | 867 ms | 867 ms | 0.98 | 1.03 | 1.03 |
| Armadillo, `ARMA_DONT_CHECK_CONFORMANCE` (round-B policy) | 800 ms | 770 ms | 752 ms | 743 ms | 1.04 | 0.98 | 0.96 |
| Armadillo, `ARMA_DONT_USE_BLAS` (**probe only**, not a config option for BELFEM) | 377 ms | 377 ms | 403 ms | 279 ms | 1.00 | 1.07 | 0.74 |

`sizeof`: Blaze `StaticMatrix<double,3,3,columnMajor>` 96 B (padded columns, `BLAZE_USE_PADDING 1`,
`config/Optimizations.h:62`), `DynamicMatrix` 40 B; Armadillo `fixed<3,3>` 208 B vs `Mat` 192 B.

**Findings (confidence):**

1. **The wrapper is not in the picture** — W/D within ±4 % on every row (high; round A's result on a
   second kernel).
2. **Blaze: fixed sizes help modestly.** Fixed `J`/`invJ` alone: 12 % of the chain; the chain fully
   static: 38 % (high on the numbers, medium on what they mean for a real assembly — the chain is one
   part of per-element work, share unmeasured).
3. **Armadillo: fixed sizes do nothing, and the chain is 3.1× slower than Blaze** (pinned, §1.2; 5.3×
   unpinned). Cause, source-traced: assignment of a product goes `glue_times::apply`
   (`glue_times_meat.hpp:557-565`) → `glue_times_redirect` → `glue_times::apply` (`:737-761`), which
   takes `gemv` for a 1-row or 1-column operand (`:759-760`) and otherwise `gemm<>::apply`, whose
   `double` overload is `apply_blas_type` (`mul_gemm.hpp:534-547`);
   the inline emulation is taken only for **square** matrices of order ≤ 4 with a square right operand
   (`mul_gemm.hpp:413-417`); everything else goes to BLAS (`:454-460`, `ARMA_USE_BLAS` is set in
   `armadillo.hpp:20`). `dNdXi·X` (3×8 · 8×3) and `invJ·dNdXi` (3×3 · 3×8) are not square, so both are
   `dgemm` calls into OpenBLAS per integration point. `fixed` changes storage, not dispatch. With the
   BLAS path removed (probe only) the chain drops to 377 ms — 2.2× faster — and is still 2.4× Blaze
   (high).
4. **Where the Armadillo time goes is the finding O3 was looking for**, and it is not the finding O3
   expected. Blaze reaches its number with its own kernels: under `BELFEM_NETLIB` it has no BLAS at all
   (`BLAZE_BLAS_MODE 0`, `blaze_config.hpp:61-62`, which every non-MKL build is, `config_mkl.cmake:66-67`),
   and under MKL the column-major threshold `BLAZE_TDMATTDMATMULT_THRESHOLD 4900` target elements
   (`Thresholds.h:292-293`, tested at `TDMatTDMatMultExpr.h:563`) keeps products this small off BLAS. A
   hand-written loop matches that on **Armadillo**; on Blaze it does not beat Blaze's kernel (§1.1).

### 1.1 Second probe — option K against the backend product (same flags, 2·10⁶ reps, best of 5)

`bench_k.cpp`: the same chain, wrapper types throughout, with the two products either as the backend
expression (`J = dN*X; dNdX = iJ*dN`) or as backend-free `(i,j)` loops (K) with `dim` a template
parameter and `n` runtime, exactly the pattern `inv2`/`inv3` use.

| element | Blaze: backend / K | Armadillo checks on: backend / K | Armadillo checks off (round B): backend / K |
|---|---|---|---|
| quad4 (2-D, 2×4·4×2) | 70 / **33** ms (0.48) | 604 / 120 ms (0.20) | 522 / **65** ms (**0.13**) |
| hex8 (3-D, 3×8·8×3) | 151 / 156 ms (1.03) | 693 / 275 ms (0.40) | 679 / **205** ms (**0.30**) |
| hex27 (3-D, 3×27·27×3) | 320 / 480 ms (**1.50**) | 1088 / 850 ms (0.78) | 989 / 607 ms (0.61) |

5. **K pays on Armadillo for every element measured (1.6×–8.0× unpinned; 1.4×–6.9× pinned, §1.2), and
   is a wash or a loss on Blaze** (par at hex8, 1.4–2.1× gain at quad4, 1.3–1.5× loss at hex27, where
   Blaze's vectorized kernel beats a scalar loop). So K is an **Armadillo-side** kernel behind the existing backend split (`src/linalg/armadillo/`
   vs `src/linalg/blaze/`), with the Blaze side forwarding to the backend product — not a generic loop.
   (high on the numbers; the Blaze quad4 gain is real but small in absolute terms.)
6. **The round-B check policy matters for K too**: with Armadillo's bounds checks on, the `(i,j)` loop
   pays a check per element read (hex8 275 → 205 ms, quad4 120 → 65 ms once the checks are off). K is
   sequenced after round B for that reason, or ships with the loops reading `memptr()`-based column
   pointers — a design point for the jury.
7. **2-D is the stronger case** on Armadillo (0.14 pinned): less arithmetic per `dgemm` call, same call
   overhead. The pre-registration's expectation (and the former O1) is confirmed.

### 1.2 Third probe — pinned threads, every output observed (Codex round-1 finding 2)

`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`; two `X` entries in two columns perturbed per iteration; every
element of `J` and `dNdX` summed into the accumulator. Same flags, 2·10⁶ reps, best of 5.

| | Blaze | Armadillo, checks off |
|---|---|---|
| W / D / F / A (hex8) | 178 / 181 / 159 / 127 ms → F/D **0.88**, A/D **0.70** | 570 / 567 / 556 / 556 ms → F/D 0.98, A/D 0.98 |
| K vs backend, quad4 / hex8 / hex27 | 48/69 = **0.70**, 186/184 = 1.01, 583/434 = **1.34** | 60/415 = **0.14**, 219/571 = **0.38**, 728/1011 = **0.72** |
| Blaze D vs Armadillo D (hex8) | 181 ms | 567 ms → **3.1×** |

The first probe's Armadillo times (826/842 ms) were inflated by the unpinned OpenBLAS pool — the
thread-pinning risk both auditors raised, confirmed. All ratios of §1 and §1.1 survive within a few
percent; v2's headline numbers are these.

## 2. Options

| option | Blaze gain on the chain | Armadillo gain on the chain | radius | verdict |
|---|---|---|---|---|
| **K — small-product kernels behind the backend split**: Armadillo side as explicit `(i,j)` loops like `inv3`/`inv2`; Blaze side forwards to the backend product (its kernel wins at hex27) | 0 by construction (Blaze keeps its product); optional: K for 2-D only on Blaze (0.48) | **2.6× hex8, 6.9× quad4, 1.4× hex27** on the chain, pinned, checks off (§1.2); more with checks on today's release (2.5× at hex8, §1.1) | `src/linalg/{armadillo,blaze}/fn_*` + the dispatch header + 3–4 call sites in `cl_FEM_Calculator.hpp` | **recommended, after measuring the chain's share (R1)** |
| **S — `SmallMatrix<T,R,C>` type** (design-C style second template) for `J`/`invJ` in `calculator::MatrixData` | 12 % of the chain | 0 | new class, `MatrixData` becomes a template or variant; direct `J()`/`Jm()` consumers are the five `const Matrix<real>& J = this->Jm(…)` lines in `cl_FEM_Calculator.cpp` (`:2304`, `:2418`, `:2533`, `:2591`, `:2698`), `cl_Maxwell_TMatrix.cpp:119` and five `det( this->J(…) )` calls (`.hpp:2116-2256`); the 52 IWG lines call `B()`/`dV()`, which S would not retype. There is **no** zero-copy conversion from a fixed type to `const Matrix&` (the wrapper owns a dynamic backend member), so S means a second API, not an adapter | not now: K first; S only if R1 shows the chain dominates on Blaze |
| **A — static element types** (templating the calculator on `nnodes`) | 38 % of the chain | with K: same as K | the calculator's whole design | no — a rewrite for an unmeasured share |
| **B — `ARMA_DONT_USE_BLAS`** | — | 2.2× on the chain | kills BLAS for every large product in Armadillo | no |

## 3. Ordered steps

- [ ] **R1 — measure the chain's share (Christian).** One 2-D and one 3-D assembly-dominated deck per
  backend (the quad4/hex27 spread is 6.9× vs 1.4×, one deck cannot see it), release, `OMP_NUM_THREADS=1`
  **and** `OPENBLAS_NUM_THREADS=1` (`parallel_execution.md:57-66`: SCLS OpenBLAS follows `OMP_NUM_THREADS`,
  pthread builds the other), BLAS implementation, revision, hardware and full commands recorded. The
  accessors are `inline` (`cl_FEM_Calculator.hpp:1718`, `:1886`, `:1939`, `:1967`) and nest
  (`Bscalar → invJ → J`), so a flat profile mis-attributes: use `perf record --call-graph dwarf` and take
  the **union** of samples attributed to the chain (never a sum of inclusive percentages), or a
  `Timer`-instrumented probe build around `J()`/`invJ()`/`B()`; under `ARMA_USE_WRAPPER` the BLAS entry
  is `arma_dgemm` (`def_blas.hpp:54`), not `dgemm_`. Decision rule: Armadillo chain share < 5 % of
  assembly on both decks → K is hygiene, behind round D; otherwise K is next after round A.
- [ ] **R2.0 — setup sizes first (Codex round 1, P1).** `mdN` is created 2×(2n) for plane stress and
  3×(3n) for elasticity (`cl_FEM_Calculator.cpp:1280`, `:1291`) but receives the dim×n product
  `invJ * dNdXi` (`cl_FEM_Calculator.hpp:1948`, `:1977`); today the backend silently resizes it on the first
  assignment (`glue_times_meat.hpp:753`, one reallocation). A size-asserting kernel would fail there.
  Inventory every product site and correct the `create_matrix` sizes (`coding_philosophy.md:613`:
  scratch sized at setup), with a test per calculator path (scalar, plane stress, elasticity,
  master/slave).
- [ ] **R2 — kernels (after R1 and R2.0; plan jury done, code jury owed).** Dispatch header
  `src/linalg/fn_small_times.hpp` (name to be settled) declaring `jacobian( aDNdXi, aX, aJ )`
  (dim×n · n×dim → dim×dim) and `gradient( aInvJ, aDNdXi, aDNdX )` (dim×dim · dim×n → dim×n) — `a`
  prefixes as `fn_inv3.hpp:47`. Armadillo implementation `armadillo/fn_AR_small_times.hpp`: `(i,j)`
  loops with `dim` a template parameter (§1.1); Blaze implementation `blaze/fn_BZ_small_times.hpp`
  forwarding to the backend product through `matrix_data()` (the `fn_BZ_dot.hpp:27-29` pattern). Contract:
  `BELFEM_ASSERT` on all four sizes, output must not alias an input (asserted in debug), `dim ∈ {2,3}`.
  **Dispatch:** not through a function pointer like `mFunInvertJ` (`:734`, `:1784`) — the probe's gain is
  the inlined loop (Grok round 1); branch on the stored dimension inside the inline accessor
  (`if( mNumberOfDimensions == 3 ) jacobian<3>(…) else jacobian<2>(…)`), measured in the gate.
  **Sites:** `cl_FEM_Calculator.hpp:1724-1726` (curved/linear), `:1740-1741`, `:1755-1756` (J, Jm, Js);
  `:1894-1895` (Bscalar), `:1915` (Bm), `:1931` (Bs), `:1948` (Bplanestress), `:1977` (Bvoigt); outside the
  calculator `fn_Mesh_compute_volume.cpp:119`, `cl_Pipette.cpp:312`, `:366`, `:487` (same product, lower
  call rate — second pass). `N * X` (`:2598`, `:2606`) is the `gemv` arm and is not touched.
  **Gate and acceptance:** numerical tests against the backend product on both backends (both dimensions,
  curved and linear, every element type in the tree) to the scope-guard tolerance; the §1.2 probe re-run
  by hand on both backends with the kernels as landed — **accept only if** Armadillo is faster on every
  measured element and Blaze is unchanged (it forwards); R1's profile repeated on the same decks;
  `make check` both backends. A checked-in benchmark that does not assert on time is documentation,
  not a gate, and is not claimed as one.
- [ ] **R3 — `SmallMatrix` for `J`/`invJ` (conditional on R1 showing Blaze's chain share matters).**
  Design C in the report: a second template with the same API subset. There is no zero-copy conversion
  to `const Matrix&` (the wrapper owns a dynamic backend member), so S is a second API with its own
  storage, lifetime and copy costs to specify — not an adapter. Decision deferred to after R2's numbers.
- [ ] **R4 — close out:** O3 in the investigation report marked RESOLVED with the verdict; devlog;
  debt row DR-200 widened (same subsystem, policy 2026-08-29), not a new row.

## 4. Open questions

- **O1 — 2-D.** RESOLVED 2026-11-07 by probe (§1.1, §1.2): quad4 is the stronger case on Armadillo (0.14).
  Unmeasured: n between 8 and 27 on Blaze (hex20, tet10, quad8/9, wedge) — Blaze forwards
  unconditionally, which cannot regress.
- **O2 — Blaze's `BLAZE_BLAS_MODE`.** Under `BELFEM_NETLIB` (every non-MKL build) it is 0: Blaze never
  calls BLAS, even for large products. Under MKL it is 1 with the 4900-target-element threshold. The
  `BELFEM_ACCELLERATE` branch (`blaze_config.hpp:63`) is set by no CMake file — dead. Not an O3 item.
- **O3 — do the kernels belong in `src/linalg` or in the calculator?** `inv2`/`inv3` live in `src/linalg`
  (`fn_inv2.hpp`, `fn_inv3.hpp`); the products are the same kind of thing. Default: `src/linalg`.
- **O4 — K on Blaze for 2-D only?** 0.48 at quad4, 1.03 at hex8, 1.50 at hex27. Default: no — one
  backend-specific kernel is enough to reason about; revisit if R1 shows 2-D Blaze assemblies dominated
  by the chain.
- **O5 — element access in the Armadillo kernel.** `(i,j)` only (`coding_philosophy.md:324`: never
  compute offsets from `data()`), as `inv2`/`inv3` do. Armadillo's own bounds check costs ~25 % of the
  loop until round B lands (275 vs 205 ms at hex8); K still wins 2.5× before that, so K is **not** gated
  on round B — round B adds to it. RESOLVED 2026-11-07 (plan jury).

## 5. Audit trail

- Probes and numbers: `tmp/ai_exchange/review_linalg_small_matrices.md` (pre-registration carries
  `bench.cpp`; the addendum `bench_k.cpp`; the round-1 verification the pinned, fully observed re-run).
  The timed inverse is a local stand-in, not production `inv3` (`fn_inv3.hpp:77`); ratios isolate the
  products because it sits on both sides.
- Plan jury round 1, 2026-11-07 18:44–18:56 CST: Codex 5 findings (P1: `mdN` sizing), Grok 12 challenged
  claims + open risks; all CONFIRMED except one REFUTED (the MORIS typedef file exists). Verification and
  reconciliation in the thread. Plan → v2. Compile line: the tree's flags with `-DBELFEM_BLAZE -I src/linalg/blaze` or
  `-DBELFEM_ARMADILLO -I src/linalg/armadillo -larmadillo`, `-lopenblas`, `-O2 -DNDEBUG`.
- Plan jury: Codex `gpt-6-astra`/high + Grok `grok-4.7`/high (critical-framework row: the calculator).
