# linalg Wrapper Design: Alias the Backend Type, or Keep the Member Wrapper?

**Date:** 2026-11-07
**Purpose:** Answer three questions about `src/linalg`: (a) is a `Matrix` that aliases the backend
object directly, with no wrapped `mMatrix` member and no `matrix_data()`, cheaper at runtime than
BELFEM's current wrapper; (b) can BELFEM adopt that design, and is it worth it; (c) what is the
impact radius, and how would an implementation plan look. The MORIS linalg (`tmp/moris/projects/LINALG`)
was named as the model that had "developed further" in this direction.
**Module:** `src/linalg` (`armadillo/cl_AR_Matrix.hpp`, `blaze/cl_BZ_Matrix.hpp`,
`armadillo/cl_AR_Vector.hpp`, `blaze/cl_BZ_Vector.hpp`, `operators/`, `fn_*.hpp`, `lapack/`);
consumers in every module
**AIs involved:** Claude Fable 5.1 (investigation, probes, draft); Codex `gpt-6-astra`/high (round-1
jury, six findings, all confirmed and folded in); **Grok `grok-4.7`/high: seat FAILED** (the `grok` CLI
is not signed in on this machine and no `XAI_API_KEY` is set) — the reserve seat Copilot (`auto` →
`gpt-5.6-luna`, protocol §9.2) substituted; round 2 at `xhigh` on this revised document; **Grok `grok-4.7`/xhigh, round 3** (seat restored 2026-11-07): ten findings, all CONFIRMED against source — one premise (that Apple builds use Clang) refuted for this machine — folded into v6
**Status:** v6, final for discussion — revised 2026-11-07 after round 3 (Grok `xhigh`, the first non-OpenAI voice: the Blaze `Matrix` expression `+=`/`-=`/`*=` exist only under `BELFEM_GCC`, so the overload set is compiler-dependent (P1); the Armadillo `Vector` has the `Matrix`'s constructor hole; `Vector::operator() const` is not virtual; the backend's compound assignment allocates on alias; two §9 rows lowered) — the implementation plan is `todo/linalg_wrapper_hardening_plan.md`; after round 1 (Codex: the "only runtime cost" claim and the
Armadillo benchmark column refuted; Copilot: R1 is a downstream API break) and round 2 (Codex `xhigh`:
a third overload gap in construction, a regression in the first R2a recipe, static call-site counts are
not executed allocation counts, six count and scope corrections; Copilot: the `final` break needs an
explicit decision (O5), `linalg_usage_guide.md:98` names the wrong Armadillo macro, `-=` coverage and a
per-site inventory before R2b). Round 2 was owed under the stop condition; a round 3 is not — the next
step is the executable gate. Decision pending Christian. **No source modified; read-only session.**
Every number below was measured or counted in this session; the grep scope is given with each count.
**The Grok voice is missing from this round**; a Grok pass should run once the CLI is signed in
(`grok login --device-code`) before the design is settled.

> **Scope guards:**
> - No source edits. Artifacts: this file, `tmp/ai_exchange/review_linalg_wrapper_design.md`, a devlog.
> - The cohomology core is closed to AI edits (protocol §7.1). Its one escape-hatch site is counted,
>   not proposed for change.
> - `nonfree/` enters only as aggregate counts; no nonfree content is quoted here or in the exchange.
> - Literature-first is not required (protocol §5: no algorithm, formulation, or physics is at stake).
>   The coding shelf was consulted for the wrapper-design rationale; see §8.

---

## 0. Summary

1. **The premise does not hold.** MORIS has not dropped the container. Every MORIS branch (33 remotes
   in the local clone, HEAD 2026-10-12) and GitHub `main` (fetched 2026-11-07) still define
   `Matrix< arma::Mat<Type> >` with a private `arma::Mat<Type> mMatrix` and a public `matrix_data()`.
   The last commit touching MORIS LINALG is from 2024-08-29. What differs from BELFEM is the template
   parameter (the backend type, not the scalar) and the absence of a separate `Vector` class. §1.
2. **A container-less wrapper is not cheaper at runtime than a complete, non-virtual member wrapper.**
   A class with one data member and inline forwarding methods compiles to the member's own code. The
   wrapper costs something only where it adds what the member does not have, or lacks an overload the
   member has. §2.
3. **BELFEM's wrapper does add two things the member does not have.** First, two `virtual` members:
   `Matrix` and `Vector` carry a virtual destructor, and `Matrix` a virtual const `operator()`. Nothing
   derives from either class. Measured: +8 bytes per object under Blaze (+16 under Armadillo), a vtable
   check on every `Matrix::operator() const` call through a reference whose dynamic type is unknown
   (the `Vector` accessor is not virtual, `cl_AR_Vector.hpp:285`, `cl_BZ_Vector.hpp:340`), and a 2.6× (Blaze) / 1.3× (Armadillo) slowdown of a 3×3 kernel
   called through `const Matrix&`/`const Vector&`, which drops to parity on both backends when the
   `virtual`s are removed. §3. Second, **gaps in the overload sets relative to the backend's
   expression overloads**, demonstrated in the tested operations: under Armadillo `A += B * dt` and
   `A -= B * dt` on a `Matrix` materialize an `arma::Mat` temporary (heap above 16 elements, local
   storage below); under Blaze every expression `+=`/`-=` on a `Vector` (`y += x * dt`, `y += x + z`,
   `y += A.col(0)`, `y += A.row(0)`) materializes a temporary `Vector` (one `posix_memalign` + `free`,
   every size); and `Matrix C( B * dt )` allocates twice on both backends where the backend allocates
   once. The backend's own operators allocate only when the right-hand side aliases the target
   (`Mat_meat.hpp:5531-5533`, `DynamicVector.h:1412-1414`); none of the shapes above alias. §3b. The `cl_IWG_Timestep.cpp` escape hatches
   are workarounds for exactly the `+=` case.
4. **Pure aliasing is incompatible with the two-backend contract and the wrapper policy (L-21).**
   Inheritance works under Blaze and **does not compile under Armadillo 15.2**, whose expression
   templates are gated on exact-type traits; specializing those traits did not rescue it. §2, §4.
5. **Recommendation:** do not alias and do not inherit. Two steps are recommended, one conditional.
   (R1) Remove the six `virtual` specifiers and mark the four classes `final` (ten declaration lines in
   four headers, zero in-tree call-site changes, a deliberate API break for any out-of-tree class that
   derives from them — none known), gated by `make check` on both backends and a before/after kernel
   timing. (R2a) Audit the four classes' overload sets against the backends' expression overloads
   (constructor, `=`, `+=`, `-=`, `*=`) and close the gaps with *constrained* overloads, gated by
   **executed** allocation counters, not assembly inspection — the first recipe (restore the generic
   Blaze `+=`) breaks `y += A.row(0)` and is withdrawn. (R2b, conditional on R2a) Retire the 109
   escape-hatch lines in 12 editable files on the live paths, each with a named numerical test; the
   Crank–Nicolson/Galerkin lines sit behind a deliberate `BELFEM_ERROR` and stay. §5, §6.

---

## 1. Premise check: what MORIS actually does

Local clone: `tmp/moris`, HEAD `49c74333e` (2026-10-12), `origin = https://github.com/kkmaute/moris.git`.

| claim | evidence | confidence |
|---|---|---|
| The generic `Matrix` still wraps a member | `projects/LINALG/src/cl_Matrix.hpp:31-36`: `template< typename Matrix_Type > class Matrix { Matrix_Type mMatrix; ...`; `matrix_data()` at `:281-285` | high |
| The Armadillo specialization still wraps a member | `Arma_Impl/cl_Matrix_Arma_Dynamic.hpp:46-51`: `class Matrix< arma::Mat< Type > > { private: arma::Mat< Type > mMatrix;`; `matrix_data()` at `:527-540` returns `mMatrix` (both overloads) | high |
| The Eigen specialization too | `Eigen_Impl/cl_Matrix_Eigen_Dynamic.hpp:22-25` | high |
| No branch has a different design | `git grep -c` for the `mMatrix` member over all 32 `origin/*` branches (plus the symbolic `origin/HEAD`) in both files: every branch hits both files | high (enumeration printed in session) |
| GitHub `main` matches the clone | raw `cl_Matrix.hpp` and `Arma_Impl/cl_Matrix_Arma_Dynamic.hpp` fetched 2026-11-07: `arma::Mat< Type > mMatrix`, `matrix_data()` present, no base class | high |
| LINALG on `main` has not changed since 2024 | `git log -- projects/LINALG/src` on `main`: newest `ffa74e15d` 2024-08-29 (clang-tidy refactor). Across all refs (`git log --all`) there are later touches — `f52d91dd6` 2025-07-14 adds `fn_join_vert.hpp`, `8bd44d821` 2026-10-14 on `origin/cmake-refactor` edits `CMakeLists.txt` — none of them to the storage design | high |
| MORIS leans harder on the escape hatch than BELFEM | `matrix_data()` appears in 154 MORIS files, 92 of them inside LINALG and **62 outside**, versus 13 BELFEM files outside `src/linalg` (§5) | high |

Three things in MORIS can read as "aliasing" and are not:

- **The template parameter is the backend type.** `Matrix< DDRMat >` with `typedef arma::Mat< real > DDRMat`
  (`linalg_typedefs.hpp:56-69`), one partial specialization per backend, and the `fn_*` headers
  dispatch on `Matrix_Type` through `linalg_internal::` (`fn_trans.hpp:41-49`, `Arma_Impl/fn_trans_Arma.hpp`).
  Every specialization still owns an `mMatrix`. The parameter choice is organizational, not a
  storage choice.
- **`moris::Cell` was renamed `moris::Vector`** on 2024-02-07 (`47eb1116d`). That `Vector` is the
  `std::vector` wrapper, BELFEM's `Cell`, not a linear-algebra vector.
- **MORIS has no linear-algebra `Vector` class**; vectors are n×1 `Matrix`. BELFEM's `Vector<T>` is
  `arma::Mat<T>` n×1 under Armadillo (`cl_AR_Vector.hpp:38`) and `blaze::DynamicVector` under Blaze
  (`cl_BZ_Vector.hpp:46`).

**Consequence.** The question is not "follow MORIS" but "would a container-less `Matrix` be better
for BELFEM". §2 evaluates the three ways to build one.

---

## 2. Three container-less designs and what each would cost

**Design A — pure alias.** `template< typename T > using Matrix = arma::Mat< T >;` (resp.
`blaze::DynamicMatrix< T, columnMajor >`). No code, no overhead, and no BELFEM API: `n_rows()`
becomes `.n_rows` under Armadillo and `.rows()` under Blaze; `set_size` becomes `set_size` resp.
`resize(..., false)`; `length()` becomes `.n_rows` resp. `.size()`; `data()` becomes `memptr()` resp.
`data()`; the `spacing()`/`capacity()` transfer contract of `cl_Matrix.hpp:16-27` has no home; the
`BELFEM_ASSERT` bounds messages vanish; and every one of the 552 consumer files (§5) names the vendor
type directly — that is what an alias changes — which L-21 forbids as policy (`scripts/check_wrapper_policy.py` has no Armadillo/Blaze
family, so this boundary is enforced by review, not by the script). Both backends are live (Blaze is
the Apple default, Armadillo the Linux default, `CLAUDE.md` "Build Configuration"), so A requires a
free-function adaptor layer to keep one source tree compiling on both, and that layer is a wrapper
again. **Verdict: not applicable while BELFEM keeps two backends.** Confidence: high.

**Design B — inheritance.** `template< typename T > class Matrix : public Backend< T >` with
`using Base::Base; using Base::operator=;` and the BELFEM accessors added on top. Probed this session
(`inherit_probe.cpp` in the session scratchpad; both backends from `/opt/scls/include`,
`-O2 -std=gnu++17`):

| backend | result | detail |
|---|---|---|
| Blaze 3.9 | **compiles and runs** | `C = A * B`, `D = trans(A) * B + C`, `+= -= *=`, `y = A * x`, `z = 2.0 * x + y`, `std::vector< Matrix >`, move; `sizeof` equals the base (40 / 24); non-polymorphic. Blaze's free functions take `DenseMatrix< MT, SO > &` and deduce `MT` through the derived-to-base conversion. Codex re-ran the probe in round 2 and confirmed both backends' results. |
| Armadillo 15.2 | **does not compile** | `trans( A )`: "no matching function"; `A * x`: "no match for operator*". Armadillo's expression templates are gated on `is_arma_type< T1 >`, built from exact-type specializations (`armadillo_bits/traits.hpp:96-116`: `is_Mat< Mat<eT> >`, `is_Mat< Col<eT> >`, …), which are false for a derived class. |
| Armadillo + `is_Mat`/`is_Col` specialized for the derived type | **still does not compile** | same two errors — because the gate is `is_arma_type2< const T1 >` (`traits.hpp:542-571`), so the unqualified specialization never reaches it. Codex's const-qualified specialization passed the gate and then failed on missing `Proxy<>` operations (Codex round 1, not re-run here). Going further means vendor-internal `Proxy<>`/trait surgery that every Armadillo release can break. |

The alternative that does compile under Armadillo is a forwarding layer that casts to the base
inside every free function and operator, i.e. today's design with `static_cast< Base & >( *this )`
in place of `mMatrix`. Its runtime is identical to a non-virtual member wrapper (§3 shows that
parity directly). What B adds is the vendor API leaking through the public type: in the Blaze probe
`A.capacity()`, `blaze::isSquare( A )` compile on a `belfem::Matrix`, and under Armadillo `A.n_elem`,
`A.is_square()` would — the uncontrolled form of what `matrix_data()` provides today under a name
you can grep for. It also changes overload resolution: BELFEM's dimension-checked `operator*`
(`operators/op_MatrixTimes.hpp:22-35`) and the backend's own `operator*` now both match a
`Matrix * Matrix`. **Verdict, stated at the width the probe supports:** *transparent* inheritance — a derived class that
Armadillo's expression templates discover on their own — fails under Armadillo 15.2 and would need
vendor-internal surgery; it works under Blaze 3.9; a cast-forwarding derived class compiles on both and is
today's design by another spelling. No variant has a runtime gain over a non-virtual member wrapper.
Confidence: high on the compile results (probe), medium (~70%) that no clean Armadillo workaround exists
(two trait hacks tried, Codex's and mine, not an exhaustive search). Probe caveat: under Armadillo
`A = 1.0` resizes to 1×1 (`Mat_meat.hpp:793-800`) rather than filling, so the probe's "scalar fill" line is
Blaze-only; it did not affect the compile result.

**Design C — MORIS-style template parameter.** `Matrix< DDRMat >` instead of `Matrix< real >`.
Zero runtime difference; it renames the type in every consumer file and would let a future
fixed-size specialization (`Matrix< F33RMat >`) coexist with the dynamic one. **Verdict: no benefit
for BELFEM as it stands.** Confidence: high.

**The cost model in one sentence.** A `class Matrix { Backend mMatrix; }` with inline forwarding
compiles to the same instructions as `Backend` itself (shown for the const-read loop in §3 under GCC
15 at `-O2`, both backends; a claim about the tested shapes, not every access pattern); the indirection
is a compile-time name, not a runtime step. The wrapper costs runtime where it introduces something
the member does not have (§3) or lacks an overload the member has, so that an expression converts
through a temporary (§3b).

---

## 3. Where BELFEM's wrapper pays today, first cost: the two `virtual`s

**The code.**

| header | line | member |
|---|---|---|
| `src/linalg/armadillo/cl_AR_Matrix.hpp` | 106 | `virtual ~Matrix() = default;` |
| `src/linalg/armadillo/cl_AR_Matrix.hpp` | 276 | `virtual inline const T & operator()( const size_t, const size_t ) const` |
| `src/linalg/blaze/cl_BZ_Matrix.hpp` | 139 | `virtual ~Matrix() = default;` |
| `src/linalg/blaze/cl_BZ_Matrix.hpp` | 309 | `virtual inline const T & operator()( const size_t, const size_t ) const` |
| `src/linalg/armadillo/cl_AR_Vector.hpp` | 154 | `virtual ~Vector() = default;` |
| `src/linalg/blaze/cl_BZ_Vector.hpp` | 208 | `virtual ~Vector() = default;` |

Nothing uses the polymorphism: no class derives from `Matrix` or `Vector` (grep for
`public Matrix<`, `public belfem::Matrix`, `public Vector<`, `public belfem::Vector` over `src`,
`tests`, `examples`, `nonfree`: 0 hits), no `dynamic_cast` or `typeid` on them (0 hits), and the
non-const `operator()` is already non-virtual. The git history cannot say why they are there: the
headers enter at the squashed `f81dd84` "ASC 2026 release" commit.

**The mechanism.** A class with a virtual member carries a vtable pointer, and a call to a virtual
member through a reference whose dynamic type the compiler cannot prove is an indirect call. GCC at
`-O2` applies speculative devirtualization: it compares the vtable slot with the expected address
and inlines the expected body, with the indirect call as fallback. The assembly of a
`sum( const Matrix<double> & )` loop in this session's probe shows exactly that per element:
load the vptr, load slot 2, `cmpq` against `_ZNK6belfem6MatrixIdEclEmm`, branch, else `call *%rax`
(both backends; the backend-reference version of the same loop has no indirect call). The guard
also blocks interprocedural purity analysis of any function that reads a `const Matrix &`.

**The numbers.** Probe `bench2.cpp` (session scratchpad), `/opt/scls/bin/g++` (GCC 15, SCLS) with
`-O2 -DNDEBUG -std=gnu++17 -m64` and the build's defines (`BELFEM_GCC BELFEM_MPI OMP BELFEM_NETLIB`
plus the backend selector), **without** the build's `-fopenmp -rdynamic -fPIE -fextended-identifiers`
(`cmake-build-debug/tests/linalg/CMakeFiles/test_linalg.dir/flags.make:9`); i9-9880H; best of 7
repeats; single session, not rerun by an auditor; "devirt" = the same four headers with the word
`virtual` removed, nothing else changed.

| quantity | Blaze as-is | Blaze devirt | Armadillo as-is | Armadillo devirt |
|---|---|---|---|---|
| `sizeof( Matrix<double> )` / backend | 48 / 40 | 40 / 40 | 208 / 192 | 192 / 192 |
| `sizeof( Vector<double> )` / backend | 32 / 24 | 24 / 24 | 208 / 192 | 192 / 192 |
| `std::is_polymorphic` | yes | no | yes | no |
| indirect calls in the const-read loop | 1 | 0 | 1 | 0 |
| 1e6-element const read, wrapper / backend | 1.08 / 1.01 ms (1.08×) | 1.01 / 1.04 ms (0.97×) | 1.17 / 1.19 ms (0.99×) | 1.16 / 1.11 ms (1.04×) |
| 3×3 mat-vec through `const Matrix&, const Vector&`, 2e7 calls, wrapper / backend | 300 / 114 ms (**2.64×**) | 119 / 120 ms (**0.99×**) | 435 / 341 ms (**1.28×**, two runs 1.28 / 1.28) | 381 / 365 ms (**1.04×**, two runs 1.04 / 1.07) |

The Armadillo column is from `bench3.cpp`, which uses the bounds-checked `operator()` on **both** sides.
A first version compared the wrapper's `x( j )` with the backend's unchecked `x[ j ]` (`Mat_meat.hpp:6914`
vs `:6937`) and reported 2.00× / 1.66×; Codex caught it in round 1 and the numbers above replace it.
With `-DARMA_DONT_CHECK_CONFORMANCE` (the macro that actually disables Armadillo's bounds checks,
`config.hpp:352-357`, `debug.hpp:1383-1398`) the devirtualized pair measures 399 / 398 ms (1.00×): the
checks cost nothing measurable in this kernel.

**Reading.** On a streaming read the guard costs 0–8 %; on the small-kernel shape of a FEM hot path
(Jacobians, shape-function products, through `const` references, which 325 function signatures in
the tree take — grep for `const Matrix< ... > &` outside `src/linalg`) it dominates, because the
guard and the blocked purity analysis prevent the compiler from treating the 9 reads as it treats
the backend's own. This is a microbenchmark, not a kernel; the magnitude on a real assembly loop is
**unmeasured** and is the first gate in §6. Confidence: high on the mechanism (assembly), medium
(~60%) that a real kernel shows a double-digit percentage.

**Side observations, not part of the verdict.**

- Under Armadillo a `Vector<T>` is a full `arma::Mat` (192 bytes, including the
  `ARMA_MAT_PREALLOC` local store that keeps ≤16-element matrices off the heap); under Blaze a
  `DynamicVector` (24 bytes, always heap). The two backends have different small-object behaviour
  that no wrapper design changes.
- Neither `ARMA_NO_DEBUG` nor `ARMA_DONT_CHECK_CONFORMANCE` is defined anywhere in the build
  (`CMakeLists.txt`, `config/`, `src/linalg`: 0 hits; `linalg_usage_guide.md:98` documents the former as
  deliberate), so Armadillo bounds-checks every `()` in release. Note that `ARMA_NO_DEBUG` would **not**
  turn those checks off — only `ARMA_DONT_CHECK_CONFORMANCE` does (`config.hpp:333-357`). The probe
  measured no cost from them in the 3×3 kernel. Separate decision (O2).

---

## 3b. Second cost: gaps in the overload sets, demonstrated in the tested operations

Raised by Codex in round 1 for Armadillo `+=`; the Blaze `Vector` and the constructor cases were found in
verification and in round 2. **Method and its limit:** the numbers below are allocation and free **call
sites in the function body** of a `-O2 -DNDEBUG -S` compile with the de-virtualized headers (so the two
costs are separated), counted over `call` and tail-call `jmp` to `malloc`/`posix_memalign`/`free`. A
call site is static evidence that a temporary *can* be built; it is not an executed count — a site on
an aliasing branch may never run, and Armadillo temporaries of **≤ 16 elements use local storage**
(`arma_config.hpp:29` `mat_prealloc = 16`, `Mat_meat.hpp:320-324`), so for a 3×3 the Armadillo rows are
a copy, not a `malloc`. Blaze has no such preallocation. R2a's gate is therefore an *executed* counter
(§6), and this table is what tells it where to look.

| expression | Armadillo wrapper / backend (alloc + free sites) | Blaze wrapper / backend | cause |
|---|---|---|---|
| `A += B * dt` (Matrix) | **1 + 2 / 0** | 0 / 0 **under GCC only** — the expression `+=` sits inside `#ifdef BELFEM_GCC` (`cl_BZ_Matrix.hpp:593-601`); a `BELFEM_CLANG`/`BELFEM_INTEL` build converts through the `blaze::Matrix<MT,SO>` constructor (`:80-82`) into a temporary `Matrix` (Grok round 3) | `cl_AR_Matrix.hpp:488-546` has `operator+=` for `const T&`, `T&`, `const Matrix<T>&`, `const MatrixType&` and no `template< ET >` overload; the `eOp` converts to an `arma::Mat` temporary (`Mat_meat.hpp:5461`) instead of reaching `Mat::operator+=( const eOp& )` (`:5525`) |
| `A -= B * dt` (Matrix) | **1 + 2 / 0** | 0 / 0 **under GCC only** (`:625-642`, which also guards the plain `MatrixType` `-=`; `*=` likewise at `:666-683`) | same |
| `y += x * dt`, `y += x + z`, `y += A.col(0)`, `y += A.row(0)` (Vector) | 0 / 0 | **1 + 1 / 0** each | `cl_BZ_Vector.hpp:458-466`: the generic `template< ET > operator+=( const ET & )` is commented out; `:468-470` (and the `-=` twin `:519-521`) accept only `DMatScalarMultExpr` (a *matrix* expression); every vector expression converts to a temporary `Vector` through the constructor at `:128`/`:132` |
| `Matrix C( B * dt )` (construction) | **2 + 3 + `memcpy` / 1 + 1** | **2 + 2 / 1 + 1** | Armadillo: `cl_AR_Matrix.hpp:75-83` constructs from `const MatrixType &` or `arma::Op<ET,OP>` only, so an `eOp` materializes through `Mat(const eOp&)` and is copied into `mMatrix`; the Armadillo `Vector` has the same hole (`cl_AR_Vector.hpp:105-113`: `VectorType` and `arma::Op` only; Grok round 3). Blaze: counts reproduced, overload path not traced — the forwarding constructor `:80-82` → `DynamicMatrix(const Matrix&)` (`DynamicMatrix.h:991-1000`) size-constructs once and assigns, which does not by itself explain a second allocation |
| `Matrix C( trans( B ) )` | 1 + 1 / 1 + 1 | 1 + 1 / 1 + 1 | no gap: `arma::Op` has a wrapper constructor; Blaze's `blaze::Matrix< MT, SO >` constructor matches |
| `A = B * dt + C` (Matrix) | 0 / 0 | 1 + 1 / 1 + 1, on the `canAlias` branch only (`DynamicMatrix.h:4966-4970`) | not a wrapper cost |

**What this means.** The wrapper's second cost is not one hole but a *pattern*: wherever its overload
set lacks the expression overload the backend has, the expression converts through the backend-type
constructor or the `const MatrixType &` overload and a temporary is built. The three gaps above are the
ones demonstrated; R2a audits the full set (constructor, `=`, `+=`, `-=`, `*=`, both classes, both
backends) rather than fixing these three. These are hidden allocations in exactly the shapes the time
integrator writes once per element per step, and they violate the no-temporary rule of
`doc/coding_philosophy.md` §"No temporary `Vector`/`Matrix`" wherever the size exceeds 16 (Armadillo) or
always (Blaze). **They are also why the escape hatches exist:** `cl_IWG_Timestep.cpp:851-866` writes
`aJ.matrix_data() += mTimeStepMatrices->K().matrix_data() * mDeltaTime` — the hatch reaches the backend's
allocation-free `operator+=`, which the wrapper denies it. So the 109 hatch lines of §5 are not only
friction; a share of them is a deliberate performance workaround, and removing them without first
closing the gaps would *add* allocations. **The first R2a recipe — restore the commented-out generic
Blaze `+=`/`-=` verbatim — is withdrawn:** probed in round 2, it makes `y += A.row(0)` fail to compile
(a row expression forwarded into a column `DynamicVector` without the orientation change the constructor
performs). The replacement is constrained overloads: `blaze::DenseVector< VT, false >` for column
expressions and a transposing overload for `DenseVector< VT, true >`, with compatibility tests for row
and column views, scalars, wrapper operands and expressions.

Confidence: high on the call-site rows (assembly, re-derived by Codex for the Armadillo rows); the
executed cost at the sizes BELFEM uses is **unmeasured** (R2a's gate); the share of hatch sites that are
performance workarounds versus convenience is not counted (medium, ~60% that `cl_IWG_Timestep.cpp`'s
79 lines are mostly the former).

**Compiler dependence (Grok round 3, P1).** The table's Blaze `Matrix` rows describe a `BELFEM_GCC` build,
which is what the §3 probe passed and what every tree on this machine is (`cmake-build-debug/CMakeFiles/belfem_lib.dir/flags.make`:
`-DBELFEM_GCC`, SCLS `g++` on Apple). BELFEM also configures `BELFEM_CLANG` (`config/compiler/config_gcc.cmake:4-8`)
and `BELFEM_INTEL` (`config_icc.cmake:3`); on those the three guarded blocks in `cl_BZ_Matrix.hpp` vanish and
`A += B * dt` goes through a temporary `Matrix`. The guard's reason is lost — it arrived in the squashed
2026-09-05 release commit. Probed 2026-11-07 (`-fsyntax-only`, Apple clang 17, `-DBELFEM_CLANG`, Blaze):
the unguarded overload set compiles; Intel is untested here. R2a therefore includes lifting the guards.

---

## 4. Verdict on (a) and (b)

**(a) Is aliasing cheaper?** No. Relative to a non-virtual, complete member wrapper it is
runtime-identical. Relative to BELFEM's *current* wrapper there are two demonstrated differences, the
`virtual`s (§3) and the overload gaps (§3b), and both are defects of this wrapper, not of wrapping:
ten declaration lines remove the first, an audited overload set removes the second, neither touches
the design. Whether those two are the *only* differences is not established (§9). Confidence: high on
the two; medium on exhaustiveness.

**(b) Can BELFEM apply it, and is it worth it?** Pure aliasing (A): no, it dissolves the two-backend
contract and the wrapper policy. Inheritance (B): transparent inheritance works under Blaze and fails
under Armadillo 15.2, the Linux default; the cast-forwarding variant that compiles on both is today's
design renamed. Backend-type parameter (C): possible, pointless. None of the
three is worth a change that touches the type named in 552 consumer files. **What is worth doing is
de-virtualizing the four headers and completing their operator sets.** Confidence: high.

---

## 5. Impact radius (c)

All counts are `grep` over `src/` and `nonfree/` (`*.hpp`, `*.cpp`), excluding `src/linalg/` unless
stated. "Lines" are matching lines, not call sites; a line can hold two.

| what | count | scope note |
|---|---|---|
| source files in `src` + `nonfree` | 1058 | `find` |
| files outside `src/linalg` mentioning `Matrix<` or `Vector<` | **552** (498 `src`, 54 `nonfree`); 64 more inside `src/linalg` | by module: fem 180, physics 85, numerics 81, mesh 51, math 38, sparse 24, homology 20, io 5, circuit 4, visualizer 3, executables 3, containers 2, comm 2 (sum 498) |
| `Matrix<` lines / `Vector<` lines | 1406 / 3596 | matching lines, not receivers |
| scalar types instantiated (occurrences outside linalg) | `Vector<real>` 2525, `Matrix<real>` 1238, `Vector<id_t>` 480, `Vector<index_t>` 362, `Matrix<int>` 96, `Vector<uint>` 81, `Vector<int_t>` 43, `Matrix<index_t>` 35, `Vector<metis_t>` 21, `Vector<scotch_t>` 15, `Matrix<id_t>` 13, `Vector<cplx>` 3, … | textual occurrences, not template instantiations; every integer and complex instantiation must survive any redesign |
| `const Matrix< … > &` lines (whitespace-tolerant) | 343 | not all are parameters; the hot-path signature shape of §3 |
| `.n_rows()`/`.n_cols()` lines / `.length()` lines | 244 / 717 | `length()` includes `Cell`-free contexts only where the name is linalg's; not separated |
| `.row(` / `.col(` / `.submat(` lines (backend views) | 161 | views are backend types already (`auto`); 3 `auto` holders found |
| `.set_row(` / `.set_col(` lines | 162 | |
| `.data()` lines (all containers) | 308 | `Cell::data()` not separated |
| `.spacing()` outside linalg | 4 | all in `src/comm/commtools.hpp` (`:1534, :1593, :1751, :1815`) |
| `.spacing()`/`.data()` inside `src/linalg/lapack` | 55 | the LAPACK wrappers (`lapacktools.hpp`, `fn_gesv.hpp`, …) consume the BELFEM API, not the backend |
| **linalg escape-hatch lines** `matrix_data()`/`vector_data()` on `Matrix`/`Vector` | **110 in 13 files** (109 in the 12 editable files + 1 in the closed module) | `cl_IWG_Timestep.cpp` 79, `cl_FEM_Calculator.hpp` 9, `cl_Gas.cpp` 9, `cl_GM_EoS_Cubic.cpp` 3, `cl_Material_UserDefined.{hpp,cpp}` 3, `cl_Gradient.cpp` 1, `fn_Mesh_compute_edge_lengths.cpp` 1, `cl_FEM_DofMgr_EigenValues.cpp` 1, `cl_Cohomology.cpp` 1 (closed module), `fn_ODE_DOP853.cpp` 1, `fn_ODE_RK45.cpp` 1, one nonfree file 1 |
| other `vector_data()`/`matrix_data()` lines outside linalg (`Cell`'s `std::vector` accessor, `Spline`, `ProtoMesh`) | 45 (= 155 total − 110) | **not** linalg; listed so the two are not confused |
| linalg module itself | 4 class headers, 8 operator headers, 26 `fn_*` dispatch headers + 18 backend `fn_*` + 6 backend `op_*`, 10 LAPACK files, 2 guides (1002 + 320 lines) | `src/linalg`, `src/linalg/doc` |
| tests | `tests/linalg` 4184 lines (`test_Matrix` 920, `test_Vector` 845, `test_LinalgSolvers` 943, `test_lapack` 1434), label `fast` | `MatrixDataExposesBackend` (`test_Matrix.cpp:192`) pins the escape hatch |

The escape-hatch sites are instructive. `cl_IWG_Timestep.cpp:755-790, 850-1162` writes
`aSwap.vector_data() = mBeta(0) * mCalc->qold(0).vector_data()` and
`aJ.matrix_data() += M().matrix_data()` where the wrapper operators should suffice;
`cl_FEM_Calculator.hpp:1724-1757` assigns a backend expression to `matrix().matrix_data()` instead
of to the `Matrix`; `cl_Gradient.cpp:196` chains `inv( mJ.matrix_data() ) * dNdXi.matrix_data() * mPhi.vector_data()`.
The first version of this report read these as overload friction with no runtime effect. Round 1
refuted that: the `+=`/`-=` forms reach the backend's allocation-free compound assignment, which the
wrapper's overload set denies them (§3b). The `cl_Gradient.cpp` and `cl_FEM_Calculator.hpp` forms are
plain assignments through the wrapper's `template< ET > operator=`; whether each is temporary-free
was **not** probed (the Blaze `inv()` in the Gradient line is evaluated on purpose, `fn_BZ_inv.hpp:24`,
and stays that way), so the split between workaround and convenience across the 110 lines is not
counted here.

Three radii follow, one per option in §6:

- **R1, de-virtualize + `final`:** 10 declaration lines in 4 headers (six `virtual` removed, four
  `final` added); **zero in-tree call-site changes**, but a deliberate **downstream API break** for any
  out-of-tree class that derives from `Matrix` or `Vector` (none known) and an **ABI change** for every
  translation unit that sees the layout: −8 bytes per object under Blaze, −16 under Armadillo, so the
  whole tree is rebuilt. Plugins are **not** in that set by contract: the plugin API is backend-free
  (`UserLibraryTemplate.cmake:42-46`, `cl_Material.hpp:24`, enforced by the DR-37 gate in
  `tests/physics/backendfree/`), so a conforming plugin never sees a `Matrix`; a plugin that violates
  the contract is outside the supported surface already. `final` turns any downstream derivation into
  a compile error at the first rebuild; a `static_assert( !is_polymorphic )` would not (it inspects only
  the base). `tests/linalg` + `make check` on both backends is the gate.
- **R2a, audit and complete the overload sets:** 4 class headers; the gaps demonstrated so far are
  `cl_AR_Matrix.hpp` constructor-from-`eOp`/`+=`/`-=` and `cl_BZ_Vector.hpp` `+=`/`-=`, plus a Blaze
  construction case whose path is untraced. Gate: **executed** allocation counters (not assembly
  sites) per expression shape, at sizes below and above 16 elements, with and without aliasing, both
  backends — wrapper form must allocate exactly what the backend form does — plus `make check`.
- **R2b, retire the hatches (after R2a):** 12 editable files, 109 lines, plus 1 line in the closed
  cohomology core (report to Gregory, do not touch). Gate per replaced expression: the R2a counter
  **and** a named numerical test of the path. Today's coverage of the timestep paths is thin:
  `tests/fem/test_BdfTimestepMethod.cpp:25-27` records that nothing exercises BDF2–5 end to end;
  `test_QoldTable.cpp` is the only other file naming the class. The Crank–Nicolson and Galerkin
  branches are **deliberately disabled** by `BELFEM_ERROR` (`cl_IWG_Timestep.cpp:84-94`, `:162-169`: the
  Newton tangent is exact for the BDF family only), so their hatch lines are dead code behind a
  rejection the `:313` test protects — they are not an R2b prerequisite and are left as they are.
- **R3, inheritance (not recommended):** `src/linalg` only *if* the public API stays
  source-compatible — every signature, including `matrix_data()` (which can return `*this` by
  reference) — then the 552 consumer files, the LAPACK wrappers and `commtools.hpp` compile unchanged;
  the object layout changes regardless (an ABI change, full rebuild); otherwise the radius is the whole
  tree. Plus `tests/linalg`, the two guides, `doc/coding_philosophy.md` §"Container
  Selection" and `CLAUDE.md` §"Matrix Memory Layout". And it is Blaze-only (§2).

---

## 6. Implementation plan (conditional on Christian's decision)

Every step is read-only until approved. R1 and R2a are recommended; R2b is conditional on R2a; R0 is
the measurement that sizes R1's benefit and does not gate R2a.

- [ ] **R0 — measure first (no change to the tree).** In an isolated build tree, build as-is and with
  the §3 header edit applied, and time one representative assembly (`belfem` on an `examples/` deck
  with element-matrix work dominating, both backends). Any layout-dependent code loaded into that run
  must be built from the same headers (the loader has no version handshake; conforming plugins are
  backend-free and unaffected). Record wall time and the `check-fast` result. Falsifier: if the real
  kernel shows < 2 % difference, R1 is still correct (it removes a vptr and a guard that buy nothing)
  but stops being a performance item and becomes a hygiene item.
- [ ] **R1 — de-virtualize `Matrix` and `Vector` and mark them `final`** (plan jury + code jury per
  `code-with-jury`; critical-framework row, Codex `gpt-6-astra`). Remove the six `virtual` specifiers at the lines in
  §3; add `final` to the four class heads (O4). Add to `tests/linalg/test_Matrix.cpp` and
  `test_Vector.cpp` a `static_assert( !std::is_polymorphic< belfem::Matrix< belfem::real > >::value, ... )`
  (and `Vector`) as a regression tripwire for the base itself — knowing it does not see a downstream
  derived class; `final` covers that. Gate: `make check` on Blaze (Apple) and Armadillo (Linux CI),
  plus an assembly probe showing no `call *` in a `const Matrix &` read loop (the §3 probe, checked in
  under `tests/linalg` as a compile-only TU, or kept as the R0 artifact). Document in
  `linalg_usage_guide.md` that the classes are not base classes. Rebuild every dlopen'ed plugin.
- [ ] **R2a — audit and complete the overload sets** (same jury row). For each of the four classes,
  list every backend expression overload of the constructor, `=`, `+=`, `-=`, `*=` and the wrapper
  overload that reaches it without a temporary; close each gap with a **constrained** overload.
  Known gaps: `cl_AR_Matrix.hpp` needs `template< ET > operator+=/-=( const ET & )` (mirroring the
  Armadillo `Vector`) and a constructor for `eOp`/`eGlue`/`Glue` expressions; `cl_AR_Vector.hpp:105-113`
  needs the same constructor; `cl_BZ_Matrix.hpp:593-683` must lose its `#ifdef BELFEM_GCC` guards so
  Clang and Intel builds get the same overload set (Grok round 3; Apple clang 17 compiles the unguarded
  set, probed); `cl_BZ_Vector.hpp` needs
  `operator+=/-=( const blaze::DenseVector< VT, false > & )` for column expressions and a transposing
  overload for `DenseVector< VT, true >` — **not** the commented-out generic `template< ET >` at
  `:458-466`, which makes `y += A.row(0)` fail to compile (probed in round 2). Gate: an **executed**
  allocation counter (a `malloc`/`operator new` hook or `malloc_zone` statistics around the expression,
  setup excluded) per shape in §3b, below and above 16 elements, with and without aliasing, both
  backends, **for `+=` and `-=` alike** (the §3b table measured `+=`; the `-=` twins at
  `cl_BZ_Vector.hpp:509-521` have the same shape and are not measured); compatibility tests for row and
  column views, scalars, wrapper operands and expressions; `make check`.
- [ ] **R2b — retire the escape hatches** (after R2a). First an inventory of the 109 lines: operation,
  operand shapes, live or dead path, workaround or convenience — the split §3b leaves uncounted. Then,
  for each live site, replace the hatch with the plain wrapper expression; gate each with the R2a counter and a
  named numerical test of the path — for `cl_IWG_Timestep.cpp` that means first adding the BDF2–5
  assembly tests that `test_BdfTimestepMethod.cpp:25-27` says do not exist. The Crank–Nicolson and
  Galerkin hatch lines sit behind a deliberate `BELFEM_ERROR` and are left alone. The cohomology site
  is reported to Gregory only.
- [ ] ~~**R3 — inheritance prototype behind a CMake flag.**~~ Struck in this report: fails to compile
  under Armadillo 15.2 (§2). Un-strike only if O1 is decided as "Blaze only".

**Open questions for Christian (logged, not decided):**

- **O1** Does BELFEM keep both backends? The Armadillo failure in §2 is what rules out design B; a
  Blaze-only BELFEM could revisit it (and would then lose Armadillo's ≤16-element preallocation and
  SuperLU/ARPACK glue, `config/linalg/config_matrix.cmake:15`).
- **O2** `ARMA_DONT_CHECK_CONFORMANCE` in release builds (not `ARMA_NO_DEBUG`, which leaves the
  checks on): it removes the bounds checks on `()` **and** the same-size, multiplication-size and
  BLAS-size conformance checks (`debug.hpp:1383-1398`), i.e. the Armadillo-side detection of a shape
  mismatch in release. The probe measured no 3×3-kernel effect (1.00× with the checks off, §3). Low
  priority; a timing alone cannot decide it — inventory what each removed check would have caught
  and what BELFEM's own `BELFEM_ASSERT`s (compiled out in release) leave uncovered.
- **O3** **INVESTIGATED 2026-11-07 → `todo/linalg_small_matrix_kernels.md`**: fixed sizes buy 12 % of the Jacobian
  chain on Blaze and nothing on Armadillo, whose chain is 5.3× slower because its non-square small products
  are BLAS calls (`mul_gemm.hpp:413`); the plan is backend-free small-product kernels, not a second type.
  Original text: fixed-size small matrices (`arma::Mat<T>::fixed<3,3>` / `blaze::StaticMatrix<T,3,3>`) for
  2×2 and 3×3 Jacobians — the one MORIS idea (Eigen `F33RMat`, `linalg_typedefs.hpp:49`) that could
  pay in BELFEM's calculator. Out of scope here; would need the §2 design-C parameterization or a
  second template.
- **O4** `final` on `Matrix` and `Vector` — **proposed resolution: yes, as part of R1.** It is the
  only mechanism that makes a downstream derived class a compile error (Codex round 1); it costs
  nothing and forbids derivation, which nothing in `src`, `tests`, `examples`, `nonfree` does today.
- **O5** The `final` API break is a deliberate compatibility decision, not a side effect: R1 ships only
  with Christian's explicit acceptance of it, recorded with the plan approval. The one compatibility
  fact that bounds it: the plugin API is backend-free by contract (`UserLibraryTemplate.cmake:42-46`,
  `cl_Material.hpp:24`, DR-37 gate), so no conforming plugin can be affected; anything else out of tree
  that derives from the classes is unknown to this report. The layout change reaches wider than
  derivation: any out-of-tree translation unit that holds a `Matrix`/`Vector` **by value** must be rebuilt,
  derived or not (Grok round 3) — `final` does not diagnose that; only a clean rebuild does.
- **Doc fix, no jury owed:** `src/linalg/doc/linalg_usage_guide.md:98` says the Armadillo release
  bounds checks stay on "because BELFEM leaves its conformance checks enabled (`ARMA_NO_DEBUG` is not
  defined)". The parenthetical names the wrong switch: `ARMA_NO_DEBUG` leaves the checks on
  (`config.hpp:333-336`); `ARMA_DONT_CHECK_CONFORMANCE` is what would turn them off (`:352-357`).

---

## 7. What would change the verdict

- A MORIS commit or branch that actually removes `mMatrix` (none found; §1 enumerates the search).
  Then re-evaluate against that design, not against the idea of one.
- An Armadillo mechanism that lets a derived class pass `is_arma_type` without vendor-internal
  surgery. Then design B becomes backend-symmetric and the trade is API leakage vs. nothing.
- R0 showing no measurable difference on a real kernel. Then R1 is hygiene, not performance, and
  this report's §3 magnitude claim is wrong in the direction of overstatement. R2a's case does not
  depend on R0: a hidden temporary per element per step is a hot-path allocation by the philosophy's
  own rule, whatever it costs — with the caveat that below 16 elements Armadillo's is a copy into
  local storage, so the executed counter decides which shapes are heap allocations.
- The executed allocation counter (R2a's gate) showing that the §3b call sites never run at
  BELFEM's element sizes. Then the gaps are a latent hazard, not a present cost, and R2a drops in
  priority behind R1.
- A Clang or Intel build in which the Blaze `Matrix` `+=` cell of §3b is still 0 / 0 — it cannot be,
  the overload is absent there (`cl_BZ_Matrix.hpp:593`); conversely, if lifting the guard breaks a
  compiler BELFEM supports, R2a needs a per-compiler path rather than one overload set.
- A consumer that relies on `Matrix`/`Vector` polymorphism (none found in `src`, `tests`,
  `examples`, `nonfree`; a downstream user plugin could). R1's `final` turns that into a compile
  error at the plugin's first rebuild; the `static_assert` in `tests/linalg` would not see it.

---

## 8. Literature

Not a literature-first task (protocol §5). Two sources were checked for the wrapper rationale:
Oliveira & Stewart 2006 §7.6 on library-interface design and scratch-space policy (the argument for
a thin, uniform interface over vendor calls) and the passage noting that expression-template
libraries exist to avoid temporaries (`oliveira2006.txt:6983`); Rouson, Xia & Xu 2011 on wrapper
interfaces being both safer and simpler than direct legacy calls (`rouson2011.txt:1502-1503`).
Neither book addresses `virtual` in value classes; that is C++ mechanics, settled here by the probe.

---

## 9. Confidence summary

| claim | confidence | evidence rung (protocol §11) |
|---|---|---|
| MORIS still wraps a member; no branch differs | high | static source trace, enumerated |
| The `virtual`s cost a vcall guard per const read and up to 2.6× on a small kernel through const refs | high on the mechanism (assembly, re-derived by Codex); medium on the ratio (single-session timing, GCC 15, -O2, one CPU) | numeric probe |
| The overload gaps build a temporary (the §3b rows) | high that the call sites exist (assembly, Armadillo rows re-derived by Codex); **unmeasured** whether they execute at BELFEM's sizes | static call-site count, below the probe rung for execution |
| Those two are the *only* wrapper runtime costs | low | three gaps found in two rounds after one was claimed; the question is open until R2a's audit |
| Magnitude on a real BELFEM kernel is double-digit % | medium (~60%) | microbenchmark only; R0 is the gate |
| Alias (A) incompatible with two backends + L-21 | high | source trace of the API differences |
| Inheritance (B) fails under Armadillo 15.2 | high | compile probe |
| No clean Armadillo workaround for B exists | medium (~70%) | two trait hacks tried (`is_arma_type`, `Proxy`); Grok round 3 read `Base`, `Proxy`, `unwrap` and found no supported extension point |
| The demonstrated `+=`/`-=` hatch lines (`cl_IWG_Timestep.cpp:851`, `:866`, `:921`) are performance workarounds for §3b | high | probe + source trace; the class as a whole is the next row (lowered in round 3) |
| The Blaze `Matrix` compound-assignment set is complete | high **under GCC only**; absent under `BELFEM_CLANG`/`BELFEM_INTEL` | source (`#ifdef BELFEM_GCC`, `cl_BZ_Matrix.hpp:593-683`); Apple clang 17 compiles it unguarded (probe) |
| The share of hatch sites that are workarounds vs convenience | low–medium | not counted |
| Impact-radius counts | high — recounted in round 2 with the stated scopes (six corrections folded in); still textual matches, not receiver types | grep, scopes stated |

---

## 10. Audit trail

- Exchange: `tmp/ai_exchange/review_linalg_wrapper_design.md` — Claude pre-registration, round-1
  jury (Codex `gpt-6-astra`/high, Grok `grok-4.7`/high), verification, reconciliation; round 2 on
  the revised document at `xhigh`.
- Probes (session scratchpad, not in the tree): `probe_vtable.cpp` (sizes, `is_polymorphic`,
  assembly), `bench2.cpp` / `bench3.cpp` (timings, both backends, as-is vs devirt headers; `bench3`
  uses identical accessors), `temp_probe.cpp` (allocation counts of §3b), `inherit_probe.cpp` and
  `inherit_probe_traits.cpp` (design B, both backends). `bench2` and `inherit_probe` are reproduced in
  the exchange pre-registration, `temp_probe.cpp` and the `bench3` accessor change in the round-1
  reconciliation entry, so every number here can be re-run from the exchange alone. R1 and R2a check
  the two gate probes into `tests/linalg`, which retires this caveat.
- Round-1 record: Codex `gpt-6-astra`/high, six findings, all CONFIRMED and folded in; Grok seat
  failed (CLI not signed in); Copilot reserve seat substituted (`auto` → `gpt-5.6-luna`), four findings,
  all CONFIRMED. Round-2 record: Codex `gpt-6-astra`/xhigh on v3, nine findings, all CONFIRMED and folded
  in (one withdrawn recipe, one new gap, six count corrections); Copilot reserve seat on v3, six findings,
  all CONFIRMED (two duplicating Codex), folded in as O5, the guide doc fix, the `-=` gate and the R2b
  inventory. Both rounds had two OpenAI-hosted voices and no Grok, so the vendor independence the jury is
  bought for was reduced; a Grok pass is owed once the CLI is signed in.
- Round-3 record (2026-11-07 15:03 CST): Grok `grok-4.7`/xhigh, blind, on v5 — ten findings, all
  CONFIRMED against source in the verification entry; one premise (the Apple gate being a Clang build)
  REFUTED for this machine (SCLS `g++`, `-DBELFEM_GCC`), the compiler-dependence itself stands. P1:
  `cl_BZ_Matrix.hpp:593-683`. Probe added: `clang_probe/probe.cpp` (four configurations, all compile).
  The jury has now had a non-OpenAI voice; the vendor-independence caveat of rounds 1–2 is lifted for v6.
- Probes added in round 2: `row_probe.cpp` (Blaze `Vector` `+=` with views, scalars, sums; and the
  withdrawn generic overload), `ctor_probe.cpp` (construction from an expression, both backends); their
  text is in the round-2 verification entry.
