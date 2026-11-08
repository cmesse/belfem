# linalg Wrapper Hardening: De-virtualize, Check Policy, Overload Sets, Hatch Retirement

**Date:** 2026-11-07
**Purpose:** Remove the two measured costs of BELFEM's `Matrix`/`Vector` wrapper found in
`todo/linalg_wrapper_design_investigation.md` — six `virtual` members nothing derives from, and
overload-set gaps that force temporaries — and align the backend size/bounds checks with the build
mode (on in debug, off in release). Mechanism: rounds A–C edit the four backend class headers, one backend configuration header and
`tests/linalg`; round D edits the consumers that bypass the wrapper. **Four implementation rounds**
(virtuals → check policy → overload sets → consumers), each with its own code jury and executable
gate, on **both** backends.
**Module:** `src/linalg` (rounds A–C); `src/fem/iwg`, `src/fem/kernel`, `src/fem/postproc`,
`src/physics`, `src/numerics`, `src/mesh` (round D, consumers)
**AIs involved:** Claude (plan), Codex `gpt-6-astra` (plan + code audits), Grok `grok-4.7` (third
voice; restored 2026-11-07 after two rounds on the report without it)
**Status:** IN PROGRESS — **rounds A, B and C VERIFIED (Christian's `make check`, Blaze and Armadillo, 2026-11-08) and committed; round D next (R4.0 hatch inventory, read-only, then per-file sub-rounds).** O-decisions as recorded below. Plan jury record: **round 2 at `xhigh` done** (Codex: 5 findings, 3 P1 — eigen status discarded, complex transpose
semantics, stale §6 counter sketch; Grok: 4 P1 + 8 P2 — one-TU tripwire, API wording, row-view add is a transpose-add not UB, …; all
CONFIRMED, two by run probe) and folded in (gaps #36–#37, R2.0, R3.0, O12–O13). O2a, O2b, O4, O5 **decided 2026-11-07 by Christian**
(below). Drafted 2026-11-07, **plan jury round 1 done the same day** (Codex `gpt-6-astra`/high:
8 findings; Grok `grok-4.7`/high: 13 findings + open risks; all CONFIRMED except one REFUTED by probe; three P1s
folded in as gaps #32–#34 and the rewritten #18). **No source modified**; source edits begin only after round 2 is
reconciled and Christian approves. Grok's blind audit of the investigation report (round 3 of that thread, `grok-4.7`/xhigh, 2026-11-07)
is **folded in**: gaps #16, #27–#30 and step R3.7 come from it.

> **Scope guards (from the task brief, 2026-11-07, Christian):**
> - **Both backends.** Every step is specified for Armadillo (`src/linalg/armadillo/`) and Blaze
>   (`src/linalg/blaze/`) and is gated on both: Blaze is the Apple default, Armadillo the Linux
>   default (`CMakeLists.txt:92-97`; `config/linalg/config_matrix.cmake` only rejects both-or-neither), so one workstation never exercises both — the
>   Linux CI matrix (debug/release × Armadillo/Blaze, `.gitlab-ci.yml`) is part of every gate.
> - **Order:** the virtuals first, the consumers last. Rounds B and C sit between because round D
>   needs the overloads of round C, and round C's gate needs the check policy of round B settled.
> - **Decided 2026-11-07, Christian:** size-conformance checks **on in debug, off in release**
>   (resolves report O2 — see O2a/O2b below for the two residuals this plan cannot decide alone).
> - **Out of scope:** the MORIS-style container-less designs (report §4: alias, inheritance,
>   backend-type parameter — all rejected or struck); fixed-size small matrices (report O3);
>   the `moris::Vector`/`Cell` question; the review-script defects (`ask_copilot.sh` bash 3.2,
>   `cross_review.sh` Grok-auth reporting) — tooling scripts owe their own jury round.
> - **Compatibility:** two API changes, both deliberate. (1) Both classes become `final` (report
>   O4/O5, a downstream API break with no in-tree or `nonfree` user —
>   grep today: no class derives from either, no `dynamic_cast`/`typeid`/`is_polymorphic` on
>   either). Object layout changes (−8 B Blaze, −16 B Armadillo per object): full rebuild of every
>   tree and every dlopen'ed plugin; conforming plugins are backend-free by contract
>   (`src/physics/materials/UserLibraryTemplate.cmake:42-46`, `cl_Material.hpp:24`,
>   `tests/physics/backendfree/`). (2) Round C **adds** constructors and compound-assignment overloads;
>   every existing expression call site rebinds to them (e.g. the four `aJ += K().matrix_data()*dt`
>   lines and `aRHS -= K().matrix_data() * qold().vector_data()` in `cl_IWG_Timestep.cpp`), which is
>   a behaviour change across the tree, gated by the round-C tests and `make check` on both backends.
>   Rebuild rule: every translation unit that includes the linalg headers, in every tree; backend-free
>   plugins are unaffected by contract.
> - **Closed to AI:** the one hatch line in `src/homology/cl_Cohomology.cpp` is reported to Gregory,
>   never edited. The one `nonfree` site is handled in a `nonfree/` session with no exchange tier.
> - **Christian runs builds and binaries.** Every gate below names what he runs; Claude syntax-gates
>   (`-fsyntax-only`) and says "reviewed", never "verified", until the gate has run.

---

## 1. Current Behaviour and How It Fails

Findings of the investigation report (v5, two jury rounds, 25 findings all confirmed), re-anchored
today:

| Failure | Mechanism | Evidence |
|---|---|---|
| Hot-path indirection | Six `virtual` members: the destructor in all four classes and `operator()(i,j) const` in both `Matrix` classes. GCC speculatively devirtualizes through `const Matrix &` (vptr load, slot compare, `call *` fallback) and the class carries a vptr | `cl_AR_Matrix.hpp:106,:276`, `cl_AR_Vector.hpp:154`, `cl_BZ_Matrix.hpp:139,:309`, `cl_BZ_Vector.hpp:208`; report §3: 3×3 mat-vec 2.64× (Blaze) / 1.28× (Armadillo) slower than the backend, parity (0.99× / 1.04×) with the `virtual`s removed |
| Temporaries where the backend has none | Overload gaps: Armadillo `Matrix` has `+=`/`-=`/`*=` only for `T`, `Matrix<T>` and `MatrixType` (`cl_AR_Matrix.hpp:489-572`), so `A += B*dt` converts the expression to a `MatrixType` temporary; the Armadillo `Vector` has the `template<ET>` forms (`cl_AR_Vector.hpp:399-401,:438-440,:491-493`). Blaze `Vector` has its generic `+=`/`-=` commented out (`cl_BZ_Vector.hpp:441,:460,:501,:511`) and keeps only a `DMatScalarMultExpr` overload (`:468-470,:519-521`), so `y += x*dt` builds a temporary. Direct-initialization from an expression (`Matrix<real> C( B*dt )`; copy-init `C = B*dt` is rejected today, two user-defined conversions) goes through the `MatrixType` converting constructor on Armadillo (`cl_AR_Matrix.hpp:75`; only `arma::Op<ET,OP>` has a direct ctor, `:80-81`) | report §3b: Armadillo `A += B*dt` 3 allocation sites (wrapper) vs 0 (backend); Blaze `y += x*dt` 2 vs 0; static call-site counts, executed counters are round C's gate |
| Consumers bypass the wrapper to get the backend's path | 110 `matrix_data()`/`vector_data()` lines in 13 files, 79 of them in `cl_IWG_Timestep.cpp` (`collect_qhist` `:741-794`, `bdf1`–`bdf5`/`*_nok` `:909-1122`, …), written as `aJ.matrix_data() += M().matrix_data()` where the wrapper operator should have sufficed | report §5; the shapes are 58 `+=`, 7 `-=`, 8 `*=`, 6 `=` in that file (grep today) |
| Check policy differs by backend | Blaze: bounds checks in `()` are `BLAZE_USER_ASSERT` (`DynamicMatrix.h:1055-1056`, `DynamicVector.h:944`), switched on only without `NDEBUG` by `blaze_config.hpp:20-37`; size-conformance is an always-on `BLAZE_THROW_INVALID_ARGUMENT` (`DynamicMatrix.h:1667,:1707`, `DMatDMatAddExpr.h:1092-1093`, `DVecDVecAddExpr.h:1096-1097`). Armadillo: bounds **and** conformance both hang on `ARMA_CHECK_CONFORMANCE`, defined by default (`armadillo_bits/config.hpp:180-181`, `debug.hpp:1383-1405`) and never touched by `armadillo.hpp:20-32`, so a release build checks every `A(i,j)` and every expression on Armadillo, while release Blaze checks shapes (exceptions) but not bounds | `ARMA_NO_DEBUG` would not change this (`config.hpp:333-336`); `ARMA_DONT_CHECK_CONFORMANCE` would (`:352-357`), and `ARMA_DEBUG` re-enables everything (`:338-350`). BELFEM's own debug bounds checks exist on both backends (`BELFEM_ASSERT` in `operator()`, e.g. `cl_AR_Matrix.hpp:258-263`, `cl_BZ_Matrix.hpp:291-296`) |
| Guide states the wrong switch | `linalg_usage_guide.md:97-98` names `ARMA_NO_DEBUG` as the reason bounds checks stay on in release | `config.hpp:333-336` |

**Bottom line:** the wrapper is the right design and costs nothing where it is complete and
non-virtual; it costs where it is polymorphic (six declarations) and where its overload set is
thinner than the backend's (two classes), and consumers have been routing around the second defect
by hand for years. All three are fixable inside `src/linalg`; only the clean-up of the consumers
reaches outside it.

## 2. Architecture: Why the Member Wrapper Stays the Spine

The one-member, non-virtual, inline wrapper compiles to the backend's own code (report §3, probes on
both backends). The three container-less alternatives were weighed and rejected in the report: a
pure alias dissolves the two-backend contract and the L-21 wrapper policy; transparent inheritance
compiles under Blaze 3.9 but not under Armadillo 15.2, whose expression templates key on exact
types (`traits.hpp:542-571`); MORIS's backend-type template parameter changes nothing at runtime.
Nothing here is new design — each round removes a defect from the design BELFEM already has.

## 3. Gap Table

Class: **(a)** mechanical from the citation, **(b)** open question (On), **(c)** needs explicit
handling or a probe.

| # | State | Needed for | Handled today? | Class | Citation / rationale |
|---|---|---|---|---|---|
| 1 | `virtual ~Matrix()` ×2, `virtual ~Vector()` ×2 | A | present | (a) | `cl_AR_Matrix.hpp:106`, `cl_AR_Vector.hpp:154`, `cl_BZ_Matrix.hpp:139`, `cl_BZ_Vector.hpp:208` |
| 2 | `virtual inline const T & operator()(i,j) const` ×2 | A | present | (a) | `cl_AR_Matrix.hpp:276`, `cl_BZ_Matrix.hpp:309` — the non-const twin (`:256`, `:288`) is not virtual |
| 3 | `final` on the four class heads | A | absent | (a) — O4/O5 resolved 2026-11-07 | `cl_AR_Matrix.hpp:28`, `cl_AR_Vector.hpp:32`, `cl_BZ_Matrix.hpp:27`, `cl_BZ_Vector.hpp:40` |
| 4 | Regression tripwire: `static_assert(!std::is_polymorphic<…>)` for both classes, both scalar families | A | absent | (a) | lands in `tests/linalg/test_Matrix.cpp`/`test_Vector.cpp`; guards the base only — `final` guards derivation |
| 5 | `sizeof` tripwire: `sizeof(Matrix<real>) == sizeof(backend type)` | A | absent | (a) | the report's layout numbers (48→40 Blaze, 208→192 Armadillo) become a compile-time check |
| 6 | Armadillo conformance/bounds off in release | B | **off never** | (a) | `armadillo.hpp` gains the tree's debug predicate `#if !defined(NDEBUG) \|\| defined(DEBUG)` (as `blaze_config.hpp:20` and `assert.hpp:22-30`) and defines `ARMA_DONT_CHECK_CONFORMANCE` on the release side, before `#include <armadillo>`. `DEBUG` **is** defined under `USE_DEBUG` (`config_gcc.cmake:23-24`, `config_icc.cmake:17-18`, `finalize_compiler.cmake:48-49`) and `NDEBUG` otherwise — the plan's first draft said `DEBUG` was unused; Grok round 1 refuted that |
| 7 | Armadillo prints `note: #pragma message: WARNING: conformance checks disabled` once per TU when #6 is set | B | — | (b) **O2b** | `config.hpp:352-355`: fires whenever `ARMA_CHECK_CONFORMANCE` was defined (always, `:180-181`) and `ARMA_WARN_LEVEL >= 2` (default 2, `:20-21`). Probed today with `-fsyntax-only`: note present at level 2, silent at level 1. GCC cannot suppress `#pragma message` by flag |
| 8 | Blaze size-conformance in release | B | always on (exceptions) | (b) **O2a** | not an assertion; `BLAZE_THROW` can be redefined (`Exception.h:93-95`) but a no-op throw lets the code run on with mismatched sizes — rejected (Appendix B) |
| 9 | Compile-time tripwire for the policy | B | absent | (a) | a test TU, **one branch per backend** (a TU is `BELFEM_ARMADILLO` or `BELFEM_BLAZE`, never both): under `BELFEM_ARMADILLO`, `#if BELFEM_ASSERTIONS_ACTIVE` → `#ifndef ARMA_CHECK_CONFORMANCE #error` (it is an **empty** define, `config.hpp:181`, so `#ifdef`, not `#if`), `#else` → `#ifdef ARMA_CHECK_CONFORMANCE #error`; under `BELFEM_BLAZE`, `#if BELFEM_ASSERTIONS_ACTIVE != BLAZE_USER_ASSERTION #error` (0/1 value). Guarded on the exported predicate (`assert.hpp:28`), never on a private `#ifndef NDEBUG` |
| 10 | Debug-mode tests for a shape mismatch and an out-of-bounds `()` | B | partly (`test_Matrix.cpp:446-466` covers wrapper `()` via `EXPECT_THROW(std::runtime_error)`, guarded `#ifndef NDEBUG` — pre-existing guard drift, out of scope) | (c) | the binary runs with `set_throw_on_error(true)` (`test_linalg_main.cpp:25-27`), so tests are `EXPECT_THROW`, not death tests. Shape mismatch must go through a path with no `BELFEM_ASSERT` (`operator+=`, `cl_AR_Matrix.hpp:513-517`): Armadillo throws `std::logic_error` (`debug.hpp:134`), Blaze `std::invalid_argument` (always). Backend bounds via `matrix_data()(i,j)`: Armadillo `std::out_of_range` (`debug.hpp:164`); Blaze asserts (no test). All guarded on `BELFEM_ASSERTIONS_ACTIVE` |
| 11 | Guide `:97-98` wrong switch; no "checks by build mode" table | B | wrong | (a) | `src/linalg/doc/linalg_usage_guide.md`; Codex language sweep owed (user-facing) |
| 12 | Armadillo `Matrix`: constructor from any expression | C | only `MatrixType` and `arma::Op<ET,OP>` (`:75`, `:80-81`) | (c) | add `template<typename ET> Matrix(const arma::Base<T,ET>&)`; probe that `arma::Base` is the right constraint for `eOp`, `eGlue`, `Glue`, `mtOp`, `subview` and that it does not ambiguate the `MatrixType` and copy ctors |
| 13 | Armadillo `Matrix`: `+=`/`-=`/`*=` from an expression | C | `MatrixType` only (`:514`, `:543`, `:572`) → temporary | (c) | add `template<typename ET> operator+=(const arma::Base<T,ET>&)` etc. — **constrained**, not the unconstrained `const ET&` the Armadillo `Vector` ships (`cl_AR_Vector.hpp:399-405`); probe no ambiguity with the `T` and `Matrix<T>` overloads. `operator=` is **not** a hole: `cl_AR_Matrix.hpp:476-482` already takes `const ET&` |
| 14 | Blaze `Vector`: `+=`/`-=` from a column expression | C | commented out (`:460`, `:511`); only `DMatScalarMultExpr<ET,T,false>` (`:468-470`, `:519-521`) | (c) | add `template<typename VT> operator+=(const blaze::DenseVector<VT,false>&)`; the generic `template<ET>` form is **refuted** (report round 2: `y += A.row(0)` fails to compile) |
| 15 | Blaze `Vector`: `+=`/`-=` from a row expression (`A.row(0)`, `trans(x.vector_data())`) | C | no direct overload; `y += A.row(0)` **compiles today** through the row ctor (`:127-136`) into a temporary `Vector` (probed) | (c) | add the transposing `template<typename VT> operator+=(const blaze::DenseVector<VT,true>&)` with a body of its own (§6) — the row **ctor** at `:127-136` is an element loop, not a `trans`, and is not a template for a `+=`. The `DMatScalarMultExpr<ET,T,false>` overloads (`:468-474`, `:519-525`) add `blaze::column( expr, 0 )` of a **row-major** matrix-scalar product (`StorageOrder.h:71`: `rowMajor = false`) and are **not** subsumed by a `DenseVector` overload; they stay unless R3.4's probe finds no caller, checked by value |
| 16 | Blaze `Matrix` `+=`/`-=`/`*=` from expressions | C | **GCC only**: `:596`, `:628`, `:669` and the plain `MatrixType` `-=`/`*=` sit inside `#ifdef BELFEM_GCC` (`cl_BZ_Matrix.hpp:593-601`, `:625-642`, `:666-683`); a `BELFEM_CLANG`/`BELFEM_INTEL` build converts through the ctor `:80-82` into a temporary | (c) | Grok round 3 P1. Lift the guards (R3.7). Probe 2026-11-07: Apple clang 17 `-fsyntax-only -DBELFEM_CLANG` accepts the unguarded set **with the bodies instantiated** (a planted error in the body is diagnosed; compile rung); the guard arrived in the squashed 2026-09-05 release commit, reason lost; Intel untested here |
| 17 | Blaze `Vector` construction from an expression | C | ctors `:121`, `:128` present | (c) | report §3b lists a construction case "untraced" on Blaze — trace it with the counter before declaring it closed |
| 18 | Executed allocation counter | C gate | absent | (c) | new `tests/linalg/test_LinalgAllocations.cpp`. **Both backends allocate through `posix_memalign`** (`memory.hpp:72`, `Memory.h:94`), never `operator new`, so the hook is a strong `posix_memalign` definition in the test binary that increments a **`thread_local`, monotonic** event counter (frees not subtracted — a heap-size delta cannot see a temporary that is freed inside the expression; `thread_local` keeps an MPI progress thread out of the window). **Positive control:** a deliberate backend temporary must count ≥ 1 or the test fails as "instrumentation not working". Sizes: the temporary must exceed Armadillo's 16-element local storage (`arma_config.hpp:27`), so matrices 8×8 and 64×64, vectors 32 and 256; n = 3 is recorded as "no heap on either side" only. Oracle: `EXPECT_EQ(wrapper count, backend count)` per shape — with aliasing the backend itself allocates (`DynamicMatrix.h:1670-1672`, `Mat_meat.hpp:5531-5533`), so the oracle is parity, not zero. Label `fast` (`Add_Test.cmake:43-45`). The interposer is `extern "C"`, forwards to the real allocator via `dlsym(RTLD_NEXT)` with the requested alignment (`memory.hpp:69`), and is initialised race-free (the binary runs MPI, `test_linalg_main.cpp:23`). It measures the test TU's **own** instantiations, which bind to the executable's symbol under static and shared `libbelfem` alike. If interposition fails on a CI platform the gate is **incomplete there**, stated as such — no fallback to heap statistics. Parity is the oracle **except** for shapes whose backend behaviour differs by design (the row operand, #32): those are tested by value in #33 and the counter compares against the explicit backend form (`strans`) |
| 19 | Hatch inventory, 110 lines / 13 files | D | none | (c) | per site: operation, operand shapes, live or dead path, workaround or convenience; Appendix D of this plan is filled as the inventory is produced |
| 20 | BDF2–5 assembly tests | D (Timestep) | **absent** | (c) | `tests/fem/test_BdfTimestepMethod.cpp:25-27`: nothing exercises BDF2–5 end to end; add element-level tests of `bdf2`–`bdf5`, `bdf*_nok`, `collect_qhist` against the hand formula before any line in `cl_IWG_Timestep.cpp` changes |
| 21 | Crank–Nicolson / Galerkin hatch lines | D | dead behind `BELFEM_ERROR` | (a) | `cl_IWG_Timestep.cpp:858`, `:872`; the rejection at `:84-92`, `:158-166` is tested at `test_BdfTimestepMethod.cpp:313` — leave as they are |
| 22 | `cl_Cohomology.cpp` site | D | — | closed | report to Gregory (protocol §7.1) |
| 23 | `nonfree` site (1 line) | D | — | (c) | `nonfree/` session, `nonfree/devlog/`, no exchange |
| 24 | `MatrixDataExposesBackend` test pins the hatch | D | present (`test_Matrix.cpp:192`) | (a) | stays: the hatch is a supported API (`commtools.hpp`, LAPACK glue), round D only removes its *misuse* |
| 25 | Documentation: "not a base class", check-policy table, overload-set contract | A–C | absent | (a) | `linalg_usage_guide.md`, `src/linalg/doc/README.md`; `doc/coding_philosophy.md` and `CLAUDE.md` only if they mention the virtuals (grep at close-out) |
| 26 | Debt-register row | approval | absent | (a) | one `[P]` row opened at plan approval, widened per round (policy 2026-08-29: extend, do not multiply) |
| 27 | Armadillo `Vector`: constructor from any expression | C | only `VectorType`, `arma::Op`, `subview_col/row` (`cl_AR_Vector.hpp:105-113`, `:120`, `:129`) | (c) | same hole as #12 (Grok round 3); same `arma::Base<T,ET>` fix |
| 28 | `belfem::trans( Vector )` does not exist | C tests | `fn_trans.hpp:72-88` unwraps `Matrix` only | (c) | round-C compatibility tests write `y += trans( x.vector_data() )`; adding a `Vector` overload is **O10**, not assumed |
| 29 | R0/R1 probe kernel must read through `const Matrix &` | 0, A | — | (a) | only `Matrix::operator() const` is virtual (`cl_AR_Matrix.hpp:276`, `cl_BZ_Matrix.hpp:309`); `Vector::operator() const` is not (`cl_AR_Vector.hpp:285`, `cl_BZ_Vector.hpp:340`) — a kernel holding non-const references or `Vector`s only will not show the §3 guard |
| 30 | Compilers, not only backends | all gates | GCC everywhere here | (c) | configured: `BELFEM_GCC` / `BELFEM_CLANG` (`config/compiler/config_gcc.cmake:4-8`) / `BELFEM_INTEL` (`config_icc.cmake:3`). Every tree on this Mac and the Linux CI are SCLS `g++`; Clang and Intel have **no** gate. Round C syntax-gates with `/usr/bin/clang++ -DBELFEM_CLANG` as well; Intel is stated untested |
| 31 | Tree names do not say the mode | all gates | — | (a) | `cmake-build-debug/` on this Mac is configured with `-DNDEBUG` (release; `CMakeFiles/belfem_lib.dir/flags.make`). Gates name the mode by `USE_DEBUG`, never by directory |
| 32 | **One orientation contract for `+=`/`-=` of a row operand, both backends** | C | Blaze: temporary via the row ctor (#15). Armadillo: `template<ET> operator+=(const ET&)` (`cl_AR_Vector.hpp:399-405`) forwards a 1×n operand into an n×1 `arma::Mat`: with checks on, `arma_assert_same_size` (`subview_meat.hpp:1720`, `debug.hpp:648-658`) → `std::logic_error`; with round B's checks off, a **row view** takes the `n_rows == 1` fast path (`subview_meat.hpp:1727-1749`) and is a silent in-bounds **transpose-add**, other mismatches are unchecked (`:1751-1756`). The ctor `:129-131` transposes (with `arma::trans`, see #37) | (c) | Codex round 1 P1, refined by Grok round 2. **Routing fact (probed by run):** a `Base<T,ET>`-constrained template never beats the existing unconstrained `const ET&` template, for any Armadillo type — so routing can only be done by (a) a **non-template** `subview_row<T>` overload (wins over both templates) using `arma::strans`, and (b) for 1×n *expressions* (`trans(x)*s`, `A.row(0)*s`) either a runtime shape branch inside the generic `+=`/`-=` (one compare per call) or "unsupported: debug throws, release adds elementwise" — **O12**. Blaze decides at compile time (`TF`). No consumer does `+= …row(` today (grep); #33 makes the contract real |
| 33 | **Numerical gate for every changed overload** | C | absent: `test_Vector.cpp:348` (`PlusEqualsScalar`) and `test_Matrix.cpp:336` (`RowViewReadsCorrectly`) cover scalar/wrapper operands and view reads only | (c) | Codex round 1 P1. Per changed overload: value **and** shape assertions, operands initialized independently, expected result computed before mutation, cases: expression, row view, scaled row, transposed column, column view, aliasing (`A += A*2`), matrix-expression `*=`, copy-init and direct-init construction (#12), `-=` twins of everything, **complex** operands with nonzero imaginary parts (#37), and the `DMatScalarMultExpr` column-0 overload's value if it stays, on both backends. Round D's consumer tests are too late to establish round C |
| 34 | Build-system registration of the new TUs | C | `tests/linalg/CMakeLists.txt:6-11` is an explicit `SOURCES` list; `Add_Test.cmake:33` builds `test_<name>_main.cpp` + `SOURCES` | (a) | add the files to `SOURCES`; the `posix_memalign` interposer needs no link flag (strong symbol in the executable; `--wrap` is **not** used, it is GNU-ld only); confirm in the gate that `make check` lists and runs the new tests |
| 35 | **Blaze runs OpenMP-parallel assignments under BELFEM's flags** | C gate, R0 | `USE_OPENMP=ON` → `-fopenmp` (`config_gcc.cmake:46-48`; present in `cmake-build-debug/…/flags.make`) → `_OPENMP` → `BLAZE_OPENMP_PARALLEL_MODE 1` (`blaze/system/SMP.h:65-66`, `config/SMP.h:58`) | (c) | Found while settling my round-2 pre-registration. Effect on this plan: the counter's `thread_local` count is complete only below Blaze's SMP thresholds — vectors < 38000 elements, dense mat-vec < 330 rows (`config/Thresholds.h:461-726`); the matrix sizes must stay below `BLAZE_SMP_DMATDMATMULT_THRESHOLD` (3025 target elements, `Thresholds.h:1632`) — 64×64 is not (D20: nondeterministic counts under the full `ctest` run on 16 cores); the test uses 8×8 and 48×48, vectors 32/256, stated in the test. R0 timings on large global vectors are OpenMP-affected and need `OMP_NUM_THREADS=1` for a like-for-like. **Side finding, out of scope:** every `Vector<real>` operation above those thresholds already runs on `OMP_NUM_THREADS` threads under Blaze — `doc/parallel_execution.md`/CLAUDE.md say BELFEM's own parallelism is MPI only; logged as **O11** for Christian, not changed here |
| 36 | **Discarded backend status becomes an unchecked index after round B** | B prerequisite | `fn_AR_eigen.hpp:33` ignores `arma::eig_gen`'s `bool`; on nonfinite input or `geev` failure the eigenvalues are `soft_reset()` (`fn_eig_gen.hpp:81-83`, `auxlib_meat.hpp:1098`, `:1131`) and `:42` indexes `tValues(k)` — today `std::out_of_range` (`Mat_meat.hpp:6901`), after round B an out-of-bounds read. `fn_BZ_eigen.hpp:34` relies on `blaze::eigen` throwing | (c) | Codex round 2 P1. **R2.0**: audit every `fn_AR_*`/`fn_BZ_*` for a discarded status before any check is removed; failure is an expected algorithmic outcome → status or `BELFEM_ERROR`, never an index into an emptied result (`coding_philosophy.md:643`) |
| 37 | **Complex transpose semantics differ between the backends today** | C | `arma::trans` is the conjugate transpose (`op_htrans`, `op_htrans_meat.hpp:277`); Blaze `trans` copies values (`DVecTransExpr.h:182-184`). `belfem::trans` already uses `arma::strans` (`src/linalg/fn_trans.hpp:28-33`), but the Armadillo row ctor `cl_AR_Vector.hpp:129-131` uses `arma::trans`: a `Vector<cplx>` built from a row is conjugated on Armadillo, not on Blaze (`cl_BZ_Vector.hpp:127-136`). Three `cplx` wrapper uses outside `src/linalg` | (c) | Codex round 2 P1. Round C: `strans` in every new overload **and** in the ctor; complex tests with nonzero imaginary parts (`test_lapack.cpp:246` already instantiates complex types) |

### 3.1 Cross-cutting findings

- **ABI, not just API.** Round A changes the object layout; every tree and every dlopen'ed plugin
  must be rebuilt from the same headers. The loader has no version handshake. Mixed builds in the
  shared `cmake-build-debug/` and `build/` trees are the one way round A can corrupt a run silently;
  the gate is a clean rebuild, not an incremental one.
- **Two backends, two machines.** No single workstation runs both; a round is not closed until the
  Linux CI matrix is green on both backends in both modes. Round B changes **release** behaviour
  only, so its gate must include a release tree on Armadillo — the debug-default workstation tree
  proves nothing about it.
- **`-fsyntax-only` is not a build.** It skips the unused-function and codegen-time diagnostics that `-Werror` turns fatal (D10); the
  pre-jury gate compiles each touched TU with `-c` under the build's own `flags.make`.
- **Static counts are not evidence.** Call-site counts in assembly were mis-read twice in the
  report's history (tail-call `jmp _free`, ≤16-element local storage). Round C's gate is the
  executed counter (#18), and round D accepts a replacement only when the counter shows the wrapper
  form at the backend form's count.
- **Order of rounds is load-bearing.** Round D's replacements rely on round C's overloads; round C's
  counter must run under round B's final check policy (a conformance check does not allocate, but
  the gate should be measured on the configuration that ships).

## 4. Ordered Steps

Each round: implement → syntax gate → **code jury** (same slug, Codex `gpt-6-astra`, Grok
`grok-4.7`; `high` for a first pass, `xhigh` for a repeat) → verification pass + reconciliation →
⛔ Christian runs the gate → close-out. Fixes raised by a jury are applied in a separate approved
step, never inside the round that raised them.

### Round 0 — baseline (optional, Christian, no change to the tree)

- [ ] **R0** In an isolated build tree (mode named by `USE_DEBUG`, gap #31), time one `examples/` deck
  whose assembly path reads element matrices through `const Matrix &` (inspect first: the transient
  Timestep path goes through `matrix_data()` and never calls the virtual accessor, so a transient deck
  can show nothing) as-is and with the round-A header edit, both backends, release. Record wall time and `make check-fast`. Falsifier: < 2 % difference means round
  A is hygiene, not performance — it still lands. Does not gate anything. The kernel must read
  `Matrix` elements through `const Matrix &` (gap #29) or it measures nothing.

### Round A — de-virtualize and seal (R1; `src/linalg` only)

- [x] **R1.1** Remove the six `virtual` specifiers (gap #1, #2). No signature changes. *(edited 2026-11-07; awaiting code jury + gate)*
- [x] **R1.2** Add `final` to the four class heads (gap #3; O4/O5 approved). *(edited 2026-11-07)*
- [x] **R1.3** Tripwires in `tests/linalg/test_Matrix.cpp` and `test_Vector.cpp` (gap #4, #5) *(edited 2026-11-07, via the public `MatrixType`/`VectorType` typedefs)*:
  `static_assert(!std::is_polymorphic<Matrix<real>>::value)` and
  `static_assert(sizeof(Matrix<real>) == sizeof(decltype(std::declval<Matrix<real>>().matrix_data())))`,
  same for `Vector`, plus one integer instantiation each (`Matrix<int>`, `Vector<index_t>`).
- [x] **R1.4** (report §3 probe re-run against the tree headers 2026-11-07; not checked in) Compile-only probe TU for the codegen claim: the §3 kernel of the report, built with
  the real flags; the acceptance is "no `call *` in the `const Matrix &` read loop" on both backends.
  Kept under `tests/linalg/probes/` as a documented manual probe, **not** in the test binary (an
  assembly grep is not a portable unit test).
- [x] **R1.5** Guide: `Matrix`/`Vector` are not base classes; `final` is a contract. *(edited 2026-11-07, Key Features item 6)*
- **Gate (Christian):** clean rebuild of both backends; `make check` (release) on Blaze/Apple and
  on the Linux CI matrix; rebuild of every plugin in use; R1.4 probe output attached to the devlog.
- **What could break:** nothing in-tree derives, casts or inspects polymorphism (grep today, gap
  scope guard). A `delete` through a base pointer cannot exist without a base. Out-of-tree
  derivation becomes a compile error at the first rebuild — the intended, approved break.

### Round B — check policy: debug on, release off (R2; `src/linalg/armadillo/armadillo.hpp`, docs)

- [x] **R2.0** (prerequisite, gap #36) Audit every `fn_AR_*.hpp`/`fn_BZ_*.hpp` for a backend status that is
  discarded and then relied upon *(inventory 2026-11-07, Appendix E; E1 `BELFEM_ERROR` on `eig_gen`, E2 `cross` bodies, E3 `BELFEM_ASSERT`s in `fn_dot`/`fn_det`/`fn_inv` applied the same day)*; fix each (`fn_AR_eigen.hpp:33`: check `arma::eig_gen`'s return, report
  through the existing `int_t` result or `BELFEM_ERROR`, never index an emptied vector) with a test that
  feeds a nonfinite matrix. Lands before R2.1 in the same round, both backends.
- [x] **R2.1** In `armadillo.hpp`, before `#include <armadillo>`: the tree's debug predicate
  (`#if !defined(NDEBUG) || defined(DEBUG)` — the same as `blaze_config.hpp:20` and `assert.hpp:26`),
  and on its `#else` side `#define ARMA_DONT_CHECK_CONFORMANCE` unless already defined (gap #6).
  `DEBUG` is defined under `USE_DEBUG` (`config_gcc.cmake:24`), so the two macros are not complements
  and the predicate, not `NDEBUG` alone, is the definition of a debug check in this tree.
- [x] **R2.2** The per-TU note (gap #7): `ARMA_WARN_LEVEL 1` on the release side of the same block
  (O2b → (ii), decided). The guide table of R2.6 states what is lost.
- [x] **R2.3** Blaze: **no change** (gap #8, O2a); the asymmetry — Blaze keeps size-conformance as
  an always-on exception — is documented, not engineered away (Appendix B).
- [x] **R2.4** Compile-time policy tripwire TU in `tests/linalg` (gap #9), guarded on
  `BELFEM_ASSERTIONS_ACTIVE`, both configurations.
- [x] **R2.5** Debug-mode `EXPECT_THROW` tests for one shape mismatch (through `operator+=`, expecting
  `std::logic_error` on Armadillo and `std::invalid_argument` on Blaze) and one backend out-of-bounds
  `matrix_data()(i,j)` (Armadillo `std::out_of_range`; Blaze has no throwing bounds check) per gap #10,
  guarded on `BELFEM_ASSERTIONS_ACTIVE`, registered per gap #34. Plus one **unguarded** test that runs in
  release too: Blaze shape mismatch throws `std::invalid_argument` in every mode (that is O2a), and the
  Armadillo twin is skipped in release (it would be unchecked).
- [x] **R2.6** *(sweep still owed)* Guide fix at `:23`, `:94` and `:97-98` (all three state that Armadillo checks bounds in
  release) and a "checks by build mode × backend" table (gap #11); Codex language sweep
  (`gpt-5.6-luna`/medium) on the touched section only.
- **Gate (Christian):** `make check` in a **release** Armadillo tree and a debug Armadillo tree
  (Linux), plus the Blaze pair; R2.4 compiles in all four; the Armadillo release build log is
  inspected for the O2b outcome (note present/absent as decided). Timing is **not** a gate here
  (report §3: 1.00× on the 3×3 kernel); the deliverable is the policy.
- **What could break:** any path that today relies on an Armadillo exception as its error handling —
  R2.0 is the inventory of exactly those. After it, a release run that today dies with an Armadillo
  `logic_error` on a shape mismatch will instead compute on garbage or segfault, exactly as it does on Blaze-bounds and on
  every `BELFEM_ASSERT` today. That is the decided policy; the debug tests of R2.5 are the net.

### Round C — complete the overload sets (R3; `src/linalg` class headers + one new test)

- [x] **R3.0** (prerequisite, gap #20, moved here from R4.1 by Codex round 2) Element-level assembly tests
  for `bdf2`–`bdf5`, `bdf*_nok`, `collect_qhist`, `derivative*`, `static_subfield` against the hand
  formula, on both backends — the baseline that round C's rebinding is measured against.
- [x] **R3.1** Inventory matrix (shape × class × backend): which overload each expression reaches
  and whether a temporary is constructed — from source trace first, then confirmed by R3.2's counter.
  Shapes: `C = B*dt`, `Matrix C(B*dt)`, `A += B*dt`, `A -= B*dt`, `A *= s`, `y = x*dt`,
  `Vector y(x*dt)`, `y += x*dt`, `y -= x*dt`, `y += A.row(0)`, `y += A.col(0)`, `y += A.col(k) / s`
  (`cl_ThinShellFactory.cpp:935`), `y += trans(x.vector_data())` (gap #28), `y += A*x`,
  `y += ( A - B*dt ) * q` (mixed operands, `cl_IWG_Timestep.cpp:852`), `y -= K*q` (`:1129`),
  `y += A.row(0) * s`, `y += trans(x.vector_data())` with complex `T`, `A += trans(B)`, `A *= B*dt`,
  each with wrapper and with raw-backend operands, aliased and not. Copy-initialization
  `Matrix<real> C = B*dt` is **rejected today** (two user-defined conversions; probed) while
  direct-initialization `Matrix<real> C( B*dt )` goes through the temporary — after R3.3 both compile
  and both are in #33.
- [x] **R3.2** The allocation-counter test (gap #18) **and** the numerical tests (gap #33), written
  **first** so the counter fails on the gaps before R3.3 closes them (the test is the reproducer) and
  the numerical tests pin the values the new overloads must reproduce. Registered per gap #34.
- [x] **R3.3** Armadillo `Matrix`: constrained expression constructor and `+=`/`-=`/`*=` (gaps #12,
  #13), and the same constructor for the Armadillo `Vector` (gap #27). Constraint `const arma::Base<T, ET> &`; probe ambiguity against `MatrixType`, `Matrix<T>`,
  `T`, `arma::Op<ET,OP>` (`:80-81` may become redundant — remove only if the counter and the
  compile of `tests/linalg` say so).
- [x] **R3.4** Blaze `Vector`: `+=`/`-=` for `blaze::DenseVector<VT,false>` and the transposing
  `DenseVector<VT,true>` (gaps #14, #15); the Armadillo `Vector` gets the **non-template**
  `subview_row<T>` overloads with `arma::strans`, the row ctor `:129-131` switches to `strans` (gap #37),
  and the 1×n-expression case follows O12 (gap #32); the `DMatScalarMultExpr` overloads
  (`:468-474`, `:519-525`) stay unless a caller search and a by-value probe show them unused. Delete the commented-out generic forms
  (`:441`, `:460`, `:501`, `:511`) — commenting guideline: no dead code.
- [x] **R3.5** Trace the Blaze construction case (gap #17) and the Blaze `Matrix` set (gap #16) with
  the counter; close or open rows accordingly.
- [x] **R3.7** Lift the three `#ifdef BELFEM_GCC` guards in `cl_BZ_Matrix.hpp` (gap #16) so the
  overload set is the same on every compiler. Syntax gate with SCLS `g++` **and** Apple clang, the
  exact command being
  `/usr/bin/clang++ -std=gnu++17 -fsyntax-only -DBELFEM_BLAZE -DBELFEM_CLANG -DBELFEM_MPI -DOMP -DBELFEM_NETLIB -DNDEBUG -I src/core -I src/comm -I src/containers -I src/linalg -I src/linalg/lapack -I src/linalg/operators -I src/linalg/blaze -I <build>/generated -isystem /opt/scls/include <tu>`
  on `tests/linalg/*.cpp` (gap #30). Clang is **compile-gated only**: linking and running the test
  binary with Clang against GCC-built SCLS libraries is not attempted. Also one R3.1 row per compiler for
  the ranking of the `MatrixType` non-template overloads against the `blaze::Matrix<MT,SO>` templates. If a
  supported compiler rejects the unguarded set, the fallback is a per-compiler path, recorded here.
- [x] **R3.8** Retire `dot( Matrix, Vector )` and `dot( Vector, Matrix )` on both backends with their `fn_dot.hpp`
  doc (O14); a compile-fails test pins their absence.
- [x] **R3.6** Guide: the overload-set contract — "for the supported expression set, the wrapper adds
  **no materialization beyond the backend's own** for ctor, `+=`, `-=`, `*=`" (`=` already forwards);
  not "temporary-free": Armadillo's plain in-place product `*=` **always** materializes
  (`glue_times_meat.hpp:586-597`), and `+=` of a product does so on alias or outer product (`:612`).
  The counter test is the proof of parity.
- **Gate (Christian):** R3.2 green on both backends, debug and release (counts differ between modes
  only if a check allocates — it must not); `make check`; `tests/linalg` label `fast` keeps
  `check-fast` under its budget.
- **What could break:** overload ambiguity is a compile error and surfaces in `tests/linalg`; a
  wrong transposing overload is a wrong result and surfaces only in the new tests of gap #33 (the
  existing ones do not reach these overloads). **Every existing call site rebinds** to the new
  overloads — including the four `aJ += K().matrix_data()*dt` lines already in `cl_IWG_Timestep.cpp`
  — so this is a behaviour change across the tree, gated by `make check` on both backends, not a
  test-only change.

### Round D — retire the escape hatches in consumers (R4; after round C)

Sub-rounds, each with its own code jury:

- [ ] **R4.0** Inventory (gap #19) → Appendix D. Classify each of the 110 lines: live/dead,
  workaround/convenience, replacement expression, test that covers the path. Lines that are
  *legitimate* uses of the hatch (handing the backend object to LAPACK or MPI glue,
  `commtools.hpp`, `hdf5_tools.hpp`) are marked "keep" and are not touched.
- [ ] **R4.1 — `cl_IWG_Timestep.cpp` (79 lines).** First the BDF2–5 tests (gap #20); then, per
  method, replace `aJ.matrix_data() += M().matrix_data()` by `aJ += M()`, `aSwap.vector_data() =
  mBeta(0) * qold(0).vector_data()` by `aSwap = mBeta(0) * qold(0)`, etc. Crank–Nicolson and
  Galerkin lines stay (gap #21). Gate: the new tests, the R3.2 counter run on each replaced shape,
  `make check`, and one transient example deck run before/after agreeing to a **stated tolerance**
  (relative 1e-14 on J and RHS; bitwise identity recorded when it occurs, not required — a fused
  expression may round differently). The class has 20 naming files, not one: every IWG deriving from
  it (`cl_IWG_Maxwell`, `cl_IWG_TransientHeatConduction`, `cl_IWG_Poisson`, `cl_IWG_StaticHeatConduction`),
  the `mt_*` matrix factories, `cl_FEM_Controller`, `cl_FEM_DofMgr_EigenValues`, `test_QoldTable.cpp`
  (derives `TwoFieldIwg` from it) — their tests are in the gate.
- [ ] **R4.2 — `src/fem/kernel` + `src/fem/postproc`** (`cl_FEM_Calculator.hpp` 9,
  `cl_FEM_DofMgr_EigenValues.cpp` 1, `cl_Gradient.cpp` 1, `fn_Mesh_compute_edge_lengths.cpp` 1, …).
  The Calculator's nine lines are three assignments `mJ->matrix().matrix_data() = …` (`:1724`,
  `:1740`, `:1755`) and six `dot( phi().vector_data(), mX.col(·) )` calls (`:2233-2298`, mixed
  wrapper/backend operands of `dot`); the former become `template<ET> operator=` through the wrapper,
  the latter need a `dot` overload inventory first — both counter-gated. The Crank–Nicolson and
  Galerkin hatch statements are `cl_IWG_Timestep.cpp:862-866` and `:876-881` (the earlier `:858`/`:872`
  were the signatures). `cl_Gradient.cpp:196` keeps
  the eager `inv()` (`fn_BZ_inv.hpp:24`, deliberate).
- [ ] **R4.3 — `src/physics`, `src/numerics`, `src/mesh`** (`cl_Gas.cpp`, `cl_GM_EoS_Cubic.cpp`,
  `cl_Material_UserDefined.*`, `fn_ODE_RK45.cpp`, `fn_ODE_DOP853.cpp`, …). Gas-model tests are in
  the `check` fan-out under `-DUSE_GASMODELS=ON` (CI since 2026-09-18).
- [ ] **R4.4** `nonfree` site — separate `nonfree/` session (gap #23).
- [ ] **R4.5** `cl_Cohomology.cpp` — report to Gregory (gap #22).
- **Gate per sub-round (Christian):** named numerical test of the path + R3.2 counter + `make check`
  both backends.
- **What could break:** a replaced expression that silently picks a *different* overload (e.g. the
  Blaze `DMatScalarMultExpr` path) and changes semantics — the counter catches the allocation, the
  numerical test catches the value. Replacements in `collect_qhist` change nothing in the math, but
  the BDF2–5 paths have **no** end-to-end test today, which is why R4.1 starts with the tests.

### Round E — close-out

- [ ] **R5.1** Devlog per round (`devlog/dlYYYYMMDD_linalg_wrapper_round_<X>.md`) + index line.
- [ ] **R5.2** Debt-register row updated per round (gap #26); plan Status line refreshed; boxes ticked
  in the same session.
- [ ] **R5.3** Guides swept (Codex language pass, user-facing sections only); `check_doc_claims.py`
  after any edit to `CLAUDE.md`/`coding_philosophy.md`; `check_wrapper_policy.py` after round D (the
  hatch retirement must not introduce a direct backend call outside `src/linalg`).
- [ ] **R5.4** Move the investigation report and this plan to `todo/closed/` with DONE summaries.

### 4.0 Implementation Progress (updated 2026-11-07, evening: round C probes)

**Round C — R3.2 counter design probed (2026-11-07, `counter/probe.cpp`, strong `extern "C" posix_memalign` + `dlsym(RTLD_NEXT)`, thread-local count,
positive control):**
- **Blaze (GCC, this Mac):** instrumentation works (positive control 1). Executed counts, wrapper / backend: `A += B*dt` 0/0, `A -= B*dt` 0/0,
  `C = B*dt + A` 0/0, `A += trans(B)` 0, `Matrix D( B*dt )` **1/1** (8×8 and 64×64); `y += x*dt` **1/0**, `y -= x*dt` 1, `y += x + z` 1,
  `y += A.col(0)` 1, `y += A*x` **1/0** (32 and 256). So **gap #17 closes** (Blaze construction has parity; the report's "2+2 vs 1+1" was a static
  count) and the Blaze `Vector` gaps #14/#15 are confirmed as executed allocations.
- **Armadillo on macOS:** the SDK defines `_POSIX_ADVISORY_INFO` as −1, so `ARMA_HAVE_POSIX_MEMALIGN` is never set here and `memory::acquire`
  falls back to `std::malloc` (`compiler_setup.hpp:106-109`, `memory.hpp:92-95`); the `posix_memalign` hook counts nothing (positive control 0),
  dyld `__interpose` sections are honoured only in dylibs, and a strong `malloc` with `dlsym` inside is unsafe. **Decision for #18:** the test hooks
  `posix_memalign`, runs its positive control first, and `GTEST_SKIP`s with a message on a platform where the control counts 0 — i.e. Armadillo
  on macOS is not a counted configuration; Armadillo on Linux (glibc, `posix_memalign`) and Blaze everywhere are. Stated in the test and here.
- **Armadillo compile fact for gap #13:** `A += trans( B )` does **not compile** today — ambiguous between `operator+=( const Matrix<T> & )`
  (via the `arma::Op` constructor `:80-81`) and `operator+=( const MatrixType & )` (via `Mat(const Op&)`). The gap is a compile error, not only a
  temporary; R3.3's `Base`-constrained overloads resolve it (exact template match beats both conversions).
- `belfem::trans` on Armadillo returns `Op<Mat, op_htrans>` for real types because `arma::strans` itself does; for complex it is the plain
  transpose (probed: `(5,6)`, same as Blaze). No finding.

**Approved 2026-11-07 (Christian):** plan v3, source edits, O12 (ii). Debt row **DR-200** opened.

**Round A — VERIFIED (both backends, Christian's gates) and committed with round B.** Four headers: `final` on the class head, `virtual` removed from the
destructor (all four) and from `operator()(i,j) const` (both `Matrix` classes); no other change. Tripwires
in `test_Matrix.cpp` / `test_Vector.cpp`: `!is_polymorphic` and `sizeof(wrapper) == sizeof(backend typedef)`
for `real` and one integer instantiation each. Guide: Key Features item 6. R1.4 (assembly probe) is the
report's §3 probe, re-run as part of the gate, not checked in. Syntax gate and code jury: see §8.

**Defect tracker (Dn)** — code jury round A, 2026-11-07 (Codex `gpt-6-astra`/high, Grok `grok-4.7`/high):
- [x] **D1** (LOW, Grok) — the tripwires pin `!is_polymorphic` and `sizeof`, not `final` itself; removing `final` would leave them green.
  Fix: add `static_assert( std::is_final< … >::value )` beside them (compiles on both backends, probed).
  **Fixed 2026-11-07** (approved by Christian; applied by Claude; syntax-gated, not yet built).
- [x] **D2** (LOW, Codex + Grok) — guide item 6 and the assert messages say "exactly its backend object"; the asserts establish size
  equality and non-polymorphism, not type identity. Fix: "adds no object-size overhead and has no vtable".
  **Fixed 2026-11-07** (approved by Christian; guide item 6 and both assert messages reworded).
- No P0/P1. Rebuild rule (every TU that includes the headers, every tree, plugins that reach into linalg) restated by both auditors.

Code jury round B, 2026-11-07 (Codex `gpt-6-astra`/high, Grok `grok-4.7`/high; A+B diff; all CONFIRMED, none refuted):
- [x] **D3** (LOW, Codex + Grok) — guide table's LAPACK row: `BELFEM_ERROR` only when `AbortOnError` (default); otherwise `info` is returned
  (`fn_gesv.hpp:214-216`). Fix: qualify the row.
- [x] **D4** (LOW, Codex) — `=` is not a shape error, backend assignment resizes (`Mat_meat.hpp:920`, `DynamicMatrix.h:4891`). Fix: drop `=` from the row.
- [x] **D5** (LOW, Grok) — element-bounds cell: debug Armadillo is `BELFEM_ASSERT` (it throws before the backend is reached); `std::out_of_range`
  only through `matrix_data()`. Fix: the cell.
- [x] **D6** (LOW, Grok) — "pins the table" → pins the switch and a sample of rows; "the one asymmetry" → also `ARMA_CHECK_NONFINITE` stays on in
  release Armadillo. Fix: the two sentences.
- [x] **D7** (**HIGH**, Grok) — `eigen`/`eigen_sym` squareness is `BELFEM_ASSERT` only; with Armadillo's own check off, a non-square release call
  reaches `geev` (`auxlib_meat.hpp:1092` is `arma_conform_check`). Fix: `BELFEM_ERROR` for squareness in all four `eigen*` implementations
  (setup-grade functions, one compare).
- [x] **D8** (LOW, Grok) — `cross( Vector<T> &, const ET & )` takes a non-const `Vector` (`fn_AR_cross.hpp:33`, `fn_BZ_cross.hpp:33`), so a const
  `Vector` falls through to the `(ET,ET)` template and fails; the two mixed overloads have no length-3 assert. Fix: `const Vector<T> &` and the
  assert, both backends.
- [x] **D9** (**HIGH**, Claude probe, Appendix E E7; Grok "untraced") — `blaze::eigen` on a non-finite matrix returns zero eigenvalues silently
  (`geev.h:139-142` throws for `info > 0` only). Fix: a finiteness pre-check with `BELFEM_ERROR` in both `eigen` implementations, so the backends
  fail the same way; the NaN test becomes backend-neutral.
- [x] **D10** (LOW, Christian's Armadillo release gate) — `test_CheckPolicy.cpp` helper `square()` unused on Armadillo release → `-Werror=unused-function`. Fixed 2026-11-07 (helper removed). **Gate lesson:** `-fsyntax-only` does not run the unused-function pass; the syntax gate now compiles with `-c` under the build's `flags.make` (gap #31 companion).
- **O14** (Grok confirms E5) — Blaze has no `dot(Matrix, Vector)`; the two documented overloads are dead there.
  **RESOLVED 2026-11-07 → retire both overloads and their doc in round C (Christian)**; step R3.8.
- D3–D9 **fixed 2026-11-07** (approved by Christian): `fn_is_finite.hpp` (new, backend-free), `BELFEM_ERROR` for squareness and finiteness in
  all four `eigen*` implementations, `const Vector<T> &` plus length asserts in the mixed `cross` overloads, guide rows D3–D6, tests
  extended (NaN and non-square `eigen`/`eigen_sym` in every build on both backends; const operands in the mixed `cross` test).

**Round B — VERIFIED 2026-11-07: Christian's `make check` on Blaze (release, VTK on) and Armadillo (after D10) green for `test_linalg`; committed.** Earlier same evening gated on this Mac 20:13 CST (Blaze, release): `test_linalg` 289 passed incl. the four every-build `CheckPolicy` cases; the 54 other failures are exactly the DR-202/203 set (HDF5Test 30, JcDatabaseTable 10, ElectricalCircuit 9, TensorMeshBSpline 3, BfmFile 1, AlloyRRR 1), none new. Armadillo and debug legs owed to CI.** `armadillo.hpp`: `ARMA_DONT_CHECK_CONFORMANCE` +
`ARMA_WARN_LEVEL 1` on the release side of `!defined(NDEBUG) || defined(DEBUG)` (the per-TU note is gone, probed). E1–E3
applied. `tests/linalg/test_CheckPolicy.cpp` (new, registered): compile-time tripwire per backend; every-build tests
(Blaze shape mismatch → `std::invalid_argument`; Armadillo failed `eigen` → `std::runtime_error`; mixed `cross` overloads by
value); debug-only tests (`dot`/`det`/`inv` asserts; Armadillo shape mismatch → `std::logic_error`, backend bounds →
`std::out_of_range`). Guide pitfall 5 rewritten with the checks-by-mode table (Codex language sweep owed). Syntax gate:
g++ Blaze/Armadillo × release/debug, Armadillo with `-DNDEBUG -DDEBUG` (behaves as debug, as the predicate says), Apple
clang Blaze — all clean. **R2.0 inventory (read-only), Appendix E.** Four findings: (E1) `fn_AR_eigen.hpp:33`
discards `arma::eig_gen`'s status — the one real hole (callers: `tests/linalg/test_LinalgSolvers.cpp` only, 14
lines, none in `src`/`nonfree`); (E2) `cross( ET, Vector )` has a body that calls `vector_data()` on the backend
operand in both backends (`fn_AR_cross.hpp:44-46`, `fn_BZ_cross.hpp:44-46`) — never instantiated, a latent
compile error; (E3) with round B, Armadillo stops checking `dot` lengths, `cross` lengths, `det`/`inv`
squareness in release, and BELFEM's own dispatch headers assert none of these except `cross`
(`fn_cross.hpp` 1 assert; `fn_dot.hpp`, `fn_det.hpp`, `fn_inv.hpp` 0) — Blaze throws on all of them in every
mode, so a `BELFEM_ASSERT` in each dispatch header keeps debug symmetric; (E4) every LAPACK wrapper checks
`info` (`fn_gees/geev/gels/gesv/gesvd/getrf/getri/posv.hpp`), nothing to do. Round B's diff therefore also
carries E1–E3.

**Round C — implemented 2026-11-07, code-juried (Codex `gpt-6-astra`/high + Grok `grok-4.7`/high), D11–D18 applied 2026-11-08 with
approval, gate not run (Christian).** Jury: two P1s — `v += 1` stopped compiling on Armadillo (the unconstrained template beat the scalar overload
and instantiated `arma::Proxy<int>`; D11 constrains on `is_arma_type`) and `Vector( trans(x) )` stayed 1×n (the old `arma::Op` constructor beat
the `Base` one; D13 deletes the `Op` constructors of both classes) — plus: my O14 slice had taken Blaze's generic `dot(ET1,ET2)` and `dot(Rows,ET)`
with it (D12, restored); a generic `Op` operand was evaluated twice by the shape branch (D14: branch only for `eOp`/`eGlue`/`mtOp`/`mtGlue`/transpose
Ops, everything else forwarded; D15: products, `arma::Row`, `submat` views are documented as not oriented, debug assert on a row result, trait
private and CamelCase); dead `tQhist` block (D16); `dlsym` null check (D17); guide `dot` paragraph (D18). Re-gated: seven linalg test TUs,
`test_BdfTimestepMethod.cpp`, twelve consumer TUs, both backends × both modes, clean; value probes pass on both backends; Blaze allocation
parity holds for every shape. Nothing built in the tree. First `make check` (Christian, Blaze release): 16/17 suites green (DR-202 has landed in that tree); `test_fem` segfaulted in the new
`BdfAssembly` fixture, which had not called `link_to_group( tBlock )` before assembling — the kernel's own loop does (`cl_FEM_DofManager.cpp:587`).
**D19 fixed 2026-11-08**, compile re-gated on both backends and modes; the run is Christian's. Second run: Armadillo green incl. `BdfAssembly`; Blaze `test_linalg` red on the 64×64 `*=` parity count, which is nondeterministic above Blaze's SMP product threshold under `-fopenmp` — **D20 fixed 2026-11-08** (48×48). Third run: `make check` green on Armadillo and on Blaze — **round C VERIFIED on both legs**, committed.
- R3.0: `tests/fem/test_BdfTimestepMethod.cpp` §7 — `SyntheticMkfIwg` (synthetic M, K, f) + `AssemblyProblem`; `Bdf1ToBdf5MatchTheHandFormula`
  and `RampRunsTheOrderTheHistorySupports` check `compute_jacobian_and_rhs` against `J = alpha M + h K`, `b = h f + M Σ(−1)^k beta_k q_k` with
  the IWG's own coefficients (§4 pins those). Compiles in both modes.
- R3.2/#18: `tests/linalg/test_LinalgAllocations.cpp` — strong `extern "C" posix_memalign` forwarding via `dlsym(RTLD_NEXT)`, thread-local event
  count, positive control → `GTEST_SKIP` where the hook is blind (Armadillo on macOS: `std::malloc`). Shapes: matrix `+=`/`-=`/`*=`/`=`/ctor/copy-init/
  aliased, vector `+=`/`-=` scaled, sum, column, mat-vec, ctor, and the row operand against the explicit transposing backend form. Probe with the
  new headers on Blaze: every former `1/0` is now `0/0`.
- #33: `tests/linalg/test_OverloadSets.cpp` — values and shapes for every new overload, both backends, incl. aliasing, `*=` product, copy-init,
  integer instantiation, complex row (not conjugated), transposed column, and the Blaze `DMatScalarMultExpr` column-0 semantics. Run probe of the
  same cases on both backends: all pass.
- R3.3/#12/#13: `cl_AR_Matrix.hpp` — `template<ET> Matrix( const arma::Base<T,ET>& )` and `+=`/`-=`/`*=` for `const arma::Base<T,ET>&`
  (`A += trans(B)` was a compile error before: ambiguous between the `Matrix` and `MatrixType` conversions).
- #27/#32/#37 (O12 ii): `cl_AR_Vector.hpp` — `Vector( const arma::Base<T,ET>& )`, non-template `subview_row` overloads for `=`/`+=`/`-=` with
  `arma::strans`, the row ctor on `strans`, and `assign/plus/minus_oriented`: products, `Mat` and views go to the backend unchanged; every other
  expression asks a `Proxy` for its shape (lazy for element-wise expressions and `trans`) and is wrapped in `strans` when it is 1×n.
- R3.4/#14/#15: `cl_BZ_Vector.hpp` — `+=`/`-=` for `DenseVector<VT,false>` (as is) and `DenseVector<VT,true>` (`blaze::trans`); the generic
  `operator=` transposes a row expression (`if constexpr IsRowVector_v`); the four commented-out blocks are gone; the `DMatScalarMultExpr`
  overloads stay (pinned by test).
- R3.7/#16: the three `#ifdef BELFEM_GCC` guards in `cl_BZ_Matrix.hpp` lifted; the expression `*=` now forms the product in a temporary and moves it,
  as Blaze's own operator does (the old body did not compile for an expression operand).
- R3.8/O14: `dot( Matrix, Vector )` and `dot( Vector, Matrix )` retired on both backends with their `fn_dot.hpp` doc; the Armadillo mixed `dot`
  takes `const Vector<T>&` now (it took a non-const reference, the D8 wart again).
- R3.5/#17: closed — the Blaze construction case has parity (1/1 at every size, executed).
- R3.6: guide Key Features item 6 (expression operands, orientation contract, no matrix `dot`). Codex sweep owed with pitfall 5.
- Gate so far: `-c` under the build's `flags.make`, both backends × both modes: the five linalg test TUs, `test_BdfTimestepMethod.cpp`, and the
  consumer TUs `cl_IWG_Timestep.cpp`, `cl_FEM_Calculator.cpp`, `mt_maxwell_phi.cpp`, `cl_ThinShellFactory.cpp`, `cl_Gradient.cpp`,
  `fn_IF_initialize_integration_points_on_facet.cpp`, `cl_TensorMeshFactory.cpp`, `cl_Maxwell_TMatrix.cpp`, `fn_ODE_DOP853.cpp` — all clean.
  Nothing built or run beyond the two scratchpad probes.

**Still unimplemented:** round D (R4.0 inventory next), round E; R0 (optional); the Codex sweep of the guide's user-facing changes (items 3, 6, pitfall 5, `dot`).

## 5. Open Design Questions (logged, not silently decided)

- **O1 — both backends kept?** **RESOLVED 2026-11-07 → yes** (Christian: "this has to work for both
  Armadillo and Blaze"). Consequence: design B (inheritance) stays struck.
- **O2 — conformance checks in release?** **RESOLVED 2026-11-07 → on in debug, off in release**
  (Christian). Two residuals:
  - **O2a — Blaze keeps its always-on size exceptions.** Proposal: accept. The check is one
    size comparison per *expression*, not per element; the only way to remove it is to redefine
    `BLAZE_THROW` (`Exception.h:93-95`), and a no-op there lets `A += B` run on with mismatched sizes
    (undefined behaviour instead of an exception). Cost of accepting: a release Blaze run still
    dies cleanly on a shape mismatch while a release Armadillo run will not — an asymmetry the
    guide table (R2.6) states. **RESOLVED 2026-11-07 → accepted (Christian).**
  - **O2b — the per-TU `#pragma message` note on Armadillo release builds.** Options:
    (i) accept the note (one `note:` per translation unit, ~1000 in a full Linux release build log,
    harmless, drowns real warnings); (ii) `#define ARMA_WARN_LEVEL 1` in release — silences the note (probed) **and** Armadillo's
    level-2 runtime warnings (`config.hpp:22-28`: poorly conditioned systems in `solve()`/`spsolve()`,
    approximate-solution fallbacks). These are **numerical diagnostics**, not invariant checks, so the
    release-assertion policy does not settle them (Codex round 1). What BELFEM can emit at level 2:
    the `arma::inv`, `arma::det`, `arma::eig_*` calls in `fn_AR_inv.hpp`, `fn_AR_det.hpp`,
    `fn_AR_eigen.hpp`; BELFEM calls no `arma::solve` (grep). Debug builds keep level 2 either way.
    **Recommendation: (ii)**, on the narrower ground that a release build log with a thousand identical
    notes hides real warnings and the three emitters are wrapped calls whose failures BELFEM already
    checks itself (`fn_AR_inv.hpp`). **RESOLVED 2026-11-07 → (ii) `ARMA_WARN_LEVEL 1` in release
    (Christian)**; the §6 block is now unconditional.
- **O3 — fixed-size small matrices** (report O3). **Investigated 2026-11-07 (Christian's request) →
  `todo/linalg_small_matrix_kernels.md`**: not a wrapper item; Armadillo's BLAS dispatch for non-square
  small products is the cost, an Armadillo-side kernel behind the backend split the fix (Blaze's own kernel
  wins at hex27). Plan-juried 2026-11-07 (v2). Sequenced after round A and its R1 profile; independent of
  B–D (round B adds ~25 % to the kernel's gain, it does not gate it).
- **O4 — `final` on both classes.** Grounds: it is the only mechanism that makes downstream derivation
  a compile error; nothing derives today. **RESOLVED 2026-11-07 → yes, in round A (Christian).**
- **O5 — the `final` API break is a deliberate decision**, recorded with the plan approval, bounded
  by the backend-free plugin contract. The layout change is wider than derivation: any out-of-tree TU
  holding a `Matrix`/`Vector` by value must be rebuilt, derived or not (Grok round 3). **RESOLVED
  2026-11-07 → the break is accepted (Christian); recorded here as the compatibility decision.**
- **O6 — jury granularity for round D.** Proposal: one code jury per sub-round (R4.1, R4.2, R4.3),
  since R4.1 alone is 79 lines plus new tests. Alternative: one jury for all of D — cheaper,
  weaker. **Recommendation: per sub-round.**
- **O7 — when R0 runs.** Proposal: in parallel with round A's jury, in an isolated tree; it never
  gates. Christian decides whether to run it at all.
- **O8 — delete `arma::Op<ET,OP>` ctor (`cl_AR_Matrix.hpp:80-81`) once the `arma::Base` ctor
  exists?** Decided in R3.3 by probe + counter, logged here when decided.
- **O10 — add `belfem::trans( const Vector<T> & )`?** `fn_trans.hpp:72-88` has `Matrix` only, so
  `trans(x)` on a wrapper vector is ill-formed today (gap #28). Not needed by rounds A–C; decided when
  round D's inventory shows whether any consumer wants it. Default: no.
- **O12 — Armadillo row *expressions* in `Vector::operator+=`/`-=` (gap #32).** Row **views**
  (`A.row(i)`) get a non-template overload and transpose on both backends. For a 1×n *expression*
  (`trans(x)*s`, `A.row(0)*s`) Armadillo cannot select by type (probed: the unconstrained template wins
  over any `Base` overload). Options: (i) **unsupported** — documented; debug throws `std::logic_error`
  as today, release adds elementwise (equal length) or runs unchecked; (ii) a runtime shape branch in the
  generic `+=`/`-=` (`if( n_rows == 1 && n_cols > 1 ) mVector += strans(expr) else …`), one compare per
  call, symmetric with Blaze. **Recommendation: (ii)** — one integer compare buys backend symmetry on a
  path that is per-element-matrix, not per-element. **RESOLVED 2026-11-07 → (ii) runtime shape branch
  (Christian).**
- **O13 — the new Armadillo expression constructor is implicit.** Copy-initialization
  `Matrix<real> C = B*dt` becomes legal on Armadillo; it is already legal on Blaze
  (`cl_BZ_Matrix.hpp:80-82` is a non-explicit template — expected, probed in R3.3). **Recommendation:
  implicit, for backend parity**; the existing `Matrix(const MatrixType&)` (`:75`) is implicit already.
  Logged, default taken unless Christian objects.
- **O14 — the flat `dot( Matrix, Vector )` overloads exist on Armadillo only** (Appendix E, E5). **RESOLVED
  2026-11-07 → retire on both backends in round C, with the doc (Christian).**
- **O11 — Blaze OpenMP mode (gap #35).** Blaze parallelizes large assignments with OpenMP whenever
  `-fopenmp` is on, which BELFEM's `USE_OPENMP=ON` guarantees. Keep (free speed on large dense vectors,
  thread count already pinned by `Allrun`'s `OMP_NUM_THREADS`), or set `BLAZE_USE_SHARED_MEMORY_PARALLELIZATION 0`
  in `blaze_config.hpp` to make the "MPI only" statement literally true? Out of this plan's scope; the
  documentation statement is the thing that is wrong today. **Needs Christian.**
- **O9 — review-script repairs** (`ask_copilot.sh` under bash 3.2, `cross_review.sh` reporting a
  Grok auth failure only as "Audit FAILED"). Out of scope; their own tooling round. Until then the
  scripts are run with `/opt/bash/bin/bash` (bash 5.3).

## 6. Interface Design — the exact edits, for audit

### Round A (identical pattern in all four headers)

```cpp
// cl_AR_Matrix.hpp:27-28  (and cl_AR_Vector.hpp:31-32, cl_BZ_Matrix.hpp:26-27, cl_BZ_Vector.hpp:39-40)
template< typename T >
class Matrix final
{
...
// :106  (and :154, :139, :208)
~Matrix() = default;            // was: virtual ~Matrix() = default;
...
// :276  (and cl_BZ_Matrix.hpp:309)
inline const T &
operator()( const size_t aRowIndex, const size_t aColIndex ) const   // was: virtual inline const T &
```

### Round B (`src/linalg/armadillo/armadillo.hpp`, before `#include <armadillo>`)

```cpp
// release: bounds and size-conformance checks off, like Blaze's BLAZE_USER_ASSERTION ( blaze_config.hpp );
// debug: on. ARMA_NO_DEBUG would not do this — it only drops ARMA_DEBUG. Same predicate as assert.hpp.
#if !defined( NDEBUG ) || defined( DEBUG )
// debug: Armadillo's default ARMA_CHECK_CONFORMANCE stays
#else
#  ifndef ARMA_DONT_CHECK_CONFORMANCE
#    define ARMA_DONT_CHECK_CONFORMANCE
#  endif
#  ifndef ARMA_WARN_LEVEL             // O2b (ii), decided 2026-11-07: silences the per-TU note;
#    define ARMA_WARN_LEVEL 1         // drops Armadillo level-2 runtime warnings in release only
#  endif
#endif
```

### Round C — Armadillo `Matrix` (`cl_AR_Matrix.hpp`)

```cpp
// constructor from any Armadillo expression ( eOp, eGlue, Glue, mtOp, subview, ... )
template< typename ET >
Matrix( const arma::Base< T, ET > & aExpression ) :
    mMatrix( aExpression.get_ref() ) {}

template< typename ET >
Matrix< T > &
operator+=( const arma::Base< T, ET > & aExpression )
{
    mMatrix += aExpression.get_ref();
    return *this;
}
// -= and *= likewise; the existing T, Matrix<T> and MatrixType overloads stay
```

### Round C — Blaze `Vector` (`cl_BZ_Vector.hpp`)

```cpp
// column expression: same orientation, no temporary
template< typename VT >
Vector< T > &
operator+=( const blaze::DenseVector< VT, false > & aExpression )
{
    mVector += *aExpression;
    return *this;
}

// row expression ( A.row(0), trans(x) ): transposed into the column
template< typename VT >
Vector< T > &
operator+=( const blaze::DenseVector< VT, true > & aExpression )
{
    mVector += blaze::trans( *aExpression );
    return *this;
}
// -= likewise; the DMatScalarMultExpr overloads stay or go per the R3.4 probe
```

`*aExpression` is Blaze 3.9's own unwrap (`DynamicMatrix.h:1666` uses `*rhs`); the member is `mVector`
(`cl_BZ_Vector.hpp:53`). The sketch stands on its own: the in-file row **ctor** (`:127-136`) is a
resize-and-copy loop and is not the pattern for a `+=`. The Armadillo twin for gap #32:

```cpp
// row view: transposed into the column. Non-template, so it beats the unconstrained template at :399.
// strans, not trans: arma::trans conjugates complex elements, Blaze's trans does not ( gap #37 )
Vector< T > &
operator+=( const arma::subview_row< T > & aRow )
{
    mVector += arma::strans( aRow );
    return *this;
}
```
The 1×n *expression* case cannot be routed by type (probed: `template<ET> f(const ET&)` beats
`template<ET> f(const Base<T,ET>&)` for every Armadillo type, `Mat` included) and follows O12. The existing
row ctor at `:129-131` changes `arma::trans` → `arma::strans` in the same round.

### Round C — the counter (`tests/linalg/test_LinalgAllocations.cpp`, sketch — matches gap #18)

```cpp
// Both backends allocate buffers through posix_memalign ( armadillo_bits/memory.hpp:72, blaze/util/Memory.h:94 ),
// never operator new. The hook is a strong, extern "C" definition in the test executable that forwards to
// the real allocator and counts events on the calling thread only. Frees are not counted: a temporary that
// is allocated and freed inside the expression must still show up.
static thread_local std::size_t gAllocEvents = 0;      // thread_local: an MPI progress thread cannot move it

extern "C" int posix_memalign( void ** aPtr, size_t aAlignment, size_t aSize )
{
    using Fn = int ( * )( void **, size_t, size_t );
    static Fn tReal = reinterpret_cast< Fn >( dlsym( RTLD_NEXT, "posix_memalign" ) );   // resolved once, on
    ++gAllocEvents;                                                                     // the first call,
    return tReal( aPtr, aAlignment, aSize );                                            // before MPI init
}

// positive control, first test in the file: a deliberate backend temporary must count >= 1,
// otherwise the instrumentation is broken and every other test in the file is meaningless
// per shape ( R3.1 list ), sizes chosen so a temporary exceeds Armadillo's 16-element local storage:
//   matrices 8x8 and 64x64, vectors 32 and 256; n = 3 is recorded, never used as an oracle
//   tBefore = gAllocEvents; <wrapper expression>;                 tWrapper = gAllocEvents - tBefore;
//   tBefore = gAllocEvents; <backend expression via matrix_data()>; tBackend = gAllocEvents - tBefore;
//   EXPECT_EQ( tWrapper, tBackend );   // parity, not zero: the backend itself allocates on alias
```

Platform notes, stated in the test: the symbol interposes for every call compiled into the executable,
which is where the test's expressions are instantiated, under static and shared `libbelfem` alike; on
macOS the executable's own definition wins for its own calls (two-level namespace). If the positive
control fails on a platform, that platform's gate is **incomplete** and the devlog says so — there is no
heap-statistics fallback (a size delta cannot see a temporary freed inside the expression). Blaze's OpenMP
assignment thresholds (gap #35) are far above these sizes, so the thread-local count is complete.

## 7. Definition-of-Done Checklist

- [ ] Every gap-table row maps to a step or an On.
- [ ] Each gap is backed by a citation re-opened on 2026-11-07, not by the report alone.
- [ ] Rounds ordered with dependencies; each has a named gate Christian runs.
- [x] Open questions logged, not silently decided — O1, O2, O2a, O2b, O4, O5 resolved by Christian on
  2026-11-07 and recorded in place; O3, O6, O7, O8, O9, O10 open.
- [ ] The plan-jury P1s (counter design #18, orientation contract #32, numerical gate #33, discarded
  status #36, complex transpose #37) are closed by tests that exist before the header they concern
  changes; the BDF baseline (R3.0) exists before any overload is added.
- [ ] End-to-end: `make check` green on Armadillo + Blaze, debug + release, after every round;
  the allocation counter green on both backends; the BDF2–5 tests exist before `cl_IWG_Timestep.cpp`
  changes.

## 8. Audit Trail

- Investigation thread: `tmp/ai_exchange/review_linalg_wrapper_design.md` (rounds 1–2 Codex
  `gpt-6-astra` high/xhigh + Copilot reserve; round 3 Grok `grok-4.7` xhigh on the report v5,
  2026-11-07 15:03 CST: ten findings, all CONFIRMED, one premise refuted for this machine; verified in
  the entry that follows it).
- Plan thread: `tmp/ai_exchange/review_linalg_wrapper_hardening.md` — pre-registration (11 findings,
  4 predictions), then **round 1, 2026-11-07 15:15–15:31 CST**: Codex `gpt-6-astra`/high (8 findings:
  3 P1 — counter blind to alloc+free, Armadillo row-operand UB after round B, no numerical gate — 5 P2),
  Grok `grok-4.7`/high (13 findings + open risks: `DEBUG` is defined, exception types, Blaze row ctor is
  a loop, `DMatScalarMultExpr` not subsumable, `=` not a hole, caller set, Calculator 3+6, …; one
  REFUTED by probe: `-fsyntax-only` does instantiate). Every citation re-checked; verification entry
  and reconciliation table in the thread; plan → v2. My own pre-registration #8 was wrong. **Round 2, 17:05–17:20 CST, `xhigh`:** Codex 5 findings (3 P1:
  `fn_AR_eigen.hpp` discards status; `arma::trans` conjugates vs Blaze; §6 counter sketch stale — plus
  `Base` overload ranking shown by run probe, BDF baseline before round C), Grok 4 P1 + 8 P2 (one-TU
  tripwire, API wording, row-view add is a transpose-add not UB, n = 3, `*=` citation, `:1129`, CN
  statement lines, release-Blaze throw test, Clang command, P3 slips) + open risks (implicit ctor → O13,
  shared-lib interposition, column-0 semantics). All CONFIRMED; none refuted. Plan → v3. Code rounds
  append to the same slug.
- Every auditor citation is re-checked against source before it enters this plan.

## Appendix A — Decision: `final`, not `static_assert`, seals the classes

`static_assert(!std::is_polymorphic<Matrix<T>>)` inspects the base only: a downstream
`class Foo : public Matrix<real> { virtual … }` is polymorphic without making `Matrix` so, and the
assert never sees it (Codex round 1, Copilot round 1 on the report). `final` makes the derivation
itself ill-formed at the first rebuild. The assert still earns its place as a tripwire against
someone re-adding a `virtual` to the base. Cost of `final`: none at runtime; it also lets the
compiler devirtualize any remaining virtual call statically, though after round A there is none.

## Appendix B — Decision: Blaze size-conformance stays on in release

Blaze's shape checks are exceptions (`BLAZE_THROW_INVALID_ARGUMENT`), not assertions
(`DynamicMatrix.h:1667`, `DMatDMatAddExpr.h:1092-1093`). The single customization point is
`BLAZE_THROW` (`Exception.h:93-95`). Redefining it to anything that returns (`(void)0`, a log) makes
the throwing branch fall through into the operation with mismatched sizes — out-of-bounds writes,
not a clean death. Redefining it to `abort()` keeps the check and only changes how it dies. So the
Christian-decided policy is reachable on Armadillo (one macro) and *not* reachable on Blaze without
patching the library; the cost of the Blaze check is one integer comparison per expression, far
below anything measurable in the report's probes. Accept, document, move on.

## Appendix C — Decision: constrained overloads, not a generic `template<ET>`

The commented-out `template<typename ET> operator+=(const ET&)` in `cl_BZ_Vector.hpp:458-466` looks
like the obvious fix and is wrong: with a row expression (`A.row(0)`) it forwards a `rowVector` into
a `columnVector` compound assignment and fails to compile (probed, report round 2). Constraining on
`blaze::DenseVector<VT, TF>` and dispatching on `TF` gives both orientations a temporary-free path
and keeps scalars and wrappers on their existing overloads. On Armadillo the analogous constraint is
`arma::Base<T, ET>`, which every expression type derives from and which `Mat<T>` itself satisfies —
the ambiguity with the copy constructor and the `MatrixType` overload is the thing R3.3's probe must
show absent (exact-match non-template overloads win over a template by the usual ranking; the probe
confirms, the plan does not assume).

## Appendix E — R2.0: backend status and conformance inventory (2026-11-07, read-only)

| # | site | today (checks on) | after round B (Armadillo release) | Blaze (every mode) | action in round B |
|---|---|---|---|---|---|
| E1 | `fn_AR_eigen.hpp:33` `arma::eig_gen( tValues, … )`, `bool` discarded; `:42` indexes `tValues(k)` | nonfinite input / `geev` failure → `soft_reset()` → `std::out_of_range` at `:42` | out-of-bounds read | `blaze::eigen` throws (`geev.h:140,:198`) | check the status; `BELFEM_ERROR( status, … )` (setup-grade call; `eigen_sym` already does this at `:77`) |
| E2 | `fn_AR_cross.hpp:44-46`, `fn_BZ_cross.hpp:44-46`: `cross( const ET & aA, const Vector<T> & aB )` body uses `aA.vector_data()` | ill-formed if instantiated; no caller does (39 `cross` sites, the mixed form unused) | same | same | body → `cross( aA, aB.vector_data() )` |
| E3a | `dot` of two `Vector`s of different length (`fn_dot.hpp`, no assert) | `std::logic_error` (`op_dot_meat.hpp:284-689`) | unchecked read over the shorter operand | `std::invalid_argument` (`DVecDVecInnerExpr.h:271`) | `BELFEM_ASSERT( aA.length() == aB.length() )` in `fn_dot.hpp` |
| E3b | `cross` with a non-3 operand | `std::logic_error` (`glue_cross_meat.hpp:37`) | unchecked | `std::invalid_argument` (`DVecDVecCrossExpr.h:806`) | already asserted (`fn_cross.hpp`, 1 assert) — confirm it covers every overload |
| E3c | `det` of a non-square `Matrix` (`fn_det.hpp`, no assert) | `std::logic_error` (`op_det_meat.hpp:54,:101,:127`) | unchecked | Blaze throws | `BELFEM_ASSERT( n_rows == n_cols )` in `fn_det.hpp` |
| E3d | `inv` of a non-square `Matrix` (`fn_inv.hpp`, no assert) | `std::logic_error` | unchecked | Blaze throws | `BELFEM_ASSERT` in `fn_inv.hpp`; singular stays a throw on both (`op_inv_gen_meat.hpp:36,:72`, `Inversion.h:153`) |
| E3e | element access `A(i,j)` | Armadillo `std::out_of_range` + `BELFEM_ASSERT` | `BELFEM_ASSERT` only (compiled out) | `BELFEM_ASSERT` only | nothing — this is the decided policy |
| E5 | `dot( Matrix, Vector )` and `dot( Vector, Matrix )`, documented in `fn_dot.hpp:38-57` as flat-sequence inner products | Armadillo: `arma::dot` accepts them | unchanged | **do not exist**: `blaze::dot` has no matrix×vector overload, the `decltype` SFINAEs both wrappers out (`fn_BZ_dot.hpp:37-49`), a call is a compile error (probed 2026-11-07); no in-tree caller, or the Apple build would not compile | round C item: either give Blaze the flat semantics explicitly (loop over `data()`… no — through `(i,j)`) or retire the two overloads and their doc on both backends. **O14** |
| E6 | `inv` of a non-square `Matrix` on **Blaze** | `std::invalid_argument` "Invalid non-square matrix provided" from the inverse expression (`DMatInvExpr.h:410`), every mode (**probed** 2026-11-07; the `BLAZE_INTERNAL_ASSERT`s in `Inversion.h` are not the only check — a first reading of the headers said otherwise and was wrong) | — | throws | nothing beyond the new `BELFEM_ASSERT`; guide table row kept as one line |
| E7 | **`eigen()` on a non-finite matrix under Blaze returns zero eigenvalues silently** (**probed**: `blaze::eigen` on `{{NaN,0},{0,1}}` → no throw, `w = (0,0)`, while LAPACK prints `** On entry to DGEEV parameter number 4 had an illegal value` — an `info < 0` the Blaze wrapper does not throw on: `geev.h:139-142` has an internal assertion for `info >= 0` and throws for `info > 0` only) | same today | — | silent garbage in every mode | **pending, post-jury fix (D-tracker):** a finiteness pre-check with `BELFEM_ERROR` in both `eigen()` implementations (setup-grade call, n² scan acceptable), so the two backends fail the same way; the Armadillo test of `test_CheckPolicy.cpp` then becomes backend-neutral |
| E4 | LAPACK wrappers `fn_gees/geev/gels/gesv/gesvd/getrf/getri/posv.hpp` | every `info` checked with `BELFEM_ERROR` | unchanged (not Armadillo checks) | — | nothing |

The symmetry rule that falls out of E3: after round B the two backends agree in **debug** (BELFEM assert or
Blaze exception, both fatal) and differ in release (Blaze throws, Armadillo computes on) — that is O2a, accepted,
and the guide table of R2.6 states it per operation.

## Appendix D — Hatch inventory (filled in R4.0)

| file | line | expression | shape | live? | kind | replacement | covering test |
|---|---|---|---|---|---|---|---|
| (to be filled) | | | | | | | |
