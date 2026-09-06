# Coding Philosophy Revision — Jury Feedback Harvest

**Date:** 2026-08-08
**Purpose:** Full source-verified revision of `doc/coding_philosophy.md` per the jury-feedback work order (Phase 1 verification V1–V11, Phase 2 edits E1–E12, Phase 3 proposals recorded only)
**Files touched:** `doc/coding_philosophy.md` (rewritten), `tmp/whitepaper/coding_philosophy_CHANGES.md` (Phase 1 ledger + Phase 3 proposals appended)

## Method

Four parallel read-only source sweeps (allocation/ownership, comm layer,
container/smart-pointer audit, build/CI), every claim resolved to `file:line`
before editing. Per-claim verdicts and citations in
`tmp/whitepaper/coding_philosophy_CHANGES.md` (2026-08-08 section). No source
code was modified; standing constraints (raw pointers stay; Phase 3 =
proposals only) honored.

## Headline verification results (confidence: high — all file:line-cited)

1. **`-fno-exceptions` was never real.** No build file sets it; release
   (`USE_DEBUG=OFF`) is `-O2` + `NDEBUG` (GCC/Clang), `-O1 -xHost` (ICC) —
   never `-O3`. `assert.hpp:99` contains a `throw` compiled into every default
   build, so a genuine `-fno-exceptions` build would not compile. Doc now
   presents exception-free production as design intent. **`CLAUDE.md:115`
   carries the same false flag claim — not edited (out of scope), needs a
   follow-up.** Decisions previously made against the abort model (e.g. the
   try/catch rejection in `todo/coreduce_…_findings.md:335`, ~24 throwing
   `std::stod` sites) should be revisited.
2. **Latent Blaze-backend MPI bug (code, not doc).** Matrix transfers ship the
   raw `capacity()` buffer by design (padding round-trips correctly), but
   Blaze `resize(m,n,false)` never shrinks `capacity_`: a shrunk persistent
   matrix reports a stale-large size and the receiver — fresh, exact-fit —
   posts receives for the sender's count with no size assertion
   (`commtools.hpp:1841-1874`). Reachable via
   `cl_FEM_DofMgr_SolverData.cpp:2026` (`send(mRhsMatrix)`) → `collect()`.
   Armadillo backend immune. Needs a fix decision + regression test.
3. **Mesh allocation reality:** per-object `new`/`delete`, hierarchical
   ownership (Mesh→Block→Element, SideSet→Facet); no pool/arena exists; the
   doc's placement-new pool example was aspirational and is gone.
   `ShiftRegister` is the one correct placement-new pool in the tree.
4. **Alignment reality:** zero `aligned_alloc`/`posix_memalign` calls in
   `src/`; alignment is backend-provided (Armadillo ≤32 B, Blaze 32 B under
   AVX2 / 64 B only under AVX-512). The doc's 64-byte story and its defective
   hand-rolled example are gone.
5. **Smart-pointer census:** 15 hits; dominant use is factory-held
   `shared_ptr<Kernel>`/`shared_ptr<Controller>` (Maxwell/thermal), zero in
   `tests/` and `src/io/` — the old "I/O and unit tests" characterization was
   wrong on both counts.
6. **Rule of Five:** `ShiftRegister`, `DynamicBitset`, `StringList` compliant.
   Open gap: `Logger` owns a `FILE*` and is copyable (double-`fclose`
   candidate, `cl_Logger.hpp:48`).
7. **No CI, no sanitizers, Valgrind manual-only; `USE_TEST` defaults OFF.**
   Build bug found: `tests/core` is registered with ctest but missing from the
   `check`/`check-fast` dependency lists (`CMakeLists.txt:298-299,315`), so
   `make check` never builds `test_core`.
8. **Tag ceiling confirmed** (~32,768 ranks at INT_MAX `MPI_TAG_UB`,
   `cl_Communicator.cpp:169-187`) and now documented as a limit; ">100,000
   cores" removed. Redesign is comm-layer-local (no public API takes a tag;
   tags recomputed symmetrically) — small blast radius.

## Document changes (E1–E12, all applied)

Memory chapter rebuilt around "allocate deliberately, not manually" with the
entity rule and ownership rule; false smart-pointer claims removed
(`unique_ptr` double-allocation, non-deterministic-RAII); `Cell` rationale
rewritten (wraps `std::vector`, interface-and-discipline argument, real API in
examples — `size()`, free-function `unique()`); ownership annotations +
sanctioned owning-`Cell<T*>`-holder pattern (SideSet verbatim example);
canonical per-object lifecycle example; error-handling chapter now
semantic-rule-primary with the frequency rule for placement and a third
"expected algorithmic failure → status/retry" category; OpenMP chapter
reframed as a scope decision with the STRUMPACK `MPI_THREAD_MULTIPLE`
exception called out; MPI-4 large-count acknowledgment; all unmeasured
performance figures deleted (surviving heuristics cited to Drepper/Agner Fog
and labeled); "Benchmarks We Owe Ourselves" section lists the measurements
that would upgrade each claim; naming rules deduplicated and the Files rule
corrected (`cl_`+CamelCase etc.); ⚠️ audit annotations moved to a Source
Reconciliation appendix; Further Reading updated (Ousterhout, Lakos, Drepper,
Agner Fog, Core Guidelines Per.*/ES.84, Acton; Clean Code borrowings scoped).

## Phase 3 proposals (NOT implemented — await sign-off)

P1 snippet-extraction CI, P2 `owner<T*>`/`observer<T*>` + clang-tidy,
P3 `length()` alias, P4 debug-asserted span view, P5 Rule-of-Five fixes
(Logger first), P6 ASan/UBSan CI (+ `tests/core` dep fix), P7 tag-scheme
redesign. Estimates in the CHANGES ledger. The V3 matrix-receive overflow fix
is recommended to jump the queue (correctness, ~2 h + test).
