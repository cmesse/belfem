# linalg wrapper design: MORIS-style aliasing weighed, two real wrapper costs found instead

**Date:** 2026-11-07
**Purpose:** Christian asked whether `src/linalg` should follow what he understood to be the newer MORIS
design — a `Matrix` that aliases the backend object directly, without the wrapped `mMatrix` member and
`matrix_data()` — and, if so, what the impact radius and plan would be. Read-only investigation with a
three-round jury; the deliverable is `todo/linalg_wrapper_design_investigation.md` (v6). The same day
carried on into the implementation plan (`todo/linalg_wrapper_hardening_plan.md`, v3, two plan juries),
round A of it (de-virtualize + `final`, code-juried, gated on this Mac), the O3 follow-up
(`todo/linalg_small_matrix_kernels.md`, v2, plan-juried) and the gate's fallout (DR-202, DR-203).
**Module:** `src/linalg` (`cl_AR_Matrix.hpp`, `cl_BZ_Matrix.hpp`, `cl_AR_Vector.hpp`, `cl_BZ_Vector.hpp`),
with the consumer census over the whole tree
**AIs involved:** Claude Fable 5.1 (investigation, probes, report); Codex `gpt-6-astra`/high (round 1) and
`gpt-6-astra`/xhigh (round 2); **Grok `grok-4.7`: seat failed both rounds** — the `grok` CLI is not signed
in on this machine (`grok login --device-code`) and no `XAI_API_KEY` is set; reserve seat Copilot
(`auto` → `gpt-5.6-luna`, effort ignored) sat in both rounds per protocol §9.2; **Grok `grok-4.7`/xhigh
round 3** after Christian restored the CLI (15:03 CST, blind, on v5). Plan juries on the hardening plan:
Codex `gpt-6-astra`/high + Grok `grok-4.7`/high (round 1), both at `xhigh` (round 2). Code jury round A:
Codex `gpt-6-astra`/high + Grok `grok-4.7`/high. Plan jury on the O3 plan: Codex `gpt-6-astra`/high +
Grok `grok-4.7`/high
**Status:** investigation COMPLETE (report v6); hardening plan v3 APPROVED by Christian (plan + source edits;
O1, O2, O2a, O2b, O4, O5, O12 decided); **round A implemented, code-juried, D1/D2 applied, gated on this Mac
(Blaze, release)**; O3 plan v2 juried, awaiting Christian's R1 profile. **The round-A source diff (7 files) is
left uncommitted until the Linux CI confirms the Armadillo and debug legs**; this commit is the paper trail.
Rounds 1–2 on the report had two OpenAI-hosted voices; every later round had Codex + Grok
**Literature References:** none required (protocol §5); Oliveira & Stewart 2006 §7.6 and Rouson et al.
2011 consulted for the wrapper-interface rationale
**Verification:** probes (scratchpad compiles against the tree's headers and `/opt/scls/include`, `-O2
-DNDEBUG`, GCC 15): object sizes, assembly, timings, allocation call sites, inheritance compile tests,
overload-ranking and instantiation probes, the O3 Jacobian-chain benchmarks. Round A: syntax gate in five
configurations (g++ Blaze/Armadillo, Apple clang Blaze, both backends `-DDEBUG`), then **Christian's `make check`
in `cmake-build-debug` (Blaze, release, clean rebuild): `linalg` passed with the new tripwires; 54 failures in
io/circuit/physics/mesh/the-nonfree-python-extension all traced to pre-existing causes (DR-201 the-nonfree-python-extension extension, DR-202 VTK's
mangled HDF5 in the test include order, DR-203 never-run TensorMeshBSpline tests), none to the diff.** The
Armadillo and debug legs are owed to the Linux CI matrix — round A is *verified on one leg*, reviewed on the
others
**Exchange:** `tmp/ai_exchange/review_linalg_wrapper_design.md` (report, 3 rounds),
`review_linalg_wrapper_hardening.md` (plan rounds 1–2, code round A, pre-registrations, probes),
`review_linalg_small_matrices.md` (O3, probes, round 1)
**Exchange:** `tmp/ai_exchange/review_linalg_wrapper_design.md` (pre-registration, two jury rounds,
two verification entries, two reconciliation tables, probe sources)

## Summary

The premise did not survive contact with the MORIS clone: every MORIS branch (32 remotes, HEAD 2026-10-12)
and GitHub `main` still hold `arma::Mat<Type> mMatrix` behind `matrix_data()`; nothing in MORIS aliases
the backend. The question was therefore answered on its merits. A complete, non-virtual member wrapper is
runtime-identical to the backend; aliasing is incompatible with the two-backend contract and L-21;
transparent inheritance compiles under Blaze 3.9 and fails under Armadillo 15.2 (exact-type traits). What
the investigation did find is that BELFEM's wrapper is neither non-virtual nor complete, and pays twice
for it.

## Key Findings

- **MORIS still wraps a member.** `tmp/moris/projects/LINALG/src/cl_Matrix.hpp:35`,
  `Arma_Impl/cl_Matrix_Arma_Dynamic.hpp:51`, all 32 `origin/*` branches, GitHub `main` fetched
  2026-11-07. Newest LINALG commit on `main`: 2024-08-29. Confidence: high.
- **Cost 1, the `virtual`s.** `cl_AR_Matrix.hpp:106,276`, `cl_BZ_Matrix.hpp:139,309`,
  `cl_AR_Vector.hpp:154`, `cl_BZ_Vector.hpp:208`. Nothing derives (0 hits in `src tests examples nonfree`).
  Probe: `sizeof` +8 (Blaze) / +16 (Armadillo); GCC emits a speculative-devirtualization guard with an
  indirect-call fallback on every const element read through `const Matrix&`; a 3×3 mat-vec through
  const refs runs 2.64× (Blaze) / 1.28× (Armadillo, equal accessors) slower than the backend, and at
  parity (0.99× / 1.04×) with the six `virtual`s removed. Confidence: high on mechanism, medium on the
  ratio (one session, one CPU); the real-kernel magnitude is the R0 gate.
- **Cost 2, overload gaps.** Found by Codex (round 1) and widened in verification and round 2: Armadillo
  `Matrix` `+=`/`-=` with an expression materializes an `arma::Mat` temporary (heap above 16 elements);
  Blaze `Vector` `+=`/`-=` with *any* expression (`x*dt`, `x+z`, `A.col(0)`, `A.row(0)`) materializes a
  temporary `Vector` (the generic overload at `cl_BZ_Vector.hpp:458-466` is commented out);
  `Matrix C( B*dt )` allocates twice on both backends where the backend allocates once. Static call-site
  counts, not executed counts — the R2a gate is an executed allocation counter. **The `cl_IWG_Timestep.cpp`
  escape hatches (79 of the 109 editable hatch lines) are workarounds for exactly the `+=` gap**, not
  friction. Confidence: high on the sites, unmeasured on execution at BELFEM's sizes.
- **Recipe withdrawn in round 2:** restoring the commented-out generic Blaze `+=` verbatim makes
  `y += A.row(0)` fail to compile (row expression into a column vector); constrained overloads are needed.
- **Armadillo macro correction:** `ARMA_NO_DEBUG` does not disable bounds checks; `ARMA_DONT_CHECK_CONFORMANCE`
  does, and takes the size-conformance checks with it (`config.hpp:333-357`, `debug.hpp:1383-1398`).
  `linalg_usage_guide.md:98` names the wrong macro — a one-line doc fix is owed. With the checks off the
  3×3 kernel measured no difference (1.00×).
- **Impact radius** (recounted in round 2): 552 consumer files outside `src/linalg` (498 `src`, 54
  `nonfree`); 110 linalg escape-hatch lines in 13 files, 1 in the closed cohomology core; `.spacing()`
  only in `commtools.hpp` (4 lines); the LAPACK wrappers consume the BELFEM API (55 lines).
- **Plan in the report:** R0 measure a real kernel (Christian runs it); R1 remove the six `virtual`s and
  add `final` (10 declaration lines, zero in-tree call-site changes, a deliberate downstream API break —
  O5 — bounded by the backend-free plugin contract); R2a audit and complete the overload sets with
  constrained overloads, gated by executed allocation counters; R2b retire the hatches on live paths
  only, after a per-site inventory and with named numerical tests (BDF2–5 assembly tests do not exist,
  `test_BdfTimestepMethod.cpp:25-27`; Crank–Nicolson/Galerkin are deliberately disabled and stay). R3
  inheritance struck.

## Jury record

| round | seat | depth | findings | outcome |
|---|---|---|---|---|
| 1 | Codex | `gpt-6-astra`/high | 6 | all CONFIRMED; refuted "the `virtual`s are the only cost" and the Armadillo benchmark column (unequal accessors `()` vs `[]`) |
| 1 | Grok | `grok-4.7`/high | — | FAILED: CLI not signed in |
| 1 | Copilot (reserve) | `auto`/n.a. | 4 | all CONFIRMED (tripwire is `final`, not the `static_assert`; R1 is an API break; 10 lines; probes into the record) |
| 2 | Codex | `gpt-6-astra`/xhigh | 9 | all CONFIRMED (R2a recipe regression; call sites ≠ executions; construction gap; six count/scope corrections; CN/Galerkin disabled; plugin API backend-free; O2 scope) |
| 2 | Copilot (reserve) | `auto`/n.a. | 6 | all CONFIRMED (two duplicate Codex; O5; guide macro; `-=` gate; per-site inventory) |
| 3 | Grok | `grok-4.7`/xhigh | 10 | all CONFIRMED against source; one premise (Apple builds are Clang) REFUTED for this machine (SCLS `g++`). P1: Blaze `Matrix` expression `+=`/`-=`/`*=` exist only under `BELFEM_GCC` (`cl_BZ_Matrix.hpp:593-683`); Armadillo `Vector` ctor hole; `Vector::operator() const` not virtual; backend `+=` allocates on alias; two §9 rows lowered. Report → v6 |

Nothing refuted in either round; nothing lifted past the probe rung. Stop condition met after round 2.

## Tooling defects found on the way (not fixed — `.claude/scripts/` edits need their own round)

- `.claude/scripts/ask_copilot.sh:167` expands `"${EFFORT_ARGS[@]}"` under `set -u`; with
  `COPILOT_MODEL=auto` the array is empty and the shebang's `/bin/bash` 3.2.57 aborts with
  `EFFORT_ARGS[@]: unbound variable`. Worked around by running it under `/opt/bash/bin/bash` 5.3. Fix
  candidate: `${EFFORT_ARGS[@]+"${EFFORT_ARGS[@]}"}` or `#!/usr/bin/env bash`.
- `scripts/cross_review.sh --jury` reports a Grok auth failure only as "Audit FAILED" in the exchange;
  the cause (`Not signed in`) is in `tmp/ai_exchange/_cross_review.log`.
- The Grok seat is unavailable on this machine until `grok login --device-code` (or `XAI_API_KEY`) is
  done; the jury this session ran with two OpenAI-hosted voices.

## Changes Made / Proposed

- New: `todo/linalg_wrapper_design_investigation.md` (v5) and its line in `todo/README.md` (Active 9 → 10).
- Proposed, pending Christian: R1, R2a (jury rounds), R2b (conditional), the `linalg_usage_guide.md:98`
  macro fix, O1–O5.
- No source file touched.

## Same day, second half: implementation plan and its first jury

Christian's ruling (afternoon): both backends stay; size-conformance checks on in debug, off in release;
implementation in rounds, virtuals first, consumers last; one jury round on the plan before any source
edit; scripts run under bash 5.3 (`/opt/bash/bin/bash`). Plan written: `todo/linalg_wrapper_hardening_plan.md`
(rounds 0/A/B/C/D/E, gap table #1–#34, O1–O10, exact edits in §6). Grok round 3 on the report (above) was
folded into the plan before the plan jury.

Plan jury round 1 (`tmp/ai_exchange/review_linalg_wrapper_hardening.md`, Codex `gpt-6-astra`/high +
Grok `grok-4.7`/high, blind): 8 + 13 findings, all CONFIRMED against source except one REFUTED by probe
(`-fsyntax-only` does instantiate used member templates — a planted error is diagnosed). Three P1s changed
the plan's gates: the allocation counter must be a monotonic `thread_local` `posix_memalign` event count
with a positive control (a heap-size delta cannot see a temporary freed inside the expression; neither
backend uses `operator new`); the Armadillo `Vector` forwards a row operand unchanged into an n×1 `Mat`,
which is `std::logic_error` today and undefined behaviour once release checks are off, so round C owes one
orientation contract on both backends; and the new overloads need value/shape tests of their own. Also
refuted in my own text: `DEBUG` **is** defined under `USE_DEBUG` (`config_gcc.cmake:24`), so the round-B
predicate is `!defined(NDEBUG) || defined(DEBUG)` as in `assert.hpp`, and tests guard on
`BELFEM_ASSERTIONS_ACTIVE`. Two probes recorded in the thread (copy-init rejection; template instantiation
under `-fsyntax-only`). Plan → v2. At the STOP Christian decided: plan round 2 at `xhigh` on v2 before approving source
edits; O2b → `ARMA_WARN_LEVEL 1` in release; O2a → the Blaze asymmetry is accepted; O4/O5 → `final` on
both classes, the downstream break accepted. Round 2 (`xhigh`, blind, 17:05–17:20 CST) found three more
P1s, all confirmed: `fn_AR_eigen.hpp` discards `arma::eig_gen`'s status and would read out of bounds once
release checks are off (new prerequisite R2.0: audit every `fn_AR_*`/`fn_BZ_*` for discarded status);
`arma::trans` conjugates complex elements while Blaze's `trans` does not, and the existing Armadillo row
constructor already diverges from Blaze for `Vector<cplx>` (fix with `strans`, complex tests); the §6 counter
sketch still described the rejected heap-delta oracle (rewritten: `extern "C"` forwarding `posix_memalign`
interposer, thread-local monotonic count, positive control). Grok corrected the plan's own text twice:
`y += A.row(0)` on Armadillo with checks off is a silent transpose-add, not undefined behaviour, and a
`Base<T,ET>`-constrained overload can never beat the unconstrained `const ET&` template (shown by a run
probe) — so the Armadillo row contract needs a non-template `subview_row` overload and a decision on row
*expressions* (O12). Plan → v3. Christian approved v3 and source edits and chose O12 (ii). **Round A implemented** the same evening:
`final` on the four class heads, `virtual` removed from four destructors and two `Matrix::operator() const`,
`!is_polymorphic` + `sizeof` tripwires in `tests/linalg`, guide item 6; syntax gate green in five configurations
(g++ Blaze/Armadillo, Apple clang Blaze, both backends with `-DDEBUG`). Code jury (Codex astra/high + Grok 4.7/high,
blind, own-files diff): no P0/P1; two P2s — pin `final` with `std::is_final` (Grok), soften "exactly its backend
object" (both) — logged as D1/D2, approved and applied the same evening (`std::is_final` asserts added;
guide and assert messages reworded; syntax gate re-run green). Debt row DR-200 opened. Nothing built or run: **reviewed,
not verified**; the gate (clean rebuild, `make check`, both backends) is Christian's.

## O3 — small matrices (same evening, Christian's request)

Three probes on the per-integration-point Jacobian chain (`J = dNdXi·X`, closed-form inverse,
`dNdX = invJ·dNdXi`; hex8, quad4, hex27; GCC 15 `-O2`, tree headers after round A). Fixed-size `J`/`invJ`
(`arma::Mat::fixed`, `blaze::StaticMatrix`) buy 12 % of the chain on Blaze and nothing on Armadillo, whose
chain is 3.1× slower than Blaze with the OpenBLAS thread count pinned (5.3× unpinned): Armadillo's `gemm`
takes its inline path only for square matrices ≤ 4×4 (`mul_gemm.hpp:413`), so every 3×8 product is a
`dgemm` call. A hand `(i,j)` loop beats that by 1.4× (hex27) to 6.9× (quad4) on Armadillo and loses to
Blaze's own kernel at hex27 — so the item is an Armadillo-side kernel behind the backend split, not a
second matrix type. Plan `todo/linalg_small_matrix_kernels.md` v2 after one plan jury (Codex astra/high,
Grok 4.7/high): Codex's P1 — the calculator's `dN` scratch is created 2×2n / 3×3n and silently resized by
the backend on first use, which a size-asserting kernel would catch — is R2.0; Grok corrected the dispatch
citation, an inverted ratio, the consumer census and the scheduling text, and one of its claims was refuted
(the MORIS typedef file exists). Next step is R1, Christian's profile of the chain's share on one 2-D and
one 3-D deck per backend, with both thread variables pinned.

## The round-A gate and what it turned up (evening)

Christian ran `make check` in `cmake-build-debug` (Blaze, release, every object newer than the header edit).
`linalg` passed. The 54 failing cases split three ways, none touching the diff. **DR-202:** every test
target compiles with `-isystem /opt/scls/include/vtk-9.6/vtkhdf5/src` before `-isystem /opt/scls/include`
(`test_io` flags, positions 28 vs 51), so `<hdf5.h>` in a test TU is VTK's bundled copy with
`vtk_hdf5_mangle.h` (`:3202` renames `H5Tget_class`), while `libbelfem.a` was compiled with the SCLS
directory first; six test objects reference `vtkhdf5_*`, two HDF5 libraries exchange type ids in one process,
and every HDF5 read fails with "not a datatype", class −1 (`hdf5_tools.hpp:246`) — 51 of the 54. Library
objects are clean; `USE_VTK` is ON in this cache since Oct 26 and OFF on CI. **DR-203:** the three
`TensorMeshBSpline` tests, added 2026-10-22 as "test not run", build their mesh in memory and fail with zero
node sources. **DR-201** (the-nonfree-python-extension, extension not importable) is another session's (belfem-5e). Both new
rows filed unowned; ids coordinated with that session.

## Open Questions

- O1 keep both backends; O2 `ARMA_DONT_CHECK_CONFORMANCE` (needs an inventory of what the checks catch,
  not a timing); O3 fixed-size small matrices; O4/O5 `final` and its API break — all Christian's.
- Whether the §3b call sites execute at BELFEM's element sizes (R2a's counter decides).
- A Grok pass on the report before the decision.

## Files Updated

- `todo/linalg_wrapper_hardening_plan.md` (new, v3, IN PROGRESS round A), `todo/linalg_small_matrix_kernels.md` (new, v2),
  `todo/debt_register.md` (DR-200 opened; DR-202, DR-203 filed), `todo/README.md` (three index lines; investigation line corrected),
  `tmp/ai_exchange/review_linalg_wrapper_hardening.md` (new thread), `tmp/ai_exchange/review_linalg_wrapper_design.md`
  (Grok round 3 + verification), `tmp/ai_exchange/review_linalg_small_matrices.md` (new), this devlog.
  **Uncommitted on purpose until CI:** `src/linalg/{armadillo,blaze}/cl_{AR,BZ}_{Matrix,Vector}.hpp`,
  `tests/linalg/test_Matrix.cpp`, `test_Vector.cpp`, `src/linalg/doc/linalg_usage_guide.md` (round A).


- todo/linalg_wrapper_design_investigation.md (new)
- todo/README.md
- devlog/dl20261107_linalg_wrapper_design_investigation.md (this file)
- devlog/README.md
