# Devlog 2026-08-31 — Documentation Currentness Sweep (Claude/Codex/Grok)

**Date:** 2026-08-31
**Purpose:** Record of the overnight three-AI sweep of all 96 markdown files in `doc/` and
`src/*/doc/`, covering currentness, readability, theory alignment against `./literature`,
magnet-design scope coverage, and the framework's community value proposition
**Module:** cross-cutting

---

## Summary

Christian asked five questions of the documentation tree: are the claims still true, is it
readable, does the theory match `./literature`, what magnet-design scope does BELFEM cover
(against `literature/books/magnets/`), and what benefit the code brings to the community.

The sweep ran as a jury: six disjoint batches covering all 96 files, each audited blind by
Codex (`gpt-5.6-terra`, `high`) and Grok (`grok-4.6`, `high`), pre-registered in
`tmp/ai_exchange/doc_currentness_sweep.md` before any auditor data existed. Every load-bearing
citation was re-derived by Claude against the tree before acceptance.

**All twelve auditor runs completed.** **The single most consequential finding: `hphirun` and `hphiTrun` are no longer built, and at
least six documents still instruct users to run them — including the root `README.md` and
`doc/getting_started.md`, the two files a newcomer reads first.**

## Key Findings

### The executables are gone from the build, not from the docs

`src/executables/CMakeLists.txt:40-47` comments out both `Add_Executable` blocks under the note
"the old executables are obsolete and will eventually be retiredI" (sic). Only `belfem` and
`electricalCircuit` are built. `cmake-build-debug/bin/`, built 2026-08-31 02:56, contains
`belfem`, `electricalCircuit`, `material`, `gas`, `db2exo`, and no `hphirun` or `hphiTrun`.

The `.cpp` sources still exist, so documentation citing `hphirun.cpp` **as source** stays valid.
Only *run* instructions break — a distinction that keeps the severity honest.

This finding also produced the round's clearest methodological lesson. Codex reported it as P0.
Grok contradicted it, calling the path "stale, not immediately fatal", on the strength of
`examples/README.md:38-40` — "`hphirun` and `hphiTrun` still exist and still work". Grok read a
*document* where the tree was available. That sentence in `examples/README.md` is itself false
and joins the fix list. Iron rule 1 paid for itself in a single round.

### The comm guide describes three non-collective primitives as collective

Found by Codex, widened by Grok's B5 round, verified here by counting the MPI calls in each
implementation. `comm_usage_guide.md:881-890` calls `share(data)` an "all-to-all broadcast" and
shows it invoked unconditionally on every rank. `share(Vector)` has **0 receives**;
`distribute(Cell<T>)` (`commtools.hpp:766-812`) has **0 receives** yet `:698-709` calls it on all
ranks and claims "data on rank i contains i*10"; `collect(Cell<T>)` (`:821-862`) has **0 sends**
yet `:807-813` presents it as a standalone gather.

CLAUDE.md states the contract directly — the rank guard is required because calling either on the
wrong rank is a deadlock — and `doc/coding_philosophy.md:537-541` teaches the correct guarded
form. The module's own usage guide teaches the hang the framework documentation exists to prevent,
and it is a systematic error rather than three typos: the guide reasons about BELFEM's asymmetric
API as if it were MPI collectives.

The same file's "Broadcasting Uninitialized Data" pitfall at `:103-120` is self-refuting — the
WRONG and CORRECT snippets are the same `broadcast(v, 0)` call with the same setup, differing only
in the comment.

### Two taught code examples do not compile, and one silently computes the wrong answer

- `src/fem/maxwell/doc/maxwell_usage_guide.md:945,951` teaches
  `material->rho( norm(b), angle, T )`. The signature is `rho( T, B, beta )`
  (`cl_Material.hpp:680`). All three parameters are `real`, so the wrong order **compiles** and
  returns a wrong resistivity — the all-`real` trap named in `lessons_learned.md` L-06. The
  sibling `materials_usage_guide.md:70,259` teaches the correct order, so the convention is not
  in doubt. Found by Claude.
- `src/fem/kernel/doc/dof_manager_usage_guide.md:2332-2337` teaches an object API
  (`comm.send(data, target, tag)`, `comm.sum_all(...)`) that does not exist. BELFEM's MPI layer
  is free functions with no tag parameter, and `sum_all` appears nowhere in `src/`. Found by
  Claude.
- `iwg_usage_guide.md:263-298` calls `create_iwg` with three arguments; the factory declares two
  (`cl_IwgFactory.hpp:53`). `maxwell_usage_guide.md:506-568` calls `new Kernel(mesh)`; the
  constructor takes `KernelParameters *` (`cl_FEM_Kernel.hpp:98`). Both found by Codex, both
  re-derived here.

### Defaults documented backwards

`USE_DEBUG` is documented as defaulting **ON** in `doc/getting_started.md:24-27`,
`doc/coding_philosophy.md:56` and `README.md:25`. It defaults **OFF**
(`find_scls_flavor.cmake:23`; `CMakeLists.txt:98`) unless `$SCLS` names the `debug` flavor. A
newcomer following the quick start gets a release build with assertions compiled out while
believing the opposite. Found independently by Codex, Grok and Claude — the only triple-agreement
of the round, and the one place where agreement was cheap because the check is mechanical.

The Coulomb-gauge penalty is documented as "off by default" in two places
(`coulomb_gauge_penalty_theory.md:11`, `src/fem/maxwell/doc/README.md:33`). It has been **on** at
`chi = 1e-4` since 2026-08-27 — `cl_IWG_Maxwell.cpp:48-51` says so in a comment written for
exactly this purpose.

### Countable claims have drifted

Layer 1 of `lessons_learned.md` is 200 lines (19→219), documented as "141 lines"
(`doc/README.md:12`) and "~140 lines" (`CLAUDE.md:53`). Layer 2 has 21 cards (L-01..L-21),
documented as 18. The incident catalogue holds 558 distinct IDs on 593 rows; the docs say 539 and
537. `check_doc_claims.py` guards none of these — it passes 37/37 and is structurally blind to
them.

Separately, `CLAUDE.md:520-523` lists 23 modules with a `doc/` directory; there are 26
(`fem/postproc`, `fem/thermal`, `visualizer` are missing). The adjacent "only one of the 23
carries a DOI" remains **true** in substance — that one is `src/homology/doc/README.md`.

### Theory: one open formulation question, one misattribution

Claude predicted a document would reproduce the pre-erratum Arsenault air-domain coupling, then
grepped for `E = 0 in air` phrasings, found none, and reported it clean. Codex found what that
grep structurally could not: `maxwell_usage_guide.md:329-369` states the air region as the static
`∇·(μ∇φ) = 0` Laplace form, while `mt_maxwell_phi.cpp:26-46` assembles a φ mass contribution. Doc
and code differ on whether a time derivative is present in the φ region. **This is a formulation
question and is routed to Christian, not decided here** (escalation trigger 2).

`maxwell_usage_guide.md:186-210` and `CLAUDE.md:699,711` attribute static condensation to
Messe et al. 2023. Read firsthand at `literature/papers/fem/messe2023.txt:509-517`, the paper
says the opposite of the attribution: *"For the sake of simplicity, we opted to couple the φ-field
and the in-plane field h_t using Lagrange multipliers… we wish to highlight that static
condensation may be the preferred method."* Static condensation is credited there to Alves et al.
in GetDP. The engineering advice in the guide is sound and matches what BELFEM does; only the
citation is wrong.

No BELFEM document cites the 2026 erratum beside its Arsenault 2023 citations
(`dof_manager_usage_guide.md:2825,3074`, `iwg_usage_guide.md:2061`), though
`doc/literature_references.md:40` records it and CLAUDE.md's routing table requires reading the
two together.

### Undocumented module

`src/numerics/opt` — an optimizer with NLOPT backing — has no `doc/` directory and appears in
neither `doc/README.md` nor CLAUDE.md's architecture section. It is live code:
`src/numerics/CMakeLists.txt:5` adds the subdirectory, and `USE_NLOPT` defaults ON
(`CMakeLists.txt:93`).

## Open Questions

- **The φ-region time derivative.** Is the static Laplace form in `maxwell_usage_guide.md`
  a stale document or a deliberate modelling choice? Physics — Christian's call.
- **The Messe attribution.** How should the static-condensation citation read, given the paper's
  own implementation used Lagrange multipliers? Christian's paper, Christian's call.
- **Should `check_doc_claims.py` grow probes** for the `USE_DEBUG` default (it cannot see through
  `${BELFEM_DEFAULT_USE_DEBUG}`), the built-executable set (it inventories `src/executables/*.cpp`
  stems, not CMake targets), and the lessons-learned counts? All three classes of drift found
  tonight are mechanically checkable.

## Files Updated

None. This sweep was read-only by construction; every finding is recorded here and in
`todo/doc_currentness_fixes.md` for a later editing session.

## Status

**Reviewed, not verified.** No executable gate ran. The `hphirun` finding rests on the contents
of a build tree, which outranks a source trace but is still short of a run.
