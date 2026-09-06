# CLAUDE.md stale-content sweep

**Date:** 2026-08-11
**Purpose:** Record the removal of the proprietary literature inventory from the open-source instruction file, the retirement of the `paperN` citation aliases, and the nine factual defects the sweep and its two auditors turned up.

---

## 1. What triggered it

The `literature/` file tree was reproduced inside `CLAUDE.md` — a proprietary, separately-versioned
repository's inventory sitting in the open-source instruction file. It had also drifted: it listed
`felippa` under `books/` (moved to `lectures/` since), and omitted `lectures/`, `papers/topology/`,
`books/gross.*` and `books/cormen.*` entirely. `literature/README.md` already carries a current
structure tree, so the copy was pure duplication of something the open repository should not state.

The `paperN` aliases (`paper1`, `paper6`, `F0`, …) were retired in the same pass: they are opaque to
anyone without the proprietary library, while author-year citations stay useful to every reader.

## 2. What was removed or rewritten

- The `literature/` structure block is gone; the section now states the copyright position and points
  at `literature/README.md` for layout and routing. Example `.txt` filenames were genericized so the
  file no longer names extraction files at all.
- Every `paperN` and `F0`-style alias is gone from the routing table, the implementation-choices
  table, the pitfalls table, the citation-format block, and the key-papers list. The only surviving
  mentions are the note that the aliases are **retired** and that `doc/literature_references.md`
  decodes the ones in older records.
- The key-papers list was brought up to date with the current library: the **Arsenault et al. 2026
  erratum** (it corrects the 2023 air-domain coupling), Denis et al. 2026, Lucchini 2025, Dular et al.
  1997/1999 (which is what the `circuit/` module rests on), and the topology set.
- The Quick Start guide's duplicated naming/container/memory/error blocks and its second copy of the
  implementation-choices table were dropped; navigation and the literature list stay.

## 3. Factual defects found — the part worth reading

The sweep was audited read-only by Codex and Grok. Nine claims in the file were false against the
tree, and each was confirmed at file:line before being touched. Four of them would have produced
wrong code if followed literally:

| Claim in the file | Reality |
|---|---|
| Blaze is the Linux default backend | `CMakeLists.txt:61-67` — **Armadillo** is the default everywhere except Apple. (This one was introduced by the sweep itself and caught by the audit.) |
| "no `share` overload exists for `Cell`" | `share( Cell<T>& )` is at `commtools.hpp:1648`, `receive` at `:549`. Only `Matrix` lacks `share`. |
| `aligned_alloc(64, …)` shown as the SIMD pattern | `src/` contains no `aligned_alloc`/`posix_memalign` call. Alignment is backend-provided and never 64 (Blaze: build SIMD width, 32 under AVX2; Armadillo: 32 above 1 KiB, never 64). |
| `BELFEM_ERROR` for "convergence failures" | `coding_philosophy.md:621` — non-convergence is a third category that returns status and triggers the retry path; aborting there is called a design error. |

And five that were merely wrong rather than dangerous: release flags were stated as
`-O3 -DNDEBUG -fno-exceptions` (actually `-O2 -DNDEBUG`, `-Og -g` in debug); `make tests` does not
exist (`make check` / `check-fast`, and `USE_TEST` defaults OFF); the executables were listed as
"dipole, helix, magnet" (actually `hphirun`, `hphiTrun`, `electricalCircuit`); `fvm/` was listed as a
normal module though `src/CMakeLists.txt` keeps it commented out, and `visualizer/` is `USE_VTK`-gated
OFF; "both static and shared library builds supported" — project libraries are STATIC only, the sole
SHARED targets being the user-plugin templates. The claim that every module README carries DOI-backed
external references held for **1 of 23**.

The audit also caught two omissions with operational weight. The **nonfree AI-exchange ban** (protocol
§2 — the exchange is routed through external vendors, so proprietary work must not touch it) was
absent from `CLAUDE.md` while its checklist told sessions to *use* the exchange. And the
"reviewed ≠ verified" evidence ladder (§11) had no presence in the bootstrap file at all.

## 4. Two edits outside `CLAUDE.md`

`doc/ai_collaboration_protocol.md:171` mandated the `paperN` citation format *and* named `CLAUDE.md`
as its definition. Since the protocol outranks `CLAUDE.md` in the precedence order, retiring the
aliases in one file alone would have left the rule inoperative — so that line and its two example
strings were moved to author-year.

`doc/coding_philosophy.md:525` carried the same false "no `share` overload for `Cell`" row that
`CLAUDE.md` had copied from it; both are corrected, with the `commtools.hpp` line numbers cited so the
next reader can check rather than trust.

## 5. Closed the same day (DR-60, DR-61)

Both gaps this sweep opened were closed on Christian's go-ahead, in the same session.

**DR-60 — the public citation list.** `doc/literature_references.md` is what a reader without the
proprietary `literature/` repository has, and it carried no entry for the seven topology papers, for
Gross & Kotiuga, for CLRS, or for the Evans 2017 lecture notes — precisely the references behind the
`homology/` module. It now has a **Computational Topology and Optimal Cuts** section split three ways
(theory and method: Gross & Kotiuga 2004, Pellikka et al. 2013, Mrozek & Batko 2009; optimality and
complexity: Dey et al. 2011, Chen & Freedman 2011, Dunfield & Hirani 2011, Costantini 1998, Haken
1961; algorithms: CLRS), each entry with a DOI or ISBN and a line on what it answers — including the
negative result, that minimal cuts over Z₂ are NP-hard to approximate, which is worth meeting before
proposing a "just minimize it" cut algorithm. Evans 2017 went into a new **Lecture Notes**
subsection, with Felippa moved beside it. **One claim in this session's own todo file was wrong and
is corrected there:** it said *two* lecture sets were missing. Felippa was already listed under
Advanced References; only Evans was absent — which is what Codex's audit had actually said.

**DR-61 — the alias conversion.** 163 further occurrences across 18 documents, on top of the 14 code
comments: `src/**/doc/*.md` (15 files, led by `dof_manager_usage_guide.md` at 37 and
`maxwell_usage_guide.md` at 23), `doc/input_file_reference.md` (5), and `doc/coding_philosophy.md:100`
(1) — that last one a worked example of the citation style, so the philosophy document had been
*teaching* the retired form. Two rules, applied by script and then read line by line in the diff:
delete the alias where an author-year already precedes it, expand a bare alias to the full citation.
**`paper5` and `paper6` are why a blind delete would have been wrong** — both read "Alves et al.
2022" in the surrounding prose, so deleting the alias would have left an ambiguous citation; they
became `2022a` and `2022b` explicitly. Seven parentheticals the substitution left reading
redundantly (`**Messe et al. 2023** - BELFEM core (Messe, SUST 2023)`) were tidied to the journal
alone. Scope ruling recorded as O1: `todo/` and `devlog/` are left alone, and the alias tables in
`doc/literature_references.md` are kept — relabelled "retired — decoder for older records" — because
they are what keep the dated entries readable. Measured after: `src/` at zero, `doc/` at 15, all of
them inside the decoder tables and the prose explaining the retirement.

**46 of those 163 were in `doc/literature_references.md` itself**, whose body turned out to be written
in alias vocabulary throughout — entry headings, the question → literature mapping, the pitfalls
table. The FVM `F0`–`F4` aliases went the same way, and they are the better argument for retirement
than `paperN` ever was: the 2026-08-07 jury devlog uses `F1` and `F2` as *finding* IDs in the same
sentences where other entries use them as papers.

**And converting that body turned up a defect nobody had flagged.** Putting the entry headings beside
the alias table showed that **ten of the eleven paper titles disagreed with the DOI-backed citations
twelve lines above them** — the body called Messe et al. 2023 "BELFEM: A Finite Element Framework for
HTS Applications" where the published title is "BELFEM: a special purpose FE code for magnetodynamic
modeling of HTS tapes", and Schnaubelt et al. 2023 "3D Thin-Shell h-φ for No-Insulation Coils" against
"Electromagnetic Simulation of No-Insulation Coils Using H-φ TSA". They read as titles, not
descriptions, and this document exists so that a reader *without* the proprietary library can find the
works — a paraphrased title is a search that fails. All ten now carry the published title. Logged as
D11; it was invisible until the conversion put the two forms side by side.

## 6. The checker (DR-62) — built, and it found bugs in itself first

Christian ruled to build the mechanical check, **against my recommendation**, which had been to wait
and fold the probes into the input-contract checker nobody has written yet. The ruling was the better
call, for a reason that only appeared during construction.

`scripts/check_doc_claims.py` verifies 21 claims across `CLAUDE.md` and `doc/coding_philosophy.md`,
each anchored on a searchable token rather than a line number — the input contract's own rule, for the
same reason. It probes the optimisation flags, the `-fno-exceptions` design-goal-versus-build-fact
distinction, make targets and the `USE_TEST` default, the `src/executables` inventory, modules
commented out of `src/CMakeLists.txt`, the backend defaults either side of `if( APPLE )`,
STATIC-versus-SHARED libraries, hand-rolled aligned allocation, every `share( … )` overload denial
against `commtools.hpp`, and the error-tier claim about convergence failure. Prose and rationale are
deliberately out of scope: those need a reader, which is what the cross-review round is for.

**A green run on documents I had just fixed would prove nothing**, so the real test was to run it
against the pre-sweep `CLAUDE.md`, recovered from `HEAD`. It reports 13 failures there and
rediscovers, unaided, `-O3`, `-fno-exceptions`, `make tests`, the `USE_TEST` default, the executable
list, `fvm/`'s build status, the backend inversion, "static and shared", the `aligned_alloc` snippet
and the convergence-failure tier. Current tree: 21/21, exit 0.

**Three bugs surfaced in the checker itself, all the same class, and it is the class that matters
here.** First it split `config_gcc.cmake` at the wrong `else()` — an unrelated `BELFEM_USE_CLANG`
block precedes the `USE_DEBUG` one — then, after that fix, a non-greedy `endif()` match was cut short
by the nested `if( APPLE )`. **Codex's audit found the third:** the parser matched only `if(`, never
the `if (` spelling, which exists at `config_gcc.cmake:66` and `CMakeLists.txt:101,231` — so an
`endif` could decrement a depth that was never incremented. Correct for `USE_DEBUG` today, wrong the
moment it is pointed anywhere else. In every case a probe returned an *empty* result, so the check
reading it passed **vacuously while appearing to run**: the flag assertions were no-ops reporting
success. Neither of the first two would have been visible without the pre-sweep regression run.

**The structural fix matters more than the three bugs.** `REQUIRED_FACTS` now names every fact a
check depends on, and an empty one is a *failing row* rather than a silently dropped check — Codex
demonstrated the general case by feeding empty facts into `checks()` and getting zero failures with
fewer rows, which is the same silent success wearing a different hat. Its six other findings are all
applied: `-fno-exceptions` was asserted in a detail string but never probed; the library-kind check
read only `Add_Library.cmake` although real `SHARED` targets exist in the two user-plugin templates,
so the probe now scans the tree and requires the document to name that exception; `forbid()`'s
negation window looked only *behind* the match, so "`-O3` is never used" passed while "`-O3` is not
used" would have been failed as an error; and `vtk_gated` was computed but consumed by nothing.

Validated four ways: 32/32 on the current tree; 15 failures on the pre-sweep `CLAUDE.md`; 32/32 on a
fixture carrying both negation wordings, so no false positive; and a fixture with the compiler config
removed fails its probe rows instead of passing.

Not wired to anything — no hook, no CI, since there is none. It is run by hand; `CLAUDE.md` says when.

## 7. Left open

- **DR-63** — 22 of 23 module READMEs carry no external reference, against the stated convention:
  populate them, or amend the convention. A direction is owed, not code.
- **O3, from Codex's audit** — `doc/coding_philosophy.md:621` makes non-convergence a status/retry
  outcome, but `src/fem/interpolation/cl_IF_ElementMapper.cpp:369` aborts on it:
  `BELFEM_ERROR( tCount < 100, "ElementMapper::evaluate_general: failed to converge." )`. This is a
  fair challenge to the *rule*, not only to the site. It may well be correct — the philosophy's rule
  is about algorithms with a retry policy above them, and a point-inversion Newton inside element
  mapping has none; a failure there means the point is not in the element, which is a caller bug, not
  a recoverable numerical state. If that reading holds, the philosophy should say so, because as
  written the rule reads as unconditional. Source untouched; Christian's ruling.
- Nothing in this session was compiled or run. Every change is documentation or comment text, but the
  corrected build, container, and MPI statements were each read out of the source they describe, and
  the bibliographic fields added in DR-60 were copied from the `literature/` navigation headers rather
  than recalled.
