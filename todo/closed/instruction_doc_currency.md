# Instruction-Document Currency: Keeping CLAUDE.md and the Convention Docs True

**Date:** 2026-08-11
**Purpose:** The 2026-08-11 `CLAUDE.md` sweep found nine claims that were false against the tree,
four of which would have produced wrong code if followed. All nine are fixed. This plan records what
was found, closes the residual gaps the sweep opened (an incomplete public citation list, legacy
citation aliases in living documents), and answers the question the sweep raised but did not settle:
whether a mechanical code → document check is worth building for the convention documents, the way
the input contract already has one for `input.conf` keys.
**Module:** meta (documentation convention)
**AIs involved:** Claude (sweep + fixes), Codex (audit), Grok (audit)
**Status:** IN PROGRESS — D1–D9 and D11 fixed (D10 retracted), **R1–R4 all DONE (2026-08-11)**,
O1 and O2 resolved. `DR-60`, `DR-61`, `DR-62` **and `DR-63`** are closed;
`scripts/check_doc_claims.py` guards a probe-derived claim set (34 as of 2026-08-14 — the count
floats with the tree: the fvm move to nonfree retired its build-status probe, and the USE_TEST
default check became bidirectional, ON and OFF, after the one-way form went silent when the
default flipped to ON). R4 was closed by *amending the
convention* rather than manufacturing 22 reference sections — `doc/documentation_guidelines.md`
now names `doc/literature_references.md` as the single source of truth for citations, says module
docs link to it rather than carry their own DOIs, and records that 1-of-23 is the expected state
(commit `7432385a`).
**The only thing left in this file is O3** — whether `BELFEM_ERROR` on non-convergence is correct
in `ElementMapper::evaluate_general`, i.e. whether `coding_philosophy.md`'s rule should be narrowed
or the site changed. That is Christian's ruling, not work. Documentation and one new script;
nothing compiled, and the checker is run by hand.

> **Scope guards:**
> - Devlogs are **out of scope** for any alias conversion — they are a dated historical record.
> - `doc/literature_references.md`'s alias index is **kept**, deliberately: it is what makes the
>   citations in older records readable.
> - No source behaviour is in scope. The only `src/` change in this campaign was comment text plus
>   one error-message string.

---

## 1. What the sweep was, and why the findings matter

`CLAUDE.md` is read at the start of every AI session and is the only document many sessions read in
full. A false statement in it does not sit inert — it is acted on. Four of the nine defects were of
exactly that kind: they instructed a session to write code that the codebase forbids elsewhere.

The document had also been reproducing the `literature/` directory's file inventory, which is a
separate, proprietary repository. That was the trigger for the sweep (Christian, 2026-08-11); the
defects below were found while removing it.

## 2. Defects (all fixed 2026-08-11)

- [x] **D1 — HIGH. The open-source instruction file reproduced the proprietary literature inventory.**
  A full directory tree of `literature/` sat in `CLAUDE.md`, in a repository where that directory is
  `.gitignore`d for copyright reasons. It had drifted besides: `felippa` was still listed under
  `books/` though it moved to `lectures/`, and `lectures/`, `papers/topology/`, `books/gross.*` and
  `books/cormen.*` were absent entirely.
  *Fixed:* tree deleted, replaced by a pointer to `literature/README.md`, which already carries a
  current one. Example extraction filenames were genericized, so the open repository now names no
  file inside `literature/` at all.

- [x] **D2 — CRITICAL. "No `share` overload exists for `Cell`" — false.**
  `share( Cell< T > & )` is at `src/comm/commtools.hpp:1648`, with `receive( Cell< T > & )` at `:549`.
  Only `Matrix` genuinely lacks `share`. A session following the table would have reached for an
  unchunked `broadcast` on exactly the large payloads that `broadcast` is documented to fail on.
  The same false row was in `doc/coding_philosophy.md:525` — `CLAUDE.md` had copied it from there,
  so this is one defect in two documents, not two findings.
  *Fixed:* both tables corrected, with the `commtools.hpp` line numbers cited in the philosophy doc
  so the next reader can check rather than trust. Found by Grok.

- [x] **D3 — HIGH. `aligned_alloc(64, …)` taught as the SIMD allocation pattern.**
  `src/` contains no `aligned_alloc` or `posix_memalign` call anywhere, and no backend delivers 64
  bytes (Blaze aligns to the build's SIMD width, 32 under AVX2; Armadillo to 32 above 1 KiB and
  deliberately never 64). `doc/coding_philosophy.md:547-555` had already identified this snippet as
  defective and removed it; `CLAUDE.md` kept teaching it.
  *Fixed:* replaced with the alignment reality — allocate through `Vector`/`Matrix`, let the backend
  align. The stale "need explicit alignment (64-byte for SIMD)" bullet in the manual-memory rules
  went with it. Found by Grok, confirmed by Codex.

- [x] **D4 — HIGH. `BELFEM_ERROR` prescribed for convergence failures.**
  `doc/coding_philosophy.md:621-623` makes non-convergence a third category that returns a status and
  triggers the retry path above it, and calls `BELFEM_ERROR("did not converge")` a design error
  because it converts a recoverable state into a run abort.
  *Fixed:* removed from the `BELFEM_ERROR` list; the third category is now stated in `CLAUDE.md`
  with the retry-path rationale. Found by Grok.

- [x] **D5 — HIGH. Matrix backend defaults inverted.** `CMakeLists.txt:61-67`: Armadillo is the
  default everywhere **except** Apple, where Blaze is. **Introduced by the sweep itself** and caught
  by the audit round — the pre-sweep text said only "Armadillo or Blaze" and was not wrong.
  *Fixed:* corrected. Recorded here because it is the sharpest argument for auditing a documentation
  rewrite the same way a code change is audited: the sweep that removed eight falsehoods added one.

- [x] **D6 — MEDIUM. Six build and architecture facts wrong.** Release flags given as
  `-O3 -DNDEBUG -fno-exceptions` (actually `-O2 -DNDEBUG`; debug is `-Og -g`,
  `config/compiler/config_gcc.cmake:15-32`); `make tests` as the test command (it is `make check` /
  `check-fast`, and `USE_TEST` defaults OFF, `CMakeLists.txt:81,288-332`); executables listed as
  "dipole, helix, magnet" (actually `hphirun`, `hphiTrun`, `electricalCircuit`); `fvm/` listed as a
  normal module though `src/CMakeLists.txt:12` keeps it commented out; `visualizer/` listed
  unconditionally though it is `USE_VTK`-gated, default OFF; "both static and shared library builds
  supported" though `config/scripts/Add_Library.cmake:27` builds every project library STATIC and the
  only SHARED targets are the user-plugin templates.
  *Fixed:* all six. `circuit/`, `fvm/`, `visualizer/`, `fem/thermal/` and `math/quaternion` were also
  missing from the architecture list. Found by Claude and Codex independently.

- [x] **D7 — MEDIUM. The protocol mandated a citation format that `CLAUDE.md` forbids.**
  `doc/ai_collaboration_protocol.md:171` required "Author et al. YYYY (paperN), Section X" *and*
  named `CLAUDE.md` as the definition of that format. Since the protocol outranks `CLAUDE.md`
  (precedence order, §"Instruction Precedence"), retiring the aliases in `CLAUDE.md` alone would have
  left the rule inoperative — an AI could not satisfy both documents.
  *Fixed:* that line and its two example strings moved to author-year. Found by Grok.

- [x] **D8 — HIGH. The nonfree AI-exchange ban was absent from `CLAUDE.md`.**
  Protocol §2 suppresses the `tmp/ai_exchange/` tier entirely for `./nonfree/` work, because the
  exchange is routed through external AI vendors. `CLAUDE.md`'s compliance checklist instructed
  sessions to *use* the exchange, with no exception — so a session that read only `CLAUDE.md` would
  have written proprietary findings into a vendor-bound channel.
  *Fixed:* the exception is now in both the key points and the checklist. Found by Grok.

- [x] **D9 — LOW. "Each module README provides external references with DOIs" — false for 22 of 23.**
  Only `src/homology/doc/README.md` carries one.
  *Fixed:* the claim is now hedged and points sessions away from expecting a literature trail there.
  The underlying gap is R4. Found by Codex.

- [x] **D11 — MEDIUM. Ten paper titles in the public reference list were paraphrases, not titles.**
  Found while converting `doc/literature_references.md` for R2: its body entries carried titles that
  disagreed with the DOI-backed citations in its own alias table twelve lines above — e.g.
  "BELFEM: A Finite Element Framework for HTS Applications" against the published "BELFEM: a special
  purpose FE code for magnetodynamic modeling of HTS tapes", and "3D Thin-Shell h-φ for No-Insulation
  Coils" against "Electromagnetic Simulation of No-Insulation Coils Using H-φ TSA". All ten Tier 1–4
  entries were affected. In a document whose stated purpose is letting a reader without the
  proprietary library find the works, a paraphrased title is a search that fails.
  *Fixed:* every body entry now carries the DOI-backed title from the alias table. Nobody flagged
  this — it surfaced only because the conversion put the two forms side by side.

- [x] **D10 — FALSE POSITIVE (retracted 2026-08-11).** Codex reported a surviving literature alias at
  `src/fem/kernel/cl_FEM_Calculator.cpp:287`. The `( F4 )` there is a review finding ID, not the
  Wheeler 2011 FVM alias — Wheeler concerns mixed FVM on hexahedra and has nothing to say about
  thin-shell facet machinery on buffer blocks. Pattern matched, meaning did not. Recorded so the
  next alias sweep does not re-flag it.

**Also landed 2026-08-11 (not defects):** the `paperN` aliases were retired across `CLAUDE.md` and
converted in all 14 `src/**.cpp/.hpp` comment sites (`cl_EdgeFunctionFactory.cpp` ×2,
`cl_FEM_Controller.cpp` ×10, `cl_FEM_Controller.hpp` ×2, `cl_FEM_DofMgr_SolverData.cpp` ×1), all of
which were `paper1` = Messe et al. 2023. One of those is a user-facing `BELFEM_ERROR` string; no
executable line changed. `CLAUDE.md` also gained the "reviewed ≠ verified" evidence-ladder line
(protocol §11) and the current key-paper set, including the **Arsenault et al. 2026 erratum** that
corrects the 2023 air-domain coupling.

## 3. Remaining steps

- [x] **R1 — Complete `doc/literature_references.md`. DONE 2026-08-11.** It is the public fallback
  for readers without the proprietary library, and it covered the FEM papers, the FVM papers and the
  FEM textbooks — but **not** the seven `papers/topology/` works, `books/gross` (Gross & Kotiuga),
  `books/cormen` (CLRS), or the Evans 2017 lecture notes. These are the references behind the
  `homology/` module, so the gap was not cosmetic.
  *Correction to this plan's own first draft:* it claimed **two** missing lecture sets. Felippa was
  already listed, under "Advanced References" — only Evans was absent. Codex's audit had said exactly
  this ("only surfaced Felippa among those lecture/topology names") and the plan overstated it.
  *Landed:* a new **Computational Topology and Optimal Cuts** section (theory and method → Gross &
  Kotiuga, Pellikka et al. 2013, Mrozek & Batko 2009; optimality and complexity → Dey et al. 2011,
  Chen & Freedman 2011, Dunfield & Hirani 2011, Costantini 1998, Haken 1961; algorithms reference →
  CLRS), each with DOI or ISBN; a **Lecture Notes** subsection holding Evans 2017 and Felippa (moved
  there from Advanced References); a fourth category in the Overview; a routing block for cuts and
  cohomology questions; and topology search tags. `CLAUDE.md`'s scoping hedge is removed.
  Bibliographic metadata came from each work's `.md` navigation header — citations, not extracted
  text.

- [x] **R2 — Alias conversion in living documents. DONE 2026-08-11** *(scope per O1)*. **163
  occurrences across 18 files.** 117 of them in `src/**/doc/*.md` (15 files — `dof_manager_usage_guide.md`
  37, `maxwell_usage_guide.md` 23, `hanging_dofs_static_condensation.md` 13,
  `nonlinear_controller_theory.md` 9, and 11 others), `doc/input_file_reference.md` (5), and
  `doc/coding_philosophy.md:100` (1) — that last one a *worked example of the citation style*, so it
  had been actively teaching the retired form.
  The remaining 46 were in **`doc/literature_references.md` itself**, which turned out to be written in
  alias vocabulary throughout its body — entry headings (`**[paper1]** …`), the question → literature
  mapping, the pitfalls table, the implementation-guideline headings. The decoder tables at the top are
  kept and now cover both alias families, so the body could be converted without orphaning any older
  citation; the FVM `F0`–`F4` aliases were converted the same way, since devlogs use those too (and
  ambiguously — `F1`/`F2` also appear as *finding* IDs in the 2026-08-07 jury entry).
  **`src/` and `src/**/doc/` are now at zero**; the 15 alias tokens left in `doc/` are the two decoder
  tables, the prose explaining the retirement, and the protocol's retirement note.
  Two conversion rules, applied by script and then reviewed line by line: where an author-year already
  preceded the alias the alias was deleted (`Messe et al. 2023 (paper1)` → `Messe et al. 2023`), and a
  bare alias became the full citation (`From paper1:` → `From Messe et al. 2023:`). **`paper5` and
  `paper6` were the reason a blind delete would not do:** both read "Alves et al. 2022" in the
  surrounding prose, so dropping the alias would have produced an ambiguous citation — they became
  `2022a` and `2022b` explicitly. Seven redundant parentheticals left behind by the substitution
  (`**Messe et al. 2023** - BELFEM core (Messe, SUST 2023)`) were tidied to the journal alone.
  `todo/` and `devlog/` untouched, per O1 and the scope guard.

- [x] **R3 — DONE 2026-08-11: `scripts/check_doc_claims.py` built** *(after: O2, decided by Christian)*.
  21 claims across `CLAUDE.md` and `doc/coding_philosophy.md`, each anchored on a searchable token
  rather than a line number, per the input-contract rule. Probes: release/debug optimisation flags
  (`config/compiler/config_gcc.cmake`), the `-fno-exceptions` design-goal-vs-build-fact distinction,
  make targets and the `USE_TEST` default, the `src/executables` inventory, modules commented out of
  `src/CMakeLists.txt`, the matrix-backend defaults either side of `if( APPLE )`, STATIC-vs-SHARED
  libraries, hand-rolled aligned allocation, every `share( … )` overload denial against
  `commtools.hpp`, and the error-tier claim about convergence failure.
  **Validated in both directions:** 21/21 on the current tree (exit 0), and **13 failures against the
  pre-sweep `CLAUDE.md` recovered from `HEAD`** (exit 1) — it independently rediscovers `-O3`,
  `-fno-exceptions`, `make tests`, the `USE_TEST` default, the executable list, `fvm/`'s build status,
  the backend inversion, "static and shared", the `aligned_alloc` snippet and the convergence-failure
  tier. A green run on documents that were just fixed proves nothing; that second run is the evidence.
  **Three bugs were found in the checker itself, all the same class:** a probe returned an empty
  result, so the check reading it passed *vacuously* while appearing to run. Two were mine, found
  while building it — splitting `config_gcc.cmake` at the wrong `else()` (an unrelated
  `BELFEM_USE_CLANG` block precedes the one wanted), then a non-greedy `endif()` cut short by the
  nested `if( APPLE )`. The third came from **Codex's audit**: the parser matched only `if(` and not
  the `if (` spelling, which exists at `config_gcc.cmake:66`, `CMakeLists.txt:101` (`endif ()`) and
  `:231` (`elseif`) — so an `endif` could decrement a depth that was never incremented. Correct for
  `USE_DEBUG` today, wrong the moment it is pointed anywhere else.
  **The structural fix matters more than the three bugs:** `REQUIRED_FACTS` now names every fact a
  check depends on, and an empty one is a **failing row**, not a skipped check. Codex demonstrated the
  general case by feeding empty facts into `checks()` — every scenario produced zero failures with
  fewer rows, i.e. silent success. A checker whose failure mode is reporting success is worse than no
  checker.
  **Codex's other findings, all applied:** `-fno-exceptions` was asserted in the detail string but
  never probed (now grepped from the CMake files); the library-kind check read only
  `Add_Library.cmake` although real `SHARED` targets exist in `UserLibraryTemplate.cmake:92` and
  `UserMaterialTemplate.cmake:96`, so the probe now scans the tree and *requires* the document to name
  that exception; `forbid()`'s negation window looked only behind the match, so "`-O3` is never used"
  passed but "`-O3` is not used" would have false-positived — both sides are checked now; and
  `vtk_gated` was a computed-but-unused probe, now consumed as a real row.
  **Validated three ways:** 32/32 on the current tree; **15 failures on the pre-sweep `CLAUDE.md`**;
  32/32 on a fixture whose philosophy document carries both negation wordings (no false positive); and
  a fixture with the compiler config removed fails its probe rows instead of passing.
  *Not done:* the checker is not wired to anything — no hook, no CI (there is none). It is run by hand;
  `CLAUDE.md` says when.

- [x] ~~**R3 (original wording) — Decide whether a mechanical doc → tree check is worth building**~~
  *(kept for the reasoning; the decision was **yes**, ruled by Christian 2026-08-11, and the entry
  above is what was built.)*
  The input contract already has this shape of check in the **code → schema** direction, and
  `CLAUDE.md` explains why: the failure mode is not neglect, it is that whoever adds the key is
  precisely the person unaware the document exists. The nine defects above are the same failure mode
  in a second place. What distinguishes them is that several are *mechanically checkable*: compiler
  flags, `make` targets, executable names, module presence in `src/CMakeLists.txt`, backend defaults,
  and whether a named symbol (`share( Cell )`, `aligned_alloc`) actually exists. A check that greps
  the claims out of the convention documents and diffs them against the tree would have caught
  D2, D3, D5 and D6 — six of the nine — with no judgement involved.

- [x] **R4 — Module README external references. DONE 2026-08-11** (`7432385a`, DR-63 closed);
  ticked in the 2026-08-11 currentness sweep after verifying the amendment in the tree. 22 of 23
  `src/*/doc/README.md` carried no DOI or external citation, against a documented convention that
  said they should. **Resolved by amending the convention, not by manufacturing 22 reference
  sections:** `doc/documentation_guidelines.md:312` now states that
  `doc/literature_references.md` is the single source of truth for full citations and DOIs, that
  module docs *link* to it, and that a "References" heading with nothing under it is worse than its
  absence — so modules implementing no published method (`core`, `containers`, `comm`, `io`) are
  expected to cite nothing and 1-of-23 is the expected state, not a gap.

## 4. Open questions

- [x] **O1 — RESOLVED 2026-08-11 → convert `src/**/doc/` and `doc/`, leave `todo/` and `devlog/`.**
  Christian delegated the ruling with the instruction to handle DR-61; the recommendation below was
  applied as written. Rationale: `src/**/doc/` and `doc/` are living reference material read by people
  who may not have the proprietary library, and `coding_philosophy.md:100` was teaching the retired
  form outright. `todo/` files are rewritten by the campaigns that own them, and `devlog/` is a dated
  record. `doc/literature_references.md`'s alias table is kept deliberately as the decoder, now
  relabelled "retired — decoder for older records".

- [x] **O2 — RESOLVED 2026-08-11 → build it. Ruled by Christian**, against Claude's recommendation,
  which had been to fold the checks into the (still unwritten) input-contract checker rather than
  build anything standalone. The ruling is the better call on the evidence: the checker found two
  vacuous-pass bugs *in itself* during construction, which is the same class of defect it exists to
  catch, and waiting on a tool nobody has scheduled would have left the class unguarded indefinitely.
  Evidence that argued for it: the sweep's own D5 shows a careful rewrite still injects errors, and a
  grep-level check is cheap. Evidence that argued against, and still stands: there is no rate estimate
  — the nine defects accumulated over the document's whole life, not per month, so the checker's value
  is unproven until it catches something nobody was looking for.

- [ ] **O3 — Is `BELFEM_ERROR` on non-convergence correct in `ElementMapper::evaluate_general`?**
  Raised by Codex while auditing the checker (2026-08-11), and it is a fair challenge to the *document*
  rule, not just to the tool: `doc/coding_philosophy.md:621` makes non-convergence a status/retry
  outcome, yet `src/fem/interpolation/cl_IF_ElementMapper.cpp:369` aborts —
  `BELFEM_ERROR( tCount < 100, "ElementMapper::evaluate_general: failed to converge." )`.
  **This may well be correct:** the philosophy's rule is about algorithms that *have a retry policy
  above them*, and a point-inversion Newton inside element mapping has none — a failure there means
  the point is not in the element, which is a caller bug rather than a recoverable numerical state.
  If that reading holds, the philosophy should say so explicitly, because as written the rule reads
  as unconditional. Claude did not touch the source. Needs Christian's ruling: narrow the rule, or
  change the site.

## 5. Verification status

**Nothing here was compiled or run.** The claims corrected in the documents were each read out of the
source they describe, at the file:line cited above — that is a *static source trace* on the evidence
ladder (protocol §11), not a verified result. The `src/` comment conversions are text-only and cannot
change behaviour, but they have not been through a compiler.

What *was* mechanically checked, for R1 and R2: every bibliographic field added in R1 was copied from
the corresponding `literature/**/*.md` navigation header rather than recalled; the R2 substitution was
run as a script and then reviewed line by line in the diff, and the residual counts were measured, not
assumed — `src/` and `src/**/doc/` at zero, `doc/` holding only the decoder table (60) and the
protocol's retirement note (1).

## 6. Threads

- `tmp/ai_exchange/claudemd_sweep.md` — Codex audit (12 findings)
- `tmp/ai_exchange/claudemd_sweep_grok.md` — Grok audit (16 findings + 12 omissions)
- `devlog/dl20260811_claudemd_stale_content_sweep.md` — session record
