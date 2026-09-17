# Commenting Guidelines {#doc_commenting_guidelines}

**Date:** 2026-09-16
**Purpose:** What a source comment in BELFEM may say, what it must not say, where it goes, and how it is written — so that the comments that carry a contract, a reason, or a citation are not buried under the ones that do not
**Module:** cross-cutting (C++ and Fortran)

**Status:** proposed; revised after one blind jury round (Codex and Grok). Consolidated from four independent drafts and the language agreements of 2026-09-15. Once adopted, the rules apply to new code at once and to existing code through the cleanup sweeps planned in `todo/`. The census figures quoted below are approximate counts from `src/` on 2026-09-16; they motivate the rules and are not part of them.

---

## 1. The one rule

**The code is the only authoritative description of what the code does. A comment may only add what the code cannot say.**

All four books this document rests on arrive at that rule from different directions. Martin says the proper use of a comment is "to compensate for our failure to express ourselves in code" and asks that the code be made to say it first. Oliveira & Stewart put it as *once and only once*: there is exactly one authoritative description, the code itself, and a comment that restates it becomes a lie the moment the code changes (`/* increment i */` above `i = i - 1`). Thomas & Hunt (Topic 9, "Duplication in Documentation") reach the same place from DRY: a comment that repeats the code states the intent twice, and "given time, we can pretty much guarantee the comment and the code will get out of step"; their Tip 13 restricts non-API comments to "*why* something is done, its purpose and its goal", because the code already shows *how*. Rouson, Xia & Xu open their style rules with "write code that comments itself" and note that a named constant "serves dual roles as declaration and effectively as comment".

The corollary is the cost model. A comment is not free because it is not compiled. It costs the reader's time now, and it costs false confidence later, when it has drifted. **A comment earns its place when what it says is not derivable from the code and will still be true after the next three edits.**

Two consequences that matter for this codebase:

- **Density hides signal.** When one comment accompanies every second statement, the reader learns to skip comments, and then skips the one that mattered. A large share of the `//` lines in `src/` today, on the order of a third, narrate the next statement. The comments that carry a real fact (`cl_FEM_Controller.cpp` on which reset undoes which shift, `commtools.hpp` on why a matrix is sent as `spacing() * n_cols`) sit inside that noise. Deleting the narration is what makes them visible.
- **Comment density follows the difficulty of the problem, not the length of the code.** Rouson et al. raise their comment frequency when the *problem* grows complex (§7.2.1), not when the file does. A 200-line assembly kernel that mirrors a published equation may need one citation. A 20-line index-space remapping may need a paragraph.

---

## 2. What a comment may say

Each kind below adds something the code cannot. The test: **could a competent reader reconstruct this from the code they are looking at?** For a comment inside a body, that means the surrounding statements. For a contract at a declaration, it means the declaration alone: a caller reads the header, not the body, so an ownership rule or a precondition that is obvious from the implementation still passes the test at the interface. If the answer is yes, the comment is not one of these kinds, and it goes.

| Kind | What it carries | Example in the tree (anchor is a greppable token, not a line) |
|---|---|---|
| **Contract** | precondition, postcondition, units, required layout (column-major, sorted, 1-based), ownership, what happens on the degenerate case, MPI collectiveness — design by contract in comment form (Thomas & Hunt, Topic 23) | `cl_Mesh_ConnectivityCalculator.hpp`: `nothing to release; all containers are owned by the mesh` — the parenthetical is the contract; the `@brief Destructor` around it is the form §5 deletes |
| **Reason** | the decision behind a choice that looks wrong or optional, including a one-line name of the alternative that was rejected | `cl_FEM_Controller.cpp`: the thermal algorithm reset marked `( leaked across retries )` |
| **Literature anchor** | which equation or section this block implements | `cl_FEM_Controller.cpp`: `Messe et al. 2023 Eq. 13`; `geometrytools.cpp`: `Eq. 8.152c` of Bronshtein (the source spells it `Bronstein`, a typo the sweep corrects) |
| **Warning** | a trap the reader cannot see: aliasing, buffer reuse, thread-unsafety, a call that is collective, a call order that must be respected, a non-obvious cost | `commtools.hpp`: `Never use capacity() here: a matrix that shrank keeps its old, larger allocation` |
| **Amplification** | "this looks inconsequential and is not" — a `+1`, a reset that must precede a shift, a `trim()` that decides a lookup | the `nan` guard in the controller watchdog (see §8 for its rewrite) |
| **Stage heading** | *what* a multi-line block accomplishes, at a higher level than its statements: `// binary search for i with x(i) <= t < x(i+1)` (Oliveira & Stewart §7.9) | acceptable above a stage in a long kernel; never repeated inside it |
| **Symbol annotation** | the mathematical symbol a statement or a variable stands for, inside a kernel written in domain notation | `// xi` above the affine map of a facet's integration points in `fn_IF_initialize_integration_points_on_facet.cpp`; units on an enumerator, `// [Ω·m]` |
| **Clarification of foreign code** | the meaning of a library argument you cannot rename | the `ICNTL`/`CNTL` slots in `cl_SolverMUMPS.cpp`; the LAPACK work-array contracts in `src/linalg/lapack/` |
| **Invariant** | a statement guaranteed true at this point | prefer `BELFEM_ASSERT` where it can be executed; a comment then explains *why* the invariant holds, not *that* it does |
| **Actionable TODO** | a job with an owner and an exit condition (§7) | — |
| **License header** | the copyright block at the top of every file | the existing block; it and `TODO(owner)` are the only places a person's name belongs |

Doxygen blocks are the contract kind applied to a public interface (§5).

---

## 3. What a comment must not say

| Anti-pattern | Why it hurts | Present in the tree (approximate, 2026-09-16) |
|---|---|---|
| **Narration** — restating the next statement (`// increment counter`, `// delete fields`, `// loop over all blocks`, `// allocate memory`, `// solve system`) | teaches the reader to skip comments; drifts into a lie when the statement changes | on the order of 3 000 lines in `src/` (a classification, not a token count) |
| **Name restatement in Doxygen** (`@brief Constructor`, `@brief Get the material type`, `@param aMesh the mesh`, `return the size of the Cell` above `size()`) | Martin's "mandated comments" (2025, Ch. 5); Thomas & Hunt's "mechanical comment writing" that leaves "two things to update" (Tip 13) | about 90 of the roughly 410 `@brief` blocks in headers |
| **History** — dates, bylines, "changed on", "removed 2026-08-29", "measured 2026-08-21 on tapestack3d", `Decision: <name>, <date>` | git holds the *when* and the *who*; a date in a comment is a journal entry that never gets a follow-up | about 130 dated comment lines in `src/`; 28 of them in `cl_FEM_Controller.cpp` |
| **Provenance of a review** — `( Codex+Grok audit )`, `( Grok R1c )`, `( three-AI round 2026-08-21 )`, `// Created by claude on 1/11/25` | says who agreed, not why the code is right; meaningless to anyone outside the session | about 30 lines in `src/` |
| **Pointers into the working record** — `tmp/ai_exchange/…`, `INC-551`, `DR-77`, `L-18`, a todo *step* ID such as `R3` | the exchange is swept; the incident IDs index a document the reader of the source does not have open; the same rule already bars `doc/` from citing devlogs and todos (Martin 2025, Ch. 5 "Nonlocal Information"). A durable `todo/<task>.md` *file path* is different and is allowed in a TODO (§7) | 3 lines in `src/` |
| **Experiment logs** — iteration counts, timestep traces, residuals from the run that motivated a rule | the rule matters; the run belongs in the devlog that recorded it | the controller watchdog, the STRUMPACK ordering notes, `powerlaws.hpp` |
| **Derivations and tutorials** — a page of the paper above three lines of arithmetic, the user-library walkthrough repeated in three headers | the reader wanted the loop; the derivation belongs in `src/<module>/doc/` or the paper | `cl_Material_UserDefined.hpp` and `cl_MaterialFactory.hpp` carry the same tutorial |
| **Nonlocal information** — a default that this function does not own, a behavior of a distant module | cannot be kept in sync from here (Martin). A member's own default initializer (`= 1e-6 ; // relative error criterion`) is local and is not this | check headers that restate `input.conf` defaults |
| **Commented-out code** | nobody dares delete it; git already remembers it (Martin 2025, Ch. 5 "Commented-Out Code") | about a hundred commented-out statements, plus the two `/* … */`-disabled `save_fields`/`load_fields` bodies in `cl_Mesh.cpp` |
| **Mumbling** — `// help constant`, `// special case`, `// important`, a symbol above a statement it does not describe | means something only to the author | scattered; a symbol annotation in a kernel (§2) is not mumbling and the sweep must tell the two apart |
| **Venting, jokes, apologies** — `// todo: old function, should be obsolete soon` | Martin's "Give me a break!"; an apology is not a plan | the bare `// todo:` remarks in `src/fem` and `src/homology` |
| **Closing-brace labels** (`} // end for`) | a symptom of a function that is too long | rare here; keep it that way |
| **HTML in comments** | unreadable in the editor, which is where comments are read; Doxygen emits its own markup | almost none; keep it that way |

**Misleading is worse than absent.** A comment that is slightly wrong sets an expectation the code will not meet. Four such comments were found in one review pass: a "FACE IDs" heading above a cell-dof loop in `cl_FEM_DofMgr_DofData.cpp`; `wait until receive is complete` above an `MPI_Wait` on a *send* request in `commtools.hpp`; a derivative contract in `cl_JcFunction_Database.hpp` that one branch violates; and a MUMPS symmetry note in `cl_FEM_DofMgr_EigenValues.cpp` that said "lower triangle" where MUMPS accepts either triangle. These are corrected before anything is shortened.

**Not comments at all.** `!$omp` and `!dir$` lines in Fortran and `#pragma` lines in C++ are compiler directives. A sweep that matches comment syntax must exclude them, and so must any census. The same holds for string literals that happen to contain a date or a name: two user-facing error messages in `cl_FEM_Controller.cpp` name the file-format cutoff date, and they are program text, not comments.

---

## 4. Placement and length

**One comment per idea, above the block it describes.** A comment is a heading for the reader. It sits on its own line above the code, not at the end of the line. The exception is a unit or symbol annotation on a declaration or an enumerator (`real tRho ; // [Ω·m]`).

**Above the block, not woven through it.** A single comment beside one surprising statement inside a stage is fine; that is where a local trap or an order dependency belongs. A comment on every statement is narration. If a loop body seems to need prose throughout, ask first whether a named variable, a `const`, or a `BELFEM_ASSERT` would carry the fact instead. If a block needs several headings, it is several blocks.

**Where each kind of information lives:**

| Information | Home |
|---|---|
| Calling contract: meaning, units, ownership, mutation, collectiveness | the declaration (Doxygen on a public interface, a plain comment otherwise) |
| Reason for a local decision, a local trap, an order dependency | beside the decision, once |
| Derivation, algorithm taxonomy, the full index-space model, a tutorial | `src/<module>/doc/*.md`; the source keeps the local conclusion and names the file |
| Project-wide convention | `doc/` |
| Full bibliographic entry | `doc/literature_references.md`; the source cites author-year and the equation or section |
| The experiment that motivated a rule; superseded attempts; who reviewed it; who ruled | the devlog; the source keeps the rule |
| An old implementation | git |

**Cite, do not reproduce.** `// Messe et al. 2023, Eq. 13` is complete. The retired `paperN` aliases stay retired.

**A warning is stated once, at the place the mistake would be made.** If two independently edited sites can each make the same mistake, such as the matrix `send`, `receive`, and `broadcast` in `commtools.hpp`, state the warning at the first site. Add a one-line pointer at later sites instead of copying it. The MUMPS symmetry note exists in three copies today (the wrapper, `cl_IWG.hpp`, and the eigenvalue caller), and one of them drifted into the "lower triangle" error above. That is the cost of the copy.

**There is no line budget, but there is a length tripwire.** When a comment runs past about five lines, read it once more for a date, a measurement log, the *derivation* of a rejected alternative, a review credit, or a paper's argument. Those parts move. What is left — the rule, the one-line name of what was rejected, and the one fact the rule depends on — is usually two or three lines. Do not shorten by compressing into private abbreviations; shorten by removing what does not belong.

**Keep the local conclusion when moving background out.** A bare `// see doc/foo.md` beside a surprising branch is not enough; the reader should not need to open a second file to understand why the branch exists. State the condition and the consequence here, and point to the document for the derivation.

---

## 5. Doxygen documents the contract, not the name

`make doc` generates the API reference from `/** */` blocks, so those blocks are read twice: in the header and on the generated page. Both readers are served by the same rule.

**Write a `/** */` block on a public member when it states at least one thing the signature does not:** ownership transfer or non-transfer, units, a required layout, a precondition the caller must satisfy, a call-order dependency, MPI collectiveness, a non-obvious cost, or a status return that encodes an expected failure.

**Do not write one otherwise.** Delete `@brief Constructor`, `@brief Get the label`, `@param aMesh the mesh`, and `@return the number of rows`. A getter that returns `mFoo` gets nothing. A private helper used in one `.cpp` gets nothing unless its invariant is surprising. When a block is *both* — `@brief Destructor (nothing to release; all containers are owned by the mesh)` — the parenthetical stays and the restated name goes: `/** Nothing to release: every container is owned by the mesh. */`.

A class-level `@brief` — one sentence saying what the class is for — is welcome on any public class; it is the one place a Doxygen block may say *what* rather than *contract*. A header that today carries a multi-page tutorial in Doxygen markup (`cl_Material_UserDefined.hpp`, the top of `powerlaws.hpp`) keeps the class-level sentence and the API contract; the walkthrough moves to the module's `doc/` and the example deck.

```cpp
// Contract, not name: the one-argument Cell constructor reserves and does
// not create elements, which a caller reading Cell<T>( n ) cannot know.
/** Create an empty Cell with capacity for at least aReserve elements. */
Cell( const std::size_t aReserve );

// Nothing: the name is the documentation.
size_t size() const;
```

`@param` lines go only on the arguments that carry a unit, a layout, an ownership rule, or a precondition. `@param[in]` / `@param[out]` beats a prose paragraph. No `@author`, no `@date`, no change lists, no HTML.

The Fortran wrappers follow the same rule: `mumpstools.f90` documents its `bind(c)` entry points with `!>` blocks, `arpacktools.f90` with plain `!` headers. The argument meanings and the status encoding are the contract a C++ caller cannot see, and they belong there.

---

## 6. Banners

The `//------` line between definitions is house style. It is a position marker, and Martin's warning applies: a banner is only useful when banners are rare enough to be noticed. Three rules keep it that way in new and revised code:

- **One banner between top-level definitions** in a `.cpp`, and at most one per logical section (an access specifier, or a group of related accessors) in a header. Do not place a banner between every member, as several headers (`cl_Cell.hpp`, `cl_Mesh.hpp`) do today.
- **No banners inside a function body.** If a body needs visual separation, it needs a stage heading (§2), or it is two functions.
- **No sub-banners** (`// - - - - -` in C++, `! - - - - -` in Fortran). About 500 exist. Remove them when their containing file is touched, and write no new ones.

Whether to thin the existing banner population tree-wide is a formatter decision, taken separately from the commenting sweeps. Banners can be noisy but cannot be wrong, so thinning them can wait.

---

## 7. TODO comments

A TODO in source is a job that cannot be completed now but must be completed. It includes enough information for someone else to complete it. It carries an owner and an exit condition, and may point at the durable task file that tracks it:

```cpp
// TODO(cm): drop the dense fallback once the sparse Schur path is validated
//           on the 12-tape coil; see todo/<task>.md
```

(`todo/<task>.md` stands for the real task file. Write the actual file name, and only one that exists.)

This is narrower than Martin, who withdrew his first-edition tolerance of TODO comments: "TODO means Don't Do", so he no longer checks them in and instead does the thing, removes the need for it, or moves it to the backlog. The owned TODO above is that backlog entry left in place as a pointer, and it is the only form allowed.

Martin's second edition no longer checks any TODO in at all: "TODO means Don't Do" — do the thing, remove the reason for it, or put it in the backlog (2025, Ch. 5 "TODO Comments"). BELFEM keeps one narrow exception to that position, the owned form above, because its backlog is `todo/` and a one-line pointer from the code to the task file is cheaper for the next reader than a search. Everything else follows Martin: a bare `// todo: this should be obsolete soon`, `// todo: we might not need this anymore`, or the same `// todo: optimize` line pasted eight times is an apology, not a job. Each of the roughly 44 existing `TODO`/`FIXME` lines has one of three outcomes: it is rewritten to the form above, moved to a `todo/*.md` file with a one-line pointer left behind, or deleted as stale. The triage does not resolve the TODO itself.

**Exception — the closed cohomology core.** Eight of those lines sit in `fn_Smith.hpp` and `fn_Smith.cpp`, the two files of one of the six units closed to AI edits (`doc/ai_collaboration_protocol.md` §7.1). The protocol treats comment-only edits there as banned. Those TODOs are reported to Gregory in a devlog and are not touched by a sweep.

---

## 8. No history, no provenance, no dates

This is the rule most often broken in recent code, so it gets its own section.

**Git holds the *when* and the *who*. A comment says *what* and *why*, and nothing else.** No dates, no "changed on", no "removed 2026-08-29", no reviewer names, no audit round labels, no session IDs, no incident or debt-register numbers, no paths into `tmp/`, and no `Decision: <name>, <date>` lines. The maintainer of record takes the blame by default, and `git blame` answers the question when it is asked. Devlogs, todo files, and `doc/*.md` keep their dates because their templates require them; source comments do not.

The reason is not tidiness. A dated comment is a journal entry, and journal entries accumulate: the reader of `cl_FEM_Controller.cpp` today meets 28 of them and has to decide, for each, whether the date still matters. It never does. The rule the date was attached to either still holds, in which case the date is noise, or no longer holds, in which case the comment is a lie with a timestamp.

A review credit is worse than a date. `( Codex+Grok audit )` tells the reader that two models agreed, which is the weakest rung of the evidence ladder (`doc/ai_collaboration_protocol.md` §11) and says nothing about *why* the code is right. If the mechanism is worth recording, record the mechanism.

**Worked rewrite — a pool-size constant.** Before, from `mumpstools.f90`:

```fortran
! FIXED at 8. The setter that used to change it was removed 2026-08-29:
! it had no caller anywhere, it did not reallocate the pool, and a call
! after the first create would have left gMaxNumSolvers and the actual
! array size disagreeing -- a scan past the end, or a deallocate with
! high slots still occupied. Raising the pool size is a reallocation
! feature; add it as one if it is ever wanted, not as a bare setter
integer( int_t ), save :: gMaxNumSolvers = 8
```

Six of the seven lines are the history of a setter that no longer exists. After:

```fortran
! The pool size is fixed. Growing it requires reallocating gSolvers.
integer( int_t ), save :: gMaxNumSolvers = 8
```

**Worked rewrite — a guard with its experiment attached.** Before, from the controller watchdog:

```cpp
// SPARE a stalled-best iterate while the line search is regrowing
// its relaxation: a no-new-best window that merely spans an
// overshoot/backtrack/recovery cycle is not a stall. Measured
// ( coarse tapestack3d, 2026-08-21 ): the AIMD recovery from a
// backtracked omega ~0.1 takes ~17-22 iterates at beta = 1.1,
// window W = 8 cut it mid-climb and CASCADED ( halving improved
// the best residual and was cut again, 2.10 -> 1.05 -> 0.53 ms ).
// A genuine stall has omega pinned or shrinking ( the hd floor
// grinds, and 383@2.10 with omega frozen at 0.053 ), so it still
// fires. The motivation parallels Chamberlain-Powell-Lemarechal's
// watchdog technique -- do not punish a nonmonotone step that is
// still contracting -- though theirs relaxes a line search and
// this spares a timestep cut. NaN prev ( cold start ) must not
// spare.
if ( ! std::isnan( tOmegaPrev ) && aOmega > tOmegaPrev )
```

The rule is four sentences; the rest is the run that found it and belongs in the devlog that recorded it. After:

```cpp
// A rising relaxation factor means the line search is recovering from a
// backtrack, not stalling, and that recovery can outlast the watchdog
// window: do not cut the timestep while omega is still climbing. A
// genuine stall has omega pinned or shrinking and still fires. A NaN
// previous omega (cold start) is not a recovery. Unlike the classical
// watchdog, which relaxes a line search, this spares a timestep cut.
if ( ! std::isnan( tOmegaPrev ) && aOmega > tOmegaPrev )
```

What went: the deck, the date, the iterate counts, the trace, and the cascade. What stayed: the mechanism, the reason the window alone cannot decide, the cold-start exception, and the one sentence that keeps a reader from mapping this guard onto the textbook watchdog. The numbers went *because* they were conditional: the recovery length depends on the growth factor and the residual, and the window itself is rescaled during retries, so "eight" and "twenty" would have been stale against the tree on the day they were written.

---

## 9. Reference implementations and disabled code

Commented-out code is deleted. Git remembers it, and the compiler does not check a commented block. It is a second specification that is guaranteed to rot.

Oliveira & Stewart (Ch. 12) recommend keeping the "simple, intelligible (unoptimized)" version of a tuned kernel in a comment, because architectures and compilers change. BELFEM pursues the same goal differently: **a reference implementation that is worth keeping is worth executing.** It lives in a test under `tests/`, compared against the tuned kernel with a stated tolerance, where `make check` keeps it honest. A small, labeled formula or pseudocode fragment above the tuned loop is fine as *explanation*; an old implementation in comment form is not.

Delete disabled diagnostics (`// print_banner();`, `// std::cout << …`). A disabled implementation next to its replacement (the old `save_fields` in `cl_Mesh.cpp`) is deleted once its replacement has a passing gate.

The ownership lines that sometimes sit above disabled code — `// facets are deleted by sideset` — are contract, not history. They stay; the loop under them goes.

---

## 10. Language

Comments are read by people who did not sit in the session that wrote them, and often by non-native speakers. The 2026-09-15 agreements for `todo/` prose also apply unchanged to source comments:

1. **Name the function or file, never a coined label.** Write "the `clean()` call in the constructor", not "the constructor clean". A term invented during a session does not appear without a definition.
2. **Define a term on first use or do not use it.** Expand an abbreviation the first time.
3. **One idea per sentence.** Split a sentence that carries a decision, a reason, and a consequence.
4. **No pipeline arrows** (`a -> b -> c`). One sentence per handoff.
5. **Say what kind of thing a term is** — a mesh property, an implementation failure, a user-visible behavior.
6. **Full sentences for a reason or a contract; fragments are fine for a citation, a unit, or a symbol.**
7. **US spelling**, no jokes, no venting, no mumbling.

A comment that needs its own explanation has failed (Martin's "unobvious connection"). If the terms in an expression need naming before the comment makes sense, name them in the code — a `const` with a unit — and the comment shrinks or disappears (Rouson Rule 1.2).

---

## 11. Extract a function, or write a comment?

Martin's first remedy for an explanatory comment is a function whose name is the comment. `doc/coding_philosophy.md` already rejects that as the default inside numerical kernels. Splitting a stencil across a dozen named methods destroys the correspondence between the code and the published equation that Mathematical Readability protects. Rouson makes the same exception for domain notation (`u`, `v`, `i`, `j` inside a Navier–Stokes kernel).

The resolution: **extract at the physics-step boundary, not at the line.** A forty-line `if` nest that needs the heading "handle the superconducting-to-normal transition" may become `is_quenching()`. The transition ODE inside it stays as one block with one citation. In `src/math` and `src/physics`, a one-line equation anchor *is* the extraction, and a symbol annotation beside a statement is how the block stays readable without being split.

The second edition of *Clean Code* prints the other side of this argument: in its Appendix, Ousterhout holds that interfaces and abstractions cannot be defined without many comments and that missing comments cost more than bad ones, while Martin holds that comments as generally practiced are a net loss. The two agree on one sentence, and it is the one this document adopts for implementation code: "implementation code only needs comments when the code is nonobvious" (Martin 2025, Appendix "Comments Summary"). For public interfaces this document sits nearer Ousterhout: the contract is written down (§5), because the caller does not read the body.

---

## 12. Comments are edited with the code

A change that alters a loop's direction, a unit, a default, an ownership rule, or a call order edits or deletes the comment in the same diff.

A reviewer who finds a stale *narration* or a stale *stage heading* deletes it rather than asking for a fix. A stale *contract* — ownership, units, collectiveness, call order — is investigated first, because the disagreement may be the code's fault: the comment may be recording the rule the implementation broke. The same holds for a scientific claim. The code establishes what executes, not what should execute; rewriting the comment to bless the implementation defeats the purpose of having written it.

If you are not willing to update a comment in the commit that changes the code, do not write the comment. Write a check instead, chosen by the error-tier policy in `doc/coding_philosophy.md` (`BELFEM_ASSERT` is compiled out in release; `BELFEM_ERROR` is not), or a test that actually exercises the condition. A comment is kept honest by nobody; a check that runs is kept honest by the run. Thomas & Hunt say the same of a contract: writing it as a comment "is a great start", and having the program check it is the greater benefit (Topic 23 "Assertions", Topic 25).

---

## 13. Review checklist

For each comment in a diff:

- [ ] Does it say something the names, types, prefixes, or checks do not?
- [ ] Is it a contract, a reason, a citation, a warning, an amplification, a stage heading, a symbol annotation, or an owned TODO?
- [ ] Is it beside the code it describes, and only about that code?
- [ ] If it is narration or a stage heading: would it still be true if the next line changed in a plausible way? (A citation, a contract, or a warning is judged by §2, not by this question.)
- [ ] Does it carry a date, a review credit, a round label, an incident or debt-register ID, a ruling with a name on it, or a `tmp/` path? Remove it.
- [ ] Does it carry a measurement log, the derivation of a rejected alternative, or a paper's argument? Move it; keep the rule and the one-line name of what was rejected.
- [ ] If it is Doxygen, is this a public interface, and does the block state a contract?
- [ ] If it is a banner, is it between definitions and not inside a body?
- [ ] If I delete it, does the file get harder to use, or only shorter?

If the last answer is "only shorter", delete it.

For a cleanup sweep, five further rules:

- **No change to any executable path.** Comment lines only. Line numbers embedded by `BELFEM_ERROR` and `BELFEM_ASSERT` through `__LINE__` will shift, so a byte-identical binary is not the gate; an unchanged `make check` and an unchanged diff of the preprocessed non-comment tokens are.
- **Never delete a comment you do not understand.** A line that matches a noise pattern but reads as a warning or a contract is rewritten, not removed.
- **Correct the misleading ones before shortening anything**, in their own commit, so a reviewer can judge accuracy separately from volume.
- **Exclude what is not a comment**: compiler directives (`!$omp`, `!dir$`, `#pragma`) and string literals, even when they contain a date or a name.
- **Exclude the closed cohomology core** (`cl_Cohomology`, `cl_Homology`, `cl_SimplicialComplex`, `cl_Chain`, `cl_Cochain`, `fn_Smith`). Findings there go to Gregory in a devlog.

---

## 14. Relation to the other conventions

- `doc/coding_philosophy.md`, "Mathematical Readability", the error-tier policy, and the note on *Clean Code*: the citation form, the check selection, and the kernel-decomposition exception this document builds on.
- `CLAUDE.md`, "Citation Format": author-year plus equation or section; `paperN` aliases retired. `CLAUDE.md`, "Documents state facts inline": a `doc/` file must not cite a devlog, a todo file, or a task ID; this document extends that rule to source comments.
- `doc/documentation_guidelines.md`: where derivations, guides, and session records live.
- `doc/ai_collaboration_protocol.md` §7.1 (the closed core) and §11 (why a review credit in a comment is not evidence).
- The 2026-09-15 language agreements (`todo/` prose for non-native readers; no dates in source comments) are folded into §8 and §10.

None of this is enforced mechanically yet. A comment census script (banners, narration patterns, dated lines, provenance markers, dead code, `@brief` restatements, with the directive, string-literal and closed-core exclusions above) is the enforcement target. It will be specified in the `todo/` plan for the first cleanup sweep, which does not exist yet.

---

## Sources

- Robert C. Martin, *Clean Code: A Handbook of Agile Software Craftsmanship*, 2nd edition, Addison-Wesley / Pearson, 2025, Ch. 5 "Comments" (pp. 89–116) — the taxonomy of good comments (legal, informative, intent, clarification, warning of consequences, amplification, public-API docs) and bad ones (mumbling, redundant, misleading, mandated, journal, noise, scary noise, TODO, position markers, attributions and bylines, commented-out code, HTML, nonlocal information, too much information, unobvious connection, function headers, non-public API docs). Two changes from the 2008 first edition matter here: TODO comments moved from the good list to the bad list (§7 records the exception), and the closing-brace category was dropped (the row in §3 stands as a house rule, not a citation).
- Suely Oliveira and David E. Stewart, *Writing Scientific Software: A Guide to Good Style*, Cambridge University Press, 2006, §7.9 "Comments" (interface contracts, invariants, warnings, *what* not *how*, once and only once), §7.10 "Documentation" (keep it small; document interfaces), Ch. 12 introduction (reference code beside tuned code).
- Damian Rouson, Jim Xia, and Xiaofeng Xu, *Scientific Software Design: The Object-Oriented Way*, Cambridge University Press, 2011, §1.7 "Elements of Style", Rules 1.1–1.3 (code that comments itself, name all constants, make constants constant; the domain-notation exception), §7.2.1 (comment intent more often as the problem grows).
- David Thomas and Andrew Hunt, *The Pragmatic Programmer*, 20th anniversary ed., Addison-Wesley, 2020 — Topic 7, Tip 13 "Build Documentation In, Don't Bolt It On" (API comments yes, mechanical per-function comments no, non-API comments discuss *why*); Topic 9 "Duplication in Documentation" (the comment and the code drift apart); Topic 23 "Design by Contract" (the contract as a comment, then as an assertion); Topic 25 (assertions check the impossible). Cited by Topic and Tip, never by page.

All four books are in `literature/books/coding/` when that directory is present (`oliveira2006`, `rouson2011`, `martin2025`, `thomas2020`, each with a navigation `.md`, routed from its `index.md`); the citations above are the published works, not the local extracts. Full bibliographic entries: `doc/literature_references.md`, "Software Craft".