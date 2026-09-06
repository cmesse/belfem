# L-21 Repair: the Wrapper-Policy Card Cites a Contract That Moved, and Cannot Enforce Itself

**Date:** 2026-08-31
**Purpose:** Repair lesson card `L-21` ("Third-party libraries are accessed only through their
dedicated wrappers") in `doc/lessons_learned.md:491-528`. Its central citation points at a stub,
the wrapper file it names is the one file that does not carry the contract it describes, and its
proposed enforcement grep is blind to the largest body of direct MPI left in the tree while
flagging type declarations and comments as violations. Mechanism: correct the three factual
claims, then replace the enforcement recipe with a call-aware sweep — once Christian rules on
whether the Fortran solver shims are an intended exception.
**Module:** `doc/` (+ `doc/coding_philosophy.md`; no `src/` change proposed)
**AIs involved:** Claude (verification + plan), Codex `gpt-5.6-sol`/xhigh, Grok `grok-4.6`/xhigh
(both raised the findings independently in the 2026-08-31 round-2 jury)
**Status:** DONE 2026-08-31 — R1–R7 applied, O1 ruled by Christian (Fortran drivers ARE the
wrapper the C++ tree uses), O2 settled as a sibling script. The R5 sweep's first run was red on
11 sites, none of them MPI, and Christian ruled all 11 by design the same day (PETSc's scope
reaches into `comm/`; the materials HDF5 database probes are sanctioned). Both encoded as owner
entries that print their reason; sweep green, falsifier re-run after the widening. **No residue.**
See `devlog/dl20260831_l21_card_repair.md`.

> **Scope guards:**
> - **Documentation only.** No `src/` edit is proposed. In particular, `mumpstools.f90` and
>   `parpacktools.f90` are **NOT** to be "fixed" by this plan — whether they are violations or an
>   intended exception is O1, and rewriting them would break the Fortran/C++ boundary the tree
>   deliberately uses.
> - The **rule itself is not in question.** Christian stated it as hard policy on 2026-08-30 and
>   the four C++ sites it closed stay closed. This plan repairs the card's *claims* and its
>   *enforcement*, not its content.
> - Building the mechanical check is in scope as a **design** (R5); wiring it into CI is not.

---

## 1. Current Behaviour and How It Fails

All five defects were raised by the round-2 auditors and **re-verified firsthand by Claude on
2026-08-31** against the named files. Confidence high on all except D5 (a judgement about the
file's own promise).

| Failure | Mechanism | Evidence |
|---|---|---|
| **D1** The card's key citation names a stub | L-21 cites `src/core/assert.cpp:294-320` for "the abort path's `MPI_Initialized`/`MPI_Finalized` guards and its deliberate use of `MPI_COMM_WORLD` over `gComm.world()`". None of that is there. `error_abort()` at `:294-305` is a five-line body whose own comment says the contract "moved verbatim into `comm_abort`… on 2026-08-30"; `:309-328` is `gThrowOnError` / `throw_on_error()` / `set_throw_on_error`, no MPI at all | `doc/lessons_learned.md:504-506`; `src/core/assert.cpp:294-305`, `:309-328`; live contract at `src/comm/cl_Communicator.hpp:234-257`, definition `src/comm/commtools.cpp:67` |
| **D2** `coding_philosophy.md` is stale the same way | It still says `error_abort()` "under MPI aborts the job with `MPI_Abort(MPI_COMM_WORLD, 1)` — guarded by `MPI_Initialized`/`MPI_Finalized`", as though the code were still in `assert.cpp`. It is the sentence L-21 was written from | `doc/coding_philosophy.md:631` |
| **D2b** …and repeats itself | The same paragraph carries the sentence "**It is a test hook, not a configuration knob**: a production run must keep the abort, because a throw that escapes `main` terminates one rank and leaves its peers blocked in a collective, and nothing in `src/` catches BELFEM errors." **twice**, near-verbatim, the second time with "the `MPI_Abort` reaction" | `doc/coding_philosophy.md:631` (same line) |
| **D3** The named wrapper file is the wrong one for the example | The Rule says "`commtools.hpp` for MPI". But `comm_abort` — the card's own worked example — is declared in `cl_Communicator.hpp` **on purpose**: `commtools.hpp` pulls `cl_Vector`/`cl_Matrix` from linalg, which is not on core's include path, and `assert.cpp` must reach it. So the one file the card names is the one that does not carry the contract it holds up | `doc/lessons_learned.md:499`; `src/comm/cl_Communicator.hpp:225-231` (the placement rationale, in-tree) |
| **D4** The enforcement recipe cannot enforce the rule | `grep -rn 'MPI_' src/ --include='*.cpp' --include='*.hpp' \| grep -v src/comm/` is **blind to Fortran**: `src/sparse/mumpstools.f90` (8 `MPI_` lines) and `src/sparse/parpacktools.f90` (25) are outside `src/comm/` and are violations *by the letter of the card*. It is simultaneously **noisy on non-calls**: `typedef int MPI_Comm ;` (a serial-build shim), `MPI_Comm` parameters and members, and comments | `doc/lessons_learned.md:525-527`; `src/sparse/mumpstools.f90`, `src/sparse/parpacktools.f90`; `src/sparse/petsctools.hpp:40,68`, `src/sparse/st_SolverPetscData.hpp:30` |
| **D5** The card cites no incident | `doc/lessons_learned.md:10` promises "Every rule cites the incidents that paid for it". L-21's Evidence paragraph narrates the 2026-08-30 sweep in prose and carries no `INC-NNN` | `doc/lessons_learned.md:509-517` vs `:10` |

**Bottom line:** the card that exists to stop stale-citation and wrong-contract errors contains
one of each — D1 is exactly the `INC-543` shape (a genuine line number attached to a statement
that is no longer true), inside the card meant to prevent it — and its enforcement clause would
pass a tree that still contains 33 direct MPI lines while failing one that contains none.

## 2. Why Repair Rather Than Rewrite

L-21's *rule* has paid for itself: four sites closed the day it was stated, one of them written
that same day by a session that knew the tree well. The failure is in the supporting apparatus,
which is repairable in place. Rewriting the card would lose the "the violations were not
ignorance, they were the default behaviour under time pressure" observation, which is the reason
it is a tripwire rather than a style note.

The alternative considered and rejected: **demote L-21 to prose in `coding_philosophy.md`.**
Rejected because the trigger ("about to write a vendor-prefixed call") is exactly the moment a
Layer 1 tripwire is for, and the philosophy document is not loaded per-session.

## 3. Gap Table

| # | Claim in the card | Needed for | True today? | Class | Citation / rationale |
|---|---|---|---|---|---|
| 1 | abort contract lives in `assert.cpp:294-320` | the "why the rule earns its keep" argument | **no** — moved 2026-08-30 | (c) explicit fix | `assert.cpp:294-305` is a stub; `cl_Communicator.hpp:234-257` is the contract |
| 2 | `commtools.hpp` is *the* MPI wrapper | telling a reader where to add a missing operation | **partly** — true for collectives, false for `comm_abort` | (c) explicit fix | `cl_Communicator.hpp:225-231` |
| 3 | vendor calls appear only in the owning module | the rule's universal form | **no** for Fortran | (b) open — see **O1** | `mumpstools.f90`, `parpacktools.f90` |
| 4 | the grep is a usable check | `Status:` line's enforcement claim | **no** | (c) explicit fix, after O1 | D4 evidence |
| 5 | every rule cites its incidents | the file's own contract | **no** for L-21 | (a) rebuildable — the sweep is recorded in `dl20260830_*` | `:509-517` |
| 6 | `coding_philosophy.md` describes the abort path | the source L-21 was written from | **no** + duplicated sentence | (c) explicit fix | `:631` |

### 3.1 Cross-cutting finding

**Both D1 and D2 descend from one edit.** The 2026-08-30 move of the contract from `assert.cpp`
into `comm_abort` updated the code and its in-place comments meticulously — `assert.cpp:296-303`
and `cl_Communicator.hpp:225-231` both explain the move — and updated neither of the two
documents that describe it. The code is self-documenting about its own history; the prose layer
did not follow. Any fix that repairs only the card leaves `coding_philosophy.md` as a live source
for re-introducing the same claim.

## 4. Ordered Steps

- [x] **R1** — Repoint the D1 citation. In `doc/lessons_learned.md:504-506`, replace
      `src/core/assert.cpp:294-320` with `src/comm/cl_Communicator.hpp:234-257` (contract) and
      `src/comm/commtools.cpp:67` (definition). Keep the sentence's claim — the guards and the
      `MPI_COMM_WORLD` choice are real, they just live elsewhere now.
- [x] **R2** — Fix `doc/coding_philosophy.md:631`: state that `error_abort()` delegates to
      `comm_abort`, and **delete the duplicated sentence** (D2b). Do not restate the contract a
      third time — point at the header.
- [x] **R3** — Fix the Rule's wrapper naming (`:499`): "the **comm module** for MPI
      (`commtools.hpp` for collectives; `cl_Communicator.hpp` declares `comm_abort`, which core
      must reach without pulling linalg)". One clause, because the exception is load-bearing.
- [x] **R4** — **Resolve O1 with Christian.** Blocks R5. Nothing else in this plan waits on it.
      **Ruled 2026-08-31, option (a):** "MUMPS and other Fortran tools are excluded from the MPI
      wrapper ban. The wrapper requirement is specific to C++. Fortran drivers for third packages
      such as MUMPS, PARDISO, ARPACK, ARE the wrappers that the C++ codebase uses. The
      communication calls stay within the confinements of the driver, such as MUMPS internal MPI
      communications are confined within MUMPS."
- [x] **R5** *(after R4)* — Replace the enforcement recipe. Requirements the new sweep must meet:
      include `*.f90`; match **calls**, not `MPI_Comm` type uses, parameters, members, `typedef`s
      or comments; carry an explicit per-module allowlist rather than a single `grep -v src/comm/`;
      and print what it skipped. Design only — wiring into CI is out of scope.
- [x] **R6** — Add `INC-NNN` citations to L-21's Evidence paragraph (D5). The 2026-08-30 sweep is
      already in the catalogue; find its rows rather than filing new ones.
      **Correction on execution:** it was *not* in the catalogue. The addendum ran INC-540…564 and
      held no row for the wrapper sweep, only the devlog `dl20260830_mpi_wrapper_contract.md`. Two
      rows filed instead: INC-565 (the four-site sweep) and INC-566 (the wrapper drafted to close
      it re-introducing the `gComm.world()` contract hours later, plus four further defects the
      blind jury caught). Header range bumped to INC-566.
- [x] **R7** — **Gate.** Run the R5 sweep and confirm three things by inspection of its output:
      (a) it reports the `mumpstools.f90` / `parpacktools.f90` sites (or excludes them by a named
      allowlist entry, per O1); (b) it does **not** report `petsctools.hpp:40`'s `typedef`,
      `st_SolverPetscData.hpp:30`'s member, or the `assert.hpp` comments; (c) it reports **zero**
      C++ violations outside the comm module, matching the 2026-08-30 closure. Falsifier: deliberately
      open-code an `MPI_Barrier` in a scratch `.cpp` under `src/sparse/` and confirm the sweep goes
      red, then remove it. A check that cannot go red is not a check (`L-01`).

## 5. Open Design Questions

- **O1 — Are the Fortran solver shims an intended exception to L-21?** `mumpstools.f90` and
  `parpacktools.f90` call `MPI_BARRIER`, `MPI_COMM_SIZE`, `MPI_ALLGATHER(V)`, `MPI_ALLREDUCE` and
  use `MPI_COMM_WORLD` directly. By the letter of the card they are violations; in practice they
  are Fortran kernels talking to MUMPS and PARPACK, both collective on `MPI_COMM_WORLD`, on the
  far side of a language boundary the C++ wrapper layer does not cross. **Nothing in the tree
  states this either way** — that is the actual defect, independent of which answer is right.
  Options: **(a)** write the exception into the card ("the C++ wrapper layer; Fortran kernels
  interface their own libraries directly") — *recommended*, it matches the tree and costs nothing;
  **(b)** require a Fortran-side wrapper module — real work, and the contracts differ
  (`MPI_INTEGER`, not `comm_type<T>`); **(c)** case-by-case, which is what exists now and is why
  the question is open. **Needs Christian — this is a policy scope decision, not a vote.**
- **O2 — Where does the R5 check live?** `scripts/check_doc_claims.py` is scoped to `CLAUDE.md`
  and `doc/coding_philosophy.md` and checks *documentation* claims; this is a *source* sweep.
  Options: extend it (cheap, wrong shape), a sibling script (`scripts/check_wrapper_policy.py`),
  or a `make` target. Recommend a sibling script; decide at R5.

## 7. Definition of Done

- [x] R1, R2, R3, R6 applied; every citation in the touched paragraphs re-opened at its file
      before the edit lands (this plan's own findings came from doing exactly that).
- [x] O1 answered by Christian and written into the card in one clause.
- [x] R5 sweep designed, R7 gate run with the deliberate-break falsifier, output pasted into the
      devlog.
- [x] L-21's `Status:` line updated to state what is now mechanically checked and what is not —
      **and the file header's "none is mechanically enforced yet" (`doc/lessons_learned.md:14`)
      re-examined, since it is already false for L-08** (`check_doc_claims.py:292-347`). That
      header claim is out of scope to fix here but must not be left contradicting a card this
      plan updates.
- [x] Devlog entry; `todo/README.md` entry updated to DONE; plan moved to `todo/closed/`.

## 8. Audit Trail

- `tmp/ai_exchange/review_cohomology_edit_ban.md` — round 2 (2026-08-31), Codex `gpt-5.6-sol`
  xhigh + Grok `grok-4.6` xhigh. D1/D4 raised by both independently; D3 and the `.f90` blindness
  by Grok; D2/D5 by Codex. Claude's verification pass re-opened every citation; D2b was found in
  that pass and is in neither auditor's entry.
- Ephemeral. Distil before the exchange is swept.
