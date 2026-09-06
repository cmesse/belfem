# L-21 Repaired: the Card Now Cites What Is True and Runs a Check That Can Go Red

**Date:** 2026-08-31
**Purpose:** Record the repair of lesson card `L-21` (wrapper policy) — three corrected claims,
Christian's ruling on the Fortran drivers, the call-aware sweep that replaces the card's
unusable grep, and the eleven live violations that sweep found on its first run
**Module:** `doc/`, `scripts/`

## What was wrong

Raised by the 2026-08-31 round-2 jury (Codex `gpt-5.6-sol`/xhigh, Grok `grok-4.6`/xhigh),
planned in `todo/closed/l21_wrapper_card_repair.md`, every citation re-opened firsthand before
the edit landed.

- **D1** The card's central citation, `src/core/assert.cpp:294-320`, pointed at a five-line stub.
  The MPI lifecycle contract it described moved into `comm_abort` on 2026-08-30 —
  `src/comm/cl_Communicator.hpp:235-258` (contract) and `src/comm/commtools.cpp:67` (definition).
  `INC-543`'s exact shape, inside the card written to prevent it.
- **D2/D2b** `doc/coding_philosophy.md:631` was stale the same way, and repeated one sentence
  ("It is a test hook, not a configuration knob…") near-verbatim.
- **D3** The Rule named `commtools.hpp` as *the* MPI wrapper, but `comm_abort` — its own worked
  example — is declared in `cl_Communicator.hpp` deliberately: `commtools.hpp` drags
  `cl_Vector`/`cl_Matrix` from linalg, which core's include path does not carry.
- **D4** The enforcement grep could not enforce the rule. On the tree as it stood it returned
  **fourteen** `MPI_` hits outside `src/comm/`, **none** of them a call: comments, a serial-build
  `typedef int MPI_Comm`, parameters and members, a Blaze macro, a PETSc error enumerator. It was
  simultaneously blind to the 33 `MPI_` lines in `src/sparse/*.f90`.
- **D5** The card carried no `INC-NNN`, against the file's own promise at `:10`.

**Both D1 and D2 descend from one edit.** The 2026-08-30 move updated the code and its in-place
comments meticulously — `assert.cpp:296-303` and `cl_Communicator.hpp:225-231` both explain it —
and updated neither document that describes it. The code was self-documenting about its own
history; the prose layer did not follow.

## Christian's ruling on the Fortran drivers (O1)

> MUMPS and other Fortran tools are excluded from the MPI wrapper ban. The wrapper requirement is
> specific to C++. Fortran drivers for third-party packages such as MUMPS, PARDISO, ARPACK **are**
> the wrappers that the C++ codebase uses. The communication calls stay within the confines of the
> driver — MUMPS's internal MPI communication is confined within MUMPS.

Written into the card's **Exceptions** in one clause, and into the sweep as a named, printed
exclusion rather than a silent `grep -v`. Nothing in `src/` was changed; the plan explicitly
forbade "fixing" the Fortran, and the ruling confirms there was nothing to fix.

## The sweep (`scripts/check_wrapper_policy.py`)

Matches **calls** — a vendor-prefixed identifier followed by an open parenthesis, after comments
and string literals are stripped — across MPI, HDF5, PETSc, MUMPS, STRUMPACK, SuperLU,
LAPACK/BLAS and ARPACK, each with its owning module. Every exclusion prints by name with its
reason. In-wrapper calls are counted and shown, not hidden, so an empty owner column is visible
as a broken pattern rather than read as a clean tree.

**Falsifier run** (R7 gate, scratch root, then removed): a deliberate `MPI_Barrier` under
`sparse/` was reported; on the same file a commented `MPI_Barrier`, a `typedef int MPI_Comm`, an
`MPI_Comm` parameter and a `"MPI_Barrier("` string literal were **not**; a call inside `comm/`
counted as in-wrapper. Exit 1 on the violation, and a second probe confirmed two calls on one
line report as two work items. A check that cannot go red is not a check (`L-01`).

**The in-wrapper column immediately earned its keep, against its own author.** The first
STRUMPACK pattern matched a `STRUMPACK_` C prefix; STRUMPACK is a C++ namespace API
(`strumpack::StrumpackSparseSolver< real, int >( … )`), so the pattern matched nothing anywhere
and would have read as "no STRUMPACK violations" forever. The zero was visible only because
in-wrapper calls are printed. Fixed, and generalised: a family whose pattern matches nothing
anywhere now prints a WARNING, unless it carries a stated reason to be empty — MUMPS and ARPACK
do, since both are driven from `*.f90` and have no direct C++ call by design, which is the O1
ruling showing up as a testable prediction and holding.

**Out of reach, stated rather than hidden:** vendor *constants* (`MPI_SUM`, `H5T_NATIVE_DOUBLE`)
crossing a module boundary leak the same token, but are not syntactically distinguishable from a
mention. That stays review work.

## What the first real run found — 11 sites, none of them MPI, and all of them by design

The MPI closure of 2026-08-30 holds: zero direct MPI calls outside `src/comm/`. The families the
old grep never swept produced eleven hits, and **Christian ruled all eleven by design the same
day** — the sweep's first output was not a defect list but a map of two undocumented ownership
facts.

| site | call | note |
|---|---|---|
| `src/comm/cl_Communicator.cpp:138` | `PetscOptionsSetValue` | PETSc initialization sits beside `MPI_Init` in the communicator. Arguably an ownership-map question rather than a violation — needs a ruling |
| `src/physics/materials/cl_JcFunction_Database.hpp` (7 sites) | `H5Fopen`, `H5Lexists`, `H5Gopen2` ×2, `H5Gclose` ×2, `H5Fclose` | carries a source comment: "raw handle rather than the HDF5 wrapper: `/source` is a SIBLING of the group Database opened, and the wrapper exposes no root handle to probe" |
| `src/physics/materials/fn_rho_database_is_current.hpp` (3 sites) | `H5Fopen`, `H5Lexists`, `H5Fclose` | the idiom the above cites as precedent |

**Ruled 2026-08-31, both by design.** PETSc's requirements are unusual enough that its scope
reaches into the communicator: `PetscInitialize` is paired with `MPI_Init` in
`Communicator::init` and cannot live behind a `Solver*` class. And the materials HDF5 probes are
sanctioned as written — they were never the "wrapper lacks the operation" work item they looked
like from outside.

Both are now **owner entries that print their reason** on every run, not silent allowlist lines:
`comm/` owns PETSc by ruling, `physics/materials/` owns HDF5 by ruling. The scoping is tight and
was proven so rather than asserted — a re-run falsifier confirms an `MPI_` call in
`physics/materials/` is still a finding and an `H5*` call in `physics/gasmodels/` is still a
finding. Widening an allowlist is exactly when a check needs re-proving it can go red, so it was
re-proved after the widening, not before.

**The sweep is green on the tree**, and the green is now meaningful: two ownership facts that
existed only in people's heads are written down and mechanically applied.

## Landed

- `doc/lessons_learned.md` — L-21: citation repointed (R1), wrapper naming corrected with the
  layering reason (R3), Enforcement clause replaced (R5), Exceptions rewritten to carry all three
  rulings (Fortran drivers, PETSc into comm, materials HDF5), `Status:` green, Evidence citing
  INC-565/566 (R6). File header `Status:` corrected — it claimed nothing was mechanically
  enforced, which was already false for L-08.
- `doc/lessons_learned_evidence.md` — new addendum block, INC-565 (the four-site sweep) and
  INC-566 (the wrapper drafted to close it re-introducing the contract it was wrapping, plus four
  further defects found by the blind jury). The plan assumed these rows existed; they did not.
- `doc/coding_philosophy.md:631` — abort path now delegates to `comm_abort` with the contract
  cited at the header; duplicated sentence removed (R2).
- `scripts/check_wrapper_policy.py` — new (R5), falsifier-gated twice (R7): once on first
  write, and again after the two by-design owners were added, since a widened allowlist is
  precisely when a green run stops being evidence.
- `CLAUDE.md` — sweep registered in the tooling block. `check_doc_claims.py`: 37/37 after.

## Owed

Nothing from this plan. The two questions it opened were ruled the same day and are encoded.

Worth naming for whoever reads the sweep next: it is green because three exceptions are declared,
and each one is a place where the tree's real structure differs from the rule's clean statement.
That is the honest shape of an enforced policy, not a weakening of it — but a fourth exception
added without a ruling would be indistinguishable from these three, so the rule that keeps this
check trustworthy is that **exceptions are ruled by Christian and print their reason**.
