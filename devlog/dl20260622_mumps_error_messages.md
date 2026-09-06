# MUMPS error/warning message overhaul (5.9.0)

**Date:** 2026-06-22
**Purpose:** Bring `MUMPS::error_message` up to MUMPS 5.9.0, fix latent bugs, add a
proper warning decoder, and make `free()` diagnostics use real INFO(2) details.
**Module:** sparse (`src/sparse/cl_SolverMUMPS.{cpp,hpp}`, `mumpstools.{f90,hpp}`)

## Context

`MUMPS::error_message` only covered INFO(1) codes `-1..-53` (transcribed from an
older manual) and had several real bugs. Task: extend/repair it against the
5.9.0 user guide (`tmp/userguide_5.9.0.txt`), keeping the C++ interface simple —
BELFEM deliberately exposes only a few ICNTLs (the Fortran glue sets just ICNTL
4,7,10,11,14,29,33,35; `mumpstools.f90:269-324`), so messages hinging on
unexposed ICNTLs were simplified to drop unactionable `ICNTL(xx)` jargon.

Two-AI workflow: Claude drafted, **Codex** audited twice (message correctness,
then the Fortran/C interop change) — both signed off high confidence.

## What changed

**1. Bug/factual fixes in the old `-1..-53` table** (wording otherwise preserved):
- `-2`: was "matrix size out of range" testing `INFO(4)`; it is purely NNZ
  out-of-range with detail in `INFO(2)` (N is the separate `-16`).
- `-5`, `-7`: now apply the "negative INFO(2) ⇒ ×1e6" size rule (was printing the
  raw negative value).
- `-9`: read `INFO(3)`; the missing-entry count is in `INFO(2)`. Reworded to
  "Main real/complex workarray S is too small ( missing N entries )".
- `-11`: the non-positive branch left a literal `%i` unsubstituted; identifier
  `LWK_USER`→`LWKUSER`.
- `-22`: corrected the pointer-array label table (`SHUR`→`SCHUR`,
  `LISTVAR_SHUR`→`LISTVAR_SCHUR`, stray `R` before `PERM_IN`,
  `IRN_loc/JCN_loc`→`IRN_loc, JCN_loc or A_loc`) and added entries 17 (`IRHS_loc`)
  / 18 (`RHS_loc`).
- Factual: `-27` `NZ_RHS+1`→`NRHS+1`; `-28` `IRHS_PTR`→`IRHS_PTR(1)`; `-30`
  `SHUR_LLD`→`SCHUR_LLD`; `-32` `NHRS`→`NRHS`; `-36` `ICNTL(26)`→`ICNTL(25)`;
  `-41` `LWK USER`→`LWKUSER`; `-44` detail is `ICNTL(31)` not `(32)`; `-45`
  `RHS`→`NRHS`; `-47`/`-48` reworded for accuracy.
- Prose typos swept (`unsufficent`, `hadle`, `Incompativle`, `resprected`,
  "but but", "values if", trailing dash, unicode `−`→`-`).

**2. New codes `-54..-90` + `-800`** added (BLR, distributed-RHS, config-file,
save/restore, SCOTCH, out-of-core, …), simplified per the ICNTL policy.
Deliberately **omitted** `-92, -100, -102, -104, -105, -106, -107`
(GPU/cuBLASXt/XKBlas/OpenMP-thread/rank-revealing) — features BELFEM never
enables; they route through an improved `default:` that prints the actual
`INFO(1)`.

**3. New `MUMPS::warning_message(const int_t*, Cell<string>&)`** — positive INFO(1)
is a **bitwise sum** of warnings (`+1/+2/+4/+8/+16`), so it decodes each set bit
into its own line. Only `+1` and `+16` carry an INFO(2) detail, and INFO(2)
belongs to whichever was raised *last*; when both are set the count is suppressed
to avoid mis-attribution.

**4. Call sites** (`free()` + both `solve()` overloads): now call `error_message`
only on the error path (`INFO(1)<0`) and loop over `warning_message` line-by-line
for warnings (`INFO(1)>0`).

**5. `free()` now returns full `INFO(1:40)`** (Codex caught that
`mumpstools_free_solver` wrote only the scalar `INFO(1)`, leaving `mInfo[1..]`
stale from the prior solve, so any free-phase message interpolating INFO(2)
printed garbage). Fixed the Fortran wrapper to `intent(out), dimension(40)` with a
separate `ierr` for `MPI_BARRIER` and a `forall` fill (mirroring
`mumpstools_solve`); C decl `int_t& → int_t*`; call site `mInfo(0) → mInfo.data()`.

## Verification

- `mpicxx -fsyntax-only` (full debug flags, `-Werror`) passes **with and without**
  `-DBELFEM_MUMPS` (latter exercises the call-site + warning-loop edits).
- Case coverage checked (no spurious `-59/-68/-82..-87`; intended omissions only);
  braces balanced 135/135.
- Codex audited messages (round 2 fixes applied) and the `free()` interop change —
  both high confidence, no blockers. Residual pre-existing note: `free()` assumes
  `gOccupied`/`aSolverID` valid before indexing, guarded in C++ by `mInitialized`.
- Fortran not standalone-compiled here (needs the MUMPS module); change mirrors the
  proven solve-wrapper convention. **Run a full `make reset && make <target>` to
  rebuild the static libs** before trusting at runtime (CMake doesn't refresh
  static libs on plain `make`).

## Status

Complete, uncommitted. No behavioral change to the solve path beyond richer/correct
diagnostics; `free()` diagnostics are now accurate.

## Addendum — convergence improved, but NOT from this work

After the `make reset && make` rebuild done to test the above, the user observed
warning **+8** ("iterative refinement exceeded ICNTL(10) steps") gone and better
convergence. Investigated whether our edits changed MUMPS behavior — they did not
(everything here is diagnostics + the JOB=−2 `free()` wrapper; nothing touches
factorization/solve settings).

Real cause: `init_defaults` sets `ICNTL(14)` (memory relaxation) `0 → 30`
(`cl_SolverMUMPS.cpp:104`) — the documented `dl20260620` fix for the `-9`/`+8`
symptoms. Git state proves it: **HEAD = 0**, **staged-in-index = 30** (from a prior
session, never committed), our message edits unstaged on top (`git status` = `MM`).
Because **plain `make` does not rebuild static libs** (known BELFEM build quirk),
that staged `=30` had never been compiled into the running binary; the `make reset`
for *this* task finally activated it. ICNTL(14)=30 → ~30% more workspace → less
dynamic pivoting → cleaner factorization → refinement converges within ICNTL(10)=9
steps → +8 disappears. The `ICNTL(14)=30` fix remains staged-but-uncommitted (not in
HEAD); flagged to the user to commit so it isn't lost.
