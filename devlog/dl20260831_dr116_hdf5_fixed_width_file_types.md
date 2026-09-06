# DR-116: fixed-width HDF5 file types, and a bool reader that overran its target

**Date:** 2026-08-31
**Purpose:** Record the DR-116 writer change, the byte-identity evidence behind it, and the reader defect it uncovered
**Module:** `src/io`

## Context

Asked one hour before a design freeze how risky DR-116 would be to fix. Backward
compatibility was explicitly waived; the standing constraint was "do not break the reader."

## What changed

`filetype< T >()` added to `hdf5_types.hpp` beside `datatype< T >()`. The integer cases
dispatch on `sizeof` and `std::is_signed`, not on the C type name: `long` is 8 bytes under
LP64 and 4 under LLP64, and plain `char` is signed on x86 but unsigned on ARM, so a
name-keyed table would bake in the very platform assumption the row exists to remove.

All four templated writers (`save_scalar_to_file`, `save_array_to_file`,
`save_vector_to_file`, `save_matrix_to_file`) plus `save_bool_to_file` now create datasets
with `H5Tcopy( filetype< T >() )` and pass the NATIVE type to `H5Dwrite` as the memory type
— which also closes the row's second, latent big-endian clause. `H5Tset_order` no longer
appears anywhere in `src/io`.

`long double` deliberately keeps the native layout: it is 16 bytes carrying 80 bits of
precision, so `H5T_IEEE_F64LE` would silently truncate it. No writer instantiates it
(`real` is `double`).

## The evidence that mattered

The change is a **no-op on disk**, and that is the point rather than a disappointment:

- For all 15 types with a `datatype< T >()` specialization, `H5Tequal( old writer type,
  filetype< T >() ) == 1` — measured against the shipped table, not a hand-written mirror.
- Files written by the old and new paths are byte-for-byte identical (`cmp`, 4368 bytes
  each; `h5diff -v` reports 0 differences).

Because the emitted bytes do not change, the reader receives an unchanged input by
construction. That is a stronger guarantee than "the reader still works."

## The gate in the register was a null gate

The row proposed "`h5dump -H` shows fixed-width file types" as its acceptance criterion.
That criterion **already passed before the fix**: `h5dump` prints the canonical name of the
datatype description, and `H5Tcopy( H5T_NATIVE_ULONG )` + `H5Tset_order( LE )` already
carries the `H5T_STD_U64LE` description. The check could not distinguish fixed from unfixed
and would have signed off on a no-op. Replaced in the row by the equality, byte-identity,
compile and end-to-end checks.

The general lesson: a gate whose pass state is identical before and after the change is not
a gate. Ask what output the criterion would produce on the unfixed tree before adopting it.

## The reader defect this uncovered

`load_bool_from_file` passed the FILE type as the `H5Dread` MEMORY type — the DR-82 pattern,
un-migrated because DR-82 covered only the four TEMPLATED loaders and this one is not a
template. HDF5 therefore copies at the FILE's element size into a one-byte `hbool_t`.

Not theoretical. Reproduced before fixing, with a canary struct around the target: a 4-byte
dataset read into a guarded `hbool_t` **clobbered 3 bytes past it**, silently, and still
returned the right value — so the corruption is invisible from the call site.

Fixed by mirroring DR-82's two rungs: `check_read_datatype< bool >()` after `H5Dget_type`,
and `datatype< bool >()` as the read memory type.

`load_string_from_file` was examined and is **not** defective — it sizes its buffer from
`H5Dget_storage_size`, so file-type-as-memory-type is self-consistent there.

## Gates run

All executed, none inferred:

1. Type-equality probe over the shipped `filetype< T >()` — 15/15.
2. Old vs new writer output — byte identical.
3. `-fsyntax-only` with the real `belfem_io` flags on `hdf5_tools.cpp`, plus forced explicit
   instantiation of all four writers and `save_array_to_file` over real/int/uint/ulong.
4. End-to-end against the real fixed code (`hdf5_tools.cpp` compiled fresh, linked against
   the prebuilt core/comm archives): bool round-trips true/false; the formerly-corrupting
   4-byte-into-bool case converts with 0 canary bytes clobbered; a double dataset read as
   bool is refused by the class guard with the intended message, exit 134.

5. `make check` — **run 2026-08-31 by Christian, passes.**

The build gate is load-bearing rather than assumed. The rebuilt `cmake-build-debug/test/test_io`
is stamped 02:01:29, later than every edited source (latest `hdf5_tools.cpp` 01:42:51), so the
binary demonstrably carries the change. And the suite reaches both changed bool functions:
`cl_HDF5.cpp:443` and `:458` wrap `save_bool_to_file` and `load_bool_from_file`, driven by
`ScalarBoolTrueRoundTrip` and `ScalarBoolFalseRoundTrip` (`tests/io/test_HDF5.cpp:320`, `:337`).

**What the suite does not do is discriminate.** Both bool cases write a one-byte bool and read a
one-byte bool — which passed before the fix as well. So `make check` is a regression guard for
this change, not evidence for the overrun repair. The discriminating evidence remains the canary
probe and the end-to-end run against the real fixed code. This is the DR-140 lesson in a
different subsystem: a green suite says the suite is green, and the separate question is whether
any case in it could have run red before.

## Left alone, deliberately

The vlen path in `cl_HDF5_Dataset.hpp` builds `H5Tvlen_create( datatype< T >() )` and uses
it as both file and memory type. Same class of issue on the writer side; its reader is
correct. Untouched and unmeasured — out of scope for a change made against a freeze
deadline.

## Working-tree note

`src/io/filetools.{cpp,hpp}` and `todo/debt_register.md` were modified by another session
concurrently with this work. `src/io` was clean when this session started; `filetools` gained
89 lines mid-session and the register was substantially rewritten (rows struck, DR-120
archived, `[W]` rows removed). Those changes are not part of this work and were not touched.
`check_doc_claims.py` reports 35/37 solely because that concurrent rewrite has not yet
updated its own `[P]`/`[W]` header counts — left for its author to finish.

## Closure

**DR-116 STRUCK and archived 2026-08-31** (Christian's ruling) under the DR-42/DR-49
exception: the design and code work is finished and the residue is a single run. The gate —
a real `make` plus the io ctest suite — survives in the archived row's gate column.
**Struck is not verified.**

One register error was made and corrected in this session: the first status-cell edit split
the row on `' | '` and wrote into the wrong cell, duplicating `P3` into the status column and
pushing the `blocking 1.0` value (`no`) off the end of the row. Caught by inspecting the row
against `git show HEAD` before archiving, and repaired as part of the strike. The lesson is
narrow and mechanical: a register row must be rebuilt cell-by-cell against a known column
count, never by positional assignment into a split whose length was only asserted as `>= 6`.
