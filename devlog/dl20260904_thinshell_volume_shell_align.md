# Thin-shell volumes: flat rank-0 exchange keyed by shell id

**Date:** 2026-09-04
**Purpose:** fix the MPI abort `Element … of thin shell block … maps to facet ordinal 5190, but the shell has only 1456 facets` on the `undulator2d` deck (proc 4, np ≥ 5); serial ran clean
**Module:** fem/kernel (`Kernel::compute_element_volumes`)
**Method:** plan + Codex/Grok audit → code + Codex/Grok audit, both rounds at `gpt-5.6-terra`/`xhigh` and `grok-4.6`/`xhigh` (safety-boundary row: MPI collectives). Exchange: `tmp/ai_exchange/thinshell_volume_shell_align.md`, restricted diff `…_shell_align.diff`

## Root cause (confirmed by both auditors, high)

`compute_element_volumes()` (`src/fem/kernel/cl_FEM_Kernel.cpp:865`) computed thin-shell
layer volumes inside `for ( tShell : tMesh->thin_shells() )`, with one `share`/`receive` of the
facet areas and one `broadcast` of the per-block first element ids **per shell**. That assumes
every rank iterates the same shell list. It does not: on workers `ProtoMesh::create_thinshells()`
(`src/mesh/cl_ProtoMesh.cpp:1491,1498`) skips a shell whose sideset or layer blocks have no element
on that rank, because empty blocks are deleted before they reach the block map. The surviving
shells keep master order, with holes.

`share`/`receive` reuse two pair tags in FIFO order and `broadcast` is a collective matched by
invocation order (`src/comm/commtools.cpp:89-109`), so a worker missing one shell consumes rank 0's
data for the *previous* shell at each later iteration. The numbers fit exactly: the six undulator
tapes have 1456 LINE2 facets each; layer elements and layer edges share one id counter, so the
first-id stride between consecutive tapes is 1456 + 2·1456 = 4368, and 5190 − 4368 = 822 is a
valid ordinal. Equal facet counts made the received vector look right until the id arithmetic
failed. Had the debug assert not fired, rank 0's leftover `MPI_Ibcast` would have met the worker's
`comm_barrier()`; in release the wrong first id indexes the area vector out of bounds.

The id-arithmetic scheme dates from 2026-08-15 (60f81aeb); the worker-side shell skip from
2025-12-02 (a7a37ec9). Earlier decks always had every tape on every partition.

## Change (one function, thin-shell section only)

Rank 0 walks **all** its shells once and builds two flat payloads: every facet area concatenated in
shell order, and a header of `[ shell id, area offset, facet count, block count,
( block id, first element id ) × blocks ]` per shell. Both go out **before** the shell loop as two
ordered `share`/`receive` pairs — header by `share` rather than `broadcast` because it is
variable-size (Codex P1, CLAUDE.md's payload table). Every worker receives both unconditionally,
including a rank with no thin shell at all. Each rank then resolves its own shells and blocks
**by id**, never by position, and applies area × thickness to the elements it owns exactly as before.

Checks, all setup-tier `BELFEM_ERROR` (runs once): duplicate shell id on rank 0; empty or
facet-count-mismatched layer block (was a debug assert, and also guarded a
`elements()( 0 )` on an empty block); shell or block absent from the header; header record and
area slice inside the received lengths. The per-element ordinal check stays `BELFEM_ASSERT`.

Kept on the auditors' insistence: the `++tNumNotMyElements` branch for non-owned layer elements — it
sizes the aura-volume gather later in the function (Grok C1, P1). Serial path: no communication,
same arithmetic, two new abort conditions that no deck that ever ran can hit.

## Jury record

- **Plan round:** both *accept with changes*. Adopted: header via ordered `share`/`receive`;
  exchanges unconditional; header bounds errors; duplicate-id rejection; areas and header from one
  rank-0 walk; keep the non-owned count. Codex corrected F8 (`index_t` is 64-bit only under
  `BELFEM_INT64`; still ≥ `id_t`) and F7 (the precedent is global ordering, not literally one
  share). Grok flagged `collect_thin_shell_facet_ids` (`cl_FEM_DofMgr_BlockData.cpp:160-181`) as a
  pattern *not* to copy: it counts shells with a selected block but dumps every shell's facets —
  a latent write past `set_size` if the two sets differ. Not touched this session.
- **Code round:** both *merge*, no P0/P1. Grok P2s declined, unreachable from this producer:
  overflow-safe form of the header bounds check and a livelock on a corrupt stride word both need a
  header no rank 0 emits; `gNoID` as the block-not-found sentinel collides only with a first id of
  `UINT_MAX`. The unsigned `id − firstId` wrap in release is pre-existing and out of scope.
- Grok did not modify the tree (git status checked after both rounds).

## Evidence

- Syntax check of `cl_FEM_Kernel.cpp` with the tree's real flags (`-std=gnu++17 -Wall -Werror`
  from `cmake-build-debug/src/fem/kernel/.../flags.make`): clean.
- **Gate ran (Christian, same day): `undulator2d` at np 8, cold start and warm restart, both
  clean.** The abort is gone at a rank count above the one that failed. Verified to the
  end-to-end deck rung; no tape-volume comparison against serial was made.
- **`make check` green (Christian, same day)** — no regression in the existing suite.

## Owed

- ~~Rebuild and rerun `undulator2d` at the rank count that aborted~~ — done, np 8, cold and warm.
- An MPI regression that puts a *later* shell on a worker with no element of an *earlier* one;
  `make check` today exercises MPI only in `tests/comm` and `sparsempi`. A zero-shell worker tests
  the hang, not the shifted mapping (Codex).
- `collect_thin_shell_facet_ids` count/fill mismatch — a separate audit item.
