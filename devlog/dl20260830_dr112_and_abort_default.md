# DR-112, an abort policy that was right for the code and wrong for the workflow, and one destroyed line

**Date:** 2026-08-30
**Purpose:** record the mandatory-BDF-triple loader contract, why the `gComm.size()`-dependent
error reaction was implemented and then reverted, and the shared-checkout incident that closed
the session
**Module:** `src/fem/kernel`, `src/core`, `src/comm`

## DR-112 — the memdump contract

`load_memdump` accepted a dump with no `bdf_h` / `bdf_step_count` / `bdf_last_dt` triple by
cold-starting the BDF order ramp, so a run configured BDF2–5 silently re-anchored at order 1 and
said so only in a log line. A *partial* triple already hard-errored; the wholly-absent case was
the hole.

Christian's ruling — "from now on, all memdumps will have the bdf variables" — turned the triple
into a format marker rather than an order-dependent need, and `save_memdump` already wrote it
unconditionally. The loader became the writer's mirror: require the magnetic triple always, the
thermal `*2` triple whenever `mEquation2 != nullptr`. That deleted the order-dependent predicate
the jury had specified, including the configured-`mOrder`-versus-live-ramp distinction that was
the easiest part to get wrong.

**The plan round earned its place twice.** Grok found the thermal check sitting *after* the
magnetic restore, so a coupled dump missing only the `*2` triple would have restored half the
integrator and then aborted — falsifying the very "no half-restored BDF state" rationale used to
justify the placement. Fixing that carried its own trap: hoisting the flag broadcasts required
*removing* the originals, because two collectives where the ranks expect one is a deadlock.

Gate, run by hand because `make check` never reaches `load_memdump` (both auditors confirmed):
a current dump resumes; the same dump with the triple deleted is refused on **4 of 4 ranks with
0 re-anchoring**; `restart : false` starts fresh. The 4-of-4 is the rank symmetry the hoist
exists for, tested rather than argued.

## The abort default — correct reasoning, wrong conclusion

Christian asked whether `BELFEM_ERROR` should call `MPI_Abort` outside debug. Release already
did; the throw is the *debug* default. I proposed making the reaction depend on `gComm.size()`:
abort in parallel, throw in serial. The motivating class is real — any `BELFEM_ERROR` inside an
`if ( rank == 0 )` block throws on one rank and leaves the others in a collective.

It was implemented, audited across two plan rounds and three code rounds, and **reverted**.

The deciding fact is not in the source. Christian debugs parallel runs with one `lldb` per rank,
launched by `mpirun`, each in its own terminal with `bt` scripted. A throw stops *that rank's*
debugger on the failure with a live backtrace and leaves the peers inspectable. `MPI_Abort`
tears the job down and shows nothing — the failure being debugged, made invisible. The change
would have broken the primary parallel debugging method outright.

Both vendors endorsed the design across four audits. They were reasoning correctly from the
code, and the code does not know how anyone debugs. That is a real limit on what an audit round
can establish, and worth remembering the next time vendor agreement feels like confirmation.

The rationale is now recorded beside the flag in `assert.hpp` as an explicit "do not make this
rank-dependent", because the argument *for* changing it is persuasive on its own terms and will
be rediscovered.

### What survived the revert

Four bugs the round exposed, none touching the policy:

- `error_abort()` used `gComm.world()` = `mComms( 0 )`, and `mComms` is EMPTY between `MPI_Init`
  and the push in `Communicator::init` — a nested assert in debug, undefined behaviour in
  release, **on the path already handling an error**. Now `MPI_COMM_WORLD` guarded by
  `MPI_Initialized` / `MPI_Finalized`, falling through to `std::abort()`. The documented form in
  `core_usage_guide.md` had been right all along; the code had drifted.
- A failed `MPI_Init` reported itself through `BELFEM_ERROR`, whose reaction may query MPI —
  undefined after a failed init. The error handler could become the second fault. Three sites now
  terminate locally. The PETSc checks that follow a *successful* init deliberately keep
  `BELFEM_ERROR`: MPI is live there.
- `finalize()` never restored `mSize` / `mCommRank`, so `rank()` and `size()` kept describing a
  communicator that no longer existed.
- `print_errorbox` would print `proc 2147483647` before init.

`check_doc_claims` also gained a `[F]` row-count check. It verified `[P]` and `[W]` but not
`[F]`, because `[F]` is stated in words — "the `[F]` list is EMPTY" — rather than the `( n rows`
digit form the existing parse greps. That count went stale twice on 2026-08-29 while the `[P]`
check passed, so the green line was actively reassuring.

## INC-558 — the last five minutes

Closing out, I saw `M devlog/README.md`, ran `git diff <file> | grep '^[+-]'`, got no visible
output, concluded "touched, no content change", and ran `git checkout --`. It was not a no-op: it
carried another session's index line for `dl20260830_physics_stats_footer.md`, and
`checkout --` has no undo.

`git diff --stat` on the same file read `1 insertion(+)`. The contradicting evidence was on
screen. I acted on the weaker signal because it agreed with what I already believed — which is
the same failure this session spent the night cataloguing in other forms: a stale `libbelfem.a`
reporting ten passing decks with the validator absent from the binary, a `make check` that
predated the change it appeared to bless, a YAML file that parsed while silently deleting three
keys from the accepted contract.

Two rules, the second from the peer session that repaired it:

- **`git stash` is the reversible form of `git checkout --`**, and the correct default in a tree
  with concurrent writers. It costs nothing when you turn out to be right.
- **`git diff --stat` is the authority on whether a file changed**, not a grep whose output can
  be swallowed by formatting.

The devlog body was untouched, so the entry was reconstructable. A peer rebuilt it from the body
rather than the header and produced a better line than my reconstruction — and correctly stripped
the apology I had put *inside* the index entry, on the grounds that an index is a navigation
artifact and every future reader should not pay for our incident. That judgement was right.

## Status

DR-112 struck. `[F]` empty, `[P]` 16, `check_doc_claims` 37/37. Commits `ceaab7f0`, `792c10c9`,
`fd8e36c9`. The abort-policy change is not in the tree and should not be re-proposed.
