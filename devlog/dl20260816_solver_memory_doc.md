# Solver Memory Doc: Two Blind Drafts, Four Corrections, One Live Defect

**Date:** 2026-08-16
**Purpose:** Record the blind-parallel documentation campaign for
`src/sparse/doc/solver_memory_and_compression.md` and the DR-79 defect it
uncovered
**Module:** sparse
**Round:** `tmp/ai_exchange/solver_memory_doc_{brief,codex,grok,verified_facts,crossaudit_grok}.md`

## The method (Christian's design)

Codex and Grok drafted the full page blind and in parallel from one
pre-registered brief; each then cross-audited the other's draft; Claude
unified under its own lead, having independently fetched and read the
MUMPS 5.7.3 Users' Guide §5.19 to settle every [VERIFY-MANUAL] marker.

The method earned its cost within the hour. Corrections that no single
pass would have caught:

- **Grok refuted the brief itself, twice:** workers do NOT carry the full
  mesh (partition + aura only — the "6 GB replicated per rank" model this
  session had repeated for two days is false, and the committed 8-proc
  RSS numbers prove it), and the "four OOM kills" count is not
  reconstructable from committed logs.
- **Codex's draft shipped the wrong knob** in its BLR example
  (`relative tolerance` — the refinement key — where the compression
  cutoff belongs; the cutoff has NO deck key today, which is now an
  explicit code-round deliverable) — caught by Grok's cross-audit.
- **Grok's draft documented a key that does not parse** (the same
  missing `compression cutoff`) — caught by Codex's cross-audit, along
  with a Messe et al. citation-anchor cleanup.
- **Claude's own reassurance was wrong:** the manual fact "CNTL(7)=0.0 =
  full precision" is true, but BELFEM never sends 0.0 —
  `MUMPS::initialize()` overwrites it with the 1e-8 cutoff under a
  mislabeled comment. Both cross-audits confirmed the corrected trace.
- **Grok corrected both drafts on the STRUMPACK estimate line:** it is a
  SUM over ranks, not per-rank; the true per-rank peak line is compiled
  out unless the library was built with flop counting. The doc's sizing
  workflow was rewritten around that.

## The defect: DR-79

The cross-audit turned a consistency itch into a live bug: every MUMPS
deck without an explicit `compression scheme : off` runs lossy BLR at
1e-8 *absolute* today ( ICNTL(35)=automatic + CNTL(7)=1e-8 ), and
`examples/garber` pairs that with a 1e-10 nonlinear tolerance — negative
headroom, the 2026-07-06 failure pattern shifted two decades. Registered
with the fix design ( AUTOMATIC ≡ off everywhere, cutoff flows only on
explicit blr, comment fixed, `compression cutoff` becomes a real key
with parser + both contract artifacts + copy ctor + synchronize ).

## What shipped

`src/sparse/doc/solver_memory_and_compression.md` — Codex's skeleton
(short rule, deck contracts, headroom-pairing table, pitfalls), Grok's
memory model (rank-0 overlay incl. the MC64 gather and global matrices;
worker submeshes; the committed imbalance table; verbosity levels per
library), the settled manual facts ( ICNTL(35) semantics, CNTL(7)
absolute-vs-relative asymmetry, no MUMPS leaf-size equivalent —
Christian's direct question, answered ), and the campaign observations
labeled as such, separated from committed measurements. Indexed in the
sparse doc README.

Next: the unification code round ( plan+audit → code+audit ) and the
MUMPS driver jury ( abs-vs-rel cutoff semantics, ICNTL(35)=2 explicit,
the headroom warning tier ).

## The unification round (same day)

Plan dual-audited; both audits substantive. Codex refuted the plan's F4
(SIX shipped decks run MUMPS compression-keyless, not one — Grok added
the per-deck headroom table: garber negative, costheta/circuit/
RLC_Circuit at zero, helix/sidecoating at two decades). Grok rewrote C1
into the STRUMPACK shape (BLR opts IN via its own case; default —
OFF/AUTOMATIC/any future enumerator — stays exact, so this defect class
cannot recur through a new enum value), hardened C3 (setter-validated
> 0 AND finite, declaration-order copy-ctor insertion — `-Werror` would
have caught the wrong spot the hard way — have-flag, synchronize
10→11 uints + 1→2 reals), and found two holes in the warning design
(library-gate it to MUMPS/STRUMPACK so a stray blr key on a petsc block
does not earn the 07-06 story; print UNGUARDED on rank 0, because `-v 0`
legally silences the info-gated banners and a correctness notice must
not vanish with them). It also exposed `set_mumps_blr` as a non-working
pre-init escape hatch (initialize() clobbers its epsilon) — documented
as such rather than repaired, zero callers.

Implemented same-session: C1/C2/C3/C4 as amended, contract artifacts
(new `compression cutoff` row + rewritten `compression scheme` row +
schema entries + the ungated-case paragraph now pointing at the
warning), tests (default pin, validation matrix incl. ±inf/NaN,
copy-preservation walk over every tunable member), doc-page banner
flipped to landed, DR-79 → implemented/reviewed. All TUs gated green
with -Wreorder; check_doc_claims 34/34.

**Open: the warning tier.** Grok formally challenged warning-only with a
three-row split ( unstated-default cutoff + blr + tight tolerance →
ERROR; stated cutoff with zero/negative headroom → ERROR; stated, thin
but positive → WARNING ), arguing the first two rows are not "a
memory-bound run with open eyes" but the 2026-07-06 failure re-armed.
Christian ruled the same day: zero-or-negative headroom errors,
thin-but-positive warns. Implemented with one unifier judgment stated
openly: the error fires on the EFFECTIVE cutoff — stated or inherited —
because erroring only on stated values would make stating a cutoff
riskier than omitting it; the message names the provenance instead.
Grok's row 1 ( unstated default, thin-positive headroom ) thus remains a
warning; its negative-headroom subcase — the re-armed 07-06 — errors.

## GMRES default + the assumption DR-35 warned about (same day)

Christian resolved docket item 5a before the jury convened: uniform
krylov-method semantics — STRUMPACK's `auto` now resolves to
factorization-preconditioned GMRES ( PREC_GMRES; on a healthy factor the
first preconditioned iterate is at machine precision, on a degraded one
GMRES minimizes where refinement provably floored ), `preonly` documented
as the strict-tolerance choice. Implemented and verified same day.

The verify's V3 corrected a claim made in the DR-76 round and repeated
when announcing this change: **STRUMPACK returns SUCCESS at the
iteration cap** — the iterative kernels return only the residual and the
sparse-solver callers discard it ( source-verified in tmp/STRUMPACK ).
The maxit cap is a cost bound, never a named failure. DR-35's closed row
had flagged exactly this at medium confidence five days earlier ( "the
source is not installed... only positive evidence is a string in the
binary" ) — the register's calibrated-uncertainty discipline did its
job; the DR-76 round's mistake was asserting past that flag. Four sites
carrying the never-true "fails fast with NO_CONVERGENCE" claim corrected
( both contract artifacts, the strumpacktools comment, the register row
annotation ), plus the session memory.
