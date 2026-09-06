# Trace: when the periodic-aware greedy cleaner left the cohomology core

**Date:** 2026-09-02
**Purpose:** Answer Christian's question whether `src/homology` changed in the last week, after the cohomology algorithm with periodic boundary conditions "seems broken"; record where `Cohomology::clean_greedy()` went and how to get it back
**Module:** src/homology (read-only trace; cohomology core closed to AI edits since 2026-08-31)
**Verification:** git history and source read only. **Reviewed, not verified** — no deck was run.

## Answer

`Cohomology::clean_greedy()` was deleted on 2026-08-25 in the DR-107 hygiene round
(`dl20260825_dr107_homology_cleanup.md`, jury Codex + Grok on Christian's condition "if Grok and
Codex agree, we can clean up") and reached the repository on 2026-08-26 in commit `0db03e32`
("backup"). The whole `mFunClean` pointer dispatch went with it; `clean()` now calls
`clean_spfa()` directly. The removal predates the 2026-08-31 AI edit ban on the core.

It had been unreachable since 2026-07-15 (`c29ef8fe`): Christian's own scaffold set
`mFunClean = &clean_spfa` in all three constructors, keeping the greedy body as the one-line
switch-back option. That option is what the deletion took away.

What survives in `clean_spfa()`: `fire_node_coboundary()` reproduces the greedy firing pattern
(master + slave quotient adjacency), and `rectify_greedy_sweeps()` runs the greedy sweeps on
feasible generators. Two deliberate deviations from the deleted routine, both recorded in
`dl20260715_spfa_clean_implementation.md`: the second adjacency pass folds boundary nodes via
`is_flagged()` instead of `is_periodic()`, and slave *faces* are no longer flagged (only nodes
and edges).

Recovery, if wanted, is a plain history checkout of the pre-deletion file
(`git show 0db03e32^:src/homology/cl_Cohomology.cpp`); the pointer scaffold and the greedy body
are both intact there. Authorization is Gregory's (protocol §7.1).

## Other homology-adjacent changes in the window (2026-08-25 to 2026-09-02)

| commit | date | change | periodic relevance |
|---|---|---|---|
| `ac73de0d` | 08-25 | `PeriodicityFactory::match_edges`: pure half-cut edge pairs are now TIED (deferred pass) instead of left untied | high — changes which periodic edges are tied before the quotient is built |
| `0db03e32` | 08-26 | DR-107 bundle: `clean_greedy`, field constructors, `SideSetFactory`, `create_sidesets_2d/3d`, `_old` transform deleted; `tTransInv` real → int; Progressbar in `clean_spfa` | the int change alters the `kernelImage` path that picks free-cut generators (integer division truncates, real `floor` does not, for negative entries) |
| `04c33664` | 08-28 | `Homology::reorient_generators()`: 2-D sign flipped to −1 (3-D unchanged); autopins (`MaxwellFactory::find_autopins`) | every 2-D cut changes sign |
| `a78ee40e` | 08-27 | `CutSet`: comment only (weights are 1 at any order) | none |
| `c7458e64` | 08-31 | `classify_periodic_pairs` / `mPairVerdict` / `CutPairVerdict` removed (rejected Step 6c residue) | none by design — the table had no consumer |
| `461d6b34` | 09-01 | docs only | none |

The homology suite passed after `0db03e32` (`dl20260827_overnight_dr_runs.md`, three green
`make check` runs), so the removal did not break `test_CohomologyPeriodic`. That test covers the
periodic bar / infinite tape only; the corc periodic deck is not gated.

Same-day observation that may be the same symptom: G7 in `todo/thinshell_conductor_normal_field.md`
(corc_solder rerun, 2026-09-02) finds the 10 A φ jump sitting on a cut that hugs the periodic
planes, with mid-cell rings not winding. Recorded there as pre-existing relative to the autopin
campaign; whether it predates 2026-08-25 is not established.

## Addendum (same day): Christian's ruling and the rollback footprint

Christian: the 2026-08-25 / 08-26 commits "should not have happened"; the greedy cleaner is
important. Handoff with the error phenomenology: `tmp/ai_exchange/corc_solder_rerun_handoff.md`
(G7: the 10 A φ jump sits on a cut hugging the periodic planes; mid-cell rings do not wind).

Anchor: `00c81576` (08-25 16:46, "backup") is the last state before either commit for every
file below; nothing between it and `0db03e32` touched `src/homology` or the periodicity factory.

| file(s) | what a checkout of `00c81576` restores | what it would lose (landed after) |
|---|---|---|
| `src/mesh/cl_Mesh_PeriodicityFactory.cpp` | untied half-cut edges (pre `ac73de0d`) | the pure half-cut tie; its own gate showed the corc cap-corner antisymmetry gone with it (`dl20260825_halfcut_tie_fix.md`), so reverting brings that artifact back |
| `src/fem/maxwell/cl_MaxwellFactory.cpp` (4 guard hunks of `ac73de0d`) | nothing needed: the guards are inert once no half-cut edge carries an edge source | file heavily rewritten since (autopins); do not checkout, hand-edit if at all |
| `cl_Cohomology.{cpp,hpp}`, `cl_Homology.{cpp,hpp}`, `cl_CutFactory.{cpp,hpp}`, `cl_SideSetFactory.{cpp,hpp}`, `CMakeLists.txt` | `clean_greedy` + `mFunClean` switch, field ctors, `create_sidesets_2d/3d`, `Matrix<real>` transform, `_old` transform | Progressbar in `clean_spfa` (Christian's hunk), two Verbose demotions, the 2-D sign flip in `reorient_generators` (`04c33664`, 08-28; a no-op in 3-D) and its comment in `cl_CutFactory.cpp` |
| `cl_CutProcessor.*`, `cl_CutSet.*` | `classify_periodic_pairs` / `mPairVerdict` (consumer-less table) | removed 08-31 (`c7458e64`), outside the ruled window; the 08-25 bisect commit `ecb4d851` is on no branch. Recommend leaving as is |

Not in the ruled window but on the cut plumbing the same days, recommend keeping: `d0f4db7d`
(08-26, `uint16_t` hanging-source counter, was wrapping on 464-cut trunks) and `c6ef67ef` (08-26,
thin-shell facet nodes kept out of master re-derivation; `ConnectivityCalculator` trusts inherited
master/slave links on facets from a `.bfm` or the distributor). `ac73de0d`'s `rel_tol` default
1e-9 → 1e-10 restores an earlier regression and should stay.

The 09-01 test additions (`tests/homology/test_Cohomology.cpp`) use only `clean()`,
`updatekGeneratorsFromHomology` and `get_Generators`, all present at `00c81576`; no test references
a deleted symbol.

Not done here: the core files are under the 2026-08-31 ban, so the checkout waits for Gregory's
authorization relayed by Christian. No run before 08-25 of the corc_solder deck is on record to
show G7 absent, so the rollback is a hypothesis test, not a confirmed fix.

## Ruling (same session, later): nothing is restored

Christian, as Gregory Giard's mentor, authorized a core edit for the rollback; the full
`00c81576` restore was staged, then withdrawn after the analysis below, and the tree is back at
HEAD (`21a2bb12`) in `src/`. His ruling:

- **Keep the `Matrix<int>` transform** (`0db03e32` D4): the cochain coefficients can only be
  integers, so the narrowing was intended. The kernel routine is exact on either type; the two
  instantiations differ only in the quotient (true floor vs truncation toward zero), which can
  change the free-cut basis, never the constrained rows, which come from the integer Smith form.
- **Keep the half-cut ties** (`ac73de0d`): reverting them brings the cap-corner antisymmetry back
  (`dl20260825_cap_corner_defect.md`, `dl20260825_halfcut_tie_fix.md`).

Facts established with the corc-solder session ("corc-solder-96", Fable) reading the rerun files:
`clean_greedy` was dead from the commit that introduced SPFA (`c29ef8fe`, 2026-07-15 20:45; the
pointer was never flipped in any commit); the greedy rectifier inside `clean_spfa` and the double
`clean()` run are intact at HEAD; so `0db03e32` removed no executed path on the cut side. The
rerun's G7 (10.25 A jump on `cut_1`, a deformed cross-section annulus on the periodic planes and
outer boundary; no through-cell θ cut; mid-cell rings not winding) matches the "pairing-correct
DUST" cuts already recorded on 08-25 morning (`dl20260825_corc_interface_contract.md`) because
`reduce_complexPellikka` never runs with suggestions on. G7 predates both ruled commits; it is a
cohomology question for Gregory. The switch-back option for the standalone greedy (`c29ef8fe^`)
remains available from history if wanted.

## "The greedy was some sort of a fallback" (Christian's recollection, checked)

Both sessions independently: the fallback is the **dense θ write-back inside `clean_spfa`**, not
the standalone `clean_greedy`. Per generator (`cl_Cohomology.cpp`, HEAD): SPFA certifies
(infeasible → `BELFEM_ERROR` with the negative-cycle certificate and `error.exo`, :532-543);
feasible → `rectify_greedy_sweeps` (:548, the greedy is the primary rectifier); greedy stall →
dense θ write-back with a re-solve guard (:551-567, "never observed"). `clean_greedy` was only
the one-line pointer switch-back (`dl20260715_spfa_clean_implementation.md:12-15`) and was never on
a runtime path (no assignment outside the three constructors at any revision). No record of the
dense fallback firing on any deck after 07-16; the one run record says all 14 corc generator
instances rectified sparse (`dl20260716_census_density_pivot.md:79`).

Visibility loss since `0db03e32`: the stall WARNING (:558) and the per-generator summary (:612)
are `Verbose`, previously `Default`, so a default-verbosity run that took the dense path is
silent. Restoring the WARNING to `Default` is a one-word edit inside the banned core.
