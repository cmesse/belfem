# Branch `Claude`: does the Smith diagonal need to be inverted for grouped terminals?

**Date:** 2026-09-19
**Purpose:** Record of the experiment Christian asked for on 2026-09-18 evening: a quarantined
branch in which the cohomology module may be edited, to test whether the k × per-crossing
current of a k-curve bracket group disappears when the terminal transform uses a proper
pseudo-inverse of the Smith form, and whether the existing decks are unaffected. The outcome
is evidence for Gregory Giard, who owns `src/homology` and decides.
**Module:** `src/homology` (`cl_Cohomology`, `cl_Homology` comment, `cl_CutProcessor`,
`cl_CutSet`, new `fn_terminal_transform`), `tests/homology`
**Status:** COMPLETE 2026-09-19 — plan and code juries folded in, all gates on the final binary 9a4e65f passed; report `dl20260919_smith_inverse_report_for_gregory.md`; branch NOT merged. **Addendum 2026-09-20:** Gregory's single-crossing convention tested on both binaries and confirmed (devel exact at 4 mm); the k × current needs a k-crossing bracket group, which the count-guard message and the input reference steer users into; the patch is now option (b) of the report's §0.

## Authorization record (protocol §7.1, Codex finding 2, Grok C9)

`cl_Cohomology`, `cl_Homology`, `cl_SimplicialComplex`, `cl_Chain`, `cl_Cochain` and `fn_Smith`
are closed to AI edits; the protocol requires Gregory's named authorization, relayed by Christian,
recorded in a devlog. That authorization does NOT exist. What exists is Christian's instruction
in chat on 2026-09-18 (evening, before leaving for the weekend): create a branch `Claude`, give
it full write access to the cohomology module, test the Rutherford deck and then one or two
shipped decks, never merge `Claude` into devel or main, and write a report for Gregory, whose
call it is whether anything is adopted. This file is the record of that scope: branch `Claude`
only, evidence only, no merge. Every closed-core change on this branch is a proposal.

## What was changed (commits on `Claude`, base devel `26bfd28`)

- `30d808a` the Rutherford example (the uncommitted devel state, so the deck is reproducible)
- `44b56f4` the plan `todo/smith_inverse_grouped_terminals.md`
- `5fd1724` and follow-ups: `fn_terminal_transform.{hpp,cpp}` (new, outside the closed core):
  `T = D R Σ̃ᵀ Q_` for the constrained rows with Σ̃ the unit-normalized Smith diagonal, the
  invariant factor of each condition read off `(D R)_i`, refusal when factors mix, when a
  condition couples to no generator, or when the row does not satisfy `A F = k D` exactly; the
  free-cut rows and the kernel completion moved over unchanged. `Cohomology::
  updatekGeneratorsFromHomology` calls it and stores the factors (`terminal_factor( c )`, 1 by
  default). `CutProcessor` turns the factors into jump weights `1/k`, logs any weight ≠ 1, and
  `CutSet::create_duplicates` writes them into the duplicate-node sources (the original node's
  slot stays 1.0). `cl_Homology.cpp:662-668` comment updated to the new pairing contract.
  `tests/homology/test_TerminalTransform.cpp`: seven cases (unit factor equals the old formula
  entry for entry including the free rows; k = 2, 5 stay unit and report k; two conductors;
  2-D single chain; two conductors of different multiplicity in one group keep a free split;
  dependent conditions refused; input linking once against output linking twice refused).

## Code jury and the final form (9a4e65f)

Code jury (`tmp/ai_exchange/review_smith_inverse_code.md`, Codex astra/high, Grok 4.6/high) on
`git diff 26bfd28..1be1d9d -- src tests`. Neither disputes the algebra. Codex found the
factor rule wrong: for two independent conditions on conductors linked 2 and 3 times the Smith
form of diag(2, 3) is diag(1, 6) with both factors mixed into each row of D R, so "all touched
factors equal" refused a solvable case. The factor is now the least common multiple of the
denominators of the exact pseudo-inverse row `Σ_j (D R)_ij / σ_j Q_(j,:)`, and that row times
the factor is the integer generator (tested: factors 2 and 3, unit rows). Codex also showed an
inherited hole, a visible output chain with zero pairing whose current the kernel completion
then leaves to a free dof; a zero pairing is now accepted only on a chain no cocycle sees, and
a free row may not pair with a chain a condition uses (tested; the periodic slave-plane shape
still passes). Grok: header contract without the transpose, the Homology Stokes sentence off by
1/k, `a`-prefixed locals, no unit test of the weight plumbing (open; the gate log line and the
30 A are the evidence). All folded in; `make check-fast` 8/8, homology test binary 19/19.

**Gate hygiene, learned the hard way:** a rerun in a directory that holds `<mesh>.bfm` loads
the cached cuts and weights of the binary that wrote it and skips the cut stage; the
final-binary gates were rerun after deleting the caches (the log must show "Creating cuts").

## Evidence so far

Probe (`tmp/ai_exchange/rutherford_figs/smith_probe.{cpp,out}`, compiled against
`fn_Smith.cpp`): for F = [[k, −k],[0, 0]], D = [1, −1] the old `D R Σᵀ Q_` gives A·F = k²·D
(1, 4, 25 for k = 1, 2, 5); `s` = 1, 0, 0. Solver runs on devel: five-curve group 5.00 × I_1 per
crossing and 25.0 × around the cable; single curve 1.000 × (record:
`dl20260918_rutherford_cable_generator.md`).

Plan jury (`tmp/ai_exchange/review_smith_inverse_plan.md`, Codex astra/high, Grok 4.6/high):
both confirm the diagnosis; findings folded into the diff and the plan (see the plan's §4.0).

## Gates (branch binary built from 1fd937e; devel binary built 18:25:34 from the tree committed as 26bfd28)

**Five tapes, 4 mm air (the mesh devel refused as a thick cut), 0.3 s ramp to 150 A.** The
run passes the thin-cut stage and logs "cut 0 : its terminal links the conductor 5 times, jump
weight 1/5". YBCO current per crossing through z = L/2 (triangle/plane intersection):

| t (s) | I_1 imposed | per crossing (A) | sum / I_1 | Ampère at r = 15 mm / I_1 | ∫Hz dz on axis (A) |
|---|---|---|---|---|---|
| 0.1 | 50 | 10.00, 9.84, 10.00, 9.83, 10.00 | 0.993 | 1.000 | 9.98 |
| 0.3 | 150 | 29.99, 29.52, 29.99, 29.50, 29.99 | 0.993 | 1.000 | 29.93 |

Devel on the same deck (1.5 mm mesh): 5 × I_1 per crossing, 25 × I_1 around the cable. The
amplitude is now the total through the group; the tape carries I_1 / 5, the deck's intended
30 A at the end of the ramp; mean J/Jc 0.48, the profile of the first field solution
(`dl20260918_rutherford_cable_generator.md`) at the current it was meant for.

**Five tapes, 1.5 mm air (the mesh of the first field solution), 0.3 s ramp to 150 A.** Step 3,
I_1 = 150 A: per crossing 29.99, 29.51, 29.99, 29.51, 29.99 A at z = L/2 (0.993 × I_1 in sum),
Ampère 1.001 × I_1, ∫Hz dz on axis 29.9 A, mean J/Jc 0.479. Same numbers as the 4 mm mesh to
three digits.

**Single tape, one pitch (unit factor), 0.3 s ramp to 30 A, devel vs branch on the same deck
and mesh.** Step 1: I_1 = 10 A, 10.002 A through z = L/2 on both, U_1 = 3.913e-6 V on both,
I_2 = −0.02606 on both; nodal fields phi, H, B, J, J/Jc agree to a worst relative difference
of 1.2e-8 (Newton tolerance 1e-9), i.e. the unit-factor path is unchanged.

**`tapestack3d` (eight soldered tapes, solder-face terminals `[42:48]`, unit factor), 76 ms,
devel vs branch, same deck and mesh.** Identical adaptive timestep history (10 steps, same Δt
sequence, no rejection); final step t = 93.02 ms: I_1 = 0.237375900693 A and U_1 = 3.202e-9 V
on both; phi 8.8e-9, Jz 3.0e-8, Hx 3.8e-8 worst relative difference on the load-bearing fields;
the 3e-5 maximum is on the transverse J components whose magnitude is 1e-7 of Jz.

**`corc_periodic_bc` (three tapes, translated periodicity, `[1,2]` = two halves of one crossing,
unit factor), STRUMPACK, 0.31 ms, devel vs branch.** Identical timestep history through the
first save; step 1 (t = 0.1 ms): I_1 = 5.025721452501 A, U_1 = 1.297063264e-3 V, I_2 =
4.929260027481 on both; worst relative field difference 7.6e-11.

Final step (t = 0.31 ms): I_1 = 15.557679027584 A, U_1 = 1.291755521e-3 V, I_2 =
15.259071165983 on both; worst relative field difference 7.5e-11; identical Δt history.

**`tapestack3d` over 1 s, devel vs branch:** 40 saved steps, 48 timesteps with identical Δt history and no rejection; I_1 and U_1 identical to all printed digits at every step; phi, Jz, H within 1.4e-6 relative at every step (`compare_all_steps.py` in `tmp/ai_exchange/rutherford_figs/`).

**`corc_periodic_bc` over 3 ms, devel vs branch (in progress).** Identical to 7e-11 through the
first five saved steps (0.5 ms). From step 63 (t ≈ 0.63 ms) the deck's Picard iteration stops
after 4 to 6 iterations at a residual of −20 to −30 dB (tolerance 1e-7 = −140 dB, `stall
tolerance : 0.001`, so the stall exit accepts an unconverged step); the two runs print the SAME
residual to two decimals at every step through step 84 (t = 0.76 ms) and then part in the last
digit. By saved step 7 (0.7 ms) the fields differ by up to 10 % on a few dozen of the 219 294
nodes, on the tapes near r = 2.8 mm, while I_1 and U_1 stay identical to nine digits and the
per-tape J/Jc maxima and means agree to three digits. That is the signature of two trajectories
of an unconverged nonlinear iteration parting on roundoff, not of a different constraint. A
devel-against-devel control run in a fresh directory is running to settle it; the outcome is
recorded below when it finishes. Finding for Christian regardless: `corc_periodic_bc` as
shipped runs its Picard loop into the stall exit from 0.63 ms on.

**Attribution of the corc seed (02:50).** The devel binary is deterministic: a second run of
`build/bin/belfem` in a fresh directory is bit-identical to the first through 0.7 ms (0.0
difference on every field at every saved step). Devel's own source (`26bfd28`) rebuilt in a
fresh tree (`../belfem-devel-rebuild`, same 36 options, same compilers) differs from the shipped
devel binary by 6.5e-11 at the first saved step, the same size as branch-vs-devel (6.4e-11) and
as branch-vs-rebuilt (6.5e-11). The seed is build-to-build floating-point roundoff, not the
terminal transform. Running the rebuilt devel binary to 1 ms to show it parts from the shipped
one at the same step and by the same amount as the branch did.

**Closing comparison (03:35).** The rebuilt devel binary against the shipped devel binary, same
source, parts at the same place and by the same amount as the branch did: 6e-11 through 0.5 ms,
5e-6 at 0.6 ms, then 12 % on phi/Jz/H and 40 % on the transverse J at 0.7 ms (branch vs shipped
devel: 2.7e-6 then 7 % and 28 %; branch vs rebuilt: 2.5e-6 then 6 % and 12 %). I_1 and U_1 stay
identical to nine digits in every pairing. Two runs of the rebuilt binary differ from each other
by 6e-11 at the second step, so ~6e-11 is the solver stack's own floor. Verdict for
`corc_periodic_bc`: the branch is within that floor for as long as the deck is reproducible at
all; from 0.63 ms on, where its Picard loop leaves through the stall exit every step, the deck
is not reproducible across two builds of the same source, and no binary comparison can say more.
The devel control was stopped at 0.7 ms and the rebuilt run at 0.7 ms once this was established;
the devel/branch pair runs on to 3 ms for the record of the globals.

**Corc pair to 3 ms, final (10:25).** Both runs: 308 timesteps, identical Δt history, none
rejected. I_1 identical to nine digits at all 30 saved steps; U_1 identical to four digits at 26
of them and different in the fourth digit at 4 (the conjugate reflects the fields, which differ
by 15 to 57 % in the stalled-Picard regime, as two builds of the same source do). Nothing
beyond 0.6 ms on this deck can distinguish binaries.

**Same for `tapestack3d`:** shipped devel vs rebuilt devel differs by 1.5e-7 to 6.6e-7 on phi,
Jz, H and 2e-5 on the noise-level transverse J over the first six saved steps, the same sizes
as branch vs devel; I_1 and U_1 identical. The branch-vs-devel differences on both regression
decks are build-to-build roundoff.

**The example as shipped** (4 mm mesh from `python/main.py`, 150 A over 10 s, 15 s): runs to
the end on the branch binary, 163 timesteps, none rejected, 30.0 A per crossing on the plateau,
J/Jc max 0.80 during the ramp; refused on devel.

**All five gates were rerun with the final binary (9a4e65f) in directories cleared of `.bfm`,
`memdump.hdf5` and `iv_results.csv`, so the cut stage ran ("Creating cuts" and, for the five-tape
decks, "jump weight 1/5" in every log).** The numbers above are from those runs; the earlier
runs with the intermediate binaries gave the same values to the printed digits. `make
check-fast` on the final build: 8/8.

## Addendum 2026-09-20: Gregory's reading, tested

Gregory (via Christian): the unit diagonal is a design assumption; on a translation-periodic
model define ONE terminal pair on ONE conductor. The deck had one condition with all five
crossings in its bracket, which is the k = 5 group. Single-crossing deck (`[1]`, `[2]`, 30 A,
4 mm mesh) on devel `26bfd28` and the branch binary, fresh cut stage, `reg3/onecross_*`: devel
passes the cut stage without the thick-cut refusal; per crossing 10.00, 9.84, 10.00, 9.83, 10.00 A
at I_1 = 10 A and 29.99, 29.52, 29.99, 29.50, 29.99 A at 30 A on both; Ampère 4.998 × I_1; I_1, U_1
identical at all four saved steps; fields within 1.2e-7. Gregory is right for the convention he
describes; the report's §0 reframes the decision as whether a k-crossing group is a legal input,
with the guard-only option (a) listed first. The deck and generator were switched to the
single-crossing convention in both checkouts. `corc_periodic_bc`'s `[1,2]` are the two halves
of one tape crossing, so every shipped deck already follows the convention.

**Retraction (2026-09-20, Christian's question "so the inverted Smith concept is wrong?").** The
report's framing "defect in the formula" is withdrawn. Under Gregory's contract (one condition on
one crossing) every invariant factor is 1 and `Σᵀ = Σ⁺`; the transform is exact by construction.
What remains is an unguarded precondition (`s` from `smithForm` discarded, non-unit factors not
refused) and the count-guard message plus `doc/input_file_reference.md` steering users into the
k-crossing group. The mechanism section (what the code does with such a group) stands; the patch
is a generalization to non-unit factors, offered as option (b), not a fix.

**Blind referee round (2026-09-20, Christian's request).** Codex astra/high and Grok 4.6/high, in a
clean detached worktree of devel 26bfd28 with the original report, patch, probe, as-run deck and
Gregory's verbatim reply, nothing of the single-crossing test or the retraction. Both match the
sealed position on every row (one condition with five summed crossings; assumption holds for
`s == t`; devel exact with one crossing; "defect" overstated, missing guard; option (a) plus the
pairing check, not the patch; the count-guard message misleads). Corrections accepted: generalized
inverse, not Moore–Penrose; the as-run README's "5 × 30 A" was a factor short; corc README's
"each tape connects to itself" is false for the 1/3-pitch translated model (fixed in the README).
Citations verified; blind tree clean after both runs. `tmp/ai_exchange/review_smith_blind.md`.

**Gregory's ruling (2026-09-21, via Christian).** Confirms the single-crossing convention works as
he suggested, and states the rule: grouping terminals is not a way to model periodicity, it
artificially increases the coefficients of the cut-to-terminal pairing matrix and makes no sense
topologically; on a periodic problem always define exactly one terminal pair per connected
component. Landed as wording in `doc/input_file_reference.md` §"Lists, ranges, groups" (both
checkouts; outside the closed core). The count-guard message in `cl_Cohomology.cpp` still says
"merge ... into a single bracketed group" and is his to change.
