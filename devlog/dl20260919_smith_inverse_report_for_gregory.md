# Report for Gregory: grouped terminals and the Smith diagonal in `updatekGeneratorsFromHomology`

**Date:** 2026-09-19
**Purpose:** Evidence pack from branch `Claude` (never merged) on one question: when a bracket
group of curves links the same conductor k times, is the k × per-crossing current a consequence
of using the transposed Smith diagonal where a pseudo-inverse is needed, and does the obvious
correction leave every existing deck unchanged? Prepared by Claude at Christian's request
(2026-09-18 evening) for Gregory's decision. Nothing here is a change request.
**Module:** `src/homology`
**Status:** gates on the final binary (9a4e65f) and both juries complete; see §5 for what is open.
**Addendum 2026-09-20:** §0. Gregory's single-crossing convention tested on both binaries and confirmed: under his contract every invariant factor is 1, the formula is exact by construction, and devel is right. What this report calls a defect is an unguarded precondition (`s` discarded, the k-crossing group not refused) that BELFEM's own messages steer users into; the branch patch is a generalization, not a fix. The sections below stay as written as the record.

## 0. Addendum 2026-09-20: Gregory's reading, tested

Gregory's response (relayed by Christian): the unit diagonal is a design assumption; a
translation-periodic model has one connected conductor, hence one cut plus the periodicity
cut, and ONE input/output terminal on any one conductor should be defined, as in his corc
1/3-twist example; he suspects five conditions were defined.

**What the deck actually had:** one condition, not five, but its bracket held all five
crossings of the conductor on the back plane, `input curves : [1,3,5,7,9]`. That is the
k = 5 group of §1, not five conditions (five conditions are refused by the count guard).

**His convention, run on both binaries** (five tapes, 4 mm air, the shipped mesh;
`input curves : [1] ; output curves : [2] ; amplitude : 30 A`; devel `26bfd28` and the final
branch binary, fresh cut stage in both):

| binary | cut stage | I_1 = 10 A, per crossing at z = L/2 | I_1 = 30 A, per crossing | Ampère around the cable |
|---|---|---|---|---|
| devel | passes (no thick-cut refusal) | 10.00, 9.84, 10.00, 9.83, 10.00 A | 29.99, 29.52, 29.99, 29.50, 29.99 A | 4.998 × I_1 (five crossings) |
| branch | passes, no jump-weight line (k = 1) | identical | identical | identical |

I_1 and U_1 identical at all four saved steps; nodal fields within 1.2e-7 relative, the
build-to-build class of §4. The per-crossing pattern at 30 A is, digit for digit, what the
branch gave for the five-crossing group at 150 A (§4).

**Consequences.**

1. Gregory's reading is confirmed: with one crossing per condition the old formula is exact,
   every crossing carries the conductor current, and the design assumption `Σ = unit` holds.
   The k² of §1 arises only when one bracket links one conductor k times. Whether that input
   is legal is the whole question; the algebra of §2 is not in dispute.
2. The thick-cut refusal at 4 mm (§2 item 5, and the reason the deck was called a thick-cut
   test case) was a consequence of the k = 5 generator, not of the mesh: the single-crossing
   deck passes the thin-cut stage on devel at 4 mm.
3. What sends a user into the k-crossing group is BELFEM itself. The count-guard message
   (`cl_Cohomology.cpp:1088-1101`) says "merge the terminals that belong to the same
   conductor into a single bracketed group", and `doc/input_file_reference.md` §"Lists,
   ranges, groups" says ids inside one `[ ]` "are electrically connected and act as a single
   terminal". For a chained periodic conductor the k crossings ARE electrically connected,
   so a user following either text writes the k-crossing group, gets no refusal, and k × the
   current.
4. `corc_periodic_bc`'s `[1,2]` is one crossing (the geo splits each tape at its midline,
   two sidesets per tape), so the shipped decks all follow the single-crossing convention
   and none is affected by any of the options below.

**Options, in order of intrusion (your call):**

- (a) Guard only: `smithForm` already returns `s`, the number of unit invariant factors;
  refusing `s < t` in `updatekGeneratorsFromHomology` with a message "define the condition on
  ONE crossing of the conductor, as in corc_periodic_bc" turns the silent k × into a setup
  error. `s = t` is necessary, not sufficient: an input linking once against an output linking
  twice gives `F = [[1,-2],[0,0]]`, a unit factor, and `A F ≠ D` without a word; the guard
  should also check that each constrained row pairs `D_i` with every chain some cocycle sees
  (0 on the invisible slave chain), the one part of the branch patch worth keeping. No formula
  change, no effect on any existing deck. Plus the message and doc corrections of item 3
  (outside the closed core).
- (b) The branch patch of §3: accept the k-crossing group with the 1/k jump. Only worth it if
  you want the group to be a legal input; it also lets the 4 mm mesh through with the group.
- (c) Documentation only.

**Blind referee round (2026-09-20).** At Christian's request Codex (gpt-6-astra, high) and Grok
(grok-4.6, high) refereed the disagreement blind: a clean worktree of devel `26bfd28`, this
report as first sent (without §0), the patch, the probe, the deck as run and Gregory's reply
verbatim; nothing of the single-crossing test. Both, independently: one condition whose bracket
sums five crossings (Gregory right in substance, wrong in letter); the unit-diagonal assumption
holds exactly when every nonzero invariant factor is 1, guaranteed by one crossing per
condition; with `[1]`/`[2]` devel gives I_1 per crossing and the k = 5 thick-cut refusal
disappears (predicted; measured above); "defect" overstated, missing guard; minimal change (a)
with the pairing check, not the patch; the count-guard message and the input reference lead a
user into the deck that was run. Corrections to this report from the round: `R Σ⁺ Q_` is a
generalized, not the Moore–Penrose, inverse (§2); the as-run README's "5 × 30 A per tape" was
a factor short (750 A per crossing at the 150 A amplitude). Record:
`tmp/ai_exchange/review_smith_blind.md`.

The deck, its README and generator have been changed to the single-crossing convention in
both checkouts (`input curves : [1] ; output curves : [2] ; amplitude : 30 A`); the branch
patch stays on the shelf.

## 1. The observation

`examples/rutherford` (new deck, on branch `Claude` as `30d808a`): five 4 mm REBCO tapes wound at
30° around a flat core, one fifth of a pitch with translated periodicity, so the five tape
sidesets are one conductor that crosses each end plane five times. The deck names the five
crossings on the back plane in one bracket, `input curves : [1,3,5,7,9]`, and imposes 150 A.

Measured on devel (`26bfd28`, binary `build/bin/belfem` of 2026-09-18 18:25), YBCO current
integrated across the plane z = L/2 per crossing, and Ampère on the nodal H around the cable:

| deck | crossings in the group | per crossing / I_1 | around the cable / I_1 |
|---|---|---|---|
| five tapes, 1.5 mm air | 5 | 5.00 | 25.0 |
| single tape, one pitch | 1 | 1.000 | 0.999 |

Same binary, same generator, same deck skeleton; the only difference is the number of curves in
the bracket. (Record: `dl20260918_rutherford_cable_generator.md`; the three-voice round that
established it: `tmp/ai_exchange/review_rutherford_saturation.md`, Codex astra xhigh and Grok 4.6
xhigh, with Codex's exact −∇φ integral −153.750 A at I_1 = 6.15 A as the independent measurement.)

## 2. The mechanism (source trace, confirmed by a probe on `fn_Smith.cpp`)

1. A thin-shell bracket group becomes ONE suggested chain, the sum of its curves'
   segments (`cl_Homology.cpp:377-386`). Each of the k curves is a loop linking the conductor
   once, so the coupling `F(j, u) = ⟨c_j, γ_u⟩` carries ±k in the conductor's class.
2. `updatekGeneratorsFromHomology` builds the constrained row as `D · R · Σᵀ · Q_`
   (`cl_Cohomology.cpp:1123, 1141-1143`; comment at `:1119`: "Compute A = D * F^{dagger} via
   Smith Normal Form"). `smithForm` maintains `Q_ F₀ R = Σ` (`fn_Smith.hpp:133-137, 166-171`),
   so `F₀ = Q Σ R⁻¹` and `F₀⁺ = R Σ⁺ Q_` is a generalized inverse (`F F⁺ F = F`; not Moore–Penrose, `R` and `Q_` are unimodular, not orthogonal). The code has `Σᵀ` where `Σ⁺` belongs. For a 0/1
   diagonal `Σᵀ = Σ⁺`, so every deck with unit invariant factors is right; for `Σ₀₀ = k` it
   multiplies by k instead of dividing, `A F₀ = k² D`. `s`, the number of unit invariant factors
   `smithForm` returns (`fn_Smith.cpp:193-195`), is discarded; only the rank `t` is checked
   (`:1131`).
3. Probe (`tmp/ai_exchange/rutherford_figs/smith_probe.cpp`, compiled against `fn_Smith.cpp`),
   F = [[k, −k],[0, 0]], D = [1, −1]:

   | k | Σ | s | A = D R Σᵀ Q_ | A·F₀ | pairing per crossing |
   |---|---|---|---|---|---|
   | 1 | diag(1, 0) | 1 | [1, 0] | [1, −1] | 1 |
   | 2 | diag(2, 0) | 0 | [2, 0] | [4, −4] | 2 |
   | 5 | diag(5, 0) | 0 | [5, 0] | [25, −25] | 5 |

4. Why no shipped deck showed it: every bracket group in `examples/` links each generator
   class at most once (one curve or sideset per conductor; several DISTINCT conductors in one
   bracket, as in `racetrack_usermat`'s 18 tape ends, give 18 entries of 1, not one entry of 18;
   `corc_periodic_bc`'s `[1,2]` are the two halves of one crossing). The Rutherford fifth-pitch
   deck is the first with one conductor crossing the plane k times inside one group.

5. Two facts the real matrices add (found on the branch, not visible in the probe):
   on a periodic deck the chain on the slave plane pairs with NO cocycle, so the solver's
   coupling is F₀ = [[0, −1],[0, 0]] for one tape and [[0, −5],[0, 0]] for the five-tape group,
   and the transformation has always projected D onto the row space of F and imposed the
   current through the master plane alone; and the coarse (4 mm air) five-tape mesh was refused
   as a thick cut because a generator scaled by 5 needs five disjoint cut sheets, which the
   two-edge periodic loops of the coarse outer air cannot host.

## 3. What was changed on branch `Claude` (commits 5fd1724, 5aa9f02, 1fd937e)

The fraction 1/k cannot live in an integer cocycle (`Cochain` coefficients are `int`;
`CutData::collect_coefficients` records ±1; `clean_spfa` certifies unit representatives), so it
goes into the jump:

- new `src/homology/fn_terminal_transform.{hpp,cpp}` (outside the closed core): with
  `Q_ F R = Σ` the exact pseudo-inverse row of condition i is `Σ_j (D R)_ij / σ_j Q_(j,:)`;
  the least common multiple k_i of its denominators is the condition's factor and k_i times the
  row its integer generator, which pairs k_i (not k_i²) with the group chain. For unit factors
  k_i = 1 and the row is the previous `D R Σᵀ Q_` entry for entry (tested against a replica of
  the old pipeline, free rows included). Refused at setup: a condition that couples to no
  generator; a pairing that is neither `k_i D` nor 0 on a chain no cocycle sees (on a periodic
  deck the slave-plane chain is such a chain, and the old formula already imposed the current
  through the master plane alone); a free row pairing with a chain a condition uses. The
  free-cut rows (left null space of F, kernel completion) moved over unchanged.
- `cl_Cohomology.cpp`: steps 5 and 6 replaced by the call; `mTerminalFactors` stored;
  `terminal_factor( c )` returns 1 for free cuts and before the transform ran. `cl_Homology.cpp:
  662-668` comment updated (pairs +k, jump 1/k).
- `cl_CutProcessor`: `mCutWeights( c ) = 1 / factor`, refused if not positive, logged when ≠ 1;
  `CutSet::create_duplicates` writes the weights into the duplicate-node sources (original node
  slot stays 1.0). The dof manager consumes them as the T-matrix entries of the duplicate φ
  (`cl_FEM_DofManager.cpp:1316`).
- `tests/homology/test_TerminalTransform.cpp`: seven cases, including entry-for-entry equality
  of the full n × n transformation with the old pipeline for unit factors, the periodic
  slave-chain shape, and the refusals.

Semantics implemented: the amplitude of a k-crossing group is the total through the group
(I_1 stays what `iv_results.csv` calls it; each crossing carries I_1 / k). The alternative,
amplitude = conductor current, is weights 1.0 with the unit generator plus a deck-contract
change; both are one decision away and neither is decided here.

## 4. Gates (branch binary from 9a4e65f, fresh cuts in every run; devel binary as above)

| gate | result |
|---|---|
| homology test binary; `make check-fast` | 19/19; 8/8 |
| five tapes, 4 mm air (devel: refused as thick cut) | passes the cut stage; at I_1 = 50 A: 10.00, 9.84, 10.00, 9.83, 10.00 A per crossing, Ampère 1.000 × I_1; at 150 A: 29.99, 29.52, 29.99, 29.50, 29.99 A |
| five tapes, 1.5 mm air (the first field solution's mesh) | at 150 A: 29.99, 29.51, 29.99, 29.51, 29.99 A per crossing, Ampère 1.001 × I_1, mean J/Jc 0.479 |
| single tape, one pitch, devel vs branch | I_1, U_1, I_2 identical; nodal fields to 8e-9 |
| `tapestack3d`, devel vs branch, 1 s (40 saved steps, 48 timesteps) | identical Δt history, no rejection on either; I_1 and U_1 identical to all printed digits at every saved step; phi, Jz, H within 1.4e-6 relative at every step; the 2e-5 maximum is on transverse J components 1e-7 of Jz; T_max 77 K on both |
| `corc_periodic_bc`, devel vs branch, 3 ms (30 saved steps, 308 timesteps, identical Δt history) | fields within 7e-11 through 0.5 ms; I_1 identical to nine digits at every saved step, U_1 to four digits at 26 of 30; from 0.63 ms the deck's Picard loop exits through its stall criterion every step (residual −20 to −30 dB against a 1e-7 tolerance) and two builds of the SAME devel source part there by the same amount (12 % of the fields at 0.7 ms) as branch and devel do; the solver stack's run-to-run floor is 6e-11 |

**The shipped deck, as written.** `examples/rutherford` regenerated with `python/main.py`
(default 4 mm air) and run on the branch binary exactly as shipped (150 A over 10 s, 15 s of
simulation, MUMPS, Newton, BDF5): 163 timesteps, none rejected, 151 saved frames. At 5 s
15.03 A per crossing, from 10 s on 30.0 A per crossing (Ampère 1.000 × I_1 throughout), J/Jc
peaks at 0.80 during the ramp and settles to 0.71, U_1 decays to 2e-10 V on the plateau. On
devel this deck stops at the cut stage.

## 4a. Two side findings for Christian

- `corc_periodic_bc` as shipped runs its Picard loop into the stall exit from t = 0.63 ms on
  (4 to 6 iterations, residual −20 to −30 dB, tolerance 1e-7, `stall tolerance : 0.001`), and
  from there the solution is not reproducible across two builds of the same source. Any
  regression gate on that deck has to stop before 0.6 ms or compare globals only.
- A rerun in a directory holding `<mesh>.bfm` loads the previous binary's cuts and weights and
  skips the cut stage.

## 5. What is still open

- O1 semantics of a k-crossing group's amplitude (total vs conductor current): Christian's and
  your call. O4 the voltage conjugate under weighted jumps: consistent by construction of the
  T-matrix, not independently gated. O2/O3 mixed multiplicities and non-unit free cuts: refused
  or left alone, no shipped deck needs them.
- Not tested: multi-rank runs of a non-unit deck (the weights travel with the hanging nodes
  through the distributor, Codex's reading of `cl_Mesh_Distributor.cpp:791`, not executed);
  voltage-driven non-unit decks; the manual-cut path (`CutProcessorManual` keeps its own unit
  weights, it carries no terminal transform).
- Authorization: protocol §7.1 was not lifted by you; the branch exists under Christian's grant
  for evidence only and will not be merged. The closed-core diff is 20 lines in
  `cl_Cohomology.cpp` (call site, member, accessor) and a comment in `cl_Homology.cpp`.

## 6. Records

`todo/smith_inverse_grouped_terminals.md` (plan, both juries folded in),
`devlog/dl20260919_smith_inverse_branch.md` (authorization record, commits, gates),
`tmp/ai_exchange/review_smith_inverse_plan.md` and `review_smith_inverse_code.md` (juries),
`tmp/ai_exchange/rutherford_figs/` (probe, figures, analysis script).
