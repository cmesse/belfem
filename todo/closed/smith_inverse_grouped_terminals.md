# Grouped Terminals With a Non-Unit Invariant Factor: Invert the Smith Diagonal

**Date:** 2026-09-19
**Purpose:** When one bracket group of curves (or terminals) links the same cohomology class k > 1
times, `updatekGeneratorsFromHomology` currently multiplies the constrained generator by k instead of
dividing, so the imposed current appears k times per crossing and k² around the conductor.
Mechanism of the fix: build the constrained generator from the **unit-normalized** Smith diagonal
(an integer unit cocycle), read the invariant factor of each condition off the Smith form, refuse a
condition that mixes different factors, and carry 1/k into the cut-set weights so the amplitude of
a group means "the total through the group". That meaning is what this campaign's own deck and
writer state (`examples/rutherford/README.md`, `python/rutherford/belfem.py`); the input contract
(`doc/input_file_reference.md:75-80`, `input_schema.yaml`) only says a bracket is one condition
and does not define how the amplitude scales with crossings (Grok plan audit C1), so O1 is open.
**Module:** `src/homology` (`cl_Cohomology` — closed core, opened for branch `Claude` only by
Christian's grant of 2026-09-18; `cl_CutProcessor`, `cl_CutSet`; new `fn_terminal_transform`)
**AIs involved:** Claude Fable (plan + code), Codex gpt-6-astra (plan and code audits), Grok grok-4.6
(third voice); Opus (build orchestration). Evidence round: `tmp/ai_exchange/review_rutherford_saturation.md`.
**Status:** CLOSED 2026-09-23 — the experiment served its purpose: the report went to Gregory, and branch `Claude` stays unmerged by design. Earlier status: ✅ COMPLETE 2026-09-19 (branch `Claude`, never merged). 2026-09-20: Gregory's single-crossing convention tested and confirmed exact on devel; this patch is option (b) of the report's addendum, a guard on `s < t` with corrected message and docs is option (a). Landed: `fn_terminal_transform`
(lcm-of-the-rational-row factor, unit generators, visibility and free-row checks), the call site
and `terminal_factor()` in `cl_Cohomology`, jump weights 1/k through `CutProcessor`/`CutSet`, 9
unit tests. Verified: five-tape Rutherford deck 30.0 A per crossing at 150 A on the 1.5 mm and
the 4 mm mesh (the latter refused on devel), single tape and `tapestack3d`/`corc_periodic_bc`
identical to devel at solver tolerance, `make check-fast` 8/8, all on the final binary 9a4e65f
with fresh cuts. Report: `devlog/dl20260919_smith_inverse_report_for_gregory.md`. Open: O1
semantics, O4 voltage conjugate, multi-rank and voltage-driven non-unit gates, a `CutSet`
unit test of the weight path (Grok 5).
into devel or main; the outcome is a report for Gregory Giard, who owns `src/homology`.

> **Scope guards:**
> - Nothing here changes how a cohomology generator is COMPUTED (SNF, reductions, `clean_spfa`);
>   only how the terminal-constrained generator is scaled and how its jump is weighted.
> - Every deck whose coupling matrix has unit invariant factors must produce a bit-identical
>   transformation matrix and identical cut-set weights (1.0). That is the regression contract.
> - Semantics of a k-crossing group: "amplitude = total through the group" (documented contract).
>   The alternative "amplitude = conductor current" is O1, decided by Christian/Gregory, and is a
>   one-line switch in R4. This plan does not decide it; it implements the documented one.
> - No merge into devel/main from branch `Claude` (Christian, 2026-09-18).

---

## 1. Current Behaviour and How It Fails

`Cohomology::updatekGeneratorsFromHomology` (`src/homology/cl_Cohomology.cpp:993-1235`):

1. Coupling `F(i, j) = ⟨c_i, γ_j⟩` between the SNF cocycle basis and the suggested chains
   (`:1000-1012`). A thin-shell bracket group of k curves is ONE chain, the sum of the k loops
   (`cl_Homology.cpp:377-386`), so a group whose curves are k crossings of one conductor pairs
   ±k with that conductor's class.
2. Condition matrix `D` (input +1, output −1; `:1105-1115`).
3. `auto [Q, Q_, R, R_, s, t] = smithForm(F)` (`:1123`) turns `F` into `Σ` in place. The
   elementary operations (`fn_Smith.hpp:133-137, 166-171`) maintain `Q_ · F₀ · R = Σ`, so
   `F₀ = Q Σ R⁻¹` and the pseudo-inverse is `F₀⁺ = R Σ⁺ Q_`.
4. The code builds `A = D · R · Σᵀ · Q_` (`:1141-1143`), i.e. `Σᵀ` where `Σ⁺` is needed. For a
   rectangular diagonal of 0/1 entries `Σᵀ = Σ⁺`, so every deck with unit invariant factors is
   correct. For `Σ₀₀ = k`, `A F₀ = k² D` on that class: the constrained generator pairs k with
   every single crossing. `s`, the count of unit invariant factors, is discarded; only
   `aNumConditions ≤ t` (rank) is checked (`:1131`).
5. `CutSet::create_duplicates` ties `φ_dup = φ_orig + Σ_c 1.0 · I_c` (`cl_CutSet.cpp:43-53`), and
   `IWG_Maxwell::set_currents` fixes abstract dof c to the amplitude (`cl_IWG_Maxwell.cpp:119`).

| Failure | Mechanism | Evidence |
|---|---|---|
| Imposed current realised k × per crossing, k² around the conductor | `Σᵀ` for `Σ⁺` at `cl_Cohomology.cpp:1141-1143` | probe `tmp/ai_exchange/rutherford_figs/smith_probe.out` (k = 1, 2, 5 → A·F₀ = 1, 4, 25); five-tape deck 5.00 × per crossing, 25.0 × Ampère; single-tape control 1.000 × |
| Silent: no guard sees it | `s` discarded, rank check only | `:1123`, `:1131` |
| Coarse meshes refuse with "thick cut" | a generator scaled by k needs k disjoint cut sheets | 4 mm Rutherford mesh refused, 1.5 mm accepted |

**Bottom line:** the transformation is a pseudo-inverse only for unit invariant factors; a
k-crossing group produces the factor k, and nothing refuses or corrects it.

## 2. Architecture: Why "Unit Generator + Weighted Jump" Is the Right Spine

The exact rational solution `A = D R Σ⁺ Q_` has the entry 1/k: no integer cocycle realises it, and
BELFEM's cuts are integer cocycles (`Cochain` with `int` coefficients, `CutData::collect_coefficients`).
The fraction therefore has to live in the jump, not in the generator. Two places can carry it:

- the cut-set weights (`tWeights` in `CutSet::create_duplicates`), consumed by the dof manager as the
  T-matrix entries of the duplicate φ node (`cl_FEM_DofManager.cpp:1316`; the writes at
  `cl_FEM_DofMgr_DofData.cpp:3088, 3493` copy node weights into edge/face sources, Grok C8), or
- the imposed value (`set_currents`, divide the amplitude by k before fixing the dof).

The weights are chosen (R4): the abstract dof then IS the group total that the deck names, `I_1` in
`iv_results.csv` keeps its documented meaning, and the conjugate `U_1` is the voltage of the whole
group. Both places give identical fields; O1 records the alternative.

The generator itself is built from `Σ̃ᵀ`, the transpose of `Σ` with every nonzero entry replaced by
1: `A_int = D R Σ̃ᵀ Q_` (u × n, so the product is c × n for a rectangular F). For F = [[k, −k],[0, 0]]
this is the row [1, 0] (probe, k = 1 row). `A_int` is an integer combination of the basis (R and Q_
unimodular) whose pairing with the terminal chains is normalized; its EDGE coefficients are not
thereby bounded, so the certificate and rectification stage (`Cohomology::clean()` after the
update, `cl_CutFactory.cpp:352`) stays as it is and whether a coarse mesh now passes it is an
experimental outcome. For unit factors `Σ̃ᵀ = Σᵀ`, so nothing changes for existing decks.

**Per-condition factor.** Row i of `A_int` comes from `v = (D R)_i`; the invariant factors it touches
are `{σ_j : v_j ≠ 0, j < t}`. If they are all equal, that is the condition's factor `k_i`; if they
differ, or none is touched, `BELFEM_ERROR` at setup (once per run: `doc/coding_philosophy.md:621`).
**Pairing check (Codex plan audit, finding 1, corrected on the real deck):** every entry of
`A_int F₀` in row i must be 0 or `k_i D_i`, and at least one must be nonzero. It is NOT `A F = k D`
entry for entry: on a periodic deck the chain on the slave plane pairs with no cocycle, so the
solver's coupling is F₀ = [[0, −1],[0, 0]] for the single tape and [[0, −5],[0, 0]] for the
five-tape group, D = [1, −1], and the old formula always projected D onto the row space and
imposed the current through the master plane alone (A F₀ = [0, −1]); the strict check refused
every periodic deck, found by the first gate run. The per-entry rule keeps that behaviour and
still refuses an input chain linking once against an output chain linking twice
(F₀ = [[1, −2],[0, 0]]: the entry −2 is neither 0 nor −1). A group whose curves cross two
DISTINCT conductors 2 and 3 times is a single summed chain with column (2, 3, 0): factor 1, the
group total is imposed and the split between the conductors is a free generator (3 c_a − 2 c_b),
which is the documented semantics; tested.

**Where the algebra lives.** A new pure function `terminal_transform( D, F ) → ( A_int, factors )` in
`src/homology/fn_terminal_transform.{hpp,cpp}` (outside the closed core), unit-tested on synthetic
matrices; `updatekGeneratorsFromHomology` calls it and stores the factors. The change inside the
closed file is then the call site plus a member and accessor.

## 3. Gap Table

| # | Gap | Where | Fix step |
|---|---|---|---|
| G1 | `Σᵀ` used as pseudo-inverse | `cl_Cohomology.cpp:1141-1143` | R2, R3 |
| G2 | invariant factors never inspected | `:1123-1131` | R2, R3 |
| G3 | cut-set weights hard-coded 1.0 | `cl_CutSet.cpp:53` | R4 |
| G4 | no regression for k > 1 | `tests/homology` | R5 |
| G5 | deck README said "5 × 30 A per tape"; measured is k × the amplitude per crossing, 750 A per tape at 150 A (Grok C5) | `examples/rutherford/README.md`, `belfem.py` | done 2026-09-19 |

## 4. Ordered Steps

- [x] **R1** — Plan jury (Codex astra/high, Grok 4.6/high) on this file; fold findings in.
- [x] **R2** — `fn_terminal_transform.{hpp,cpp}`: `TerminalTransform terminal_transform( const Matrix<int>& D, Matrix<int> F )` returning `A` (`aNumConditions × n`, integer) and `factors` (`Cell<int>`, one per condition), plus the free-cut rows exactly as today (left null space rows `t..n-1` of `Q_`, kernel completion). Refuse mixed factors with `BELFEM_ERROR`. Doxygen contract: `Q_ F₀ R = Σ`; `A_int = D R Σ̃ Q_`; `A_int F₀ = diag(k) D`.
- [x] **R3** — `Cohomology::updatekGeneratorsFromHomology` calls it; new member `Cell<int> mTerminalFactors` (size n, 1 for free cuts) with accessor `terminal_factor( c )`. No other change in the closed core. (after: R2)
- [x] **R4** — `CutProcessor` reads the factors at construction (`cl_CutProcessor.cpp:21-45`) into `Vector<real> mCutWeights` (`1.0 / factor`); `CutSet::create_duplicates( aMaxNodeID, aAbstractNodes, aCutWeights )` fills `tWeights( tCount ) = aCutWeights( c )` for each set bit c instead of the constant 1.0 (`cl_CutSet.cpp:43-62`, both the plain and the periodic branches use `tWeights`). `tWeights` is initialised to 1.0 and only the set-bit abstract slots are overwritten, the last slot (the original node) stays 1.0 (Grok D4); a non-positive factor is refused and a weight ≠ 1 is logged. Consumers: `cl_CutProcessor.cpp` is the only caller; `CutProcessorManual` hard-codes its own `tWeights(2, 1.0)` (manual cuts carry no terminal transform) and `InterfaceProcessor` copies hanging originals, so 1/k rides through (Grok C7). (after: R3)
- [x] **R5** — Tests: `tests/homology/test_TerminalTransform.cpp` with F = [[k, −k],[0, 0]], D = [1, −1] for k = 1, 2, 5 (A_int = [1, 0], factors = k, `A_int F₀ = k D`), a 2-D single-column case, a two-condition unit case (A = D F⁺ bit-identical to the old formula), and a mixed-factor case expecting the error. Registered in `tests/homology/CMakeLists.txt` (label fast). (after: R2)
- [x] **R6** — Code jury (astra/high, Grok/high) on the diff; fold findings in. (after: R4, R5)
- [x] **R7** — Gates, all with the branch binary: (a) `make check-fast`; (b) Rutherford five-tape, 1.5 mm air, `amplitude 150 A`: YBCO current per crossing = I_1 / 5 ± 1 %, Ampère = I_1, and the run log shows "jump weight 1/5" (the first executable test of the weight path, Grok D5); (c) Rutherford single tape: unchanged 1.000 ×; (d) `tapestack3d` and (e) `corc_periodic_bc` run with BOTH binaries (devel `build/bin/belfem`, built 18:25:34 from the working tree committed two minutes later as `26bfd28`, whose tier-6 change touches the calculator, not homology, Grok D6; and the branch's) for the same short interval: `iv_results.csv` rows and a checksum of the nodal `phi`, `Hx..Hz`, `Jx..Jz` fields at the last step identical to roundoff (MUMPS nondeterminism at 1e-12 tolerated); (f) the 4 mm Rutherford mesh: expected to pass the thin-cut stage now (unit generator) — record the outcome either way. README status and `belfem.py` comment corrected. (after: R6)
- [x] **R8** — Report for Gregory: `devlog/dl20260919_grouped_terminal_multiplicity.md` on branch `Claude` (probe, runs, diff, both semantics, what was not tested) and a copy for Christian. (after: R7)

### 4.0 Implementation Progress (updated 2026-09-19, gates)
R7 status: (a) `make check-fast` 8/8 on the branch build before the change, homology test binary
17/17 after; (b) five tapes 1.5 mm: 30.0 A per crossing at I_1 = 150 A, Ampère 1.001 × I_1, log
shows "jump weight 1/5"; (c) single tape: devel/branch identical to 1.2e-8; (d) `tapestack3d`:
identical timestep history, fields to 3e-8 on Jz/phi/H; (e) `corc_periodic_bc`: step 1 identical
to 7.6e-11, final step pending; (f) 4 mm five-tape mesh: passes the thin-cut stage, 30.0 A per
crossing at 150 A. Code jury (R6) on `git diff 26bfd28..1be1d9d -- src tests`: Codex 3 findings (factor rule
refused independent non-unit conditions → lcm of the exact rational row; visible chain with zero
pairing and a completion row pairing with it → visibility and free-row checks; naming/Cell),
Grok 5 (header transpose, Stokes sentence, naming, mixed-factor test superseded, weight
plumbing untested at unit level) — all folded in 9a4e65f, 19/19 homology tests; gates rerun on
that binary.

R2–R5 drafted and committed as 5fd1724 before the plan jury returned (Christian's autopilot
grant; ordering defect acknowledged, Grok D1: the jury is folded into the diff). Codex plan
audit: 7 findings; 1, 3, 4, 7 applied (exact-pairing check + regression test, safe default
factor and positivity guard, full n × n transform returned and tested against the old pipeline
entry for entry, CMake registration), 2 → O5, 5 → wording in §2, 6 → O4. Grok plan audit:
V1–V10 confirm the diagnosis; C1 (contract citation), C5 (750 A, not 5 × 30 A), C7/C8 (consumer
census, T-matrix site), C10 (Homology comment), D4 (last weight slot), D5 (weight path logged
and gated in R7), D6 (binary pinned), D7 (literature, O7), D8 (factors logged at Verbose)
folded in; C4: the probe lives in the devel checkout's `tmp/ai_exchange/rutherford_figs/`,
copied to this worktree and reproduced in the report; C9 → O5 with the grant recorded in
`devlog/dl20260919_smith_inverse_branch.md`.

## 5. Open Design Questions (not silently decided)

- **O1** — Semantics of a k-crossing group's amplitude: total through the group (documented; this
  plan) or conductor current (I per crossing; weights 1.0 and the factor applied nowhere, i.e. R4
  becomes "unit generator only"). Both are one line apart in R4. To Christian and Gregory.
- **O2** — Mixed-factor conditions: refuse (this plan) or realise with per-class weights (each cut
  gets 1/σ of its own class). Refusal first; nothing in the shipped decks needs it.
- **O3** — Whether to also refuse when `s < t` and NO condition touches a non-unit factor (a free
  cut with a non-unit factor). Free cuts carry no imposed value, so a non-unit free generator only
  affects the certificate stage; left as is, logged.
- **O4** — Voltage conjugate (Codex finding 6): with the weighted jump the abstract dof is the
  group total and its conjugate residual is `Σ_c w_c ×` (loop voltage), i.e. the voltage of the
  whole group per period; a voltage-driven k-crossing deck and an independent power check
  (`I·U` vs `∫ J·E`) are owed before the semantics are called verified. Not in this round's gates.
- **O6** — `cl_Homology.cpp:662-668` states the pairing contract; with this change it pairs +k on
  the group chain and the jump carries 1/k. The comment is updated on the branch (Grok C10).
- **O7** — Literature (Grok D7): Kaczynski et al., *Computational Homology*, for the Smith form
  (`fn_Smith.cpp:14`); Pellikka et al. 2013 eq. 4.5 for the meaning of the basis coefficients
  (integrals of h over independent 1-cycles, `literature/papers/topology/pellikka2013.txt:373`);
  `src/homology/doc/thin_cut_nonunit_rectification.md` for why a non-unit coefficient cannot be a
  thin cut. None of them treats a terminal that links its conductor k times; that scaling is
  BELFEM's own bookkeeping.
- **O5** — Authorization (Codex finding 2, Grok C9): protocol §7.1 requires Gregory's named authorization
  for the closed core; this branch works under Christian's grant of 2026-09-18 for branch
  `Claude` only, so R3 is a PROPOSAL carried as code for evidence, not a change request.
  Recorded here and in the devlog; Gregory decides.

## 6. Interface Design

```cpp
// fn_terminal_transform.hpp  ( belfem::mesh::cohomology )
struct TerminalTransform
{
    Matrix< int > T ;         // n x n: condition rows, then the free cuts ( left null space + kernel completion )
    Cell< int >   factors ;   // n entries: the factor of each row, 1 for free cuts
};
TerminalTransform terminal_transform( const Matrix< int > & aD, Matrix< int > aF );
```
`Cohomology::terminal_factor()` returns 1 when the transform never ran (empty cell) or the index
is out of range (Codex finding 3); `CutProcessor` refuses a non-positive factor before taking the
reciprocal. `fn_terminal_transform.cpp` is registered in `src/homology/CMakeLists.txt` (finding 7).
`Cohomology::terminal_factor( const uint aCut ) const` → `int` (1 for free cuts).
`CutSet::create_duplicates( id_t &, Cell< Node * > &, const Vector< real > & aCutWeights )`.

## 7. Definition-of-Done Checklist

- [x] R1–R8 ticked, both juries recorded in §8
- [x] `make check-fast` green on the branch build
- [x] Rutherford five-tape: 30 A per crossing at 150 A imposed (plane integral and Ampère)
- [x] `tapestack3d` and `corc_periodic_bc`: devel and branch outputs identical to roundoff
- [x] Report delivered; branch `Claude` NOT merged

## 8. Audit Trail

- `tmp/ai_exchange/review_rutherford_saturation.md` — the evidence round (pre-registration, astra
  xhigh, Grok xhigh, verification, reconciliation, probe).
- `tmp/ai_exchange/review_smith_inverse_plan.md` — R1 plan jury (Codex astra/high, Grok 4.6/high).
- `tmp/ai_exchange/review_smith_inverse_code.md` — R6 code jury (Codex astra/high, Grok 4.6/high), verification and reconciliation.
