# Ghost interface: the SWIP average flux as the operator of record (DR-168)

**Date:** 2026-09-29
**Purpose:** `h_ghost()` assembles the consistency and adjoint-consistency terms with a JUMP of
the two outward through-thickness slopes weighted by the harmonic mean, which is inconsistent at
every resistivity ratio but 3. Christian ruled on 2026-09-24 that the symmetric weighted interior
penalty (SWIP) average of Ern, Stephansen & Zunino 2009 is the operator of record. This plan lands
that in one code round: flip the slave slope operator, halve the weight, rewrite the two coefficient
tests to the average, add the manufactured-flux fixture (R7 of `test_gap_followups.md`), and
correct the doc and the kernel comment, including the three-way α-floor disagreement.
**Module:** `src/fem/maxwell` (kernel), `tests/fem` (battery), `src/fem/maxwell/doc` (guide)
**AIs involved:** Claude (plan + implementation), Codex + Grok (plan jury, code jury)
**Status:** IN PROGRESS — approved 2026-09-29 (Christian, quiz); R1-R5 done (code jury terra/high + Grok 4.7/high reconciled, no P0/P1; Ern et al. 2009 Eqs. 7/8/15/16 verified from the primary text, DOI 10.1093/imanum/drm050); R6 (a) private-tree gate owed by Claude, (b) `make check-fast` owed by Christian; close-out after

> **Scope guards:**
> - Only the flux operator `D−` sign and the flux weight change in code. The penalty `α`, its
>   regularization `k_reg`, the four block lines, the symmetry, the `PENTA6TS`-only restriction and
>   the ghost switch are untouched.
> - No `input.conf` key changes meaning: `eta` and `k_reg` keep their roles, so
>   `doc/input_file_reference.md` and `doc/input_schema.yaml` are not edited (O3 asks the jury to
>   confirm this reading).
> - Consumer inventory (§3.2): no shipped deck assembles the ghost, so no shipped result changes.
> - The thin-shell interface condition itself (flux continuity, Alves et al. 2022b) is not
>   re-derived here; the 2026-09-24 round established it and Christian ruled on it.

---

## 1. Current behaviour and how it fails

`h_ghost( Calculator *, TimestepMatrices * )`, `src/fem/maxwell/matrices/mt_maxwell_h.cpp:436-603`.
Per Gauss point, with `Em`/`Es` the master/slave edge interpolation on the interface and `Dm`/`Ds`
the per-layer slope operators:

| Item | Code today | Evidence |
|---|---|---|
| weight | `rho_harm = 2 ρm ρs / (ρm + ρs)`, guarded to 0 when `ρm + ρs < ε` | `mt_maxwell_h.cpp:473-474` |
| master slope | `Dm = ( Em_top q_top − Em_top q_bottom ) / hm` = `+∂z h_m` (master is the lower prism, +z is its outward normal) | `:558`, `:563`, `:567` |
| slave slope | `Ds = ( Es_bottom q_bottom − Es_bottom q_top ) / hs` = `−∂z h_s` (the slave's own outward normal) | `:559`, `:564`, `:568` |
| blocks | `Kmm += EᵀE α − EmᵀDm ρ − DmᵀEm ρ`, `Kms −= EmᵀEs α + EmᵀDs ρ − DmᵀEs ρ`, `Ksm −= EsᵀEm α − EsᵀDm ρ + DsᵀEm ρ`, `Kss += EsᵀEs α + EsᵀDs ρ + DsᵀEs ρ` | `:573-576` |

Writing `[[q]] = Em q_m − Es q_s` (interface jump) and `F(q) = Dm q_m + Ds q_s`, the assembled
residual is

```
r_m =  α Emᵀ[[q]] − ρ Emᵀ F(q) − ρ Dmᵀ [[q]]
r_s = −α Esᵀ[[q]] + ρ Esᵀ F(q) − ρ Dsᵀ [[q]]
```

so the bilinear form is `α⟨[[q]],[[v]]⟩ − ρ⟨F(q),[[v]]⟩ − ρ⟨[[q]],F(v)⟩`: the SWIP structure, with
`ρ F` in the role of the averaged flux. With `Ds = −∂z h_s`, `F = ∂z h_m − ∂z h_s` is the JUMP of the
z-slope. For an exact solution with continuous flux `φ = ρm ∂z h_m = ρs ∂z h_s` it evaluates to
`ρ_harm ( φ/ρm − φ/ρs ) = 2φ (ρs − ρm)/(ρm + ρs)` where consistency needs `φ`: consistent only on
the ray `ρs = 3ρm`, zero at equal ρ, and an O(1) interface error that in-plane refinement does not
remove (three independent derivations, `tmp/ai_exchange/review_dr168_ghost_operator.md`, 2026-09-24).

The SWIP average (Ern, Stephansen & Zunino 2009: weights `ω_m = ρs/(ρm+ρs)`, `ω_s = ρm/(ρm+ρs)`) is
`ω_m ρm ∂z h_m + ω_s ρs ∂z h_s = ρ_w ( ∂z h_m + ∂z h_s )` with `ρ_w = ρm ρs/(ρm + ρs) = ρ_harm/2`.
The 2026-03-19 devlog's own transcription of the paper's weight (`δ⁺ρ⁺ = ρ⁺ρ⁻/(ρ⁺+ρ⁻)`) is that
`ρ_w`.

Second item, the `(0,0)` penalty floor is stated three ways:

| Where | Says | Algebra with `k_pen = 2 km_reg ks_reg/(km_reg + ks_reg)`, `α = η k_pen` (`:516-518`) |
|---|---|---|
| code | `α = η k_reg` when both `k` vanish; `α → 2η( k_other + k_reg )` when one `k → ∞` | right: `k_pen = 2 k_reg²/(2 k_reg) = k_reg`; `k_pen → 2( k_other + k_reg )` |
| `src/fem/maxwell/doc/ghost_penalty_stabilization.md:70-71` | `α ~ 2·eta·k_reg` (wrong), `α ~ 2·eta·(other_k + k_reg)` (right) | |
| kernel comment `mt_maxwell_h.cpp:496-497` | `alpha ~ 2 k_reg` (wrong), `alpha ~ 2 * min( other_k, k_reg )` (wrong) | |

`GhostPenaltyCoefficient` asserts the code's floor (`tests/fem/test_InterfaceOrientation.cpp:1217-1221`).

**What the tests pin today.** `GhostConsistencyCoefficient` (`:1091-1148`) and
`GhostPenaltyCoefficient` (`:1150-1230`) drive a one-sided constant master field `C` (slave dofs
zero) and assert, with `g_m = ∫ Em_topᵀ C dS`, `g_s = ∫ Es_bottomᵀ C dS` (`ghost_interface_projections`,
`:967-1040`):

| rows | dofs | today (coded jump, weight `ρ_h`) |
|---|---|---|
| 0-2 | master bottom (far) | `+ρ_h/hm · g_m` |
| 3-5 | master top (interface) | `( α − ρ_h/hm ) g_m` |
| 6-8 | slave bottom (interface) | `−( α + ρ_h/hs ) g_s` |
| 9-11 | slave top (far) | `+ρ_h/hs · g_s` |

Rows 6-8 and 9-11 carry the slave sign through the adjoint term `Dsᵀ[[q]]`, so they go red under
the average, as the test's own comment says (`:1086-1090`).

**Bottom line:** the kernel is a symmetric interior-penalty form whose flux operator is the slope
jump instead of the weighted average; the fix is two sign lines and one weight, and the §2 table is
the authority for every test row (rows 6-8 change sign structure, not scale).

## 2. Mechanism

Flip the two `Ds` assignments so `Ds = ( Es_bottom q_top − Es_bottom q_bottom ) / hs = +∂z h_s`,
and replace the weight by `ρ_w = ρm ρs/(ρm + ρs)` with the same `ε` guard. `F = ∂z h_m + ∂z h_s`,
`ρ_w F = ω_m ρm ∂z h_m + ω_s ρs ∂z h_s` is the SWIP average, the four block lines are unchanged and
`K` stays exactly symmetric (the adjoint term flips with the consistency term).

For the one-sided constant field (`F = 0`, `[[q]] = C`) the rows become:

| rows | after (average, weight `ρ_w = ρ_h/2`) |
|---|---|
| 0-2 | `+ρ_w/hm · g_m` |
| 3-5 | `( α − ρ_w/hm ) g_m` |
| 6-8 | `−( α − ρ_w/hs ) g_s` |
| 9-11 | `−ρ_w/hs · g_s` |

(rows 6-8: `Dsᵀ C` on the slave-bottom dofs is now `−g_s/hs`, so `−α g_s − ρ_w(−g_s/hs)`; rows 9-11:
`Dsᵀ C` on the slave-top dofs is `+g_s/hs`, so `−ρ_w g_s/hs`.)

For a manufactured field that is linear in z in each layer with continuous flux (`[[q]] = 0`, so
`α` and `η` never appear): interface rows `r_m(3-5) = −ρ_w( s_m + s_s ) g_m = −φ g_m`,
`r_s(6-8) = +φ g_s`, far rows zero. The coded kernel gives `−2φ(ρs−ρm)/(ρm+ρs) g_m` (`−18/11 φ g_m` at
ratio 10, `0` at equal ρ), an average with the full `ρ_h` gives `−2φ g_m`: ratio 10 and ratio 1
both separate the three, ratio 3 does not (the coded form is accidentally exact there).

**Rejected alternative:** keep the jump as a deliberate different operator and stop citing SWIP.
Ruled out 2026-09-24 (Christian, quiz): no reading makes the jump consistent, and the doc's claim
that the terms vanish for the exact solution is false today.

## 3. Gap table

| # | Item | Class | Change | Citation |
|---|---|---|---|---|
| 1 | slave slope sign | (c) | `Ds( i, j ) = −Es( i, j )`, `Ds( i, m+j ) = +Es( i, j )` | `mt_maxwell_h.cpp:559`, `:564` |
| 2 | flux weight | (c) | `rho_harm` → `rho_w = ρm ρs/(ρm+ρs)`, same guard; `rwdS = rho_w * wdS` | `:473-474`, `:572` |
| 3 | kernel comment on the weight | (c) | "harmonic mean of resistivities" → the SWIP weight `ω_m ρm = ω_s ρs`, Ern, Stephansen & Zunino 2009 | `:473` |
| 4 | kernel comment on the α floor | (c) | three limits: both `k → 0`: `alpha = eta k_reg`; one SC, other metal `k ≫ k_reg`: `alpha → 2 eta k_reg`; one `k → ∞`: `alpha → 2 eta ( other_k + k_reg )` | `:496-497` |
| 5 | consistency test rows | (c) | rows 0-2 `+ρ_w/hm`, rows 9-11 `−ρ_w/hs`; the four-row oracle comment (`:846-858`), the section banner (`:841`), the linearity note (`:780-781`), the helper comment (`:861-863`) and the header comment (`:1086-1090`) rewritten in the same diff | `test_InterfaceOrientation.cpp:1139-1142` |
| 6 | penalty test rows | (c) | rows 3-5 `α − ρ_w/hm`, rows 6-8 `−( α − ρ_w/hs )`; `tAlphaRead = tRG/tGG + ρ_w/hm` | `:1194`, `:1196`, `:1212` |
| 7 | test helper | (c) | `ghost_rho_harm` → `ghost_rho_w` returning `ρm ρs/(ρm+ρs)` (O2) | `:865-869` |
| 8 | R7 fixture | (c) | new `TEST( InterfaceOrientation, GhostConsistencyManufacturedFlux )`, §4 R3 | new |
| 9 | doc §2 operator | (c) | state the bilinear form of §1 (`α⟨[[q]],[[v]]⟩ − ⟨ρ_w F(q),[[v]]⟩ − ⟨[[q]],ρ_w F(v)⟩`, `F = ∂z h_m + ∂z h_s`, `ρ_w F` the SWIP average `ω_m ρm ∂z h_m + ω_s ρs ∂z h_s`); replace "they vanish for the exact solution" by: the adjoint term vanishes on a continuous trace, the consistency term reproduces the single-valued flux `φ`, which is what makes the form consistent | `ghost_penalty_stabilization.md:40-44` |
| 10 | doc §2 α floor | (c) | `α ~ eta·k_reg` for both `k` tiny; add the one-SC/one-metal limit `2·eta·k_reg` (the case the regularization exists for); keep `:71` | `:70` |
| 11 | doc §2 history sentence | (a) | keep: symmetrized April 2026, ruled intentional 2026-08-28; add that the average replaced the jump (no date in source; the doc may carry it) | `:46-53` |
| 12 | `GhostElementContract` items 2-4 | (a) | unchanged: far columns still scale linearly in ρ (item 2), a continuous constant field still has `[[q]] = 0` and `F = 0` (item 3), symmetry kept (item 4) | `:698-810` |
| 13 | fixture header comment | (a) | "couple only through the rho_harm-scaled gradient terms" → `rho_w` | `:640-642` |
| 14 | debt register row DR-168, `test_gap_followups.md` R7/O2, `devlog/campaigns` | (c) | close-out | — |

### 3.1 Cross-cutting

- **Symmetry is the invariant to re-measure:** `GhostElementContract` item 4 measures
  `max|K − Kᵀ| < 1e-12 max|K|`; it must stay green after the flip (the adjoint term is built from the
  same `Ds`).
- **The SC-SC guard** (`ρm + ρs < ε → 0`) carries over unchanged to `ρ_w`; the `(0,0)` pair in both
  coefficient tests keeps its expectations (weight 0, floor `η k_reg`).
- **No thickness enters the R7 oracle** (`(q_top − q_bottom)/h` cancels the slope's `h`), so the
  single-thickness fixture `TS_TestGhostStack( 0.1, flip )` is sufficient; the register row's
  "unequal thickness" wording was superseded by the 2026-09-24 round.

### 3.2 Consumer inventory (ghost-on decks and readers)

- `grep -rn "nitsche ghost penalty" examples/*/input.conf`: 2 hits, `tapestack3d/input.conf:39` and
  `tapestack3d_quench/input.conf:60`, both inside `//` comments. Ghost is opt-in since 2026-09-01
  (`fn_FEM_ghost_switch.hpp`), so **no shipped deck assembles `h_ghost`** and no shipped result changes.
- Tests that assemble it: the `InterfaceOrientation.Ghost*` battery (`tests/fem/test_InterfaceOrientation.cpp`).
  `ThinShellSelection.NoGhostFacetAcrossABuffer` checks facet creation, not the kernel.
- Prose: `src/fem/maxwell/doc/ghost_penalty_stabilization.md` (the operator, §2); `README.md:65`
  names the kernel only; `doc/input_file_reference.md:257-322` documents the switch and the two
  parameters, not the flux weight.

## 4. Ordered steps

- [x] **R1 — kernel.** Items 1-4 of §3 in `mt_maxwell_h.cpp`. The new weight comment cites
  Ern, Stephansen & Zunino 2009 (O1 on equation numbers) and states the operator:
  `ρ_w ( ∂z h_m + ∂z h_s )` = `ω_m ρm ∂z h_m + ω_s ρs ∂z h_s`.
- [x] **R2 — coefficient tests.** Items 5-7, 13. Every rewritten expectation is derived in §2 and
  pre-registered in the code-jury entry before the round.
- [x] **R3 — manufactured-flux fixture** (R7 of `test_gap_followups.md`), `tests/fem/test_InterfaceOrientation.cpp`
  beside the two coefficient tests, reusing `TS_TestGhostStack`, `ghost_edge_dof`,
  `ghost_interface_projections`, `ghost_assemble`, `ghost_residual`, `ghost_expect_rows`:
  - field: in-plane direction `C = ( 0.7, −0.4, 0 )`, profile `1 + s z` per layer with the master
    (lower prism) on `z ∈ [−hm, 0]`, slave (upper prism) on `z ∈ [0, hs]`; `s_m = φ/ρm`,
    `s_s = φ/ρs` for a chosen `φ` (say `1e-6`, so `s h ~ O(1)` at the fixture's `ρ ~ 1e-8..1e-7`);
  - dofs: master bottom row `( 1 − s_m hm ) · ghost_edge_dof( edge, C )`, master top row and slave
    bottom row `ghost_edge_dof( edge, C )` (each side's own edge objects), slave top row
    `( 1 + s_s hs ) · ghost_edge_dof( edge, C )`;
  - the profile is written in interface-relative z; the fixture's nodes sit at `z = 0, h, 2h`
    (`cl_TS_TestStack.hpp:741-745`, master `[0,h]`, slave `[h,2h]`) and node z is never sampled —
    the dofs are built from the bullet above only;
  - premises asserted first, with the test's own lines (`ghost_interface_projections` forms the
    master trace only, `:980-1004`): `max_k |Em_k q_m − Es_k q_s|` on the interface `< 1e-12 |C|`
    (`[[q]] = 0`), `Em q_m` reproduces `C`, `mMaxEmBottom`, `mMaxEsTop` as in `ghost_coefficient_setup`;
  - pairs `( ρ, ρ )`, `( ρ, 10ρ )`, `( 10ρ, ρ )` with `ρ = 1e-8`, both `tFlip ∈ { −1, 1 }`;
  - oracle: rows 3-5 `−φ g_m`, rows 6-8 `+φ g_s`, rows 0-2 and 9-11 zero, all to
    `1e-12 · ( φ + α ) · mMaxG` absolute — `α ≈ η k_reg = 4e-3` at the fixture's `ρ` (`mPenalty =
    { 4.0, 1e-3, 1e-4 }`, `cl_IWG_Maxwell.cpp:54`) and its roundoff `ε α g` would beat a `φ`-only
    tolerance at `φ = 1e-6`;
  - discrimination stated in the test comment, per the mutation table in R6.
- [x] **R4 — doc and comment corrections** (items 9-11): land with the round, no jury on the text;
  `check_doc_claims.py` and `check_citations.py` run after; the guide is user-facing, so §2 gets
  the Codex language sweep (`gpt-5.6-luna`/`medium`) with the formulas and citations fenced off.
- [x] **R5 — code jury** on the own-files diff (`git diff HEAD -- src/fem/maxwell/matrices/mt_maxwell_h.cpp tests/fem/test_InterfaceOrientation.cpp src/fem/maxwell/doc/ghost_penalty_stabilization.md > tmp/ai_exchange/review_ghost_swip_average.diff`),
  Codex `gpt-5.6-terra`/`high`, Grok `grok-4.7`/`high` (plan-or-diff row; `src/fem/maxwell/matrices`
  is not on the astra row). Verification pass, reconciliation, `git status` after Grok.
- [◐] **R6 — gate.** (a) done 2026-09-29 by Claude on the private non-MKL tree with the final diff: `InterfaceOrientation.*` 38/38 including `GhostConsistencyManufacturedFlux`; the two mutation checks ran on the implementing agent's build of the same diff: `−0.818182`/`+0.818182` (sign-only, 1:10 / 10:1) and `−2` (weight-only), both red as predicted. (b) Christian's `make check-fast` (Blaze) owed. (a) On Claude's private non-MKL debug tree (`/opt/scls/debug`, Armadillo):
  `test_fem --gtest_filter='InterfaceOrientation.Ghost*'` green, and two single-line mutation checks
  on the fixed kernel, each recorded with the number it produced:

  | mutation | equal ρ | `( ρ, 10ρ )` | `( 10ρ, ρ )` |
  |---|---|---|---|
  | none (SWIP) | `−φ` | `−φ` | `−φ` |
  | `Ds` sign reverted, `ρ_w` kept | `0` | `−9/11 φ` | `+9/11 φ` |
  | weight reverted to `ρ_h`, sign kept | `−2φ` | `−2φ` | `−2φ` |
  | both reverted (today's kernel) | `0` | `−18/11 φ` | `+18/11 φ` |

  Equal ρ cannot tell the sign-only revert from today's kernel; ratio 10 can, and the expected red
  value for the sign-only mutation is `−9/11 φ`, not `−18/11 φ`. (b) Christian: `make check-fast` on the shared trees (Blaze). Both
  backends are needed because the battery's tolerances are relative and the fixture uses `data()`-free
  accessors only — (b) is the Blaze leg.

## 5. Open design questions

- **O1** — **RESOLVED 2026-09-29 → Eqs. 7 (weighted average), 8 (SWIP form), 15 (weights), 16 (γ_K = δ⁺δ⁻/(δ⁺+δ⁻), "2γ_K is the harmonic average"), read from the authors' preprint; DOI 10.1093/imanum/drm050 resolved.** Original question: equation numbers for Ern, Stephansen & Zunino 2009. The register row says "Eqs. 6-8 and
  15", the 2026-03-19 devlog says "eq. 2.16"; `literature/` holds no extraction of the paper and
  `doc/literature_references.md:150` lists it with "DOI not confirmed locally". Proposal: cite the
  work and the quantity (the harmonic flux weights), no equation number, until someone with the
  paper confirms one. Jury: if you have the paper, state the equation and the DOI.
- **O2** — **RESOLVED 2026-09-29 → renamed.** rename `rho_harm` → `rho_w` (kernel) and `ghost_rho_harm` → `ghost_rho_w` (test helper).
  The minimal-rename rule applies to fixes that keep a variable's meaning; here the value stops being
  the harmonic mean, so keeping the name would mislabel it. Proposal: rename, two definitions and
  their uses (`mt_maxwell_h.cpp:474`, `:572`; test `:865`, `:1133`, `:1189`, plus the comments named in gap item 5).
- **O3** — **RESOLVED 2026-09-29 → no input-contract edit (both auditors agree).** input contract. `eta` and `k_reg` keep their meaning and defaults; the assembled operator
  changes. Proposal: no edit to `doc/input_file_reference.md` / `doc/input_schema.yaml` (they
  document the switch and the parameters, not the flux weight). Jury: object if you read the
  contract rule as reaching this.

## 7. Definition of done

- [ ] R1-R4 landed as approved, R5 reconciled, R6 (a) run by Claude with the two mutation checks
  recorded, R6 (b) run by Christian.
- [ ] `todo/debt_register.md` DR-168 struck to `debt_register_closed.md` with the gate named;
  `todo/test_gap_followups.md` R7: the "unequal-thickness" wording and "no kernel edit before it"
  struck with the reason (thickness inert; the round ran), the landed fixture named, O2 resolved;
  `todo/test_gap_backlog.md:123` (T-24) `rho_harm·[h]/h` updated to the average;
  `devlog/dl20260929_ghost_swip_average.md` + `devlog/README.md`.

## 8. Audit trail

- Physics round (2026-09-24): `tmp/ai_exchange/review_dr168_ghost_operator.md` (Codex terra/high,
  Grok 4.7/high, blind, reconciled), brief `dr168_ghost_operator_brief.md`.
- This plan's jury and the code jury: `tmp/ai_exchange/review_ghost_swip_average.md`.
