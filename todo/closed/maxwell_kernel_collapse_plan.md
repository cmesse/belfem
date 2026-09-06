# Maxwell / MaxwellThermal Matrix-Kernel Collapse onto `calculator::MaxwellData`

**Date:** 2026-07-13
**Purpose:** Collapse the 38 near-identical matrix-assembly variants (25 `maxwell::h_*` in
`mt_maxwell_h.{hpp,cpp}`, 13 `T_h_*` in `mt_thermal_h.{hpp,cpp}`) into a small set of generic
kernels that delegate all per-integration-point material math to the new memoizing helper
`calculator::MaxwellData` (`aCalc->maxwell()`), then shrink the two IWG dispatch trees. The
acceptance criterion is **proven numerical equivalence, variant by variant** — these kernels
build the tangent/stiffness of the nonlinear HTS solve, where a wrong sign, angle convention,
temperature source, or missing clamp produces silently wrong physics, not a crash.
**Module:** `src/fem/maxwell`, `src/fem/thermal` (+ `src/fem/kernel` for MaxwellData extensions)
**AIs involved:** Claude/Fable (exploration + plan), Codex (audit pending), Grok (third voice pending)
**Status:** PLAN + wiring fixes in tree, kernels untouched — drafted 2026-07-13; Codex audit
round 1 **confirmed D1-D5**; Codex prose polish applied; Christian's decision rounds adopted
the Newton-split kernel pair (§6.0). D1/D5/D8/D9/D10/D11 are fixed; D2 is fixed in the
helper with the legacy thermal-Conductor dispatch residual removed at R12; D3/D6 policy nits
and D4/O1 material-policy work are queued in R2; D7 is waived. **O1-O10 are all closed**,
including O10's executed powerlaw reorder to `(normJ, T, normB, angleNxB)` with Codex+Grok
0 defects. **2026-07-14: R1 build+smoke PASSED** (first attempt caught D12 — allocate-tail
MaxwellData wiring dead for blocks; hoisted + pointer-defaulted same day) **and R2a-R2h
IMPLEMENTED + Codex-audited** ("no blocking defect"; both caveats fixed; Grok unavailable,
2x truncation). Remaining tails: R2d Controller-side diagnostic, D2 legacy
Conductor-dispatch residual (dies at R12). **2026-07-19 status pass:** R6 had landed
silently in `abb71c2e` (2026-07-15) — bulk `h_calc` + thermal `T_h_calc` live, but the
TS-alloy half hit **D13** (thin-shell sidesets never build MaxwellData → `maxwell()`
assert = the real cause of the 2026-07-16 double-corc failure, mis-blamed on D1 in the
devlog) and was reverted to legacy `h_alloy`/`h_alloy_t`. D4 closed (R2c in + audited);
T12 resolved (intpoint asserts restored). Next: §4.1 verification of the three live R6
flips (Christian runs); D13 sideset-MaxwellData design before R9; then R7+.
**2026-07-21: Christian LEAPFROGGED R6-R10** — `link_to_group`'s entire Conductor/ThinShell
material tree replaced by an algorithm-based pick over the new kernel triple
(`h_picard` / `h_newton_mu0` / `h_newton_mu`, 3-AI audited + fixed same day; Ferro
analogously `phi_ferro_picard`/`_newton`). Claude fixed three crossed ternary slots that
dispatched the mu0-hardcoded `h_newton_calc` (wrong M for non-mu0, missing dmudh for
field-mu; approved). ALL legacy `h_*` material kernels + `h_calc`/`h_newton_calc` are now
unreachable from dispatch — dead until the R12 deletion; the family-by-family R6-R10 flip
sequence is superseded, so §4.1 verification now applies to the whole leap at once
(baseline = pre-rewrite commit `ce6a0e8b`). **2026-07-21 (later): thermal half COMPLETE.**
Conductor/ThinShell now dispatch to `T_h_picard` (renamed from `T_h_calc`); Christian's
verification gate passed; all 13 legacy `T_h_*` bodies were deleted (`mt_thermal_h.cpp`
1339 → 55 lines; R11 [x], R12 thermal half done, §6.1 sign-split closed). Codex found the
Alloy dependency-bit defect D14, now fixed. `Controller::set_thermal_kernel` now shares
materials for thermal Buffer/Ferro blocks. Remaining: R12 final gate (§4.1 matrix,
≥2 ranks) + R2d/O11 tails. See R11/R12 notes + §3.2 D14.

**2026-08-09 currentness sweep — the collapse is CODE-COMPLETE; only the gate is left.**
Verified in tree at `df3d8f90`, **re-checked 2026-08-11** (counts refreshed; the structure is
unchanged):
- `mt_maxwell_h.cpp` is **495 lines** (478 at `df3d8f90`; the growth is the DR-64 signed
  history-scalar fix and its comment, not a new variant) and declares exactly `h_picard`,
  `h_newton_mu0`, `h_newton_mu`, `h_ghost` plus the two side-connector kernels
  `h_side_connector` / `h_side_connector_newton` (a later, separate campaign —
  `hex8tb_phase2_fem_wiring.md`). All 25 legacy `h_*` variants are gone.
- `mt_thermal_h.cpp` is **129 lines**: `T_h_picard` + `T_h_newton`. All 13 `T_h_*` gone.
- `cl_IWG_Maxwell.cpp` is 629 lines and `cl_IWG_MaxwellThermal.cpp` is 105; both dispatch
  trees are collapsed to algorithm-based picks with no material lookup at link time.
- D13's premise no longer holds: conductive layers are FEM **Blocks** (Christian's
  `DomainType::ThinShell` selection edit), and `Calculator::allocate` builds `MaxwellData`
  for every non-Air block present in the peer dof manager
  (`cl_FEM_Calculator.cpp:1111-1137`) — the sideset-MaxwellData design that R9 was waiting
  on was never needed.

Therefore **R6–R10 are ticked below as "superseded by the 2026-07-21 leap"**, not as
individually verified flips — the distinction matters because §4.1 was never run
family-by-family.

**R12 stays `[◐]`. Two of its sub-items are undone, not one:**
1. **The final gate** — §4.1 matrix on ≥2 MPI ranks, helix + corc. Tracked as DR-02 in
   `debt_register.md`. This is the *only* verification the whole refactor will get, precisely
   because the family flips were leapfrogged.
2. **The `src/fem/maxwell/doc/` kernel architecture note** — never written. Confirmed
   2026-08-09: no file under `src/fem/maxwell/doc/` or `src/fem/kernel/doc/` mentions
   `h_picard` / `h_newton_mu0` / `h_newton_mu` / `T_h_picard` / `T_h_newton`. The only prose
   describing the collapsed kernel set lives in devlogs, which are ephemeral and must not be
   cited as documentation. A reader of `src/` today has no account of why 38 variants became
   five, or of the §6.0 two-way algorithm pick that replaced the material tree.

(R12's other sub-items are done: the dispatch cutover, the 37 body deletions and
`get_thermal_calculator` are all confirmed gone; "strip the shadow harness" is moot since
R3 was cut and no harness was ever built.)

**2026-08-11 currentness sweep — both R12 sub-items are still open, and one now has tooling.**
The architecture note is still unwritten: no file under `src/fem/maxwell/doc/` or
`src/fem/kernel/doc/` names any of the five kernels (re-grepped this sweep; the only hits are the
sources themselves). What *has* changed is the gate's instrumentation — the DR-02 Gate A
assembly dump is in the tree and committed (`bc578b5e`), env-gated on `BELFEM_DUMP_SYSTEM` at
both `cl_FEM_Controller.cpp:921` and `:1577`, writing the assembled system before the solve.
That does not close DR-02 and does not revive the abandoned pre-collapse A/B — the legacy
kernels exist only in history, which is why Christian cut Gate A — but the dump is the mechanism
any future equivalence comparison would need, and it is no longer a "commented-out call site".

Tails: R2d's Controller-side converged-at-clamp diagnostic, and O11 kernels + validation
case. §6.1's sign-split decision (DR-01) is closed — see the R12 note.

> **Scope guards:**
> - IN scope: the 25 `h_*` and 13 `T_h_*` variants, the `mFunMKF` selection in
>   `cl_IWG_Maxwell.cpp:306-397,437-527` and `cl_IWG_MaxwellThermal.cpp:70-174`, and the
>   MaxwellData extensions required for equivalence (D1-D5 below).
> - OUT of scope: `phi_*` / `symmetry_*` / `background_*` / L2-projection kernels;
>   `MaxwellPostprocessor` rerouting (logged as O8, not part of this refactor);
>   the **thermal Newton tangent** (stays with
>   `todo/thermal_matrices_cleanup_and_newton_plan.md` §3/R4-R5 — this plan *supersedes only
>   its R2-R3 collapse mechanism*, replacing local `mt_thermal_h` helpers with MaxwellData
>   delegation); the rho/lambda argument-convention decision
>   (`todo/rho_lambda_argument_convention.md`, linked as O10).
> - Compatibility promise: bit-comparable assembly wherever the FP evaluation order is
>   unchanged; the only tolerated deviations are the two documented FP-reassociation /
>   bug-fix exceptions (§3.1-T2 mu0 placement for generic-routed materials, §3.1-T8 2D
>   defect z-coordinate). Everything else must diff to zero.
> - `h_ghost` stays a separate kernel (facet-penalty logic, no per-point material math);
>   only its `element_rho` input contract is in scope (O5).

---

## 1. Current Behaviour and Why It Needs Collapsing

Each Maxwell variant assembles, per integration point k (weights `w`, measure `dV(k)`,
edge-interp `E(k)`, curl `C(k)`, edge dofs `q`):

```
b   = -mu0 * E(k) * q            (bulk;  ts: bt = -mu0 E q, b = bt + bn via compute_bn)
j   =  C(k) * q
beta = {pi/2 | bj_angle | bn_angle}          (per family, see §3.1-T1)
rho  = mat->{rho | rho_powerlaw | rho_piecewise}( ..., T, [x,y,z,t] )
M() += trans(E)*E * ( w(k)*mu0*dV(k) )
K() += trans(C)*C * ( w(k)*rho*dV(k) )
HTS only: dKdx_times_x() += (Cᵀj)(Cᵀj)ᵀ * ( w(k)*(drho_dJ/norm_j)*dV(k) )   if norm_j > eps
element bookkeeping: rho_m += w*rho*dV ; V += w*dV ; save_resistivity(rho_m/V)
```

Thermal variants mirror this with `M += w Nᵀ(density·cp)N dV`, `K += w Bᵀ·lambda·B dV`,
`f += w Nᵀ·rho·|j|²·dV`, reading b/j/E/C from the **paired Maxwell calculator** and T from
the thermal element itself.

The 38 bodies differ only along five axes — bulk vs thin-shell, metal/alloy/HTS/UserDefined,
defect (x,y,z,t) or not, power-law vs piecewise, thermally coupled (`_t`) or not — exactly
the axes `MaxwellData`'s constructor already dispatches on via function pointers chosen once
(`cl_FEM_Calculator.cpp:78-238`). Every past convention fix (e.g. the β `std::abs` fix,
`todo/closed/beta_angle_normal_sign.md`) had to be applied at up to 9+ call sites; the
duplication is the standing defect.

| Failure mode of the status quo | Mechanism | Evidence |
|---|---|---|
| Convention drift across variants | 22+ hand-copies of the same loop | β abs() fix needed 9 sites (`closed/beta_angle_normal_sign.md`); `T_h` vs `T_h_metal` lambda routing already disagree (§3.1-T10) |
| Silent physics divergence generic vs specialized | generic `h`/`h_ts` carry UserDefined branches the specialized fns lack | `mt_maxwell_h.cpp:1999-2033` vs `:384-453` |
| Hidden per-point allocation | `Matrix<real> Coords(N(k)*X)` per k in every defect variant | `mt_maxwell_h.cpp:589,680,916,...`; also flagged as D-item in `thermal_matrices_cleanup_and_newton_plan.md` |
| 2D defect coordinates read out of bounds | `Coords(0,2)` on a 1×2 matrix (mX is n×dim) | `mt_maxwell_h.cpp:592-594` + `cl_FEM_Calculator.cpp:651-657`; see §3.1-T8 |

**Bottom line:** the material math is one algorithm with five dispatch axes; it belongs in one
place (MaxwellData, already built and memoizing), and the kernels should reduce to the pure
FEM accumulation loop shown above.

---

## 2. Architecture: Why MaxwellData Is the Right Spine

- **Dispatch cost model matches the legacy design.** Legacy selects a specialized function
  once per group (`link_to_group`); MaxwellData selects member-fn-pointers once per block at
  `link_maxwell`/`allocate` time (`cl_FEM_Calculator.cpp:78-238`). Per point the cost is one
  indirect call — the same zero-abstraction budget (`doc/coding_philosophy.md`), plus
  memoization (`mLastIndex` keyed on k, `cl_FEM_Calculator.hpp:2307-2317`) that legacy code
  *lacks*: today `T_h_*` recompute b/j that the Maxwell pass already computed.
- **Reset is already wired.** `Calculator::link(Element*)` fires
  `mFunResetDataContainer` → `MaxwellData::reset()` (`cl_FEM_Calculator.cpp:1244`,
  `cl_FEM_Calculator.hpp:2294-2305`); kernels never call `reset()`.
- **Literature precedent.** This is the "equation object dispatches by domain type"
  architecture of Messe et al. 2023 (paper1), §2.5, applied one level down: the material
  database (§2.6, power law after Rhyner — see `powerlaws.hpp:94-95`) becomes the single
  source of ρ(J,B,β,T). The Newton term (Cᵀj)(Cᵀj)ᵀ·(∂ρ/∂J)/|J| and the ε<10⁻¹¹ iteration
  strategy it feeds are paper1's Quasi-Newton scheme (Eq. 12-13 and the checkerboarding
  discussion); the T-coupling into ρ follows the magnetodynamic h-φ coupling of
  Arsenault et al. 2023 (paper3), Section II.
- **Rejected alternative — template/CRTP kernel generation:** would keep the combinatorial
  variant surface (one instantiation per axis combination), cannot share memoized b/j/T with
  the thermal kernels, and duplicates the dispatch the helper's ctor already performs.
- **Rejected alternative — collapse thermal onto local `mt_thermal_h` helpers**
  (the pre-MaxwellData plan `thermal_matrices_cleanup_and_newton_plan.md` R2): superseded —
  it would create a *second* home for the same material math the helper now owns.

---

## 3. Variant Inventory (Gap Table)

> **HISTORICAL (reconciled 2026-07-14).** Two same-day code changes postdate most citations
> in §3-§5: (a) the O10 powerlaw reorder — the T-bearing `rho_powerlaw`/`rho_piecewise`/
> `drho_*_dJ` overloads are now `( normJ, T, normB, angleNxB [, x,y,z,t ] )`; quotes of the
> old `( normJ, normB, angleNxB, T )` order in D4/T14/O1 describe pre-reorder code. Line
> numbers in `powerlaws.hpp`, `mt_maxwell_h.cpp`, `mt_thermal_h.cpp` remain valid (in-place
> reorder; line counts unchanged, spot-checked). (b) the D5/D8-D11 dispatcher work grew
> `cl_FEM_Calculator.{hpp,cpp}` — hpp citations from the original draft may sit up to ~26
> lines below their quoted position.

Axes legend: TS = thin-shell, D = defect(x,y,z,t), PW = piecewise, T = thermally coupled.
"Class": (a) = mechanically replaceable by the generic kernel + existing helper,
(b) = replaceable once an On decision lands, (c) = must be handled explicitly.
All Maxwell rows also do: `save_resistivity(mean)` (O5) and M/K accumulation.
Confidence: high on every row (all 38 bodies read this session).

### 3a. `maxwell::h_*` (mt_maxwell_h.cpp)

| Fn | Lines | Axes | β convention | T source | dKdx | Class | Notes |
|---|---|---|---|---|---|---|---|
| h_metal | 28-86 | bulk·metal | bj_angle (:72) | gTbulk | — | (a) | rho(gTbulk,B,β) :75 |
| h_ts_metal | 88-157 | TS·metal | bj_angle (:142) | gTbulk | — | (a) | compute_bn :123 |
| h_metal_t | 159-221 | bulk·metal·T | bj_angle (:208) | FEM :204 | — | (a) | peer link :181 → D5 |
| h_ts_metal_t | 223-295 | TS·metal·T | bj_angle (:283) | FEM :279 | — | (a) | |
| h_alloy | 297-335 | bulk·alloy | none | gTbulk | — | (a) | **outlier: no b/j at all** — rho(gTbulk) :325; helper laziness reproduces this (§3.1-T13) |
| h_alloy_t | 337-382 | bulk·alloy·T | none | FEM :368 | — | (a) | **also dispatched for TS alloy** (O6) |
| h_hts | 384-453 | bulk·HTS | **dummy π/2** (:427) | gTbulk | :444-449 | (a) | Newton guard norm_j>eps :445 |
| h_ts_hts | 455-538 | TS·HTS | bn_angle (:511) | gTbulk | :529-533 | (a) | |
| h_hts_defect | 540-616 | bulk·HTS·D | π/2 (:580) | gTbulk | :607-612 | (a)* | Coords :589, time :594; *2D → §3.1-T8 |
| h_ts_hts_defect | 618-707 | TS·HTS·D | bn_angle (:672) | gTbulk | :698-703 | (a)* | |
| h_hts_piecewise | 709-778 | bulk·HTS·PW | π/2 (:749) | gTbulk | :769-774 | (a) | |
| h_ts_hts_piecewise | 780-865 | TS·HTS·PW | bn_angle (:837) | gTbulk | :856-861 | (a) | |
| h_hts_defect_piecewise | 867-943 | bulk·HTS·D·PW | π/2 (:907) | gTbulk | :934-939 | (a)* | |
| h_ts_hts_defect_piecewise | 945-1037 | TS·HTS·D·PW | bn_angle (:1002) | gTbulk | :1028-1033 | (a)* | |
| h_hts_t | 1040-1114 | bulk·HTS·T | π/2 (:1082) | FEM :1091 | :1105-1110 | (a) | |
| h_ts_hts_t | 1116-1205 | TS·HTS·T | bn_angle (:1174) | FEM :1182 | :1196-1201 | (a) | |
| h_hts_defect_t | 1207-1293 | bulk·HTS·D·T | π/2 (:1254) | FEM :1263 | :1284-1289 | (a)* | extra intpoint assert :1231-1233 |
| h_ts_hts_defect_t | 1295-1391 | TS·HTS·D·T | bn_angle (:1353) | FEM :1361 | :1382-1387 | (a)* | |
| h_hts_piecewise_t | 1393-1472 | bulk·HTS·PW·T | π/2 (:1440) | FEM :1449 | :1463-1468 | (a) | |
| h_ts_hts_piecewise_t | 1474-1563 | TS·HTS·PW·T | bn_angle (:1532) | FEM :1540 | :1554-1559 | (a) | |
| h_hts_defect_piecewise_t | 1565-1651 | bulk·HTS·D·PW·T | π/2 (:1612) | FEM :1621 | :1642-1647 | (a)* | |
| h_ts_hts_defect_piecewise_t | 1653-1749 | TS·HTS·D·PW·T | bn_angle (:1711) | FEM :1719 | :1740-1745 | (a)* | |
| h_ghost | 1751-1922 | facet penalty | n/a | gTbulk / facet-avg T (:1765-1777) | — | (c) | **outlier: stays**; consumes clamped element_rho :1785-1790 (O5) |
| h (generic) | 1927-2190 | bulk·any | π/2 (HTS) / bj_angle (metal :2170) | thermal-optional :1978-1986 | :2159-2163 | (b: O1) | UserDefined reduced-overload branches :1999-2033, 2048-2083, 2100-2119, 2129-2148; **M·mu0 at end** :2183-2186 (§3.1-T2) |
| h_ts (generic) | 2192-2472 | TS·any | bn_angle (HTS :2275) / bj_angle (metal :2455) | thermal-optional :2258-2266 | :2443-2447 | (b: O1) | UserDefined branches :2283-2317, 2332-2377, 2384-2438; M·mu0 at end :2466-2469 |
| get_thermal_calculator | 2474-2492 | util | — | — | — | (c) | link + intpoint assert → absorbed by D5/link_peer |

### 3b. `T_h_*` (mt_thermal_h.cpp) — thermal builds M (density·cp), K (lambda), f (ρ|j|²)

| Fn | Lines | Axes | β for ρ-source | lambda routing | Class | Notes |
|---|---|---|---|---|---|---|
| T_h (generic) | 17-144 | bulk·any | π/2 (HTS :102-121) / bj (metal :127-129) | `depends(lambda,normB)` :81-88 | (b: O1,O2,O3,O4) | T clamp gTmin :72; rho clamp :137; density fallback :90; **no UserDefined branches** (brief correction — those are in T_h_ts only) |
| T_h_ts (generic) | 146-343 | TS·any | bn_angle (HTS :260) / bj (:328) | depends() :239-249 | (b: O1..O4) | UserDefined branches :267-319; bn precomputed once before loop :193-208 |
| T_h_metal | 345-413 | bulk·metal | bj (:399) | **hard-coded lambda(T,B,β)** :401 | (a)+(O2..O4) | clamps :395,:409; density(gTroom) :403 |
| T_h_ts_metal | 415-517 | TS·metal | bj (:503) | hard-coded :505 | (a)+(O…) | |
| T_h_alloy | 519-575 | bulk·alloy | none | lambda(T) :563 | (a)+(O…) | no b at all; **also dispatched for TS alloy** (O6) |
| T_h_hts | 577-645 | bulk·HTS | π/2 :640 | lambda(T) :632 | (a)+(O…) | |
| T_h_ts_hts | 647-751 | TS·HTS | bn_angle :737 | lambda(T) :739 | (a)+(O…) | |
| T_h_hts_defect | 753-826 | bulk·HTS·D | π/2 :819 | lambda(T) :808 | (a)*+(O…) | |
| T_h_ts_hts_defect | 828-935 | TS·HTS·D | bn_angle :916 | lambda(T) :918 | (a)*+(O…) | |
| T_h_hts_piecewise | 937-1005 | bulk·HTS·PW | π/2 :1000 | lambda(T) :992 | (a)+(O…) | |
| T_h_ts_hts_piecewise | 1007-1109 | TS·HTS·PW | bn_angle :1095 | lambda(T) :1097 | (a)+(O…) | |
| T_h_hts_defect_piecewise | 1111-1184 | bulk·HTS·D·PW | π/2 :1177 | lambda(T) :1166 | (a)*+(O…) | |
| T_h_ts_hts_defect_piecewise | 1186-1293 | TS·HTS·D·PW | bn_angle :1286 | lambda(T) :1275 | (a)*+(O…) | |

Dispatch trees today: `cl_IWG_Maxwell.cpp` Conductor `:306-397` (PureMetal/LookupAlloy/HTS
/default × thermal × defect × piecewise), ThinShell `:437-527` (same, LookupAlloy → shared
**bulk** `h_alloy`/`h_alloy_t` at `:454-464`), Ghost `:555-559`;
`cl_IWG_MaxwellThermal.cpp` Conductor `:70-118`, ThinShell `:119-174` (alloy → bulk
`T_h_alloy` at `:136-140`; the ThinShell branch dereferences the **Maxwell** calculator's
material for dispatch, `:121-128`).

### 3.1 Cross-cutting equivalence traps — verdicts

Every trap from the task brief, re-verified against the code this session. Verdict key:
**EQ** = proven equivalent as helper stands, **EQ\*** = equivalent conditional on a Dn fix,
**DIV** = divergent, needs work (Dn/On), **CHG** = accepted intentional change.

| # | Trap | Verdict | Evidence |
|---|---|---|---|
| T1 | β four-convention reproduction | **EQ\*** (needs D3) | bulk HTS π/2: `mBeta` init/reset (`cl_FEM_Calculator.hpp:144,2304`) passed as dummy (`:2590`) = `mt_maxwell_h.cpp:426-433`. bulk metal bj: `compute_rho_metal` `:2567-2582` = `:70-75`. TS HTS bn: `compute_rho_powerlaw_ts` `:2624-2645` (+n-currency assert `:2632`) = `:510-518`. TS metal bj: same `compute_rho_metal`, b already = bt+bn via `compute_b_mu0_ts` = `:140-145`. All four exact **unless** the shared β slot is poisoned (D3). |
| T2 | mu0 sign & M placement | **EQ / CHG(FP)** | `compute_b_mu0_bulk` `hpp:2493-2499` = `-mu0·E·q` (`mt_maxwell_h.cpp:65`); `compute_b_mu0_ts` `hpp:2509-2520` = bn + bt (`:132-135`); `compute_bn` writes bn *and* n (`hpp:2243-2286`). M-mu0: 22 specialized fns multiply per term (`:79`), generic h/h_ts at end (`:2183-2186,2466-2469`) — algebraically equal; collapsed kernel uses **per-term** (majority form) → materials currently routed via generic get an FP-reassociation-level change (documented tolerance exception). The `compute_b_mu_*`/`compute_b_h_*` variants are currently unreachable generalizations: base `Material` ctor sets constant mu = mu0 (`cl_Material.cpp:161`); only bh-curve materials change mu (`:777-779`) and those are Ferro → phi kernels. |
| T3 | 2D j special case | **EQ** | helper j is 1-component in 2D only when it *creates* the vector (`link_vector`, `hpp:2178-2193`); IWG pre-creates "j" at mesh dim (`cl_IWG_Maxwell.cpp:590`) but legacy `j = C*q` resizes to C's rows anyway. `bj_angle_2d` returns π/2 and sets both norms (`hpp:2164-2173`). `norm_j` and Ctj dims match legacy. |
| T4 | Newton correction HTS-only | **EQ** | non-HTS `mFundRhodJ = return_zero` (`cl_FEM_Calculator.cpp:134,168,189,223`); with the §6.0 Newton split, non-HTS materials dispatch to `h` (no dKdx code at all) and HTS to `h_newton`, whose `norm_j > BELFEM_EPSILON` guard reproduces the legacy term exactly. |
| T5 | element_rho bookkeeping | **kernel-side (O5)** | `save_resistivity` mean is unclamped, clamp only at `get_resistivity` (`mt_maxwell_h.hpp:26-37`); sole consumer h_ghost (`mt_maxwell_h.cpp:1785-1790`); hidden from file output (`cl_MaxwellFactory.cpp:699-703`). Recommendation: keep the rho_m/V accumulation in the collapsed kernel using the same memoized `compute_rho(k)` values. |
| T6 | T clamp gap | **RESOLVED → O2 [x]** | `compute_T_fem` raw (`hpp:2402-2411`); thermal kernels use `max(·,gTmin)` (`mt_thermal_h.cpp:72,236,395`; gTmin=1.0, `typedefs.hpp:65`); Maxwell `_t` kernels use **raw** T (`mt_maxwell_h.cpp:204`). Decision (O2, 2026-07-14): unified clamp inside `compute_T` on both sides + `mTClamped` + zero dT-derivatives — deliberate (documented) upgrade of the raw-T maxwell side. |
| T7 | rho clamp gap (thermal f) | **RESOLVED → O3 [x]** | thermal kernels clamp into [gRhoMin,gRhoMax] before f (`mt_thermal_h.cpp:137,337,409-411,641-643`); helper raw (`hpp:2330-2339`); Maxwell K uses raw rho. Decision (O3, 2026-07-14): clamp inside `compute_rho` + `mRhoClamped`, one rho for all consumers — inert on the K path (`Material::mRhoMin = 1e-16` == gRhoMin, `cl_Material.hpp:255`); `compute_drhodj` returns 0 while clamped (expected dKdx shadow deviation on deep-subcritical elements). |
| T8 | defect (x,y,z,t) + 2D | **EQ / CHG(2D)** | helper: `compute_x` + ctor-bound `mTime` (`hpp:2716-2745,2794-2869`; time binding `cl_FEM_Calculator.cpp:68`); legacy: per-point `Coords` temporary + controller chain (`mt_maxwell_h.cpp:589-599`). 3D exact. 2D: legacy `Coords(0,2)` is an **out-of-bounds read** — mX is n×dim (`cl_FEM_Calculator.cpp:651-657`), so `N(k)*X` is 1×2; helper's mZ=0 (F7, accepted by Christian per `tmp/ai_exchange/maxwelldata_helper.md`) is the fix. Document as intentional non-bug-compatibility. Confidence high. |
| T9 | vector-aliasing side effect | **EQ, narrower than briefed** | mB/mJ/mBn/mBt/mN are references into the Maxwell calculator's work vectors (`link_vector`, `cl_FEM_Calculator.cpp:69-74`), so delegating keeps them fresh exactly when legacy did. Grep shows the only cross-file reader is `MaxwellPostprocessor`, which overwrites "bn"/"b" itself (`cl_MaxwellPostprocessor.cpp:637-656`). Alloy paths leave b/j stale — identical to legacy `h_alloy`. Confidence: high (grep), medium that no exotic consumer exists. |
| T10 | lambda routing (F6) | **EQ generic / DIV specialized — resolved toward generic** | helper keys on `depends(lambda,normB)` (`cl_FEM_Calculator.cpp:96-97,229-238`) = generic `T_h` (`mt_thermal_h.cpp:81-88`). Specialized `T_h_metal` hard-codes field-dependent lambda (`:401`) and HTS/alloy hard-code lambda(T) (`:563,632,739,...`). For PureMetal materials `depends(lambda,normB)` is true and HTS/alloy false ⇒ helper reproduces the specialized behavior *provided* the property flags match the type — verify per material in R2; any mismatch surfaces in the shadow compare (R3). Confidence: medium-high. |
| T11 | density gap | **DIV → O4** | thermal M needs density·cp (`mt_thermal_h.cpp:90-91,403-405`); helper has no density accessor (`hpp:177-214`). Also: generic `T_h` falls back `have(density) ? density(gTroom) : ref_density()` (`:90,252`) while specialized call `density(gTroom)` unconditionally (`:403,507,565,...`) — unify on the generic (fallback) form. |
| T12 | paired-calculator link + intpoint assert | **[x] RESOLVED 2026-07-19 — D8/D9 fixed, intpoint asserts restored in both coupled dispatchers (`cl_FEM_Calculator.cpp:2462-2500`)** | legacy: `get_thermal_calculator` links the peer element and asserts matching intpoint counts (`mt_maxwell_h.cpp:2474-2492`); thermal mirrors it inline (`mt_thermal_h.cpp:22-38,350-366`). Now handled inside `Calculator::link(Element*)` → `link_element_maxwell` (`cl_FEM_Calculator.cpp:1246,2336-2360`, Christian 2026-07-13); recursion proven bounded (depth 2, see D5); still open: D8 branch selection, D9 null guard, re-adding the intpoint assert. |
| T13 | h_alloy / h_ghost not the common shape | **EQ / out-of-scope** | `compute_rho_bulk` (`hpp:2538-2542`) never touches b/j ⇒ the collapsed kernel reproduces h_alloy's skip-b/j optimization *and* its stale-vector behavior for free. h_ghost stays (facet penalty; only element_rho contract via O5). |
| T14 | UserDefined dependency overloads | **DIV → D4/O1** | helper always calls the full overloads. The 8-arg defect `rho_powerlaw` **asserts** full (normB,angleNxB,T) jc-dependence (`powerlaws.hpp:229-238`); the 4-arg has a constants-fallback (`:199-209`) but evaluates ρₙ = rho(T) (`:211`) where the legacy 3-arg path uses rho(gTbulk) (`:156`). Legacy generic h/h_ts select reduced overloads by `depends(jc,·)` (`mt_maxwell_h.cpp:1999-2033,2283-2317`). Parity does NOT hold for partial-dependency UserDefined superconductors — fix in the helper ctor (O1). |
| T15 | no thin-shell alloy variants | **EQ — explained** | alloy rho is T-only, so b (the only TS-specific quantity) is never queried; the *bulk* kernel is correct on a TS block because E/C/dV come from the TS calculator regardless (`cl_IWG_Maxwell.cpp:454-464`, `cl_IWG_MaxwellThermal.cpp:136-140`). The collapsed design reproduces this by construction (helper laziness). |
| T16 | classification mismatch (enum vs have/depends) | **first-class → O1** | see O1 for the full mapping proof and the two residual divergences (UserDefined-SC routing, D4). |

### 3.2 Defect register (current tree, found this session — block delegation until fixed)

- [x] **D1 — FIXED 2026-07-13 (Christian, verified by Claude in tree same day).**
  `allocate()` excludes Air (`domain_type() != DomainType::Air`, `cl_FEM_Calculator.cpp:1062`)
  and the `link_maxwell()` early-return now also rejects Air
  (`cl_FEM_Calculator.cpp:1083`), closing the immediate-rebuild residual.

  **Original finding, preserved for audit history:** `MaxwellFactory::create_controller` calls
  `link_maxwell` on **every** magnetic block (`cl_MaxwellFactory.cpp:730-737`);
  `Calculator::allocate()` gates only on `mMaxwellKernel != nullptr && GroupType::BLOCK`
  (`cl_FEM_Calculator.cpp:1062`). Air blocks never get a material (`assign_materials` skips
  "air", `cl_MaxwellFactory.cpp:2368-2377`; `Group::mMaterial` defaults nullptr,
  `cl_FEM_Group.hpp:97`), yet the MaxwellData ctor binds `aCalculator->group()->material()`
  (`cl_FEM_Calculator.cpp:67`) and dereferences it (`:80,:85`). Fix: gate construction on
  `group()->material() != nullptr` (or domain type ∈ {Conductor, ThinShell}). Found by Claude
  2026-07-13; **confirmed by Codex 2026-07-13** ("no later gate/default found"), with the
  refinement that the null case is **air-specific**: Buffer domains require a material in
  `Domain::Domain()` (`cl_FEM_Domain.cpp:36-45`). Confidence: high (code path); medium that
  no run has hit it yet (helper rewire flagged "pending build" in
  `tmp/ai_exchange/maxwelldata_helper.md`).
- [◐] **D2 — FIXED in the helper 2026-07-13 (Christian, verified by Claude in tree same
  day):** the MaxwellData ctor now binds
  `aMaxwellKernel->dofmgr()->block( group()->id() )->material()` (`cl_FEM_Calculator.cpp:67`)
  — always the Maxwell-side material, matching the legacy convention. **Residual (legacy,
  outside the helper):** the thermal Conductor dispatch still reads the thermal group's
  (null) material (`cl_IWG_MaxwellThermal.cpp:70-74`) — unexercised today, disappears at R12
  when the dispatch collapses; until then a volume-conductor thermal model still crashes
  there. Original finding, preserved for audit history:
  Thermal-side MaxwellData must use the **Maxwell** block material, not the thermal group.
  `Controller::set_thermal_kernel` relinks the **thermal** kernel's block calculators
  (`cl_FEM_Controller.cpp:1827-1833`), but `ThermalFactory` never assigns materials (zero
  references in `cl_ThermalFactory.cpp`) and `auto_set_materials` only copies within one
  kernel (`cl_FEM_DofManager.cpp:1053-1073`, `mIndex > 0`). Legacy thermal kernels read the
  **Maxwell** calculator's material (`mt_thermal_h.cpp:32`). Fix: bind
  `mMaxwellCalculator->group()->material()` in the ctor (after D1's gate).

  Found by Claude 2026-07-13; **confirmed by Codex 2026-07-13** with two additions: the thermal
  field is dof manager index 0 (`cl_FEM_Kernel.cpp:587-608`) so `auto_set_materials` (mIndex>0)
  can never help it, and — a **pre-existing suspect independent of this refactor** — the
  thermal Conductor dispatch already reads `aGroup->material()` to choose `T_h_*`
  (`cl_IWG_MaxwellThermal.cpp:70-74`) while the ThinShell dispatch correctly uses the Maxwell
  material (`:119-129`); verify in R1 how (or whether) the Conductor path currently works.
  Confidence: high.

  **Clarification (2026-07-13, answering Christian's "the kernel is getting materials from
  somewhere"):** it gets them from the **Maxwell side**, at both levels — every `T_h_*` body
  reads `tMaxwellCalculator->group()->material()` (`mt_thermal_h.cpp:32,360,...`), and the
  ThinShell dispatch does the same (`cl_IWG_MaxwellThermal.cpp:121-128`). The thermal groups
  themselves carry **no** material (verified: `ThermalFactory::create_thermal_kernel` read in
  full, `cl_ThermalFactory.cpp:59-252`, no material assignment; `Block` ctor none). The
  Conductor dispatch branch (`:70-74`) that reads the thermal group's material is
  **unexercised** — current thermal examples are thin-shell (`Tape_Quench` uses
  `thinshell : tape`, `tmp/examples/Tape_Quench/BuiltinMat/input.conf:123`); a
  volume-conductor thermal model would null-deref there **today**, independent of this
  refactor. Both fixes align on the same convention: material always from
  `mMaxwellCalculator->group()`.
- [x] **D3 — RESOLVED by policy 2026-07-13 (decided by Christian).** The two β definitions
  are per-material-family and never combined in one material: (b,n) for HTS/REBCO (only HTS
  materials have or need jc; their lambda is T-only), (b,j) for metals. The slot aliasing is
  therefore unreachable by policy. ~~Fix: split into two slots (betaBJ / betaBN)~~ → reduced
  action (code → **R2a**): a comment at the β slot documenting which convention each family
  uses, plus a cheap ctor guard making the policy loud —
  `BELFEM_ERROR( !( tIsHTS && tLambdaFieldDependent ), "material %s combines jc with
  field-dependent lambda — beta conventions would alias" )` — so a future UserDefined
  material that violates the policy is rejected at setup instead of silently mixing angles.

  **Original finding, preserved for audit history; superseded action noted above:**
  `compute_lambda_metal`/`compute_rho_metal` cache **bj**_angle under slot
  `MaxwellDataValue::beta` (`hpp:2434-2447,2569-2582`); `compute_rho_*_ts` cache
  **bn**_angle under the *same* slot (`hpp:2628-2645,2674-2707,2746-2792,2822-2869`);
  bulk-HTS relies on mBeta staying at the reset π/2 (`hpp:2304,2590`). For a material with
  `have(jc) && depends(lambda,normB)` (constructible via UserDefined,
  `cl_Material_UserDefined.cpp:101-103`), whichever runs first at index k poisons the other —
  legacy kept the two betas in separate locals (`mt_thermal_h.cpp:78` vs `:102`).
  Fix: split into two slots (betaBJ / betaBN) or a dedicated `mBetaLambda`. Found by Claude
  2026-07-13; **confirmed mechanically by Codex 2026-07-13**, reachability medium: no built-in
  combines the two (YBCO lambda is T-only, `cl_Material_YBCO.cpp:107-108`; Cu/Ag have
  field-dependent lambda but no jc, `cl_Material_Copper.cpp:602-608`,
  `cl_Material_Silver.cpp:487-493`), but a legal UserDefined material can
  (`cl_Material_UserDefined.cpp:93-116,128-166`). Fix anyway — the aliasing is
  order-dependent and silent. Confidence: high (code), medium (reachability).
- [x] **D4 — CLOSED 2026-07-19:** the R2c implementation is fully in tree and
  Codex-audited — `Material::jc_eval()`/`n_eval()` dependency-routing helpers
  (`cl_Material.hpp:301-304`) used at all four full-overload sites in `powerlaws.hpp`
  (`:98,108,221-241,646-714`); policy per O1 (full-signature calls, routed internally on
  `depends()`). History below (original finding kept for the overload map):
  **UserDefined reduced-dependency overload parity broken** (= trap T14).
  **Policy DECIDED 2026-07-13 (Christian): full-signature calls, routed internally on
  `Material::depends()`/`dependencies()` — see O1 (RESOLVED) for the decision record and
  the three implementation items (assert relaxation, per-JcFunction eval check, doc
  paragraph) → **R2c/R2g**.** Original finding and analysis:
  Fix inside the helper ctor: extend the HTS branch to key the rho/drhodj pointers on
  `type()==UserDefined` × `depends(jc, T/normB/angleNxB)` exactly mirroring
  `mt_maxwell_h.cpp:1999-2033` (see O1). Found by Claude 2026-07-13 (brief flagged
  "confirm parity"; verified it does NOT hold). **Codex 2026-07-13 confirmed the piecewise
  family too:** `rho_piecewise(nJ,nB,β,T)` asserts full deps and uses rho(T)
  (`powerlaws.hpp:636-655`, defect twin `:717-734`) while the reduced B/angle overloads use
  rho(gTbulk) (`:484-495,559-570`); the `drho_piecewise_dJ` full overloads skip the
  dependency assert but still use rho(T) (`:1750-1764,1825-1839`) vs rho(gTbulk) reduced
  (`:1600-1614,1675-1689`). Confidence: high.
  **Policy direction 2026-07-13 (Christian):** the code is unpublished and the material API
  is not stable — exact bug-for-bug parity with the legacy generic branches is NOT the goal;
  the requirement is (i) a clear, documented jc-dependency policy in the material API and
  (ii) all provided examples keep working (they are the acceptance test). O1 reframed
  accordingly; the overload matrix above stays as the map of what any chosen policy must
  cover. *(Argument order in this entry predates the O10 reorder — see the §3 reconciliation
  note.)*
- [x] **D5 — CLOSED 2026-07-13** with the D8/D9 fixes (Christian) plus the final tail
  (Claude, source edit approved by Christian): dispatchers folded to a single
  null-or-stale condition with the peer cached in a local, and the legacy intpoint-match
  assert restored in both coupled dispatchers
  (`cl_FEM_Calculator.cpp:2382-2418`). History below.
  **Implemented 2026-07-13 (Christian), differently from the plan's `link_peer()` proposal
  and better:** `Calculator::link(Element*)` now tail-dispatches through `mFunLinkElement`
  (`cl_FEM_Calculator.cpp:1246`); Maxwell-data blocks get `link_element_maxwell()`
  (`cl_FEM_Calculator.cpp:2336-2360`), which resets the helper, sets `mElement`, and
  cross-links the peer calculator with an element-ID guard. **Claude double-checked the
  recursion 2026-07-13 (Christian's ask): NO infinite loop** — `mElement = aElement` is set
  *before* the cross-link, and the peer's re-entrant `link_element_maxwell` compares element
  IDs against the already-updated first calculator, so recursion terminates at depth 2.
  Follow-up defects D8/D9 (found in the first version) were fixed by Christian's
  three-way dispatcher split same day; **D11** (rebuild-path selection) remains, and the
  legacy intpoint-match assert (`mt_maxwell_h.cpp:2487-2489`) was not carried over — add a
  `BELFEM_ASSERT( num_intpoints() == peer->num_intpoints(), … )` after the peer link.
  Original clarification, preserved for the record:
  **Clarification (2026-07-13, answering Christian's "the calculator is relinked to every
  element — I don't understand the concern"):** the **assembling** kernel's own calculator is
  relinked per element (`compute_mkf` → `link()`, `cl_IWG_Maxwell.cpp:228`). The missing link
  is the **other kernel's** calculator. When the thermal IWG assembles element X, the
  *Maxwell* block calculator — whose `E(k)/C(k)/q()` the helper reads for b/j — still points
  at whatever element it linked last (a different element, or last timestep's).

  Legacy kernels bridge that gap explicitly every element: `T_h` calls
  `tCalculator->link( tElement )` on the Maxwell calculator (`mt_thermal_h.cpp:34`), and the
  `_t` Maxwell kernels mirror it through `get_thermal_calculator()`
  (`mt_maxwell_h.cpp:2474-2492`). After delegation to the helper, no code performs this
  cross-link — `Calculator::link()` resets only its own MaxwellData
  (`cl_FEM_Calculator.cpp:1244`). Without it, `compute_j`/`compute_T_fem` silently evaluate
  the wrong element. Fix: `MaxwellData::link_peer()` (§6): thermal-side links
  `mMaxwellCalculator` to `dofmgr->element(id)`, Maxwell-side `_t` links
  `mThermalCalculator`; both re-assert matching intpoint counts (contract from
  `mt_maxwell_h.cpp:2487-2489`). Kernels call it once per element before the k-loop.
  Confidence: high.
- [x] **D6 — RESOLVED 2026-07-13 (decided by Christian): keep the self-assignment.**
  `mB = (this->*mFunB)(aIndex)` (`hpp:2319-2328`) is deliberately in-your-face (readability
  for the tired 2 AM human beats micro-tidiness). Action reduced to a one-line comment at the
  site noting it is a knowing self-assign. Code → **R2b**.
- [x] **D7 — RESOLVED 2026-07-13 (Christian): non-issue.** All elements of a block have the
  same node count, so `mCoords = N*X` (`hpp:2205-2220`) sizes once per block and never
  reallocates afterwards. (Still strictly better than legacy, which constructed a fresh
  `Matrix` per integration point, `mt_maxwell_h.cpp:589`.) No action.
- [x] **D8 — FIXED 2026-07-13 (Christian, verified by Claude in tree same day):** the
  dispatcher is now split three ways — `link_element_maxwell` (no thermal kernel),
  `link_element_maxwell_thermal`, `link_element_thermal_maxwell` — selected in `allocate()`
  from `mThermalKernel == nullptr` / `mGroup->parent()->iwg()->type()` with a loud
  `BELFEM_ERROR` fallback (`cl_FEM_Calculator.cpp:1062-1086`). `iwg()` is valid at that
  point (the Calculator ctor itself already reads `parent()->iwg()->delta_time()`, `:251`).
  This kills the dead `mIwgType` member path entirely and saves the per-element branch.
  **But see D11:** the selection is made only in `allocate()` — the `link_maxwell` rebuild
  path does not re-select. Original finding, preserved for audit history:
  **`mIwgType` is never set on block calculators, so
  `link_element_maxwell` branch selection is dead.** `mIwgType` is assigned only in
  `Calculator::link( Group * )` (`cl_FEM_Calculator.cpp:1108`), and the **only** caller of
  that overload is the SideSet path (`cl_FEM_SideSet.cpp:88`); blocks never call it and
  `IWG::link_to_group` doesn't either (`cl_IWG.cpp:383-…`, verified: sets `mCalc` but never
  links it to the group). Every block calculator therefore carries
  `mIwgType == IwgType::UNDEFINED` and takes the **else** branch of `link_element_maxwell`
  (`cl_FEM_Calculator.cpp:2351-2358`). Consequences: thermal-side calcs work *by accident*
  (else is their correct branch); **Maxwell-side calcs never link the thermal peer** —
  for them `mMaxwellData->maxwell()` is *itself*, whose `mElement` was set two lines above,
  so the ID guard always compares equal and the branch no-ops. `compute_T_fem` then reads a
  stale thermal element during every Maxwell `_t` assembly → **silently wrong ρ(T), no
  crash** — exactly the failure class this plan guards against. Fix (recommended): drop the
  IwgType test entirely and bind the peer once at construction —
  `mPeer = ( this == mMaxwellData->maxwell() ) ? mMaxwellData->thermal()
  : mMaxwellData->maxwell()` — identity-based, immune to type wiring; alternatively set
  `mIwgType` in `allocate()` from `mGroup->parent()->iwg()->type()`. Found by Claude
  2026-07-13 while double-checking the loop question. Confidence: high (all three code
  points read this session).
- [x] **D9 — FIXED 2026-07-13 (Christian, verified by Claude in tree same day):** both
  coupled dispatchers now check `peer->element() == nullptr` first and link+return
  (`cl_FEM_Calculator.cpp:2361-2370,2380-2389`). The extra per-**element** branch is
  well inside budget — it runs once per element, not per integration point, and is dwarfed
  by the coordinate/edge-function work `link()` already does. Style option (not required):
  fold the
  two ifs into one — `if ( tPeer->element() == nullptr || tPeer->element()->id() !=
  aElement->id() )` — with the peer cached in a local `Calculator * tPeer`, which also drops
  the repeated `mMaxwellData->thermal()` chases. Original finding: a never-linked peer has
  `mElement == nullptr` (init `cl_FEM_Calculator.hpp:353`; `element()` unchecked,
  `:2048-2051`) → null-deref on the first coupled element. Found by Claude 2026-07-13.
- [x] **D10 — FIXED 2026-07-13 (Claude, source edit approved by Christian):** dead
  `reset_data()` declaration removed; replaced by the documented
  `select_link_element_dispatcher()` declaration (`cl_FEM_Calculator.hpp:1020-1024`).
- [x] **D11 — FIXED 2026-07-13 (Claude, source edit approved by Christian):** the
  three-way selection now lives in `Calculator::select_link_element_dispatcher()`
  (`cl_FEM_Calculator.cpp:1076-1095`), called from `allocate()` (`:1066`) **and** from
  `link_maxwell`'s already-allocated rebuild branch (`:1126`) — so a thermal-side
  calculator that allocated before `set_thermal_kernel` gets its dispatcher upgraded when
  the relink fires. Original finding, preserved for audit history:
  **`link_maxwell`'s rebuild path does not (re)select `mFunLinkElement`.** The dispatcher choice lives only in `allocate()`
  (`cl_FEM_Calculator.cpp:1062-1086`); the `link_maxwell` already-allocated branch
  (`:1108-1115`) rebuilds `mMaxwellData` but leaves `mFunLinkElement` at whatever
  `allocate()` chose. The thermal kernel initializes **eagerly inside its factory**
  (`mThermalField->initialize()`, `cl_ThermalFactory.cpp:242`) — *before*
  `Controller::set_thermal_kernel` runs (`hphiTrun.cpp:83`) — so its block calculators
  allocate with `mMaxwellKernel == nullptr` → `mFunLinkElement = link_element_default`.
  The later relink builds their MaxwellData but per-element linking stays `default`:
  **`MaxwellData::reset()` never fires between elements**, so k-indexed memoized values
  (`mLastIndex` matches) from element A are silently served for element B — wrong rho/T/b/j
  on the whole thermal side, no crash. (The magnetic kernel escapes only because its
  `initialize()` is deferred to Controller runtime, after `link_maxwell`.) Fix: extract the
  three-way selection into one private helper (e.g. `select_link_element_dispatcher()`)
  called from both `allocate()` and `link_maxwell`'s rebuild branch. Found by Claude
  2026-07-13 while verifying D8/D9. Confidence: high (ordering verified:
  `cl_ThermalFactory.cpp:242` vs `hphiTrun.cpp:83`; `mIsAllocated` set at `allocate()`
  entry, `cl_FEM_Calculator.cpp:634`).
- [x] **D12 (CRITICAL) — MaxwellData construction + dispatcher selection were DEAD CODE for
  blocks: placed below `allocate()`'s "done if this is a block" early return.**
  Caught by Christian's first R1 smoke run (SIGSEGV, jump to 0x0 from
  `Calculator::link` — `mFunLinkElement` null on the first assembled block element).
  `allocate()` returns early for every BLOCK group (`cl_FEM_Calculator.cpp:833,940`,
  comment "done if this is a block"), so the tail section that built MaxwellData and
  selected `mFunLinkElement` never executed for blocks — the entire allocate-path wiring
  (recommended in the original F1 audit as "construct at the END of allocate()") was
  unreachable for exactly the groups that need it; nothing noticed until `link(Element*)`
  started dispatching through the pointer. **FIXED 2026-07-14 (Claude):** the section is
  hoisted to immediately after `mIsAllocated = true` (`cl_FEM_Calculator.cpp:661-677`) —
  work vectors are no precondition since `link_vector()` creates missing ones on demand
  (F5 design); belt-and-braces: `mFunLinkElement` now has a default initializer
  (`= & Calculator::link_element_default`, hpp declaration) so no allocate() path can ever
  leave it null again. Found by smoke test 2026-07-14; diagnosed + fixed same day.
  Confidence: high (backtrace matches; both early returns audited — the `:653` one exits
  before `mIsAllocated` and is safe).
- [x] **D13 (CRITICAL) — ThinShell groups never build a MaxwellData: `h_calc` on the TS·LookupAlloy branch fires the `maxwell()` assert.**
  **Fixed / dissolved — confirmed 2026-08-09.** Christian's `DomainType::ThinShell`
  addition to the magnetic block selection (`cl_MaxwellFactory.cpp`) made conductive
  layers real FEM Blocks on the magnetic side too, and `Calculator::allocate` now
  constructs `MaxwellData` for every non-Air block that exists in the peer dof manager
  (`cl_FEM_Calculator.cpp:1111-1137`). The "sideset MaxwellData design" this defect was
  waiting for is therefore not needed, and `h_picard`/`h_newton_*` serve Conductor and
  ThinShell from the same branch (`cl_IWG_Maxwell.cpp`). `debt_register.md` DR-04 closed
  with this box.
  **Discovery.** Found by Claude 2026-07-19 during the R6 status assessment.

  **Wrong prior attribution.** This is the actual cause of the 2026-07-16
  double-corc failure recorded in `devlog/dl20260716_periodic_seam_phase_two.md`,
  not D1. The devlog's "D1, null-material air blocks" attribution is wrong: the
  D1 gates are in tree and functional.

  **Mechanism.** (Precise two-level statement, corrected 2026-07-20 after
  Christian's objection to "thin shells are FEM sidesets".) At the **mesh**
  level thin shells ARE blocks: `ThinShellFactory` creates one block per layer
  and stamps `DomainType::ThinShell` on each (`cl_ThinShellFactory.cpp:257-271`).
  The **thermal** kernel selects those layer blocks, so there they become FEM
  Block groups (`cl_ThermalFactory.cpp:74-136`). The **magnetic** kernel does
  not: its block selection admits only Air/Buffer/Conductor/Ferro — the
  ThinShell layer blocks fall through the `default:` and are skipped
  (`cl_MaxwellFactory.cpp:513-532`). Instead the factory stamps the tape
  *sidesets* with `DomainType::ThinShell` (`cl_MaxwellFactory.cpp:2571-2577`),
  and those become FEM SideSet groups (`cl_FEM_DofMgr_SideSetData.cpp:95-170`).
  `IWG_Maxwell::link_to_group` has a single `domain_type()` switch
  (`cl_IWG_Maxwell.cpp:255`), so on the magnetic dofmgr its ThinShell case fires
  on the sideset group's calculator. MaxwellData construction is gated on
  `GroupType::BLOCK` in both `allocate()` (`cl_FEM_Calculator.cpp:675-682`) and
  `link_maxwell()` (`:1146-1151`), and the factory only calls `link_maxwell` on
  `dofmgr()->blocks()` (`cl_MaxwellFactory.cpp:736`). The R6 TS-alloy flip
  therefore dispatched `maxwell::h_calc` onto a calculator whose `maxwell()`
  accessor asserts (`cl_FEM_Calculator.hpp:3023`).

  **Why corc hits it.** corc reaches the branch because hastelloy is a
  LookupAlloy tape layer (`examples/corc/input.conf`).

  **Thermal side — believed unaffected, with an OPEN verification item (added
  2026-07-20).** Thermal shell layers are FEM blocks on the thermal dofmgr, and
  `Controller::set_thermal_kernel` relinks both kernels' blocks
  (`cl_FEM_Controller.cpp:1820-1833`), so the thermal `T_h_calc` flips do build
  a MaxwellData. BUT: the MaxwellData ctor binds its Maxwell-side calculator and
  material via `aMaxwellKernel->dofmgr()->block( aCalculator->group()->id() )`
  (`cl_FEM_Calculator.cpp:65-67`), and `DofManager::block(id)` returns the
  shared **empty block** (null material) for ids the magnetic dofmgr does not
  carry (`cl_FEM_DofMgr_BlockData.hpp:130-142`). If shell-layer block ids really
  have no magnetic FEM block, that bind is null — same pattern as the legacy
  thermal ThinShell dispatch read (`cl_IWG_MaxwellThermal.cpp:123-129`), which
  predates R6 and presumably worked in coupled runs, so either a magnetic-side
  block under those ids exists after all (mechanism not yet found) or the
  coupled TS path was never exercised. Must be traced before trusting the
  thermal TS `T_h_calc` flip. Confidence on "thermal unaffected": downgraded to
  medium.

  **Mitigated 2026-07-19 (Claude, approved).** The TS·LookupAlloy dispatch
  (`cl_IWG_Maxwell.cpp:450-465`) reverted to legacy `h_alloy`/`h_alloy_t`
  routing. Bulk-Conductor `h_calc` and both thermal `T_h_calc` flips stay.

  **Residual design work (R9 prerequisite).** MaxwellData construction + material
  binding for thin-shell sidesets still need design work. The ctor's
  `dofmgr()->block( group()->id() )` material lookup has no meaning for a sideset
  id, so R9/R10 TS flips need a designed sideset path: extend the gates + add a
  TS-aware material bind. This cannot be fixed by just removing the gate.

  **Confidence.** high on mechanism; medium-high that the TS branch was corc's
  exact crash site (a backtrace would settle it).

  **2026-07-20 follow-up (Christian's block-selection change + deeper trace).**
  - Christian added `DomainType::ThinShell` to the magnetic block selection
    (`cl_MaxwellFactory.cpp:522`), so conductive layer blocks now become FEM
    Blocks. Facts verified around it: `ThinShellFactory` stamps ThinShell on
    every layer block (`cl_ThinShellFactory.cpp:264`) but `create_buffers`
    demotes rho-less layers to Buffer (`:1730-1736`) — corc's magnesia layer was
    therefore ALREADY selected before the change (Buffer was in the switch).
    The dof tables already anticipate ThinShell blocks: `collect_block_dofs`
    routes ThinShell → Conductor dof table (`cl_Maxwell_FieldList.cpp:256-263`).
    Per-block materials come from `mMaterialBlockAssignment`, which includes
    layer blocks (`cl_MaxwellFactory.cpp:2361-2385`). Both `link_maxwell` loops
    now reach the layer blocks, so per-layer MaxwellData construction works
    without touching the BLOCK gates, and the empty-block hazard on the thermal
    side closes (ids resolve to real magnetic blocks).
  - Open consequences of the change: `block_activation_mode(ThinShell)` is not
    in the activation map, so it defaults to FULL `GeometryAndDofs`
    (`cl_IWG.cpp:2135-2146`) — layer blocks will be assembled by the DofManager
    block loops, and `link_to_group`'s ThinShell case currently dispatches the
    legacy facet-based TS kernels for them (wrong for volume blocks). Needs
    either a `mBlockActivationModes[ThinShell] = GeometryOnly` line in
    `IWG_Maxwell::init_activation_maps` (`cl_IWG_Maxwell.cpp:109-118`, blocks as
    MaxwellData carriers only) or a GroupType split in the ThinShell dispatch
    case (blocks → collapsed bulk kernels). Christian's intent decides.
  - VESTIGE FOUND: `Group::number_of_thin_shell_layers()` has NO override
    anywhere in the current tree (the historical `Mesh` accessor is gone); the
    base always `BELFEM_ERROR`s, yet `IWG::link_to_group`'s BlockSpecific
    SIDESET branch calls it for ThinShell sidesets (`cl_IWG.cpp:459-466`). A TS
    sideset reaching that branch dies with the "sidesets or shells" message —
    NOT the h_calc assert. Since the 2026-07-16 corc run died at the h_calc
    assert, the TS *sideset* branch was likely never the crash site; confidence
    on "TS branch = corc crash site" drops to LOW. How thin-shell assembly is
    actually driven today (and which group fired `h_calc` on 2026-07-16) needs
    Christian's knowledge or a backtrace.

  **2026-07-20 Codex audit (thread `tmp/ai_exchange/thinshell_assembly_tree.md`;
  citations spot-checked by Claude; confidence high except Q3 medium).**
  - **Assembly driver:** with Christian's selection edit, the layer BLOCKS
    assemble in the normal transient block loop
    (`cl_FEM_DofManager.cpp:633-665`) and are the thin-shell physics carrier —
    the `h_ts_*` kernels are **volume-style** (integration-point loop with
    `compute_bn`, `mt_maxwell_h.cpp:89-151`), NOT facet kernels as Claude
    assumed. Ghost sidesets (DomainType::Ghost, active by default) carry the
    facet-based `h_ghost` interface kernel (`mt_maxwell_h.cpp:1752-1787`). Tape
    ThinShell sidesets are GeometryOnly and skipped by the transient guard
    before `link_to_group` (`cl_FEM_DofManager.cpp:690-696`) — no duplicate
    physics from the edit. Corollary: BEFORE the edit, neither blocks nor tape
    sidesets assembled the TS conductor physics — the edit completes the
    migration started when the sideset-layered TS machinery was dismantled
    (override removed at `e5f533ae`; it existed at `c54027b3`).
  - **Vestige verdict:** the `number_of_thin_shell_layers()` call is dead in
    the transient path (guard above). Latent trap: the non-transient
    `compute_jacobian()` sideset loop has NO `is_active()` guard
    (`cl_FEM_DofManager.cpp:496-510`) — would die if ever used on a TS model.
    Cleanup candidate for R12.
  - **Crash attribution REFUTED:** at `abb71c2e`/`ce6a0e8b` the TS·alloy
    `h_calc` dispatch was unreachable (blocks unselected, sidesets skipped),
    and `examples/corc/input.conf` has NO bulk Conductor block (air blocks 1:2
    only; hastelloy = tape layer + built-in LookupAlloy,
    `cl_Material_HastelloyC276.cpp:41-42`). The 2026-07-16 double-corc `h_calc`
    assert is NOT explained by the pristine example deck + normal transient
    path — the campaign deck may differ (bulk LookupAlloy block?) or a
    MaxwellData gate was bypassed on that run path. Needs the actual deck or a
    backtrace (Christian).
  - **Revert parity note:** `h_alloy` (the D13 revert target) is bulk-style
    without `compute_bn`, unlike `h_ts_metal` — that is the pre-existing T15
    "TS-alloy-uses-bulk-kernel" property, i.e. exact legacy parity, not a
    regression. R9 owns bn-handling when the collapsed kernel takes over TS.
  - **Updated fix direction:** with layer blocks selected and carrying
    MaxwellData, the collapsed kernels become viable for TS layers via the
    BLOCK route — the sideset-MaxwellData design this entry originally called
    for may be unnecessary. Pending Christian: confirm the
    layer-blocks-as-volume-carrier design (then GeometryAndDofs default is
    intended, layer edges gain conductor h-dofs — a deliberate dof-structure
    change) and supply the crash backtrace.

- [x] **D14 — `material::Alloy` field overrides unreachable through dependency-based
  dispatch — FIXED 2026-07-21 (Claude, approved by Christian; found by Codex during the
  thermal-collapse audit).** `Alloy` is `SplineLookupTable( MaterialType::PureMetal )`
  (`cl_Material_Alloy.cpp:37,49-53`) with genuine field-aware overrides — 3D
  ρ(T, log₁₀B, β) database lookup and WF-quotient λ(T,B,β)
  (`cl_Material_Alloy.hpp:141-184`) — but its splines register T-only dependency
  (`cl_Material_SplineLookupTable.cpp:47`) and no normB/angleBxJ bits were ever set
  (only Copper/Silver set them, e.g. `cl_Material_Copper.cpp:602-608`). Since MaxwellData
  keys metal routing on `depends(rho/lambda, normB)` (`cl_FEM_Calculator.cpp:94,100-101`),
  Alloy silently fell to the T-only `compute_rho_bulk`/`compute_lambda_bulk` — on BOTH
  kernels once the collapses landed, while legacy dispatch (keyed on `MaterialType`)
  had called the field-aware overloads via `T_h_metal`/`h_metal`. Fix: set the
  rho/lambda T+normB+angleBxJ dependency bits at the end of
  `Alloy::populate_rho_database()` (`cl_Material_Alloy.cpp:484-499`) — the single funnel
  through which the field database becomes valid (compute and file-load paths), so
  decks that never build the database keep T-only routing and never hit the
  `mRhoData` assert. Not corc-relevant: hastelloy is `Metal(…, LookupAlloy)`
  (`cl_Material_HastelloyC276.cpp:41`), legacy-thermally T-only `T_h_alloy` anyway.
  Confidence: high (all sites read this session).

---

## 4. Ordered Steps

Sequencing rule: helper first (it must be provably correct before any kernel delegates),
then one dispatch family per step, each independently buildable and shadow-verified.
**Christian runs all builds/examples** (standing rule); AI prepares diffs and analysis.

> **Checkbox semantics (added 2026-07-14 after this bit Christian):** boxes in §3.2 and §5
> track the *finding or decision* — `[x]` there means "we know exactly what to do (or the
> fix is in)". Boxes here in §4 track *code + verification*. So a ticked D/O upstream with
> an open R-item here simply means **decided, not yet implemented**. To keep the two
> lifecycles apart, R2's sub-items carry their own IDs (R2a…R2g) and name the decision they
> implement.

- [x] **R1 — Build/smoke gate for completed wiring. PASSED 2026-07-14 (Christian).**
  Wiring code in tree (2026-07-13/14): D1 Air gates, D2 ctor material rebind, D5/D8/D9
  three-way `link_element_*` dispatchers (Christian); D10/D11 shared selector + intpoint
  asserts, R2 batch (Claude, approved). The gate earned its keep: the first smoke attempt
  SIGSEGV'd on hphirun → **D12** (MaxwellData wiring was dead code for blocks, below
  allocate()'s block early-return — a hole three read-only audit passes missed); fixed by
  hoisting + pointer default, rebuild, **smoke passes**. Any §4.1 case not yet in the
  smoke set gets exercised again at R4 baselines and every family flip.
- [ ] **R2 — Helper hardening to equivalence-grade.** Every design decision is already
  made (§3.2/§5); this step is pure implementation of them. Boxes = code landed:
  - [x] ~~wiring items (D1 residual gate, D8, D9, D10, D11, intpoint asserts)~~ — landed
    2026-07-13, details in §3.2; kept so R2's original scope stays visible.
  - [x] ~~D5 `link_peer()`~~ — obsolete: superseded by the `link_element_*` dispatchers
    (D5 closed in §3.2); nothing left to implement.
  - [x] **R2a** (implements D3): β-convention comment + ctor guard
    `BELFEM_ERROR( !( tIsHTS && tLambdaFieldDependent ), … )`; `beta_dummy()` accessor with
    the π/2 assert at all 8 bulk-HTS sites. *(Claude 2026-07-14, audit pending)*
  - [x] **R2b** (implements D6): "deliberate self-assignment" comment at `compute_b`.
  - [x] **R2c** (implements D4/O1): `Material::jc_eval()`/`n_eval()` dependency-routing
    helpers; all 8 T-bearing powerlaw/piecewise rho+drho overloads delegate; full-dep
    asserts removed (reduced-overload asserts kept). Per-subclass eval verification PASSED
    (ModifiedKim 2-arg, Database 3-arg, UserDefined both — devlog table).
  - [x] **R2d** (implements O2): clamp in `compute_T_fem` into `[gTmin, material T_max]` +
    `mTClamped` + zero dT-derivatives while clamped + `T_clamped()` accessor — DONE;
    `compute_T_const` clamps + flags identically (added after Codex's audit caveat);
    ~~**deferred:** the once-per-step converged-at-clamp diagnostic (needs Controller-side
    convergence context; accessors ready)~~ — **CLOSED as won't-do (2026-08-14, Christian's
    ruling against Claude's proposal, on the evidence of the tapestack3d printout):** the
    runaway path into clamping goes through degraded convergence, which
    `print_thermal_stall_warning` already voices; a genuinely converged-at-ceiling state
    shows as a `T_max` plateau in the exodus T field, which quench analysis always
    inspects; and the stabilized Controller is not worth perturbing for a diagnostic whose
    triggering scenario has two existing observability channels. The `T_clamped()` /
    `rho_clamped()` accessors stay as-is (zero consumers, trivially small) in case a
    future campaign wants them.
  - [x] **R2e** (implements O3): clamp inside `compute_rho` + `mRhoClamped`;
    `compute_drhodj` forces rho current and returns 0 while clamped. **Includes a latent-bug
    fix found en route:** `gRhoMin`/`gRhoMax` had no defaults in
    `Communicator::set_globals()` (zero-initialized → the new clamp would have pinned every
    rho to [0,0] outside hphirun/hphiTrun); defaults now `0 / BELFEM_REAL_MAX`.
  - [x] **R2f** (implements O4): `mDensity` ctor-cached (`density(gTroom)` →
    `ref_density()` → NaN + asserting getter) + physics-trap comment.
  - [x] **R2g** (docs for O1+O4): "Assembly Contract: the Full-Signature Policy" and
    "Density and the Undeformed Mesh (Physics Trap!)" added to
    `materials_usage_guide.md` (+ revision row).
  - [x] **R2h** (T10 audit): table in `devlog/dl20260714_maxwell_helper_r2_hardening.md` —
    built-in HTS/alloy/Cu/Ag rows match; TWO ⚠ rows for the R3 shadow harness to
    adjudicate: non-Kohler PureMetals (helper follows deps → rho(T); legacy type-dispatch →
    rho(T,B,β) via the Kohler pointer) and user Alloys typed `MaterialType::PureMetal`
    (`cl_Material_Alloy.cpp:38,53` — intentional? question for Christian).
  *(implemented 2026-07-14 by Claude/Fable, Christian's go-ahead. **Audit round done same
  day:** Codex — "no blocking R2a-R2h defect found", high confidence; its two caveats were
  fixed on the spot (`compute_T_const` now clamps + flags like the FEM path, and the two
  reduced derivative overloads gained the missing symmetric `!depends_on(T)` asserts —
  a pre-existing asymmetry, not an R2 regression). Grok third voice UNAVAILABLE this round
  (two truncated runs; its four audit questions answered by direct code check instead —
  see the exchange). R2 closes when the R2d converged-at-clamp diagnostic finds its
  Controller-side home and Christian's rebuild+smoke passes.)*
- [ ] ~~**R3 — Shadow-compare harness.**~~ **CUT 2026-07-14 (decided by Christian):** not
  worth the effort — the tree is committed before each flip, so "revert on breakage"
  replaces matrix-level A/B comparison; verification is run-based (see §4.1, revised).
  Original design kept below for the record in case a flip ever misbehaves subtly enough
  to want it back: Debug-only compile switch (e.g.
  `BELFEM_SHADOW_MAXWELL_KERNELS`) in `IWG_Maxwell::compute_mkf` /
  `IWG_MaxwellThermal::compute_mkf` (`cl_IWG_Maxwell.cpp:223-234`,
  `cl_IWG_MaxwellThermal.cpp:37-49`): run the NEW kernel, snapshot M/K/f/dKdx, reset, run the
  LEGACY kernel, assert `max|Δ| ≤ 1e-12·‖·‖∞` per matrix (0 exactly where §3.1 predicts
  bit-equality), leaving the *legacy* results in place (legacy-last so shared work-vector
  state — including alloy's stale b/j — matches the old code for any downstream reader;
  **ordering confirmed by Codex 2026-07-13**). Between the two passes, explicitly reset BOTH
  `TimestepMatrices` (its `reset()` only zero-fills storage, `cl_TimestepMatrices.cpp:63-107`)
  AND the helper caches (`aCalc->maxwell()->reset()`) — `Calculator::link()` resets the helper
  only once, before the harness (`cl_FEM_Calculator.cpp:1244`), so the new kernel's memoized
  β/norm/T would otherwise leak into the legacy pass (Codex refinement). This is the
  per-variant equivalence proof a reviewer replays. *(after: R2)*
- [ ] ~~**R4 — Freeze baselines.**~~ **CUT 2026-07-14 (decided by Christian):** the git
  commit IS the baseline — any reference output can be regenerated on demand by checking
  out the pre-flip commit and rerunning. Original step kept for the record: With the pre-refactor tree (R1/R2 helper fixes included but
  kernels untouched), archive for each regression case (§4.1): exodus/HDF5 output at fixed
  step counts, per-timestep nonlinear iteration counts, and the `element_rho` field.
  *(after: R2; independent of R3.* Note: strictly only R1 is required — R2 touches only the
  helper, which legacy kernels never call, and its powerlaws changes are assert-only — but
  freezing after R2 pins a single tree state for the whole campaign; pull earlier only if
  R2 stalls.)*
- [x] **R5 — Write the collapsed kernel pair. DONE — build confirmed by Christian
  2026-07-14 ("looks really nice").** Code in tree (Claude):
  `maxwell::h_calc` / `maxwell::h_newton_calc` (`mt_maxwell_h.{hpp,cpp}`) and
  `fem::T_h_calc` (`mt_thermal_h.{hpp,cpp}`), declared+documented beside the legacy
  functions, zero call sites (dispatch 100% legacy). Bodies mirror the legacy expression
  forms literally (`w(k)·rho·dV(k)` groupings, per-term mu0, `norm_j > BELFEM_EPSILON`
  Newton guard, `rho·norm_j·norm_j` source). Remaining: Christian's build. (new functions beside the old:
  `maxwell::h_calc` / `maxwell::h_newton_calc` and `maxwell::T_h_calc` /
  `maxwell::T_h_newton_calc` working names; final naming takes over `h`/`h_newton` and
  `T_h`/`T_h_newton` at R11/R12). Bodies per §6, mirroring the legacy expression forms
  literally (FP-order preserved — doubly important with R3 cut). `T_h_newton_calc` is NOT
  part of R5 (new physics, owned by `thermal_matrices_cleanup_and_newton_plan.md` R4-R5).
  Builds with dispatch still routing 100% legacy. *(after: R2)*
- [x] **R6 — Flip LookupAlloy** (bulk + TS + `_t`; Maxwell + thermal trees). Simplest family:
  exercises T routing, `compute_rho_bulk`, `density`, memoization, `save_resistivity` — and
  the TS-alloy-uses-bulk-kernel property (T15). Verify per §4.1. *(after: R5)*
  **Status 2026-07-19:** the flip LANDED silently in `abb71c2e` (2026-07-15 "backup",
  unrecorded here) at all 4 sites, but the TS-alloy half was defective (D13 — sidesets
  have no MaxwellData) and was reverted to legacy `h_alloy`/`h_alloy_t` on 2026-07-19.
  Currently in tree: bulk-Conductor `h_calc` (`cl_IWG_Maxwell.cpp:329`) + both thermal
  `T_h_calc` sites (`cl_IWG_MaxwellThermal.cpp:83,140`); TS-Maxwell legacy. **§4.1
  verification has NOT been run** — no baseline commit, no iteration-count parity check.
  Remaining: verify the three live flips, then re-flip TS-alloy once D13's sideset
  design lands (R9 prerequisite).
- [x] ~~**R7 — Flip PureMetal**~~ SUPERSEDED by the 2026-07-21 leap (bulk/TS/`_t`): bj_angle β, metal rho/lambda, and no Newton
  kernel dispatch under the §6.0 `have(jc)` split. *(after: R6)*
- [x] ~~**R8 — Flip HTS bulk**~~ SUPERSEDED by the 2026-07-21 leap (± defect ± piecewise ± `_t`): π/2 dummy β, route to
  `h_newton`, defect x/y/z/time, peer-T. Iteration counts must match baseline exactly
  (§4.1 criterion — the Newton path is the sensitive one, paper1 Eq. 13/ε<10⁻¹¹).
  *(after: R7)*
- [x] ~~**R9 — Flip HTS thin-shell**~~ SUPERSEDED by the 2026-07-21 leap (± defect ± piecewise ± `_t`): bn_angle β, route to
  `h_newton`, compute_bn, n-vector currency assert. Verify corc. *(after: R8)*
- [x] ~~**R10 — Flip default/UserDefined**~~ SUPERSEDED by the 2026-07-21 leap (retire generic `h`/`h_ts` routing): exercises
  the O1 full-signature policy; UserDefined superconductors route by `have(jc)` to
  `h_newton`. Verify Tape_Quench/CustomMat (piecewise + defect + `.so` material). Expected
  FP-tolerance deviation: M-mu0 placement only (§3.1-T2). *(after: R9)*
- [x] **R11 — Thermal tree flips** mirroring R6→R10 order for existing `T_h_*` bodies into
  `T_h`; `T_h_newton` is the separate thermal-Newton target from
  `thermal_matrices_cleanup_and_newton_plan.md` §3/R4-R5. *(after: R10 — or interleaved per
  family with R6-R10 if the shadow harness is in; note the interleave in this file either way)*
  **2026-07-21 dispatch collapse LANDED (Fable session per
  `todo/closed/handoff_thermal_iwg_collapse_session.md`, Codex-audited, Christian approved):**
  `IWG_MaxwellThermal::link_to_group` Conductor+ThinShell cases collapsed to an
  unconditional `T_h_picard` pick (renamed from `T_h_calc` — Christian's ruling: B/β are
  frozen per thermal solve, not thermal dofs, so the collapsed kernel IS the Picard kernel
  and `T_h_newton` needs only the dcp/dT, dλ/dT, dρ/dT channels); Buffer/Ferro/Air keep
  `T_phi`; dead Controller/Kernel includes dropped. D2's legacy Conductor-dispatch residual
  is GONE. All 13 legacy `T_h_*` bodies are now unreachable from dispatch — deletion
  pending the §4.1-style verification gate (handoff S4/T5, baseline `ce6a0e8b`).
  Codex-audit deltas (thread `tmp/ai_exchange/thermal_iwg_collapse.md`): T upper clamp
  and UserDefined full-signature routing = closed O2/D4 policy decisions; NEW cross-cutting
  defect D14 (Alloy dependency bits) found and fixed same session.
  **2026-07-21 (later): verification gate PASSED (Christian) and the 13 legacy bodies
  DELETED** — `mt_thermal_h.cpp` 1339 → 55 lines, hpp 77 → 40 (keep-set: `T_h_picard`);
  straggler sweep clean. Christian also reported the thermal Buffer/Ferro material gap
  (`T_phi` null-deref); `Controller::set_thermal_kernel` now shares the maxwell-side
  material pointer for thermal blocks lacking one (`cl_FEM_Controller.cpp:1818-1828`).
  R11 is DONE.
- [◐] **R12 — Cutover + cleanup.** Point dispatch cases at the final kernels
  (Conductor/ThinShell reduce to the §6.0 two-way picks), delete the 37 retired bodies
  (24 `h_*` — all but `h_ghost` — plus all 13 `T_h_*`) + `get_thermal_calculator`, strip
  the shadow harness or leave it compiled-out, update
  `src/fem/maxwell/doc/` (kernel architecture note), devlog, move this file toward closed.
  **2026-07-21 Maxwell half DONE (Claude, approved):** all legacy `h_*` bodies +
  `h`/`h_ts`/`h_calc`/`h_newton_calc` + `get_thermal_calculator` + the dev pragmas
  deleted from `mt_maxwell_h.{hpp,cpp}` (2724 → 325 lines; keep-set: `h_picard`,
  `h_newton_mu0`, `h_newton_mu`, `h_ghost`, `save_resistivity`, `get_resistivity`);
  the legacy `compute_bn` wrapper removed from `cl_FEM_Calculator.hpp` (zero callers).
  Straggler sweep clean. **Open Newton-tangent follow-up (attribution corrected by
  Christian same day):** with constant jc/n the powerlaw tangent is exact; the missing
  J-channel is the METAL/ALLOY layers' ρ(T,|B|,β) field derivative (Kohler/table,
  `drhodj = return_zero`). Consumers for the new Metal/Database dB/dβ machinery still
  needed: MaxwellData `compute_drhodb/dbeta` dispatchers + E-channel block in
  `h_newton_mu0`/`h_newton_mu`. When jc/n become lookup tables, add
  `drho_powerlaw_dB/_dbeta` (∂ρ_PL/∂jc chain + parallel-combination factor). **Consequence for the §6.1 sign ruling:** the TS normal-metal
  `−E·q` kernels are gone — the pending document-vs-patch decision now applies ONLY to
  the `mt_thermal_h.cpp` mirrors, which die with the thermal-tree collapse (R11
  remainder: 13 `T_h_*` + thermal dispatch still legacy).
  **2026-07-21 thermal half DONE (Fable session, Christian approved + gate passed):**
  thermal dispatch now collapses to `T_h_picard`; all 13 `T_h_*` bodies, hpp
  declarations, and dead includes were deleted (`mt_thermal_h.cpp` 1339 → 55 lines);
  straggler sweep clean.
  **§6.1 sign-split CLOSED:** the last legacy `−E·q` mirrors are gone — the
  document-vs-patch decision is moot; the corrected convention (h = +E·q via
  MaxwellData) is now the only implementation in tree.
  Final gate: full §4.1 matrix on ≥2 MPI ranks for helix + corc. *(after: R11)*

### 4.0 Implementation Progress (updated 2026-07-14)

**Implemented and audited (pre-R1, wiring layer only — kernels untouched):**
- MaxwellData wiring is delegation-ready: D1 Air gates, D2 Maxwell-side material binding,
  D5/D8/D9 three-way `link_element_*` per-element dispatchers with null-first peer guards,
  D10/D11 shared `select_link_element_dispatcher()` from both `allocate()` and
  `link_maxwell()`, intpoint asserts restored. Defect history with citations: §3.2.
- O10 executed in parallel (Opus session, Codex+Grok 0 defects): the powerlaw family is
  now `( normJ, T, normB, angleNxB [, x,y,z,t ] )` across 8 decls / 8 impls / 16 helper
  sites / 48 legacy callers — the helper and the legacy kernels already share the new
  argument order, shrinking the R5 diff.
- O1-O10 are closed; the Newton-split kernel pair is adopted (§6.0).

**R1 smoke PASSED, R2 implemented + audited, R3/R4 cut, R5 built and confirmed
(all 2026-07-14).** Also landed en route: B7/qhist consistent-tangent fix (Codex: no
blocking defect) and the O11 ferro-conductor dispatch semantics. **Still pending:** the
pre-flip commit, then the R6-R11 dispatch flips (one commit per family), R12 cutover;
tails: R2d Controller diagnostic, rank-1 phi_ferro tangent commit, O11 kernels +
validation case.

### 4.1 Regression matrix & per-step verification

Cases (all present in-tree, checked 2026-07-13):

| Case | Path | Exercises |
|---|---|---|
| helix | `examples/helix/` (has `topology{…periodic…}`, input.conf:54-67) | bulk conductor, cross-block periodic edges (memory: regression for collect_edges fallback), MUMPS |
| corc | `examples/corc/` (`thinshell : tape`, input.conf:91) | thin-shell HTS, ghost facets (h_ghost/element_rho), cuts |
| Tape_Quench/BuiltinMat | `tmp/examples/Tape_Quench/BuiltinMat/` | Maxwell+thermal (`hphiTrun` path), metal/HTS `_t` variants, T_h family |
| Tape_Quench/CustomMat | `tmp/examples/Tape_Quench/CustomMat/` (`custom` + `piecewise` + `defect`, input.conf:64-73) | UserDefined superconductor, D4/O1 branches, defect+time |
| Validation set | `tmp/examples/Validation/{Superconducting_Slab_2D,Superconducting_Wire_2D,Wire_2D}` | 2D bulk HTS/metal, bj_angle_2d, 2D-j path, T8 |

**Verification (revised 2026-07-14 after the R3/R4 cut — run-based, git as baseline):**
per flipped family, run every §4.1 case that reaches the family and check:
1. the solve **converges with the same per-timestep iteration counts** as the pre-flip
   commit (the most sensitive cheap indicator — the Newton tangent is the checkerboarding
   guard, Messe et al. 2023 (paper1), Eq. 10-13); expected exceptions: the documented
   T2/T8 FP notes and O2/O3 clamp deviations (log when a clamp flag fires);
2. output fields and `element_rho` look physically sane (spot check in ParaView);
3. on any suspicion: checkout the pre-flip commit, rerun the case, and diff the outputs
   directly — the commit replaces the frozen R4 baseline;
4. helix additionally re-run on 2 and 4 ranks (periodic partner ghosting, see memory notes).

One flip per commit, so `git revert` cleanly undoes exactly one family.

---

## 5. Open Design Questions

- [x] **O1 — jc-dependency policy for the material API.
  RESOLVED 2026-07-13 → decision (Christian): full-signature policy, routed internally on
  the existing dependency machinery.** The assembly path always calls the full
  `( normJ, normB, angleNxB, T [, x,y,z,t ] )` interface; **inside** the material layer the
  routing keys on `Material::depends()` / `Material::dependencies()`
  (`cl_Material.hpp:354,357-362`; per-property bitset storage `:249`) — which is already the
  class's stated design intent ("Automatic property function dispatching based on
  dependencies", `cl_Material.hpp:202`). Correct dependency flags — from users and from the
  built-in library alike — are the contract that makes this clean. **Also decided: the
  reduced rho-overload set on the base Material class stays.** Christian weighed collapsing
  it (the many rho functions are the remaining annoyance) and ruled that a rarely-used API
  is less trouble than making new-material implementation annoyingly difficult —
  user-friendliness for material implementers wins; do not re-propose slimming the overload
  zoo. Implementation items (→ R2): relax the two 8-arg defect full-dependency asserts to
  dependency-keyed routing (`powerlaws.hpp:229-238`, `:717-734`, fallback pattern
  `:199-209`); verify `JcFunction::eval( normB, angleNxB, T )` per subclass for
  partial-dependency functions (Codex check item); policy paragraph into
  `src/physics/materials/doc/`. *(Reframed history and the original parity analysis stay
  below as the coverage map.)*

  **Proposal (Claude, 2026-07-13 — ACCEPTED same day with the dependency-machinery
  refinement recorded in the resolution above): the "full-signature policy".**
  *One sentence: the assembly path always calls the full HTS interface —
  `rho_powerlaw` / `rho_piecewise` `( normJ, normB, angleNxB, T [, x, y, z, t ] )` and their
  `drho_*_dJ` twins — and the material consumes the arguments its jc/n functions depend on
  and ignores the rest; ρₙ (the normal-matrix branch of the parallel combination) is always
  evaluated at the local T.*

  Why this shape: it is what MaxwellData already does (zero helper changes); the caller
  should not know a material's dependency set — the material should (information hiding,
  same once-per-setup fn-ptr dispatch pattern as MaxwellData itself); and it upgrades the
  legacy quirk where a partial-dependency material got ρₙ(gTbulk) even in a thermal problem
  (`powerlaws.hpp:156` vs `:211`) — physically, the matrix resistivity should track local
  temperature. In magnetics-only problems T == gTbulk, so nothing changes there.

  Mechanics (small, all in `src/physics/materials`):
  1. Relax the full-dependency asserts in the 8-arg defect overloads
     (`powerlaws.hpp:229-238` powerlaw; `:717-734` piecewise twin) to the same
     constants/partial fallback the 4-arg already has (`:199-209`).
  2. Guarantee `JcFunction::eval( normB, angleNxB, T )` is well-defined for
     partial-dependency functions — either eval ignores unused arguments, or
     `set_jc_function` binds an internal dispatch pointer once at setup (verify per
     JcFunction subclass; check item, not assumed).
  3. State that user superconductors provide jc via a `JcFunction` or constants; the
     `jc_custom(T)` path stays builtin-only (served by the `(normJ,T)` overloads).
  4. Keep the reduced overloads public for postprocessing/tests, documented as
     non-assembly; the legacy generic-kernel dependency branches die with the kernels
     at R10.

  Behavior deltas vs legacy, to disclose (uncoupled runs are bit-identical since
  T == gTbulk): (i) partial-dep UserDefined superconductors in **coupled** problems get
  ρₙ(T_local) instead of ρₙ(gTbulk) — deliberate improvement, small result change;
  (ii) partial-dep **defect** materials no longer assert-fail in debug. Acceptance per
  Christian's framing: Tape_Quench/CustomMat and IV_Characteristic run and pass the
  shadow compare with exactly these two documented exemptions; the policy paragraph lands
  in `src/physics/materials/doc/`. Rejected alternative: dependency-keyed fn-ptr matrix in
  the MaxwellData ctor (up to 24 extra targets — reproduces inside the helper the
  combinatorial explosion this plan exists to kill, and preserves the gTbulk quirk).
  Mapping proof (verified): the trees branch on `MaterialType` (`cl_Material.hpp:79-88`);
  the helper classifies by `have(jc)` (set by `set_jc_function`, `cl_Material.cpp:651`, and
  by UserDefined SC config, `cl_Material_UserDefined.cpp:101,195`) and
  `depends(rho,normB)` (`cl_FEM_Calculator.cpp:85-90`). Hence: `HTS` ⇒ have(jc) ✓;
  `PureMetal` ⇒ ¬jc ∧ depends(rho,normB) ✓ (Kohler magnetoresistance);
  `LookupAlloy` ⇒ ¬jc ∧ ¬depends ⇒ bulk ✓. The partitions **coincide for the built-in types
  but not for UserDefined**: a UserDefined superconductor is `default:`→generic in the trees
  but HTS-family in the helper — which is *fine* iff the helper reproduces the generic
  UserDefined branches, which today it does not (D4). Decision needed: (i) exact-parity ctor
  branches keyed on `type()==UserDefined` + `depends(jc,·)` [recommended — provable, ugly],
  vs (ii) key *all* HTS materials on jc-dependencies [cleaner, but changes behavior for
  built-in HTS with partial-dep Jc functions, cf. `set_jc_function` setting only the deps the
  function has, `cl_Material.cpp:637-649`]. Recommendation: (i) now, revisit (ii) after R12.
- [x] **O2 — T-clamp ownership and its consistent derivative.
  RESOLVED 2026-07-14 → decision (Christian confirmed Claude's proposal below): unified
  clamp inside `MaxwellData::compute_T` (both sides), `mTClamped` flag lives in MaxwellData,
  all dT-derivative accessors return 0 while clamped, converged-at-clamp diagnostic;
  code → **R2d**.**
  Analysis record (2026-07-13, physics rationale Christian, analysis Claude):

  *Why the clamp exists (Christian):* highly nonlinear iterates — especially in the low-cryo
  region — can transiently swing to negative T (and past the upper table edge). The math
  allows it; the physics doesn't, and the material table raises a physics error. So property
  evaluation must clamp on **both** edges. The clamp guards table **inputs**; it is an
  iterate-transient guard, not part of the converged physics.

  *Is the derivative zero at the clamp? — Yes, and it has a name.* This is the
  **algorithmic (consistent) tangent** principle: Newton must linearize the algorithm the
  residual actually evaluates, not the underlying smooth function. The canonical instance is
  return-mapping plasticity, where the tangent of the *projected* stress update — not the
  continuum modulus — is what delivers quadratic convergence (Belytschko, §5.9, "algorithmic
  moduli consistent with the underlying stress update scheme"; Zienkiewicz & Taylor Vol 2,
  §4.4 — "convergence is very rapid, especially when consistent tangent moduli are
  available" — citing Simo & Taylor 1985). Property evaluation through a clamp is exactly
  such a projection: `p̃(T) = p( clamp(T) )`, so `dp̃/dT = p'·(dclamp/dT)` = **0 whenever
  clamped**. Using the smooth derivative instead (an *inconsistent* tangent) predicts
  residual changes the clamped residual will not deliver → oscillation near the bound.
  Common-practice notes: (i) the zero-slope plateau does not strand Newton, because the
  residual itself still sees the unclamped T dof (only property evaluation is clamped) and
  the relaxation/line-search layer (paper1, Eq. 12-13) handles the pull-back; (ii) the
  alternative — linear table *extension* beyond the range — keeps slope information but
  evaluates unphysical values and must extend value+derivative together; unbounded, not
  recommended as default; (iii) a C¹ "soft clamp" is the fallback only if kink chatter is
  ever observed; (iv) a **converged** solution sitting at a clamp is a modeling error —
  emit a once-per-step diagnostic instead of silently accepting it.

  *Proposed decision:* the clamp lives **inside `MaxwellData::compute_T`** — unified for
  BOTH maxwell-side and thermal-side instances. (Christian's table-error rationale applies
  identically to the Maxwell `rho(T,…)` call; the legacy raw-T maxwell `_t` kernels
  (`mt_maxwell_h.cpp:204`) were exposed to the same table error — an oversight, not a
  contract.) `compute_T_fem` memoizes T **plus an `mTClamped` flag**, and every
  dT-derivative accessor — `compute_dcpdT`, `compute_dlambdadT`, and the future `drhodT`
  for `T_h_newton` — returns **0 when `mTClamped`**: value and derivative are owned by the
  same algorithm in one place. Bounds: lower `gTmin` (`typedefs.hpp:65`); upper from the
  material's `T_max` property where defined (`MaterialProperty::T_max` exists; there is no
  global gTmax) — cache the window at construction. Behavior change vs legacy (maxwell-side
  clamping) is documented; it fires only in pathological transients.
- [x] **O3 — rho clamp ownership. RESOLVED 2026-07-14 → decision (Christian): NO separate
  `compute_rho_clamped` accessor — mirror the O2 pattern instead: `compute_rho(k)` clamps
  internally into `[gRhoMin, gRhoMax]` and sets an `mRhoClamped` flag; both flags
  (`mTClamped`, `mRhoClamped`) live inside the MaxwellData container; code → **R2e**.** One rho value for
  every consumer (Maxwell K, Joule source, `element_rho` mean).

  *Behavior-neutrality check (Claude, 2026-07-14):* the universal clamp is inert on the
  Maxwell K path in practice — the power law's own floor is `Material::mRhoMin = 1e-16`
  (`cl_Material.hpp:255`), numerically identical to `gRhoMin`, and the upper bound is
  unreachable (the parallel combination is bounded by ρₙ ~ 1e-6; metals/alloys sit decades
  below 1e10). **One instructive subtlety:** the parallel combination
  `1/(1/ρₙ + 1/ρ_PL)` is strictly *below* `min(ρₙ, ρ_PL)`, so when ρ_PL sits at its floor
  the result dips a hair under 1e-16 and the clamp **does** fire in deep-subcritical
  states → `mRhoClamped = true` there.

  *Consistent-tangent consequence (deliberate, better than legacy):* while `mRhoClamped`,
  `compute_drhodj` returns **0** — same principle as O2. Legacy computed `drho_*_dJ`
  unconditionally, including at the floor where the true derivative of the clamped value is
  zero; the new scheme zeroes the Newton term exactly where ρ is pinned (deep-subcritical
  elements are effectively linear there, and this removes the derivative swings that
  destabilize Newton at extreme ρ). Expected, documented shadow-compare deviation:
  `dKdx_times_x` on deep-subcritical elements. The future `drhodT` (T_h_newton) returns 0
  while `mRhoClamped || mTClamped`. Clarification record (2026-07-13):

  *What this clamp is:* `gRhoMin = 1e-16`, `gRhoMax = 1e10`, set in the executables
  (`hphirun.cpp:56-57`, `hphiTrun.cpp:48-49`). Two consumers only: (1) the thermal Joule
  source `f += w·Nᵀ·ρ·|j|²·dV` — every `T_h_*` clamps ρ first (`mt_thermal_h.cpp:137,337,
  409,…`); (2) the ghost-penalty resistivity read `get_resistivity`
  (`mt_maxwell_h.hpp:32-37`). The Maxwell **K** path deliberately uses raw ρ. (The power
  law's internal `mRhoMin` floor, `powerlaws.hpp:108`, is a separate material-level
  regularization, not this clamp.)

  *Role, distinct from O2:* the T-clamp guards table **inputs** (physical validity of T);
  the ρ-clamp guards the power-law **output** against J-driven swings — HTS ρ spans ~15
  decades, and transient iterates can push the source term into overflow/conditioning
  trouble. Both are iterate-transient guards; O2's converged-at-clamp diagnostic idea
  applies here too.

  ~~*Sharpened proposal (superseded by the resolution above):* separate
  `compute_rho_clamped(k)` accessor for the source term, raw `compute_rho` for the Maxwell
  K~~ — rejected 2026-07-14 (Christian): two rho values for one quantity is the more
  complicated design; the neutrality check above shows the raw/clamped distinction never
  matters in practice, so one internally-clamped `compute_rho` + `mRhoClamped` wins.
  `element_rho`/`get_resistivity` keeps its own read-side clamp unchanged (O5).
- [x] **O4 — density accessor. RESOLVED 2026-07-13 → decision (Christian): density at the
  undeformed-mesh reference temperature is a PHYSICS REQUIREMENT, not a convenience.**
  Although real density changes with temperature, all calculations run on the **undeformed
  mesh**, and the transport properties are already **corrected for thermal expansion** —
  so the mass matrix must use the density at which the mesh is not deformed, which is by
  default room temperature. Evaluating `density(T)` per point would double-count expansion:
  this is a **dangerous physics trap** for anyone touching these kernels. Implementation
  unchanged from the earlier recommendation: ctor-time
  `mDensity = have(density) ? density(gTroom) : ref_density()` + trivial getter
  (`mt_thermal_h.cpp:90,403` legacy forms), with a trap-warning comment at the member.
  **Doc item (→ R2g; accessor → R2f):** record this circumstance in
  `src/physics/materials/doc/materials_usage_guide.md` — density is evaluated at the
  undeformed-mesh reference temperature; transport properties are expansion-corrected.
  The fallback-form note stands: materials lacking `density` that would assert in the
  specialized legacy variants now use `ref_density()` (benign extension, flag in R2 notes).
- [x] **O5 — element_rho mean home. RESOLVED 2026-07-14 (Christian: "why not store
  element_rho_mean onto the mesh?" — it already is, and it stays there).** The mean lives in
  the mesh field `"element_rho"` today (`save_resistivity` writes
  `mesh()->field_data("element_rho")(element index)`, `mt_maxwell_h.hpp:26-30`; registered
  NonDof, `cl_IWG_Maxwell.cpp:77`; hidden from file output, `cl_MaxwellFactory.cpp:699-703`).
  The only question O5 ever asked was who *accumulates* the mean after the collapse — answer:
  the kernel, unchanged (`rho_m += w·ρ·dV; V += w·dV; save_resistivity(rho_m/V)`), using the
  helper's memoized `compute_rho(k)` values. h_ghost keeps `get_resistivity`'s read-side
  clamp (`mt_maxwell_h.hpp:32-37`). Helper-side accumulation rejected: the helper is
  per-point, the mean is per-element state.
- [x] **O6 — thin-shell alloy. CONFIRMED OBSOLETE 2026-07-14 (Christian).** No decision was
  ever needed — T15 shows bulk==TS for T-only rho; the collapse makes the question
  disappear. Kept only so audits confirm rather than assume.
- [x] **O7 — peer-linking API shape. RESOLVED 2026-07-13 by implementation (Christian):**
  folded into `Calculator::link(Element*)` via the three-way `mFunLinkElement` dispatcher —
  the opposite of Claude's recommendation, and the implementation addressed both stated
  objections: groups that never touch T get `link_element_maxwell`/`link_element_default`
  (zero peer logic — no fire-on-every-group cost), and the mutual-recursion hazard was
  disproved (depth-2 bound, D5) then hardened (D8/D9/D11 fixes). Kernels need no per-element
  peer call at all — strictly better than `link_peer()`.
- [x] **O8 — postproc scope. DIRECTION DECIDED 2026-07-14 (Christian), work deferred:**
  address the postprocessors *after* the main kernels are fixed, and when doing so their
  behavior and design must **mirror the main kernels** in terms of matrix computation —
  i.e. route through MaxwellData with the assembly conventions. That mirroring will also
  settle the convention discrepancies found this session:
  `MaxwellPostprocessor::compute_superconductor{,_ts}` currently recompute h/b/j/T/jc with
  their own conventions — `mH = E·q` (no −mu0), `Temp = norm(N·T)` instead of `dot`, a
  signed (non-abs) β (`cl_MaxwellPostprocessor.cpp:530-699`, β at `:661`) — and unifying on
  the helper changes plotted output (a fix, per this decision). Follow-up todo to be opened
  after R12; out of this plan's scope.
- [x] **O9 — rollback/coexistence. RESOLVED 2026-07-14 (Christian): coexistence is not an
  issue — soft cutover confirmed.** Legacy variant functions stay compiled during migration
  and get deleted once the new system works (R12). Each Rn family flip is a
  one-line-per-case dispatch change, so rollback = revert that line. No runtime switch; the
  shadow harness (R3) covers A/B needs.
- [x] **O11 — domain type vs material μ: dispatch semantics for ferro-conductors.
  DECIDED 2026-07-14 (Christian).** Within BELFEM naming logic, **"Ferro" means
  field-dependent μ but NO current** — a topological role (φ-region; cohomology cuts may
  penetrate a Ferro block). Setting an iron material on a **Conductor** domain is legal:
  it declares a current-carrying region (h-dofs, cuts terminate) whose μ happens to be
  field-dependent. Therefore the dispatcher, when linking the calculation function to a
  block, must check the MATERIAL — `is_constant(mu) && constant_property(mu) == μ0` →
  the existing `h_calc`/`h_newton_calc`; field-dependent μ → new variable-μ h-kernels
  (mirroring the fixed `phi_ferro`: per-point `μ(H)` mass + `dMdx_times_x`/`dMdx_times_h`
  tangent blocks; exact rank-1 form `u uᵀ·μ′/H`, `u = Eᵀ(Eq)`, available at the cost of
  the existing HTS Newton term). This is the O1 property-over-type principle applied to
  the μ axis, and MaxwellData's ctor already computes the needed flags
  (`cl_FEM_Calculator.cpp:116-129`). **Corollary (silent gap today):** iron-on-Conductor
  currently runs with μ0 in the mass term — the legacy `h_*` kernels never inspect μ.
  Remaining before implementation: a validation problem (no legacy oracle — new physics;
  Christian to pick), and the Q2 history fix is prerequisite (landed 2026-07-14, see
  `closed/bdf_nonlinear_mass_verification.md` B7).
- [x] **O10 — rho/lambda argument convention. RESOLVED + EXECUTED 2026-07-14** (decided by
  Christian — option (c) unify, normJ-first placement; executed by a parallel Opus session;
  Codex + Grok read-only audit: **0 code defects**; tick approved by Christian after
  Claude's verification). The last un-unified family (Jc/powerlaw) moved:
  `rho_powerlaw`/`rho_piecewise`/`drho_powerlaw_dJ`/`drho_piecewise_dJ` T-bearing overloads
  are now `( normJ, T, normB, angleNxB [, x,y,z,t ] )` (verified in tree:
  `cl_Material.hpp:603-606`; caller spot check `mt_maxwell_h.cpp:433`
  `rho_powerlaw( norm_j, gTbulk, norm_b, beta )`). Reduced overloads unchanged;
  `JcFunction::eval( normB, angleNxB, T )` internal order deliberately kept (separate
  contract). Full record: `todo/rho_lambda_argument_convention.md` (Status updated there) and
  the O10 entry in this plan's exchange thread. **Remaining tail (tracked in that todo, not
  here):** the materials doc pass — guides still show the old powerlaw order.

---

## 6. Target Design (concrete)

### 6.0 Kernel granularity — Newton split (adopted 2026-07-13, Christian)

Instead of one kernel with a runtime `drho != 0` guard, provide **two versions of each
kernel, with and without the Newton part**:

- `maxwell::h` (non-Newton: metals/alloys and any non-HTS material) and
  `maxwell::h_newton` (adds the `dKdx_times_x` block) — dispatch on `have(jc)`;
- `T_h` and `T_h_newton` — where **`T_h_newton` is new physics, not a port**: the thermal
  Newton contribution has never been implemented (the IWG ctor sets the derivative-matrix
  flags with a `//todo`, `cl_IWG_MaxwellThermal.cpp:29-32`) and becomes buildable now with
  the new material API (`compute_dcpdT`/`compute_dlambdadT` exist,
  `cl_FEM_Calculator.hpp:2369-2400`; a dρ/dT sibling is the missing piece). Term recipes and
  math live in `thermal_matrices_cleanup_and_newton_plan.md` §3/R4-R5, which rebases onto
  `T_h_newton`.

This supersedes Appendix A's single-kernel answer **in part**: the material/geometry axes
still collapse into the helper; only the Newton axis stays a kernel-level split (it changes
which matrices are written, not how material values are computed). Benefit: no dead
`compute_drhodj` call on metal/alloy blocks, and the Newton term reads as its own kernel.
Cost: one more dispatch line per tree — dispatch stays a two-way pick.

### 6.1 Collapsed Maxwell kernels (replaces all non-ghost `h_*`)

Sketch below shows the `h_newton` body; `h` is the same without the `dKdx` block and without
the `compute_drhodj` call (§6.0):

```cpp
// mt_maxwell_h.{hpp,cpp} — bulk+TS+all-materials kernel, Newton variant
void
h_newton( Calculator * aCalc, TimestepMatrices * aMatrices )
{
    const Vector< real > & w = aCalc->integration()->weights();
    calculator::MaxwellData * mx = aCalc->maxwell();

    // peer element already linked: Calculator::link() dispatched through
    // link_element_maxwell (cl_FEM_Calculator.cpp:2336-2360; D5, needs D8+D9)

    real rho_m = 0.0;
    real V     = 0.0;

    for ( uint k = 0; k < aCalc->num_intpoints(); ++k )
    {
        const Matrix< real > & E = aCalc->E( k );
        const Matrix< real > & C = aCalc->C( k );

        // keep shared work vectors ("b","j",…) fresh for downstream readers (T9)
        // — compute_rho pulls b/j lazily itself; rho path decides what it needs
        real rho  = mx->compute_rho( k );      // full β/T/defect dispatch inside
        real drho = mx->compute_drhodj( k );   // h_newton is reached only for have(jc)

        real wdV = w( k ) * aCalc->dV( k );
        rho_m += rho * wdV;
        V     += wdV;

        aMatrices->M() += trans( E ) * E * ( constant::mu0 * wdV );   // per-term mu0 (T2)
        aMatrices->K() += trans( C ) * C * ( rho * wdV );

        if ( mx->norm_j( k ) > BELFEM_EPSILON )        // legacy guard; h_newton is HTS-only (§6.0)
        {
            const Vector< real > & j = mx->compute_j( k );
            Ctj = trans( C ) * j ;                     // "Ctj" workspace (e x 1), bound
            aMatrices->dKdx_times_x() += Ctj * trans( Ctj )
                                       * ( ( drho / mx->norm_j( k ) ) * wdV );
        }
    }
    save_resistivity( aCalc, rho_m / V );                             // O5
}
```

Notes: identical accumulation order and factor grouping as the specialized legacy bodies
(`w(k)*rho*dV(k)` association preserved via `wdV` — verify bit-equality in R3; if the
`rho*wdV` regrouping vs legacy `w(k)*rho*dV(k)` shows one-ulp diffs, revert to the literal
legacy grouping). `Ctj` binds the preallocated workspace by reference OUTSIDE the loop
(`Matrix< real > & Ctj = aCalc->matrix( "Ctj" );` — created e×1 in
`create_custom_vectors_and_matrices`); per-point `Matrix` construction is forbidden
(Christian 2026-07-21, hidden-allocation rule — the earlier "legacy-identical temporary"
note is superseded).

> **Sign-convention ruling (Christian 2026-07-21, supersedes the F3 open question):**
> the operator matrices N, B, E, C are physics-agnostic. Field reconstruction rules:
> **h = −B·φ** (gradient operator), **h = +E·q** (edge interpolation), **j = C·q**.
> Legacy kernels computing `h/b = −µ·E·q` are WRONG; MaxwellData's H-first paths and the
> postprocessor implement the correct convention. Cancellation sweep (Claude 2026-07-21,
> Codex adversarial verification same day): the legacy sign error is silent for **bulk
> kernels** (consumers quadratic/norm-based/abs-protected — full-b flips masked by
> `abs(dot)` in `bj_angle`, `cl_FEM_Calculator.hpp:2250-2266`), for **TS HTS kernels**
> (`bn_angle(b,n)`: bn ∥ n correct-signed, bt ⊥ n — `|b|` and `dot(n,b)` invariant),
> and for ghost/φ/interface/BC/constraint paths (no signed consumer). **NOT provably
> silent for TS normal-metal paths** (Codex counterexample, verified): `h_ts_metal`
> (`mt_maxwell_h.cpp:249-268`) and thermal mirrors mix `bt` (wrong sign) + `bn`
> (correct sign) into `bj_angle(b,j)` — `abs` protects a FULL flip, not a partial one:
> `abs((−bt+bn)·j) ≠ abs((bt+bn)·j)` whenever both `bt·j` and `bn·j` are nonzero, and
> TS `j = C·q` has no local guarantee of zero normal component. Affects β-dependent
> metal rho (Kohler magnetoresistance) on TS layers only. Decision pending (Christian):
> (a) document the restriction and live with it until R12, or (b) flip `bt = +µ0·E·q`
> in the TS normal-metal kernels + thermal mirrors now (aligns with the ruling; shifts
> TS-metal β baselines slightly toward correct physics). Either way the R12 cleanup
> deletes the convention split and the legacy `compute_bn` wrapper (bn = µ0·hn).

### 6.2 Collapsed thermal kernel (replaces all 13 `T_h_*`)

```cpp
void
T_h( Calculator * aCalc, TimestepMatrices * aMatrices )
{
    const Vector< real > & w = aCalc->integration()->weights();
    calculator::MaxwellData * mx = aCalc->maxwell();    // thermal-side instance

    // Maxwell peer element already linked via link_element_maxwell (D5, needs D8+D9)

    for ( uint k = 0; k < aCalc->num_intpoints(); ++k )
    {
        const Matrix< real > & B = aCalc->B( k );
        const Matrix< real > & N = aCalc->N( k );

        real cp     = mx->compute_cp( k );              // T via clamped compute_T (O2)
        real lambda = mx->compute_lambda( k );          // depends(lambda,normB) routing (T10)
        real rho    = mx->compute_rho( k );             // clamped internally, mRhoClamped (O3)
        real nj     = mx->norm_j( k );

        real wdV = w( k ) * aCalc->dV( k );
        aMatrices->M() += trans( N ) * ( mx->density() * cp ) * N * wdV;      // O4
        aMatrices->K() += trans( B ) * lambda * B * wdV;
        aMatrices->f() += trans( N ) * ( rho * nj * nj ) * wdV;
    }
}
```

The thermal Newton matrices (`dKdX_times_x` etc., flags set with a `//todo` in the IWG ctor,
`cl_IWG_MaxwellThermal.cpp:29-32`) are **not** filled by `T_h`; they belong in the separate
`T_h_newton` body described by `thermal_matrices_cleanup_and_newton_plan.md` R4-R5.

### 6.3 MaxwellData extensions (all in `cl_FEM_Calculator.{hpp,cpp}`)

| Addition | Purpose | Trap/Defect |
|---|---|---|
| construction gate: material non-null (+ domain check) | air/buffer blocks | D1 |
| material source: `mMaxwellCalculator->group()->material()` | thermal-side binding | D2 |
| β-convention comment + ctor guard (no slot split — policy per D3 resolution) | convention aliasing | D3 |
| ctor UserDefined×depends(jc,·) rho/drho branches (8 targets: {powerlaw,piecewise}×{TnBβ, nBβ, T-only}×defect±) | generic-branch parity | D4/O1/T14 |
| ~~`link_peer()`~~ → landed as `Calculator::link_element_maxwell` via `mFunLinkElement` (Christian 2026-07-13); pending: identity-based peer binding (D8), null guard (D9), intpoint assert | peer currency | D5/T12, D8, D9 |
| unified clamp inside `compute_T` (both sides) + `mTClamped` flag; dT-derivative accessors return 0 when clamped (consistent tangent) | table-input guard | O2/T6 |
| clamp inside `compute_rho` + `mRhoClamped` flag (one rho for all consumers); `compute_drhodj` returns 0 while clamped | power-law-output guard | O3/T7 |
| `density()` (ctor-memoized at gTroom, fallback form; undeformed-mesh physics — see O4) | thermal M | O4/T11 |

### 6.4 Before → after counts

| | before | after |
|---|---|---|
| `mt_maxwell_h.cpp` kernels | 25 + `get_thermal_calculator` | 3 (`h`, `h_newton`, `h_ghost`) |
| `mt_thermal_h.cpp` kernels | 13 | 2 (`T_h`, `T_h_newton` — the latter NEW physics, §6.0) |
| `cl_IWG_Maxwell.cpp` Conductor case | ~92 lines (`:306-397`) | ~3 lines |
| `cl_IWG_Maxwell.cpp` ThinShell case | ~91 lines (`:437-527`) | ~3 lines |
| `cl_IWG_MaxwellThermal.cpp` Conductor + ThinShell | ~105 lines (`:70-174`) | ~8 lines (ThinShell keeps no material lookup at all) |
| β-convention implementations | 4 conventions × ~9 sites | 1 site each (helper) |

---

## 7. Definition-of-Done Checklist

- [ ] Every §3 row mapped to an Rn step or an On question (T-column verdicts all EQ/CHG).
- [ ] D1-D5 fixed and re-audited; D6/D7 closed or waived in writing.
- [ ] O1-O5, O7, O9 decided by Christian and annotated in place (O6 confirmed, O8/O10 spun out).
- [ ] Shadow harness passed on the full §4.1 matrix for every family flip.
- [ ] Baseline diffs + iteration-count parity archived in the devlog for R6-R11.
- [ ] helix on 2/4 ranks + corc pass post-cutover (R12).
- [x] `thermal_matrices_cleanup_and_newton_plan.md` updated to rebase its R4-R5 on the new
  `T_h_picard`/`T_h_newton` pair (2026-07-21, incl. Christian's T-derivatives-only ruling).
- [ ] Codex prose+technical audit of this file; findings folded in with attribution.

## 8. Audit Trail

- Exchange threads: `tmp/ai_exchange/maxwell_kernel_collapse.md` (this plan's audit thread,
  opened 2026-07-13); `tmp/ai_exchange/maxwelldata_helper.md` (helper history, F1-F7
  conventions — F7 = 2D mZ=0 accepted; F6 = lambda routing hoist; both load-bearing here).
- Scouting facts from the task brief were Codex-cross-audited beforehand; every load-bearing
  claim was independently re-verified against the tree this session (all §3 line citations
  re-read, not inherited). Brief corrections found: generic UserDefined branches exist in
  `T_h_ts` only, not bulk `T_h`; the vector-aliasing consumer set is narrower than stated
  (§3.1-T9); `Coords(0,2)` in 2D is out-of-bounds, making F7 a bug fix (§3.1-T8).
- D1-D5 are new findings by Claude/Fable 2026-07-13; **Codex audit round 1 (2026-07-13,
  confidence high) confirmed all five** — refinements folded into §3.2/R3 in place: D1
  air-specific (Buffer requires material, `cl_FEM_Domain.cpp:36-45`); D2 plus the
  pre-existing thermal-Conductor dispatch suspect (`cl_IWG_MaxwellThermal.cpp:70-74`); D4
  extended to the piecewise family with citations; R3 needs an explicit
  `maxwell()->reset()` between shadow passes; D3 reachability bounded (no built-in material
  reaches the alias; UserDefined can). Grok third voice + Codex prose polish still pending.

## Appendix A — Decision: one generic kernel vs per-family kernels

> **SUPERSEDED IN PART (2026-07-13, Christian):** the Newton axis is split back out into
> `h`/`h_newton` and `T_h`/`T_h_newton` — see §6.0. The material/geometry collapse argument
> below stands unchanged; the closing "single kernel means alloy blocks pay dead calls" risk
> is resolved by the split rather than by the runtime guard. Read "single kernel" below as
> "single non-Newton Maxwell kernel plus the HTS Newton companion."

**Question:** collapse to a single `h` (helper dispatches everything) or keep one thin kernel
per material family (h_metal-like shells around helper calls)?
**Answer:** single kernel. The per-family differences are *entirely* inside quantities the
helper already dispatches (β, T, rho, drho, defect coords); a per-family shell would retain
the dispatch tree and re-open convention drift. The tempting middle ground — separate bulk
vs TS kernels — is also unnecessary: `mFunB` (`cl_FEM_Calculator.cpp:116-129,171-184`) keys
bulk/TS at construction from `domain_type()`, and T15 shows even the alloy asymmetry
disappears. Precedent: the helper's own ctor was already structured around exactly these
axes. What stays separate: `h_ghost` (different integrand entirely) and everything
phi-/symmetry-side. **Risk accepted:** a single kernel means alloy blocks pay two dead
virtual-ish calls per point (`compute_drhodj`→`return_zero`, `norm_j` guard short-circuits
before evaluating) — measured cost ≈ one indirect call per point; if profiling ever shows it,
split alloy back out *after* equivalence is proven, not before.
