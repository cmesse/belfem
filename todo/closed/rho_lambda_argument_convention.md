# rho / lambda Argument-Order Convention

**Date:** 2026-07-11
**Purpose:** Decide how to handle the deliberate but footgun-prone opposite
argument conventions of the magneto-dependent `rho` and `lambda` families.
**Module:** `src/physics/materials`
**Scope:** Discuss and resolve **within** the ongoing `mt_maxwell_h.*pp` /
`mt_thermal_h.*pp` refactoring — not a standalone change.
> **CONTESTED (2026-08-28).** This file's central decision — option (c), unify on **T-first**
> — is under review for **reversal**. Christian decided 2026-08-27 to move the family to
> uniform field-first (`rho( B, angle, T )`, `lambda( B, beta, T )`); the plan and its
> inventory live in [rho_lambda_b_first_reorder.md](rho_lambda_b_first_reorder.md), audit
> dispatched, nothing implemented. Read that file before acting on anything below.
> **Two corrections to the Tasks section that hold either way:** the materials-module half of
> the doc pass was completed on 2026-08-25 by the three-AI materials documentation overhaul
> (all seven `src/physics/materials/doc/` files rewritten from `( B, angle, T )` to
> `( T, B, beta )`), so the inventory table below overstates what is open by six files; and
> the surviving sites are `src/fem/kernel/doc/dof_manager_usage_guide.md:3118,3154,3205,3269`,
> `src/fem/maxwell/doc/maxwell_usage_guide.md:841,952,958` and `cl_Material.hpp:1200` — whose
> line numbers have all drifted from the ones recorded here. If the reversal lands, those
> seven sites are already correct and the seven materials documents become wrong instead.

**Status:** CODE RESOLVED, **DOC PASS STILL OPEN** — this is the only thing left in this
file, and it is the last surviving item of the whole rho/lambda unification (re-verified
2026-08-09; **inventory re-baselined again 2026-08-11**, see Tasks below — it grew, it did not
shrink). The code side was re-confirmed in the 2026-08-11 sweep at the declaration itself:
`cl_Material.hpp:573/:589` are `rho( const real T )` and `rho( const real T, const real B,
const real beta )`, so every remaining `rho(B, angle, T)` in the guides is teaching an order
the compiler no longer accepts in that sense — it silently means something else.
Originally: RESOLVED 2026-07-11 — option (c) chosen (unify on **T-first**,
`(T, B, beta)`, for rho/lambda/drhodT/dlambdadT). Reorder applied and a 3-AI
audit closed ~7 stale old-order sites (see audit trail below). **Update
2026-07-13:** the residual split (private `mFunctionRhoKohler` family + user
callback still T-last) was eliminated — the entire Kohler pointer family and the
`mUserRhoFunction`/`mUserLambdaFunction` callbacks are now T-first too; 3-AI
unanimous PASS (`tmp/ai_exchange/kohler_pointer_unify_audit.md`,
`devlog/dl20260713_kohler_pointer_order_unify.md`). **Update 2026-07-14:** the
final family — Jc/**powerlaw** — was unified too (option (c), decided by Christian
via the maxwell_kernel_collapse O10 discussion). `rho_powerlaw`/`rho_piecewise`/
`drho_powerlaw_dJ`/`drho_piecewise_dJ` moved from `(normJ, normB, angleNxB, T
[,x,y,z,t])` to `(normJ, T, normB, angleNxB [,x,y,z,t])` — normJ stays first
(matches the existing `(normJ,T)` overload); only the T-bearing overloads changed,
reduced forms untouched. 8 decls + 8 impls + 16 helper sites + 48 legacy callers;
JcFunction eval internal order `(normB,angleNxB,T)` deliberately kept. Codex+Grok
read-only audit: 0 code defects (`tmp/ai_exchange/maxwell_kernel_collapse.md`).
Remaining: the documentation pass (now also needs the powerlaw order updated —
both auditors flagged stale guides/README/contracts).

## Background

The magneto-dependent material accessors intentionally put the temperature
argument in opposite positions:

| Family | Signature | Rationale |
|--------|-----------|-----------|
| `rho`  / `drhodT`    | `(normB, angle, T)` — **T last**  | The Maxwell solve's unknown is `H` (hence `B`); the field is the primary argument. |
| `lambda` / `dlambdadT` | `(T, normB, angle)` — **T first** | The thermal solve's unknown is `T`; temperature is the primary argument. `cp`/`lambda` are thermal properties parameterised by `T`. |

Base decls: `cl_Material.hpp:511,514` (lambda/dlambdadT), `:558,561` (rho/drhodT).
This is a considered design choice (Christian), not an accident.

## The hazard

All arguments are `real`, so a call with the wrong order **compiles silently**
and swaps `T <-> B <-> angle`. This already bit the materials driver (fixed
2026-07-11, see below). Every new call site is a chance to repeat it, and the
two families sit side-by-side in the thermal assembly.

## Options to weigh (decide during the refactor)

- [ ] **(a) Keep + document.** Add a prominent comment on both signature
      clusters stating the convention and why; treat it as house style. Cheapest;
      relies on caller discipline.
- [ ] **(b) Type-safe wrappers.** Wrap the field/temperature in tiny distinct
      tagged structs (e.g. `Temperature{}`, `FluxDensity{}`) so a swapped call
      fails to compile. Safest; touches the accessor surface and all call sites.
- [x] **(c) Unify on one order.** CHOSEN — unify on **T-first** `(T, B, beta)`.
      Public `rho`/`lambda`/`drhodT`/`dlambdadT` are all T-first; the private
      Kohler helpers (`rho_kohler`/`rho_table`) stay in the `mFunctionRhoKohler`
      pointer-contract order `(normB, angle, T)`.

## Tasks

- [x] Reorder all magneto callers to `(T, B, beta)` — 3-AI audit (Claude+Grok+
      Codex) confirmed ~7 stale sites; all fixed 2026-07-11: `mt_thermal_h.cpp`
      rho callers (129/329/409/513); `mt_maxwell_h.cpp` :75 (B/beta swap),
      :2172/:2456; Metal `rho_table`/`drhodT_table` reverted to pointer-contract
      order; Metal/Alloy `lambda`/`dlambdadT` internal calls; base `Material::rho`
      /`drhodT` param names; `main.cpp` :88/:379.
- [x] **Reorder the Jc/powerlaw family to T-first** (2026-07-14, Claude; option (c)
      per Christian, maxwell_kernel_collapse O10). `rho_powerlaw`/`rho_piecewise`/
      `drho_powerlaw_dJ`/`drho_piecewise_dJ`: `(normJ, normB, angleNxB, T [,x,y,z,t])`
      → `(normJ, T, normB, angleNxB [,x,y,z,t])`. Only T-bearing overloads changed;
      2/3/6/7-arg reduced forms untouched; JcFunction `eval(normB,angleNxB,T)`
      internal order kept. `cl_Material.hpp` (8 decls), `powerlaws.hpp` (8 impls),
      `cl_FEM_Calculator.hpp` (16 helper sites), `mt_maxwell_h.cpp`/`mt_thermal_h.cpp`
      (48 legacy callers). Codex+Grok read-only audit: 0 code defects.
- [ ] **Doc pass — STILL OPEN, and larger than first scoped. Line references re-baselined
      2026-08-09** (the 2026-07-14 numbers had all drifted; the earlier list also missed two
      files). The guides still teach the pre-2026-07-11 orders — both `rho/lambda(B, angle, T)`
      and the pre-2026-07-14 powerlaw `(normJ, normB, angleNxB, T)`. Current inventory:

      **Inventory re-baselined again 2026-08-11** (currentness sweep — re-grepped, not copied;
      the 2026-08-09 numbers had drifted in three of the six files and were short by ~15 sites):

      | file | lines |
      |---|---|
      | `src/physics/materials/doc/materials_usage_guide.md` | 46, 68, 73, 76, 202, 203, 270, 400, 404, 455, 472, 479, 481, 517, 729, 730, 927, 973, 1235, 1274, 1286, 1287 |
      | `src/physics/materials/doc/materials_contracts_and_invariants.md` | 200, **314** *(a full `Metal::rho(real B, real beta, real T)` signature — the most misleading site in the set)*, 330, 344, 370, 420, 514, 517 |
      | `src/physics/materials/doc/README.md` | 123, 232 *(both powerlaw; the earlier 87/118/227 no longer point at old-order text)* |
      | `src/fem/kernel/doc/dof_manager_usage_guide.md` | 3192, 3228, 3279, 3343 |
      | `src/fem/maxwell/doc/maxwell_usage_guide.md` | 791, 910, 916 |
      | `src/physics/materials/doc/alloy_transport_mixing.md` | 37 *(the `lambda` line at :73 reads `lambda(B, beta, T)` as **physics notation** in a derivation, not a call — decide per the note below)* |

      **One source site belongs with them** (found in the 2026-08-11 sweep, not fixed — source
      edits are out of this pass's scope): `src/physics/materials/cl_Material.hpp:981` still
      describes "three-parameter properties like `rho(B,angle,T)`" in a doc-comment, while the
      declaration 400 lines above it is `rho( const real T, const real B, const real beta )`
      (`:589`). It is comment text only, so it cannot change behaviour — but it is the one stale
      site a reader meets *inside the header they are calling*.

      Note the two distinct wrong orders: metal/alloy sites need `(T, B, beta)`, powerlaw
      sites need `(normJ, T, normB, angleNxB)`. `alloy_transport_mixing.md` is arguably
      *mathematical* notation rather than an API signature — decide per line whether to
      reorder or to mark it explicitly as "physics notation, not the call signature".
      Deferred from the 2026-07-11 commit (prose-entangled; needs a coherent pass).
- [ ] MaterialData deferred design items (HTS lambda dispatch, drho/dJ dispatcher,
      peer relink, memoization invalidation) — see audit trail; tracked in
      `thermal_matrices_cleanup_and_newton_plan.md`.

## Cross-references

- Refactoring home: [thermal_matrices_cleanup_and_newton_plan.md](thermal_matrices_cleanup_and_newton_plan.md)
  (R4 Material dT-derivative helper, O1 FD-vs-analytic-vs-API derivatives) and the
  planned `MaterialData` consolidation in `cl_FEM_Calculator.*pp`.
- Second-audit thread (derivative math verified, this footgun flagged by all
  three voices): `tmp/ai_exchange/material_derivatives_audit.md`.
- Reorder + MaterialData audit thread (3-AI, the ~7 stale-site findings):
  `tmp/ai_exchange/material_reorder_and_materialdata_audit.md`.
- Related fix applied same day: `is_constant()` NaN comparison
  (`cl_Material.hpp:1391-1393`) — was `!= BELFEM_QUIET_NAN` (always true),
  now `!std::isnan(...)`; load-bearing for the future `MaterialData` dispatch.
