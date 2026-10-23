# IWG_TransientHeatConduction: Make It Usable End to End

**Date:** 2026-10-22
**Purpose:** Replace the placeholder coefficients of `IWG_TransientHeatConduction` with material data
and make the IWG assemble in production.
**Module:** `src/fem/iwg`
**AIs involved:** Codex `gpt-5.6-terra`/high and Grok `grok-4.7`/high found the gaps (plan jury of
`todo/closed/guide_sweep_fixes_modules.md`, item B11)
**Status:** OPEN — deferred by Christian on 2026-10-22 ("leave it for now"). Nothing implemented.

---

## 1. Current Behaviour

- `compute_mkf` uses a unit mass coefficient and `111.0` as conductivity
  (`cl_IWG_TransientHeatConduction.cpp:101`, `:105`) and fills no load into `f` (`:79`).
- The IWG sets no `MatrixFlag`, so `TimestepMatrices::initialize` never sizes `M` and `K`
  (`cl_TimestepMatrices.cpp:22-42`); `IWG_MaxwellThermal` shows the working pattern
  (`cl_IWG_MaxwellThermal.cpp:27-36`, `:119`).
- `tests/fem/test_Pipette.cpp:408-439` calls `compute_mkf` with no material and asserts
  `Σ M = area`, which holds only for the unit coefficient.

## 2. Steps

- [ ] **R1** — Christian's ruling of 2026-10-22: take ρ·c_p(T) and λ(T) from the material, through
  `Calculator::material()` (its null assert, `cl_FEM_Calculator.cpp:2990`), T interpolated from
  `node_data( "T" )` as in `IWG_StaticHeatConduction` (`:96-106`); ρ by the MaxwellData rule
  (density at `gTroom`, else `ref_density`, times `density_correction`).
- [ ] **R2** — set the `M`, `K`, `F` flags and initialize the matrices from the group.
- [ ] **R3** — give the pipette fixture a material and a temperature field, and change the mass
  oracle to ρ·c_p·area.
- [ ] **R4** — plan jury, code jury (`/code-with-jury`).
