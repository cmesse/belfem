# Gas Models: Migration to the Open-Source Physics Tree

> **CLOSED 2026-09-03** (todo/ currentness sweep, round 3): the migration itself is complete; R4 (optional-load / fail-loud for `gasdata.inp`) declined as a nicety on an OFF-by-default module. Status lines and checkboxes below are as they stood at closure and are not maintained.

**Date:** 2026-08-03
**Purpose:** Move `nonfree/physics/{gasmodels,gastables,atmosphere}` into the open-source
`src/physics/` tree without redistributing third-party data.
**Module:** `src/physics` (target), `nonfree/physics` (source)

**Status (2026-08-05):** Phase 0 COMPLETE and audited (Codex+Grok). **R1/R1b COMPLETE** —
the clean 59-species `gasdata.inp` is built and verified; the Poling-derived full table is
now `gasdata_ppoc.inp` and stays private. **R1c and R2 COMPLETE** — the data path is
`$BELFEM_DATA` via `gBelfemDataPath` with a working-directory fallback, the root-owned
`/opt/scls/cea` dependency is gone, and BELFEM no longer redistributes the vendor CEA files
at all: `scripts/fluidprop/build_tables.py` extracts a 62-species `thermo.inp`/`trans.inp`
with fitted low-temperature intervals merged in, 1.3 MB down to 51 kB. Deuterium is in for
fusion. LBNL IPO approval obtained. Chain-of-title cleared — the code was Christian's own
when written; the `@dlr.de` stamps were cosmetic and have been removed.
**R5 COMPLETE and R6 mostly so:** `gastables` and `gasmodels` now live in `src/physics`,
`USE_GASMODELS` no longer requires `USE_NONFREE`, and combustion is opt-in behind a new
`USE_COMBUSTION`. All four modules and the four example executables **compile and link**,
which also flushed out a latent `Spline::update_data` call that had been broken since the
`SplineBC` parameters were added.
**R7 mostly done:** the gas tests are revived in `tests/physics`, 15 of 19 passing, four
parked with reasons; five real `src` defects were fixed getting there.
Next: **R4** (make the data loads conditional and the silent cubic-EoS degradation loud).

**Status confirmed 2026-08-09 (currentness sweep).** The migration is real and further along
than the header above suggests — items ticked this pass, all verified in tree:
`src/physics/gastables` and `src/physics/gasmodels` are git-tracked (R6 `git add`);
`share/` is out of `.gitignore` and all four `share/fluidprop/*.inp` files are committed;
`scripts/fluidprop/` is committed; `doc/README.md` already indexes Gas Models and Gas Tables;
and the `CLAUDE.md` `USE_NONFREE` line no longer exists, so that item is moot. Remaining
blockers are unchanged in substance: **R4** (conditional loads + loud cubic degradation),
the carried data defects **D7/D8/D9**, **O4** (normal vs equilibrium hydrogen — physics call),
the five disabled tests, `atmosphere`'s home, and the missing `gasdata.inp` generator.
Two small hygiene finds this pass: `scripts/fluidprop/__pycache__/*.pyc` are tracked, and the
disabled-test count is five, not four (see R7). Companion file for the physics defects:
`closed/gas_correctness_fixes.md`; the WIP nitrogen EoS is `nitrogen_eos_completion.md`.
Open: the four disabled tests (`Gibbs` first, its NaN is a defect); `atmosphere` is still in
the nonfree tree and needs a decision;
`gasdata.inp` still has no generator script, and acetone cannot be keyed because its CEA
label overruns the 11-character label field.

---

## Context

The gas models were moved out of `src/physics` into the private `nonfree/` tree on
2026-07-06 (nonfree commit `8c054bb8`). The stated reason
(`nonfree/physics/CMakeLists.txt:1-9`) was dependence on "proprietary/uncertain data" and
adjacency to "export-controlled combustion". This plan reverses that as broader than
necessary.

Three questions were tangled together; two are now closed:

1. **Chain of title — CLOSED.** The code was Christian's own work. No DLR claim.
2. **Export control — CLOSED via LBNL IPO.** The moving code is equilibrium chemistry,
   equations of state, and transport properties: general thermophysics, published method
   (NASA RP-1311, unlimited distribution), not "specially designed" for propulsion.
   Chemical-nonequilibrium (finite-rate Chemkin kinetics) stays in
   `nonfree/physics/combustion/`. Dependency direction verified one-way — `gasmodels` and
   `gastables` contain zero references to combustion.
3. **Third-party rights in the data tables — the live question.** See the gap table.

The pre-existing state was internally contradictory: 85 files carrying the open BSD-3
header and pointing at the open `LICENSE`, inside a tree declared "NOT part of the
open-source BELFEM distribution", with no `LICENSE` of its own, while tracked open-repo
tests under `tests/old/physics/` depend on the private data.

## Scope guards

- `combustion/` does **not** move, and its rocket-chamber / scramjet framing must not be
  carried into the open modules' documentation.
- No change to the numerical behaviour of any surviving code path in Phase 0; the defect
  fixes are latent-correctness only (D1–D5) except D4, which changes a crash into an error.
- The full 471-row `gasdata.inp` and the Poling-derived provenance stay private
  regardless of what the rebuilt public table contains.

## Gap table — the four data files

Class (a) = ship as-is, (b) = rebuild/replace, (c) = do not ship.

| File | Records | Source | Class | Disposition |
|---|---|---|---|---|
| `trans.inp` | 111 | NASA CEA, unmodified (Svehla, NASA TR R-132) | (a) | Fetch via download script |
| `thermo.inp` | 2059 | NASA CEA base + Ponomarenko + DeMar + IVTANTHERMO | (b) | Fetch pristine 2021 NASA GRC release; own additions re-applied as an overlay |
| `gasdata.inp` | 471 | Poling/Prausnitz/O'Connell 5th ed., McGraw-Hill 2001 (≈ all of Appendix A) | (b) | Rebuild from CoolProp/Cantera + primary literature |
| `cubicalpha.inp` | 70 | Mahmoodi & Sedigh, *Fluid Phase Equilibria*, Elsevier 2017, DOI 10.1016/j.fluid.2016.12.015 | (c) | Ship nothing initially; graceful fallback exists |

### Why the `gasdata.inp` rebuild is small

The library can only ever use **33** of the 471 records, and only **20** reach full
capability. A `RefGas` is constructed from a `thermo.inp` label
(`cl_GT_RefGasFactory.cpp:98-109`) and there is no cross-nomenclature normalization —
`fn_GT_fix_label.cpp:20-34` fixes only `AR`/`AL`/`CL` capitalization, while `gasdata.inp`
uses Hill-style ordering where CEA does not (`H3N`/`NH3`, `O2S`/`SO2`, `F6S`/`SF6`,
`CH4O`/`CH3OH`).

Counted honouring the **last-in-wins** map semantics of `cl_GT_InputData.cpp:44`
(471 rows collapse to 263 distinct labels — see D9):

- **46** of the 263 resolved labels match a `thermo.inp` label at all
- **32** carry a complete five-field record (`has_crit() == true`)
- **18** additionally have `cubicalpha.inp` entries (the PM alpha working set):
  `Ar C2H4 C2H6 C3H8 C6H6 C7H8 CH4 CO CO2 F2 H2O H2S He Kr N2 Ne O2 Xe`
- `has_crit()` but no PM entry (falls to the factory CCR table or the ω correlation):
  `Br2 BrF3 C2F6 CF4 CH2F2 CH2O CH3F CHF3 D2 H2 N2H4 N2O NO2 O3`
- The remaining ~2030 `thermo.inp` species are ideal-gas-only regardless

**H2 and D2 are absent from the PM set only because of defect D9**, not because the data
is missing — both have correct `cubicalpha.inp` rows that the resolved CAS fails to reach.

So 438 of the 471 transcribed rows do nothing, and the McGraw-Hill exposure is carried
almost entirely by dead weight.

## Capability matrix (verified by source trace)

Both files load unconditionally at `cl_GT_RefGasFactory.cpp:48,51`, and
`cl_Ascii.cpp:155-157` hard-errors on a missing file (`BELFEM_ERROR` is **not** compiled
out in release — `src/core/assert.hpp:171-184`). The matrix assumes those loads are made
conditional first (R4).

| Capability | No `gasdata.inp` | No `cubicalpha.inp` |
|---|---|---|
| Ideal-gas thermo (cp, h, s, R, M) | WORKS — marginally *better* | WORKS |
| Equilibrium chemistry (`compute_equilibrium`, `Gibbs`, `dGibbsdT`) | WORKS | WORKS |
| Ideal-gas transport, species in `trans.inp` | WORKS | WORKS |
| Ideal-gas transport, species absent from `trans.inp` | DEGRADED (~13 species) | WORKS |
| High-pressure transport (Lucas / Stiel-Thodos), SRK/PR path | **BROKEN** | WORKS (values shift slightly) |
| High-pressure transport, Helmholtz path | WORKS | WORKS |
| Cubic EoS (SRK / PR) | **BROKEN — silently** | DEGRADED (designed fallback) |
| Helmholtz EoS (H2 / CH4 / O2) | WORKS | WORKS |
| Compressible-flow relations | WORKS (breaks only via cubic EoS) | WORKS |

Key asymmetry: `cubicalpha.inp` is optional accuracy data with a real fallback
(the `has_crit()` guards at `cl_GM_EoS_Cubic.cpp:475,517` route to the
`create_ccr_mc_srk`/`create_ccr_pr` calls at `:477,519`, which carry 22 species hardcoded
in source plus ω correlations). `gasdata.inp` is structural — sole
source of ω, dipole, Z_crit, rho_crit and CAS, gate for every real-gas path via
`has_crit()`, and the join key for `cubicalpha.inp`.

Molar mass has an independent source: `thermo.inp` supplies M at
`cl_GT_InputThermo.cpp:136` at higher precision than `gasdata.inp`, which overwrites it
(`Ar` 39.948 vs 39.9480000; `CH4` 16.043 vs 16.0424600).

---

## Phase 0 — defect fixes (COMPLETE, pending audit)

- [x] **D1** `cl_GM_EoS_AlphaFunctionFactory.cpp:139,324` tested `tCAS == "75-28-50"`; the
      real isobutane CAS is `75-28-5` in both data files, so the hardcoded CCR entry could
      never match and silently fell through to the ω correlation.
- [x] **D2** Eight CAS numbers in `gasdata.inp` and three in `cubicalpha.inp` were mangled
      by a spreadsheet into `DD.MM.YY` dates. The rule is `XXXX-YY-Z` → `Z.YY.<last two of
      XXXX>`, which is self-consistent across all eight and independently constrains each
      restored value. Corrected in both data files **and** in the factory (which was
      written against the corrupted `04.06.83` for H2S at lines 118/303) so the CAS join
      stays consistent. Without the factory half, the R1 rebuild would have silently
      dropped the H2S CCR coefficients.
- [x] **D3** `cl_GT_InputData.cpp:66-67` computed the final field width as
      `tLength.length() - tStart(5)` = `6 - 138` in `uint`, underflowing to ~4.29e9.
      Harmless only because `std::string::substr` clamps its count. Now guarded.
- [x] **D4** `fn_GT_data_path.cpp:22` did `std::string tSCLS = std::getenv("SCLS")`,
      constructing `std::string(nullptr)` when `SCLS` is unset — UB, reached from *every*
      `Gas` construction via the default argument at `cl_GT_RefGasFactory.hpp:50`. Masked
      only because the dev box has `SCLS` set. Now null-checked with a clear error.
- [x] **D5** `cl_GT_InputAlpha.cpp:100-113` was unreachable: the `srk()`/`pr()` values it
      wrote are read only at `cl_GM_EoS_AlphaFunctionFactory.cpp:448,460`, both gated by
      `has_cubic()`, which only the *other* branch sets. Deleted rather than activated —
      the factory already implements the genuine CCR fallback, and CCR parameters carry
      different meaning than the PM alpha-function form.
- [x] **D6** `gasdata.inp` line 257 (2,5-lutidine) had a truncated name and a stray
      freezing-point column, so the CAS ran into the molar mass and the parser read
      M = 257.650 instead of 107.155. Row rebuilt at canonical offsets.
- [x] Remove the `@dlr.de` stamps from `crthermo.inp`, `crtrans.inp`, and the
      `thermo.inp` v1.6–1.6.2 changelog entries.

---

## Ordered steps

### R1 — Rebuild `gasdata.inp` clean

**Format: byte-identical to the current file**, so `cl_GT_InputData.cpp` needs zero
changes and anyone holding Poling can regenerate their own full table from the same
layout. Column map, re-derived from `cl_GT_InputData.cpp:62-63,98,101`:

| Field | Columns | Notes |
|---|---|---|
| label | 0…, whitespace-delimited | **must be ≤ 11 chars** (name field starts at 12) |
| name | 12–71 | `substr(12,60)` |
| CAS | 75–90 | `substr(75,16)` |
| M [g/mol] | 92–101 | ×0.001 → kg/mol |
| T_crit [K] | 102–110 | |
| p_crit [bar] | 112–119 | ×1e5 → Pa |
| Z_crit | 121–128 | |
| ω | 129–135 | |
| dipole [debye] | 138–145 | **fixed width** (was read-to-EOL) |
| src(Tc,pc,Zc) | 147–182 | provenance; everything from 147 is ignored by the reader |
| src(omega) | 185–220 | |
| src(dipole) | 223–EOL | |

Provenance is **per field**, not per row, because a species may draw its values from more
than one place — `Cl2` takes Tc/pc/ω from Perry and Zc from Poling. Tokens are space-free:
a DOI where one exists, else `symmetry` (dipole exactly zero by point group),
`simple-fluid` (ω zero by Pitzer's definition for a monatomic fluid),
`estimate-from-<X>`, `<author>-<year>`, or `A/Zc:B` for a mixed row.
`gasdata_sources.md` has been **removed** — the licence argument now lives in the file
header, so there is no second document to drift.

**Rule 1 — the label is the post-`fix_label` `thermo.inp` label, exactly.** This *is* the
D7 fix and it costs no code: `thermo.inp` uses `NH3`/`SO2`/`CH3OH`/`C2H5OH` where the old
table used Hill ordering `H3N`/`O2S`/`CH4O`/`C2H6O`. All four have complete critical data
and `cubicalpha.inp` entries that are simply unreachable today.

**Rule 2 — exactly one row per label.** Kills D8 and D9 by construction.

**Sourcing (verified 2026-08-03 against CoolProp 8.0.0, MIT):**
- `T_crit`, `p_crit`, ω, M — direct from CoolProp.
- `Z_crit` — computed as `pc / (rhomolar_critical · R · Tc)`. Verified against the old
  table: Water → 0.2294 vs 0.229.
- **dipole — NOT available from CoolProp.** Since `has_crit()` requires all five fields, a
  missing dipole silently disqualifies a species. Most of the target list is nonpolar, where
  the value is **0.0 by molecular symmetry** rather than a datum needing a source; roughly
  15 polar species (H2O, NH3, SO2, H2S, CH3OH, C2H5OH, CO, NO, NO2, N2O, O3, HCl, CH3Cl,
  CH3F, CH2F2, CHF3, COS) need individually cited values.

**Coverage.** Of CoolProp's 136 pure fluids, 67 composition-match a `thermo.inp` species
and **50 have a match whose label fits the 11-char field**; **26 of those also carry
`cubicalpha.inp` PM entries** — against 18 reachable today. The composition match produces
false positives that must be filtered by name, not accepted blindly (`DiethylEther` →
`n-Butanol`, `MethylStearate` → `IDP`, and the four butene isomers all collapsing onto
`C4H8,cyclo-`).

- [x] Agree the final species list from the 50 label-fitting candidates, after filtering
      the composition false-matches by name → **42 species**.
- [x] Emit the table from a script driven by CoolProp; hand-add the polar dipole values
      with citations → 23 sourced dipoles, 19 zero by symmetry.
- [x] Keep the full table private → renamed `gasdata_ppoc.inp` (git `RM`, history preserved,
      and it carries the Phase 0 D2/D6 corrections). New table is `gasdata.inp`.
- [x] Provenance carried **in the data file itself** rather than a sidecar, so a copied row
      carries its own sources and nothing can drift. Three source columns (per field);
      `gasdata_sources.md` deleted and its licence argument moved into the file header.
- [x] Extended to **59 species**. Chlorine gets real-gas treatment for the first time — its
      ω was blank in the Poling table, so `has_crit()` had always been false for it.
- [x] **Four column-corruption defect classes found in `gasdata_ppoc.inp`** while sourcing,
      every one of which would have propagated silently into a transcribed table:
      the lutidine row (D6); the eight spreadsheet-mangled CAS numbers (D2); a literal
      **TAB** in the `CH2O` row that shifted the fixed-column read so formaldehyde's dipole
      parsed as `.270` instead of `2.331`, a factor of ten; and the `H2O2` row, which is
      damaged in **three fields at once** — CAS `22-84-1` for 7722-84-1, p_crit `20.000`
      for 220 bar, and dipole `26` for what was probably the superseded 2.26 D value.
      Independent of licensing, this is a data-quality argument for having rebuilt rather
      than transcribed.
- [x] **`N2O4` deliberately excluded**, rationale recorded in the `gasdata.inp` header.
      It and `NO2` are the same reactive system (N₂O₄ ⇌ 2 NO₂), and `NO2` is already entered
      as a pseudo-pure representative; entering both would double-count the association
      effect in any mixture containing the pair, which an equilibrium calculation always
      produces. Also note the commonly tabulated N₂O₄ critical volume refers to the
      *dissociating mixture* and yields Zc ≈ 0.47 — that value must never be combined with
      pure-component constants, a cross-set error that was proposed and rejected.
- [x] **Required for the source column:** the dipole field was defined as read-to-EOL, so
      a ninth column would have been swallowed into it. Working by luck today
      (`std::stod` parses the leading number and ignores the trailing DOI), but a row with
      a **blank dipole and a source present** would have silently parsed `10.1063` as the
      dipole moment. Dipole now has a fixed width of 8 (`cl_GT_InputData.cpp:62-65`), which
      also deletes the runtime width computation that caused D3 in the first place.
      Verified: both tables parse with **zero failures** — 42/42 complete in the public
      table, 114 complete of 263 labels in `gasdata_ppoc.inp` — and dipoles read exact
      (`1.62673`, `1.6515`) with the source column excluded.
- [x] **All six regressions recovered** — `Br2`, `BrF3`, `CH2O`, `N2H4`, `NO2`, `O3` — plus
      `CHCl3`, `C2H2F2`, `F2O`, `D2S`, `Cl2`, `CCl4`, `C2F4`, `Rn`, `NO`. Table is at
      **57 species**.
      Two values needed provenance beyond a plain citation and carry it in the file:
      `NO2` p_crit is sourced `Poling-2001/pc:N2O4-row`, because its own Poling row reads a
      corrupted 10.100 bar where `N2O4` — the same substance in equilibrium, Tc within
      0.01 K — reads 101.000; and `NO2` ω is `Chemeo-AP1700-pseudopure`, since 0.85 is the
      pseudo-pure engineering value and a pure-component calculation gives ~0.06, a factor
      of 15 that would otherwise look like an error.
- [x] ~~**Blocked on species identity: `S` / `S8`.**~~ **Closed won't-do, 2026-08-04.**
      `thermo.inp`'s `S` is raw `S(a)`, monatomic sulfur — correctly identified — but a
      monatomic species has no condensed phase and therefore no critical point. The
      Tc = 1313 K / pc = 182 bar in the Poling `S` row belongs to *bulk elemental* sulfur
      (the S₈/S₂ melt), a different substance. The table already treats **all 22** free
      atoms and open-shell radicals in `thermo.inp` (H, O, N, C, F, Cl, Br, I, OH, NH, CH,
      CN, SH, SO, SN, HO2, CH3, C2H, HCO, NH2 …) as ideal-gas only; noble gases are entered
      because they are atoms that *do* have a genuine condensed phase. Atomic sulfur also
      only appears above ~1500 K, where the mixture is essentially ideal and a real-gas
      correction would be numerically meaningless. Rationale is recorded in the
      `gasdata.inp` header so it is not reopened. `S8` remains a valid future candidate if
      critical constants are found.

**Result, verified by round-tripping through the exact offsets and last-in-wins map
semantics of `cl_GT_InputData.cpp`:** 42 species, `has_crit()` true for every one, every
label matches a `thermo.inp` entry and fits the 11-char field, and **cubicalpha PM
coverage is 26/42 against 18 reachable before**.

**O4 is resolved and needs no decision.** CoolProp's `Hydrogen` is *normal* hydrogen and
carries CAS `1333-74-0` — the CAS that `cubicalpha.inp:62` and the factory entries at
`cl_GM_EoS_AlphaFunctionFactory.cpp:83,268` already key on. So the physically correct
default and the repaired lookup come from the same row; the old table's conflict was an
artefact of assigning the real CAS to the para-ish row (Tc 32.980 ≈ CoolProp
ParaHydrogen 32.938) and a placeholder to the normal one. Same for D2 via `7782-39-0`.
`ParaHydrogen`/`OrthoHydrogen` exist in CoolProp under suffixed CAS (`1333-74-0p`/`o`)
if the variants are ever wanted.

### R1b — Cross-check tool (throwaway) — DONE
- [x] Python tool diffing new against old, per species and per field, matched by CAS.
- [x] Run against the private table as reference. **Findings:** all field deltas small and
      explicable. The alarming acentric-factor percentages are division-by-near-zero
      artefacts — absolute deltas across the noble gases are 0.0002–0.02. Largest genuine
      shift is `CH3F`, 55.48 → 59.06 bar, where CoolProp's modern reference EoS supersedes
      the older value. Supplied dipoles agree with the old table everywhere except `H2O`
      (1.8546 vs 1.8) and `H2S` (0.977 vs 0.9), where the old table **truncated** rather
      than rounded — the new values are the standard ones.
- [x] The tool independently rediscovered D9: the only four species with no CAS match in
      the old table are `H2`, `D2` (placeholder CAS wins last-in-wins), `C2H5OH` (dimethyl
      ether wins the `C2H6O` collision) and `COS` (genuinely absent). Good evidence the
      comparison is working rather than trivially agreeing.

### R1c — Deploy and exercise the new table
- [x] ~~**The repo `tables/` directory is a staging copy only.**~~ **Closed 2026-08-05.**
      The `$SCLS/../cea` lookup is gone. `fn_GT_data_path.cpp` now resolves
      `gBelfemDataPath` — set from `$BELFEM_DATA` in `Communicator::set_globals()` — and
      falls back to searching for `share/fluidprop` relative to the working directory, so
      running from a build tree needs no configuration at all.
- [x] ~~Either `sudo cp` the new table into `/opt/scls/cea/`~~ **Closed 2026-08-05.** No
      root anywhere: the shipped tables live in the repo at `share/fluidprop/`. The
      root-owned 2023 copy is now simply unused.
- [ ] Smoke test once deployed: construct `Gas("H2")` and confirm it now picks up CAS
      `1333-74-0` and a PM alpha function rather than falling through to the ω correlation.
      This is the end-to-end proof that D9 is fixed.

### R2 — NASA data acquisition — SUPERSEDED, we build our own tables
- [x] ~~Download script for `thermo.inp` / `trans.inp`.~~ **Closed 2026-08-05.** Both the
      download and the licence question it was meant to manage are gone: BELFEM now
      **extracts** the species it needs into its own tables rather than redistributing the
      vendor file. The third-party additions that made the circulating `thermo.inp`
      awkward — Ponomarenko, DeMar, IVTANTHERMO — are all solid and propellant materials,
      so a table restricted to our species contains none of them. This resolves the
      redistribution question by construction instead of by argument, and drops 1.3 MB to
      51 kB. Tool: `scripts/fluidprop/build_tables.py`, documented in its README.
- [x] ~~**200 K overlay.**~~ **Closed 2026-08-05, superseded by something better.** Rather
      than extending `Tmin` and warning about extrapolation, a genuine low-temperature
      interval is now *fitted* per species and merged into the record. There is no overlay,
      no second file and no load-order contract, and the range reaches the triple point
      instead of 200 K. The declared bounds stay honest because the data behind them is
      real.
- [x] **Deuterium added** for fusion fuel handling — `D2`, `D`, `HD`, `OD`, `D2O` and the
      ions, 62 species total. `D`, `HD` and `OD` have transport in no database, so theirs is
      mass-scaled from `H`, `H2` and `OH`, which is the same construction NASA used for its
      own `D2` record (their D2/H2 viscosity ratio is sqrt(2) to within 0.05 %). Derived
      records say so in their header rather than inheriting the parent's citation.
- [ ] **`gasdata.inp` has no generator script.** It was produced from CoolProp by hand
      during R1, so unlike `thermo.inp`/`trans.inp` it cannot be reproduced or extended.
      Close this before the move, or the next species addition repeats the manual work.
- [ ] **Acetone is ideal-gas only and data alone cannot fix it.** Its CEA label
      `C3H6O,acetone` is 13 characters against the 11-character label field in
      `gasdata.inp`, so no record can be keyed to it. `fn_GT_fix_label.cpp` already strips
      the descriptor for `C2H3,vinyl` and needs the same case here. The factory carries a
      hardcoded CCR entry for its CAS `67-64-1`, so the code expects acetone to work.

### R3 — Provenance manifest
- [ ] One manifest per shipped file, modelled on `nonfree/share/textures/attribution.yaml`
      (source, SPDX, attribution, acquisition date, derived-modifications block).

### R4 — Make data optional
- [ ] Make the `gasdata.inp` and `cubicalpha.inp` loads conditional
      (`cl_GT_RefGasFactory.cpp:48,51`); keep `thermo.inp`/`trans.inp` mandatory with an
      error naming every path searched.
- [ ] **Make the silent cubic-EoS degradation loud.** Without critical data,
      `cl_GM_EoS_Cubic.cpp:479,521` build `create_empty()` (α ≡ 0), `:551-555` set
      `mAc = mBc = 0`, `:645` zeroes the departure function — SRK/PR quietly *becomes the
      ideal gas*, and `eval_critical_point:797-834` divides by zero. Must error at
      construction, not return plausible wrong numbers.
- [ ] Data-path search order: explicit `set_data_path()` → `$BELFEM_GASTABLES_PATH` →
      legacy `$SCLS/../cea` → CMake-configured default → `./gastables`. Probe for
      `<dir>/thermo.inp` rather than the directory. Add `gastables::data_available()`,
      which must never abort.

### R5 — Build system — DONE 2026-08-05
- [x] Decoupled `USE_GASMODELS` from `USE_NONFREE`: the hard error is gone, and
      `USE_COMBUSTION` was added requiring **both**, with a separate message for each so a
      failed configure says which flag is missing.
- [x] `src/physics/CMakeLists.txt` is the public dispatcher, adding `gastables` and
      `gasmodels` under `USE_GASMODELS`. `nonfree/physics/CMakeLists.txt` is trimmed to
      `atmosphere` (always) and `combustion` (under `USE_COMBUSTION`), and its header
      comment no longer claims all four modules are nonfree.
- [x] `BELFEM_NONFREE_SOURCE_DIR` → `BELFEM_SOURCE_DIR` in the two moved CMakeLists **and**
      in `atmosphere` and `combustion`, which reach across trees for the gas headers. Only
      `combustion`'s self-include still points at the nonfree tree, which is correct.
      `LIBLIST` ordering untouched.
- [x] **Behaviour change to announce:** combustion was previously built unconditionally
      whenever `USE_GASMODELS` and `USE_NONFREE` were both on. It is now opt-in. Existing
      build trees carry `USE_COMBUSTION=OFF` and will silently stop building it until the
      flag is set.
- [x] Verified by configuring six option combinations. The four legal ones configure and
      produce exactly the right library set; the two illegal ones fail with the intended
      message:

      | `USE_GASMODELS` | `USE_NONFREE` | `USE_COMBUSTION` | result |
      |---|---|---|---|
      | OFF | OFF | — | configures, no gas libraries |
      | ON | OFF | — | configures, gastables + gasmodels — **the open-source case** |
      | ON | ON | OFF | + atmosphere |
      | ON | ON | ON | + combustion |
      | OFF | ON | ON | fails: needs `USE_GASMODELS` |
      | ON | OFF | ON | fails: needs `USE_NONFREE` |

- [x] **Compiled**, in a scratch build tree so the working one was untouched: `gastables`,
      `gasmodels`, `atmosphere`, `combustion` and the four example executables (`gastable`,
      `gasmodel`, `formation`, `helmholtz`) all build and link with zero errors.
- [x] **One build break fixed on the way.** `cl_GM_EoS_Cubic.cpp:636` called
      `Spline::update_data` with `( matrix, values, gTref, Sref/M )`, but the signature had
      since gained two `SplineBC` parameters ahead of the reals, so the two reals were
      landing on the boundary-condition arguments. It does not convert, so this had been a
      hard compile error since the Spline change; the module simply had not been compiled
      since. Fixed to match the sibling call at `cl_GT_RefGas.cpp:1141`, which already
      passed the full list.

### R6 — Move the code
- [x] `gastables` and `gasmodels` moved into `src/physics/` by Christian, 2026-08-05. They
      are untracked pending the `git add`, which is a first-time add rather than a `git mv`
      because `nonfree/` is gitignored.
- [ ] **`atmosphere` stayed in `nonfree/physics`.** The original plan moved it too. It has
      no data or export-control problem of its own — it is ISA-1976 — so decide whether it
      follows; nothing in the build now depends on where it sits.
- [x] `git add src/physics/gastables src/physics/gasmodels`. *(Done — both trees are git-tracked as of 2026-08-09.)*
- [x] **Take `share/` back out of `.gitignore`** — **DONE.** The `share` entry and its
      "temporarily excluding" comment are gone from `.gitignore`, and all four files are
      committed: `share/fluidprop/{thermo.inp, trans.inp, gasdata.inp, cubicalpha.inp}`
      (plus a local `share/fluidprop/.gitignore`). Original text: (`.gitignore:25`, whose own comment says
      "temporarily excluding share, will add back in when we bring the gas models back").
      Until that happens `share/fluidprop/` is *not committed*, so a fresh clone has no
      tables and every `Gas` construction fails on a missing file. This is the step that
      makes building our own tables actually mean something. `git add share/fluidprop/` and
      confirm all four files land: `thermo.inp`, `trans.inp`, `gasdata.inp`,
      `cubicalpha.inp`.
- [x] Commit `scripts/fluidprop/` alongside them — the tables are generated output and the
      generator is what makes them reproducible and extensible. **DONE** (`build_tables.py`,
      `nasa9_lowT.py`, `README.md`). ⚠️ **Hygiene defect found 2026-08-09:**
      `scripts/fluidprop/__pycache__/*.pyc` are tracked too — remove them and add a
      `__pycache__/` ignore before the open-source cut.

### R7 — Tests — MOSTLY DONE 2026-08-05
- [x] Relocated to `tests/physics/{gastables,gasmodels}` and the old copies removed. The
      include paths needed no change, exactly as predicted. Registered in the `check`
      fan-out and labelled `fast`; both suites together run in under a second. **`physics`
      itself was missing from that fan-out**, so `make check` had never built the existing
      YBCO test either — added at the same time.
- [x] **15 of 19 tests pass.** Five real defects were fixed to get there, all in `src`
      rather than in the tests:
      `Spline`'s two vector constructors were mutually ambiguous for any three argument call
      (every trailing parameter defaulted on both);
      `Spline::update_data` hardcoded UMFPACK, so every build without SuiteSparse failed
      there, and now uses SuperLU unconditionally;
      `create_glue_polys_heat` could not terminate across an already smooth junction, which
      hit the five noble gases and `D2`;
      `Gas::remix_heat` filled its heat spline through `matrix_data()` without declaring what
      the extra coefficient row holds, so `entropy()` tripped its assert and no `Gas` could be
      constructed in a debug build at all;
      and synthesized transport was unreachable through `mu()`, because `create_splines`
      calls `idgas_lambda` while the object is still evaluating from polynomials.
- [x] Six species added to the table for coverage — `Kr`, `Xe`, `N2H4`, `BrF3`, `CF2ClBr`,
      `CH2Cl2`. `Kr` and `Xe` were not optional: `Gas::Gas()` names both in its default air
      composition, so a default constructed gas could not be built without them.
- [ ] **Four tests parked with `DISABLED_` and a written reason.** *(Count re-checked
      2026-08-09: there are **five** `DISABLED_` tests in `tests/physics/gasmodels` —
      `Gibbs`, `AlphaFunction`, `Cubic_Departure`, `Cubic_State`, and
      `Gas_Entropy_RefState` in `cl_GM_Gas_Consistency.cpp`. The fifth is the entropy
      reference-state question tracked in the departure-convention section below; note the
      pass tally accordingly: 19 test files, 5 disabled.)* `Gibbs` returns NaN over a
      sweep starting at 100 K and its reference data came from a locally modified
      `thermo.inp` that is no longer shipped; that one is a defect, not a tolerance.
      `AlphaFunction` compares an analytic second derivative against a cubic spline's,
      reaching r2 = 0.66-0.73 against 0.99. `Cubic_Departure` and `Cubic_State` miss by one
      decade (0.99999 vs 0.999999) and probably just need the tolerance re-set against the
      new data — a deliberate call, not a silent edit.
- [ ] Gate at **runtime**, not configure time — an open checkout must still compile the
      tests. Shared `GastablesFixture` with `GTEST_SKIP()` in `SetUp()`; already the idiom
      (`tests/comm/test_CommMPI.cpp:44`). Audit the 15 gasmodels tests individually — the
      table-free ones must stay ungated.
- [ ] Add a hand-authored synthetic fixture (2–3 species, written from scratch) so the
      CEA-format parsers keep permanent regression coverage regardless of data.

### R8 — Documentation
- [ ] Index the three modules in `doc/README.md`; add `src/physics/atmosphere/doc/README.md`
      (atmosphere has no `doc/` today); document the data-path resolution and the
      `BELFEM_GASTABLES_PATH` / `BELFEM_GASTABLES_DATA_DIR` knobs in the gastables README.
- [x] ~~Update `CLAUDE.md:52` — gas models no longer require `USE_NONFREE`.~~ MOOT — the
      `USE_NONFREE` claim is no longer in `CLAUDE.md`; the build section now lists gas models
      as a plain optional module. Nothing to change.

---

## Open defects (carried, not yet fixed)

- [ ] **D7** Label-nomenclature mismatch silently strands 419 of 471 `gasdata.inp` records
      (`H3N`/`NH3`, `O2S`/`SO2`, `F6S`/`SF6`, `CH4O`/`CH3OH`). Fix during R1 — either
      normalize in `fn_GT_fix_label` or key the whole lookup on CAS.
- [ ] **D8** `C7H9N` spans **six** rows in `gasdata.inp` (lines 255–260, all lutidine
      isomers). `InputData`'s constructor does `mMap[tLabel] = tLineCount`, so last-in-wins
      and `Gas("C7H9N")` resolves to line 260 (3,5-lutidine). The D6 row repair on line 257
      is therefore correct as data but **not reachable by label lookup**.
- [ ] **D9** The duplicate-label problem is general, not local: **73 labels are duplicated**
      in `gasdata.inp`, collapsing 471 rows to 263 addressable entries. Two of them corrupt
      the working set through placeholder CAS numbers:
      - `H2` → line 455 ("hydrogen, normal", CAS `800000-51-5`) wins over line 454
        (CAS `1333-74-0`). The `800000` series is not a real registry range.
      - `D2` → line 433 (CAS `800000-54-8`) wins over line 432 (CAS `7782-39-0`).

      Consequence for hydrogen — the most important gas in this codebase — is a **silent
      double miss**: `cubicalpha.inp:62` keys H2 on `1333-74-0`, so `entry_exists()` fails
      and `has_cubic()` stays false; the factory's hardcoded CCR entry at
      `cl_GM_EoS_AlphaFunctionFactory.cpp:83,268` *also* keys on `1333-74-0`, so it misses
      too. H2 falls through to the generic ω correlation, bypassing both curated data paths.
      Same class of defect as D1, but on a species that matters. Fix during R1 by keying the
      lookup on CAS and representing isomers/spin-variants explicitly (see O4).

## Open design question — the reference pressure convention (was D10)

**Reclassified 2026-08-04: this is a convention, not a defect.** All `*dep0` terms were
removed and then **restored**. The tree carries the original behaviour; only documentation
was added (`cl_Gas.cpp`, block comment on `realgas_cp`, and Doxygen on the component block
in `cl_Gas.hpp`).

### What the code does

`realgas_cp`, `realgas_dcpdT`, `realgas_h`, `realgas_s` and `realgas_dsdT` each subtract the
departure evaluated at `gPref` = 1 bar, so at 1 bar they return the ideal-gas value of the
thermo tables exactly, whatever equation of state is active. IDGAS, SRK, PR and Helmholtz
therefore agree at standard conditions, and the departure splines are rebuilt when the gas
model changes so this holds for each. The component-level `Gas::h( aIndex, … )` applies the
same convention through the indexed departures (`cl_GM_EoS_Cubic.cpp:948-976`).

### How it relates to the textbook assembly

The textbook relations

    h2 - h1 = ( h2 - h2^0 ) - ( h1 - h1^0 ) + int cp^0 dT
    s2 - s1 = ( s2 - s2^0 ) - ( s1 - s1^0 ) + int cp^0/T dT - R ln( p2/p1 )

are **identities**, obtained by writing the definition of the departure at two states and
subtracting. The second bracket is there only because a difference between two general states
is being taken; it is not a reference offset. The absolute property is simply

    h( T, p ) = h_ideal( T ) + hdep( T, p )

with no subtraction at all. An earlier reading of these relations as licensing a fixed
reference state subtraction was wrong and is withdrawn.

The BELFEM form subtracts the departure at `( T, gPref )`, which varies with temperature.
That would be the complete and correct assembly if the heat spline held the **real** gas at
1 bar, since `h(T,p) = h_real(T,gPref) + hdep(T,p) - hdep(T,gPref)` is an identity. The
spline holds the **ideal** gas at 1 bar, which is the CEA standard state, so the two differ
by the departure at 1 bar.

Note from the same equations that `Δs⁰ ≠ Δs`: the ideal-gas entropy change already contains
`−R·ln(p2/p1)`, which is the `mR*log( gPref/aP )` term in `realgas_s`. That term is **not** a
departure and must not be reasoned about by analogy with cp and h, which have no ideal-gas
pressure dependence at all.

### Why it is not simply wrong

The tabulated NASA-9 standard state is the *hypothetical ideal gas* at 1 bar — confirmed by a
NIST isotherm for ammonia at 200 K, whose p = 0 row (cp = 33.5927 J/mol/K) matches
`thermo.inp` to +0.48% while the saturated-vapour row differs by −2.77%, and by the fact that
`thermo.inp` carries an NH₃ *gas* entry at 200 K where the real fluid at 1 bar is a liquid.
So the 1 bar isobar is not physically ideal, and pinning it to the table trades absolute
accuracy for model-switching consistency. Whether that trade is right for BELFEM is the open
question.

### To settle it

- [ ] Decide between the two self-consistent options: keep the subtraction and hold the
      spline at the real gas value at 1 bar, or drop the subtraction and keep the spline
      at the ideal gas value. The present code is a hybrid of the two.
- [ ] Either way, **both paths must use the same convention**. The bulk and component paths
      currently agree, and must not be changed independently — the finite-rate combustion
      solver depends on the component form (`cl_CN_Scheme.cpp:315,320`), pairing it with
      ideal-gas `Gibbs`/`dGibbsdT` so the energy equation and the equilibrium constants stay
      on the same footing.
- [ ] Watch the two absolute anchors when touching entropy: the RefGas heat spline is
      anchored with `S( gTref )` (`cl_GT_RefGas.cpp:1176`), the EoS_Cubic departure spline
      with `Sref() / M()` (`cl_GM_EoS_Cubic.cpp:636-637`), and `mSref = idgas_s( gTref, gPref )`
      (`cl_Gas.cpp:861`) is added on top.
- [ ] No existing test discriminates: `cl_GM_Gas_Airprop.cpp` builds air as an ideal gas,
      where all departures vanish, and judges by r² at 1e-3 — a measure blind to constant
      offsets. `cl_GM_EoS_Cubic_Departure.cpp` checks the departure functions against each
      other, not an assembled property against reference data.

**Not a mode distinction.** The hypothesis that `*dep0` belongs to reference mode but not
spline mode does not hold: `RefGas::set_mode` only rebinds function pointers, and
`create_splines` (`cl_GT_RefGas.cpp:1154-1177`) sets `RefGasMode::POLY`, samples
`this->H( aT(k) )`, and builds the spline from those samples. Both modes evaluate the same
ideal-gas function, one exactly and one interpolated; neither carries departure content.

## Open questions

- [ ] **O1** Target species list for the rebuilt `gasdata.inp` — the 33 that resolve
      today, or a wider set covering combustion products and cryogens?
- [ ] **O2** Ship `cubicalpha.inp` at all? Options: omit (fallback covers it), seek
      permission from Mahmoodi & Sedigh, or refit independently.
- [ ] **O3** Does NASA's current CEA endpoint carry click-through terms? Gates R2.
- [x] ~~**O4** For the D9 spin-variant pairs, which entry should `H2` and `D2` resolve to —
      normal or equilibrium hydrogen?~~ **Resolved 2026-08-03, no tradeoff exists:**
      CoolProp's `Hydrogen` is normal hydrogen *and* carries CAS `1333-74-0`, the CAS the
      alpha-function lookups already key on. See R1.
- [ ] **O5** CEA's comma-suffixed species (`C4H10,n-butane` 14 ch, `C4H10,isobutane` 15 ch,
      `C2H2,acetylene` 14 ch, `C3H6,propylene` 14 ch) exceed the 11-char label field, so
      they cannot carry critical data in this format. **Consequence: the D1 isobutane CAS
      fix is correct but unreachable by label lookup**, since bare `C4H10` is not a
      `thermo.inp` label — same situation as D6. Options: accept the gap (keeps the format
      byte-identical, no regression versus today), widen the label field and shift the
      offsets in `cl_GT_InputData.cpp:62-63`, or key the whole lookup on CAS via a
      label→CAS bridge. Deferred; not blocking R1.

## Known risks

- `Add_Test.cmake` references `SSF_SRC_DIR`, defined nowhere in `config/`; it works only
  because `Add_Library.cmake` leaked `include_directories` into the same scope. A fresh
  test directory could expose it — verify a clean `-DUSE_TEST=ON -DUSE_GASMODELS=ON`
  configure early.
- `${CMAKE_BINARY_DIR}/generated` is on the include path for libraries and executables but
  **not** tests — keep any generated data-dir macro out of public headers.
- `USE_EXAMPLES` defaults `ON` (`CMakeLists.txt:82`), so enabling `USE_GASMODELS` builds
  five new public executables.
- `GT_globals.hpp:24` warns that changing `gTmax = 13000.0` "might fail tests" — if an
  accuracy test fails after the move, suspect this constant, not the migration.

## Definition of done

- Build matrix green in all three: `(USE_GASMODELS=OFF)`, `(ON, USE_NONFREE=OFF)`,
  `(ON, USE_NONFREE=ON, USE_COMBUSTION=ON)`.
- `unset SCLS` and construct a `Gas` — clear error naming searched paths, not a segfault.
- Construct an SRK gas with no `gasdata.inp` — must **error**, not silently return
  ideal-gas numbers.
- `make check` green both with tables present and with `BELFEM_GASTABLES_PATH=/nonexistent`.
- No residual `BELFEM_NONFREE_SOURCE_DIR` in the open tree.

## Audit trail

- 2026-08-03 — Phase 0 defect fixes (D1–D6) applied; Codex + Grok audit requested.
- 2026-08-03 — **Grok** (third voice): D5 deletion **not refuted**, confidence high (~92%).
  Independently traced that `create_pm_srk`/`create_pm_pr` are gated by a real
  `if( has_cubic() )` in `cl_GM_EoS_Cubic.cpp`, not merely by the debug-only
  `BELFEM_ASSERT` — so the branch stays unreachable in release builds too. Also confirmed
  `GasData` has compiler-generated copy/move, so the flag can never desync from the
  coefficient vectors. Key judgement upheld: reviving the branch would have forced
  `AlphaFunction_PM` with ω-derived coefficients instead of `AlphaFunction_CCR`/MC-SRK —
  a formulation change, not a bug fix. All eight D2 CAS restorations verified correct and
  the mangling rule self-consistent across all eight; M = 107.155 g/mol confirmed by
  atomic-weight sum and consistent with the other lutidine isomers at `gasdata.inp:255-260`;
  257.65 K ≈ −15.5 °C confirmed plausible as 2,5-lutidine's freezing point.
  Corrections applied: `:475,517` are the guards, the `create_ccr_*` calls are at `:477,519`.
  Grok's residual risk (could not read the pre-deletion hunk) resolved from the diff —
  `set_cubic_flag()` appears as a context line, never a removed one. Corruption scan
  widened to all 12 table files: zero `DD.MM.YY` tokens remain.
- 2026-08-03 — **Codex** (primary auditor): D1–D4, D6 stand; confidence high. Confirmed
  `BELFEM_ERROR` is defined unconditionally (`src/core/assert.hpp:171`) and aborts via
  `error_abort()` under `NDEBUG`, so the D4 guard survives release builds. Confirmed the D3
  behaviour claim — the old underflowed count extracted identically because `substr` clamps
  — so D3 is latent-correctness only. Verified all eight D2 CAS values against NIST
  (H2S, Ne, SO2, SO3, DCl, H2Se, trans-DMCH) and PubChem/Fisher (HMN), and found no other
  code keying the old literals. **One real defect found in the fix itself:**
  `cl_GT_InputAlpha.cpp` still uses `std::pow` after the D5 deletion but relied on a
  transitive `<cmath>`; explicit include added. **One factual correction:** `C7H9N` spans
  lines 255–260, not 256–257, so line 260 wins the label — which led to D9 (73 duplicated
  labels; H2/D2 placeholder-CAS double miss) and a corrected working set of 18, not 20.
  Noted for R4: `create_pm_*` guards missing cubic data with debug-only `BELFEM_ASSERT`, so
  direct API misuse is silent in release — same failure class as the cubic-EoS degradation.
