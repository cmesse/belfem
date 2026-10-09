# DR-193: Extend and Source the Quench Deck's Property Tables

**Date:** 2026-10-08
**Purpose:** `examples/tape_quench_usermat` reaches 1770 K at its hot spot while its user-material
tables end at 199 K (`k_Hast`), 282 K (`cp_HTS`), 300 K (`cp_Cu`, `k_HTS`), 303 K (`cp_Hast`) and
390 K (`rho_Hast`), with the endpoint held above; none of the thirteen tables carries a provenance
line. Mechanism: replace or extend each short table with published data to at least 1300 K and give
every table a header that names its source with a DOI beside each tabulated range, so the hot-spot
numbers are reported against the model's real validity envelope instead of the 199 K one.
**Module:** `examples/tape_quench_usermat/data`, `examples/tape_quench_usermat/README.md`; no C++
**AIs involved:** Claude (plan, literature routing, table assembly); Codex (number-by-number audit of the tables against the cited sources, terra/high); no implementation code, so no Opus
**Status:** PLAN — written 2026-10-08, not yet started; a literature session with web access

> **Scope guards.** In scope: the six short tables and their provenance; a provenance header for
> the seven others (`cp_Ag`, `cp_Buffer`, `k_Ag`, `k_Buffer`, `k_Cu`, `rho_Ag`, `I_vs_t`); the README's
> table section; a rerun of the T5 deck to its first 50 ms to show the envelope no longer binds.
> Out of scope: the plugin code (`src/matlib.cpp` reads two-column files and holds the endpoint —
> unchanged), the builtin materials database (`src/physics/materials`, HTS-scoped by design), any
> change to the deck's physics.

## 1. Current Behaviour and How It Fails

| Table | rows | range (K) | needed |
|---|---|---|---|
| `k_Hast.txt` | 47 | 3.7–199 | to ≥ 1300 |
| `cp_HTS.txt` | 34 | 10–282 | to ≥ 1300 |
| `cp_Cu.txt` | 13 | 10–300 | to ≥ 1300 (13 rows is also too coarse below 300) |
| `k_HTS.txt` | 13 | 10–300 | to ≥ 1300 |
| `cp_Hast.txt` | 43 | 2.2–303 | to ≥ 1300 |
| `rho_Hast.txt` | 39 | 14.5–390 | to ≥ 1300 |
| `k_Cu.txt` | 31 | 1–1200 | header only |
| `cp_Ag`, `k_Ag`, `rho_Ag`, `cp_Buffer`, `k_Buffer` | — | to 500–1000 | header; extend if the source reaches |

Above its last row every table holds the endpoint (`usermat::Table`), so the 1770 K result is the
model outside its data (DR-105 residue 3, ruled O3 on 2026-10-03: report against the envelope).

**Bottom line:** the deck's thermal model is unsourced and stops at 199 K; the fix is data, with a
DOI beside every number (Christian's standing rule).

## 2. Sources to route through (candidates; each number is cited at the table, not here)

- **Copper, cp 300–1300 K:** NIST-JANAF thermochemical tables (Cu, cr); below 300 K the NIST
  Cryogenic Materials Properties database (OFHC copper). `k_Cu` already reaches 1200 K: find and
  cite its source (RRR-dependent; the deck's RRR must match).
- **Hastelloy C-276, cp, k, ρ:** cryogenic range Lu, Choi & Yang 2008 (Cryogenics 48, thermal
  properties of Hastelloy C-276 at cryogenic temperatures); 300–1300 K the Haynes International
  C-276 datasheet (cp, k, electrical resistivity to 1200 °C). Verify the DOI of the 2008 paper and
  cite the datasheet by title, revision and year.
- **REBCO / YBCO, cp and k:** cp of YBa₂Cu₃O₇ above 300 K from the calorimetric literature (the
  Dulong–Petit limit is the sanity check); k of a coated conductor above 300 K from Bonura &
  Senatore 2015 (Supercond. Sci. Technol., thermal conductivity of REBCO tapes) and the references
  therein. State which layer the table stands for (the deck's `k_HTS` is the ceramic, not the tape).
- **Silver and the buffer:** confirm the existing sources; extend `cp_Ag`/`k_Ag` from NIST-JANAF
  (Ag, cr) if the current rows stop short.

## 3. Ordered Steps

- [ ] **R1** Literature pass (web): for every table, the source, its DOI or datasheet id, and its
  valid range; write the provenance header lines (`# source: … DOI … range … K; units …`).
- [ ] **R2** Assemble the extended tables; keep the existing rows where the new source agrees within
  its stated uncertainty, else replace the whole table from one source and say so in the header;
  monotone in T where physics requires (cp), no kinks at the join (plot each).
- [ ] **R3** Codex terra/high number audit: every row against the cited source; Codex reads the
  source text it is given, not the web.
- [ ] **R4** README table section: one row per table with source, DOI, range.
- [ ] **R5** Gate: the T5 deck to 50 ms on a debug build; the envelope line in the log now reports
  the new minimum; DR-193 struck with the gate named.

## 4. Open Questions

- **O1** Whether the usermat plugin should refuse a temperature above the table instead of holding
  the endpoint. Christian's O3 ruling (2026-10-03) kept the hold; not reopened here.

## 5. Audit Trail

- Exchange thread: `tmp/ai_exchange/review_dr193_tables.md` (R3).
