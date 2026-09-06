# Campaign: Gas Models Open-Source Migration

**Status seed:** 2026-08-05 from `todo/gasmodels_open_source_migration.md` — every claim
`[seeded — confirm]` until Christian's correction pass.
**Branch:** `sideconnectors` (staged edits) · **Plan:** `todo/gasmodels_open_source_migration.md`

## Current accepted design `[seeded — confirm]`

Move `nonfree/physics/{gasmodels,gastables,atmosphere}` into open-source `src/physics/`,
finite-rate `combustion/` stays private; `share/` excluded for now. LBNL IPO approval
obtained; chain of title clear (code is Christian's own; `@dlr.de` stamps were cosmetic,
removed). Data strategy: the 471-record McGraw-Hill-derived `gasdata.inp` collapses to
263 labels / 32 complete / 18 fully capable — clean rebuild from CoolProp/Cantera is
small; R1/R1b staged a clean 42-species `gasdata.inp` with per-field provenance.

## Last passing reproducer `[seeded — confirm]`

Phase 0 defect fixes (D1–D6, D9 partial) audited Codex+Grok; no dedicated gasmodels
test target yet — reproducer = gastables lookup probes (DR-43).

## Open P0/P1 `[seeded — confirm]`

- D7 label-nomenclature mismatch stranding records; D8 six-way `C7H9N` collision;
  D9 placeholder-CAS last-in-wins hitting H2/D2 (misses `cubicalpha` PM row AND the
  factory CCR entry → silent generic-ω fallback) — DR-43.
- O4: normal vs equilibrium hydrogen — physics call (Christian).
- R2+: the actual tree move + CMake wiring.

## Superseded approaches `[seeded — confirm]`

Keeping the full 471-record transcription (copyright + quality); migrating `share/`.

## Dated entries

todo/gasmodels_open_source_migration.md (living plan) · f6babfb0 (staged R1/R1b)
