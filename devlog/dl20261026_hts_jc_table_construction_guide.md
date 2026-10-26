# HTS jc / n table construction guide

**Date:** 2026-10-26
**Purpose:** Christian asked for an instruction on how to build an HDF5 jc / n table for an HTS
cable, with emphasis on making the B-spline mapping work
**Module:** `src/physics/materials/doc` (docs only; no code changed)
**AIs involved:** Claude Opus 5.5 (primary); Codex `gpt-5.6-luna`/medium (language sweep)
**Status:** COMPLETE — guide written, swept, indexed in the module README
**Verification:** probe level, numpy route only — a synthetic anisotropic Kim table built with the
guide's minimal writer fails `verify_tables.py` on C1 axes 0/1/2 and the θ seam raw, and passes
every gate after `project_table.py` (face jumps ≤ 2.2e-12 % of the reference slope). The
`make_extension.py` → `mattest` route and a BELFEM load of the result were **not** run
**Exchange:** `tmp/ai_exchange/hts_jc_table_guide.md` (Codex sweep)

## Summary

New user guide `src/physics/materials/doc/hts_jc_table_construction.md`: the `jc`/`n` group
format (axes T / log₁₀ B / θ from the tape normal, log₁₀ values, C order θ fastest — matches
`TensorMeshConfig::node_index` = `k + nA·( j + nB·i )`), the grid rules the order-2 B-spline
layer imposes (odd node counts, θ exactly [0, π], seam column closed), the projection pipeline
from the `python` group embedded in every shipped REBCO table (`bspline.py`, `make_extension.py`,
`mattest`, `project_table.py`, `verify_tables.py`), why contracts go on the control net rather
than the node values, what `JcFunctionDatabase::open_checked` refuses at load (and that it never
checks C¹), the optional self-field anchor (`source/columns`+`points`, calibrated to
`meta/Icw_77p5K_sf_A_per_m / t_eff_m`), and cable-specific points.

## Findings worth keeping

- **No tool or doc existed for a user-built table.** The pipeline lived only inside the shipped
  HDF5 files' `python` groups and in `share/material/README.md` §"Using an older or modified
  table"; `build.py` is specific to the HTSDB/Robinson export layout. High confidence.
- **`layers` thickness must equal `meta/t_eff_m`, unchecked.** `t_eff_m` is read only for the
  self-field calibration; the deck's layer thickness scales every `Ic`. `examples/rutherford`
  satisfies it (310 µm). High confidence.
- **With `file` set, `ec` is not read** (`cl_MaterialFactory.cpp`, HTS branch) → criterion fixed
  at 1e-4 V/m; data at another criterion must be converted via `( E₀/E₁ )^(1/n)`. High.
- **A cable-level self-field Ic should not anchor a tape table** — it would count the cable's own
  field twice. Physics judgement, medium confidence; stated as such in the guide.
- Codex sweep flagged one real ambiguity, fixed: `interp_curve` (calibration) holds the raw curve
  at its end values, while `self_field` (evaluation) continues linearly from the end knots.

## Open

- `mattest` route and BELFEM load of a user-built table not run for this guide.
- `verify_tables.py`'s third `howtoread.py` fallback is a hard-coded `/home/christian/...` path in
  the embedded copy; harmless (it tries `$BELFEM_DATA` first) but not portable.
- Homogenized whole-cable tables are declared out of scope / untested.

## Files Updated

- `src/physics/materials/doc/hts_jc_table_construction.md` (new)
- `src/physics/materials/doc/README.md` (index row)
