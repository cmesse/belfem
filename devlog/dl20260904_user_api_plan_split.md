# User API Header Install: Plan Split, Cheap Items Landed, Spine Deferred

**Date:** 2026-09-04
**Purpose:** Resolve the status of `todo/user_api_header_install.md` (revision 2, 2026-08-29)
**Module:** `CMakeLists.txt`, `src/fem/kernel`, `src/numerics/sources`, `src/physics/materials`
**AIs involved:** Claude (Fable); no audit round — comment, doc, guard-rename and one install rule

## Status found

The plan was still at revision 2 with its re-audit owed and its core untouched: the whole-tree
header glob at `CMakeLists.txt:367`, the umbrella in `src/fem/kernel/` under
`BELFEM_BELFEM_USER_API_HPP`. What had landed was by-catch from 2026-08-31 (R5, D2, D3: templates
and example sources install, all four example decks build against a prefix) plus three items the
tracker did not record: the dead `Vector` forward declaration (D8) was gone from `cl_Material.hpp`,
the umbrella carried the LBNL block, and `UserLibraryTemplate.cmake` had made `BELFEM_BACKEND`
optional and off by default.

Two things had changed since the plan was written. An installed prefix now builds plugins, which
was the pain the plan answered; and R1 would have to rewrite the include-path blocks in four decks
and two templates that the 08-31 repairs had just fixed. Recommended a split; Christian agreed.

## Landed (all layout-independent)

| Item | Change |
|---|---|
| R7 / D4 | `install( FILES ${CMAKE_SOURCE_DIR}/LICENSE DESTINATION . )` after the templates rule — every shipped header says "see the top-level LICENSE file" |
| D1 (half) | umbrella guard `BELFEM_BELFEM_USER_API_HPP` → `BELFEM_USER_API_HPP`; prohibition note ("NEVER include cl_Vector.hpp …") and the `gTbulk` resolved-at-load sentence added |
| D9 | `cl_SourceFunction.hpp` guard `CL_FUNCTION_HPP` → `BELFEM_CL_SOURCEFUNCTION_HPP` (house style, safe in a flat public directory) |
| D10 | doxygen compile line in `cl_Material_UserDefined.hpp` now gives the three `-I` directories the module layout needs (`core`, `containers`, `physics/materials`) and points at the installed `UserMaterialTemplate.cmake` |
| tracker | D8 ticked (already in tree), boxes and Status updated, dated banner, file moved to `todo/deferred/` |
| docs | `CLAUDE.md` install sentence names `LICENSE` and the plan's new `closed/` path; `todo/README.md` active count 6 → 5 |

## Gate

- A TU including only `belfem_user_api.hpp` compiles with `g++ -std=gnu++17 -Wall -Werror
  -fsyntax-only`, no backend define, under `-DDEBUG` and `-DNDEBUG`, against the five in-tree
  module directories. A second TU including the umbrella twice plus `cl_SourceFunction.hpp`
  directly compiles too, so both renamed guards guard.
- `scripts/check_doc_claims.py`: 37/37.
- No cmake configure and no `make` were run (user runs builds). The install rule is one line and
  follows the pattern of the rule above it. **Reviewed, not verified**: the `LICENSE` landing at
  the prefix root has not been seen by a real `make install DESTDIR=`.

## Deferred

R1 (flat `install( FILES )`), the umbrella move (rest of R2), R3 (isolated staged gate, which
still compiles only `example_user_material.cpp` against inherited include paths), R4 (template
rewrite), R6 (documents R1 would make false), R8 (end-to-end gate). O1 (merge the templates) and
O5 (R8 depth) remain Christian's call. The revision-2 re-audit is still owed before any spine code.
D10's new compile line goes stale again the day R1 lands — it is written for the layout that ships.
