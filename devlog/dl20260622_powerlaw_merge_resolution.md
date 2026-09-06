# Devlog 2026-06-22 — Power-law `cl_Material` / `powerlaws.hpp` merge resolution

**Date:** 2026-06-22
**Topic:** Resolved the `devel` ← `periodic_new` merge conflict in the HTS power-law resistivity functions
**Module:** physics/materials
**AIs involved:** Claude Opus (Christian directing)
**Claude Confidence:** high (per-function 3-way diff + syntax-only compile)
**Codex Audit Confidence:** N/A (not run)
**Literature References:** N/A (mechanical merge; math unchanged — power law per Messe et al. 2023, paper1)

## Summary

The in-progress `periodic_new` (1491a8d8) → `devel` merge had a botched resolution of the
power-law resistivity functions. The two branches had diverged in opposite directions:

| | `devel` (HEAD) | `periodic_new` (MERGE_HEAD) |
|---|---|---|
| `powerlaws.hpp` | did **not** exist | **created** it; moved all power-law bodies there |
| `cl_Material.hpp` | bodies **inline**, no `#include` | bodies removed, ends with `#include "powerlaws.hpp"` |
| material logic | got the `cf4d8c76` fix | got assert/doxygen cleanup + tweaks |

The working-tree resolution had kept `devel`'s inline bodies in `cl_Material.hpp` **and**
dropped the `#include "powerlaws.hpp"` line — leaving the new `powerlaws.hpp` orphaned
(included nowhere, hence no duplicate-definition error, but all of `periodic_new`'s
refactor + cleanups silently discarded).

Resolution: **keep `periodic_new`'s structure** (power laws live in `powerlaws.hpp`,
included once from `cl_Material.hpp`) and port **only** the genuine `devel`-side change onto it.

## Key Findings

- A per-function, whitespace-insensitive 3-way comparison (ancestor `d2e6306d` vs `devel`
  vs `periodic_new`) showed the **math is identical** in every shared function. The
  apparent "34 changed" was noise: `periodic_new`'s copies are merely *cleaner* (added
  `BELFEM_ASSERT` guards, fixed a copy-paste assert message — `devel`'s `n` guard read
  `"...constant jc parameter..."`). So letting `cl_Material.hpp` "supersede" wholesale
  would have **re-introduced that bug** and deleted the doxygen.
- The **only** genuine `devel`-side logic change is one commit, `cf4d8c76`
  ("Minor fix to handle user-defined Jc(B,theta)", G. Giard), and it is exactly:
  - **(A)** three new overloads `periodic_new` lacked: `n(normB,angleNxB)`,
    `jc(B,angleNxB)`, `jc(B,angleNxB,x,y,z,t)`;
  - **(B)** in the **8** field-dependent (`normB,angleNxB`, no-`T`) functions
    (`rho_powerlaw`, `rho_piecewise`, `drho_powerlaw_dJ`, `drho_piecewise_dJ`, each plain
    + defect `x,y,z,t` form): `n = constant_property(MaterialProperty::n)` →
    `n = mNFunction->eval(normB, angleNxB)`, so `n` tracks `Jc(B,θ)` like `jc` already did.

## Changes Made

- **`cl_Material.hpp`**: reconstructed from `periodic_new`'s version (proven equal to the
  correct merge = `periodic_new` + the 3 declarations, since `cf4d8c76` is `devel`'s only
  post-ancestor change and its non-body part is just those declarations). Removed all
  inline power-law bodies; restored `#include "powerlaws.hpp"`; kept the 3 new overload
  declarations.
- **`powerlaws.hpp`** (base = `periodic_new`'s cleaner version): applied the `cf4d8c76`
  fix — the 8 `n`-evaluation edits + 3 new overload bodies.
- **Null-safety hardening** (Christian's call — each new overload mirrors its
  temperature-bearing sibling rather than `cf4d8c76`'s as-written body):
  - `jc(B,θ)` → `mJcFunction ? eval(B,θ) : constant_property(jc)`. This removes a **real
    release-mode null dereference**: `have(jc)` is `true` even for a plain-constant `jc`
    (then `mJcFunction == nullptr`), and the as-ported `assert(have(jc))` is compiled out
    in release, so the old body would have dereferenced null.
  - `jc(B,θ,x,y,z,t)` → `assert(mJcFunction != nullptr)` (matches its T-sibling).
  - `n(B,θ)` → `assert(mNFunction != nullptr)` (matches its T-sibling; was unguarded).

## Verification

- No conflict markers; braces balanced (cl_Material 68/68, powerlaws 148/148).
- Declaration↔definition parity 38/38 (the only diff is the pre-existing `normB`-vs-`B`
  parameter-naming convention — type-identical).
- `mpicxx -fsyntax-only` on a TU including `cl_Material.hpp` → **exit 0** (only an unrelated
  Armadillo OpenMP warning). Re-run clean after the null-safety edits.
- The clangd "main file cannot be included recursively" diagnostic on `powerlaws.hpp` is a
  false positive — it is a mid-file fragment, not a standalone TU.

`git diff --stat HEAD`: `cl_Material.hpp` −1411 / `powerlaws.hpp` +1895 (net move).

## Status / Next

- Merge is **still open and uncommitted** (`All conflicts fixed but you are still merging`).
  Conclude with `git commit` when ready.
- Build caveat (project-wide): a plain `make` will **not** pick up these header changes —
  use `make reset && make <target> -j20`.
