# Nédélec Edge Function Defects: EF_TET4 (confirmed) and EF_TET10 (suspected)

**Date:** 2026-08-14
**Purpose:** Track two edge-function defects surfaced as by-catch of the Nédélec LaTeX
extraction audit round (`closed/nedelec_tex_extraction_plan.md`, §4.0/§4.0b). Both are in the
`E` interpolation operators; the `C` curl operators are correct where checked.
**Module:** `src/fem/interpolation/nedelec`
**AIs involved:** Grok (found both), Claude (independently confirmed D1 by hand), Codex
(confirmed the D1 callout in the docs round)
**Status:** **CLOSED 2026-09-03 by the tree-first `todo/` sweep. Both defects are fixed in the
tree and the regression battery exists.** D1's fix is in `cl_EF_TET4.cpp` — edge 2 now builds
`η∇ξ − ξ∇η` through `mNablaEta`, not `mNablaZeta`. D2's TET10 fix landed the same day. The
circulation battery is `tests/fem/test_EdgeFunctions.cpp` with fixture
`tests/fem/support/cl_EF_TestVolume.hpp`, wired into `tests/fem/CMakeLists.txt`; it ran 22/22
green against a fresh `cl_EF_TET10.o` and reproduces the pre-fix signature against the stale
library as a negative control. **Residue: a run gate only** — Christian's rebuild, `make check`
in a `USE_TEST=ON` tree, and the 3D bulk-conductor run. Struck is not verified.

> **The status line this replaces was false, and is recorded here because the failure mode
> matters more than the file:** it read *"OPEN — awaiting Christian's ruling. **No source
> modified.** All findings are static hand algebra; no numeric probe has run."* Every clause of
> that was untrue by the end of the same day it was written — R1 and R2 are `[x]` in this very
> file, ticked with Christian's own in-tree edits and a sympy probe. A header that contradicts
> its own checked boxes is the exact reason this sweep stopped reading status lines and started
> reading the tree.

> **Scope guards:**
> - Fixes are one-line-class but touch production element code; they need an explicit go.
> - A regression gate must run after any fix (see R3); "it compiles" proves nothing here.

---

## D1 — `EF_TET4::E()` edge 2 uses `∇ζ` where `∇η` belongs (CONFIRMED, HIGH)

- [x] **Fixed 2026-08-14 (Christian, in-tree edit of `cl_EF_TET4.cpp:194-196`), after his ruling
  "I agree with that D1" and an exact sympy probe** (`scratchpad/tet4_circulation_probe.py`,
  session scratchpad; result distilled below). The probe, on a random rational-vertex tet,
  reproduced the defect exactly (circulation matrix = identity except column 2: +1/2 own edge,
  −1/2 on edge 0) and showed the corrected basis gives the exact identity, restores
  `mC == curl(E)` in all six columns, and equals the Whitney basis validated against
  DefElement's published degree-0 N1 tetrahedron functions (defelement.org, different
  parameter space and edge enumeration, mapped via barycentrics). Formula-level verification;
  the compiled binary has not yet re-run a 3D case (R3 still owed).

**The defect.** `cl_EF_TET4.cpp:193-195`: edge 2 (Exodus nodes 3→1, barycentric η→ξ under the
implemented map node1↔ξ, node2↔ζ, node3↔η from `cl_IF_TET4.hpp:28-31`) is coded as
`η∇ξ − ξ∇ζ`. The Whitney form is `η∇ξ − ξ∇η`. The code's own comment (`:192`) and the curl
operator `mC(:,2)` (`:122-124`, `= 2∇η×∇ξ`) both carry the correct form, so `E` and `C`
disagree within one class.

**Consequences (hand-derived, Claude + Grok independently):** the coded basis has circulation
1/2 instead of 1 on its own edge and −1/2 instead of 0 on edge 0 (ξ→ζ). Tangential conformity
of the interpolated field is broken on every TET4 wherever `E` is consumed: the `μ E^T E` mass
term and the L2 projections. The stiffness path through `C` is unaffected, which is presumably
why bulk-3D TET4 results have been plausible enough to survive.

**The fix:** `mNablaZeta` → `mNablaEta` in the three component lines `:193-195`.

**No test coverage:** the circulation battery planned in `todo/falsification_tooling.md` (D1
there) would catch exactly this and is unimplemented.

## D2 — `EF_TET10` scalar tables transcribed in the notes' labels, not the Exodus map (CONFIRMED, HIGH for any TET10 use)

- [x] **FIXED 2026-08-14 (Claude, MATLAB-first workflow on Christian's direction).** Root cause
  traced to one line: `tmp/tet10/tet10_generate.m` used the naive `lambda = [xi; eta; zeta; tau]`
  while its sibling `defelement.m` carried the correct swapped map. Fix chain, each stage gated:
  (1) generator line corrected; (2) pinned tables in `tet10_function.m` / `tet10_derivatives.m`
  re-pinned, all residuals exactly zero under Octave/SymPy **and independently in Christian's
  MATLAB run**; (3) tables ported into `precompute()` (`cl_EF_TET10.cpp`, 14 tables); (4) final
  gate parsed the **edited C++ file itself** back into sympy and passed all four checks: scalar
  tables == generator truth, all nine derivative tables == exact derivatives, exact
  block-identity circulation + zero face leakage on a random tet, and the `C()` chain-rule
  assembly == symbolic curl for all 24 dofs (`scratchpad/tet10_port_gate.py`). The assembly
  functions (`compute_edge_functions/derivatives`, `compute_face_functions/derivatives`, `C()`)
  were reviewed and confirmed correct throughout; only `precompute()` tables changed.
  `g++ -fsyntax-only` with the build tree's own flags: clean. Formula-level verification;
  binary gate = R3.

**Confirmed 2026-08-14 by the exact sympy probe** (`tet10_circulation_probe.py`, session
scratchpad), upgrading Grok's suspicion. Verbatim transcription of `precompute()`
(`cl_EF_TET10.cpp:247-310`) plus `compute_edge_functions`/`compute_face_functions`
(`:599-733`) on a random rational tet:

- **Coded set: 24 edge-dof conformity violations.** Base-edge pairs and the ζ→τ / η→τ spokes
  leak circulation onto foreign edges (entries up to ±1; own-edge values 2/3, −1 instead of
  1, 1). Face candidates 3, 4, 9, 10 carry circulations ±4/3, ±8/3 on edges where they must
  vanish. Only the (ξ,τ) pair (edge 3, dofs 6-7) is clean — its scalars contain no η or ζ.
- **η↔ζ-swapped scalars (gradients untouched): zero violations.** Exact block-identity edge
  circulation matrix, all face candidates zero on all edges, and the swapped set equals the
  notes' polynomials remapped to the implemented pairing (`nedelec_derivation.md` §4.1),
  symbol for symbol.

**Root cause:** the scalar tables `mG, mH, mU, mV, mW` are in the notes' generic
`(ξ,η),(η,ζ),(ζ,ξ),…` labels; the gradient assignments in the assembly functions are in the
Exodus map (node1↔ζ, node2↔η). Same slip class as D1, table-wide.

**Fix scope (NOT a three-liner):** swap η↔ζ in every `mG/mH/mU/mV/mW` entry, re-derive the
nine hand-written derivative tables (`mGxi…mWzeta`, noting the η/ζ derivative tables also
exchange contents under the swap) which feed `compute_edge_derivatives`,
`compute_face_derivatives` and thus `C` on both straight and curved paths. `mNxi/mNeta/mNzeta`
(Jacobian shape derivatives) are already correct and must not be touched. Proposed workflow:
extend the probe to transcribe and check the derivative tables and assembled `C` columns
before and after the edit.

## D3 — `EF_TET10` shape-derivative tables `mNeta`/`mNzeta`, rows 8-9 crossed (CONFIRMED, FIXED 2026-08-14)

- [x] **Found by the new test battery, fixed same day, runtime-verified.** The
  `mNxi/mNeta/mNzeta` tables in `precompute()` build the Jacobian; rows 8 and 9 of `mNeta` and
  `mNzeta` (the ζτ and ητ midside nodes) carried each other's values — the same naive-map trap,
  four entries. Consequence: the η- and ζ-columns of the TET10 Jacobian were **identical**, J
  singular, `1/detJ` amplifying the nablas to ~1e14 wherever the Jacobian was evaluated at a
  generic point. **Masking lesson:** `link()` evaluates the Jacobian at evaluation point 0; with
  point 0 on edge 0 (η = 0) the crossed entries coincide with the correct ones, so the straight
  path looked healthy and only the curved path (per-point Jacobian) exposed it. The battery now
  puts an interior point at column 0 so both paths see a generic Jacobian. `mNxi` and the `mD`
  second-derivative table were checked row-by-row against the pinned shape set in
  `tmp/tet10/tet10_lagrange.m` (whose node table is correct) and are clean. Diagnosed via
  valgrind (clean → computed garbage), then a private-access probe printing the Jacobian
  entries (b=c, e=f, h=i → singular).

## Ordered steps (after ruling)

- [x] **R1** — Apply the D1 three-line fix. *(Christian, 2026-08-14, in-tree; probe-verified)*
- [x] **R2** — Sweep `EF_TET10` against `nedelec_derivation.md` §4.1. *(Done 2026-08-14:
  diagnosis, MATLAB-source fix, C++ port, four-gate acceptance — see D2 above.)*
- [◐] **R3** — Gate. **The circulation battery now exists**:
  `tests/fem/test_EdgeFunctions.cpp` + fixture `tests/fem/support/cl_EF_TestVolume.hpp`
  (TRI3/TRI6/TET4/TET10; circulation identity on reference + distorted elements, C-vs-FD-curl,
  detJ > 0 EXODUS guard, edge-flip sign, TET10 curved-path equivalence; wired into
  `tests/fem/CMakeLists.txt`). Probe run against the prebuilt libraries with the fresh
  `cl_EF_TET10.o`: **22/22 green**, and together with the pre-existing Lagrange suite 46/46.
  Negative control: against the stale (pre-fix) library the battery reproduces the D2
  circulation signature (2/3, −1) exactly. Remaining: Christian's rebuild, `make check` in a
  `USE_TEST=ON` tree, and the 3D bulk-conductor run.
