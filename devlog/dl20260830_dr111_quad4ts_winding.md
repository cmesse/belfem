# DR-111: The 2-D Thin-Shell Winding Is a Convention, Not a Defect

**Date:** 2026-08-30
**Purpose:** Record the DR-111 fix round, in which the register row's own root cause and its
severity claim were both refuted, the fix moved from the factory to the consumer on Christian's
ruling, and an implementation that would have shipped silently inert was caught by an auditor.
**Module:** fem/kernel (+ fem/interpolation/nedelec, fem/thermal)
**Session type:** investigation → plan ×2 audit rounds → implementation → code audit → strike

## Summary

DR-111 was filed as "the 2-D thin-shell thermal assembly aborts on a negative layer-element
Jacobian", with a root cause naming `ThinShellFactory` and a remedy that would have renumbered the
QUAD4TS corners. Both halves were wrong.

The winding is **unconditional**, not orientation-dependent: `process_nodes_line2` builds the
normal as `n = ( t_y, -t_x )` (`cl_ThinShellFactory.cpp:923-924`), the tangent rotated clockwise,
so `n` is derived from `t` and reversing the facet reverses both. Replicating the four functions
in isolation gives `detJ = -5.000000e-12` for `+x`, `+y`, 37° **and the reversed facet**, matching
the row's gdb value to every printed digit. The 3-D branch uses the right-hand normal (`:1101`)
and comes out `+2.0e-16` for CW, CCW and tilted triangles alike — which is why only the 2-D
thermal path ever aborted.

The winding is also **forced**. Local edge node order must follow the facet edge direction,
because `mS` is read from node order (`cl_FEM_Element.cpp:1337-1400`, called at `:181-183`) while
`∇ξ` is read from the facet (`cl_EF_QUAD4TS.cpp:38-42`); combined with the edge-index contract
(`cl_Element_QUAD4TS.hpp:104-134`, positional linking at `cl_ThinShellFactory.cpp:1926-1958`)
exactly one corner order is legal, and it is left-handed. **No permutation is both CCW and
sign-neutral.** So the factory implements a consistent convention and the defect is downstream:
the thermal path used the *signed* determinant as a volume weight.

Christian ruled the factory must not be touched and all changes stay inside QUAD4TS. Fix:
`Calculator::dV_quad4ts = std::abs( dV_ts( aIndex ) )` at the single QUAD4TS dispatch, and
`Pipette::measure_quad4ts` likewise. `N` is reference-space and `B = J⁻¹·dN/dξ` is correct under an
orientation-reversing map, so `|det J|` is the **exact** change-of-variables weight, not a clamp.
Both shared bodies stayed byte-identical, so PENTA6TS, HEX8TS and the 3-D Pipette branch are
untouched by construction.

## What the rounds actually caught

- **Grok, round 2 — the fix would have been silently inert.** `invJ2D3D`
  (`cl_FEM_Calculator.hpp:1713-1727`) writes the same `mDetJ`/`mDetJIndex` cache `dV_ts` reads, and
  the thermal kernels call `B(k)` before `dV(k)`. QUAD4TS is linear, so `tIndex = 0` and k=0 is a
  **cache hit**: abs-on-the-assignment would have returned the signed value and never run, at
  exactly the Gauss point where the guard measures. Abs on the *return* is immune. This is the
  INC-549 class — a fix that executes and changes nothing.
- **Codex, round 2 — `EF_QUAD4TS::det_J()` is not unconditionally positive**
  (`0.25*mThickness*mLength`, thickness unchecked at parse), which made the wrapper strictly better
  than the copied-body alternative.
- **Both, code round — the same false sentence.** My comment said the magnetic solve "takes the
  other branch". False as control flow: dispatch is by element type, so magnetic *does* enter the
  wrapper. The safety is a **value** claim — `0.25·h·L > 0` makes the abs an identity. Corrected.
- **My own errors, conceded:** six fabricated `file:line` anchors (read through `sed` ranges, then
  written from memory — both auditors caught it); an area identity `-h|t|²` that should be `-h|t|`;
  and a claim that nothing rejects a negative thickness, refuted by Grok pointing at
  `cl_FEM_Kernel.cpp:1027-1036`.

## The severity correction

The row claimed a release run "would integrate NEGATIVE mass/volume silently — the exact DR-06
localized-cold-spot mechanism". For a thermal domain that is entirely thin-shell with no
Dirichlet or flux term — which is what `tape_quench_usermat` is, topology `thinshell` + `air`,
air carrying no thermal dofs — `M`, `K` and `f` all accumulate `× dV` and are **uniformly**
negated. The timestep system is homogeneous of degree one in them, so both sides scale by −1 and
the temperatures are unchanged. The real exposure is narrower: a thermal domain **mixing**
thin-shell with volume elements, or carrying a positive-measure BC term, where the negation is
non-uniform. That is why the defect only ever surfaced as a debug abort and never as a wrong
answer.

## Evidence, and what is still owed

`make check` green on the rebuilt binary; the fix confirmed present **by symbol** (`nm -C` shows
`dV_quad4ts` and `measure_quad4ts`) rather than assumed from the build; a coupled 2-D thin-shell
quench deck ran 750+ timesteps with the thermal thin-shell path live and converging.

**That run cannot discriminate, and the row says so.** It was a release binary
(`USE_DEBUG:BOOL=OFF`, `-O2 -DNDEBUG`; `strings` finds zero occurrences of "Negative Jacobian
determinant"), so the guard cannot fire — and by the severity correction above, the pre-fix binary
would have produced the *same* temperatures on that deck. The discriminating gate is the same deck
in an assertions-ON build whose binary is first confirmed to contain the guard string. Magnetic
parity was pre-registered as bitwise and is recorded **unverifiable**, not passed: the only
pre-change log predates ~520 lines of concurrent `SolverMUMPS` work and a modified `input.conf`.

Struck on Christian's ruling under the DR-42/DR-49 pattern. **Struck is not verified.**

## By-catch filed

DR-151 (no positivity check on layer thickness at parse; the guard exists but is late and generic),
DR-152 (`MeshChecker`'s swap table would feed `EF_QUAD4TS` the tape length as its thickness —
unreachable, but advertised), DR-153 (`src/fem/postproc/cl_MeshChecker.{hpp,cpp}` is an orphaned
duplicate in no `CMakeLists.txt`, which already misled one vendor in this round).

## Files

- `src/fem/kernel/cl_FEM_Calculator.{hpp,cpp}` — `dV_quad4ts` + dispatch
- `src/fem/kernel/cl_Pipette.{hpp,cpp}` — `measure_quad4ts` + dispatch, `<cmath>`
- `src/fem/interpolation/doc/nedelec_thinshell.md` — §2a, the convention and its two traps
- `todo/dr111_quad4ts_winding.md`, `todo/debt_register.md`
