# Gas Modules: Doxygen Completion, Comment Defects, Dead Code

**Date:** 2026-08-06
**Purpose:** Document the `gastables` / `gasmodels` API, fix comments that
contradict the code, and remove dead code. Scans by Grok (comment defects) and
Codex (dead code); fixes by Claude.
**Module:** physics/gastables, physics/gasmodels
**Status:** Class docs and comment defects done; dead code partly done.
Nothing compiled — build is Christian's.

## Baseline

gasmodels carried 53 doxygen tags across ~627 declarations; **gastables carried
three across ~250**. Of 35 classes, 15 had no description at all — essentially
the whole gastables layer plus `Statevals`, `AlphaFunctionFactory` and
`EoS_Hydrogen`.

A first survey was wrong and was redone: `grep -m1` matched *forward*
declarations, so it attributed classes to the wrong files (`RefGas` to
`cl_Gas.hpp`). Matching real definitions gave the 20/15 split above.

## Class documentation

All 15 now have one, written to carry what cannot be recovered from the code:
ownership, the units actually returned, and the non-obvious contracts. Examples:
`Statevals` records why a never-written slot reads as a default rather than an
error and why `remix` must reset it; `TransportPoly` records that `mScale` is
why viscosity and conductivity cannot share an object (micropoise vs
microwatt per cm K); `HeatPoly` records that the two integration constants
belong to the reference state, not the polynomial, so an interval lifted out of
its parent is not self contained. `RefGas`'s ~25 property accessors were grouped
under three `@name` blocks rather than annotated one by one — the useful content
is the convention, stated once.

## Comment defects (Grok scan, 20+ findings, all verified before fixing)

Three where a comment actively contradicted correct code:

- `cl_GM_EoS.hpp` — the block above `dsdepdp` was a verbatim copy of its
  sibling's, calling the pressure derivative a temperature derivative.
- `cl_GT_RefGas.hpp` — `mHaveConductivity` documented as "has viscosity data",
  the same text as the flag above it.
- `cl_Gas.cpp` — `// scale unit to kJ/kg` over `tData /= mM`, which produces
  J/kg; a factor of 1000 in the comment.

Also: `bar->atm` on a value that is in Pa; "isotropic" where the code is
isentropic (twice); `mKind` credited for a selection made by `mType`; the
nitrogen DOI carrying an extra digit; `thermo.imp` for `thermo.inp`; Leachman's
Eq. (31) pasted onto both derivative bodies; and three `@param` blocks
documenting arguments their functions do not take (`prandtl_meyer` listed
`aIndex`, `phi0` and `omega` document parameters on nullary functions).

**Two of the defects were mine, made earlier the same day.** The `TransportPoly`
class doc I had just written credited `mKind`; and the batch-fix script matched
the first of two identical comments, rewriting a *correct* one on
`mDeltaPowD` while leaving the wrong one on `mTauPowT`. Both repaired, and both
stamps then verified against the code that sets them.

**Corrected a substantive error in my own `RefGas` doc.** It claimed the factory
anchors so that `H( 298.15 )` equals the assigned enthalpy. It does not:
`cl_GT_RefGas.cpp:227` uses `Href + Hf`, so H( 298.15 ) is the formation
enthalpy **plus** the sensible term, and it is H( 0 ) that equals the formation
enthalpy at 298.15 K. See the reference-enthalpy resolution in
`dl20260806_gas_literature_verification.md`.

## Dead code (Codex scan) — removed

Each verified repo-wide before deletion:

- `Gas::mCommRank`, `Gas::mTemperatureSteps` (filled once, never read; the
  factory's own `create_temperature_steps` stays live), `HelmholtzTransport::
  mRefgas`, `EoS_Oxygen::mTau0`, `RefGas::mMode` (the mode is carried by the
  function pointers `set_mode` installs).
- `RefGas::set_href_flag()` — **declared with no definition anywhere**, a link
  error waiting for its first caller.
- **`EoS_Cubic::mUseAlwaysCardano` and the Newton branch behind it.** The flag
  was initialized true and never written, so the branch was unreachable — and it
  had not been kept in step: it carried neither the liquid root selection nor
  the `Z > B` filter added earlier the same day, so enabling it would have
  returned the gas root for liquid states. 33 lines removed and the live Cardano
  body dedented; brace balance verified across all 33 edited files.

## Not removed, deliberately

- **Public virtuals with no in-tree callers** — `EoS::hvap`, `EoS::pi`,
  `EoS::dpidp`, `Gas::hd`, `data_available()`, `deriv_test()`. These are API
  surface or hooks for the planned tablegas, not dead weight;
  `data_available()` in particular is called for by R4 of the migration plan.
- **The four unused `Statevals` slots** (HVAP, PI, DMDT, DMDP). Renumbering the
  macro block to reclaim four doubles per Gas is more risk than it is worth, and
  three of them are explicitly reserved for tablegas.
- Commented-out code blocks (~15 sites) — cosmetic, not yet swept.

## The transport glue abort, second round (Fable)

Christian's rebuild still aborted in `create_glue_polys_transport` — same
assert, reached **through** the morning's smoothness guard. The guard was
necessary but not sufficient: the regression has a second mechanism the heat
case never had. Simulating the exact C++ search over the shipped `trans.inp`
reproduced the abort and located it — CO conductivity at the 162.5 K junction
(first in air processing order, matching the trace), with CO2 C@250 K and
Xe C@201.2 K behind it. The aborting junctions are NOT separable by kink size
(CO's aborter has a smaller slope jump than a neighboring junction that
succeeds); what distinguishes them is geometry: the two branches curve in
opposite directions, so a glue with sign-constant curvature cannot exist at any
window width — the quintic matches both end curvatures and a sign change is
forced. The old condition demanded the impossible.

Fix: two-phase acceptance. Phase 0 is the historical zero-sign-change search,
untouched, so every junction it accepted is glued byte-identically. Phase 1
runs only on exhaustion and accepts exactly one curvature sign change — the
inflection the data demands, still rejecting wiggly glues. Full-table
simulation (34 species, `glue_final.py`): zero aborts, four phase-1 rescues
(CO, CO2, Xe, and D2 C@106.3 K — the deuterium junction that would have hit
the first fusion-fuel run), all at the minimal dT = 5 window.

Also cleared in the same round: two stale-looking diagnostics after the
dead-code removal were verified harmless (`mTau0` fully gone; the
`mUseAlwaysCardano` mention lives only in a comment), and one real orphan the
first removal pass missed — `HelmholtzTransport::mRefgas` had lost its
initializer but kept its declaration, an uninitialized-reference compile error
Christian's build caught; declaration removed, and the remaining removals
re-checked to 0 references each.

## First-run failures of the rebuilt tables, layer three (Fable)

With the two-phase glue in, Christian's next run reached
`Gas::create_reference_gases`: the formation table anchors each element on its
standard state, and the gas-only table rebuild had shipped none of the
condensed references. Closure over the shipped compositions gave the complete
reachable set — C(gr) (any carbon species, so default air), Cl2, and the
bromine reference — rather than fixing one abort per rebuild. All three added
to `species.txt` and both tables regenerated with the README's canonical
invocation; verified zero pre-existing lines altered, three thermo records
(+ Cl2 transport) appended. Two side findings: `Gas::reference_element`
mapped bromine to Br2(cr) where the 298.15 K standard state is the liquid —
fixed to Br2(L) with rationale; and element D falls through to monatomic D as
its own reference — no abort, routed to Christian as D15.

A process error of mine worth the record: the first two regeneration runs used
default flags and silently dropped the NIST-sourced low temperature transport
intervals of CO, D2, Kr, Ne and Xe. The canonical invocation requires
`--nist-transport`. The damage was caught by diffing the output against the
committed table before shipping — which is the check that should precede any
regenerated data file, and the reason the tables were committed first.

## First full run of the gas demo (Fable)

With D13 and D14 in, `bin/gas` runs to completion and its air table is
physically sensible (mu = 18.55 uPa s and lambda = 26.1 mW/m K at 300 K; the
h column shows the 0 K anchoring including CO2's negative formation share).
Two observations from Christian's run resolved:

- **"The header is gone"** — `bin/gas` (gasmodels demo main) never printed one
  and ignores `-g N2`: it hard-codes the default air sweep and writes
  `gasdata.h5`. The species-selecting executable with the header is
  `bin/gastable`. Confirmed from the output itself: the viscosity column is
  air's, not nitrogen's.
- **HDF5 error at exit** — `main` computes `return gComm.finalize()` before
  the file object's destructor runs; `PetscFinalize()` closes the HDF5
  library when PETSc carries HDF5 support, so the destructor then closes a
  dead handle. Fixed with an explicit `tFile.close()` before the finalize,
  with the rationale in a comment. Almost certainly as old as the demo main;
  today is simply the first time it has run since the migration. The other
  session's `cl_HDF5` change (a guarded `create_group` overload) was checked
  and is unrelated.

## Open

- [ ] Codex prose pass over the new documentation
- [ ] Commented-out code sweep
- [ ] Build + run: two-phase glue (D13) and element references (D14) both
      await Christian's executable gate
- [ ] D15: deuterium reference state (D vs D2) — Christian
