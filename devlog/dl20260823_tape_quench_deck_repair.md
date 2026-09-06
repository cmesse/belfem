# tape_quench deck repair — DR-02's coverage-gap exerciser, three plugin bugs fixed

**Date:** 2026-08-23
**Purpose:** Get `build/tape_quench` (Christian's restored deck, coupled magnetic+thermal,
custom `hts` material with `resistivity type : piecewise`) running, to close DR-02's
named coverage gap — no HEAD-runnable deck exercising the thermal `T_h`/piecewise family.
**Module:** physics/materials (user plugin), fem/kernel (Controller)

## Starting state

`build/tape_quench` had `input.conf`, `MatData/matlib.cpp` (source only), a mesh, and
nothing else — no build scaffolding, no `lib/` (current/defect plugins), and
`input.conf` carried a literal unresolved git merge conflict in the `nonlinear` solver
block (commit hash `8822651de...` not present anywhere in this repo's history).

**Christian's ruling:** resolve the conflict using the clean reference deck's settings
(`tmp/examples/Tape_Quench/CustomMat/input.conf`: `algorithm : Newton`, `max relaxation
: 1.0`, `min relaxation : 0.1`, `tolerance switch : 1e-18`). Confirmed the source
(`matlib.cpp`) is byte-identical to that reference, so copied its `CMakeLists.txt`,
data files, and `lib/` (current.cpp, defect.cpp) into `build/tape_quench`, fixed the
hardcoded `BELFEM_DIR` in both CMakeLists (`gregorygiard`'s path → this machine's),
built both plugins clean.

## Three real defects found in sequence, each surfaced by the previous fix

**1. Reserved-word collision — DR-65's exact mechanism, on a second deck.**
The tape's real buffer layer material was labelled `buffer`, colliding with BELFEM's
internal thin-shell reserved name (`cl_MaterialFactory.cpp:173`, `Material name
'buffer' is reserved...`). This confirms DR-65 is not a one-off dead deck — the
`CustomMat`-family deck hits the identical collision. **Christian's ruling: rename to
`custom_buffer`** (his direct instruction overrode an initial AskUserQuestion pick).
Four sites: `input.conf`'s material section name, its `label :`, the `layers : tape`
stack entry, and `matlib.cpp`'s `buffer_init` → `custom_buffer_init`.

**2. Missing `T_crit`.** `rho_piecewise`'s `T > constant_property(T_crit)` guard
(`powerlaws.hpp:752`) asserted — `hts_init` never called `set_constant(T_crit, ...)`.
Added `92.5 K`, matching `cl_Material_YBCO.cpp:87`'s literature-backed builtin value —
not a guess, the in-tree canonical number for this exact material family.

**3. jc/n registered through the wrong dispatch.** `hts_init` used the 1-dependency
`set_user_defined_function(Property, MaterialDependency::T, &fn)` overload for `jc`
and `n`. Traced via `addr2line` on the debug binary (two successive aborts at
different offsets inside the same inlined `rho_piecewise`/`jc_eval` region) and
confirmed by reading `cl_Material_UserDefined.cpp`: the 1-dependency overload
(`:55-65`) unconditionally stores into the generic `mUserFunctions` table regardless
of `Property` — it never touches `mJcFunction`/`mNFunction`. But `jc_eval()`/`n_eval()`
(`powerlaws.hpp:62-79`) only ever read `mJcFunction`/`mNFunction`, falling back to
`constant_property()` when null — which was also never set, hence the second assert.
The correct path, `JcFunctionUserDefined` via the 3-dependency overload, is the exact
mechanism `mu(H,T)` already uses at `cl_Material_UserDefined.cpp:90-100`. Fixed by
reshaping `hts_jc`/`hts_n` to `MatFunc3` (`normB`, `angle`, `T`; first two unused —
these are genuinely T-only polynomial fits) and re-registering through the
3-dependency overload. `JcFunction`'s base-class derivative defaults (`deval_dB` etc.
return 0.0) are documented as the correct fallback for exactly this case
(pre-DR-07 tangent contract), so no Jacobian gap was introduced.

Each fix was verified by execution — it unblocked the next distinct abort in the
sequence, never the same one recurring.

## Where it stands

With all three plugin bugs fixed, the deck reaches real coupled Picard/Newton
iteration for the first time — genuine progress on DR-02's coverage gap. But it does
not converge: magnetic residual oscillates 5–14 dB across 13 Picard steps plus one
Newton escalation, the timestep collapses to the 1e-6 s floor, and
`Controller::reset_timestep` aborts on floor-retry exhaustion (`cl_FEM_Controller.cpp:2267`,
20 retries exhausted). Ruled out: an oversized first-step current ramp — the deck's
`userdefined` current (table lookup, `I_vs_t_regular_smooth.txt`) is ≈ −2.4 A at
t = 0.1 ms, against a jc order 10¹⁰–10¹¹ A/m² — deep superconducting regime.

**Open question for Christian, not resolved this session:** is the convergence wall
deck-specific tuning (the abort message itself suggests `compute conditioning : true`
to check for a conditioning floor, or an `absolute tolerance` in the nonlinear
section), or a genuine defect the kernel-collapse gate exists to catch? Deliberately
not tolerance-loosened past this point — that would defeat the point of the gate.

## Register maintenance

DR-02's row amended in-session with the full repair sequence and the open convergence
question, per protocol §11.
