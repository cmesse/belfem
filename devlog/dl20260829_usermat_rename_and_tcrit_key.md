# `critical temperature` Input Key, `custom` → `usermat`, and the tape_quench Plugin Rebuild

**Date:** 2026-08-29
**Purpose:** record a session that revived the DR-65-family tape_quench example, added one
input key, renamed one input shape, and inventoried the silent-ignore class behind both
**Module:** `src/physics/materials`, `doc/` input contract, `cmake-build-debug/tape_quench_usermat`
**AIs involved:** Claude (exploration, code, runs), Codex + Grok (four audit rounds: T_crit
plan, T_crit code, plugin restructure, materials input logic)

## 1. What started it

Christian dropped a custom-material tape_quench deck into `cmake-build-debug/tape_quench_usermat`
and asked what it takes to run. Four blockers, in the order they fired:

1. **Unresolved git merge markers** in `input.conf` straddling the `nonlinear` block.
2. **`buffer` is a reserved material name** — DR-65's exact mechanism, second sighting.
   Renamed the section; the plugin label inside stays `buffer`. Behaviour-neutral, but NOT
   for the reason `cl_MaterialFactory.cpp:176-182` gives: that comment says ThinShellFactory
   tags the block by NAME. It does not — `cl_ThinShellFactory.cpp:2056` tags by
   `!have( MaterialProperty::rho )`. The builtin `buffer` works because it resolves to
   Magnesia, which has no rho.
3. **jc/n registered through the wrong overload.** `set_user_defined_function(jc, T, &f)`
   takes the `MatFunc1` path, which stores the pointer and calls `set_custom()` — and
   `set_custom` writes NaN into `mConstantProperties(jc)` while setting the have-bit. So
   `jc_eval` (`powerlaws.hpp:63`) took its null-`mJcFunction` branch and returned NaN. The
   fits were loaded and unreachable. Rewired to the three-dependency
   `(normB, angleNxB, T)` overload, which is the only one that builds a
   `JcFunctionUserDefined`.
4. **`T_crit` never set** — the actual first abort, `"Property is not constant"` at
   `cl_Material.hpp:1719`, from `rho_piecewise`. In a RELEASE build this would have been
   silent: `T > NaN` is false, so the superconducting branch runs unconditionally.

## 2. Why 90.0 K is not a judgement call

Evaluated the plugin's own polynomials: the jc quartic has a root at **90.0371 K** and is
negative above it; the n quintic crosses 1 at **91.8784 K** and zero at 92.0791 K. The YBCO
builtin's 92.5 K sits past BOTH — which is DR-105's late-run `n > 1` abort at the quench
front, reached only after a 9.8 h serial run. Single-sourced as `hts_Tcrit` used by both the
fit cutoffs and the registration, and the cutoff made INCLUSIVE (`T <= hts_Tcrit`): all
sixteen assert pairs in `powerlaws.hpp` are gated by a strict `T > T_crit`, so at exactly
`T == T_crit` the solver still enters the superconducting branch and must find polynomials
there, not the normal-state fallbacks.

## 3. The new key

`critical temperature : 90.0 K ;` in a material section, read once in the shared tail of
`MaterialFactory` (before `density correction`), overriding what the plugin's `_init` set.

**Fatal on a builtin** — Christian's ruling, taken after the plan audits and stronger than
what either auditor proposed. Grok had found that YBCO's constructor samples its Callaway
lambda spline (`cl_Material_YBCO.cpp:122`) and that sampling copies `T_crit` into
`params[8]` (`:507`), so a late override moves the resistivity gate and leaves the thermal
conductivity built around the old value. Grok recommended documenting the trap; refusing the
case removes it. `! tIsBuiltin` is the discriminator — the factory's own local, covering all
three builtin-detection forms — checked BEFORE `have( jc )`, because builtin YBCO *has* jc
by the time the tail runs and would otherwise pass.

Codex's one P1 was live and confirmed by execution: `to_real` is `strtod` and only rejects
when nothing is consumed (`stringtools.cpp:316-329`), so `inf K` parsed, survived the
parser's `!isnan` gate, and passed `> 0.0` — leaving a permanently superconducting material.
Guard is `std::isfinite( tTcrit ) && tTcrit > 0.0`.

Both audits also corrected a claim of ours that had gone into the plan: **`T_crit` does not
gate "every power law".** `rho_powerlaw` (`powerlaws.hpp:181-349`) has no temperature gate at
all, and `PowerLaw` is the DEFAULT law. Only `piecewise` and `riva` read it. Both contract
documents say so now.

## 4. Test results — 11 of 12, and one incident

Ran against Christian's build: neutral (90 K = the plugin's own value) bit-identical to
no-key over 128 Picard lines; 50 K differs and decays at a flat 0.9 ratio, i.e. a LINEAR
problem with the power law switched off; all three builtin forms abort; a `usermat` material
without jc reaches and trips the second guard; no unit / `inf K` / `92.5K` / `-5 K` / `0 K`
all abort with the right message. Not run: the key on a `curve` ferromagnet (no ferro
material in this deck).

**The incident worth keeping.** Tests 1 and 8 initially DISAGREED, reproducibly, against a
solver proven bit-deterministic by a repeat run. Bisecting the value showed only *exactly*
90.0 differed — 89.9, 90.1 and 200 all matched the no-key case, which is backwards for any
threshold and therefore could not be physics. It was a **stale `usermat.so`**: after
`touch matlib.cpp && make` with no source change, the neutral test passed bit-identically.
**The `.so` mtime was NEWER than the source**, so mtime was not a usable staleness signal —
the failure mode in `project_user_material_plugin_abi_fragile`, with its usual detector
removed.

## 5. `custom` → `usermat`

Christian's call, after his own observation that the HTS cluster (`usermat`, `defect`,
`resistivity type`) reads inconsistently, and his rationale: fool-proof naming.

Both auditors had recommended KEEPING `custom` under the present contract, and both rejected
`resistivity` — it names about a sixth of what the subsection supplies (the plugin's
`hts_init` registers jc, n, ec, rho, cp and lambda) and collides with the existing
`resistivity type` key, which selects the E-J law. Codex and Grok both offered `plugin` as
the rename target IF a per-property-sourcing grammar is ever built.

What decided it: **no shipped deck uses the shape at all.** All eleven `examples/*/input.conf`
were grepped — the only "custom" hits are `libcustom.so` file paths. The only live user is the
unshipped tape_quench deck. Renaming before the release is free; after it, it is breaking
forever. Hard rename, no alias.

Consistency rule adopted for the three plugin grammars: **a "user" marker appears where a
builtin alternative exists; a role name suffices where it does not.** So `builtin : ybco` vs
`usermat { }` (fork, both marked), `defect { }` unchanged (no builtin defects exist), and
`type : userdefined` unchanged (it is one value of a waveform enum, implemented in three
factories, and shares the stem).

Also replaced the S7 message: `"Custom materials created through input file not implemented
yet"` described an unimplemented feature rather than the user's actual mistake, which is a
missing `file` key. New message leads with the mistake and keeps the feature statement —
inline properties genuinely are unimplemented.

## 6. Plugin restructure

Sources to `./src`, data to `./data` (Christian), then: eleven near-identical 68-line table
loaders collapsed into one shared `usermat::Table` (`matlib.cpp` 1027 → 403 lines,
`current.cpp` 112 → 51), `throw std::runtime_error` → `BELFEM_ERROR`, two CMakeLists → one
building `usermat.so` and `userdefect.so`.

The `__FILE__` trick the old loaders used to find their data had ALREADY broken when the
sources moved — it bakes the source directory into the binary, and the data had moved
elsewhere. Replaced with a compiled-in data directory, which the plugin audit then flagged
as non-relocatable; Christian's ruling is that `$BELFEM_DATA` is for BELFEM's own hdf5
databases and user txt tables belong to the deck, so the compiled default stays for now.

**The interpolation rewrite was proven, not assumed.** `lower_bound` → `upper_bound` over
1.2M probes across all 13 data files, every node exactly and at ±1e-12 / ±1e-7: 11
disagreements, all at exact node hits, max 1.1e-16 — and in every one the NEW code returns
the exact tabulated value while the old computed `Y[i-1] + 1.0*(Y[i]-Y[i-1])`, one ulp off.

## 7. What was found and NOT fixed

Both filed rather than folded in, per the no-scope-creep rule:

- **DR-146 `[F]`** — thirteen silent-ignore cases in the `materials` section. One is live on
  this very deck: `file : sst-1.hdf5` beside a `usermat { }` subsection is never opened.
  Grok refuted that liveness reading `cmake-build-claude/tape_quench_dr105`, where the line
  is commented out — the wrong deck; the claim stands against
  `cmake-build-debug/tape_quench_usermat`. Blast radius MEASURED independently by Claude and
  Codex with identical results: three `examples/` decks fail a raw sweep, ONE fails with a
  losing-selector exemption. Both auditors rejected the simple consumed-set sweep in favour
  of an explicit shape-aware validator, on Grok's structural argument: `RRR` on ybco IS read
  and then discarded by `create_material`, so "was it read" and "did it do anything" are
  different predicates and only the second is the contract.
- **DR-147 `[W]`** — the three plugin `.so` paths resolve differently; the source-function one
  bypasses `material::data_file()`. The one-line fix is blocked by layering (numerics sits
  below physics), which is the whole content of the row.

Also noted, unfiled: the tape_quench deck's `file : sst-1.hdf5` is asking for something
reasonable and currently impossible — jc/n from a measured table, everything else from a
plugin. Per-property sourcing is a design question, not a bug, and it is the one condition
under which both auditors would revisit the naming.

## 8. Status

Code, both contract documents and the deck landed; `check_doc_claims.py` 36/36;
`-fsyntax-only` clean. The T_crit key is **verified** by execution (11 of 12 tests).
The rename is **reviewed, not verified** — it needs a build and a deck rerun.
Nothing committed; the tree carries a concurrent session's in-flight work throughout.
