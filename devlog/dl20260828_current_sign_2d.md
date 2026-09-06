# Devlog 2026-08-28 — 2-D current boundary condition ran with an inverted sign

**Date:** 2026-08-28
**Topic:** A declared positive `current` amplitude on a 2-D deck imposed the
correct magnitude with the sign reversed, on both bulk and thin-shell
conductors. 3-D was unaffected. Root-caused, fixed, and gated.
**AIs involved:** Claude, Codex (plan + code), Grok (code round only —
its plan round died on the turn budget)
**Claude Confidence:** high
**Codex Audit Confidence:** high (plan and code rounds)
**Literature References:** none required — this is a discrete orientation
convention, not a formulation choice. h = -grad(phi) and J = curl(h) follow
Arsenault et al. 2023 as already implemented.
**Verification:** end-to-end reproducer, four gates on the rebuilt `belfem`
(cmake-build-claude, fresh cut generation) — cos-theta bulk Ampere ratio
-0.988 -> +1.00; gantry mixed polarity 464/464 tapes at their declared sign;
corc 3-D no-op (-0.6851 -> -0.6852); `make check-fast` 10/10

## Summary

Christian noticed that the gantry model reported `J/Jc` around -1.2 while every
tape had been given a positive current. The postprocessor was suspected. It was
not at fault: the inversion sat in the cut orientation that the current
boundary condition rides on, and it had been there for every 2-D model.

## Key findings

**1. The postprocessor is internally consistent — measured, not argued.**
On the pre-existing `gantry.exo` (464 tapes, +340 A declared, t = 13.5 s),
integrating the written `Jz` over all 9280 `EF_QUAD4TS` layer elements gives
exactly **-340.000 A on every tape**, 464 of 464 negative with no scatter.
Counter-clockwise Ampere rectangles drawn in air around single 58-tape coils,
using the written `H`, give -19717 / -19723 / -19722 A against +19720 expected.
`Jz` comes from `C(k) * dofs` on the layer element and air-side `H` from
`-grad(phi)`: independent paths agreeing to four digits. So the field itself
carried -340 A per tape; nothing in postprocessing flipped it.

**2. 3-D was already correct.** On `corc` (+200 A declared) the input terminal
sideset sits at z = 0.0377 and the output at z ~ 0, and the CCW circulation
about +z is -765 A: current flows from the input terminal to the output one,
which is the intended 3-D convention.

**3. 2-D bulk conductors were inverted too.** Answered without a new run from
`costhetatest` (68 bulk `coil` blocks, +169.96 A): CCW Ampere loops on the
air-side `H` converge to -11414 A against +11557 expected, ratio **-0.988**
over three nested boxes, uniform. Coil blocks carry no written `H`, so the
loop uses air-side field only; free-current Ampere is unaffected by the iron it
crosses. This ruled out a thin-shell-only repair.

**4. Root cause.** `Homology::reorient_generators` multiplied every suggested
1-generator by -1 in 3-D and by +1 in 2-D, carrying the comments
*"it should be +1.0 here, but there must be another sign error elsewhere in the
code in 3-D"* and *"somehow the sign is good in 2-D"*. Downstream,
`updatekGeneratorsFromHomology` normalises the cut cochain to
`< c_i, gamma_i > = +1`, a duplicate node ties as `phi_dup = phi_orig + I`
(`CutSet::create_duplicates`, weight +1), and the air uses `h = -grad(phi)`, so
a loop crossing the cut from the original to the duplicate side integrates to
`-I`. The generator must therefore run clockwise for a positive declared
current to come out along +z. The 2-D branch left it counter-clockwise.

The dimension split was never a property of the machinery below it. The one
dimension-dependent step there, the `mIs2D` flip in
`CutProcessor::check_edge`, was compensated by the split in
`reorient_generators`, so the composite from the *built* generator to the
duplicated side was the same map in both dimensions. What actually differs is
how `suggest_Homology` builds the generator: 3-D gets the difference of the
input and output terminal boundaries, 2-D the input loop alone, and the two
come out with opposite sense.

**5. Practical consequence.** Because the inversion was uniform, magnitudes,
losses and temperatures were never wrong and relative polarity within a model
was preserved — a single-polarity 2-D model was a correct solution of the
mirrored problem. Every signed current and field direction an older 2-D run
reported was reversed.

## Changes made

- `src/homology/cl_Homology.cpp` — `reorient_generators()` now multiplies every
  suggested 1-generator by `-1` unconditionally; the dimension branch is gone.
  3-D is a no-op (it already had -1). The comment records the measured evidence
  in place of the two self-doubting notes.
- `doc/input_file_reference.md` — new "Sign of a declared current" paragraph in
  the boundary-condition section: 2-D positive = +z out of plane, 3-D positive
  = input terminal to output terminal, per-condition polarity for mixed
  models, and the dated history. The convention had never been written down.
- `doc/input_schema.yaml` — the same contract under the `amplitude` key as
  `sign_convention_current` with `since` and `history`, anchored on
  `"reorient_generators"`.
- `todo/current_sign_2d_fix.md` (+ `todo/README.md`) — plan and gate tracker.

## Gates

All four ran on a rebuilt `belfem` in `cmake-build-claude`, each in a fresh
directory carrying only the deck, the mesh and the material files. **A cached
`.bfm` skips the whole cut stage**, so a rerun over an existing `.bfm` cannot
see this fix at all — anyone repeating these gates must regenerate.

| gate | before | after |
|---|---|---|
| **O1** cos-theta bulk, Ampere ratio, three nested boxes | -0.9878 / -0.9876 / -0.9881 | **+1.0006 / +0.9984 / +0.9918** |
| **O2** gantry, curves 1:232 declared -340 A | (all +) | **232/232 exact**, -3.740 A @ t=0.1 s and -7.140 A @ t=0.2 s |
| **O2** gantry, curves 233:464 declared +340 A | (all +) | **232/232 exact**, +3.740 / +7.140 A, **0 sign errors** |
| **O3** corc 3-D, Ampere ratio to declared current | -0.6851 | **-0.6852** (no-op) |
| **O4** `make check-fast` | — | **10/10**, homology 6.08 s |

O2 is worth more than a sign check. Each unbracketed curve id becomes its own
condition, so the gantry deck builds 464 separate current BCs, 464 cuts and 464
abstract dofs, and `IWG_Maxwell::set_currents` assigns values to them
positionally on a contract the code says it cannot re-check
(`cl_IWG_Maxwell.cpp:97`). Under the old all-positive deck that ordering was
untestable. A mixed-polarity deck tests it, and it holds 464/464 — at both
saved timesteps, with zero spread within each group (every tape reports its
group's declared ramp value to the printed digit), and the second step reached
it through a Newton solve (17 iterates to -73 dB) rather than Picard alone.

O3 pins the no-op numerically rather than by inspection: the ratio is the
cosine of the CORC lay angle, and it agrees to four significant figures at
operating points three decades apart (186 A pre-fix, 0.48 A post-fix).

O4 was expected to certify nothing about the sign — no test in the tree asserts
a physical current direction — and it did not. It rules out collateral damage.

## Method note

The plan was pre-registered on the blackboard before any edit, with a named
discriminator (does 2-D *bulk* invert too?), a decision rule mapping each
outcome to a specific one-line repair, and a falsifier (mixed per-conductor
signs would have meant an orientation-determinism defect, not a convention
one, and would have stopped the repair). The discriminator returned uniform
inversion, so the falsifier did not fire and Fix B was selected over Fix A.

Codex audited the plan and then the diff. Its accepted findings: the doc
wording "only mirrored, not wrong" understated a reversed field and was
rewritten; the replacement comment was over-long and one claim stronger than
its static evidence, and was cut; `operator*(-1.0)` became `operator*(-1)` to
match the `int` parameter. One finding was rejected with reason — the
"unrelated prose changes" it flagged in `doc/input_file_reference.md` predate
this session and belong to an earlier prose sweep still in the working tree.

Grok returned no verdict on the plan round: it exhausted its 30-turn budget
($0.27, 32k output tokens) without emitting a final message. Retried on a
deliberately narrowed code audit — four questions, named files, an explicit
turn budget — it delivered the round's most valuable finding by **refuting
Claude's explanation of the voltage behaviour while confirming the conclusion**
(see Open questions). It also found the missing "+V drives +z current"
statement in both input-contract artifacts, two stale comments in
`cl_CutFactory.cpp` still claiming the generators are reoriented "so that
terminals point inwards", and a pre-existing false comment at
`cl_FEM_Controller.cpp:501` that states the current/voltage ordering
backwards — the loop under it skips *current* conditions, not voltage ones.
All four were fixed; the last is by-catch, unrelated to this defect, repaired
because it sits on the exact path the voltage argument runs through.

## Open questions

- **G-V1 (owed).** Claude first argued that the change is the basis reversal
  `e_i -> -e_i` on an assembled system, so the reported I and V columns both
  negate and `I*V` is invariant by co-flipping. Codex refused that on
  assertion and **Grok refuted the mechanism outright** — correctly. The system
  is *reassembled* from the new embedding, and the declared scalars are applied
  in the new coordinates with no extra minus, so the reported columns do **not**
  flip: for the same deck `save_IV` prints the same `I_coord` and `V_coord` as
  before, and it is their geometric meaning that reverses. Coordinate power
  equals physical power because `(sigma I)(sigma V)` with `sigma = ±1`, so
  invariance holds *because the reported columns stay put*, not because they
  co-flip. The operational conclusion survives the correction: a positive
  declared voltage still drives a `+z` current, so there is no sign split
  between the voltage and current BC families, and none between the FEM branch
  and a netlist in a circuit-coupled 2-D model. Right product, wrong reason —
  the durable wording is *port redefinition*, not basis change on a frozen
  system. Still evidence level 5; a 2-D voltage-driven run is owed.
- **`orient_terminal_curves_2D` is geometrically vacuous.** It builds
  `tB = rot90(facet tangent)` and then `tR = z x tB`, which is identically
  `-(facet tangent)`. The "does this curve run counter-clockwise" test carries
  no information about which side the loop encloses; it reduces to "orient the
  curve antiparallel to the tape facet's node ordering", a mesh convention.
  The 3-D sibling uses two independent normals and is a real test. Uniform in
  practice here because facet master/slave is BFS-deterministic, but it is
  debt. Not touched in this change.
- No test in the tree asserts a physical current sign, so `make check-fast`
  guards against collateral damage but cannot certify this fix.

## Files updated

- src/homology/cl_Homology.cpp
- src/homology/cl_CutFactory.cpp (two stale comments)
- src/fem/kernel/cl_FEM_Controller.cpp (one false comment, by-catch)
- doc/input_file_reference.md
- doc/input_schema.yaml
- todo/current_sign_2d_fix.md
- todo/README.md
