# Rutherford cable generator: a tape wound around a flat core

**Date:** 2026-09-18
**Purpose:** New example deck `examples/rutherford` with a Python mesh generator for a REBCO tape wound helically around a flat Rutherford-type core, modeled on the `corc_twolayer` generator.
**Status:** generator built and gated by its own checks (four configurations); the BELFEM deck has not been run.

---

## Task

Christian asked for a generator like `examples/corc_twolayer/python/main.py` for the cable of
Emily Romancew's thesis (`tmp/Development_and_Characterizati.pdf`, NIU 2026): a REBCO tape of
2, 4 or 6 mm width wound at a lay angle of 30 or 45 degrees around a flat Rutherford-type
copper core with rounded edges. His geometry description: a stadium cross-section of width
w, a bending radius r at the edges, a flat length ws = w − 2r, the lay angle measured from
the winding axis, periodic planes normal to the axis.

## What the thesis fixes and what it does not

Read from the text and Figures 3.3, 3.5, 3.8 (high confidence):

- one tape on a copper Rutherford-type core, not a cable of REBCO strands; tape widths 2, 4,
  6 mm; lay angles 30° and 45° **from the cable axis**; the thesis pitch rule
  `L = 2 w / tan θ` (Figure 3.8) uses the visible face width and neglects the round edges
- the 30 to 40 mm "bending diameters" are the mandrel the finished cable is bent around, not
  the edge radius
- **not in the thesis text at all:** core width, core thickness, edge radius. The ANSYS core
  is a "rectangular midsurface". The 1.8 mm thickness (r = 0.9) comes from Christian's
  memory of Emily's talk; the 12 mm width is a placeholder he approved

## Geometry as built

Cable axis z (the CORC convention: periodic planes z = 0 and z = L, parallel to x-y).
Stadium in x-y, width along x, thickness 2r along y, perimeter `P = 2 ws + 2π r`. The
extruded core is developable, so a constant-lay-angle tape is a straight band in the unrolled
(s, z) plane with `ds/dz = tan θ`; pitch `L = P / tan θ`; the tape occupies `wt / cos θ` of
perimeter at fixed z; no self-overlap iff `wt ≤ P cos θ`. On the round edges the tape crosses
the curvature obliquely: normal curvature `sin² θ / r`, geodesic torsion `sin θ cos θ / r`
(Euler's theorem). At 30° the effective bending radius is 3.6 mm for r = 0.9, at 45° it is
1.8 mm; the generator prints both.

Decisions taken with Christian in the session:

- **copper core not modeled**, its section is air block 1 (switching it to a conductor later
  is a retag of that block)
- **one layer for now**; multi-layer support is designed in: layer l sits on the stadium
  with r + l·gap, counter-wound (criss-cross), **at the common pitch**, so its lay angle
  differs (`tan θ_l = P_l / L`) exactly the way the CORC generator lets the lay angle vary
  per layer. Without this there is no common period
- **periodicity is a pure translation** (every tape returns to its arc length after one
  pitch). Fractional pitches close when `numPitches · numTapes` is whole, with the tapes
  permuted, the CORC `tape_shift` rule. The half-pitch closure by a 180° rotation, which the
  stadium's point symmetry would allow, is refused with a message, not implemented
- **rows at constant z.** A strip quad has one side along the perimeter and one along the
  tape, so the cut planes are mesh rows and never slice an element. With the along-tape step
  equal to the across-tape slot, the parallelogram angle is 90° − θ and the short-diagonal
  split gives equilateral triangles at 30° (45/67/67 at 45°). The flat/round junctions cross
  these elements; the chord error is micrometers
- **caps as polylines** through the exact ring node positions (one gmsh `Point` per slot,
  two-node `Line`s), instead of the CORC's transfinite circle arcs, because the slots do not
  land on the flat/round junctions. Domain circle unchanged (quarter arcs)

## Package

`examples/rutherford/python/rutherford/`: `stadium.py` (map, normals, signed distance),
`winding.py` (schedule, closure, counts, guards), `shells.py`, `caps.py`, `assembler.py`
(CORC's, scratch names changed), `postprocess.py` (CORC's minus the frame; the vertex
transform must come out as the translation `(0, 0, L)`), `belfem.py` (topology writer,
translation only, no solder branch), `cable.py`, `builder.py`/`mshio.py` (CORC's). The
`python/mesh/` package is copied unchanged. Physical tags, the ten type-15 vertices and the
`$Periodic` block follow the CORC output exactly, so the deck's `topology` block is the
CORC one with the tape and plane ids changed.

One defect on the first run: the strips came out inward-facing because the arc length runs
clockwise (seen from +z) where the CORC's angle runs counter-clockwise; the two orientation
branches were swapped and the outwardness assertion now passes.

## Gates run (in the scratchpad, gmsh 4.15.2)

| configuration | nodes | tets | γ min / p5 / median | cap pairs, residual | wall |
|---|---|---|---|---|---|
| default: 4 mm, 30°, 1 tape, 1 pitch | 24 286 | 130 826 | 0.179 / 0.572 / 0.844 | 1643, 3.6e-15 | 31 s |
| 45° | 16 767 | 97 747 | reported OK | 1613 | 22 s |
| 3 tapes, 1/3 pitch (permutation closure) | 14 966 | 87 408 | reported OK | 2008 | 18 s |
| 2 layers, 0.2 mm, criss-cross (layer 1 at 31.18°, sense −1) | 97 645 | 562 066 | 0.184 / 0.572 / 0.816 | 2463, 4.0e-15 | 147 s |

Unit tests `test_stadium.py` (landmarks, arc-length parameterization to 1e-6, outward unit
normals, wrap, round-core degenerate case) and `test_winding.py` (30 checks: pitch, spans,
closure rule, overlap and geometry guards, layer-specific lay angle, thin-slab rule) pass.
The tape geometry was inspected in projection: the band crosses the top face at the lay
angle, wraps the round edge, returns on the bottom face and lands on its start after one
pitch.

**Not done:** the BELFEM deck was not run (Christian runs the binaries); no Codex/Grok round
on the code (Christian did not ask for one this session); no language sweep on the README.

## Default changed to five tapes

On Christian's call the shipped default is five 4 mm tapes at 30° over a fifth of a pitch,
the most tapes that fit (`n ≤ P cos θ / wt` = 5.6), 88.6 % coverage, 0.59 mm of bare core
between neighbors. Every tape continues into its neighbor across the periodic planes. The
thesis single-tape sample is `numTapes = 1`, kept as an end-to-end variant.

| configuration | nodes | tets | γ min / p5 / median | cap pairs, residual | wall |
|---|---|---|---|---|---|
| default: 5 × 4 mm, 30°, 1/5 pitch, 9.0 mm long | 15 073 | 74 234 | 0.342 / 0.668 / 0.854 | 2446, 3.6e-15 | 20 s |

## First deck run: "5 conditions, 2 cohomology generators"

Christian ran the five-tape deck in `build/rutherford` and BELFEM refused at setup: five
current conditions against two generators. Correct refusal. With a fifth of a pitch tape k
leaves the back plane where tape k + 1 enters the front plane, so the five sidesets are one
tape of five times the length, one conductor, one condition (the second generator is the
free axial loop of the periodic model). The deck had inherited the CORC writer's one bracket
per tape. Fix: `Winding.conductor_groups()` splits a layer into gcd(q, n) conductors
(q = tape positions advanced, n = tapes) and the topology writer puts each conductor's
curves into one bracket, `[1,3,5,7,9]`; the amplitude is the total through the group, the
tape current times the crossings, so 30 A per tape is written as 150 A. `corc_periodic_bc`
(three tapes, a third of a pitch) uses the same single-bracket form. Deck and README
updated, tests extended. Whether the run then proceeds is Christian's gate; not seen here.

Observation, low confidence, not acted on: `corc_twolayer/input.conf` has six separate
brackets although its own `main.py` says every tape's far end is the next tape's near end,
which would be the same refusal; either the twisted closure counts generators differently
or that deck has not been run in this form.

## Second deck run: generator 0 is a thick cut

With the current block fixed BELFEM refused at the cut stage: "generator 0 is a thick cut,
no representative with only −1, 0, 1 coefficients on this mesh", naming a loop of two edges
between nodes 3518 → 10516 → 5813 near element 47306. In the mesh those nodes sit at a
radius of 40 mm, in the outer air far from the tapes (r ≤ 6 mm): 3518/5813 are a periodic
pair on the end planes and 10516 an interior node bridging the 9 mm period in two 4 mm
edges (738 such loops, all at r = 26 to 47 mm). My reading: the SNF representative of
generator 0 carries a multiple of the axial class, and the short period makes the shortest
axial loop too short to realize it. Test: the same deck with `domainResolution = 1.5` has no
two-edge periodic loop. **Christian ran it: still refused.** The reading is therefore
refuted as the whole story, or the multiple is larger than the loop length at 1.5 mm; the
evidence for the location of the flagged loop stands, the explanation does not. Christian is
working on thick cuts on the solver side and keeps the deck as the test case. Deck default
left at 4 mm.

## First field solution (Christian's thick-cut work in progress)

Christian ran the fine-air mesh through his work-in-progress thick-cut support and showed
the J/Jc plot (`build/rutherford/rutherford.png`), asking whether the fields of low current
are realistic. Read from the eleven Exodus steps (t = 0.1 to 1.1 s), YBCO layer, Jz
integrated across all five crossings on the plane z = L/2:

| t (s) | I_1 imposed (A) | YBCO current through the plane (A) | ratio | J/Jc min / mean / max |
|---|---|---|---|---|
| 0.1 | 1.65 | 41.0 | 24.85 | 0.00 / 0.13 / 0.33 |
| 0.4 | 6.15 | 152.9 | 24.86 | 0.00 / 0.49 / 0.92 |
| 0.8 | 12.15 | 301.8 | 24.84 | 0.69 / 0.95 / 0.99 |
| 1.1 | 16.65 | 398.9 | 23.96 | 1.23 / 1.25 / 1.26 |

Every crossing carries 5 × I_1, the whole tape 25 × I_1: the imposed generator current is
realized with the winding multiplicity of the thick generator (five crossings of one
conductor). The reported voltage shows the same factor: at 1.1 s the power law gives
E = 1e-4 (1.252)^20 = 8.9e-3 V/m over 52 mm of tape per period, 4.6e-4 V, and U_1 = 2.4e-3 V
is five times that. Consistent picture; the mechanism sits in the cut, not the deck.

The plot itself matches step 4 (I/Ic ≈ 0.48 per crossing with Ic = 64 A): Jc at the tape
edges and about 0.3 Jc in the tape centers is the Brandt–Indenbom / Norris transport profile
at that current (center value (2/π) atan(sqrt(w²/a² − 1)) = 0.31 Jc for a/w = 0.88), so the
low-current fields are the flux-free strip core, shaped by the neighbouring tapes and the
bends. From t = 0.9 s the tape is above Ic everywhere (J/Jc 1.06 → 1.25 uniform) because
80 A per crossing exceeds Ic. Interface nodes of the Exodus layer stack carry the neighbour's
value (silver shows J/Jc up to 1.26), so per-layer currents must be read from the YBCO block
alone.

**Retraction, same session.** Christian found the saturation pattern non-intuitive. Unrolling
the YBCO layer into the (s, z) plane (step 4, 30.8 A per crossing) shows why: every tape has
the same ONE-SIDED profile, J/Jc ≈ 0.4 at its low-s edge rising monotonically to 0.9 at its
high-s edge, where the Norris strip solution is symmetric with both edges at Jc. The
handedness is the same all around the perimeter (mid top face, mid bottom face, bends) and
uniform in z. A one-sided profile needs a perpendicular field; the cable's own perpendicular
component (≈ 3 mT at a quarter width against the tape's critical-state scale of 6 mT) flips
sign across the centerline and between faces, and a one-signed normal field on a closed
surface would carry net flux (div B ≠ 0). The helix's only chiral field is the solenoidal
B_z (≈ 4 mT inside the core, µT return outside), in-plane for the tapes. So the "realistic"
verdict above is withdrawn: the unevenness is not explained by the physics of the model and
points at something chiral in the discretization (thick-cut realization, or the one-sided
layer stack of the thin shell on a helix). The cut files (`cut_0.vtk`, `cut_1.vtk`) are
whole-domain face sets (150k faces, spread evenly across the tape bands), not ribbons, so
they do not locate it. Discriminating run prepared for Christian: one tape, one full pitch,
fine air (scratchpad `onetape/`): no neighbours, weak cable field; if the one-sided profile
persists it is intrinsic to the model.

## Three-voice investigation of the first field solution (evening, Christian away)

Christian's three questions: (1) a physical explanation for the saturation pattern, (2) are
the boundary conditions right, (3) is the post-processing right. Round: Claude pre-registered
(`tmp/ai_exchange/review_rutherford_saturation.md`), then a blind jury on a facts-only brief
(`rutherford_saturation_brief.md`): Codex `gpt-6-astra`/xhigh, Grok `grok-4.6`/xhigh, plus
two reproducer runs by Claude in the scratchpad with the devel binary. Figures and the
analysis script: `tmp/ai_exchange/rutherford_figs/`.

**Provenance corrected.** The run in `build/rutherford/` used the 1.5 mm-air mesh and the
devel binary (`build/bin/belfem`, 18:25): thin cuts were built and eleven steps saved. The
README's "refinement does not change that" was wrong and is corrected; the 4 mm default is
still refused. Christian restarted from `memdump.hdf5` at 22:26, which truncated
`iv_results.csv` to one row (t = 1.2 s) and overwrote step 11.

**Q2 — P0, unanimous, reproducer-level.** Each of the five crossings carries 5 × the imposed
I_1 and the cable 25 × (Ampère on the nodal H, Codex's exact −∇φ integral −153.750 A at
I_1 = 6.15 A, the YBCO plane integral, the solenoid field µ0 × 30.75 A / 9.03 mm = 4.27 mT
inside the core, all agree). Mechanism: the five curves of one bracket become ONE suggested
homology chain, the sum of five loops (`cl_Homology.cpp:377-386`); the coupling matrix then
has the invariant factor 5, and `updatekGeneratorsFromHomology` builds the constrained
generator as `D R Σᵀ Q⁻¹` (`cl_Cohomology.cpp:1123, 1141-1143`), which multiplies by a
non-unit factor where an inverse would divide; the unit-factor count `s` that `smithForm`
returns is discarded, only the rank `t` is checked. Codex reproduced the algebra
(F = [[5, −5],[0, 0]], D = [1, −1] → A = [5, 0], AF = [25, −25]). Control: the single-tape
deck (one curve per group) gives 1.000 × I_1 per crossing on every step, 0.33 A to 30 A.
The same mechanism explains the 4 mm refusal: a generator scaled by five needs five
disjoint cut sheets, which a two-edge periodic loop cannot host. Deck topology, translated
periodic triple, sign convention and the natural outer boundary are all correct. **Design
ruling owed by Christian:** what a k-crossing group's amplitude means; the fix is in the
closed core (Gregory). Until then the five-tape deck drives 5 × 30 A per tape.

**Q1 — physics, all three voices agree, not voted.** The uneven J/Jc is the cable-scale
critical state at the MEASURED I/Ic ≈ 0.48: every tape saturates toward the nearest narrow
side, the across-tape asymmetry flips sign across the centerline on both faces
(±0.31 at |x| > 3 mm, 0 at the center, 0 on the bends), and the mean level follows a Norris
strip of the full 12 mm width. Isolated-tape control at 30 A (0.3 s ramp): symmetric to 0.03,
center 0.35 vs Norris 0.31, edges 0.75 to 0.78 against the analytic plateau (the saturated
zone is 0.24 mm wide against 0.4 mm nodes). Grok's qualifications adopted: "saturation" is
loose at step 4 (max 0.92, n = 20), the peaks sit on the flats beside the bends (bend bin
0.66), and the whole-cable strip is a foil for a hollow sheath, not a prediction. My two
earlier verdicts in this file ("realistic" before the unroll, then "fixed handedness, not
physical") are both superseded by the binned asymmetry table and the control run.

**Q3 — consistent within documented limits, three findings.** YBCO J integrates to the
field's current to 1 % and its across-width shape matches the H jump across the shell;
J/Jc = |J| / 1e10 (constant deck jc). (i) Nodal `J`, not only `J/Jc`, is overwritten at the
shared layer-interface nodes: the normal-conductor postprocessor runs first, the
superconductor one second, both into one global field (`cl_MaxwellFactory.cpp:2528-2545`,
`cl_FEM_Postprocessor.cpp:163`); Codex's probe shows 95 A "in" the 2 µm silver. Read
currents from the ybco block only. Severity P1 (Codex P0), to Christian. (ii) `save_IV`
differentiates with `(x_n − x_{n−1})/Δt` while the solve is BDF5 (`cl_FEM_Controller.cpp:4123`
vs `cl_IWG_Timestep.cpp:995-1011`); magnitude open. (iii) `iv_results.csv` truncates on the
first save of a process (`:4155`). Also: `I_1`/`U_1` are the abstract-dof value and its
conjugate, positional names; Codex reconstructed `I_2` as the axial circulation coordinate.
The brief's own defects (literature pointer, line range, Exodus vs `.msh` counts, "milliamps"
at step 11, tautological B = µ0 H) are corrected in the exchange record.

**Runs launched by Claude this session** (in the scratchpad, not in `build/`): single tape,
one pitch, 1.5 mm air, 30 A over 10 s to 0.35 s; the same with a 0.3 s ramp to 30 A. Both
with `build/bin/belfem` of 18:25. Christian did not ask for these runs; "as much as you can
in autopilot" was read as covering a reproducer that decides the P0.

## Follow-up on branch `Claude` (2026-09-19, night)

On Christian's instruction the fix was tested on a quarantined branch `Claude` (worktree
`../belfem-claude`, never merged): the constrained generator is built from the exact rational
pseudo-inverse row scaled by the least common multiple of its denominators (unit generator,
factor k kept) and the cut jump is weighted by 1/k. Five-tape deck: 30.0 A per crossing at 150 A
on the 1.5 mm and the 4 mm mesh (the latter refused on devel); single tape, `tapestack3d` and
`corc_periodic_bc` identical to devel at solver tolerance; `make check-fast` 8/8; plan and code
juries (astra/high, Grok/high) folded in. Records on the branch:
`todo/smith_inverse_grouped_terminals.md`, `devlog/dl20260919_smith_inverse_branch.md`, and the
report `devlog/dl20260919_smith_inverse_report_for_gregory.md`. Gregory decides.

## Follow-up 2: Gregory's reading, tested (2026-09-20)

Gregory (relayed by Christian): the unit Smith diagonal is a design assumption; with translation
periodicity all tapes form one connected conductor, so ONE input/output terminal on any one
conductor is to be defined, as in his corc 1/3-twist example. Tested on the five-tape 4 mm deck
with `input curves : [1] ; output curves : [2] ; amplitude : 30 A` on devel `26bfd28` and on the
branch binary, fresh cut stage: devel passes the cut stage (no thick-cut refusal), no jump-weight
line on the branch (k = 1); per crossing at z = L/2 10.00, 9.84, 10.00, 9.83, 10.00 A at I_1 = 10 A
and 29.99, 29.52, 29.99, 29.50, 29.99 A at 30 A on both, Ampère 4.998 × I_1 around the cable
(five crossings), fields within 1.2e-7. Gregory's reading holds: with one crossing per condition
the formula is exact; the 5 × of Q2 came from the deck's one bracket holding all five crossings,
and the 4 mm thick-cut refusal was a consequence of that k = 5 generator, not of the mesh. What
led the deck there is the count-guard message (`cl_Cohomology.cpp:1088-1101`, "merge the
terminals that belong to the same conductor into a single bracketed group") and
`doc/input_file_reference.md` §"Lists, ranges, groups" ("electrically connected and act as a
single terminal"), both true of the five crossings of a chained periodic conductor. The deck,
its README, `belfem.py` and the e2e test now use the single-crossing convention; the branch
patch stays on the shelf as option (b) in the report's addendum (§0), option (a) being a guard
on `s < t` plus corrected message and docs. Gregory decides.

## Registrations

`examples/README.md` (count and feature-matrix row), `doc/getting_started.md`,
`doc/doxygen_nav.dox`, `CHANGELOG.md`.
