# corc_solder: strip diagonal fix in the mesh generator; D2 gate read on the first save

**Date:** 2026-09-02
**Purpose:** Record (1) the cause of the "degenerate tets" failure of the CORC mesh generator and its fix, (2) the executable read of the D2 gate (`dl20260901_corc_solder_hn_recovery.md`, "D2 applied") on the first save of the rerun.
**Module:** `examples/corc_solder/python` (generator, mirrored in `cmake-build-debug/corc_solder/python`); run files in `cmake-build-debug/corc_solder`

## 1. Generator: the long diagonal, not the element width

`main.py` (pitch 3, two turns, 0.22 mm inter-layer bed, `tapeResolution` 0.22) failed
`Postprocessor._check_tets` with `gamma < 1e-3`. Classified from the gmsh output of that run:
all 136 degenerate tets in Volume 2 (the inter-layer slab), every one with ≥ 3 nodes on a
tape/gap strip. Cause: `ShellGenerator.build` split every strip quad along a–d. The quad is a
parallelogram sheared by ω·ds along the helix, so a–d is the long diagonal for both winding
senses, and the strip triangles measured 21°/128° (b–c diagonal: 40°/86°). A 0.22 mm slab
bounded on both faces by 128° triangles cannot be filled with sound tets; the widths
themselves (0.21 mm angular, 0.26 mm helical) were fine.

Fix (`corc/shells.py`, six lines): split along b–c, orientation kept (`_assert_outward`
passes). **Verified** on the deck's parameters: gamma min/p5/median 0.000/0.219/0.668 →
0.212/0.655/0.855, "no ill-shaped tets" in all three volumes, Netgen illegal tets in the slab
15341 → 0, validation OK. All six `corc/tests/*.py` pass; `test_multilayer_e2e` passes its
validation for the first time (its docstring blames a gmsh 4.13 limitation; the else-branch
asks to be promoted to a hard pass — not done). The `main.py` comments recording 200 µm as the
"meshable floor" and pitch 2.15 as unmeshable were measured on the sheared triangulation and
deserve a retry (untested).

## 2. D2 gate on `corc.e-s.00001` (t = 0.02 s, mesh of 15:28, deck with terminals 14 → 13)

Read with scipy/netcdf from the save (191325 nodes, three `cut_*` sidesets). Coordinates in m,
cable length 25.1 mm (`numTurns` 4/3 in the regenerated deck).

| gate item | expected | measured |
|---|---|---|
| I₁ = I | ramp(0.02 s) = 0.0721 A | I₁ = 0.0721, I₂ = −0.0581, I₃ ≈ 0 (`iv_results.csv`) |
| φ winding at r = 5 mm = I at every z | 0.072 A | 2πr⟨H_θ⟩ = 0.000–0.002 A at 7 stations (200–244 air nodes each); r = 10 mm same; r = 30 mm −0.11 A at z = 2 mm, ≤ 0.02 elsewhere |
| bore winding 0 | 0 | 0 (φ ≡ 0 in the bore: pinned) |
| cut footprint runs the cell | one connected θ-cut through the cell, faces spread evenly over z | cut_1: 485 faces in air_03; faces per z-quarter **284 / 40 / 15 / 146**, 30 % on the periodic planes; coordinate-keyed connectivity 57 components, main body 191 faces; θ scattered in the middle slabs. cut_2: 717 faces, 439 / 35 / 20 / 223, 60 components. cut_3 (bore): 759 faces, 74 components, jump 0 |
| J/Jc not rim-only | current along the whole tape | YBCO \|J\| max per z-octant: 6.3e7, 2.2e6, 2.9e7, 2.3e6, 4e5, 5e5, 1.1e6, 3.8e7 A/m² — two orders of magnitude between the cap octants and mid-cell |
| φ jump on the cuts | I on the θ-cut | duplicate-node pairs on cut faces: cut_1 jumps 0.014 (991 pairs), 0.072 (295), 0.058 (41); cut_2 0.014 (1224), 0.058 (769), 0.072 (120); cut_3 0.0 (all). 0.014 = I₁ − I₂: cuts 1 and 2 overlap |

**Verdict:** the gate fails with D2 in place. The current condition is honoured as a φ jump
(I₁ = I to four digits) but on a θ-cut that does not run through the cell — 430 of its 485 faces
sit in the two end quarters next to the periodic planes — so Ampère's law does not hold in the field: H does not wind around
the cable anywhere along the cell, and the tape current sits at the caps. `cut.png` (patches
clustered near the caps) and `current.png` (hot spots at both caps only) are pictures of exactly this. Per the
D2 ruling in yesterday's devlog, the cause is on the cut-realization side (cohomology core,
Gregory), not the deck or the mesh. Evidence level: executable data read from the save; no
run launched here.

Not established: whether the dust cuts predate 2026-08-25 (no earlier corc_solder save on
record), and whether I₂ = −0.058 A is a physical azimuthal MMF or a by-product of the mixed
jumps.

Files: `tmp` scratch probes only; generator edit in both copies of `corc/shells.py`; run
files untouched.

## 3. Control: tapestack3d (`cmake-build-debug/tapestack3d/tapestack3d.e-s.00001`, t = 0.025 s)

Christian reports the deck healthy; the same probes agree. It also carries `periodic` and solder
end-face terminals, so neither feature alone breaks the cut. 2πr⟨H_θ⟩ around the stack =
0.21–0.23 A at r = 4, 8, 16 mm and five z stations against I₁ = 0.2094 A; I₂ = 1.6e-10. cut_1:
1150 faces, **one** connected sheet (1146 + 4) inside a 60° sector, 264–324 faces in every
z-quarter, 1 % on the z-planes. cut_2: 4480 faces, 32 components (largest 78 %), all θ,
carries I₂ ≈ 0.

**Method correction:** connectivity keyed on exodus node ids fragments *every* cut (285 and 877
"components" on the healthy deck), because cut duplicates carry different ids on the two sides
of a face edge. Key on coordinates. The §2 table was corrected accordingly; the earlier
"147/210/241 components, dust" reading in this file's first version was that artefact. The
z-distribution of the cut faces and the ring winding are the discriminators, and both still fail
on corc.

What differs on corc: the θ-cut has to leave the outer-air side of the layer-1 strips (tape
shells + bare solder gaps), the bore is a second air region with its own generator, and the
periodic map carries a 120° twist (`numTurns` 4/3); the 09-01 rerun (two turns, no twist)
failed the same way (G7), so the twist alone is not it.


## 4. Non-periodic control (`cmake-build-debug/corc_solder_control`, `corc.bfm` of 16:48)

Same mesh and deck with the `periodic` block deleted and the terminal ids bracketed. Read from
the `.bfm` (facets with physical 0 are the cut; face nodes = master ∩ slave element nodes after
folding `nodes/duplicates`):

- **Generator count is right:** one cut (`cut_1`, 242 facets), 8 pins (outer air, bore, six
  buffer sheets), 203 duplicate pairs. Outer air is an annular cylinder (H₁ = ℤ), the bore a
  solid cylinder (H₁ = 0) — one θ-generator, as found.
- **The cut is not a θ-cut.** Area 164 mm² against ~1450 mm² for any surface from the layer-1
  surface to the outer boundary; 22 coordinate-keyed components (largest 93); faces per
  z-quarter 15 / 38 / 72 / 117; θ scattered; 91 of its 94 free edges lie on the layer-1 surface,
  3 on a cap plane, **none on the outer boundary**; no cut node on the outer boundary. A
  relative 2-chain whose whole boundary sits on the conductor bounds a pocket of air against
  the conductor: it is a relative boundary, cohomologically trivial. A φ jump of I across it is a
  gauge shift of the pocket, which is exactly the periodic run's signature (I₁ = I honoured on
  the duplicates, no winding of H). Periodicity is thereby excluded as the cause, as Christian
  found from the picture.

Control on the healthy side: the tapestack3d cut (`tapestack3d.e-s.00001`) is one sheet in a 60°
sector with 264–324 faces in every z-quarter and reaches the box boundary.

### Conceptual differences, corc_solder_control vs tapestack3d (from the two decks and meshes)

| | tapestack3d | corc_solder_control |
|---|---|---|
| conductor cross-section | simply connected bar (7 solder blocks + 8 tape sheets + edge coatings) | annulus (one solder block); the tapes are its inner and outer skins |
| air | one region, box minus bar | two disconnected regions: outer annular cylinder and the bore inside the conductor (8 pins vs 9, the bore is its own φ component) |
| terminal sideset boundary | one loop, the stack envelope | two loops: the outer circle at r₁ and the inner circle at r₀ (the bore's rim) |
| deck curves | tape ∩ solder end face for the six inner tapes (edges interior to the terminal union), tape ∩ air end face for the two outer ones | all twelve are tape ∩ air cap (`12 @ k`, `11 @ k`) and lie on the terminal's boundary circles |
| tape/air interfaces | top and bottom sheets only, edges wrapped by coating blocks | every tape is a conductor/air shell, no coating; between tapes the solder itself faces the air (gaps), the tape edges lie in-plane with the exposed solder |
| geometry | straight sheets, cut plane meets them at right angles | helical sheets; the tape end curves are helical arcs on the caps |

`Homology::suggest_Homology` (`cl_Homology.cpp:204-330`, bulk branch) builds the suggested
generator as −(every curve that touches a flagged terminal edge) − ∂(terminal faces), then the
same for the output terminal. On corc every curve touches the terminal boundary circles, on
tapestack only the two outer tapes' curves do. Whether that changes the suggested cycle is a sign
question in the cohomology core (Gregory); not established here.

Discriminating experiments, cheapest first (decks only, none run here): (E1) control deck without
the `curves` block; (E2) pitch → large (straight tapes, same topology); (E3) bore filled with
solder (`conductor { blocks : 1, 2 }` plus the cap disks emitted as their own terminal
sidesets — a generator change), which makes the topology tapestack-like.

## 5. The thick cut is fine, the thin cut is eaten: `collect_facets` on a non-tight representative

Christian enabled `CutProcessor::save_debug_meshes()` (DEBUG-only call at `cl_CutFactory.cpp:470-472`)
on the non-periodic control: `cut_0.vtk` holds the thick cut (14447 cohomology edges) and the thin cut
(242 faces). Read with `tmp/ai_exchange/corc_thin_cut_collapse_scripts/*.py` against `corc.bfm`:

- Every crossed air tet is a 3-edge vertex cap (15082) or a 4-edge 2|2 diagonal (6230): the thick cut
  is a unit cocycle and passes `Cohomology::check()` (δθ = 0 on faces) and `determine_cut_case_3d`.
- It is not tight: ~21k crossed tets against ~2–5k for a half-plane; 4541 edges in the r = 2.3–2.7 mm
  shell over all 12 θ sectors, along the whole length; 1797 thick-cut nodes on the conductor surface.
  The band's vertices admit no consistent side labelling (2-colouring: 10604 conflicts).
- `collect_facets` (3D branch `cl_CutProcessor.cpp:411-599`) emits the conjugated face of positive
  vertex caps only (:429-434), self-cancels (:478-488), drops φ-boundary faces (:489-497), then runs a
  fixed-point loop (:525-570) that deletes any face with a non-exempt edge carrying < 2 flagged faces.
  For a tight cut the positive faces are exactly ∂U, U = tets touching the plus side, a closed set
  the loop leaves alone. Here they are not closed, and the loop eats inward from every hole for ~70
  iterations. Emulation: 7090 → 6797 → 6479 → 181 faces (other sign choice 7306 → 429, 176 of them
  among the real 242). The real survivors' 94 free edges: 91 on tape sidesets, 3 on a cap, 0 on the
  outer boundary — faces bounded by exempt edges only. Ampère then fails (§4).
- `Cohomology::check()` is not insufficient for what it claims; tightness is not a cohomological
  property. The unwritten contract is the consumer's: the representative must be locally a
  coboundary. Nothing tests it.
- The loop is old (Christian, 2025-02-09 `7c315189`; replaced by `manifold_filter_3d` in `6369747d`
  2025-08-07; restored in `34c44150` 2025-08-27 "to fix the physics"); documented as dangling-face
  trim (`todo/closed/2d_thinshell_gap_analysis.md:368-370`). Tapestack3d survives it.

**Jury round** (`tmp/ai_exchange/corc_thin_cut_collapse.md`; Codex gpt-5.6-sol high, Grok grok-4.6
high, blind, pre-registration in `_preregistration.md`; every citation re-opened, all confirmed; my
own prompt line numbers were 10–20 lines stale). Agreement on the mechanism. Corrections taken:
(1) the loss is spread over positive-only emission, self-cancel, boundary removal and pruning, with
pruning taking ~95 % of the remainder; (2) "closed sheath + sheet" is one reading of the census, a
multi-turn helicoid fits it too; (3) blaming `clean_spfa` is wrong — it rectifies to unit and preserves
support (`cl_Cohomology.cpp:342-358`), and `remove_cut_pockets` refuses any pocket touching a sideset
node by rule O1 (:668-682, :843-846), so a conductor-hugging support is kept **by design**; where it
is born (SNF basis, recombination `:1187-1198`, dense fallback) is not decidable statically;
(4) the guard must count holes (non-exempt edges of flagged-face degree 1 on `tTempCut` before the
loop), not a loss fraction — some trimming is in-spec; (5) disabling the loop is not a fix: a closed
sheath survives it as a T-junction bubble (Grok). By-catch (Codex): `cl_CutData.cpp:143` dereferences
`mThinCutFaces(0)` before any size check — an empty thin cut is a second failure mode.

Gates owed: dump tapestack3d's thick cut and run the same census (expected 0 conflicts, ~0 loss);
log the generator support size after each pass on corc; pairing of the dumped cochain against the
suggested cycle. Patch proposal (outside the core, not applied): hole count + `BELFEM_ERROR` with the
four sizes in `collect_facets`; size check before `cl_CutData.cpp:143`. Core question for Gregory:
which pass leaves the conductor-hugging support. Reviewed, not verified: no BELFEM run in this session.

## 6. Guard landed: `collect_facets` now aborts on a collapsed thin cut (approved by Christian)

Christian's ruling: the tight-representative problem is not solvable now (a Poisson-based node
reordering instead of RCM helped only a little); at least the solver must not run on a cut that is not
there. Edit in `src/homology/cl_CutProcessor.cpp` (outside the AI-closed core), 3D branch of
`collect_facets`: the number of faces entering the dangling-face loop and the number of passes are
recorded; after the loop `BELFEM_ERROR` fires if the thin cut is empty or if the loop removed more than
a quarter of the emitted faces, naming the cut, the counts and the pass count and pointing at the
thick/thin-cut document. The 2D branch gets the empty-cut error before its `set_size`. Threshold
rationale (jury §5): a healthy cut loses a few dangling faces in one or two passes, the corc collapse
removes ~95 % in ~70; a loss fraction is the only cheap measure that separates the two without
false-positiving the trim the loop exists for. Doc sentence updated in
`thick_thin_cuts_and_conjugate_edges.md`. Syntax-checked with the tree's compiler and flags; not
built, not run — the corc_solder_control rerun is the gate (expected: abort with "collapsed: removed
~6800 of ~7000 faces"). Diff audit (Codex terra/medium, Grok grok-4.6/medium, `tmp/ai_exchange/corc_thin_cut_guard.md`): pass
on tier, casts, placement and counter semantics; taken: wrap-safe comparison (`pruned <= emitted / 4`),
empty-cut message no longer blames pruning when nothing was emitted, doc wording "protective heuristic"
with the 2D empty abort named. Recorded risks: a small cut (a handful of faces) losing two danglers
would trip the quarter rule; no test reaches `CutProcessor`, so the decks are the gate. Christian's
same-evening `CutFactory::compute_poisson_problem()` (node reordering by a Poisson solve instead of
RCM, flag `mUsePoissonInsteadOfRcm`, default off) is in the working tree and untouched here.

## 7. Plan: `todo/cut_representative_options.md`

Christian asked for a todo explaining the three ways out and, in particular, how the cut-free
formulation (cocycle as edge basis function, Pellikka 2013 eq. 4.5) would look in BELFEM. Written with
Codex terra/high (consumer inventory of the thin cut, node duplication, abstract nodes, bfm, distributor,
postprocessing; finding: no mechanism attaches a carrier-free global dof to a volume element, the facet
lambda dof is the precedent) and Grok high (element Gram blocks, curl ψ = W₂(δz) = 0, circulation =
pairing, interface relation h_e = φ_{n0} − φ_{n1} + I z_e, both-sides shell risk, option 3 reduces to
option 2). Grok's option-1 counterexample (a closed cylinder as minimiser) refuted: pairing 0, not in
the class. Registered in `todo/README.md`. Exchange: `tmp/ai_exchange/cut_representative_options.md`.

## 8. Gregory's review of the plan (2026-09-03)

Gregory Giard, who owns `src/homology`, added two motivations for the cut-free formulation that the
draft had missed, both verified in the tree and recorded as §2.1 of the plan:

- **Integer coefficients.** The thin cut can only carry ±1 (`CutData` ± bitsets; the seven-pattern table
  in `cl_CutData.cpp:443-663`), which is why `clean_spfa` must rectify and, failing that, aborts with
  "no unit-coefficient representative … refine the mesh there and rerun" (`cl_Cohomology.cpp:534-543`).
  ψ = Σ z_e w_e is curl-free and correctly-circulating for any integers, so option 2 retires that abort
  and the open work in `todo/thin_cut_nonunit_rectification_implementation.md`.
- **Element shapes.** The cohomology engine is shape-agnostic (Gregory has run pyramids); the thin-cut
  conversion is tet-only and, worse, derives ONE element type for the whole mesh from dimension and max
  order (`cl_CutProcessor.cpp:33-35`, `:1166-1169`) with four cases and a "This should not happen"
  default (`:683`, `:802`, `:1220`) — so a pyramid or hex air element reaches the tet bitset routines by
  accident rather than meeting a clean check. Option 2 needs only Whitney edge functions; its reach would
  be bounded instead by the Calculator's Nédélec dispatch, which lacks PYRA5
  (`cl_FEM_Calculator.cpp:1100-1140`). Filed as O6.

**Standing decision (Gregory, agreed earlier with Frederic):** BELFEM should eventually implement the
thick cut instead, but not now — the majority of problems work and the priority is that the limits are
known, not that they are lifted. The thin-cut implementation was still worth having: the coefficient
problem had not been hit by anyone before, and both constructions are now understood well enough to
choose between them. The plan is a design record, not a scheduled item; the 2026-09-02 collapse guard is
what keeps the limit visible in the meantime.
