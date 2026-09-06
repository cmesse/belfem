# CORC Solder: Cut Forensics, the Sideset-Only Interface Contract, and Gap Physicals

**Date:** 2026-08-25 (second arc; deck conversion is dl20260825_corc_solder_deck.md)
**Purpose:** Record the corc_solder homology/cut investigation, the three-AI jury
on partial thin-shell coverage, and the mesh-generator fix.
**Modules:** src/homology (probe only), corc_solder python generator
**Exchange threads:** `tmp/ai_exchange/corc_solder_cuts.md`,
`tmp/ai_exchange/corc_solder_coverage_jury.md`

## Constraint wiring: verified correct, both homology lineages

Christian's main-vs-claude cut overlay (neither touching the cable) triggered a
full forensic round. Probe (temporary printf in `updatekGeneratorsFromHomology`
of both lineages + one-shot in `IWG_Maxwell::set_currents`; still in the tree,
strip when campaign closes):

- Coupling matrix on the terminals deck: input cap (gmsh $Periodic slave) pairs
  zero — folded out of the complex; output cap pairs +1 with generator 2;
  D = [1,−1], Smith rank 1; transform puts −generator2 into slot 0;
  `set_currents` fixes exactly that dof. Census 672172/3/217099
  (= 1 current λ + 2 bearings; free-cut λs float, per the ruling).
- Both implementations (`homology` ↔ `homology_stable`, swapped via symlink)
  produce BIT-IDENTICAL matrices, transform, and census on this deck. The
  earlier 535563/7/257218 divergence belongs to the retired curves-BC variant
  (secondary defect hunt, downgraded).
- iv_results decode: I_0 = exact sigmoid value (echo of the imposed λ);
  I_1 = −0.761·I_0 on a free dof; I_2 ≈ 0. (Auditors: the r_mid/pitch reading
  of I_1 is NOT a proven identification — field-side Ampère gate still owed.)
- bfm cut measurement: three cuts anchored exactly at the solder radii
  (2.4486/2.1486 mm) — the "floating cuts" impression was a rendering-scale
  artifact — BUT connected-component analysis shows they are pairing-correct
  DUST: 622 facets/163 components, 1364/339, 16/8. The air-only corc baseline
  is worse (7 cuts, 11k–17k facets, thousands of components each).

## The cut-density mechanism (literature-anchored)

`mSuggestHomologies` is hardcoded true and `compute_cohomologies` skips
`reduce_complexPellikka` whenever suggestions are on ⇒ reduction NEVER runs.
Pellikka et al. 2013: "It is due to the reduction algorithms that the computed
bases have somewhat small support." Coreduce-only SNF ⇒ correct classes,
sprawling support. Same disease the 2026-07-16 census measured (~226k-edge
generator support, ~60% of the domain). Fix directions (Christian's call):
carry suggestion chains through the reduction via its chain maps; or pair
pre-reduction and transform; optimal-cut post-processing stays optional
(TU-condition caveats). The dof gap between lineages = correct seam-duplication
semantics (claude) × dense representatives — fixing density claws back most of
it without giving up seam correctness.

## Partial-coverage jury (Claude pre-registered + Codex + Grok, blind)

Unanimous root cause: **the sideset-only interface contract.** InterfaceCondAir
is only a label on existing mesh sidesets (`detect_sideset_types`,
cl_Topology.cpp:213-248); sidesets exist only for 2-D physicals
(GmshReader:919-986); h-φ hanging runs only over such sidesets
(MaxwellFactory:1404-1488); CutProcessor's φ boundary comes from the same list.
The corc generator dropped the gap strips ⇒ ~11% of the solder pipe's boundary
had NO facets: silently unwired H_t (spurious free jump), φ boundary
"tape-only", perturbed cut topology. The helix example works because its .geo
has no physicals at all — gmsh emits every surface and interfaces
auto-materialize.

Also established (refuting the session's own caution): conductor-mounted thin
shells are a fully coded path — conductor is always master by domain rank;
`hang_thinshell_edges_on_edges_bottom` / `_on_nodes_top`
(MaxwellFactory:1543-1611); "air on both faces" in the docs motivates signed
sidesets, not a precondition. Christian's partial-coverage hypothesis: right
in effect (the geometry exposes the contract), wrong in mechanism (no
shell-confusion branch). Workaround "cover gaps with a solder monolayer shell"
is fallback only (adds a real conductive film; double-counts conductance).

## Generator fix (landed)

`postprocess.py` now emits gap strips as per-layer physicals `Gap_1..Gap_nL`,
appended after DomainBack so all existing ids stay stable; `belfem.py` warns
not to declare them (auto-typing does the wiring); tool CLAUDE.md updated.
Parameter archaeology: the proven mesh was built with tapeResolution 0.25
(unique from its 5760 = 180×16 tape facet counts); the checked-in 0.2
deterministically produces a γ=1e-3 sliver under gmsh 4.15.2-git. Regenerated
mesh (341,275 tets, γ_min 0.040) validates; physicals 1-12 exactly match the
proven mesh, Gap_1 = 2160 / Gap_2 = 5400 facets match the winding arithmetic.
Installed in the deck 15:43.

## Open gates

- Solder-only deck on the new mesh (no thinshell/layers/curves blocks) — the
  jury's discriminating experiment for the interface contract.
- Tape deck rerun — watch the shell/gap junction at the helical side curves
  (closed-curve non-duplication; edge-orientation abort MaxwellFactory:1859).
- Field-side Ampère gate (∮H·dl or ∫J_z dA vs I(t)) — the correctness
  certificate for the current constraint.
- ~~Probe removal (three files) once the campaign closes.~~ Done same day
  (16:20): all PROBE blocks stripped from homology_claude, homology_stable
  and cl_IWG_Maxwell.cpp; syntax-checked against the build flags. The
  new-mesh tape-deck run (16:10 out.txt) had already banked the probe's
  final evidence: gap wiring moved ~18k dofs into the condensed pool
  (659774/3/235491) with the constraint census and anchoring unchanged.
- Grok tooling: ask_grok.sh errexit bug kills its retry loop and deletes the
  stderr evidence; two salvages this session ran the CLI directly.
