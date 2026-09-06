# CORC Cap Anomaly: From "Noise?" to the Corner-Pinch Mechanism

**Date:** 2026-08-25 (evening arc; continues dl20260825_corc_interface_contract.md)
**Purpose:** Record the identification of the periodic-cap current-density
anomaly's mechanism, and the elimination chain that got there.
**Modules:** src/homology, src/fem/maxwell, src/fem/kernel (probes only — no
behavioral source change)
**Exchange thread:** `tmp/ai_exchange/corc_cap_current_jury.md`

## The anomaly

Christian flagged depressed |J|/jc at the periodic surface of the corc_solder
run. Quantified: a ~1.6 mm two-sided band at the identified plane, ratio to
interior frozen to 4 decimals across a 68x current rise (stationary linear
artifact, not noise — Christian's "might disappear over time" refuted by
measurement). Sharpened to: ANTISYMMETRIC corner pattern — per tape, one
corner suppressed and one elevated (immediate cap row 2.4x / 0.23x), mirrored
between the caps, polarity following the winding sense, all three tapes of a
layer identical, decaying over ~5-6 element rows. Interior corners perfectly
symmetric (healthy Norris-type edge peaking). Christian called the t-matrix
suspicion; a false "gap surface current physics" reading (mine) was retracted
when fine binning showed the "gap elevation" was corner peaking measured at
corner-adjacent nodes.

## Elimination chain (each step measured or gate-level)

1. Terminal-curve orientation convention (orient_terminal_curves — cap
   normals mirror by construction, layer windings differ): decision table
   confirmed the mirroring, but an env-gated perturbation
   (BELFEM_PROBE_FLIP_BACK_TERMINALS) inverting back-cap decisions produced
   a bit-identical corner table => EXONERATED (also excludes anything keyed
   to curve start/end).
2. Edge-on-NODE constraint rows: 171,468 fused rows censused geometrically —
   max hang mismatch 47.8 µm = half the stack (expected); zero half-cut
   fused rows (policy holds); the 14 both-cap "mismatched" rows are cut
   lambda-crossings (sector-local, not 3-fold) => clean.
3. Edge-on-EDGE constraint rows: new probe (probe_edge_edge_row, 3 sites in
   DofData) dumped 53,355 ties; orientation-vs-weight census: ZERO
   violations => clean.
4. Periodic tie signs (+1 on aligned pairs), PART 2 bottom-edge tie
   (sign-aware, hard abort), suggestion algebra (pure signed permutation —
   probe-measured), ThinShellFactory cap pairing (layer-l to layer-l),
   postprocessor (anomaly decays over rows, not row-0-only) => all clean.

## Mechanism: THE CORNER PINCH (measured)

Side-curve policy leaves tape corner nodes unduplicated while the terminal
curve between them is duplicated => the sheet's phi-jump is pinched to zero
at the corner points. Pinching a uniform jump at both ends of a segment
yields two opposite local circulations: antisymmetric, cap-mirrored (sheet
current arrives vs leaves), sigma-selected, layer-scaled, lineage-shared,
immune to orientation probes. Confirmed in the phi field via the bfm
duplicate map: plateau jump +1.0426e-1 uniform (sheet healthy); corner
values collapse and flip sign, cap rows strongest (cap corner+ -5.2e-2 vs
interior corner+ +4.7e-2). Sampling sparse (n=1-4/bin) but structure
unambiguous.

Severity: bounded, characterized, few-percent in slab averages but 2.4x at
the immediate cap row — with n = 20 a potential spurious quench-localization
seed at the caps once current sharing starts; discount cap-adjacent quench
onset until fixed.

## Fix direction (awaiting Christian's ruling + full plan/audit round)

Extend the sheet duplication to the terminal corner nodes (let the tear
reach the corners) — revisiting the "closed side curves are not duplicated"
policy at side-curve x terminal-curve junctions. Hazard partner: the
unpaired one-sided-jump duplicate class (create_facet_map "Node not
flagged" / 5b-exclusion territory) — the periodic machinery must receive
the new corner duplicates. Gates: corner-ratio table symmetric, jump
plateau reaching the corners, interior untouched, iv_results identical,
check-fast, facet-map health check.

## Probes currently in tree (env-gated, inert by default; strip after fix)

- cl_CutFactory.cpp orient_terminal_curves: decision print
  (BELFEM_PROBE_TERMINAL_ORIENT) + back-cap flip
  (BELFEM_PROBE_FLIP_BACK_TERMINALS).
- cl_FEM_DofMgr_DofData.cpp: probe_edge_edge_row (3 call sites) extending
  the existing BELFEM_PROBE_FUSED_ROWS hook with edge-on-edge dumps.

## Process note

New standing rule recorded (memory: announce-output-analysis): declare
needed output files BEFORE a run; two analyses lost inputs to directory
cycling this evening before the rule.
