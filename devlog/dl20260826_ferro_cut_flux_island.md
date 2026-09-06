# Ferro-cut flux island found on the first completed gantry run

**Date:** 2026-08-26
**Topic:** DR-108 registered — iron-piercing cuts leave a flux-starved yoke arm; evidence
measured on the serial gantry run enabled by the 2026-08-25 thin-shell facet fixes.

## What happened

The first gantry run to ever complete (serial, 12 saves to t = 1.2 s, I = 42.011 A = 12%
of the 347.2 A ramp) produced a field plot with a dead upper yoke arm. Christian's read:
"as if the jump condition was not imposed along both sides of the cut."

Before that, a false alarm was resolved: the missing lower half of the magnet in the plot
is the y<0 mirror half that the deck deliberately omits (midplane symmetry, declared in
gantry.geo and input.conf; the mesh has zero nodes below y=0). The geometry is complete —
all 928 tape sidesets, all 932 surfaces, all 464 current BCs active. Global field sanity
passed: no NaN/Inf, Bz ≡ 0 (2D), far field decays to mT, aperture field ≈ 30 mT consistent
with the 12% ramp point. One watch item for Christian: 25% of iron nodes already above
1.5 T at 12% excitation (peak 2.97 T).

## The real finding (measured on hphi_results.e-s.00012)

- The cut path fragments the iron into six solver components (5109/1161/375/23/8/6 nodes);
  the small ones are cut-dust slivers, the same signature as the corc cut forensics.
- The dead component (375 nodes, upper arm): |B| ≡ 0.000 exactly, φ constant −38.0
  (spread 0.01) — a tied but flux-starved φ island, not shielding (machine zero, not mT
  leakage).
- The jump machinery works: every duplicate-pair interface carries a perfectly uniform
  Δφ, and the main one (194 pairs) carries 9747 A = exactly 4 stacks × 58 tapes × 42.011 A.
- The defect is the VALUES: the other fragment interfaces carry 8818 / 8780 / 10250 A —
  not clean tape-current multiples — i.e., inconsistent per-segment enclosed-current MMFs.
  Mixed jump signs along the trace (162/141) point at the probe-4b in-plane
  sign-coherence family.

Register row DR-108 (P1); investigation brief with the full evidence tables, root
questions, and literature routing: `todo/ferro_cut_flux_island.md`.

## Christian's mitigation idea (recorded, not implemented)

Input flag: may cuts run through a ferromagnetic domain? If forbidden, treat the yoke like
a coil during cut generation — its cut is always condensed out and the physical cuts route
around the iron. This geometry provably admits an iron-free layout (Christian solved it
manually two years ago). When implemented, the two-artifact input-contract rule applies.

## Session state carried over from 2026-08-25

DR-100: the "strays" were identified as the postprocessor field-count handshake (value 2,
from `FieldData::distribute`) arriving early from racing workers; the destructive probe
variant deadlocked the np=10 run and was reverted; handshake-tracing probes are in
`FieldData::distribute` and `Postprocessor::synch_node_indices`, awaiting the next np=10
run. Gates still owed: fresh-build gantry np=10 (both 08-25 fixes), check-fast.
